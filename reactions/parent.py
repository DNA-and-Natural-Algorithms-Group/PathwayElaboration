from __future__ import division
import warnings
import numpy as np
import ConfigParser
import math
import cPickle as pickle
import os
import sys
sys.path.insert(0,os.path.realpath('../learndnakinetics'))
import myenums
import timeit
import copy
from  scipy.sparse.linalg import *
from multistrand.options import Options, Literals
from multistrand.builder import Builder, BuilderRate ,  transitiontype,  localtype ,  hybridizationString, dissociationString, threewaybmString
from multistrand.concurrent import  FirstStepRate, FirstPassageRate, Bootstrap, MergeSim
from multistrand.system import SimSystem
from pathways import  PathwayHelix, PathwayHairpin ,  PathwayThreewayStrandDisplacement , PathwayBubbleClosing
from multistrand.objects import  StopCondition
from subprocess import Popen, PIPE, call


configParser = ConfigParser.ConfigParser()
configParser.readfp(open(r'../learndnakinetics/config_file.txt'))
CONFIG_NAME = 'parent'
NUPACK_bin= configParser.get(CONFIG_NAME, 'NUPACK_bin')
build_truncatedCTMC= bool(configParser.getint(CONFIG_NAME, 'build_truncatedCTMC'))
build_truncatedCTMC_pathwayelaboration = bool(configParser.getint(CONFIG_NAME,'build_truncatedCTMC_pathwayelaboration'))
build_truncatedCTMC_GillespieSSA= bool(configParser.getint(CONFIG_NAME, 'build_truncatedCTMC_GillespieSSA'))
use_Gillespie_MFPT = bool(configParser.getint(CONFIG_NAME, 'use_Gillespie_MFPT'))
use_Gillespie_Multithread = bool(configParser.getint(CONFIG_NAME, 'use_Gillespie_Multithread'))


CONFIG_NAMEplots = 'plots'

bimolecular_scaling = "bimolecular_scaling"
unimolecular_scaling = "unimolecular_scaling"

"""output object returned to learndnakinetics.py"""
class Output(object):
	
	def __init__(self, **kwargs):
		for k, v in kwargs.items():
			setattr(self, k, v)
	
	def set_specific ( self,  **kwargs) :
		for k, v in kwargs.items():
			setattr(self, k, v)
			
class ParentComplex(object):
	#Contains function and variables that different type of reaction have in common
	def __init__(self , dataset) :
		self.rate_method= dataset.rate_method
		if self.rate_method == Literals.arrhenius :
			theta = dataset.theta_simulation
			self.kinetic_parameters_simulation = {localtype.stack: (theta[0] ,theta[1]) ,
			                                  localtype.loop: (theta[2] ,theta[3]),
			                                  localtype.end: (theta[4] ,theta[5]),
			                                  localtype.stackloop: (theta[6] ,theta[7]),
			                                  localtype.stackend: (theta[8] ,theta[9]),
			                                  localtype.loopend: (theta[10] ,theta[11]),
			                                  localtype.stackstack:  (theta[12] ,theta[13]),
			                                  bimolecular_scaling : (theta[14]) }
		elif self.rate_method == Literals.metropolis :
			theta = dataset.theta_simulation
			self.kinetic_parameters_simulation ={ unimolecular_scaling :theta[0] , bimolecular_scaling :theta[1] }
		else:
			raise ValueError('Error: Please specify rate_method to be Arrhenius or Metropolis in the configuration file!')
		
		self.pathwayelaboration_N= dataset.pathwayelaboration_N
		self.pathwayelaboration_beta = dataset.pathwayelaboration_beta
		self.pathwayelaboration_K= dataset.pathwayelaboration_K
		self.pathwayelaboration_kappa= dataset.pathwayelaboration_kappa
		self.pathwayelaboration_use_elaboration = dataset.pathwayelaboration_use_elaboration
		self.simulation_mode = dataset.simulation_mode
		self.do_inference   = dataset.do_inference
		self.maxsolvetime = dataset.maxsolvetime
		self.strands_list = dataset.strands_list
		self.maxiter=  dataset.maxiter
		self.solveToggle = dataset.solveToggle
		self.deltaPruning= dataset.deltaPruning #either None or equal to a value
		self.dataset_type = dataset.dataset_type
		self.reaction_type =dataset.reaction_type
		self.dataset_name = dataset.dataset_name
		
		self.startStates = None
		self.statespacesize = None
		self.fatten_time=0
		self.build_time = 0
		self.load = dataset.load
		self.save =  dataset.save
		self.temperature_change  =dataset.temperature_change
		self.join_concentration_change = dataset.concentration_change
		#self.flur_position = dataset.flur_position
		self.real_rate = dataset.real_rate
		self.reaction_type  =dataset.reaction_type
		self.bimolecular_reaction = dataset.bimolecular_reaction
		self.dataset_path = dataset.dataset_path
		self.docID = dataset.docID
		self.temperature = dataset.temperature
		self.join_concentration = float(dataset.concentration ) # concentration has to be a float in Multistrand
		self.sodium = dataset.sodium
		self.magnesium = dataset.magnesium
		self.multistrand_time = 0
		self.rates ={}
		self.statespace = []
		self.repeated_states=  [ ]
		self.matrixsolve_time = 0
		
		self.num_simulations_Gillespie =  dataset.num_simulations_Gillespie
		self.simulation_time_Gillespie = dataset.simulation_time_Gillespie
		
		#some states about the building the truncated CTMC
		self.multistrand_simulation_time= 0
		self.loading_prototransitions_time  = 0
		self.loading_protoinitialstates_time = 0
		self.loading_protofinalstates_time  = 0
		self.loading_protospace_time  = 0
		self.multistrand_simulation_fatten_time= 0
		self.loading_prototransitions_fatten_time  = 0
		self.loading_protoinitialstates_fatten_time = 0
		self.loading_protofinalstates_fatten_time  = 0
		self.loading_protospace_fatten_time  = 0
		
		
		self.builderpath= self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+"/"+myenums.Permanent_Folder.MYBUILDER.value+str(self.docID)
		
		if  self.reaction_type == myenums.ReactionType.HELIXDISSOCIATION.value :
			pathway= PathwayHelix( False , self.strands_list , self.reaction_type, self.dataset_name, self.dataset_type,cutoff = dataset.cutoff )
		elif  self.reaction_type == myenums.ReactionType.HELIXASSOCIATION.value :
			pathway = PathwayHelix( True, self.strands_list , self.reaction_type, self.dataset_name, self.dataset_type,cutoff =dataset.cutoff )
		elif  self.reaction_type == myenums.ReactionType.HAIRPINCLOSING.value :
			pathway = PathwayHairpin( True , self.strands_list , self.reaction_type, self.dataset_name, self.dataset_type ,cutoff =dataset.cutoff)
		elif self.reaction_type == myenums.ReactionType.HAIRPINOPENING.value :
			pathway = PathwayHairpin( False , self.strands_list , self.reaction_type, self.dataset_name, self.dataset_type , cutoff =dataset.cutoff)
		elif  self.reaction_type == myenums.ReactionType.FOURWAYBRANCHMIGRATION.value :
			pathway = PathwayFourwayStrandExchange (  self.strands_list , self.reaction_type, self.dataset_name, self.dataset_type , lenx = dataset.lenx, lenm= dataset.lenm, lenn = dataset.lenn , ms = dataset.ms , ns = dataset.ns , loopbases= dataset.loopbases )
		elif  self.reaction_type == myenums.ReactionType.THREEWAYDISPLACEMENT.value :
			pathway= PathwayThreewayStrandDisplacement(  self.strands_list , self.reaction_type, self.dataset_name, self.dataset_type,  toehold_length = dataset.toehold_length , cutoff = dataset.cutoff )
		elif self.reaction_type == myenums.ReactionType.BUBBLECLOSING.value :
			pathway = PathwayBubbleClosing(  self.strands_list , self.reaction_type, self.dataset_name, self.dataset_type  , dataset.flur_position, cutoff =dataset.cutoff)
	
		self.startStates = pathway.get_pathway()
		self.uniqueIDtoLength, self.auto_strand_id_list = pathway.get_uniqueIDtoLengthandOrder()
		self.boltzmann_sample_initialstates_list= pathway.boltzmann_sample_initialstates(self.pathwayelaboration_N)
		#self.startStates =  pathway.get_intial_final()
	
	def find_meanfirstpassagetime(self):
		attributes_file  = self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+ "/atrributes" +str(self.docID)+".txt"
		attributes_file = open(attributes_file, 'a')
		attributes_file.write("concentration_change: "  + str(self.join_concentration_change) + "      concentration: " + str(self.join_concentration )+   "      temperature: " + str(self.temperature )+  "      temprature_change: " + str(self.temperature_change)+ "   sodium: " +str( self.sodium)+  "    magnesium: "+ str(self.magnesium) +"\n")
	
		if build_truncatedCTMC == True and use_Gillespie_MFPT == False  :
			bounds = None
			return self.find_meanfirstpassagetime_truncatedCTMC(  ), bounds
		elif build_truncatedCTMC == False and use_Gillespie_MFPT == True :
			return self.find_meanfirstpassagetime_Gillespie()
		else:
			raise  ValueError( 'In config_file.txt,  set build_truncatedCTMC=1 and use_Gillespie_MFPT=0 or vice versa' )
	
	def print_trajectory(self,o):
		print o.full_trajectory[0][0][3]   # the strand sequence  #if you get  list index out of range,  set this options.output_interval = 1
		print o.start_state[0].structure   # the starting structure
		for i in range(len(o.full_trajectory)):
			time = o.full_trajectory_times[i]
			state = o.full_trajectory[i][0]
			struct = state[4]
			sequence = state[3]
			dG = state[5]
			print struct + ' t=%11.9f seconds, dG=%6.2f kcal/mol' % (time, dG)
	
	"""Estimating the Mean First Passage Time with Gillespie SSA"""
	def find_meanfirstpassagetime_Gillespie(self) :
		if use_Gillespie_Multithread == True :
			myMS = MergeSim()
			myMS.setNumOfThreads(1)
			myMS.setOptionsFactory2( doReaction2,self.num_simulations_Gillespie , [self.num_simulations_Gillespie,self.pathwayelaboration_kappa, self.reaction_type, self.dataset_type,    self.strands_list  , self.sodium, self.magnesium, self.kinetic_parameters_simulation, self.bimolecular_reaction, self.temperature , self.temperature_change,  self.join_concentration_change , self.join_concentration, self.rate_method ,self.startStates, self.simulation_mode])
			myMS.setPassageMode( )
			myMS.setTerminationCriteria(terminationCount=self.num_simulations_Gillespie)
			myMS.run()
			myRates = myMS.results
			del myMS
		else:
			options = doReaction([self.num_simulations_Gillespie,self.simulation_time_Gillespie, self.reaction_type, self.dataset_type,    self.strands_list  , self.sodium, self.magnesium, self.kinetic_parameters_simulation, self.bimolecular_reaction, self.temperature , self.temperature_change,  self.join_concentration_change , self.join_concentration, self.rate_method ,  self.startStates , self.simulation_mode])
			options.output_interval = 1 #to store all the transitions!!!!!!!!!!!!!!!!!!
			s = SimSystem(options )
			s.start()
			if self.simulation_mode ==  "first passage" or self.simulation_mode =="trajectory" :
				myRates = FirstPassageRate( options.interface.results)
			else:
				raise ValueError('Should use other simulation mode?! ' )
			del s
		#print self.print_trajectory(options)
		if self.bimolecular_reaction  == True :
			k  = myRates.kEff(self.join_concentration)
			bootstrap = Bootstrap(myRates, concentration=options.join_concentration, N=1000, computek1=False )
		#print("Estimated 95% confidence interval: [","{:.2e}".format(bounds[0]),",","{:.2e}".format(bounds[1]),"] ", sep="")
		else :
			k = myRates.k1()
			bootstrap = Bootstrap(myRates, concentration=options.join_concentration, N=1000, computek1=True)
		
		#print("Estimated 95% confidence interval: [","{:.2e}".format(bounds[0]),",","{:.2e}".format(bounds[1]),"] ", sep="")
		klog= math.log10 (k )
		bounds = bootstrap.ninetyFivePercentilesLog()
		#print self.bimolecular_reaction, klog,  bounds
		del myRates
		return  klog , bounds
	
	
	def call_Builder(self) :
		if build_truncatedCTMC_pathwayelaboration == True and build_truncatedCTMC_GillespieSSA == False   :
			"""Building Truncated CTMC with Pathway Elaboration"""
			myBuilder = Builder(doReaction, [ self.pathwayelaboration_K, self.pathwayelaboration_kappa, self.reaction_type, self.dataset_type,    self.strands_list  , self.sodium, self.magnesium, self.kinetic_parameters_simulation, self.bimolecular_reaction,self.temperature , self.temperature_change,  self.join_concentration_change , self.join_concentration, self.rate_method ,  self.startStates ,self.simulation_mode])
		elif  build_truncatedCTMC_pathwayelaboration == False and build_truncatedCTMC_GillespieSSA == True  :
			"""Building Truncated CTMC with Gillespie"""
			myBuilder = Builder(doReaction, [ self.num_simulations_Gillespie, self.simulation_time_Gillespie, self.reaction_type, self.dataset_type,    self.strands_list  , self.sodium, self.magnesium, self.kinetic_parameters_simulation, self.bimolecular_reaction,self.temperature , self.temperature_change,  self.join_concentration_change , self.join_concentration, self.rate_method ,self.startStates ,self.simulation_mode])
		else:
			raise  ValueError( 'In config_file.txt, set build_truncatedCTMC_pathwayelaboration=1 and build_truncatedCTMC_GillespieSSA=0 or vice versa' )
		if self.pathwayelaboration_beta <0 or self.pathwayelaboration_beta > 1 :
			raise ValueError('pathwayelaboration_beta should be between 0 and 1 ')
		
		myBuilder.outfile = self.builderpath
		attributes_file  = self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+ "/atrributes" +str(self.docID)+".txt"
		attributes_file = open(attributes_file, 'a')
		
		start_time = timeit.default_timer()
		
	
		if self.pathwayelaboration_K > 0 :
		
			build_starttime = timeit.default_timer()
			fatten_regardless = True
			if build_truncatedCTMC_GillespieSSA == False   :
				reactionsargs = [ self.pathwayelaboration_K, self.pathwayelaboration_kappa, self.reaction_type, self.dataset_type,    self.strands_list  , self.sodium, self.magnesium, self.kinetic_parameters_simulation, self.bimolecular_reaction,self.temperature , self.temperature_change,  self.join_concentration_change , self.join_concentration, self.rate_method ,  self.startStates , self.simulation_mode]
			elif build_truncatedCTMC_GillespieSSA ==True :
				reactionsargs = [ self.num_simulations_Gillespie, self.simulation_time_Gillespie, self.reaction_type, self.dataset_type,    self.strands_list  , self.sodium, self.magnesium, self.kinetic_parameters_simulation, self.bimolecular_reaction,self.temperature , self.temperature_change,  self.join_concentration_change , self.join_concentration, self.rate_method ,  self.startStates , self.simulation_mode]
			
			ID=   self.docID
			reactionsargspath= self.builderpath + "reactionargs"+ ID
			startstatespath= self.builderpath + "startstatespath"+ ID
			uniqueIDtoLength_path = self.builderpath + "uniqueIDtoLength" + ID
			auto_strand_id_list_path = self.builderpath + "auto_strand_id_list" + ID
			boltzmann_path = self.builderpath + "boltzmann"   +  ID
			with open (uniqueIDtoLength_path, "wb") as  p :
				pickle.dump (self.uniqueIDtoLength, p)
			
			with open (auto_strand_id_list_path, "wb") as  p :
				pickle.dump (self.auto_strand_id_list, p)
			
			with open (reactionsargspath, "wb") as  p :
				pickle.dump (reactionsargs, p)
			
			with open (startstatespath, "wb") as  p :
				pickle.dump (self.startStates, p)
			
			with open(boltzmann_path ,"wb") as p:
				pickle.dump( self.boltzmann_sample_initialstates_list, p )
			printMeanTime = True
			builderpicklepath = self.builderpath + "builderpickle" + ID
			breaklooppath= self.builderpath + "breakloop" + ID
			len_protoSpacepath = self.builderpath + "len_protoSpace" + ID
			
			convergence_crit= "-10000" #not used
			command  = ["python", "helper_doreaction.py"  ,  self.builderpath, str(convergence_crit ), str(int(printMeanTime)), str( int(fatten_regardless) ),ID, reactionsargspath, startstatespath, builderpicklepath, breaklooppath, len_protoSpacepath  , str(self.pathwayelaboration_N)  , str(self.pathwayelaboration_beta)  , uniqueIDtoLength_path , auto_strand_id_list_path, boltzmann_path, str(int(self.pathwayelaboration_use_elaboration))]
			shell = call( command )
			del shell
			self.build_time= timeit.default_timer()- build_starttime
			with open(builderpicklepath, "rb") as p :
				myBuilder= pickle.load(p)
			self.multistrand_simulation_time= myBuilder.multistrand_simulation_time
			self.loading_prototransitions_time  = myBuilder.loading_prototransitions_time
			self.loading_protoinitialstates_time = myBuilder.loading_protoinitialstates_time
			self.loading_protofinalstates_time  = myBuilder.loading_protofinalstates_time
			self.loading_protospace_time  = myBuilder.loading_protospace_time
			with open(len_protoSpacepath, "rb") as p :
				len_protoSpace= pickle.load(p)
			os.remove(uniqueIDtoLength_path)
			os.remove(auto_strand_id_list_path)
			os.remove(builderpicklepath)
			os.remove(breaklooppath)
			os.remove(len_protoSpacepath)
			
			os.remove(boltzmann_path)
			os.remove(reactionsargspath)
			os.remove(startstatespath)
			start  = 0
			##removedfattening comment below
			myBuilder.protoSpacebackup= copy.deepcopy(myBuilder.protoSpace)
			myBuilder.protoinitialpathwayuniqueIDsbackup= copy.deepcopy(myBuilder.protoinitialpathwayuniqueIDs)
		
			pathprotoinitialpathwayuniqueIDs=   self.builderpath + "protoinitialpathwayuniqueIDsbackup"
			with open(pathprotoinitialpathwayuniqueIDs   , "wb" ) as p :
				pickle.dump(myBuilder.protoinitialpathwayuniqueIDsbackup,  p )
			pathspace=   self.builderpath + "protoSpacebackup"
			with open(pathspace   , "wb" ) as p :
				pickle.dump(myBuilder.protoSpacebackup,  p )
			pathsequences= self.builderpath + "protoSequences"
			with open(pathsequences , "wb" ) as p :
				pickle.dump(myBuilder.protoSequences,  p )
			pathoptions = self.builderpath + "pathoptions"
			with open(pathoptions , "wb" ) as p :
				pickle.dump(myBuilder.optionsArgs,  p )
			fatten_starttime = timeit.default_timer()
			if self.pathwayelaboration_use_elaboration  == True:
				onlyfatteninitialpathway = 0
			elif self.pathwayelaboration_use_elaboration == False :
				onlyfatteninitialpathway  = 1  # because  they are not connected so you need to fatten the initial pathway to make them connected
			batchsize = 1000
			lenp = len( myBuilder.protoSpace)
			print "Adding missing transitions in the truncated CTMC"
			while start < lenp :
				end = min(lenp ,  start+batchsize)
				
				command = ["python", "helper_fatten.py"  , str(start), str(end) , self.builderpath, pathspace, pathsequences, pathoptions , pathprotoinitialpathwayuniqueIDs, str(onlyfatteninitialpathway) ]
				shell = call(command )
				del shell
				print "progress --------------------  " ,  str(end)  , " / ", str(lenp)
				st  = timeit.default_timer()
				
				with open( self.builderpath + "pt" + str(start)+"-"  +str(end) , "rb" ) as p:
					protoTransitions  = pickle.load (  p)
				os.remove(self.builderpath + "pt" + str(start)+"-"  +str(end))
				myBuilder.transitionMerge_2(protoTransitions)
				with open( self.builderpath + "pt" + str(start)+"-"  +str(end) + "-fattentimes", "rb" ) as myfile:
					for line in myfile:
						
						line = line.split()
						self.multistrand_simulation_fatten_time += float( line[0  ]  )
						self.loading_protospace_fatten_time  += float( line[1 ]  )
						self.loading_prototransitions_fatten_time +=float( line[2]  )
						self.loading_protoinitialstates_fatten_time += float( line[3 ]  )
						self.loading_protofinalstates_fatten_time  +=  float( line[4  ]  )
				os.remove(self.builderpath + "pt" + str(start)+"-"  +str(end) + "-fattentimes")
				start  =  end
			os.remove(pathprotoinitialpathwayuniqueIDs)
			os.remove(pathspace)
			os.remove(pathsequences)
			os.remove(pathoptions)
			del myBuilder.protoSpacebackup
			##removedfattening comment above
			try :
				self.fatten_time= timeit.default_timer()- fatten_starttime
			except :
				self.fatten_time = 0
			print "Finished adding missing transitions in the truncated CTMC"
		else:
			raise ValueError("pathwayelaboration_K should be set to greater than 0, even if not elaborating")
		
		attributes_file.close()
		self.multistrand_time =timeit.default_timer() - start_time
		
		return myBuilder
	
	def doDeltaPruning(self, myBuilder ):
	
		print "Going in doDeltaPruning"
		deltapruningstatespacesize_file_name =  self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+ "/deltapruningstatespacesize_file_name"  +str(self.docID)+".txt" +"numberofnucleotides"
		deltapruningstatespacesize_file = open(deltapruningstatespacesize_file_name, 'wb')
		deltapruningmatrixtime_file_name = self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+ "/deltapruningmatrixtime_file_name"  +str(self.docID)+".txt" +"numberofnucleotides"
		deltapruningmatrixtime_file = open(deltapruningmatrixtime_file_name , 'wb')
		deltapruninglog10k_file_name = self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+ "/deltapruninglog10k_file_name"  +str(self.docID)+".txt" +"numberofnucleotides"
		deltapruninglog10k_file = open(deltapruninglog10k_file_name, 'wb')
		
		deltapruningstatespacesize_file.close()
		deltapruningmatrixtime_file.close()
		deltapruninglog10k_file.close()
		
		#for deltaValue in [ 0, 10 ** (-7), 10 ** (-6) , 10 ** (-5), 10 ** (-4), 10 ** (-3), 10 ** (-2),   0.05, 10 ** (-1) ,  0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 ]:
		if self.do_inference == False:
			valuelist= [ 0,  0.01, 0.05 ,  0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8  ]
		elif self.do_inference == True:
			#valuelist= [ 0.01, 0.05 ,  0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 ]  in bara 200 iter estefade shod vali kheili kond bood
			valuelist= [  0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 ]  #in bara 100 estefade mishe iter
		
		times_donotchange = None
		mfpt_donotchange= None
		for deltaValue in valuelist:
		
			deltapruningstatespacesize_file_name = self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+ "/deltapruningstatespacesize_file_name"  +str(self.docID)+".txt" +"numberofnucleotides"
			deltapruningstatespacesize_file = open(deltapruningstatespacesize_file_name, 'a')
			deltapruningmatrixtime_file_name = self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+ "/deltapruningmatrixtime_file_name"  +str(self.docID)+".txt" +"numberofnucleotides"
			deltapruningmatrixtime_file = open(deltapruningmatrixtime_file_name , 'a')
			deltapruninglog10k_file_name = self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+ "/deltapruninglog10k_file_name"  +str(self.docID)+".txt" +"numberofnucleotides"
			deltapruninglog10k_file = open(deltapruninglog10k_file_name, 'a')
			print "Going  to delta prune and solve matrix with value ", deltaValue
			myBuilder.pathwayelaboration_use_elaboration= self.pathwayelaboration_use_elaboration
			newbuilder = copy.deepcopy(myBuilder)
			newbuilder.deltaPruningValue = deltaValue
			
			times_donotchange ,mfpt_donotchange = 	newbuilder.deltaPruning(newbuilder.deltaPruningValue, times_donotchange ,mfpt_donotchange , printCount=True, maxiter= self.maxiter)
			#newbuilder.deltaPruning(newbuilder.deltaPruningValue, printCount=True)
			sizeofstatespace =  len(newbuilder.protoSpace) - len(newbuilder.protoFinalStates)
			newbuilderRate = BuilderRate(newbuilder)
			solvestart= timeit.default_timer( )
			times2 = newbuilderRate.averageTime(maxiter = self.maxiter)
			solvefinish = timeit.default_timer()
			matrixtime  = solvefinish - solvestart
			
			solutions_builderRate2 = newbuilderRate.weightedPassageTime(times2 ,bimolecular=  False )
			meanfirstpassagetime = solutions_builderRate2
			if self.bimolecular_reaction == True :
				predicted_rate= 1.0 / (meanfirstpassagetime * self.join_concentration)
			else :
				predicted_rate= 1.0  / meanfirstpassagetime
			
			predicted_log_10_rate =np.log10(predicted_rate)
			
			deltapruningstatespacesize_file.write(str(sizeofstatespace) + ",")
			deltapruningmatrixtime_file.write(str(matrixtime) + ",")
			deltapruninglog10k_file.write(str(predicted_log_10_rate) +",")
			deltapruningstatespacesize_file.close()
			deltapruningmatrixtime_file.close()
			deltapruninglog10k_file.close()
			if self.do_inference == True:
				if matrixtime <self.maxsolvetime:
					return  newbuilder
			if self.do_inference == True and deltaValue != valuelist[-1]:
				del newbuilder
				del newbuilderRate
		if self.do_inference == True:
			return newbuilder
	
	"""Estimating the Mean First Passage Time with a Truncated CTMC"""
	def find_meanfirstpassagetime_truncatedCTMC(self  ) :
		attributes_file  = self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+ "/atrributes" +str(self.docID)+".txt"
		attributes_file = open(attributes_file, 'a')
		if self.load == True  :
			self.builderpath = self.builderpath
			if os.path.exists(self.builderpath ):
				start_time = timeit.default_timer()
				attributes_file.write( "Starting to load Builder   ")
				with open( self.builderpath   , "rb" ) as builderpickle:
					myBuilder  = pickle.load( builderpickle)
				attributes_file.write("\n"+"finished loading Builder took :  " +str(timeit.default_timer()  -start_time) + "\n")
				
				
				if self.rate_method == Literals.metropolis:
					set_Metropolis_params(myBuilder.options, self.kinetic_parameters_simulation)
				elif  self.rate_method == Literals.arrhenius:
					set_Arrhenius_params(myBuilder.options,self.kinetic_parameters_simulation)
				else:
					raise ValueError(' Parameter method not set correctly ' )
		
		
		successful_simulation = False
		while successful_simulation == False :
			#print "in while loop until successful_simulation  == True ",  self.docID
			try:
				#print "self.load" , self.load
				if self.load == False :
					print "Going to build truncated CTMC from scratch"
					myBuilder= self.call_Builder( )
					print "Finished building truncated CTMC from scratch"
				
				#Resetting concentration and temperature to correct values , becauese might have used different values to accelerate simulations
				myBuilder.options.join_concentration = self.join_concentration
				myBuilder.options.temperature= self.temperature
				if self.save == True :
					attributes_file.write("Started saving builder")
					start_time = timeit.default_timer()
					builderpickle =  open( self.builderpath  , "wb" )
					pickle.dump(myBuilder,  builderpickle, -1)
					builderpickle.close()
					attributes_file.write("Finished saving builder" +str(timeit.default_timer() - start_time)+ "\n")
				start_time = timeit.default_timer()
				attributes_file.write("Starting to make builderRate from  Builder")
				print "Starting to make builderrate from builder"
				solutions_builderRate = None				
				builderRate = BuilderRate(myBuilder)
				print "Finished making builderrate from builder"
			
				attributes_file.write( "finished making builderRate from Builder :  "  +str(timeit.default_timer()  - start_time )+ "\n")
				builderRate.solveToggle =self.solveToggle
				
				attributes_file.write(str(myBuilder)+ "\n")
			
				start_time  = timeit.default_timer()
				attributes_file.write("started solving matrix iteratively ")
				#times= builderRate.averageTime(x0= prev_times, maxiter = self.maxiter)
				
				print "Going to solve matrix "
				times= builderRate.averageTime( maxiter = self.maxiter)
				print "Finished solving matrix.  "
				matrixtime  = timeit.default_timer()- start_time
				attributes_file.write("Finished solving matrix. "  + " Matrix solving time is: " +  str( matrixtime)+ "\n")
				start_time2  = timeit.default_timer()
				solutions_builderRate = builderRate.weightedPassageTime(times ,bimolecular=  False )
				attributes_file.write("Weighted passage time is: " +  str( timeit.default_timer()- start_time2 )+ "\n")
				attributes_file.write(str(myBuilder)+"\n")
				attributes_file.write("Solve toggle is:" + str(builderRate.solveToggle) + "     Max iteration is: " +str(self.maxiter) +  "      Solution is: " + str( solutions_builderRate ) + "\n")

				self.matrixsolve_time = builderRate.matrixTime
				self.statespacesize  = len(builderRate.build.protoSpace)
			
				#del builderRate
				#del prev_times
				successful_simulation =  True
			except  Exception as e:
				print "exception #! " , e,  e.message, e.args , self.builderpath
				
				if "file" in e  or "directory" in e :
					print "Got no such file or directory error!"
				else :
					print "State didn't exist!"
				print "Continuing after exception, starting all over again :("
				
				self.load = False
				self.save  = True
		
		attributes_file.close()
		
	
	
		if  self.deltaPruning  == True :
			if self.do_inference== False:
				self.doDeltaPruning(myBuilder)
			elif self.do_inference == True and matrixtime > self.maxsolvetime  :
				#for inference only do delta pruning if matrixtime is large
				newbuilder = self.doDeltaPruning(myBuilder)
				#if self.save  == True:
				builderpickle =  open( self.builderpath  , "wb" )
				pickle.dump(newbuilder,  builderpickle, -1)
				builderpickle.close()
				
		return solutions_builderRate
	
	""" Calculates Mean First Passage Time, from a truncated CTMC or Gillespie SSA, and then returns the log reaction rate constant"""
	def find_answers(self):
		#return   Output( error = 4,  predicted_log_10_rate = 4 , real_log_10_rate= 4 )
		start_time = timeit.default_timer()
		meanfirstpassagetime , bounds = self.find_meanfirstpassagetime()
		
		attributes_file = open(self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+ "/atrributes" +str(self.docID)+".txt", 'a')
		
		if build_truncatedCTMC == True and use_Gillespie_MFPT  == False :
			if meanfirstpassagetime == None or  meanfirstpassagetime  <= 0:
				#meanfirstpassagetime should be greater then 0!
				attributes_file.write( "None or Negative MFPT is: "+self.docID  +str( meanfirstpassagetime  )+ "\n")
				return   Output( error = np.inf ,  predicted_log_10_rate = np.inf , real_log_10_rate= np.inf )
			if self.bimolecular_reaction == True :
				predicted_rate= 1.0 / (meanfirstpassagetime * self.join_concentration)
			else :
				predicted_rate= 1.0 / meanfirstpassagetime
		
		elif build_truncatedCTMC == False and use_Gillespie_MFPT  == True:
			predicted_log_10_rate= meanfirstpassagetime
		else:
			raise  ValueError( 'In config_file.txt, set build_truncatedCTMC=1 and use_Gillespie_MFPT=0 or vice versa' )
		
		warnings.filterwarnings('error')
		try :
			if build_truncatedCTMC  == True :
				predicted_log_10_rate =np.log10(predicted_rate)
			real_log_10_rate = np.log10(self.real_rate)
			error  = math.pow( real_log_10_rate - predicted_log_10_rate, 2)
		except  Exception as e :
			print (  " Exception occurred .. " ), e,  e.args
			return   Output( error = np.inf ,  predicted_log_10_rate = np.inf , real_log_10_rate= np.inf )
		
		lenstrands =  0
		for strand  in self.strands_list :
			lenstrands += len(strand)
		uniqueidtowrite= self.dataset_path+"/"+myenums.Permanent_Folder.MYBUILDER.value+ "/mystats"  +str(self.docID)
		total_time = timeit.default_timer() - start_time
		
		#Saveing some statistics to file
		with open(uniqueidtowrite+".txt", 'a') as p:
			p.write( str( predicted_log_10_rate )+" "  + str(real_log_10_rate) + " "  +str( total_time) + " "  + str(self.pathwayelaboration_K)+ " "  + str(self.num_simulations_Gillespie) + " "  +str(self.temperature)+" "  + str(self.temperature_change) +" "  + str(self.join_concentration )  + " "  + str(self.join_concentration_change) + " "  +  str(self.statespacesize)  + " " +  str(self.matrixsolve_time +  self.build_time  + self.fatten_time ) + " "+  str(self.matrixsolve_time)  +  " "  + str(self.build_time) +  " "+  str( self.fatten_time) +  " "+  str(0)+  " "+  str(0) + " " + str(0) + " " +str(0)  +"\n" )
		with open(uniqueidtowrite+"_numberofnucleotides.txt", 'w')  as p :
			p.write( str(lenstrands))
		with open(uniqueidtowrite+"_matrixsolvetime.txt", 'w') as p:
			p.write( str(self.matrixsolve_time))
		with open(uniqueidtowrite+"_onlypathway.txt", 'a') as p:
			p.write( str(    self.multistrand_simulation_time)   +    " " + 		 str(self.loading_protospace_time)  + " " + str(self.loading_prototransitions_time)  +    " " +  		 str(self.loading_protoinitialstates_time)  +    " " + 		 str(self.loading_protofinalstates_time) +  " " + str(self.multistrand_simulation_fatten_time)   +    " " + 		 str(self.loading_protospace_fatten_time)  +    " " + str(self.loading_prototransitions_fatten_time)  +    " " +  		 str(self.loading_protoinitialstates_fatten_time)  +    " " +  str(self.loading_protofinalstates_fatten_time)  +   " " + str(0) + " " + str(0)+ " " + str(0)+ " " + str(0)+ " " + str(0) + " " +str(0) +"\n" )
		with open(uniqueidtowrite+"_statespacesize.txt", 'w') as p:
			p.write( str(self.statespacesize) )
	
		#print "----------------------------------------------------- " , self.multistrand_simulation_time, self.loading_protospace_time,self.loading_prototransitions_time,  self.loading_protoinitialstates_time,  self.loading_protofinalstates_time
		#tb =  self.multistrand_simulation_time + self.loading_protospace_time+  self.loading_prototransitions_time +   self.loading_protoinitialstates_time + self.loading_protofinalstates_time
		#print "----------------------------------------------------- " , self.build_time, tb , "             diff: ", self.build_time - tb
		#print "+++++++++++++++++++++++++++++++++++++++++++++++++++++ " , self.multistrand_simulation_fatten_time, self.loading_protospace_fatten_time,self.loading_prototransitions_fatten_time,  self.loading_protoinitialstates_fatten_time,  self.loading_protofinalstates_fatten_time
		#tb =  self.multistrand_simulation_fatten_time +self.loading_protospace_fatten_time + self.loading_prototransitions_fatten_time +   self.loading_protoinitialstates_fatten_time + self.loading_protofinalstates_fatten_time
		#print "+++++++++++++++++++++++++++++++++++++++++++++++++++++ " , self.fatten_time, tb ,  "             diff: ", self.fatten_time - tb
		
		print  "Error: " + str(error) +  "         real_log_10_rate: "  +str (real_log_10_rate)  +  "         predicted_log_10_rate: "  +str (predicted_log_10_rate)  + "\n\n\n"
		attributes_file.write(  "Error: " + str(error) +  "         real_log_10_rate: "  +str (real_log_10_rate)  +  "         predicted_log_10_rate: "  +str (predicted_log_10_rate)  + "\n" )
		attributes_file.write("Next iteration "+"\n")
		attributes_file.close()
		
		return   Output( error = error ,  predicted_log_10_rate =  predicted_log_10_rate , real_log_10_rate= real_log_10_rate )
		
def doReaction(arguments ) :
	#options = Options()
	options = Options(trials = arguments[0])
	#options.output_interval = 1
	# Important note: 	# the first argument is always the number of paths
	options.num_simulations=  int ( arguments[0] )
	options.simulation_time =  float( arguments[1] )
	options.sodium = arguments[5]
	options.magnesium = arguments[6]
	options.temperature = arguments[9] +arguments[10]
	options.join_concentration =    arguments[12] * arguments[11]
	if arguments[13] == Literals.metropolis :
		set_Metropolis_params(options, arguments[7] )
	elif arguments[13] == Literals.arrhenius :
		set_Arrhenius_params(options, arguments[7]  )
	if arguments[15]  == "trajectory":
		options.simulation_mode = Literals.trajectory
	elif arguments [15] =="first_step":
		options.simulation_mode =  Literals.first_step
	else:
		raise ValueError(' set options.simulation_mode in parent.py  ')
	
	endComplex1 = arguments[14][-1][0] #Consider the last initial state, using it's dissociated strand as complex
	#endComplex2 = arguments[15][-1][1] #Consider the last initial state, using it's dissociated strand as complex
	#stopSuccess = StopCondition(Literals.success, [(endComplex1, Literals.exact_macrostate, 0), (endComplex2, Literals.exact_macrostate, 0)])
	stopSuccess = StopCondition(Literals.success, [(endComplex1, Literals.exact_macrostate, 0)])
	options.stop_conditions = [stopSuccess]
	if build_truncatedCTMC_GillespieSSA == True or use_Gillespie_MFPT ==True  :
		options.start_state = arguments[14][0]
	return options

def doReaction2(n_trials, arguments ) :
	arguments[0] = n_trials
	return doReaction(arguments)

"""Setting parameters for the Arrhenius kinetic model"""
def set_Arrhenius_params (options, params):
	LNA_INDEX= 0
	E_INDEX = 1
	options.rate_method = Literals.arrhenius
	options.lnAStack = float( params[localtype.stack][LNA_INDEX] )
	options.EStack = float(  params[localtype.stack][E_INDEX] )
	options.lnALoop = float( params[localtype.loop][LNA_INDEX] )
	options.ELoop = float( params[localtype.loop][E_INDEX] )
	options.lnAEnd = float( params[localtype.end][LNA_INDEX] )
	options.EEnd = float(params[localtype.end][E_INDEX] )
	options.lnAStackLoop =  float(params[localtype.stackloop][LNA_INDEX])
	options.EStackLoop =   float( params[localtype.stackloop][E_INDEX] )
	options.lnAStackEnd = float(params[localtype.stackend][LNA_INDEX] )
	options.EStackEnd = float(params[localtype.stackend][E_INDEX] )
	options.lnALoopEnd = float(params[localtype.loopend][LNA_INDEX] )
	options.ELoopEnd =float( params[localtype.loopend][E_INDEX] )
	options.lnAStackStack = float(params[localtype.stackstack][LNA_INDEX] )
	options.EStackStack =float( params[localtype.stackstack][E_INDEX]  )
	options.bimolecular_scaling = float(params[bimolecular_scaling] )

"""Setting parameters for the Metropolis kinetic model"""
def set_Metropolis_params(options, params):
	options.rate_method = Literals.metropolis
	options.unimolecular_scaling =  float(params [unimolecular_scaling])
	options.bimolecular_scaling  =float( params[bimolecular_scaling])

def main(complex ):
	return complex.find_answers( )
