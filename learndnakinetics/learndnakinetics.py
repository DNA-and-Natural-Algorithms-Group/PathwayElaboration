from __future__ import division
import ConfigParser
import numpy as np
import csv
import timeit
import os
import multiprocessing
import sys
import math
sys.path.insert(0,os.path.realpath('../reactions'))
import parent
from parent import *
from datasets import *
import myenums
import string
import copy

DATASET_PATH = '../dataset'

iter = 0

def Complement(str):
	NUCLEOTIDES = "ACTG"
	TRANSLATION_TABLE = string.maketrans(NUCLEOTIDES, "TGAC")
	return ''.join(list(reversed(str))).translate(TRANSLATION_TABLE)

"""creating required directories"""
def initconf(my_name , directories ) :
	dataset_path =  parameter_folderround +my_name
	document = DATASET_PATH + my_name + '.csv'
	directories +=  [   dataset_path + "/" +  myenums.Permanent_Folder.MYBUILDER.value   ]
	if not os.path.exists(dataset_path):
		os.makedirs(dataset_path)
	row =  open_document(document)
	
	return dataset_path, document , row

"""A reaction object, general conditions are set here, reaction specific conditions are set in function read_dataset in this file and in datasets.py """
class Reaction(object):
	
	def __init__(self, **kwargs):
		#All attributes could be overrided in datasets.py
		self.counter_celllist= None # do not change this line
		self.pathwayelaboration_N= pathwayelaboration_N
		self.pathwayelaboration_beta=  pathwayelaboration_beta
		self.pathwayelaboration_K = pathwayelaboration_K
		self.pathwayelaboration_kappa =   pathwayelaboration_kappa
		self.pathwayelaboration_use_elaboration = pathwayelaboration_use_elaboration
		self.theta_simulation  = theta_simulation
		self.rate_method= rate_method
		self.load  = load
		self.save = save
		self.maxiter = 1000000
		self.maxsolvetime =  120 #will delta-prune if using more than maxsolvetime
		self.solveToggle =2000 #will default to using exact solver
		self.simulation_mode= "trajectory"
		self.do_inference= do_inference
		self.deltaPruning = deltaPruning
		self.temperature_change=  0
		self.concentration_change =1
		
		#These two parameters are used if you build the truncated CTMC with Gillespie SSA.
		self.num_simulations_Gillespie = num_simulations_Gillespie
		self.simulation_time_Gillespie=  simulation_time_Gillespie
	
		for k, v in kwargs.items():
			setattr(self, k, v)
	
	def set_specific ( self,  **kwargs) :
		for k, v in kwargs.items():
			setattr(self, k, v)
	
	def __str__(self):
		output =""
		output  += " num_simulations: " + str(self.num_simulations)
		output  += " simulation_time: " + str(self.simulation_time)
		output  += " temperature: " + str(self.temperature) + " 1000/T:  " + str(1000/self.temperature)
		output  += " concentration: " + str(self.concentration)
		output  += " temperature_change: " + str(self.temperature_change)
		output  += " concentration_change: " + str(self.concentration_change)
		output  += " sodium: " + str(self.sodium)
		output  += " magnesium " + str(self.magnesium)
		output  += " real_rate  " + str(self. real_rate ) + " log10 k:  "  + str( self.real_rate )
		return output

"""reading reactions from database and setting their conditions """
def read_dataset( done_queue ,  dataset ) :
	temperature_column= 1
	toeholdlength_column = [16, 17 ]
	rate_column =[ 4, 18]
	flur_position  =  ""
	toehold_length = ""
	sodium_column = 2
	magnesium_column = 3
	timescale_column = 5
	concentration_column = [6,7]
	strand_column = [8,9, 10, 11, 12, 13, 14, 15]
	to_column =[20, 21 ]
	mismatchposition_column = [22]   #[22,23]
	
	try :
		temperature = float(dataset.row[dataset.counter_cell][temperature_column])
	except :
		temperature = 25
	if dataset.temperature_type ==myenums.DatasetSpecifications.TEMPKELVININV.value  :
		temperature = 1000/ temperature - 273.15
	elif dataset.temperature_type ==myenums.DatasetSpecifications.TEMPKELVIN.value  :
		temperature =  temperature- 273.15
		
	try :
		sodium = float(dataset.row[dataset.counter_cell][sodium_column])
	except:
		sodium = 0.01
	if sodium < 0.01 :
		sodium = 0.01
	try :
		magnesium = float(dataset.row[dataset.counter_cell][magnesium_column])
	except:
		magnesium  = 0
	
	try :
		concentrationsize = 2
		concentration = np.max ( (float(dataset.row[dataset.counter_cell][concentration_column[0]]) , float(dataset.row[dataset.counter_cell][concentration_column[1]]) ))
	except :
		concentrationsize = 1
	
	if concentrationsize == 1:
		try :
			concentration =float(dataset.row[dataset.counter_cell][concentration_column[0]])
		except:
			concentration = 0.0000000001
	try:
		real_rate = float(dataset.row [dataset.counter_cell][rate_column[0]])
		if dataset.rate_type ==myenums.DatasetSpecifications.LOG10RATECONSTANT.value :
			
			real_rate =  math.pow(10, real_rate)
		if dataset.rate_type == myenums.DatasetSpecifications.RATECONSTANT10POW5.value:
			real_rate = real_rate *  ( 10 **5 )
	except :
		timescale = float (dataset.row [dataset.counter_cell ][timescale_column])
		if dataset.rate_type == myenums.DatasetSpecifications.TIMESCALE.value  :
			if dataset.bimolecular_reaction == True :
				real_rate  = 1 /  ( timescale  * concentration)
			elif dataset.bimolecular_reaction == False  :
				
				real_rate  = 1 / timescale
	
	strands_list = []
	for i in strand_column:
		try :
			if dataset.row[dataset.counter_cell][i] != '' :
				strands_list.append(dataset.row[dataset.counter_cell][i].rstrip() )
		except   :
			pass
			#print "no more strings to read"
	
	if dataset.dataset_name == myenums.DatasetName.RAUZAN2013.value:
		strands_list = [strands_list[0], strands_list[1][::-1] ]
	
	if dataset.dataset_name == myenums.DatasetName.DUPIUS2013.value:
		strand1= strands_list[0 ][len(strands_list[0] )- len(strands_list[2]):]
		strand2= strands_list[2]
		botleftdangle= strands_list[3]
		topleftdangle = strands_list[0][0: len(strands_list[0] )- len(strands_list[2])]
		strands_list = [strand1,strand2 , botleftdangle, topleftdangle ]
	
	if dataset.dataset_name == myenums.DatasetName.REYNALDODISSOCIATE.value  or dataset.dataset_name == myenums.DatasetName.REYNALDOSEQUENTIAL.value :
		strand1 = strands_list[1]
		strand2 = Complement(strand1)
		dangle = strands_list[2]
		botleftdangle = "GAA"
		botrightdangle = dangle[len(botleftdangle) + len(strand1):   ]
		strands_list= [strand1 ,strand2,  botleftdangle, botrightdangle ]
		toehold_length =  0
	if dataset.dataset_name == myenums.DatasetName.SUYAMA.value or dataset.dataset_name == myenums.DatasetName.ZHANG_hybridization.value :
		strand1 = strands_list[0 ]
		strand2 = Complement(strand1)
		strands_list= [strand1 ,strand2]
	
	if dataset.dataset_name ==myenums.DatasetName.ZHANG.value:
		toehold_length =  int( dataset.row[dataset.counter_cell][toeholdlength_column [0]]   )
		toehold= strands_list[2][:toehold_length]
		displacement  = strands_list[1]
		topdangle= strands_list[0]
		strands_list= [toehold, displacement,  topdangle]
	
	if dataset.dataset_name ==myenums.DatasetName.ZHANG_stranddisplacement.value:
		toehold_length =  len (strands_list[1] )  -  len(strands_list[0])
		attacker = strands_list[1] #target
		substrate   = Complement( strands_list[1])
		incumbent = attacker[  : len(attacker) - toehold_length  ]
		strands_list= [attacker, substrate, incumbent]
	
	if dataset.dataset_name ==myenums.DatasetName.BROADWATER2016.value:
		incumbent = strands_list[0][::-1]
		attacker=strands_list[1][::-1]
		substrate = strands_list[2][::-1]
		strands_list= [attacker, substrate, incumbent]
	
	if dataset.dataset_name ==myenums.DatasetName.SHERRYCHEN2016.value:
		incumbent = strands_list[0]
		attacker=strands_list[1]
		substrate = strands_list[2][::-1]
		strands_list= [attacker, substrate, incumbent]
	
	if dataset.dataset_name ==myenums.DatasetName.MACHINEK.value:
		toehold_length =  int( dataset.row[dataset.counter_cell][toeholdlength_column [0]]   )
		incumbentDangle = strands_list[0]
		incumbent = strands_list[1]
		targetMachinek = strands_list[2]
		attacker = strands_list[3]
		aa =  len (targetMachinek) - len( attacker)
		substrate =  targetMachinek [aa : ]
		substrateDangle = targetMachinek[0: aa]
		strands_list = [incumbentDangle, incumbent, substrateDangle,  substrate, attacker]
	
	if dataset.dataset_name == myenums.DatasetName.DABBY.value:
		
		bimolecular_rate = real_rate
		unimolecular_rate =  float(dataset.row [dataset.counter_cell][rate_column[1]])
		t1 = 1 / (bimolecular_rate * concentration)
		t2 = 1 /unimolecular_rate
		real_rate = 1 / ((t1 + t2) * concentration)
		Xstr = strands_list[0]
		toeholdstrm = strands_list[1]
		toeholdstrn = strands_list[2]
		lenm =  int( dataset.row[dataset.counter_cell][toeholdlength_column [0]]   )
		lenn =  int( dataset.row[dataset.counter_cell][toeholdlength_column [1]]   )
		ms = int  ( dataset.row[dataset.counter_cell][to_column[0]] )
		ns = int (dataset.row[dataset.counter_cell][to_column[1]] )
		complex1  = toeholdstrm[len(toeholdstrm)- lenm: ] + Xstr
		FullToehold = toeholdstrm[len(toeholdstrm)- ms : ] + Xstr
		reporter1  = Complement(FullToehold)
		complex2 = Complement ( Xstr) + toeholdstrn[ :lenn ]
		FullToehold = Complement(Xstr ) + toeholdstrn[ :  ns  ]
		reporter2 = Complement(FullToehold)
		strands_list = [complex1, reporter1, complex2, reporter2 ]
		lenx = len(Xstr)
		dataset.set_specific  ( lenx= lenx, lenm= lenm, lenn = lenn , ms = ms , ns = ns  )
	
	if dataset.dataset_name == myenums.DatasetName.GROVES2015.value:
		Xstr = strands_list[0]
		toeholdstrm = strands_list[1]
		toeholdstrn = strands_list[2]
		lenm =  int( dataset.row[dataset.counter_cell][toeholdlength_column [0]]   )
		lenn =  int( dataset.row[dataset.counter_cell][toeholdlength_column [1]]   )
		ms =lenm
		ns = lenn
		complex1  = toeholdstrm[len(toeholdstrm)- lenm: ] + Xstr
		FullToehold = toeholdstrm[len(toeholdstrm)- ms : ] + Xstr
		reporter1  = Complement(FullToehold)
		complex2 = Complement ( Xstr) + toeholdstrn[ :lenn ]
		FullToehold = Complement(Xstr ) + toeholdstrn[ :  ns  ]
		reporter2 = Complement(FullToehold)
		strands_list = [complex1, reporter1, complex2, reporter2 ]
		lenx = len(Xstr)
		dataset.set_specific  ( lenx= lenx, lenm= lenm, lenn = lenn , ms = ms , ns = ns  )
	
	if dataset.dataset_name == myenums.DatasetName.ALTANBONNET.value :
		flur_position= 17 	#set flur position for Altan - bonnet !
		dataset.set_specific(flur_position =flur_position)
	
	
	if dataset.dataset_name == myenums.DatasetName.GAO2006.value:
		if dataset.counter_cell == 3 or dataset.counter_cell == 7 or dataset.counter_cell ==  11:
			
			dataset.set_specific(dataset_type = myenums.DatasetType.MISMATCH.value)
		
		else:
			dataset.set_specific(dataset_type = myenums.DatasetType.NODANGLE.value)
	
	dataset.set_specific(cutoff =1 )
	if 1<= dataset.counter_cell <= 5 and dataset.dataset_name == myenums.DatasetName.GROVES2015.value:
		dataset.set_specific(loopbases=4)
	elif 10<= dataset.counter_cell <= 12 and dataset.dataset_name == myenums.DatasetName.GROVES2015.value:
		dataset.set_specific(loopbases=8)
	elif dataset.dataset_name == myenums.DatasetName.GROVES2015.value or dataset.dataset_name == myenums.DatasetName.DABBY.value  :
		dataset.set_specific(loopbases= 6)
	
	dataset.set_specific(  temperature =  temperature , concentration =concentration, sodium = sodium, magnesium =  magnesium,  real_rate = real_rate, strands_list = strands_list, toehold_length = toehold_length, flur_position = flur_position)
	complex = ParentComplex(dataset)
	output = parent.main( complex )
	done_queue.put( ( dataset.dataset_name , output.error  , dataset.counter_cell, dataset.document,    output.predicted_log_10_rate, output.real_log_10_rate)  )

"""These values are unacceptable for the Arrhenius and Metropolis models"""
def filter_undefined_parameter_set(theta_simulation , alpha, sigma):
	if alpha  <= 0  or sigma <= 0  or (rate_method ==Literals.metropolis and  ( theta_simulation[0] <= 0 or theta_simulation[1]  <= 0 )  ) :
		return True
	
""" We bound the rates to be between 10** 4 and 10 ** 9 as follows) """
def filter_parameter_set(theta_simulation , alpha ):
	global simulationTimeAVG
	simulationTimeAVG =  0
	mybool = False
	unacceptable_high_unimolecularrate = math.pow (10, 9)
	unacceptable_low_unimolecularrate = math.pow (10, 4)
	if rate_method ==Literals.metropolis :
		#simulationTimeAVG =   (1./theta_simulation  [0] )
		if theta_simulation [0 ] > unacceptable_high_unimolecularrate   or theta_simulation [0 ] < unacceptable_low_unimolecularrate or theta_simulation [1 ] > unacceptable_high_unimolecularrate   or theta_simulation [1] < unacceptable_low_unimolecularrate :
			mybool = True
			print "The rates are unacceptable! "
	R= 0.0019872036  # kcal / K mol
	RT = R * (23 + 273.15)
	if rate_method == Literals.arrhenius :
		for birate in [1, theta_simulation[-1] ]  :
			n_local_contexts = 7
			count = 0
			for i in range(n_local_contexts):
				for j in range(i, n_local_contexts):
					lnA = theta_simulation[2* i   ] + theta_simulation [2* j  ]
					E = theta_simulation[2* i  +1 ] + theta_simulation [2* j  +1]
					rate = np.e ** (lnA - E / RT) * birate
					simulationTimeAVG +=   ( 1./rate)
					if rate > unacceptable_high_unimolecularrate   or rate < unacceptable_low_unimolecularrate  :
						print "The rates are unacceptable!!!!!" , rate
						mybool = True
					count+=1
			simulationTimeAVG = simulationTimeAVG / count
	
	if mybool == True :
		print  "Skipped a parameter set "
	if filter_smallandlarge_rates  ==  False  :
		mybool = False
	return mybool

def objective_function_auxilary( done_queue, dataset_list, n ,dataset):
	if dataset.counter_celllist == None :
		row = open_document(dataset.document)
		dataset.counter_celllist = [i for i in range( 1,  n_csv_rows(row))]
	for counter_cell in dataset.counter_celllist:
		datasetx = copy.deepcopy(dataset)
		#if datasetx.dataset_name ==   myenums.DatasetName.KIMFIG5.value and counter_cell ==5 :			#these already appear in kimtable1
		#continue
		docID = datasetx.docID + str(counter_cell )
		datasetx.set_specific (docID=  docID, counter_cell = counter_cell)
		if use_multiprocess == True:
			dataset_list.append( ForMultiProcess( read_dataset, (done_queue, datasetx)))
		else :
			read_dataset(done_queue,  datasetx)
		n +=1
	return  (dataset_list, done_queue, n )

"""This function returns the error of the predictions on the datasets. In map.py, the optimizer calls this function"""
def objective_function(theta_simulationp ):
	global iter, theta_simulation, load, save
	print  "Parameter set is " , theta_simulationp
	print "Starting iteration " ,iter
	print "\n"
	
	if do_inference == False :
		if load_existing_override == True :
			load= True
			save = False
		else:
			load = False
			save = True
	
	elif do_inference == True :
		#In the first round only going to build the trucnated CTMCs. In the next rounds load the trucnated CTCMs from round 0
		if iter ==0 :
			load = False
			save = True
	
		else:
			load = True
			save = False
		
		
	start_time = timeit.default_timer()
	theta_simulation =[]
	
	if use_regularizer == True:
		if (rate_method == Literals.arrhenius and len(theta_simulationp)!= 16 ) or (rate_method == Literals.metropolis and len(theta_simulationp)!=3):
			raise ValueError('[parameter set not initialized correctly!, probably forgot to set sigma initialization in map.py')
	else:
		if (rate_method == Literals.arrhenius and len(theta_simulationp)!= 15 ) or (rate_method == Literals.metropolis and len(theta_simulationp)!=2):
			raise ValueError('[parameter set not initialized correctly! ')
	
	for x in theta_simulationp :
		theta_simulation.append(x)
	if rate_method == Literals.arrhenius:
		theta_simulation = [theta_simulationp[0] , theta_simulationp[1] , theta_simulationp[2], theta_simulationp[3] , theta_simulationp[4] , theta_simulationp[5], theta_simulationp[6] ,theta_simulationp[7], theta_simulationp[8], theta_simulationp[9], theta_simulationp[10] , theta_simulationp[11],  theta_simulationp[12] , theta_simulationp[13], theta_simulationp[14] ]
		alpha = theta_simulation [14]
	elif rate_method == Literals.metropolis :
		theta_simulation = [theta_simulationp[0] , theta_simulationp[1]]
		alpha =1
	else:
		raise ValueError('Error: Please specify rate_method to be Arrhenius or Metropolis!')
	sigma = theta_simulationp[len(theta_simulationp)-1]
	
	if filter_undefined_parameter_set(theta_simulation , alpha, sigma) == True or filter_parameter_set(theta_simulation, alpha) == True :
		#if METHOD ==myenums.MethodName.MAPNAME.value:
		return np.inf
	
	parameter_file = open(parameter_file_name, 'a')
	parameter_file.write("Iteration " + str(iter) +" "+str(theta_simulation) +  " " + str(sigma) + '\n')
	parameter_file.close()
	itertime_file = open(itertime_file_name, 'a')
	overalltime_file = open(overalltime_file_name, 'a')
	error = 0
	n = 0
	done_queue =  multiprocessing.Manager().Queue()
	dataset_list = []
	directories =[]
	
	(dataset_list, done_queue, n ) =  read_Hata2017 (done_queue  =done_queue,dataset_list =dataset_list, n = n,  directories = directories )
	(dataset_list, done_queue, n ) =  read_Cisse2012(done_queue  =done_queue,dataset_list =dataset_list, n = n,  directories = directories )
	(dataset_list, done_queue, n ) =  read_Machinek2014(done_queue  =done_queue,dataset_list =dataset_list, n = n,  directories = directories )
	(dataset_list, done_queue, n ) =  read_Bonnet1998(done_queue  =done_queue,dataset_list =dataset_list, n = n,  directories = directories )
	# dataset not released (dataset_list, done_queue, n ) =  read_Zhang2018_hybridization(done_queue  =done_queue,dataset_list =dataset_list, n = n,  directories = directories )
	
	
	
	
	check_directories (directories)
	dataset_error = dict()
	dataset_count = dict()
	terror   = multi_process(done_queue , dataset_list   , iter , dataset_error , dataset_count   )
	error += terror

	regularizer = 0
	if rate_method == Literals.arrhenius :
		for i in range( len(theta_simulation) ) :
			regularizer+=  (theta_simulation[i] * theta_simulation [i] )
	elif rate_method == Literals.metropolis:
		for i in range(0, 2 ) :
			param = math.log (theta_simulation[i] )
			regularizer += (  param * param )
		for i in range( 2, len(theta_simulation) ):
			regularizer += (theta_simulation[i] * theta_simulation[i] )
	
	else:
		raise ValueError('Error: Specify rate_method to be Arrhenius or Metropolis!')
	LAMBDA =50
	regularizer = regularizer/ (2 * LAMBDA)
	if use_regularizer== True:
		lnprob = -(n + 1  )*np.log(sigma) - (error  /( 2 *(sigma ** 2) )) -  regularizer
	else:
		lnprob = -error
	negativelnprob = -lnprob
	elapsed = timeit.default_timer() - start_time
	itertime_file.write(str(elapsed) + ",")
	overalltime_file.write(str(timeit.default_timer() - OVERALL_STARTTIME) + ",")
	
	#writing information to file
	parameter_file = open(parameter_file_name, 'a')
	parameter_file.write( "Iteration:" + str(iter) +  "		 SE:" + str( error) +   "    Iteration time:" + str(format(elapsed, "10.4E")) +  "     Overalltime: " + str( format(timeit.default_timer() - OVERALL_STARTTIME ,"10.4E"))+ "      Size of dataset: "+str(n)+ '\n')
	parameter_file.write( "Iteration:" + str(iter) +  "		 MSE:" + str(error/n) + "   Overalltime/number_of_iterations: " + str(format(  (timeit.default_timer() - OVERALL_STARTTIME)/ (iter + 1) , "10.4E" ))+ "        Size of dataset: "+str(n)+ '\n')
	print "Iteration:" + str(iter) +  "   Mean Squared Error (MSE) of dataset: " + str(error/n) + "    Overalltime/number_of_iterations: " + str(format(  (timeit.default_timer() - OVERALL_STARTTIME)/ (iter + 1) , "10.4E" ))+ "   Size of dataset: "+str(n)+ '\n\n\n\n'
	parameter_file.write( "SE:  ")
	for cs in dataset_error:
		parameter_file.write(cs + ": " + str(dataset_error[cs])  + "(size:" + str(dataset_count[cs])+")" +",    ")
	newcounter = 0
	parameter_file.write( "\n")
	parameter_file.write( "MSE:  ")
	for cs in dataset_error:
		newcounter += dataset_count[cs]
		if dataset_count[cs] != 0 :
			parameter_file.write(cs + ": " + str(dataset_error[cs]/ dataset_count[cs])  + "(size:" + str(dataset_count[cs])+")" + ",    ")
			msedata = open(dataset_file_name +cs+".csv", "a")
			msedata.write(  str(dataset_error[cs]/ dataset_count[cs]) + ",")
			msedata.close()
	parameter_file.write( "\n\n")
	parameter_file.close()
	
	iter += 1   # Do not move this line or you'll get errors later
	with open(parameters_mse_file+"_mse.csv", "a") as p :
		p.write( str(error/n) +", ")
	if rate_method ==Literals.metropolis :
		print "Saving results for Metropolis model! "
		with open(parameters_mse_file+"_kuni.csv", "a") as p :
			p.write( str(theta_simulation[0]) +", ")
		with open(parameters_mse_file+"_kbi.csv", "a") as p :
			p.write(str(theta_simulation[1])+", ")
		if use_regularizer == True :
			with open(parameters_mse_file+"_sigma.csv", "a") as p :
				p.write(str(theta_simulation[2])+", ")
	
	if np.isnan(error)  or error == np.inf:
		negativelnprob = np.inf
		
	#if METHOD == myenums.MethodName.MAPNAME.value:
	return negativelnprob
	
"""This class used for multiprocessing"""
class ForMultiProcess(object):
	def __init__(self, function_name, arguments) :
		self.function_name = function_name
		self.arguments = arguments

"""multi processing function """
def multi_process(done_queue ,  dataset_list ,  iter , dataset_error ,  dataset_count) :
	
	global   predicted_logreactionrateconstants , experimental_logreactionrateconstants
	error = 0
	
	if use_multiprocess== True  :
		pool = multiprocessing.Pool( processes = n_processors)
		for ds in dataset_list:
			compile_error = pool.apply_async( ds.function_name ,  ds.arguments )
		print  ( "Errors: " + str(compile_error.get()) )
		pool.close( )
		pool.join ()
	
	while not done_queue.empty():
		(dataset_name, s, counter_cell, document  , predicted_log_10_rate ,real_log_10_rate ) = done_queue.get()
		
		if s== np.inf:
			print "np.inf occurred!", counter_cell, document
		error += s
		
		if dataset_name in dataset_error :
			dataset_error [dataset_name ] += s
			dataset_count[dataset_name] += 1
		else :
			dataset_error[dataset_name ] = s
			dataset_count[dataset_name] = 1
		
	return error

def check_directories (directories) :
	for dir in directories:
		if not os.path.exists(dir):
			os.makedirs(dir)

"""open a csv file"""
def open_document(document) :
	my_CSV = list(csv.reader(open(document, 'rb')))
	return my_CSV

"""return the number of rows in a csv file"""
def n_csv_rows(csv) :
	count = 0
	row = csv[count]
	while row [0] != '' :
		count += 1
		if count >= len (csv)  :
			break
		row= csv [count]
	return count

"""creating directory for round mew round of of inference"""
def set_folderfiles() :
	global parameter_folderround,  dataset_file_name,  overalltime_file_name, itertime_file_name, parameter_file_name,parameters_mse_file
	
	parameters_mse_file  = parameter_folder+ "/parameters_mse"
	if rate_method ==Literals.metropolis :
		with open(parameters_mse_file+"_mse.csv", "a") as p :
			p.write("*, ")
	parameter_folderround = parameter_folder+ "/round" +str(map_i)
	if not os.path.exists(parameter_folderround):
		os.makedirs(parameter_folderround)
	dataset_file_name= parameter_folder + "/MSE"
	itertime_file_name= parameter_folder + "/ITERTIME" + ".csv"
	overalltime_file_name= parameter_folder + "/OVERALLTIME"  + ".csv"
	parameter_file_name  = parameter_folderround + "/parametersandmse"

"""reading config_file.txt and setting filenames"""
def set_configuration():
	
	global parameter_folder ,   n_processors, use_multiprocess, filter_smallandlarge_rates, OVERALL_STARTTIME,  use_regularizer, deltaPruning ,do_inference , load_existing_override, pathwayelaboration_K , pathwayelaboration_beta, pathwayelaboration_N , pathwayelaboration_kappa, pathwayelaboration_use_elaboration,num_simulations_Gillespie,simulation_time_Gillespie,rate_method
	OVERALL_STARTTIME = timeit.default_timer()
	configParser = ConfigParser.ConfigParser()
	configParser.readfp(open(r'config_file.txt'))
	CONFIG_NAME = 'learndnakinetics'
	
	rate_method= configParser.getint(CONFIG_NAME, 'rate_method')

	n_processors =    configParser.getint(CONFIG_NAME, 'n_processors')
	use_multiprocess =  bool(configParser.getint(CONFIG_NAME, 'use_multiprocess') )
	filter_smallandlarge_rates =  bool(configParser.getint(CONFIG_NAME, 'filter_smallandlarge_rates') )
	use_regularizer  = bool(configParser.getint(CONFIG_NAME, 'use_regularizer') )
	do_inference =  bool(configParser.getint(CONFIG_NAME, 'do_inference') )
	deltaPruning =  bool(configParser.getint(CONFIG_NAME, 'deltaPruning') )
	load_existing_override =bool(configParser.getint(CONFIG_NAME, 'load_existing_override') )
	
	#These 4 parameters are used for the building truncated CTMCs with pathway elaboration
	pathwayelaboration_N= int (configParser.getint(CONFIG_NAME, 'pathwayelaboration_N'))
	pathwayelaboration_beta=   float (configParser.getfloat(CONFIG_NAME, 'pathwayelaboration_beta'))
	pathwayelaboration_K=    int (configParser.getfloat(CONFIG_NAME, 'pathwayelaboration_K'))
	pathwayelaboration_kappa=   float (configParser.getfloat(CONFIG_NAME, 'pathwayelaboration_kappa'))
	pathwayelaboration_use_elaboration =  bool(configParser.getint(CONFIG_NAME, 'pathwayelaboration_use_elaboration'))
	
	#These two parameters are used if you build the truncated CTMC with Gillespie SSA.
	num_simulations_Gillespie= int (configParser.getint(CONFIG_NAME, 'num_simulations_Gillespie'))
	simulation_time_Gillespie=  float (configParser.getfloat(CONFIG_NAME, 'simulation_time_Gillespie'))


	print "Using pathway elaboration to predict the kinetics of nucleic acid kinetics ... "
	
	if do_inference == True:
		if rate_method == 1 :
			print "\nStarting to estimate parameters for the Metropolis kinetic model"
		elif rate_method ==3 :
			print "\nStarting to estimate parameters for the Arrhenius kinetic model"
		else :
			raise ValueError('In config_file.txt, rate_method should be either 1 or 3')
	elif do_inference == False :
		if rate_method == 1 :
			print "\nUsing the Metropolis kinetic model."
		elif rate_method ==3 :
			print "\nUsing the Arrhenius kinetic model"
		else :
			raise ValueError('In config_file.txt, rate_method should be either 1 or 3 in the configuration file')
	else:
		raise ValueError('In config_file.txt, do_inference should be either 1 or 0' )
	
	CONFIG_NAME2 = "parent"
	
	build_truncatedCTMC= bool(configParser.getint(CONFIG_NAME2, 'build_truncatedCTMC'))
	build_truncatedCTMC_pathwayelaboration = bool(configParser.getint(CONFIG_NAME2,'build_truncatedCTMC_pathwayelaboration'))
	build_truncatedCTMC_GillespieSSA= bool(configParser.getint(CONFIG_NAME2, 'build_truncatedCTMC_GillespieSSA'))
	use_Gillespie_MFPT = bool(configParser.getint(CONFIG_NAME2, 'use_Gillespie_MFPT'))
	
	if build_truncatedCTMC == True and build_truncatedCTMC_pathwayelaboration  == True and  build_truncatedCTMC_GillespieSSA == False  and use_Gillespie_MFPT == False:
		parameter_folder = configParser.get(CONFIG_NAME, 'parameter_folder') + "_pathwayelaboration" + "_useelaboration:" + str(pathwayelaboration_use_elaboration) +"_N:"+str(pathwayelaboration_N) +"_beta:"+str(pathwayelaboration_beta)+"_K:"+str(pathwayelaboration_K)+"_kappa:"+str(pathwayelaboration_kappa)
	elif  build_truncatedCTMC == False and build_truncatedCTMC_pathwayelaboration  == False and  build_truncatedCTMC_GillespieSSA == False  and use_Gillespie_MFPT == True:
		parameter_folder = configParser.get(CONFIG_NAME, 'parameter_folder') +"_GillespieSSA:"+"_numberofsimulations:"+str(num_simulations_Gillespie)+ "_simulationtime:"+str(simulation_time_Gillespie)
	else:
		raise ValueError("In config_file.txt, predict the MFPT and reaction rate constant by building a truncated CTMC using pathway elaboration (set build_truncatedCTMC=1,build_truncatedCTMC_pathwayelaboration=1,build_truncatedCTMC_GillespieSSA=0,use_Gillespie_MFPT=0)!  Or predict the MFPT and reaction rate constant with the standard Gillespie SSA (set build_truncatedCTMC=0,build_truncatedCTMC_pathwayelaboration=0,build_truncatedCTMC_GillespieSSA=0,use_Gillespie_MFPT=1)! " )
	
	if not os.path.exists(parameter_folder):
		os.makedirs(parameter_folder)
