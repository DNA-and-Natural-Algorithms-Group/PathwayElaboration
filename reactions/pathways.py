from __future__ import division
import myenums
from multistrand.experiment import standardOptions, makeComplex
from nupack import  *
import string

""" This file generated manually designed truncated CTMCS BUT WE DO NOT USE THEM!!!!!!
 We  only use this file to generate initial states and final states """

NUCLEOTIDES = "ACTG"
TRANSLATION_TABLE = string.maketrans(NUCLEOTIDES, "TGAC")
def Complement(str):
	return ''.join(list(reversed(str))).translate(TRANSLATION_TABLE)



class Pathway (object) :
	def __init__(self ,  strands_list , reaction_type,  dataset_type, dataset_name   )  :
		self.dataset_type = dataset_type
		self.reaction_type = reaction_type
		self.dataset_name =  dataset_name
		self.strands_list  = strands_list
	
	def get_uniqueIDtoLengthandOrder(self) :
		uniqueIDtoLength= dict()
		
		for i in range(len(self.auto_strand_id)) :
			uniqueIDtoLength [self.auto_strand_id[i] ] =  len(self.auto_sequence_list[i ])
		return  uniqueIDtoLength, self.auto_strand_id
	
	def generate_statespace( self   ):
		#Generates the state space
		state = self.initial_final_state_config()[0]
		self.statespace = []
		self.fast_access= dict()
		self.statespace.append(state)
		self.fast_access[state] = len(self.statespace)- 1
		color= dict()
		color [state] = True
		head =  0
		tail = 1
		while head < tail    :
			state = self.statespace[head ]
			pstates = self.possible_states( state )
			for s in pstates :
				if s not in color  :
					color[s] = True
					self.statespace.append(s )
					self.fast_access[s] = len(self.statespace)- 1
					tail +=1
			head += 1
		return self.statespace
	
	def myFilter(self, state):
		dotparens= dict()
		#if self.mismatch_position != None or self.mismatch_position != "" or self.mismatch_position != -1:
		dotparen_list = self.checkformismatchbasepair(state)
		#else:
		#	dotparen_list= self.dot_paren(state)
		if str(dotparen_list) in dotparens:
			return False
		else :
			dotparens[str(dotparen_list)]  = True
		
		if self.checkfordisconnectedstrand(dotparen_list) == True  :
			return False
		sequence_list= self.sequence(state)
		return sequence_list, dotparen_list
	
	def checkfordisconnectedstrand(self, dotparenlist):
		
		for dotparen in dotparenlist:
			dotparens = dotparen.split('+')
			if len(dotparens) ==1 :
				continue
			else:
				for dp in dotparens:
					counter = 0
					for i in range(len(dp)):
						if dp[i]!= '.' :
							counter +=1
					if counter == 0 :
						return True
		return False
	def checkformismatchbasepair(self, state):
		
		
		def match(sq, sq2):
			if sq == 'A' and sq2 == 'T' or sq == 'T' and sq2 == 'A':
				return True
			
			elif  sq == 'C' and sq2 == 'G' or sq == 'G' and sq2 == 'C':
				return True
			
			else:
				return False
		
		
		
		
		dotparens = self.dot_paren(state)
		sequences=  self.sequence(state , formismatch = True  )
		mystacksq =  []
		mystacki =  []
		newdotparens= []
		
		
		for dp , sq in zip ( dotparens, sequences ) :
			
			dplist= list(dp)
			
			for  i in range (len(dp)):
				
				if dp[i] == "(" :
					mystacksq.append (sq[i])
					mystacki.append(i)
				elif dp[i] == ")" :
					sq2 =  mystacksq.pop( )
					ii = mystacki.pop()
					
					if match(sq[i], sq2)  == False :
						
						dplist[ii]  = '.'
						dplist[i]  = '.'
			dplist= "".join(dplist)
			newdotparens.append(dplist)
		
		return newdotparens
	
	def placeInitialFinal(self, pathway_list):
		#initial state must be in first position of pathway_list  and final state must be in last position of pathway_list
		initial, final= self.initial_final_state_config ()
		i  = self.statespace.index(initial)
		f = self.statespace.index(final)
		if i !=  0 :
			tempi = pathway_list[0 ]
			pathway_list[0 ] = pathway_list[i]
			pathway_list[ i ] = tempi
		size = len(self.statespace) -1
		if f !=  size :
			
			tempf = pathway_list[size ]
			pathway_list[size] =  pathway_list[f]
			pathway_list[ f ] = tempf
	
	def get_pathway(self) :
		pathway_list  = self.get_startStates()
		self.placeInitialFinal(pathway_list)
		
		return pathway_list
	
	def get_initial_final(self):
		pathway_list = self.get_pathway()
		initial_final = []
		initial_final.append(pathway_list[0])
		initial_final.append(pathway_list[len(pathway_list) - 1])
		return initial_final

class PathwayHelix(Pathway):
	def __init__(self ,association ,  strands_list , reaction_type , dataset_name ,  dataset_type ,  cutoff = 1. )  :
		Pathway.__init__(self,   strands_list , reaction_type ,  dataset_type, dataset_name  )
		self.cutoff = cutoff
		self.association = association
		self.strand1 = self.strands_list [0]
		self.strand2 = self.strands_list[1]
		
		self.L = len(self.strand1)
		if self.dataset_name == myenums.DatasetName.REYNALDODISSOCIATE.value :
			self.botleftdangle = strands_list[2]
			self.botrightdangle = strands_list[3]
			self.topleftdangle= ""
			self.toprightdangle= ""
		
		elif self.dataset_name == myenums.DatasetName.DUPIUS2013.value :
			self.botleftdangle = strands_list[2]
			self.botrightdangle = ""
			self.toprightdangle =  ""
			self.topleftdangle= strands_list[3]
		else:
			self.botleftdangle = ""
			self.botrightdangle = ""
			self.topleftdangle= ""
			self.toprightdangle= ""
	
	def get_startStates(self):
		
		self.initial_final_state_config()
		self.generate_statespace()
		self.generate_unaligned()
		
		pathway_list = []
		strand_ids = [100,200]
		counter = -1
		self.UNALIGN_1 = False
		self.UNALIGN_2 = False
		self.UNALIGN_3 = False
		unique_dict  = dict()
		useall  =1
		if useall ==1 :
			listofstates = [self.statespace, self.UNALIGNED_states_list_1, self.UNALIGNED_states_list_2 , self.UNALIGNED_states_list_3]
		#listofstates = [ self.statespace,self.UNALIGNED_states_list_2 ]
		else:
			listofstates = [self.statespace]
		
		for mylist in  listofstates :
			
			counter +=1
			
			
			for state in mylist:
				if counter == 1:
					self.UNALIGN_1 = True
					self.UNALIGN_2 =  False
					self.UNALIGN_3 =  False
				
				elif counter == 2 :
					self.UNALIGN_1 = False
					self.UNALIGN_2 = True
					self.UNALIGN_3 =  False
				elif counter == 3 :
					self.UNALIGN_1 =  False
					self.UNALIGN_2 =  False
					self.UNALIGN_3 = True
				
				output = self.myFilter(state)
				
				if output == False :
					continue
				else :
					sequence_list, dotparen_list  = output
				
				if len(dotparen_list)   == 1 :
					
					#print sequence_list, dotparen_list
					new_complex = makeComplex (sequence_list, dotparen_list[0], strand_ids  )
					
					if str(sequence_list) + "-" +  str(  dotparen_list[0]) + "-" +  str(  strand_ids) not in unique_dict:
						
						pathway_list.append([new_complex] )
						unique_dict [str(sequence_list) + "-" +  str(  dotparen_list[0]) + "-" +  str(  strand_ids)]  = True
				
				
				
				elif len(dotparen_list ) == 2:
					new_complex0 = makeComplex ( [ sequence_list[0] ] , dotparen_list[0], [strand_ids[0]]  )
					new_complex1 = makeComplex ( [ sequence_list[1] ], dotparen_list[1] , [strand_ids[1]]  )
					if   str  ( [ sequence_list[0] ] ) + "-" +  str( dotparen_list[0]) + "-" +  str( [strand_ids[0]] ) + "-" +  str( [ sequence_list[1] ]) + "-" +  str( dotparen_list[1] ) + "-" +  str( [strand_ids[1]] ) not in unique_dict:
						pathway_list.append( [new_complex1, new_complex0] )
						unique_dict[str  ( [ sequence_list[0] ] ) + "-" +  str( dotparen_list[0]) + "-" +  str( [strand_ids[0]] ) + "-" +  str( [ sequence_list[1] ]) + "-" +  str( dotparen_list[1] ) + "-" +  str( [strand_ids[1]] )]  =True
			
			#the code in parent.py assumes the last element in  pathway_list is the final state, so removing it now to add it to the end of the list later on
			if counter== 0 and useall == 1  :
				finalstate=pathway_list[-1]
				pathway_list.remove(finalstate)
		
		
		if useall== 1 :
			pathway_list.append(finalstate )
		
		strand1 = self.topleftdangle +self.strand1 + self.toprightdangle
		strand2= self.botleftdangle +self.strand2 + self.botrightdangle
		
		
		
		self.auto_sequence_list =  [strand2, strand1]
		self.auto_strand_id =   strand_ids
		
		
		
		return pathway_list
	
	def boltzmann_sample_initialstates(self, K )  :
		
		if self.association == False :
			return None
		
		print "Boltzman sampling initial complexes for helix association"
		boltzman_sample_initialstates_list = []
		strand_ids = [100,200]
		strand1 = self.topleftdangle +self.strand1 + self.toprightdangle
		strand2= self.botleftdangle +self.strand2 + self.botrightdangle
		dot_parent1 =sample([strand1], K , material = 'dna')
		dot_parent2 =sample([strand2], K , material = 'dna')
		
		for i in range(len(dot_parent1)):
			new_complex0 = makeComplex ( [strand1] , dot_parent1[i], [strand_ids[0]]  )
			new_complex1 = makeComplex ( [ strand2 ], dot_parent2[i] , [strand_ids[1]]  )
			boltzman_sample_initialstates_list.append( [new_complex1, new_complex0] )
		
		return boltzman_sample_initialstates_list
	
	def generate_unaligned(self):
		self.UNALIGNED_states_list_1= []
		self.UNALIGNED_states_list_2= []
		self.UNALIGNED_states_list_3= []
		#generating unaligned
		for i in range (self.L):
			for j in range(self.L):
				#if self.L- i-1 ==j :
				#	#if i ==j :
				#	continue
				#if self.allowed_s tate_unaligned((i,j)) ==True:
				self.UNALIGNED_states_list_1.append( (i,j))
		
		for i in range (self.L-1):
			for j in range(self.L-1):
				#if self.L- i-1 ==j :
				#	#if i ==j :
				#	continue
				#if self.allowed_state_unaligned((i,j)) ==True:
				self.UNALIGNED_states_list_2.append( (i,j))
		
		for i in range (self.L-2):
			for j in range(self.L-2):
				#if self.L- i-1 ==j :
				#	#if i ==j :
				#	continue
				#if self.allowed_state_unaligned((i,j)) ==True:
				self.UNALIGNED_states_list_3.append( (i,j))
	
	
	
	"""
	def allowed_state_unaligned(self, state):
		i,j = state
		allow= False
		if self.strand1[i]  ==  'A' and self.strand2[j] == 'T':
			allow = True

		if self.strand1[i]  ==  'T' and self.strand2[j] == 'A':
			allow = True

		if self.strand1[i]  ==  'G' and self.strand2[j] == 'C':
			allow = True

		if self.strand1[i]  ==  'C' and self.strand2[j] == 'G':
			allow = True
		return allow
	"""
	
	def possible_states(self, state):
		"""Returns the neighbors of state"""
		i, j= state
		if (i == j):
			states = [(n, n + 1) for n in range(0, self.L)]
		else:
			states = [(i - 1, j),
			          (i + 1, j),
			          (i, j - 1),
			          (i, j + 1)]
		removed = False
		removeList = []
		for s in states :
			if s[0] == s[1] and 0 < s[0] <= self.L  :
				removeList.append((s[0],s[1]))
				removed= True
		for s in removeList:
			states.remove(s )
		if removed == True :
			states.append((0,0))
		return filter(self.allowed_state, states)
	
	def initial_final_state_config(self ):
		"""sets the initial and final state for helix association (association == True) and helix dissociation (association == False ) """
		if self.association == True :
			initialStateConfig = (0, 0 )
			finalStateConfig  = (0, self.L)
		if self.association == False:
			initialStateConfig = (0, self.L )
			finalStateConfig  = (0, 0)
		return [initialStateConfig, finalStateConfig]
	
	def allowed_state(self, state):
		"""Check that a state is allowed."""
		i, j= state
		
		
		allow = 0 <= i <= j <=  self.L
		
		
		
		if  ( j - i ) / self.L > self.cutoff :
			allow =  False
		
		
		return allow
	
	def dot_paren(self, state ):
		"""Returns the structure of the complex in dot-paren notation"""
		L = self.L
		i, j = state
		if self.UNALIGN_1== False and  self.UNALIGN_2== False and  self.UNALIGN_3== False:
			dotparen2 = '.' * len(self.botleftdangle)+'.' * (L - j) + '(' * (j - i) + '.' * i+ '.' * len(self.botrightdangle)
			dotparen1 =   '.' * len(self.topleftdangle) +   '.' * i       + ')' * (j - i) + '.' * (L-j) + '.' * len(self.toprightdangle)
			if i == j:
				return  [  dotparen2  , dotparen1]
			else:
				return  [ dotparen2 + '+' +  dotparen1 ]
		
		elif self.UNALIGN_1== True and  self.UNALIGN_2== False and  self.UNALIGN_3== False:
			
			dotparen2 = '.' * len(self.botleftdangle)+'.' * ( i) + '(' + '.' *  (L-i-1)  + '.' * len(self.botrightdangle)
			dotparen1 =   '.' * len(self.topleftdangle) +   '.' * (j)      + ')'  + '.' * (L-j-1) + '.' * len(self.toprightdangle)
			
			return  [ dotparen2 + '+' +  dotparen1 ]
		
		elif self.UNALIGN_1== False and  self.UNALIGN_2== True and  self.UNALIGN_3== False:
			
			dotparen2 = '.' * len(self.botleftdangle)+'.' * ( i) + '((' + '.' *  (L-i-2)  + '.' * len(self.botrightdangle)
			dotparen1 =   '.' * len(self.topleftdangle) +   '.' * (j)      + '))'  + '.' * (L-j-2) + '.' * len(self.toprightdangle)
			
			return  [ dotparen2 + '+' +  dotparen1 ]
		
		
		elif self.UNALIGN_1== False and  self.UNALIGN_2== False and  self.UNALIGN_3== True:
			
			dotparen2 = '.' * len(self.botleftdangle)+'.' * ( i) + '(((' + '.' *  (L-i-3)  + '.' * len(self.botrightdangle)
			dotparen1 =   '.' * len(self.topleftdangle) +   '.' * (j)      + ')))'  + '.' * (L-j-3) + '.' * len(self.toprightdangle)
			
			return  [ dotparen2 + '+' +  dotparen1 ]
	
	def sequence(self, state , formismatch= False):
		i, j  =state
		strand1 = self.topleftdangle +self.strand1 + self.toprightdangle
		strand2= self.botleftdangle +self.strand2 + self.botrightdangle
		if self.UNALIGN_1==False and self.UNALIGN_2 == False and self.UNALIGN_3 == False :
			
			if i ==j :
				
				return [strand2 , strand1 ]
			else:
				if formismatch == True :
					return [strand2 +'+'+  strand1 ]
				
				else:
					return [strand2 ,  strand1  ]
		
		
		elif self.UNALIGN_1 == True or self.UNALIGN_2 == True or self.UNALIGN_3 == True :
			
			if formismatch == True :
				return [strand2 +'+'+  strand1 ]
			
			else:
				return [strand2 ,  strand1  ]





class PathwayThreewayStrandDisplacement(Pathway):
	def __init__(self, strands_list , reaction_type , dataset_name ,  dataset_type,toehold_length , cutoff = 1.)  :
		Pathway.__init__(self,   strands_list , reaction_type ,  dataset_type, dataset_name  )
		self.toehold_length = toehold_length
		self.cutoff =self.toehold_length * cutoff
		
		self.dangleft= ""
		self.dangleright =""
		self.substrateDangele = ""
		self.hasincumbentDangle = False
		self.incumbentDangle = ""
		
		
		if self.dataset_name ==myenums.DatasetName.ZHANG.value:
			
			self.attacker = strands_list[1] + strands_list[0][:toehold_length]
			self.substrate = Complement(self.attacker)
			self.incumbent =self.attacker[:-toehold_length]
			self.incumbentDangle = strands_list[2]
			self.hasincumbentDangle  = True
		
		elif self.dataset_name ==myenums.DatasetName.MACHINEK.value:
			
			self.incumbentDangle = self.strands_list[0]
			self.incumbent = self.strands_list[1][16:]
			self.substrateDangle = self.strands_list[2]
			
			self.substrate = self.strands_list[3]
			self.attacker = self.strands_list[4]
			self.hasincumbentDangle = True
		
		
		elif self.dataset_name == myenums.DatasetName.REYNALDOSEQUENTIAL.value:
			
			self.attacker= self.strands_list[0]
			self.incumbent= self.strands_list[0]
			self.substrate = Complement(self.attacker)
			self.botleftdangle = self.strands_list[2]
			self.botrightdangle = self.strands_list[3]
		else:
			self.attacker =  strands_list[0]
			self.substrate= strands_list[1]
			self.incumbent= strands_list[2]
		
		
		self.L = len(self.substrate)
		self.n = len(self.incumbent)
		self.m = self.L - self.n
	
	
	def get_startStates(self):
		
		self.initial_final_state_config()
		self.generate_statespace()
		
		pathway_list = []
		#strand ids realattacker=1,  realsubstrate  = 2, realincumbent =3
		for state in self.statespace:
			output = self.myFilter(state)
			if output == False :
				continue
			else :
				sequence_list, dotparen_list  = output
			
			i, j, k , l = state
			
			if ( i == j  ):
				new_complex1 =  makeComplex(sequence_list[1], dotparen_list[1], [100])
				new_complex2 =  makeComplex(sequence_list[0], dotparen_list[0], [200,300])
				pathway_list.append([ new_complex1, new_complex2] )
			elif  k == l  :
				new_complex1 =  makeComplex(sequence_list[1], dotparen_list[1], [300])
				new_complex2 =  makeComplex(sequence_list[0], dotparen_list[0], [100,200])
				pathway_list.append([ new_complex1, new_complex2] )
			else:
				new_complex =makeComplex (sequence_list, dotparen_list[0], [100,200,300] )
				pathway_list.append([new_complex])
		
		realincumbent= self.incumbent
		if self.hasincumbentDangle == True :
			
			realincumbent = self.incumbentDangle + realincumbent
		realsubstrate = self.substrate
		if self.dataset_name == myenums.DatasetName.MACHINEK.value :
			realsubstrate =self.substrateDangle +  realsubstrate
		
		if self.dataset_name == myenums.DatasetName.REYNALDOSEQUENTIAL.value :
			realsubstrate= self.botleftdangle+realsubstrate+ self.botrightdangle
		
		
		self.auto_sequence_list =  [self.attacker,  realsubstrate ,   realincumbent]
		self.auto_strand_id =  [100,200, 300]
		return pathway_list
	
	
	def boltzmann_sample_initialstates(self, K )  :
		print "Boltzman sampling initial complexes for threeway strand displacement"
		boltzman_sample_initialstates_list = []
		strand_ids = [100,200, 300]
		
		realincumbent= self.incumbent
		if self.hasincumbentDangle == True :
			
			realincumbent = self.incumbentDangle + realincumbent
		realsubstrate = self.substrate
		if self.dataset_name == myenums.DatasetName.MACHINEK.value :
			realsubstrate =self.substrateDangle +  realsubstrate
		
		if self.dataset_name == myenums.DatasetName.REYNALDOSEQUENTIAL.value :
			realsubstrate= self.botleftdangle+realsubstrate+ self.botrightdangle
		
		dot_parent1 =sample([self.attacker], K , material = 'dna')
		dot_parent2 =sample([realsubstrate, realincumbent], K , material = 'dna')
		
		for i in range(len(dot_parent1)):
			new_complex1 =  makeComplex([self.attacker], dot_parent1[i], [100])
			new_complex2 =  makeComplex([realsubstrate, realincumbent], dot_parent2[i], [200,300])
			boltzman_sample_initialstates_list.append( [new_complex1, new_complex2] )
		
		return boltzman_sample_initialstates_list
	
	def allowed_state(self, state):
		
		"""Checks that a state is allowed."""
		
		
		i, j, k, l = state
		allow =  0 <= i <= j <= k <= l <= self.L and k >= self.m
		
		
		
		#Further prune the statespace to make computations tractable
		
		# this filter resulted in 0.4 to 0.6 error in pathway prediciton from SSA
		"""if  ( k > self.m  ) and ( i != 0 or j < k- 2) :
			allow =False
		if  (i == j and k == l) :
			allow = False """
		
		#if  ( k > self.m + 20) and ( i != 0 or j < k- 40) :
		#	allow =False
		if  ( k > self.m + 10) and ( i != 0 or j < k- 10) :
			allow =False
		if  (i == j and k == l) :
			allow = False
		
		
		return allow
	
	def possible_states( self, state):
		"""Returns the neighbors of state"""
		
		i, j, k, l = state
		if (i == j):
			states = [ ]
			if self.dataset_name ==myenums.DatasetName.REYNALDOSEQUENTIAL.value :
				states += [(i-1, j, k, l) ,  (i+1 , j,  k , l) , (i, j-1 ,  k , l), (i, j+1 ,  k , l), (i, j,  k- 1 , l) ,(i, j,  k+ 1 , l) , (i, j,  k, l +1)   ]
			else:
				states += [(n, n + 1, k, l) for n in range(0, k)]
		
		else:
			
			states = [(i - 1, j, k, l),
			          (i + 1, j, k, l),
			          (i, j - 1, k, l),
			          (i, j + 1, k, l),
			          (i, j, k - 1, l),
			          (i, j, k + 1, l),
			          (i, j, k, l + 1)]   #(i, j, k, l - 1)]
		
		removed = False
		removeList = []
		for s in states :
			if s[2] == s[3] and  0 <= s[2] < self.L :
				removeList.append((s[0],s[1], s[2] , s[3] ))
				removed= True
		for s in removeList:
			states.remove(s )
		if removed == True :
			states.append((s[0], s[1], self.L , self.L  ))
		
		
		removed = False
		removeList = []
		for s in states :
			if s[0] == s[1] and  0 < s[2] <= self.L :
				removeList.append((s[0],s[1], s[2] , s[3] ))
				removed= True
		for s in removeList:
			states.remove(s )
		if removed == True :
			states.append(( 0 , 0 , s[2],s[3]  ))
		
		
		return   filter(self.allowed_state, states)
	
	def initial_final_state_config(self ):
		"""sets the initial and final state for three-way strand displacement """
		initialStateConfig = (0, 0 , self.m, self.L )
		finalStateConfig = (0, self.L , self.L , self.L )
		return [initialStateConfig, finalStateConfig]
	
	def dot_paren( self, state):
		"""Returns the structure of the complex in dot-paren notation"""
		# Note that if i == j, then the attacker strand is
		# not present. If k == l, the incumbent is not present.
		m, L = self.m, self.L
		i, j, k, l = state
		
		attacker = '.' * (L - j) + '(' * (j - i) + '.' * i
		
		if self.dataset_name ==myenums.DatasetName.REYNALDOSEQUENTIAL.value:
			substrate =  (  '.'*len(self.botleftdangle )+ '.' * i       + ')' * (j - i) + '.' * (k - j) +
			                '(' * (l - k) + '.' * (L - l)  + '.'*len(self.botrightdangle))
		else:
			substrate =  ('.' * i       + ')' * (j - i) + '.' * (k - j) +
			              '(' * (l - k) + '.' * (L - l))
		if self.hasincumbentDangle == True  :
			
			incumbent = '.'*  len(self.incumbentDangle) +  '.' * (L - l) + ')' * (l - k) + '.' * (k - m)
		else :
			incumbent = '.' * (L - l) + ')' * (l - k) + '.' * (k - m)
		
		if self.dataset_name == myenums.DatasetName.MACHINEK.value :
			
			substrate =  len (self.substrateDangle ) * '.' + substrate
		
		if i == j:
			
			return  [ substrate + '+' + incumbent, attacker]
		elif k == l:
			return [  attacker + '+' + substrate, incumbent ]
		else:
			return [ attacker + '+' + substrate + '+' + incumbent ]
	
	
	def sequence(self,  state , formismatch= False ):
		"""Returns the sequence of the complex as NUPACK expects. The
		   first line is the number of independent strands, and the last
		   line determines how many times each strand appears."""
		# Note that if i == j, then the attacker strand is not present.
		# If k == l, the incumbent is not present.
		i, j, k, l = state
		realincumbent= self.incumbent
		if self.hasincumbentDangle == True :
			
			realincumbent = self.incumbentDangle + realincumbent
		realsubstrate = self.substrate
		if self.dataset_name == myenums.DatasetName.MACHINEK.value :
			realsubstrate =self.substrateDangle +  realsubstrate
		
		if self.dataset_name == myenums.DatasetName.REYNALDOSEQUENTIAL.value :
			realsubstrate= self.botleftdangle+realsubstrate+ self.botrightdangle
		
		if i == j:
			if formismatch == True :
				return [   realsubstrate + '+' +  realincumbent  , self.attacker ]
			return [  [  realsubstrate ,   realincumbent ] , [self.attacker] ]
		elif k == l:
			if formismatch == True:
				return   [    self.attacker  + '+' +     realsubstrate  , realincumbent ]
			return   [  [  self.attacker  ,   realsubstrate  ] ,[realincumbent] ]
		else:
			if formismatch == True :
				return [self.attacker+ '+' +  realsubstrate + '+' +    realincumbent]
			return [self.attacker,  realsubstrate ,   realincumbent]

class PathwayHairpin(Pathway):
	
	def __init__(self ,association ,  strands_list , reaction_type , dataset_name ,  dataset_type, cutoff =1. )  :
		Pathway.__init__(self,   strands_list , reaction_type ,  dataset_type, dataset_name  )
		self.cutoff = cutoff
		self.association = association
		self.hairpin = strands_list [0] +strands_list [1]  +Complement(strands_list [0])
		self.m = len(strands_list[0])
		self.L = len(self.hairpin)
	
	def boltzmann_sample_initialstates(self, K )  :
		return None
	
	def get_startStates(self):
		
		self.initial_final_state_config()
		self.generate_statespace()
		
		pathway_list = []
		strand_ids = [1]
		for state in self.statespace:
			output = self.myFilter(state)
			if output == False :
				continue
			else :
				sequence_list, dotparen_list  = output
			
			new_complex = makeComplex (sequence_list, dotparen_list[0],  strand_ids  )
			pathway_list.append([new_complex] )
		
		self.auto_sequence_list =  [self.hairpin ]
		self.auto_strand_id =   strand_ids
		return pathway_list
	
	def possible_states(self, state ):
		"""Returns the neighbors of state"""
		i, j = state
		if (i == j):
			states = [(n, n + 1) for n in range(0, self.m)]
		else:
			states = [(i - 1, j),
			          (i + 1, j),
			          (i, j - 1),
			          (i, j + 1)]
		removed = False
		removeList = []
		for s in states :
			if s[0] == s[1] and 0  < s[0]  <=   self.m  :
				removeList.append((s[0],s[1]))
				removed= True
		for s in removeList:
			states.remove(s )
		if removed == True :
			states.append((0,0))
		return filter(self.allowed_state, states)
	
	def initial_final_state_config(self ):
		"""sets the initial and final state for hairpin closing (association == True) and hairpin opening (association == False ) """
		if self.association == True :
			initialStateConfig = (0, 0 )
			finalStateConfig  = (0, self.m)
		if self.association == False:
			initialStateConfig = (0, self.m )
			finalStateConfig  = (0, 0)
		return [initialStateConfig, finalStateConfig]
	
	
	def allowed_state(self, state):
		"""Check that a state is allowed."""
		i, j = state
		
		allow =0 <= i <= j <= self.m
		
		
		
		if  ( j - i ) / self.m > self.cutoff :
			allow = False
		
		
		return allow
	
	def dot_paren(self, state):
		"""Returns the structure of the complex in dot-paren notation"""
		m, L = self.m, self.L
		i, j = state
		dotPar = '.' * i + '('* (j-i) + '.' * (m - j ) + '.' *  (L - 2*m)  + '.' *( m - j ) + ')'* (j-i) + '.' * i
		return  [dotPar]
	
	def sequence(self  ,state , formismatch= False):
		return [self.hairpin ]


class PathwayBubbleClosing(Pathway):
	def __init__(self , strands_list , reaction_type , dataset_name ,  dataset_type , flurPosition ,  cutoff= 1. )  :
		Pathway.__init__(self,   strands_list , reaction_type ,  dataset_type, dataset_name    )
		self.cutoff = cutoff
		self.flurPosition = flurPosition
		
		self.strand1 = strands_list[0]
		self.loop = strands_list[1]
		self.strand2= Complement(self.strand1)
		self.L = len(self.strand1)
		
		self.i = 0  # the i is always fixed and equal to 0
		self.l = self.L # the l is always fixed  and equal to L
		[self.initialj, self.initialk] = [self.flurPosition, self.flurPosition + 1]
		self.j = self.initialj
		self.k = self.initialk
	
	def boltzmann_sample_initialstates(self, K )  :
		return None
	
	def get_startStates(self):
		
		self.initial_final_state_config()
		self.generate_statespace()
		
		pathway_list = []
		strand_ids = [1]
		for state in self.statespace:
			
			output = self.myFilter(state)
			if output == False :
				continue
			else :
				sequence_list, dotparen_list  = output
			
			new_complex = makeComplex (sequence_list, dotparen_list[0],  strand_ids  )
			pathway_list.append([new_complex] )
		
		self.auto_sequence_list =  [self.strand2 +  self.loop+  self.strand1]
		self.auto_strand_id =   strand_ids
		return pathway_list
	
	
	def possible_states(self, state):
		"""Returns the neighbors of state"""
		
		i, j, k, l = state
		states = [ (i, j - 1, k, l),
		           (i, j + 1, k, l),
		           ( i , j, k - 1 , l),
		           (i , j, k + 1 , l)]
		removed = False
		removeList = []
		for s in states :
			if s[1] == s[2] and s[1] != self.flurPosition :
				removeList.append((s[0],s[1], s[2], s[3]))
				removed= True
		for s in removeList:
			states.remove(s )
		if removed == True :
			states.append((0,self.flurPosition, self.flurPosition, self.L))
		return filter(self.allowed_state, states)
	
	def initial_final_state_config(self ):
		"""sets the initial and final state for bubble closing """
		initialStateConfig =(0,self.initialj, self.initialk, self.L )
		finalStateConfig  = (0,  self.flurPosition , self.flurPosition,  self.L)
		
		return [initialStateConfig, finalStateConfig]
	
	def allowed_state(self, state):
		"""Check that a state is allowed."""
		i, j, k, l = state
		allow = (0 <= i <= j <= self.flurPosition <= k <= l <= self.L ) and (0 != j  and  k != self.L)
		
		
		if  ( j - k ) / self.L > self.cutoff :
			allow  =  False
		
		return allow
	
	def dot_paren(self, state):
		#return the dot paranthesis notation
		L = self.L
		i, j, k , l = state
		strand2 = '.' * (i ) + '(' * (j - i) + '.' * (k - j )+ '(' * (l - k ) + '.' * (L - l )
		strand1 =   '.' * (L-l )       + ')' * (l - k ) + '.' * (k-j) +  ')' * (j-i)  + '.' * i
		loop  = '.' * len(self.loop)
		return   [strand2 + loop  + strand1 ]
	
	def sequence(self, state , formismatch= False):
		"""Returns the sequence of the complex as NUPACK expects. The
		   first line is the number of independent strands, and the last
		   line determines how many times each strand appears."""
		return  [self.strand2 +  self.loop+  self.strand1] 