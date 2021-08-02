from __future__ import division
''''
Created on Oct 2, 2017

    Frits Dannenberg, June 1th, 2018.

    The builder code has two classes: Builder and BuilderRate.
    A Builder will accept references to a function that returns a multistrand options object.
    The sampled states and transitions will be collected.

    The BuilderRate accepts a Builder object and will compute the
    mean first passage time from the starting states until hitting an end state.

'''
import time, copy, os, sys

from multistrand.system import SimSystem
from multistrand.utils import uniqueStateID, seqComplement
from multistrand.options import Options, Literals
from multistrand.experiment import standardOptions, makeComplex

from scipy.sparse import csr_matrix, coo_matrix, csc_matrix
from scipy.sparse.linalg import spsolve, bicg, bicgstab, cg, cgs, gmres, lgmres, qmr, inv

import numpy as np
import zlib
import cPickle as pickle
from subprocess import Popen, PIPE, call

#floatT = np.longdouble
floatT= np.float64

from ast import literal_eval
import random
import hashlib
"""
    Constructs a surface of  N *  ( N -1 ) / 2 points along the reaction frontier.
    This method returns a list of starting states that should be used as initial states for the string method.
    A starting state is a list of complexes.
"""

def hybridizationString(seq):
    
    cutoff = 0.75
    ids = [65, 66]
    
    N = len(seq)
    output = []
    
    """ start with the seperated, zero-bp state """
    dotparen1 = "."*N + "+" + "."*N
    
    complex0 = makeComplex([seq], "."*N, [ids[0]])
    complex1 = makeComplex([seqComplement(seq)], "."*N, [ids[1]])
    seperated = [complex0, complex1]
    
    output.append(seperated)
    
    """
        bps is the number of basepairs between the strands
        offset is the offset of where the pairing starts
    """
    
    dotparen0 = [seq, seqComplement(seq)]
    
    for bps in range(1, N + 1):
        
        for offset in range(N):
            
            dotparen1 = list("."*N + "+" + "."*N)
            
            if (bps + offset < (N + 1)) and bps < np.ceil(cutoff * N):
                #             if (bps + offset < (N + 1)):
                
                for i in range(bps):
                    
                    dotparen1[offset + i] = "("
                    dotparen1[ 2 * N - offset - i ] = ")"
                
                dotparen1 = "".join(dotparen1)
                complex = makeComplex(dotparen0, dotparen1, ids)
                
                output.append([complex])
    
    """ We always require the final state to be the success state"""
    
    complex0 = makeComplex(dotparen0, "("*N + "+" + ")"*N, ids)
    output.append([complex0])
    
    return output


def dissociationString(seq):
    
    ''' Exactly like the association, but swap the initial and final states '''
    
    myList = hybridizationString(seq)
    
    first = myList[0]
    last = myList[-1]
    
    myList[0] = last
    myList[-1] = first
    
    return myList


''' returns a list of strings pairs that represent the hybridization steps
    for toeholds / domains during pathway elaboration method    '''


def weave(M):
    
    output = []
    
    for bps in range(1, M + 1):
        
        for offset in range(M):
            
            leftStr = list("."*M)
            rightStr = list("."*M)
            
            if (bps + offset < (M + 1)):
                
                for i in range(bps):
                    
                    leftStr[offset + i] = "("
                    rightStr[ M - offset - i - 1 ] = ")"
                
                leftStr = "".join(leftStr)
                rightStr = "".join(rightStr)
                output.append([leftStr, rightStr])
    
    return output


def threewaybmString(lefttoe, displace, righttoe):
    
    ''' toehold switch around in the substrate '''
    N = len(displace)
    rT = len(lefttoe)
    lT = len(righttoe)
    
    output = list()
    
    invaderSq = lefttoe + displace + righttoe
    incumbentSq = displace
    substrateSq = seqComplement(invaderSq)
    
    invaderID = 65
    incumbentID = 66
    substrateID = 67
    
    """ start with the separated, zero-bp state """
    complex0 = makeComplex([invaderSq], "."*(lT + N + rT), [invaderID])
    complex1 = makeComplex([incumbentSq, substrateSq], "(" * N + "+" + "."* lT + ")" * N + "." *rT, [incumbentID, substrateID])
    seperated = [complex0, complex1]
    
    output.append(seperated)
    
    ''' left invasion toehold '''
    seqsL = [invaderSq, substrateSq, incumbentSq]
    idsL = [invaderID, substrateID, incumbentID]
    weaving = weave(lT)
    
    for pair in weaving:
        dotparen = "."*(N + rT) + pair[0] + "+" + pair[1] + "("*N + "."*rT + "+" + ")"*N
        output.append([makeComplex(seqsL, dotparen, idsL)])
    
    if lT == 0:
        dotparen = "."*(N + rT) + "" + "+" + "" + "("*N + "."*rT + "+" + ")"*N
    
    ''' the dotparen of the weave is the fully hybridized toehold.
        Toggle the basepairs one step at a time.
    '''
    for invasion in range(N - 1):
        parenList = list(dotparen)
        parenList[rT + N - 1 - invasion] = "("
        parenList[lT + N + rT + 1 + lT + invasion ] = ")"
        parenList[-invasion - 1] = "."
        dotparen = "".join(parenList)
        output.append([makeComplex(seqsL, dotparen, idsL)])
    
    ''' right invasion toehold '''
    seqsR = [invaderSq, incumbentSq, substrateSq]
    idsR = [invaderID, incumbentID, substrateID]
    weaving = weave(rT)
    
    for pair in weaving:
        dotparen = pair[0] + "."*(N + lT) + "+" + "("*N + "+" + "."*lT + ")"*N + pair[1]
        output.append([makeComplex(seqsR, dotparen, idsR)])
    
    if rT == 0:
        dotparen = "" + "."*(N + lT) + "+" + "("*N + "+" + "."*lT + ")"*N + ""
    
    '''
        Toggle the basepairs one step at a time.
    '''
    for invasion in range(N - 1):
        parenList = list(dotparen)
        parenList[rT + invasion] = "("
        parenList[lT + N + rT + 1 + invasion ] = "."
        dotparen = "".join(parenList)
        output.append([makeComplex(seqsR, dotparen, idsR)])
    
    ''' Do not forget to set the final state.
        This is just the displaced strand floating freely.
        '''
    output.append([makeComplex([incumbentSq], "."*N, [incumbentID])])
    
    return output


class ConvergeCrit(object):
    
    period = 4  # average out over past X increases.
    
    def __init__(self):
        
        self.maximumIterations = 1000
        self.minimumStateIncrement = 4
        
        self.currIteration = 0
        self.currStates = -99
        
        self.precision = 0.05
        self.array = [-99.0] * self.period
    
    def converged(self, rateIn=None, statespace_size=None):  # also saves the rateIn
        
        if self.currIteration > self.maximumIterations:
            return True
        
        if self.precision < 1.0 :
            averageL = (sum(self.array) / self.period) * (1.0 - self.precision)
            averageH = (sum(self.array) / self.period) * (1.0 + self.precision)
            
            conv1 = rateIn > averageL
            conv2 = rateIn < averageH
            
            self.array[self.currIteration % self.period] = rateIn
            self.currIteration += 1;
            return (conv1 and conv2)
        
        else:
            
            if statespace_size - self.currStates < self.minimumStateIncrement:
                return True
            
            self.currStates = statespace_size;
            
            return statespace_size >= self.precision
    
    def __str__(self):
        
        return  str(self.array)


class transitiontype(object):
    
    unimolecular = "uni"
    bimolecularIn = "bi-in"
    bimolecularOut = "bi-out"
    
    array = [unimolecular, bimolecularIn, bimolecularOut ]


class localtype(object):
    
    end = "End"
    loop = "Loop"
    stack = "Stack"
    stackstack = "StackStack"
    loopend = "LoopEnd"
    stackend = "StackEnd"
    stackloop = "StackLoop"
    
    array = [end, loop, stack, stackstack, loopend, stackend, stackloop]


class Energy(object):
    
    GAS_CONSTANT = floatT(0.0019872036)  # kcal / K mol
    
    dH = 0.0;
    dS = 0.0;
    
    def __init__(self, dG, dH, temp, concentration, n_complexes, n_strands):
        
        assert(200 < temp < 400)
        
        RT = self.GAS_CONSTANT * floatT(temp);
        dG_volume = RT * (n_strands - n_complexes) * np.log(1.0 / concentration)
        
        self.dH = floatT(dH)
        #print "dG_volume " , dG_volume ,  "n_strands" , n_strands , "n_complexes", n_complexes  , "concentration" , concentration  , "dG" , dG
        self.dS = floatT(-(dG  - dG_volume - dH) / temp)
    
    def dG(self, temp):
        
        return self.dH - floatT(temp) * self.dS
    
    def dG_wihtoutanyvolume(self, temp ):
        dS = floatT(-(dG  - dH) / temp)
        return self.dH - floatT(temp) * self.dS
    
    def __str__(self):
        
        return "dH= " + str(self.dH) + " dS= " + str(self.dS)
    
    def __eq__(self, that):
        if isinstance(that, Energy):
            return self.dH == that.dH and self.dS == that.dS
        else:
            raise Exception("Not an acceptible argument", "__eq__", that)


def codeToDesc(code):
    
    output = []
    primes = [3, 5, 7, 11, 13, 17, 19]
    primesq = [9, 25, 49, 121, 169, 289, 361]
    
    for i in range(len(primes)):
        
        if code % primes[i] == 0:
            output.append(localtype.array[i])
            
            if code % primesq[i] == 0:  # code is a square number
                output.append(localtype.array[i])
                return output
    
    return output


class InitCountFlux(object):
    
    def __init__(self):
        self.count = 0;
        self.flux = 0.0;
    
    def __add__(self, x):
        self.count + x.count
    
    def __str__(self):
        return "Count: " + str(self.count) + "  flux: " + str(self.flux)
    
    def  __repr__(self):
        return str(self)


class Builder(object):
    
    verbosity = False
    
    # input function returns the multistrand options object for which to build the statespace
    def __init__(self, inputFunction, arguments):
        
        # initial argument has to be the number of trials
        self.optionsFunction = inputFunction
        self.optionsArgs = copy.deepcopy(arguments)
        
        self.doMultiprocessing = False
        self.printTimer = True
        self.numOfThreads = 8
        
        self.protoinitialpathwayuniqueIDs = dict()
        self.protoSpace = dict()  # key: states. Value: Energy
        self.protoTransitions = dict()  # key: transitions. Value: ArrheniusType (negative if it is a bimolecular transition)
        self.protoInitialStates = dict()  # key: states. Value: a InitCountFlux object that tells how many times the state has been the initial state and the join flux (rate)
        self.protoFinalStates = dict()  # key: states: Value: the result of this final state can be SUCCES or FAILURE
        self.protoSequences = dict()  # key: name of strand. Value: sequence
        
        self.firstStepMode = True
        self.startTime = time.time()
        
        self.mergingCounter = 0  # counts how many transitions from the merging have been found
        
        # save a copy for later processing -- note this copy will not have results attached to it,
        # so it won't have a large memory footprint.
        self.options = self.optionsFunction(self.optionsArgs)
        
        self.the_dir = "p_statespace/"
        self.deltaPruningValue = 0
        self.matrixsolvetime = -1
        
        self.fatten  =  False
        self.multistrand_simulation_time= 0
        self.loading_protospace_time  = 0
        self.loading_prototransitions_time  = 0
        self.loading_protoinitialstates_time = 0
        self.loading_protofinalstates_time  = 0
        self.multistrand_simulation_fatten_time= 0
        self.loading_protospace_fatten_time  = 0
        self.loading_prototransitions_fatten_time  = 0
        self.loading_protoinitialstates_fatten_time = 0
        self.loading_protofinalstates_fatten_time  = 0
    
    
    
    def __str__(self):
        
        output = "states / transitions / initS / finalS / merged     \n "
        output += str(len(self.protoSpace))
        output += "   -    " + str(len(self.protoTransitions))
        output += "   -    " + str(len(self.protoInitialStates))
        output += "   -    " + str(len(self.protoFinalStates))
        output += "   -    " + str(self.mergingCounter)
        
        return output
    
    def reset(self):
        
        self.protoSpace.clear()
        self.protoTransitions.clear()
        self.protoInitialStates.clear()
        self.protoFinalStates.clear()
        self.protoSequences.clear()
    
    def printOverlap(self, other):
        
        N = len(self.protoSpace)
        
        overlap = 0;
        
        for state in self.protoSpace:
            if state in other.protoSpace:
                overlap += 1
        
        print "The overlap is " + str(100.0 * overlap / N) + " percent. "
    
    def mergeSet(self, this, that):
        
        for key, val in that.iteritems():
            
            if not key in this:
                this[key] = val;
    
    
    def mergeSet_2(self, this, that):
        
        for key, val in that.iteritems():
            
            this[key] = val;
    
    def mergeBuilder(self, other):
        self.mergeSet ( self.protoinitialpathwayuniqueIDs, other.protoinitialpathwayuniqueIDs )
        self.mergeSet(self.protoSpace, other.protoSpace)
        self.mergeSet(self.protoTransitions, other.protoTransitions)
        self.mergeSet(self.protoInitialStates, other.protoInitialStates)
        self.mergeSet(self.protoFinalStates, other.protoFinalStates)
        self.mergeSet(self.protoSequences, other.protoSequences)
    
    ''' Merges if both source and the target exist '''
    def transitionMerge_2(self, protoTransitions):
        
        for key, value in protoTransitions.iteritems():
            
            sFrom = key[0]
            sTo = key[1]
            
            if sFrom in self.protoSpace and sTo in self.protoSpace:
                
                if not key in self.protoTransitions:
                    self.protoTransitions[key] = value
                    self.mergingCounter += 1
    
    def transitionMerge(self, other):
        
        for key, value in other.protoTransitions.iteritems():
            
            sFrom = key[0]
            sTo = key[1]
            
            if sFrom in self.protoSpace and sTo in self.protoSpace:
                
                if not key in self.protoTransitions:
                    self.protoTransitions[key] = value
                    self.mergingCounter += 1
    
    def parseState(self, line, simulatedTemperature, simulatedConc):
        
        mywords = line.split()
        
        n_complexes = int(mywords[0])
        n_strands = 0
        
        ids = []
        sequences = []
        structs = []
        
        for i in range(n_complexes):
            
            ids.append(mywords[1 + i])
            sequences.append(mywords[1 + n_complexes + i])
            structs.append(mywords[1 + 2 * n_complexes + i])
            n_strands += len(mywords[1 + 2 * n_complexes + i].split('+'))
            
            uniqueID2 = uniqueStateID(ids, structs)
            uniqueID = ""
            
            for iddd in uniqueID2 :
                uniqueID += str(iddd)
        
        dG = float(mywords[1 + 3 * n_complexes])
        dH = float(mywords[1 + 3 * n_complexes + 1])
        
        energyvals = Energy(dG, dH, simulatedTemperature, simulatedConc, n_complexes, n_strands)
        
        return uniqueID, energyvals, (sequences, ids, structs)
    
    """ Runs genAndSavePathsFile until convergence is reached"""
    
    def genUntilConvergence(self, precision):
        
        crit = ConvergeCrit()
        crit.precision = precision
        
        currTime = -1.0
        
        while not crit.converged(currTime, len(self.protoSpace)):
            self.genAndSavePathsFile()
            
            if self.verbosity:
                print "Size     = %i " % len(self.protoSpace)
            
            if precision < 1.0:
                builderRate = BuilderRate(self)
                currTime = builderRate.averageTimeFromInitial()
        
        self.fattenStateSpace()
        
        if self.verbosity:
            print "Size     = %i " % len(self.protoSpace)
    
    """ Runs genAndSavePathsFile until convergence is reached,
        given a list of initial states"""
    
    def bfs(self,  mymethod, initialstate,  pathwaystates,myhash  ,  finalstate, keyoffinalstate, structureoffinalstate , extrastates ) :
        
        state = initialstate
        color= dict()
        color [str(initialstate)] = True
        while True    :
            if  state ==None :
                break
            foundequal  = False
            mykey = ""
            for p1 in range(len( state )) :
                mykey = mykey +"-" +state[ p1 ].sequence +"-"+ state[ p1 ].structure + "-"
            if mykey in myhash:
                foundequal  = True
            myhash[mykey] = True
            if foundequal == False :
                pathwaystates.append(state)
            pstate, extrastates_temp  = self.possible_states( state , finalstate,    findkeyandstructureforfinal = False , keyoffinalstate = keyoffinalstate, structureoffinalstate= structureoffinalstate)
            state= pstate
    
    def generateInitialPathwayAutomatically(self, initialstate, finalstate,  printMeanTime = False, pathwayelaboration_N = 1 , pathwayelaboration_beta = 0 , uniqueIDtoLength = "mustset" , auto_strand_id_list = "mustset" , boltzmann_sample_initialstates_list  =  None , pathwayelaboration_use_elaboration =  True  ) :
        print "Generating biased paths",  "N is: ", pathwayelaboration_N, "  beta is: ", pathwayelaboration_beta
        self.use_extrastates = False # don't change this
        self.pathwayelaboration_beta= pathwayelaboration_beta # if you put self.pathwayelaboration_beta = 0 , then misaligned base pairs will not be allowed in initial random pathways
        self.uniqueIDtoLength= uniqueIDtoLength
        self.auto_strand_id_list= auto_strand_id_list
        self.pathwayelaboration_use_elaboration = pathwayelaboration_use_elaboration
        keyoffinalstate, structureoffinalstate = self.possible_states(finalstate, finalstate ,  findkeyandstructureforfinal = True )
        breakloop = "not even important "
        pathwaystates= []
        starttime =time.time( )
        myhash = dict()
        mymethod = 2
        if boltzmann_sample_initialstates_list == None :
            boltzmann_sample_initialstates_list = []
            for repeat in range(pathwayelaboration_N) :
                boltzmann_sample_initialstates_list.append(initialstate)
        listofinitialstates = []
        
        for initialstate in boltzmann_sample_initialstates_list :
            mykey = ""
            foundequal  = False
            for p1 in range(len( initialstate  )) :
                mykey = mykey +"-" +initialstate[ p1 ].sequence +"-"+ initialstate[ p1 ].structure + "-"
            if mykey in myhash:
                foundequal  = True
            myhash[mykey] = True
            if foundequal == False :
                pathwaystates.append(initialstate)
                listofinitialstates.append(initialstate)
        
        #pathwaystates.append(finalstate)  #do not remove
        mykey = ""
        for p1 in range(len( finalstate  )) :
            mykey = mykey +"-" +finalstate[ p1 ].sequence +"-"+ finalstate[ p1 ].structure + "-"
        myhash[mykey] = True
        pathwaystates.append(finalstate)
        
        extrastates= []
        
        for initialstate in boltzmann_sample_initialstates_list :
            self.chosenkeys_in_path= dict( )
            avoidMemoryLeak= True
            if avoidMemoryLeak== False:
                self.bfs(  mymethod , initialstate,  pathwaystates,myhash  , finalstate, keyoffinalstate, structureoffinalstate , extrastates )
            else:
                buildername =  self.outfile.split("/")[-1]
                pathfile = self.outfile
                path= ["","","","","","","","",""]
                for i in range(7):
                    path[i] = pathfile +"-"+ str(i)+"-" + buildername
                
                with open(path[0], "wb") as p :
                    pickle.dump( initialstate,p )
                with open(path[1], "wb") as p :
                    pickle.dump( pathwaystates, p)
                with open(path[2], "wb") as p :
                    pickle.dump( myhash,p)
                with open(path[3], "wb") as p :
                    pickle.dump( finalstate,p)
                with open(path[4], "wb") as p :
                    pickle.dump( keyoffinalstate ,p)
                with open(path[5], "wb") as p :
                    pickle.dump( structureoffinalstate ,p)
                with open(path[6], "wb") as p :
                    pickle.dump( extrastates ,p)
                builderpickletemp = pathfile + "temp"+ buildername
                with open(builderpickletemp, "wb" ) as p :
                    newbuilder=   Builder(self.optionsFunction, self.optionsArgs)
                    newbuilder.pathwayelaboration_beta = self.pathwayelaboration_beta
                    newbuilder.uniqueIDtoLength= self.uniqueIDtoLength
                    newbuilder.auto_strand_id_list = self.auto_strand_id_list
                    newbuilder.chosenkeys_in_path = self.chosenkeys_in_path
                    newbuilder.use_extrastates= self.use_extrastates
                    
                    pickle.dump(  newbuilder, p)
                
                command = ["python", "helper_generatebiasedpaths.py" ,  str(mymethod), path[0], path[1], path[2], path[3], path[4], path[5], path[6],  builderpickletemp ]
                shell = call(command )
                del shell
                
                
                with open(path[1], "rb") as p :
                    pathwaystates=   pickle.load(p)
                with open(path[2], "rb") as p :
                    myhash=   pickle.load(p)
                
                with open(path[6], "rb") as p :
                    extrastates =   pickle.load(p)
        
        
        
        #if you set avoidmemoryleak to True then the complex objects in listofinitialstates will be different pathwaystate!
        if avoidMemoryLeak == True:
            for i in range(7):
                os.remove(path[i])
            os.remove(builderpickletemp)
            
            
            lenlistofinitialstates= len(listofinitialstates)
            for i in range(lenlistofinitialstates):
                initstate = listofinitialstates[i]
                mykeyinit= ""
                for p1 in range(len( initstate )) :
                    mykeyinit = mykeyinit +"-" +initstate[ p1 ].sequence +"-"+ initstate[ p1 ].structure + "-"
                for ps2 in pathwaystates :
                    mykey = ""
                    for p1 in range(len( ps2)) :
                        mykey = mykey +"-" +ps2[ p1 ].sequence +"-"+ ps2[ p1 ].structure + "-"
                    if mykeyinit ==mykey :
                        listofinitialstates[i] = ps2
                        break
        #print "total number of states from biased paths ", len(pathwaystates )
        
        #The following lines are not used becaues self.use_extrastates is set to False. Do not change it
        extrastates2 = []
        if self.use_extrastates  == True :
            for state in extrastates:
                foundequal = False
                for ps2 in pathwaystates:
                    if len(state )!= len(ps2):
                        continue
                    for p1 in range(len( state )) :
                        if state[ p1 ].sequence == ps2[p1 ].sequence and state[ p1 ].structure == ps2 [p1 ].structure  :
                            if p1 == len(state)-1 :
                                foundequal = True
                                break
                if foundequal == False :
                    extrastates2.append(state)
        
        print "Finished generating biased paths"
        print "Going to elaborate state space"
        self.genAndSavePathsFromString(pathwaystates, printMeanTime=printMeanTime, extrastates=extrastates2 ,  listofinitialstates = listofinitialstates)
        print "Finished elaborating state space"
        return breakloop, len(self.protoSpace )
    
    
    def possible_states(self, currentState , finalstate , findkeyandstructureforfinal = False , keyoffinalstate = None , structureoffinalstate=  None) :
        #print "in posisble _states "
        #print currentState
        #for state in currentState:
        #   print state
        #print "+++++++++"
        ogVerb = Builder.verbosity
        Builder.verbosity = False
        counter = 0
        
        def inspectionSim(inputs):
            
            o1 = standardOptions()
            
            #added by nz, these lines should be here or the energy of the states will be different
            o1.sodium = self.options.sodium
            o1.magnesium = self.options.magnesium
            o1.temperature = self.options.temperature
            
            #o1.join_concentration = self.options.join_concentration  #do not set the concentration, becuase for low concentration this line will make the generation of initial pathway very slow
            ## addby nz
            
            o1.rate_method = self.options.rate_method
            o1.start_state = inputs[0]
            
            return o1
        
        #nz you may have to change this line
        myB = Builder(inspectionSim, [ currentState])
        myB.fatten = True # do not remove this
        myB.genAndSavePathsFile(inspecting=True)
        #self.transitionMerge(myB)
        #self.mergeSet_2(self.protoTransitions, myB.protoTransitions)
        
        #print "len is " , len(myB.protoSpace )
        if len(myB.options.start_state) == 3:
            initialstatestructure1 = [  myB.options.start_state[0].structure,myB.options.start_state[1].structure, myB.options.start_state[2].structure ]
            initialstatestructure2 = [  myB.options.start_state[0].structure,myB.options.start_state[2].structure, myB.options.start_state[1].structure ]
            initialstatestructure3 = [  myB.options.start_state[1].structure,myB.options.start_state[0].structure, myB.options.start_state[2].structure ]
            initialstatestructure4 = [  myB.options.start_state[1].structure,myB.options.start_state[2].structure, myB.options.start_state[0].structure ]
            initialstatestructure5 = [  myB.options.start_state[2].structure,myB.options.start_state[1].structure, myB.options.start_state[0].structure ]
            initialstatestructure6 = [  myB.options.start_state[2].structure,myB.options.start_state[0].structure, myB.options.start_state[1].structure ]
        
        elif len(myB.options.start_state) == 2:
            initialstatestructure1 = [  myB.options.start_state[0].structure,myB.options.start_state[1].structure ]
            initialstatestructure2 = [  myB.options.start_state[1].structure,myB.options.start_state[0].structure ]
            initialstatestructure3  = None
            initialstatestructure4  = None
            initialstatestructure5  = None
            initialstatestructure6  = None
        
        elif len(myB.options.start_state) == 1:
            initialstatestructure1  = [myB.options.start_state[0].structure]
            initialstatestructure2  = None
            initialstatestructure3  = None
            initialstatestructure4  = None
            initialstatestructure5  = None
            initialstatestructure6  = None
        else:
            raise ValueError("myB.options.start_state  is something you did not think of ")
        
        keys = dict()
        
        structures= dict()
        for key , value   in myB.protoSpace.iteritems():
            
            try :
                newkey = literal_eval(key)
                newkey = (newkey, )
            except :
                newkey = ""
                for i in range (len(key)-1):
                    if key[i] ==")" and key[i+1] == "(":
                        newkey += key[i] +","
                    else:
                        newkey += key[i]
                newkey+= key[-1]
                newkey = literal_eval(newkey)
            keys[ key ] = newkey
            myT = myB.options._temperature_kelvin
            dG1 = myB.protoSpace[key].dG(myT)
            (seqs, ids, structs) = myB.protoSequences[key]
            structures [key ] =  structs
            if structs  == initialstatestructure1 or structs  == initialstatestructure2 or structs  == initialstatestructure3 or structs  == initialstatestructure4 or structs  == initialstatestructure5 or structs  == initialstatestructure6:
                energyofstate  = dG1
                keyofstate=   key
        
        
        if findkeyandstructureforfinal  == True   :
            return keys[keyofstate],structures[keyofstate]
        
        
        if keys [keyofstate] == keyoffinalstate :
            statechoice =  None
            extrastates = []
            return statechoice , extrastates
        
        statespace_lower = []
        statespace_greater =  []
        
        statespace_greater_extra = [ ]
        
        
        builderrate= BuilderRate(myB, ignoremessage = True )
        
        allneighbors = []
        for key , value   in myB.protoSpace.iteritems():
            if key == keyofstate :
                continue
            
            myT = myB.options._temperature_kelvin
            
            dG1 = myB.protoSpace[key].dG(myT)
            (seqs, ids, structs) = myB.protoSequences[key]
            
            myState = []
            for seq, id, struct in zip(seqs, ids, structs):
                seqs = seq.split('+')
                
                newid_temp = id.split(',')
                newid= []
                # Make sure the ids are the same ids are the initial -- should be int becaues in pathways you gave int and should aslso be split as below
                for i in newid_temp:
                    newid.append(int(i.split("_")[-1]))
                
                myC = makeComplex(seqs, struct, newid)
                myState.append(myC)
            
            allneighbors.append(myState)
            myrate=  builderrate.get_rate( keyofstate, key)
            statespace_lower.append( ( keys[key], structures[key], myState, dG1, myrate ))
        
        statespace_lower_modified = []
        statespace_greater_modified = []
        probabilities_lower_modified = [ ]
        probabilities_greater_modified = [ ]
        
        probability_lower_sum  = 0
        probability_greater_sum  = 0
        timetochooseunaligned = False
        if self.pathwayelaboration_beta >0 :
            myrandomnumber= random.random()
            if myrandomnumber  < self.pathwayelaboration_beta:
                timetochooseunaligned = True
        
        for state in statespace_lower:
            if  self.allowed_state(state[0], state[1],  keys[keyofstate],structures[keyofstate],   keyoffinalstate , structureoffinalstate , timetochooseunaligned  )   == True :
                statespace_lower_modified.append(state)
                probabilities_lower_modified.append(state[4][0])
                probability_lower_sum+= state[4][0]
        
        for state in statespace_greater:
            if  self.allowed_state(state[0], state[1 ]  ,keys[keyofstate],structures[keyofstate],   keyoffinalstate  , structureoffinalstate   , timetochooseunaligned ) == True :
                statespace_greater_modified.append(state)
                probabilities_greater_modified.append(state[4][0])
                probability_greater_sum+= state[4][0]
                statespace_greater_extra.append( state[2] )
        
        
        
        if len ( statespace_lower_modified)!=  0 :
            
            probabilities_lower_modified = [i/probability_lower_sum for i in probabilities_lower_modified]
            choice= np.random.choice(range(len(probabilities_lower_modified)), p= probabilities_lower_modified)
            choice= statespace_lower_modified[choice]
            statechoice = choice[2]
            extrastates = []
        elif len (statespace_greater_modified)!= 0 :
            probabilities_greater_modified = [i/probability_greater_sum for i in probabilities_greater_modified]
            choice= np.random.choice(range(len(probabilities_greater_modified)), p = probabilities_greater_modified)
            choice= statespace_greater_modified[choice]
            statechoice = choice[2]
            extrastates=  allneighbors
            #extrastates  =  statespace_greater_extra
        else :
            print " \n \n in here only once \n\n    "
            statechoice =  None
            extrastates = []
        
        del myB
        if statechoice !=  None:
            self.chosenkeys_in_path [choice[0] ]  =True
        return statechoice , extrastates
    
    
    def allowed_state(self, key1, struct1 ,    key2,  struct2,  finalkey , finalstruct ,  timetochooseunaligned = False , reportdistance= False  , distancebasedonbasepairswithsubstrate = False ) :
        
        def hamming_distance(chaine1, chaine2):
            d1 = dict ( )
            d2  = dict ( )
            count = 0
            for chain   in [chaine1,chaine2 ]:
                count +=1
                for i in range(len(chain )  ) :
                    
                    if count ==1:
                        if chain[i]!= 0 : #because we are only considering base pairs
                            d1[ i ] = chain[i ]
                    elif count ==2:
                        if chain[i] != 0 : #because we are only considering base pairs
                            d2[ i ]  = chain[i ]
            #d1_keys = set(d1.keys())
            #d2_keys = set(d2.keys())
            numofunpaired = 0
            for i in range(len( chaine1)) :
                if chaine1[i] == 0 :
                    numofunpaired +=1
            set1 = set(d1.items())
            set2 = set(d2.items())
            return len( set1 ^ set2 ) , numofunpaired
        
        def dbps( key ) :
            #only works for machinek handcoded with ids 100,200,300    - This function is used in saveStatisticforContourPlot function !!!!
            key = changeOrder(key)
            attackerlength = self.uniqueIDtoLength[int(100)]
            substratelength = self.uniqueIDtoLength[int(200)]
            attackerbasepair= 0
            incumbentbasepair= 0
            #print "key is " , key  , len(key ), attackerlength , attackerlength+ substratelength
            for i in range(attackerlength, attackerlength+ substratelength  ):
                if 1<= key[i] <= attackerlength  :
                    attackerbasepair+=1
                elif key[i] >= attackerlength+ substratelength:
                    incumbentbasepair+=1
            #print   "attackerbasepair , incumbentbasepair  ", attackerbasepair , incumbentbasepair
            return attackerbasepair , incumbentbasepair
        
        def changeOrder(key ) :
            newkey2 = []
            newkey = dict()
            for i in key :
                newkey[i ] = []
                uniqueids = i[0].split("a_")[1:]
                point = 0
                for id in uniqueids :
                    lengthid  = self.uniqueIDtoLength[int(id)]
                    substract = point
                    addition = 0
                    for orderid  in self.auto_strand_id_list:
                        
                        if int(orderid) == int(id):
                            break
                        addition += self.uniqueIDtoLength[int(orderid) ]
                    for j in range(point, point + lengthid ):
                        newkey[i].append( j - substract + addition +1)
                    point= point + lengthid
                for j in  range(len(newkey[i])):
                    if i[1][j]  == 0 :
                        newkey2.append(0 )
                    else:
                        newkey2.append(newkey[i][ i[1][j] -1])
            
            return tuple( newkey2 )
        
        def distance2( key, finalkey  ) :
            changedkey  = changeOrder(key)
            changedfinalkey = changeOrder(finalkey )
            return hamming_distance(changedkey , changedfinalkey  )
        
        
        if distancebasedonbasepairswithsubstrate == True :
            
            return dbps(key1 )
        
        mydistance1, numofzeros1=  distance2( key1,  finalkey)
        
        if reportdistance == True :
            
            return mydistance1
        
        mydistance2, numofzeros2= distance2(key2,  finalkey)
        if  timetochooseunaligned == True and self.pathwayelaboration_beta > 0 :
            return True
        else:
            if mydistance1 <  mydistance2 :
                return True
            else:
                return False
    
    
    def genUntilConvergenceWithInitialState(self, precision, initialStates, printMeanTime=False):
        
        crit = ConvergeCrit()
        crit.precision = precision
        
        currTime = -1.0
        
        while not crit.converged(currTime, len(self.protoSpace)) :
            
            self.genAndSavePathsFromString(initialStates, printMeanTime=printMeanTime)
            
            if precision < 1.0:
                builderRate = BuilderRate(self)
                currTime = builderRate.averageTimeFromInitial()
            
            if printMeanTime:
                print "Mean first passage time = %.2E" % currTime
        
        self.fattenStateSpace()
        
        if self.verbosity:
            print "Size     = %i " % len(self.protoSpace)
    
    ''' A single iteration of the pathway elaboration method '''
    
    def genAndSavePathsFromString(self, pathway, printMeanTime=False, extrastates  =[] ,listofinitialstates = None):
        
        """ Only the first state will count towards the set of initial states """
        for state in extrastates:
            pathway.append(state)
        for state in pathway:
            if state in listofinitialstates :
                ignoreInitial = False
            
            else :
                ignoreInitial = True
            
            startTime = time.time()
            otherBuilder = Builder(self.optionsFunction, self.optionsArgs)
            #print "Number of simulations per state is: ", otherBuilder.options.num_simulations
            #print "Simulation time per state is: ",  otherBuilder.options.simulation_time
            otherBuilder.genAndSavePathsFile(supplyInitialState=state, ignoreInitialState=ignoreInitial)
            #ignoreInitial = True
            self.mergeBuilder(otherBuilder)
            self.multistrand_simulation_time += otherBuilder.multistrand_simulation_time
            self.loading_protospace_time  += otherBuilder.loading_protospace_time
            self.loading_prototransitions_time  += otherBuilder.loading_prototransitions_time
            self.loading_protoinitialstates_time += otherBuilder.loading_protoinitialstates_time
            self.loading_protofinalstates_time  += otherBuilder.loading_protofinalstates_time
            del otherBuilder
        
        if self.verbosity or printMeanTime:
            print "Size     = %i    ---  bytesize = %i " % (len(self.protoSpace), sys.getsizeof (self.protoSpace))
            print "Size T   = %i    ---  bytesize = %i " % (len(self.protoTransitions), sys.getsizeof (self.protoTransitions))
            print "Time = %f" % (time.time() - startTime)
    
    
    """Using this funciton to add missing transition between states of the pathway"""
    def fattenStateSpace(self, start  =   None , end = None , onlyfatteninitialpathway = 0  ):
        print "In fattenStateSpace.py"
        ogVerb = Builder.verbosity
        Builder.verbosity = False
        counter = 0
        
        def inspectionSim(inputs):
            
            o1 = standardOptions()
            
            #added by nz, these lines should be here or the energy of the states will be different
            o1.sodium = self.options.sodium
            o1.magnesium = self.options.magnesium
            o1.temperature = self.options.temperature
            ## added by nz
            
            o1.rate_method = self.options.rate_method
            o1.start_state = inputs[0]
            
            return o1
        self.tempstatespace = dict()
        self.temptransitions = dict()
        
        #for key, value in self.protoSpacebackup.iteritems():
        if onlyfatteninitialpathway  ==1 :
            iterlist= self.protoinitialpathwayuniqueIDsbackup
        else:
            iterlist= self.protoSpacebackup
        for key, value in iterlist.iteritems():
            #if onlyfatteninitialpathway == True :
            #    if key not in self.protoinitialpathwayuniqueIDsbackup:
            #       continue
            
            if start != None and end != None:
                if counter < start :
                    counter += 1
                    continue
                if counter >= end:
                    counter +=1
                    break
            ogVerb = True
            #if ogVerb and ((counter % 100) == 0):
            #print "Searching for missing transitions. Progress " + str(counter) +   " " + str(end) +  " / " + str(len(self.protoSpace))
            
            (seqs, ids, structs) = self.protoSequences[key]
            
            myState = []
            
            for seq, id, struct in zip(seqs, ids, structs):
                
                seqs = seq.split('+')
                ids = id.split(',')
                myC = makeComplex(seqs, struct, ids)
                
                myState.append(myC)
            
            ''' post: myState is the state we want to explore transitions for. '''
            
            myB = Builder(inspectionSim, [myState])
            myB.fatten = True # do not remove this
            myB.genAndSavePathsFile(inspecting=True)
            
            self.multistrand_simulation_fatten_time += myB.multistrand_simulation_fatten_time
            self.loading_protospace_fatten_time  +=  myB.loading_protospace_fatten_time
            self.loading_prototransitions_fatten_time  += myB.loading_prototransitions_fatten_time
            self.loading_protoinitialstates_fatten_time += myB.loading_protoinitialstates_fatten_time
            self.loading_protofinalstates_fatten_time  +=  myB.loading_protofinalstates_fatten_time
            
            #self.transitionMerge(myB)
            self.mergeSet_2(self.protoTransitions, myB.protoTransitions)
            
            counter += 1
            #ft  = time.time()
            #if counter % 100 == 0 :
            #    print "time for 100" , ft-st
        
        print "Finished fattenStateSpace.py"
        Builder.verbosity = ogVerb
    
    """ Computes the mean first pasasage times,     then selects states that are delta-close     to the set of final states.     Those states are then added to the set of final states.     This reduces the size of the matrix that is constructed.  """
    def deltaPruning(self, delta=0.01,  times = None  ,averagedMFPT = None  , printCount=False , maxiter= None):
        builderRate = BuilderRate(self)
        if averagedMFPT== None :
            times = builderRate.averageTime(maxiter = maxiter)
            sumTime = 0.0
            sumStart = 0.0
            
            #myT = self.build.options._temperature_kelvin # ino dorost kon to nzzf
            myT = self.options._temperature_kelvin # ino dorost kon to nzzf
            RT = Energy.GAS_CONSTANT * myT
            
            for state in  builderRate.initial_states:
                
                dG = self.protoSpace[state].dG(myT)
                stateprob = np.e ** ( -dG / RT )
                
                stateindex = builderRate.stateIndex[state]
                sumTime +=    ( times[stateindex]  * stateprob  )
                sumStart +=  (   stateprob )
            
            averagedMFPT=  sumTime / sumStart
        else:
            
            builderRate.setMatrix()  #this lines should be here all will get error on self.stateIndex on next lines
        
        if printCount:
            beforeN = len(self.protoFinalStates)
        
        #Now add states that are delta close to the set of final states
        
        for state in builderRate.statespace:
            
            #if state not in self.protoFinalStates:
            if state not in self.protoFinalStates and state not in self.protoInitialStates:
                
                stateindex = builderRate.stateIndex[state]
                
                if times[stateindex] < delta * averagedMFPT:
                    
                    self.protoFinalStates[state] = Literals.success
        
        if printCount:
            print "Number of final states was %i but now is %i" % (beforeN, len(self.protoFinalStates))
        
        return times, averagedMFPT
    
    
    """
    supplyInitialState : a Complex that serves as initial state for that simulation
    ignoreIntiialState: The initial state is not added to the set of initial states
    """
    def genAndSavePathsFile(self, ignoreInitialState=False, supplyInitialState=None, inspecting=False):
        
        self.startTime = time.time()
        
        space = dict()
        transitions = dict()
        initStates = dict()
        finalStates = dict()
        sequences = dict()
        initialpathwayuniqueIDs   = dict( )
        # the first argument is always the number of paths
        inputArgs = copy.deepcopy(self.optionsArgs)
        
        def runPaths(optionsF, optionsArgs, space, transitions, initStates, finalStates, sequences):
            
            myOptions = optionsF(optionsArgs)
            myOptions.activestatespace = True
            myOptions.output_interval = 1
            
            if not supplyInitialState == None:
                myOptions.start_state = supplyInitialState
            
            """ Set longer searching time for the initial state. """
            if not ignoreInitialState:
                
                myOptions.simulation_time = myOptions.simulation_time * 10
            
            simTime = time.time()
            
            s = SimSystem(myOptions)
            
            ''' a specialized routine that ends after taking all transitions'''
            if inspecting == True:
                s.localTransitions()
            else:
                s.start()  # after this line, the computation is finished.
            
            #if self.verbosity:
            
            if self.fatten ==True :
                self.multistrand_simulation_fatten_time +=  ( time.time() - simTime )
            else:
                self.multistrand_simulation_time +=  ( time.time() - simTime )
            #print "Multistrand simulation is now done,      time = %.2f" % (time.time() - simTime)
            readingTime = time.time()
            """ load the space """
            count = 0
            if len(myOptions.start_state) == 3:
                initialstatestructure1 = [  myOptions.start_state[0].structure,myOptions.start_state[1].structure, myOptions.start_state[2].structure ]
                initialstatestructure2 = [  myOptions.start_state[0].structure,myOptions.start_state[2].structure, myOptions.start_state[1].structure ]
                initialstatestructure3 = [  myOptions.start_state[1].structure,myOptions.start_state[0].structure, myOptions.start_state[2].structure ]
                initialstatestructure4 = [  myOptions.start_state[1].structure,myOptions.start_state[2].structure, myOptions.start_state[0].structure ]
                initialstatestructure5 = [  myOptions.start_state[2].structure,myOptions.start_state[1].structure, myOptions.start_state[0].structure ]
                initialstatestructure6 = [  myOptions.start_state[2].structure,myOptions.start_state[0].structure, myOptions.start_state[1].structure ]
            
            elif len(myOptions.start_state) == 2:
                initialstatestructure1 = [  myOptions.start_state[0].structure,myOptions.start_state[1].structure ]
                initialstatestructure2 = [  myOptions.start_state[1].structure,myOptions.start_state[0].structure ]
                initialstatestructure3  = None
                initialstatestructure4  = None
                initialstatestructure5  = None
                initialstatestructure6  = None
            
            elif len(myOptions.start_state) == 1:
                initialstatestructure1  = [myOptions.start_state[0].structure]
                initialstatestructure2  = None
                initialstatestructure3  = None
                initialstatestructure4  = None
                initialstatestructure5  = None
                initialstatestructure6  = None
            else:
                raise ValueError("myOptions.start_state  is something you did not think of ")
            with  open(self.the_dir + str(myOptions.interface.current_seed) + "/protospace.txt", "r") as myFile:
                
                for line in myFile:
                    
                    uniqueID, energyvals, seqs = self.parseState(line, myOptions._temperature_kelvin, myOptions.join_concentration)
                    
                    if not uniqueID in sequences:
                        sequences[uniqueID] = seqs
                        count +=1
                        #print seqs
                    if not uniqueID in space:
                        
                        space[uniqueID] = energyvals
                        if seqs[2] == initialstatestructure1 or seqs[2] == initialstatestructure2 or seqs[2] == initialstatestructure3 or seqs[2] == initialstatestructure4 or seqs[2] == initialstatestructure5 or seqs[2] == initialstatestructure6:
                            initialpathwayuniqueIDs[uniqueID]  = energyvals
                    
                    elif not space[uniqueID] == energyvals:
                        
                        print "My hashmap contains " + str(uniqueID) + " with Energy " + str(space[uniqueID]) + " but found: " + str(energyvals)
                        print "Line = " + line
            
            #print "count",count
            if self.fatten ==True :
                
                self.loading_protospace_fatten_time +=  ( time.time() - readingTime )
            
            else :
                self.loading_protospace_time +=  ( time.time() - readingTime )
            #print "loading the space  is now done,      time = %.2f" % (time.time() - readingTime)
            readingTime = time.time()
            
            
            """ load the transitions """
            with open(self.the_dir + str(myOptions.interface.current_seed) + "/prototransitions.txt", "r") as myFile:
                
                index = 0
                go_on = True
                
                myLines = []
                
                for line in myFile:
                    myLines.append(line)
                
                while go_on:
                    
                    line1 = myLines[index];
                    line2 = myLines[index + 1];
                    line3 = myLines[index + 2];
                    
                    index = index + 4  # note the whitespace
                    
                    go_on = len(myLines) > index
                    
                    uID1, ev1, seq1 = self.parseState(line2, myOptions._temperature_kelvin, myOptions.join_concentration)
                    uID2, ev2, seq2 = self.parseState(line3, myOptions._temperature_kelvin, myOptions.join_concentration)
                    
                    transitionPair = (uID1, uID2)
                    
                    if not transitionPair in transitions:
                        
                        transitionList = list()
                        
                        n_complex1 = int(line2.split()[0])
                        n_complex2 = int(line3.split()[0])
                        
                        if n_complex1 == n_complex2:
                            transitionList.append(transitiontype.unimolecular)
                        
                        if n_complex1 > n_complex2:
                            transitionList.append(transitiontype.bimolecularIn)
                        
                        if n_complex2 > n_complex1:
                            transitionList.append(transitiontype.bimolecularOut)
                        
                        if myOptions.rate_method == Literals.arrhenius:
                            # decode the transition and add it
                            transitionList.extend(codeToDesc(int(float(line1))))
                        
                        transitions[transitionPair] = transitionList
            
            
            if self.fatten ==True :
                self.loading_prototransitions_fatten_time +=  ( time.time() - readingTime )
            else :
                self.loading_prototransitions_time +=  ( time.time() - readingTime )
                #print  time.time() - readingTime
                
                #print "reading transitions is now done,      time = %.2f" % (time.time() - readingTime)
            readingTime = time.time()
            
            """ load the initial states """
            with open(self.the_dir + str(myOptions.interface.current_seed) + "/protoinitialstates.txt", "r") as myFile:
                
                myLines = []
                
                for line in myFile:
                    myLines.append(line)
                
                index = 0
                go_on = True
                
                if len(myLines) == 0:
                    print "No initial states found!"
                
                while go_on:
                    
                    line1 = myLines[index];
                    line2 = myLines[index + 1];
                    
                    index = index + 2  # note the whitespace
                    go_on = len(myLines) > index
                    
                    uID1, ev1, seq1 = self.parseState(line2, myOptions._temperature_kelvin, myOptions.join_concentration)
                    count = int(line1.split()[0])
                    
                    if not uID1 in initStates:
                        
                        newEntry = InitCountFlux()
                        newEntry.count = count
                        newEntry.flux = 777777  # arrType is the flux, and is unique to the initial state
                        
                        initStates[uID1] = newEntry
            
            
            if self.fatten ==True :
                self.loading_protoinitialstates_fatten_time +=  ( time.time() - readingTime )
            else :
                self.loading_protoinitialstates_time +=  ( time.time() - readingTime )
            
            
            #print "reading initials is now done,      time = %.2f" % (time.time() - readingTime)
            readingTime = time.time()
            
            """ load the final states """
            with open(self.the_dir + str(myOptions.interface.current_seed) + "/protofinalstates.txt", "r") as myFile:
                
                myLines = []
                
                for line in myFile:
                    myLines.append(line)
                
                index = 0
                go_on = True
                
                if len(myLines) == 0:
                    #                 raise ValueError("No succesful final states found -- mean first passage time would be infinite ")
                    go_on = False
                
                while go_on:
                    
                    line1 = myLines[index];
                    line2 = myLines[index + 1];
                    index = index + 2
                    
                    go_on = len(myLines) > (index + 1)
                    
                    uID1, ev1, seq1 = self.parseState(line1, myOptions._temperature_kelvin, myOptions.join_concentration)
                    tag = line2.split()[0]
                    
                    if not uID1 in finalStates:
                        finalStates[uID1] = tag
            
            
            if self.fatten ==True :
                self.loading_protofinalstates_fatten_time +=  ( time.time() - readingTime )
            else :
                self.loading_protofinalstates_time +=  ( time.time() - readingTime )
            #print "reading final states is now done,      time = %.2f" % (time.time() - readingTime)
            
            """ Now delete the files as they can get quite large """
            
            os.remove(self.the_dir + str(myOptions.interface.current_seed) + "/protospace.txt")
            os.remove(self.the_dir + str(myOptions.interface.current_seed) + "/prototransitions.txt")
            os.remove(self.the_dir + str(myOptions.interface.current_seed) + "/protoinitialstates.txt")
            os.remove(self.the_dir + str(myOptions.interface.current_seed) + "/protofinalstates.txt")
            os.rmdir(self.the_dir + str(myOptions.interface.current_seed))
        
        runPaths(self.optionsFunction, inputArgs, space, transitions, initStates, finalStates, sequences)
        # do not forget to merge the objects back
        self.mergeSet(self.protoinitialpathwayuniqueIDs , initialpathwayuniqueIDs)
        self.mergeSet(self.protoSpace, space)
        self.mergeSet(self.protoTransitions, transitions)
        self.mergeSet(self.protoFinalStates, finalStates)
        self.mergeSet(self.protoSequences, sequences)
        
        if not ignoreInitialState:
            
            for key, val in initStates.iteritems():
                if not key in self.protoInitialStates:
                    self.protoInitialStates[key] = val;
                else:
                    self.protoInitialStates[key].count += val.count


# FD: This class is in progress.
# this class takes a builder object and computes the average time between
# starting in an initial state and reaching a final state.
class BuilderRate(object):
    
    solveToggle = 2000
    
    # input function returns the multistrand options object for which to build the statespace
    def __init__(self, builderIn, ignoremessage = False ):
        
        self.build = builderIn
        #self.rateLimit = 1e-5
        self.rateLimit = 0.
        
        
        if ignoremessage == False :
            #self.build.pathwayelaboration_use_elaboration=  True
            if self.build.pathwayelaboration_use_elaboration  ==  False :
                print "Removing all states and transitions  that were found by pathway elaboration "
                #print  len(  self.build.protoSpace ) , "len initial statses ",  len(  self.build.protoInitialStates ) , len(  self.build.protoFinalStates )
                for state in    list(self.build.protoSpace.keys()):
                    if state not in   self.build.protoinitialpathwayuniqueIDsbackup:
                        del self.build.protoSpace[state]
                
                for state in list(self.build.protoFinalStates.keys()):
                    if state not in   self.build.protoinitialpathwayuniqueIDsbackup:
                        del self.build.protoFinalStates[state]
                for state in list(self.build.protoInitialStates.keys()):
                    if state not in   self.build.protoinitialpathwayuniqueIDsbackup:
                        del self.build.protoInitialStates[state]
                
                #print len( self.build.protoSpace ), "len initial states ",  len(  self.build.protoInitialStates), len(  self.build.protoFinalStates )
            
            
            
            if len(self.build.protoFinalStates) == 0 :
                raise ValueError('No final states found.')
            
            self.processStates()  # prunes statespace and creates objects that can be used to create the rate matrx
            #self.setMatrix()  # generates the matrix for the current temperature
    
    
    
    """ Generates the state space by traversing from the final states """
    
    def findConnectedStates(self, final_states, neighborsdict):
        
        new_statespace = set()
        unexplored = set()
        
        for state in final_states:
            unexplored.add(state)
        
        while unexplored:  # empty dicts evaluate to false in python
            
            new_unexplored = set()
            
            for state in unexplored:
                
                new_statespace.add(state)
                
                # Make sure all neighbors are explored whenever we add a state.
                neighbors = neighborsdict[state]
                
                for s in neighbors :
                    
                    if s not in new_statespace :
                        
                        new_unexplored.add(s)
            
            unexplored = new_unexplored
        
        return new_statespace
    
    """ generates a list of neighbors for each state -- each transition goes both ways """
    
    def genNeighborsTransitive(self, transitions, statespace):
        
        output = dict()  # key: state, value: a set of states that are neighboring
        
        # each state starts out with no neighbors
        for state in statespace:
            output[state] = list()
        
        for transition in transitions:
            if self.build.pathwayelaboration_use_elaboration  == False :
                if transition[0] not in self.build.protoSpace  or transition[1]  not in  self.build.protoSpace :
                    continue
            output[transition[0]].append(transition[1])
            output[transition[1]].append(transition[0])
        
        return output
    
    """ generates a list of neighbors so that for each transition only one direction is added
        This is important when we build the matrix - so that we do not doubly add transitions """
    
    def genNeighbors(self, transitions, statespace):
        
        output = dict()  # key: state, value: a set of states that are neighboring
        
        # each state starts out with no neighbors
        for state in statespace:
            output[state] = list()
        
        for transition in transitions:
            if self.build.pathwayelaboration_use_elaboration  ==  False :
                if transition[0] not in self.build.protoSpace  or transition[1]  not in  self.build.protoSpace :
                    continue
            
            # ensure the state is connected and add only one direction of a transition
            if (transition[0] in statespace) and not (transition[0] in output[transition[1]]) :
                
                
                output[transition[0]].append(transition[1])
        
        return output
    
    def processStates(self):
        
        self.statespace = set()
        self.initial_states = dict()  # key: state, value: number of times started in this state
        self.final_states = set()
        
        # Storing the final states
        for state in self.build.protoFinalStates:
            if self.build.protoFinalStates[state] == Literals.success:
                self.final_states.add(state)
        
        if len(self.final_states) == 0:
            raise ValueError("No final states found!")
        
        startT = time.time()
        # generate a list of neighbors for each state  -- output is placed in self.neighbors
        self.neighbors = self.genNeighborsTransitive(self.build.protoTransitions, self.build.protoSpace)
        
        startT = time.time()
        # now prune the statespace to only include states that can reach the final state
        self.statespace = self.findConnectedStates(self.final_states, self.neighbors)
        
        startT = time.time()
        # Now re-generated the neighbors, but only for states in the statespace -- and only "forward" transitions instead of both ways.
        self.neighbors = self.genNeighbors(self.build.protoTransitions, self.build.protoSpace)
        
        # Update the initial states to only include those connected states
        for state in self.build.protoInitialStates:
            if state in self.statespace:
                self.initial_states[state] = self.build.protoInitialStates[state]
    
    """Use the builder options object to determine which rates to compute """
    
    def halfcontext_parameter(self, localContext):
        
        if localContext == localtype.stack:
            return 	self.build.options.lnAStack , self.build.options.EStack
        
        elif localContext == localtype.loop :
            return self.build.options.lnALoop, self.build.options.ELoop
        
        elif localContext == localtype.end:
            return self.build.options.lnAEnd, self.build.options.EEnd
        
        elif localContext == localtype.stackloop:
            return self.build.options.lnAStackLoop, self.build.options.EStackLoop
        
        elif localContext == localtype.stackend :
            return self.build.options.lnAStackEnd, self.build.options.EStackEnd
        
        elif localContext == localtype.loopend :
            return self.build.options.lnALoopEnd, self.build.options.ELoopEnd
        
        elif localContext == localtype.stackstack:
            return self.build.options.lnAStackStack, self.build.options.EStackStack
        
        else :
            raise ValueError('The transition code name is unexpected ')
    
    
    
    """    Returns the transition rate from state1 to state2 and then also the reverse rate, from state2 to state1 """
    def get_rate(self, state1, state2):
        
        transitionlist = self.build.protoTransitions[(state1, state2)]
        if self.build.options.rate_method == Literals.arrhenius:
            return self.arrhenius_rate(state1, state2, transitionlist)
        else :
            return self.metropolis_rate(state1, state2, transitionlist)
    
    """Uses the Metropolis kinetic model to calculate transition rates.     Returns the transition rate from state1 to state2  and then also the reverse rate, from state2 to state1 """
    def metropolis_rate(self, state1, state2, transitionlist):
        
        myT = self.build.options._temperature_kelvin
        RT = Energy.GAS_CONSTANT * myT
        
        dG1 = self.build.protoSpace[state1].dG(myT)
        dG2 = self.build.protoSpace[state2].dG(myT)
        
        if transitionlist[0] == transitiontype.unimolecular:
            
            if dG1 > dG2 :  # state2 is more stable (negative), dG1 - dG2 is positive
                
                rate1 = self.build.options.unimolecular_scaling
                rate2 = self.build.options.unimolecular_scaling * np.exp(-(dG1 - dG2) / RT)
            
            else:  # state2 is less or equally stable (negative), dG1 - dG2 is negative
                
                rate1 = self.build.options.unimolecular_scaling * np.exp((dG1 - dG2) / RT)
                rate2 = self.build.options.unimolecular_scaling
            
            return rate1, rate2
        
        else:  # bimolecular rate
            
            collisionRate = self.build.options.join_concentration * self.build.options.bimolecular_scaling
            
            if transitionlist[0] == transitiontype.bimolecularIn:
                
                outR = self.build.options.bimolecular_scaling * np.e ** (-(dG1 - dG2) / RT)
                #print "bimolecularIn", collisionRate, outR  , dG1, dG2
                return collisionRate, outR
            
            elif transitionlist[0] == transitiontype.bimolecularOut:
                
                outR = self.build.options.bimolecular_scaling * np.e ** ((dG1 - dG2) / RT)
                #print "bimolecularOut", collisionRate , outR,  dG1 , dG2
                return outR, collisionRate
            
            else:
                raise ValueError('The transition code is unexpected ')
    
    
    """     Uses the Arrhenius kinetic model to calculate transition rates.     Returns the transition rate from state1 to state2    and then also the reverse rate, from state2 to state1     """
    def arrhenius_rate(self, state1, state2, transitionlist):
        
        lnA_left, E_left = self.halfcontext_parameter(transitionlist[1])
        lnA_right, E_right = self.halfcontext_parameter(transitionlist[2])
        
        lnA = lnA_left + lnA_right
        E = E_left + E_right
        
        bimolecular_scaling = self.build.options.bimolecular_scaling
        concentration = self.build.options.join_concentration
        
        myT = self.build.options._temperature_kelvin
        RT = Energy.GAS_CONSTANT * myT
        dG1 = self.build.protoSpace[state1].dG(myT)
        dG2 = self.build.protoSpace[state2].dG(myT)
        
        DeltaG = dG2 - dG1
        DeltaG2 = -DeltaG
        
        if transitionlist[0] == transitiontype.unimolecular:
            if DeltaG > 0.0:
                rate1 = np.e ** (lnA - (DeltaG + E) / RT)
                rate2 = np.e ** (lnA - E / RT)
            else:
                rate1 = np.e ** (lnA - E / RT)
                rate2 = np.e ** (lnA - (DeltaG2 + E) / RT)
        
        
        
        elif transitionlist[0] == transitiontype.bimolecularIn:
            
            rate1 = bimolecular_scaling * concentration * (np.e ** (lnA - E / RT))
            rate2 = bimolecular_scaling * (np.e ** (lnA - (DeltaG2 + E) / RT))
        
        elif transitionlist[0] == transitiontype.bimolecularOut:
            
            rate1 = bimolecular_scaling * (np.e ** (lnA - (DeltaG + E) / RT))
            rate2 = bimolecular_scaling * concentration * (np.e ** (lnA - E / RT))
        else :
            raise ValueError('Exception in transition rate calculations.')
        
        return rate1, rate2
    
    """Set the rate matrix for this transition. If the target state is a final state, only subtract the outgoing rate from the diagonal. """
    
    def addTransition(self, state, neighbor, rate, rates, iArray, jArray, stateIndex):
        
        if not state in self.final_states:
            
            # now add the negative outgoing rate on the diagonal
            rates[stateIndex[state]] += -1.0 * rate
            
            if not neighbor in self.final_states:
                
                # set the transition
                rates.append(rate)
                iArray.append(stateIndex[state])
                jArray.append(stateIndex[neighbor])
    
    """ given the statespace, the initial states and final states, and the original builder object,
    build rate matrix for the current temperature. """
    
    def setMatrix(self):
        
        # give every state an explicit index
        self.stateIndex = dict()
        N = 0
        
        # the index works like this: rate[i] is in position iArray[i], jArray[j]
        # This is then interperted as an NxN matrix
        rates = list()
        iArray = list()
        jArray = list()
        
        
        for state in self.statespace:
            
            if not state in self.final_states:
                
                # give every state an explicit index
                self.stateIndex[state] = N
                
                # populate the diagonal so we can manipulate this as we go along
                rates.append(0.0)
                iArray.append(N)
                jArray.append(N)
                
                N += 1
        
        # post: N is the size of the statespace
        # post: rates[stateIndex[state]] is the diagonal for that state
        
        # first, set all transitions
        
        for state in self.statespace:
            
            for neighbor in self.neighbors[state]:
                
                myRate, revRate = self.get_rate(state, neighbor)
                
                if myRate < 0.0 or revRate < 0.0:
                    ValueError('Negative transition rate found.')
                
                # This handles either state being a final state (in which case, subtract from the non-final state diagonal,
                # but do not add the transiion rate.
                
                if myRate > self.rateLimit:
                    self.addTransition(state, neighbor, myRate, rates, iArray, jArray, self.stateIndex)
                if revRate > self.rateLimit:
                    self.addTransition(neighbor, state, revRate, rates, iArray, jArray, self.stateIndex)
        
        rate_matrix_coo = coo_matrix((rates, (iArray, jArray)), dtype=floatT)
        
        self.rate_matrix_csc = csc_matrix(rate_matrix_coo)
        
        #print "finished making ratematirx", self.build.outfile, len(rates ) ,min(rates) , max(rates)
        
        self.b = -1 * np.ones(N)
        self.b = np.array(self.b)
        #         # FD: pre-compute the matrix diagonal for preconditioning
        diagons = [ ((1.0 / x), i, j) for x, i, j in zip(rates, iArray, jArray) if i == j ]
        diagons0 = [x[0] for x in diagons]
        diagons1 = [x[1] for x in diagons]
        diagons2 = [x[2] for x in diagons]
        
        diagonal_matrix_coo = coo_matrix((diagons0, (diagons1, diagons2)), shape=(N, N), dtype=floatT)
        self.rate_matrix_inverse = csr_matrix(diagonal_matrix_coo, dtype=floatT)
        
        # save two counts for later interest
        self.n_states = N
        self.n_transitions = len(rates)
    
    """
        Computes the first passage times
    """
    
    def averageTime(self, x0=None, maxiter=None):
        
        st = time.time()
        self.setMatrix()
        startTime = time.time()
        
        if self.solveToggle == 1:
            firstpassagetimes, info = bicg(self.rate_matrix_csc, self.b, x0=x0, maxiter=maxiter,atol=1e-05)
        
        elif self.solveToggle == 2:
            
            firstpassagetimes, info = bicg(self.rate_matrix_csc, self.b, M=self.rate_matrix_inverse, x0=x0, maxiter=maxiter,atol=1e-05)
        
        elif self.solveToggle == 3:
            firstpassagetimes, info = bicgstab(self.rate_matrix_csc, self.b, x0=x0, maxiter=maxiter)
        
        elif self.solveToggle == 4:
            firstpassagetimes, info = bicgstab(self.rate_matrix_csc, self.b, M=self.rate_matrix_inverse, x0=x0, maxiter=maxiter)
        
        elif self.solveToggle == 5:
            firstpassagetimes, info = gmres(self.rate_matrix_csc, self.b, x0=x0, maxiter=maxiter)
        
        elif self.solveToggle == 6:
            firstpassagetimes, info = gmres(self.rate_matrix_csc, self.b, M=self.rate_matrix_inverse, x0=x0, maxiter=maxiter)
        
        elif self.solveToggle == 7:
            firstpassagetimes, info = cg(self.rate_matrix_csc, self.b, x0=x0, maxiter=maxiter)
        
        elif self.solveToggle == 8:
            firstpassagetimes, info = cg(self.rate_matrix_csc, self.b, M=self.rate_matrix_inverse, x0=x0, maxiter=maxiter)
        
        elif self.solveToggle == 9:
            firstpassagetimes, info = lgmres(self.rate_matrix_csc, self.b, x0=x0, maxiter=maxiter)
        
        elif self.solveToggle == 10:
            firstpassagetimes, info = lgmres(self.rate_matrix_csc, self.b, M=self.rate_matrix_inverse, x0=x0, maxiter=maxiter)
        
        else:
            
            try :
                st =time.time()
                firstpassagetimes = spsolve(self.rate_matrix_csc, self.b ,  use_umfpack = True )
            except RuntimeWarning as  w:
                s = str(w)
                print w, self.build.outfile
                if 'overflow' in s :
                    print "Overflow warning :( "
                    self.matrixTime = time.time() - startTime
                    return [np.inf for i in range(len(self.b)) ]
                if 'underflow' in s :
                    print "Underflow warning :( "
                    self.matrixTime = time.time() - startTime
                    return [np.inf for i in range(len(self.b)) ]
                print "in builder.py runtimewarning , dont know what happend"
                self.matrixTime = time.time() - startTime
                return [np.inf for i in range(len(self.b)) ]
            
            except Exception as w :
                startTime = time.time()
                print  w,  self.build.outfile
                self.matrixTime = time.time() - startTime
                print "Got a memory error solving with the director solver, spsolve. Now using an iterative solver."
                firstpassagetimes, info = bicg(self.rate_matrix_csc, self.b, M=self.rate_matrix_inverse, x0=x0, maxiter=maxiter,atol=1e-05)
            except :
                print  "Exception - Don't know what happened  "
        
        
        
        self.matrixTime = time.time() - startTime
        
        return firstpassagetimes
    
    def averageTimeFromInitial(self, bimolecular=False, printMeanTime=False):
        
        times = self.averageTime()
        
        if self.build.verbosity or printMeanTime:
            print "Solving matrix took %.2f s" % self.matrixTime
        
        mfpt = self.weightedPassageTime(times=times, bimolecular=bimolecular)
        
        return mfpt
    
    
    
    """Weight the solution vector by the initial probability distribution """
    def weightedPassageTime(self, times, bimolecular=False):
        
        sumTime = 0.0
        sumStart = 0.0
        #print "len initial state ", len(self.initial_states)
        if len(self.initial_states) == 0:
            raise ValueError("The number of initial states connected to a final state is zero.")
        
        myT = self.build.options._temperature_kelvin
        RT = Energy.GAS_CONSTANT * myT
        
        for state in self.initial_states:
            
            dG = self.build.protoSpace[state].dG(myT) # no need to worry about the volume term because all initial states are  the same complexes and the relative energy will be ok
            stateprob = np.e ** ( -dG / RT )
            
            stateindex = self.stateIndex[state]
            sumTime +=    ( times[stateindex]  * stateprob  )
            sumStart +=  (   stateprob )
        return sumTime / sumStart
    
    
    """This only works  if there is one initial state and one final state and only for threeway strand displacement"""
    """
    def saveStatisticforContourPlot(self, times , MFPT ):
        def convertKey(key, convertedkeys):
            if key in convertedkeys :
                return convertedkeys[key ]
            try :
                newkey = literal_eval(key)
                newkey = (newkey, )
            except :
                newkey = ""
                for i in range (len(key)-1):
                    if key[i] ==")" and key[i+1] == "(":
                        newkey += key[i] +","
                    else:
                        newkey += key[i]
                newkey+= key[-1]
                newkey = literal_eval(newkey)
            convertedkeys [key ] = newkey
            return newkey
        convertedkeys = dict()
        print  " in save statstics for contour plot "
        contour_dict = dict()
        contour_dict["distancefrominitial"] = []
        contour_dict["distancefromfinal"] = []
        contour_dict["energy"] = []
        contour_dict["mfpt"] = []
        contour_dict["basemfpt"] = MFPT
        contour_dict["attackerbasepair"] = []
        contour_dict["incumbentbasepair"] = []
        
        #print "len of self.initial_states" , len(self.initial_states)
        
        for state1 in self.initial_states:
            (seqs_a, ids_a, structs_a) = self.build.protoSequences[state1]
            myT = self.build.options._temperature_kelvin
            #print "initial state$",  seqs_a, structs_a   , "\n"
            newinitialkey = convertKey(state1, convertedkeys)
            
         
        for state1 in self.final_states:
            (seqs_a, ids_a, structs_a) = self.build.protoSequences[state1]
            myT = self.build.options._temperature_kelvin
            newfinalkey = convertKey(state1 , convertedkeys)
            #print "final state$",  seqs_a, structs_a   , "\n"
            
           
        doaverage = False
        
        for state in  self.statespace :
            newkey  = convertKey(state, convertedkeys)
            distanceinitial = 0
            distancefinal=  0
            
            distanceinitial  += self.build.allowed_state( newkey , "" ,    "", "",  newinitialkey , "" ,  timetochooseunaligned = False , reportdistance= True   )
            distancefinal +=self.build.allowed_state( newkey , "" ,    "", "",  newfinalkey , "" ,  timetochooseunaligned = False , reportdistance= True  )
            
            attackerbasepair, incumbentbasepair = self.build.allowed_state(newkey  ,   "" ,    "", "",  "" , "" ,  "", "" ,distancebasedonbasepairswithsubstrate = True)
            if doaverage == True :
                
                distanceinitial = distanceinitial / len(self.initial_states)
                distancefinal = distancefinal / len(self.final_states)
            
            myT = self.build.options._temperature_kelvin
            dG =self.build.protoSpace[state].dG(myT)
            GAS_CONSTANT = floatT(0.0019872036)
            RT =  GAS_CONSTANT  * myT
            
            n_strands= 3
            concentration = 5 * 10 ** (  -9 )
            n_complexes =    len ( newkey )
            print "n_complexes " , n_complexes , newkey
            dG_volume = RT * (n_strands - n_complexes) * np.log(1.0 / concentration)
            #dG_new  = dG+  2 * dG_volume
            dG_new  = dG+   dG_volume
            
            (seqs_a, ids_a, structs_a) = self.build.protoSequences[state]
            
            if state not in self.final_states:
                stateindex = self.stateIndex[state]
                time = times [stateindex]
            else:
                time = 0
            
            print "\n$state$" , ", x: ", distanceinitial / 2 , ", y:",distancefinal /2 , ", delta: ", time /  MFPT , ", dG: " , dG , "dG_new: " , dG_new,     "  attackerbasepair: ", attackerbasepair , "   incumbentbasepair: " ,incumbentbasepair  , seqs_a, structs_a   , "\n"
            
            contour_dict["distancefrominitial"].append(distanceinitial)
            contour_dict["distancefromfinal"].append(distancefinal)
            contour_dict["energy"].append(dG_new)
            contour_dict["mfpt"].append(time )
            contour_dict["attackerbasepair"].append(attackerbasepair)
            contour_dict["incumbentbasepair"].append(incumbentbasepair)
        
        myT = self.build.options._temperature_kelvin
        RT = Energy.GAS_CONSTANT * myT
        contour_dict ["RT"]  = RT
        contour_dict ["R"] =  Energy.GAS_CONSTANT
        contour_dict["T"] = myT
        pathtodict="contourdictionary.pkl"
        with open(pathtodict, "wb") as p :
            pickle.dump(contour_dict, p)
    """
    
    def __str__(self):
        
        output = "This is a builder-rate object. Printing all states and transitions \n"
        
        for state in self.statespace:
            output += str(state) + " \n"
        
        output += "Transitions \n"
        
        # first, set all transitions
        for state in self.statespace:
            
            for neighbor in self.neighbors[state]:
                
                myRate, revRate = self.get_rate(state, neighbor)
                
                transitionlist = self.build.protoTransitions[(state, neighbor)]
                
                output += "s1 = " + str(state) + "  s2 = " + str(neighbor) + " for = " + "%.2E" % myRate + " back = " + "%.2E" % revRate + "   tlist = " + str(transitionlist) + "\n"
        
        output += str(self.rate_matrix_csc.toarray())
        
        return output
    
    def numOfTransitions(self):
        
        count = 0
        
        # first, count all transitions
        for state in self.statespace:
            for neighbor in self.neighbors[state]:
                
                transitionlist = self.build.protoTransitions[(state, neighbor)]
                count += len(transitionlist)
        
        return count
    
    def statsInfo(self):
        
        output = str(len(self.statespace)) + "   -    " + str(self.numOfTransitions())
        output += "   -    " + str(len(self.initial_states)) + "   -    " + str(len(self.final_states)) + "   -   (builderRate Object)"
        
        return output
