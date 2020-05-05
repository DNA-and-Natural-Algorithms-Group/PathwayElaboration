import cPickle as pickle
import sys
import os
from multistrand.builder import Builder
sys.path.insert(0,os.path.realpath('reactions'))

import parent
from parent import *

builderpath = sys.argv[1]
convergence_crit= int( sys.argv[2])
printMeanTime= bool(int( sys.argv[3] ))
fatten_regardless= bool (int( sys.argv[4] ))
ID=  sys.argv[5]
reactionargspath  = sys.argv[6]
startstatespath = sys.argv[7]
builderpicklepath =  sys.argv[8]
breaklooppath= sys.argv[9]
len_protoSpacepath = sys.argv[10]
pathwayelaboration_N = int( sys.argv[11 ] )
pathwayelaboration_beta=  float(sys.argv[12])# if you put self.pathwayelaboration_beta = 0 , then misaligned base pairs will not be allowed in initial random pathways
uniqueIDtoLength_path = sys.argv[13]
auto_strand_id_list_path = sys.argv[14]
boltzmann_path = sys.argv[15]
pathwayelaboration_use_elaboration = bool(int( sys.argv[16] ))

with open(reactionargspath, "rb") as p :
	reactionargs= pickle.load(p )

with open(startstatespath, "rb") as p :
	startStates= pickle.load(p )

with open(uniqueIDtoLength_path, "rb") as p :
	uniqueIDtoLength= pickle.load(p )


with open(auto_strand_id_list_path, "rb") as p :
	auto_strand_id_list= pickle.load(p )

with open (boltzmann_path, "rb") as p:
	boltzmann_sample_initialstates_list = pickle.load(p)

myBuilder = Builder(parent.doReaction, reactionargs)
myBuilder.outfile = builderpath


if parent.build_truncatedCTMC == True and parent.build_truncatedCTMC_GillespieSSA == False and parent.build_truncatedCTMC_pathwayelaboration == True:
	myBuilder.generateInitialPathwayAutomatically( startStates[0], startStates[-1],  printMeanTime = False , pathwayelaboration_N=pathwayelaboration_N , pathwayelaboration_beta = pathwayelaboration_beta, uniqueIDtoLength=uniqueIDtoLength , auto_strand_id_list = auto_strand_id_list , boltzmann_sample_initialstates_list  = boltzmann_sample_initialstates_list  , pathwayelaboration_use_elaboration = pathwayelaboration_use_elaboration)
	breakloop = True
	len_protoSpace = len(myBuilder.protoSpace)
#elif parent.build_truncatedCTMC == True and parent.build_truncatedCTMC_GillespieSSA == False and parent.build_truncatedCTMC_pathwayelaboration == False   :
#	breakloop , len_protoSpace = myBuilder.genUntilConvergenceWithInitialState_doublesim(convergence_crit,  startStates[:(len(startStates) )], printMeanTime=printMeanTime, fatten_regardless = fatten_regardless)
elif parent.build_truncatedCTMC == True and parent.build_truncatedCTMC_GillespieSSA== True :
	
	myBuilder.genAndSavePathsFile()
	len_protoSpace= len(myBuilder.protoSpace)
	breakloop = True


with open(builderpicklepath, "wb") as p :
	pickle.dump(myBuilder, p)

with open(breaklooppath, "wb") as p :
	pickle.dump(breakloop, p)

with open(len_protoSpacepath, "wb") as p :
	pickle.dump(len_protoSpace, p)
