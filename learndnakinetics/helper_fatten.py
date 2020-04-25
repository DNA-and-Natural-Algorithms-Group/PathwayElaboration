import sys
import cPickle as pickle
import os
import gc
sys.path.insert(0,os.path.realpath('../reactions'))

from multistrand.builder import Builder
import parent
from parent import *

import time

start = sys.argv[1]
end = sys.argv[2]
pathbuilder = sys.argv[3]
pathspace= sys.argv[4 ]
pathsequences = sys.argv[5]
pathoptions= sys.argv[6]
pathprotoinitialpathwayuniqueIDs = sys.argv[7]
onlyfatteninitialpathway = sys.argv[8]

#mytime = open ("times.txt", "a")
#mytime.write( pathbuilder + "   start " + str( start) + "end "  + str(end) + "\n" )
st = time.time()

with open(pathoptions  , "rb" ) as p:
	optionsArg  = pickle.load( p)

myBuilder=   Builder(parent.doReaction,   optionsArg )

with open(pathprotoinitialpathwayuniqueIDs  , "rb" ) as p:
	myBuilder.protoinitialpathwayuniqueIDsbackup  = pickle.load( p)

with open(pathspace  , "rb" ) as p:
	myBuilder.protoSpacebackup  = pickle.load( p)

with open( pathsequences  , "rb" ) as p :
	myBuilder.protoSequences  = pickle.load( p)

#mytime.write( "load time " + str( time.time()  - st )+"\n")
st  = time.time ( )
myBuilder.fattenStateSpace(start = int(start)  , end= int( end), onlyfatteninitialpathway = int( onlyfatteninitialpathway)  )
#mytime.write( "fatten time time " + str( time.time()   - st ) +"\n")
st = time.time()
with open( pathbuilder + "pt" + str(start)+"-"  +str(end) , "wb" ) as p:
	pickle.dump(myBuilder.protoTransitions,  p)
with open( pathbuilder + "pt" + str(start)+"-"  +str(end)  + "-fattentimes", "wb" ) as p:
	p.write(str( myBuilder.multistrand_simulation_fatten_time ) + " " +str(myBuilder.loading_protospace_fatten_time)+ " "+str( myBuilder.loading_prototransitions_fatten_time ) + " "+str( myBuilder.loading_protoinitialstates_fatten_time ) + " "+str( myBuilder.loading_protofinalstates_fatten_time ) )

#mytime.write( "save time " + str( time.time()  - st ) +"\n")
#mytime.close()

del myBuilder
#for i in range(2):
#	n = gc.collect()
