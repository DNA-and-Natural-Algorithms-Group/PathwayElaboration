import cPickle as pickle
import sys
import os
sys.path.insert(0,os.path.realpath('reactions'))
import parent
from parent import *



mymethod =  int( sys.argv[1] )
path= ["","","","","","","","",""]
path[0]=  sys.argv[2]
path[1]=  sys.argv[3]
path[2]=  sys.argv[4]
path[3]=sys.argv[5]
path[4]= sys.argv[6]
path[5]=  sys.argv[7]
path[6]= sys.argv[8]

builderpickletemp = sys.argv[9]


with open(builderpickletemp   , "rb") as p :
    builder=  pickle.load(p)

with open(path[0], "rb") as p :
    initialstate=   pickle.load(p)
with open(path[1], "rb") as p :
    pathwaystates=   pickle.load(p)
with open(path[2], "rb") as p :
    myhash=   pickle.load(p)
with open(path[3], "rb") as p :
    finalstate=   pickle.load(p)
with open(path[4], "rb") as p :
    keyoffinalstate =   pickle.load(p)
with open(path[5], "rb") as p :
    structureoffinalstate =   pickle.load(p)
with open(path[6], "rb") as p :
    extrastates =   pickle.load(p)
    
builder.bfs(  mymethod , initialstate,  pathwaystates,myhash  ,  finalstate, keyoffinalstate, structureoffinalstate , extrastates )

with open(path[1], "wb") as p :
    pickle.dump( pathwaystates, p)
with open(path[2], "wb") as p :
    pickle.dump( myhash,p)
with open(path[6], "wb") as p :
    pickle.dump( extrastates ,p)