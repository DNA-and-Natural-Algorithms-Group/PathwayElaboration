from __future__ import division
from scipy.optimize import minimize
import learndnakinetics
import timeit
from multistrand.options import Literals

def initialize_x0 () :
	if learndnakinetics.rate_method == Literals.arrhenius :
		if learndnakinetics.use_regularizer == True:
			theta = [  13.0580, 3.,  13.0580, 3.,   13.0580, 3.,  13.0580 , 3.,   13.0580, 3.,  13.0580,  3.,   13.0580 , 3.,    0.0402,1  ]
		else:
			theta=  [  13.0580, 3.,  13.0580, 3.,   13.0580, 3.,  13.0580 , 3.,   13.0580, 3.,  13.0580,  3.,   13.0580 , 3.,    0.0402  ]
	elif learndnakinetics.rate_method == Literals.metropolis:
		#theta = [8.2 *  (10 **6), 3.3  * (10**5) ,1 ]
		if learndnakinetics.use_regularizer == True:
			#theta =  [5 * 10000, 5 * 10000,1]
			theta =  [ 2.41686715e+06,   8.01171383e+05,1]
		else:
			#theta= [5 * 10000, 5 * 10000]
			theta =  [ 2.41686715e+06,   8.01171383e+05]
	else:
		raise ValueError('Error: Please specify rate_method to be Arrhenius or Metropolis!')

	return theta

def map ( ):
	""" This function runs the MAP approach with  the Nelder-Mead optimization technique!"""
	global use_Gillespie_MFPT
	learndnakinetics.OVERALLTIME= timeit.default_timer()

	#Initializing the parameters for Nelder-Mead
	options = dict()
	options['disp'] = True
	options['maxfev'] = 200
	options['maxiter']  =200
	learndnakinetics.set_configuration() #setting some configuations
	use_ensemble=False
	if learndnakinetics.do_inference == True :
		numberofiter= 20
	else:
		numberofiter= 1
	for learndnakinetics.map_i in range (0,numberofiter):
		#if learndnakinetics.map_i  == 1 :
		#if learndnakinetics.map_i  == 0 :
		#	learndnakinetics.iter= 1
		learndnakinetics.set_folderfiles()
		if learndnakinetics.map_i ==  0 :
			if use_ensemble == False :
				theta = initialize_x0()
		learndnakinetics.objective_function(theta)
		
		#optimizing the paramters
		if learndnakinetics.do_inference == True :
			thetaOptimized = minimize(learndnakinetics.objective_function, theta, method='nelder-mead',options=options )
			theta= thetaOptimized.x
		
		learndnakinetics.iter= 0
	
if __name__ == "__main__":
	map()
