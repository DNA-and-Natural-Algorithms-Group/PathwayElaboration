import myenums
import learndnakinetics

TEST = 0 

"""Cisse, Ibrahim I., Hajin Kim, and Taekjip Ha. "A rule of seven in Watson-Crick base-pairing of mismatched sequences." Nature structural & molecular biology 19.6 (2012): 623."""
def read_Cisse2012(done_queue, dataset_list,n,  directories):
	print "Read data from read_Cisse2012."
	listt = zip([ False, False ,False,  False ,False, False ],   ["Fig1_k_",  "Fig2b_Fig3c_",  "sup_Fig4c_", "sup_Fig4g_" ,  "sup_Fig4m_", "sup_Fig4j_"] ) #using subset of Cisse!
	#listt = zip([True ,False, True, False , True ,False, True, False , True ,False, True, False ],   ["Fig1_i_", "Fig1_k_", "Fig2c_Fig3d_", "Fig2b_Fig3c_", "sup_Fig4d_", "sup_Fig4c_","sup_Fig4h_", "sup_Fig4g_" , "sup_Fig4n_", "sup_Fig4m_", "sup_Fig4k_", "sup_Fig4j_"] )
	for association, j  in listt :
		bimolecular_reaction = association
		if association == True :
			reaction_type = myenums.ReactionType.HELIXASSOCIATION.value
		else :
			reaction_type = myenums.ReactionType.HELIXDISSOCIATION.value
		dataset_type = myenums.DatasetType.MISMATCH.value
		my_name = '/helix4_cisse2012/'+ j  + str(int(association))
		dataset_path, document , row = learndnakinetics.initconf(my_name , directories )
		temperature_type =myenums.DatasetSpecifications.TEMPCELCIUS.value
		rate_type = myenums.DatasetSpecifications.RATECONSTANT.value
		dataset_name= 	myenums.DatasetName.CISSE2012.value
		learndnakinetics.check_directories (directories)
		docID= dataset_name+ str(bimolecular_reaction)
		dataset= learndnakinetics.Reaction (temperature_type=temperature_type, rate_type=rate_type, dataset_name =dataset_name, docID= docID, document=document,  reaction_type=reaction_type, dataset_type=dataset_type, bimolecular_reaction=bimolecular_reaction, row=row, dataset_path=dataset_path)
		(dataset_list, done_queue,n)  =learndnakinetics.objective_function_auxilary( done_queue, dataset_list,n,dataset)
	return (dataset_list, done_queue,n)



"""Hata, Hiroaki, Tetsuro Kitajima, and Akira Suyama. "Influence of thermodynamically unfavorable secondary structures on DNA hybridization kinetics." Nucleic Acids Research (2017)."""
def read_Hata2017(done_queue, dataset_list,n,  directories):
	print "Read data from read_Hata2017."
	for association in [True   ]:
		bimolecular_reaction = association
		if association == True :
			reaction_type = myenums.ReactionType.HELIXASSOCIATION.value
			dataset_type = myenums.DatasetType.NODANGLE.value
		else :
			reaction_type = myenums.ReactionType.HELIXDISSOCIATION.value
			dataset_type = myenums.DatasetType.NODANGLE.value
		my_name = '/helix2_hata2017/Table' +str(int(association ))
		dataset_path, document , row = learndnakinetics.initconf(my_name , directories  )
		temperature_type = myenums.DatasetSpecifications.TEMPCELCIUS.value
		rate_type = myenums.DatasetSpecifications.RATECONSTANT10POW5.value
		dataset_name= 	myenums.DatasetName.SUYAMA.value
		learndnakinetics.check_directories (directories)
		docID= dataset_name+ str(bimolecular_reaction)
		if TEST == 0:
			counter_celllist =[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 42, 43, 47]
		elif TEST==1 :
			counter_celllist =[41,44, 45, 46]
		dataset= learndnakinetics.Reaction (counter_celllist=counter_celllist, temperature_type=temperature_type, rate_type=rate_type, dataset_name =dataset_name, docID= docID, document=document,  reaction_type=reaction_type, dataset_type=dataset_type, bimolecular_reaction=bimolecular_reaction, row=row, dataset_path=dataset_path)
		(dataset_list, done_queue,n)  =learndnakinetics.objective_function_auxilary( done_queue, dataset_list,n,dataset)
	return (dataset_list, done_queue,n)



"""Bonnet, Gragoire, Oleg Krichevsky, and Albert Libchaber. "Kinetics of conformational fluctuations in DNA hairpin-loops." Proceedings of the National Academy of Sciences 95.15 (1998): 8602-8606."""
def read_Bonnet1998(done_queue, dataset_list,n,  directories)   :
	print "Read data from read_Bonnet1998."
	bimolecular_reaction = False
	for hairpinclosing in [ False , True ]:
		for j in [4,6] :
			my_name = '/hairpin_bonnet1998/Fig'+str(j)  + '_' + str(int(hairpinclosing))
			dataset_path, document , row = learndnakinetics.initconf(my_name , directories)
			learndnakinetics.check_directories (directories)
			dataset_type = myenums.DatasetType.NODANGLE.value
			if  hairpinclosing == True :
				reaction_type = myenums.ReactionType.HAIRPINCLOSING.value
			else :
				reaction_type = myenums.ReactionType.HAIRPINOPENING.value
			dataset_name = myenums.DatasetName.BONNET.value + str(hairpinclosing)
			docID = dataset_name+ str(j) + str(hairpinclosing)
			temperature_type =myenums.DatasetSpecifications.TEMPKELVININV.value
			rate_type =  myenums.DatasetSpecifications.RATECONSTANT.value
			dataset= learndnakinetics.Reaction (temperature_type=temperature_type, rate_type=rate_type, dataset_name =dataset_name, docID= docID, document=document,  reaction_type=reaction_type, dataset_type=dataset_type, bimolecular_reaction=bimolecular_reaction, row=row, dataset_path=dataset_path)
			(dataset_list, done_queue,n)  =learndnakinetics.objective_function_auxilary( done_queue, dataset_list,n,dataset)
	return (dataset_list, done_queue,n)


"""Machinek, R.R., Ouldridge, T.E., Haley, N.E., Bath, J., TurberField, A.J.: Programmable energy landscapes for kinetic control of DNA strand displacement. Nature Communications 5 (2014)"""
def read_Machinek2014(done_queue, dataset_list,n,  directories):
	print "Read data from read_Machinek2014."
	my_name = '/threeway_stranddisplacement2_machinek2014/Fig2'
	if TEST == 0:
		counter_celllist = [25,2,5,9, 10,   11,12,13,24,1]
	elif TEST==1 :
		counter_celllist =[3,4,6, 7, 8, 14, 15, 16, 17, 18 , 19,  20, 21, 22, 23,26, 27, 28, 29, 30,31,  32, 33, 34, 35, 36]
	dataset_path, document , row = learndnakinetics.initconf(my_name , directories)
	learndnakinetics.check_directories (directories)
	dataset_name =myenums.DatasetName.MACHINEK.value
	docID= dataset_name
	bimolecular_reaction = True
	reaction_type = myenums.ReactionType.THREEWAYDISPLACEMENT.value
	dataset_type =   myenums.DatasetType.MISMATCH.value
	rate_type =    myenums.DatasetSpecifications.RATECONSTANT.value
	temperature_type = myenums.DatasetSpecifications.TEMPCELCIUS.value
	dataset= learndnakinetics.Reaction (counter_celllist=counter_celllist, temperature_type=temperature_type, rate_type=rate_type, dataset_name =dataset_name, docID= docID, document=document,  reaction_type=reaction_type, dataset_type=dataset_type, bimolecular_reaction=bimolecular_reaction, row=row, dataset_path=dataset_path)
	(dataset_list, done_queue,n)  =learndnakinetics.objective_function_auxilary( done_queue, dataset_list,n,dataset)
	return (dataset_list, done_queue,n)



"""Zhang, Jinny X., et al. "Predicting DNA Hybridization Kinetics from Sequence." bioRxiv (2017): 149427."""
def read_Zhang2018_hybridization (done_queue, dataset_list,n,  directories):
	print "Read data from read_Zhang2018."
	for association in [True     ]:
		bimolecular_reaction = association
		reaction_type = myenums.ReactionType.HELIXASSOCIATION.value
		dataset_type = myenums.DatasetType.NODANGLE.value
		my_name = '/helix3_zhang2018/Table' +str(int(association ))
		dataset_path, document , row = learndnakinetics.initconf(my_name , directories  )
		temperature_type = myenums.DatasetSpecifications.TEMPCELCIUS.value
		rate_type =  myenums.DatasetSpecifications.LOG10RATECONSTANT.value
		dataset_name= 	myenums.DatasetName.ZHANG_hybridization.value
		learndnakinetics.check_directories (directories)
		docID= dataset_name+ str(bimolecular_reaction)
		counter_celllist= [1, 3, 22, 34, 40, 65, 75, 86, 98, 110, 115, 119, 123, 140, 150, 162, 170, 179, 200, 209]
		dataset= learndnakinetics.Reaction (counter_celllist=counter_celllist, temperature_type=temperature_type, rate_type=rate_type, dataset_name =dataset_name, docID= docID, document=document,  reaction_type=reaction_type, dataset_type=dataset_type, bimolecular_reaction=bimolecular_reaction, row=row, dataset_path=dataset_path)
		(dataset_list, done_queue,n)  =learndnakinetics.objective_function_auxilary( done_queue, dataset_list,n,dataset)
	return (dataset_list, done_queue,n)