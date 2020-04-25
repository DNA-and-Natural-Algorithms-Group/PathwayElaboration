from enum import Enum

"""used for saving computations to save truncated CTMCs and results"""
class Permanent_Folder(Enum):
	MYBUILDER = "MyBuilder"

""" Dataset name"""
class DatasetName(Enum):
	ALTANBONNET= "bubble_altanbonnet2003_"
	MORRISON = "helix_morrison2003_"
	REYNALDODISSOCIATE = "helix1_reynaldo2000_"
	REYNALDOSEQUENTIAL= "threeway_stranddisplacement1_reynaldo2000_"
	ZHANG= "threeway_stranddisplacement_zhang2009_"
	MACHINEK= "threeway_stranddisplacement2_machinek2014_"
	BONNET= "hairpin_bonnet1998_"
	GODDARD= "hairpin1_goddard2000_"
	KIMTABLE1 = "hairpin4_kim2006_Table1_"
	KIMFIG5 ="hairpin4_kim2006_Fig5_"
	DABBY= "fourway_strandexchange_dabby2013_"
	ZHANG_hybridization= "helix3_zhang2018_"
	ZHANG_stranddisplacement= "threeway_stranddisplacement3_zhangnotpublishedyet_"
	SUYAMA = "helix2_hata2017_"
	BROADWATER2016 = "threeway_stranddisplacement5_broadwater2016_"
	SHERRYCHEN2016= "threeway_stranddisplacement4_sherrychen2016_"
	GAO2006=  "helix5_gao2006_"
	DUPIUS2013= "helix7_dupius2013_"
	RAUZAN2013 = "helix8_rauzan2013_"
	CISSE2012  = "helix4_cisse2012_"
	WALLACE2001 =  "hairpin3_wallace2001_"
	AALBERTS2003  = "hairpin5_aalberts2003_"
	GROVES2015="fourway_strandexchange_groves2015_"

"""Type of reaction"""
class ReactionType(Enum):
	HELIXASSOCIATION = "HelixAssociation"
	HELIXDISSOCIATION ="HelixDissociation"
	HAIRPINCLOSING = "HairpinClosing"
	HAIRPINOPENING = "HairpinOpening"
	THREEWAYDISPLACEMENT = "ThreeWayDisplacement"
	FOURWAYBRANCHMIGRATION = "FourWayBranchMigration"
	BUBBLECLOSING = "BubbleClosing"
	
class DatasetType(Enum):
	NODANGLE  = "NoDangle"
	DANGLETWOSUBSTRATE = "DangleTwoSubstrate"
	INCUMBENTDANGLE= "IncumbentDangle"
	MISMATCH= "Mismatch"

"""Type of how temperature, reaction rate constant, timescale  is presented in database"""
class 	DatasetSpecifications(Enum) :
	TEMPKELVININV = "TempKevlinInverse"
	TEMPCELCIUS = "TempCelcius"
	TEMPKELVIN  = "TempKelvin"
	RATECONSTANT= "RateConstant"
	LOG10RATECONSTANT= "Log10RateConstant"
	RATECONSTANT10POW5 = "RateConstant10Pow5"
	TIMESCALE = "TimeScale"

