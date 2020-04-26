# Introduction

This software implements  the **pathway elaboration** method for (multiple) interacting nucleic acid strands using the  Multistrand simulator [2]. It could be used for estimating the mean first passage time estimation (MFPT) (and reaction rate constant) of rare events and the rapid evaluation of perturbed parameters. 
 
 This software uses and modifies code from the <a href="https://github.com/DNA-and-Natural-Algorithms-Group/multistrand">Multistrand repository</a> [2]. 

For more information regarding the methods employed in our software, see our paper [1] (under review). 

Please contact us at nasimzf@cs.ubc.ca for questions regarding the software.

Disclaimer: This software is free for use only for research purposes. If you make use of this software, please cite our paper

[1]  Sedigheh Zolaktaf*, Frits Dannenberg*, Mark Schmidt, Anne Condon, and Erik Winfree, The Pathway Elaboration Method for Mean First Passage Time Estimation in Large Continuous-Time Markov Chains with Applications to Nucleic Acid Kinetics, under review.

# Requirements 

| Dependency       | Notes          
| ------------- |:-------------:|
|Python 2 |2.7.14 |
|c++11|gcc 4.8.5+|
|make|4.0+ | 
|<a href="http://www.nupack.org/">NUPACK</a> |3.0.4+| 
|Multistrand|We have modified <a href="https://github.com/DNA-and-Natural-Algorithms-Group/multistrand">Multistrand</a>.  For this project,     <br />  use  the   modified version  provided in this directory| 
|numpy|1.16.2+| 
|SciPy | 1.2.0+ |
|scikits.umfpack|0.3.2+|
| enum34 | 1.1.6+ | 
| ConfigParser | 3.5.0+ | 
|matplotlib| 2.0.2+|



# Setup
For linux:
- Download and install the requirements. 
- Set the system environment variables to point NUPACKHOME to the folder where NUPACK is installed. For example, in openSUSE, you can run : ' export NUPACKHOME=path/to/NUPACK3.0.4 '. To verify that the NUPACKHOME is correctly set, run ' echo $NUPACKHOME '. 
- Clone our software directory into your workspace. 
- In the  multistrand_modified directory,
  - Build Multistrand by running ' make clean ' and then ' make '. Upon successful  building of  Multistrand you should get a message  "Multistrand is now built. ... ".
  - Export  Multistrand as a Python library by running  ' make install '. 
  - Set the  PYTHONPATH evironment variable to point to the multistrand_modified  directory. For example, in openSUSE, you can use the following command: ' export PYTHONPATH="${PYTHONPATH}:path/to/multistrand_modified" '.  

Note that this software only works with the modified Multistrand code provided here, and will not work with the original version of Multistrand.

For efficiently solving system of linear equations: 
- Don't forget to install scikits.umfpack.  
- And turn of multithreading for the matrix solvers as follows: 
  - export MKL_NUM_THREADS=1
  - export NUMEXPR_NUM_THREADS=1
  - export OMP_NUM_THREADS=1

 

# Package Tree
This software contains 4  directories, namely, learndnakinetics, reactions, dataset, and multistrand_modified. 
- learndnakinetics and reactions:  these directories contain code to run the pathway elaboration method. 
- dataset: this directory contains reaction rate constants and timescales that we compiled from published literature
  - hairpin_bonnet1998:  These reactions are originally from [3]. 
  - helix4_cisse2012: These reactions are originally from  [4].
  - helix2_hata2017: These reactions are originally from  [5]. 
  - threeway_stranddisplacement2_machinek2014: Threse reaction are originally from [6].
- multistrand_modfied: The modified Multistrand code adapted from the branch 'fdann-devel' from  <a href="https://github.com/DNA-and-Natural-Algorithms-Group/multistrand">  Multistrand repository</a> [2]. 


# Software Usage 
To use pathway elaboration with the dataset provided:
- In ' learndnakinetics/config_file.txt ': 
- Set *n_processors* to be the number of processors for multiprocessing the computation of the objective function. Set *use_multiprocess =1* to multiprocesses the computations (each reaction will use a distinct processor). Set *use_multiprocess=0* to turnoff multiprocessing. 
  - Set *rate_method=1* to use the Metropolis kinetic model. Set *rate_method=3* to use the Arrhenius kinetic model. 
  - Set *parameter_folder* to be the path to a directory to save results. 
  - Set  *do_inference=1* to run parameter estimation. Set *do_inference=0* to turn off parameter estimation. 
  - Set  *use_regularizer=1* to use a  regularizer (for parameter estimation). 
  - Set *filter_smallandlarge_rates=1* to filter parameters that are really slow or high.  (See code in learndnakinetics.py for values)
  - Set *load_existing_override=1* to reuse truncated CTMCs (if it exists) instead of building from scratch. 
  - Set *deltaPruning=1* to do delta-pruning.  (See code in parent.py and learndnaknietics.py to change functionality)
  - Set *pathwayelaboration_N*  corresponds to N in the pathway elaboration method [1]. 
  - Set *pathwayelaboration_beta* corresponds to beta in the pathway elaboration method [1] 
  - Set *pathwayelaboration_K* Corresponds to K in the pathwaye elaboration method [1]. 
  - Set *pathwayelaboration_beta* Corresponds to beta in the pathwaye elaboration method [1]. 
  - Set *pathwayelaboration_use_elaboration=1* to include the elaboration step [1]. 
  
 - In ' map.py ': 
   - Set the initial parameter set. 
 - Run  ' learndnakinetics/map.py '. 
 
  If you have successfully installaled this software,  you will get a message such as "Starting to estimate parameters for the Arrhenius kinetic model". 

# References 

[1]  Sedigheh Zolaktaf*, Frits Dannenberg*, Mark Schmidt, Anne Condon, and Erik Winfree, The Pathway Elaboration Method for Mean First Passage Time Estimation in Large Continuous-Time Markov Chains with Applications to Nucleic Acid Kinetics, under review.

[2] Schaeffer JM, Thachuk C, Winfree E. Stochastic simulation of the kinetics of multiple interacting nucleic acid strands. InInternational Workshop on DNA-Based Computers 2015 Aug 17 (pp. 194-211). Springer, Cham.

[3] Bonnet G, Krichevsky O, Libchaber A. Kinetics of conformational fluctuations in DNA hairpin-loops. Proceedings of the National Academy of Sciences. 1998 Jul 21;95(15):8602-6.
 
[4] Cisse II, Kim H, Ha T. A rule of seven in Watson-Crick base-pairing of mismatched sequences. Nature structural & molecular biology. 2012 Jun;19(6):623.
 
[5] Hata H, Kitajima T, Suyama A. Influence of thermodynamically unfavorable secondary structures on DNA hybridization kinetics. Nucleic acids research. 2018 Jan 25;46(2):782-91.

[6] Machinek RR, Ouldridge TE, Haley NE, Bath J, Turberfield AJ. Programmable energy landscapes for kinetic control of DNA strand displacement. Nature communications. 2014 Nov 10;5(1):1-9.

[7] Zhang JX, Fang JZ, Duan W, Wu LR, Zhang AW, Dalchau N, Yordanov B, Petersen R, Phillips A, Zhang DY. Predicting DNA hybridization kinetics from sequence. Nature chemistry. 2018 Jan;10(1):91-8.
