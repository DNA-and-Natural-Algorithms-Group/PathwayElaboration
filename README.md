# Introduction

This software implements  the following two algorithms [1] for estimating parameters for nucleic acid reactions in the full state space of the Multistrand nucleic acid kinetic simulator [2]. 
 
- The Fixed Pathway Ensemble Inference (FPEI).
- The Stochastic Simulation Algrorithm Inference (SSAI).

 This software uses and modifies code from the <a href="https://github.com/DNA-and-Natural-Algorithms-Group/multistrand">Multistrand repository</a> [2]. 

For more information regarding the methods employed in our software, see our paper [1] (to appear). 

Please contact us at nasimzf@cs.ubc.ca for questions regarding the software.

Disclaimer: This software is free for use only for research purposes. If you make use of this software, please cite our paper

[1] Zolaktaf, Sedigheh, Frits Dannenberg, Erik Winfree, Alexandre Bouchard-Côté, Mark Schmidt, and Anne Condon. "Efficient Parameter Estimation for DNA Kinetics Modeled as Continuous-Time Markov Chains." In International Conference on DNA Computing and Molecular Programming, pp. 80-99. Springer, Cham, 2019.
# Requirements 

| Dependency       | Notes          
| ------------- |:-------------:|
|Python 2 |2.7.14 |
|c++11|gcc 4.8.5+|
|make|4.0+ | 
|<a href="http://www.nupack.org/">NUPACK</a> |3.2.1| 
|Multistrand|We have modified <a href="https://github.com/DNA-and-Natural-Algorithms-Group/multistrand">Multistrand</a>.  For this project,     <br />  use  the   modified version  provided in this directory| 
|numpy|1.9.3+| 
|SciPy | 1.1.0+ |
| enum34 | 1.1.6+ | 
| ConfigParser | 3.5.0+ | 
|matplotlib| 2.0.0+|
|seaborn|0.7.1+|



# Setup
For linux:
- Download and install the requirements. 
- Set the system environment variables to point NUPACKHOME to the folder where NUPACK is installed. For example, in openSUSE, you can run : ' export NUPACKHOME=path/to/NUPACK3.0.4 '. To verify that the NUPACKHOME is correctly set, run ' echo $NUPACKHOME '. 
- Clone our software directory into your workspace. 
- In the  multistrand_modified directory,
  - Build Multistrand by running ' make clean ' and then ' make '. Upon successful  building of  Multistrand you should get a message such as "Multistrand is now built. ... ".
  - Export  Multistrand as a Python library by running  ' make install '. 
  - Set the  PYTHONPATH evironment variable to point to the multistrand_modified  directory. For example, in openSUSE, you can use the following command: ' export PYTHONPATH="${PYTHONPATH}:path/to/multistrand_modified" '.  

  Note that this software only works with the modified Multistrand code provided here, and will not work with the original version of Multistrand.

# Package Tree
This software contains 5  directories, namely, multistrand_modified, learndnakinetics, reactions, dataset, plot. 
- learndnakinetics:  this  directory contains code to run the FPEI and SSAI parameter estimation algorithms on the following dataset. 
- dataset: this directory contains reaction rate constants and timescales that we compiled from published literature
  - hairpin_bonnet98: reactions 1-10 in Table 1 from [1]. These reactions are originally from [3]. 
  - helix4_cisse2012: reactions 11-15 in Table 1 from [1]. These reactions are originally from  [4].
  - helix2_hata2017: reactions 16-19  in Table 1 from [1]. These reactions are originally from  [5]. 
- multistrand_modfied: The modified Multistrand code adapted from the branch 'fdann-devel' from  <a href="https://github.com/DNA-and-Natural-Algorithms-Group/multistrand">  Multistrand repository</a> [2]. 
- plot: contains code to plot the results from parameter estimation, similar to Fig. 5 from [1].  

# Software Usage 
To use FPEI or SSAI with the dataset provided:
- In ' learndnakinetics/config_file.txt ': 
  - Set *use_FPEI_MFPT = 1* and *use_Gillespie_MFPT=0* to use FPEI. Set *use_FPEI_MFPT = 0* and *use_Gillespie_MFPT=1* to use SSAI. 
  - Set *rate_method = 3* to estimate parameters for the  Arrhenius kinetic model. Set *rate_method = 1* to estimate parameters for the Metropolis kinetic model. 
  - Set *parameter_folder* to be the path to a directory to save results. In this folder, the MSE of the parameters over the entire reactions of the dataset is logged. In addition, the parameter set of each iteration of the optimization will be logged. Also, the  fixed paths in FPEI will be saved.   
  - Set *n_processors* to be the number of processors for multiprocessing the computation of the objective function. Set *use_multiprocess =1* to multiprocess the computations. 
 - In ' map.py ': 
   - Set the initial parameter set. 
 - Run  ' learndnakinetics/map.py '. 
 
  If you have successfully installaled this software,  you will get a message such as "Starting to estimate parameters for the Arrhenius kinetic model". 

To plot the results, similar to Fig. 5 from [1], from parameter estimation using FPEI vs SSAI: 
 - in ' plot/plot_ssavsftei.py '
   - Provide the path to the *parameter_folder* directory that the FPEI results are saved.
   - Provide the path to the *parameter_folder* direcotory that the SSAI results are saved.
 - Run ' plot/plot_ssavsftei.py '
 

# References 

[1] Zolaktaf, Sedigheh, Frits Dannenberg, Erik Winfree, Alexandre Bouchard-Côté, Mark Schmidt, and Anne Condon. "Efficient Parameter Estimation for DNA Kinetics Modeled as Continuous-Time Markov Chains." In International Conference on DNA Computing and Molecular Programming, pp. 80-99. Springer, Cham, 2019.

[2] Schaeffer, J.M., Thachuk, C., Winfree, E.: Stochastic simulation of the kinetics of multiple interacting nucleic acid strands. In: Proceedings of the 21st International Conference on DNA Computing and Molecular Programming-Volume 9211 (2015)

 [3]  Bonnet, Grégoire, Oleg Krichevsky, and Albert Libchaber. "Kinetics of conformational fluctuations in DNA hairpin-loops." Proceedings of the National Academy of Sciences 95.15 (1998): 8602-8606.
 
 [4] Cisse, Ibrahim I., Hajin Kim, and Taekjip Ha. "A rule of seven in Watson-Crick base-pairing of mismatched sequences." Nature structural & molecular biology 19.6 (2012): 623.
 
  [5] Hata, Hiroaki, Tetsuro Kitajima, and Akira Suyama. "Influence of thermodynamically unfavorable secondary structures on DNA hybridization kinetics." Nucleic Acids Research (2017).
