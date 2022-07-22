# RPC Project
This repo contains the code for the controller estimation of a study regarding 
Myoelectric Control Architectures to Drive Upper Limb Exoskeletons.

The entry point of the project is the project_myo,m file.

In this file the following actions are performed:

1) Data import from the logs files, 
Chose data of the experiments acquired using a certain control architecture and decoder
Controller estimations from single and multiples log files were tried

* List of control architecture used for data collection
```    
COMP                NONE (MULTI8)
FORCE
FORCE_INT           MULTI2/MULTI8
ADMITTANCE
IMPEDANCE           MULTI2/MULTI8
VELOCITY/POS        MULTI2/MULTI8
```

* List of control architecture used for data collection
```
NONE (MULTI8)
MULTI8
PLAIN (with different gains)
```

2) Data Cleaning: temporally align position and torque signals to have only the actual motion considered, remove outliers by get rod of experiments which outliers do not end up in zero +- a given theshold

3) Controller estimation techniques tried:
* Estimate the controller using a single experiment
* Estimate the controller concatenating the experiment of the forward motions
* Estimate the controller using the concatenation method but flip one experiment every two to partially remove discontinuities
* Estimated the whole model W and then extract the controller exploting algebraic operations

* Assumed models for the estimation
```
C                   W = CG/(1+CG)

1 poles             3 poles
1 zeros, 1 poles    1 zeros, 3 poles
1 zeros, 2 poles    1 zeros, 4 poles
2 zeros, 2 poles    2 zeros, 4 poles 

G assumed to be of the following form 1/(Js^2+ds)
```   

* From the cleaned data, 80% of the experiments was used for the estimation, 20% for the testing
* Save the best controller C_best
 
3) Test the C_best estimated with the data corresponding to the pairs Architecture/Decoder with the Architectures listed at point 1 and save the step responses

4) The goal of the procedure is to try to estimate a controller which adapts well to all given architectures

Link to external resources: 
Davide Constanzi's PhD Thesis: https://iris.univr.it/handle/11562/1061781

N.B. This repository is not working by itself but it requires the myo_tools_testing
repository and  myo_interface log data, the first of which is present on GitLab, 
the access to this repository require access to be granted by the owner Davide Costanzi
myo_testing_tools repo: https://gitlab.com/davide.costanzi/myo_tools_testing
myo_interafce data: https://univr-my.sharepoint.com/:f:/g/personal/davide_costanzi_univr_it/EjxulHf-4XBNh5BEWGRDnK4B22zxYxDER_DAsW-j8O49UA?e=efy68k

To work the folder containig this repository must be placed into the myo_testing_tools folder,
the add to path with subfolders both myo_interface and myo_testing_tools folders
