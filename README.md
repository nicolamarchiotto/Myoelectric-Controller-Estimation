# RPC Project
This repo contains the code for the controller estimation of a study regarding 
Myoelectric Control Architectures to Drive Upper Limb Exoskeletons.
 
**THE GOAL OF THE PROJECT IS TRY TO ESTIMATE A GENERAL CONTROLLER WHICH ADAPTS WELL TO ALL CONTROL ARCHITECTURES USED IN THE STUDY**

In all the scripts of this repo the following action are performed

1) Data import from the logs files, 
Chose data of the experiments acquired using a certain control architecture and decoder
Controller estimations from single and multiples log files were tried

* List of control architecture used for data collection
```    
COMP               
FORCE
FORCE_INT           
ADM
FIX_IMP
POS_V           
```

* List of control architecture used for data collection
```
NONE (MULTI8)
MULTI8
PLAIN (with different gains)
```

2) Data Cleaning: temporally align position and torque signals to have only the actual motion considered, remove outliers by get rod of experiments which outliers do not end up in zero +- a given theshold

# concat_estimation.m
The script compares different estimation methods, actions performed in this script

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

The estimation from the single experiment was choosen to advance in the project

# main.m
Script were main work for estimation was performed, actions performed in this script

* The experiments were extracted from logs_eval_gains_2019_10_09_10_18__2020_06_26_30, session containing multiple sessions of experiments

* The controller C and whole model W were estimated using the single experiment technique explained in point 3)
* The data were cleaned in the same way explained at point 2)

* The best controller and whole model W were chosen testing on all the experiment set

* From the best model W, the controller was extracted using algebraic operations

* The best Controller and the controller extracted from W were tested on the simulink architectures

# External resources and important guidelines
Davide Constanzi's PhD Thesis: https://iris.univr.it/handle/11562/1061781

N.B. This repository is not working by itself but it requires the myo_tools_testing
repository and  myo_interface log data, the first of which is present on GitLab, 
the access to this repository require access to be granted by the owner Davide Costanzi
myo_testing_tools repo: https://gitlab.com/davide.costanzi/myo_tools_testing
myo_interafce data: https://univr-my.sharepoint.com/:f:/g/personal/davide_costanzi_univr_it/EjxulHf-4XBNh5BEWGRDnK4B22zxYxDER_DAsW-j8O49UA?e=efy68k

To work the folder containig this repository must be placed into the myo_testing_tools folder,
the add to path with subfolders both myo_interface and myo_testing_tools folders
