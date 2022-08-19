# RPC Project
This repo contains the code for the controller estimation of a study regarding 
Myoelectric Control Architectures to Drive Upper Limb Exoskeletons.
 
**THE GOAL OF THE PROJECT IS TRY TO ESTIMATE A GENERAL CONTROLLER WHICH ADAPTS WELL TO ALL CONTROL ARCHITECTURES USED IN THE STUDY**

List of control architecture and adapters used to collect the data used for the estimation

```    
COMP        -   NONE               
FORCE       -   PLAIN_P
POS_V       -   PLAIN_P
FIX_IMP     -   PLAIN_P
ADM         -   PLAIN_P
FORCE_INT   -   PLAIN_P
FORCE       -   MULTICH8
POS_V       -   MULTICH8
FIX_IMP     -   MULTICH8
ADM         -   MULTICH8
FORCE_INT   -   MULTICH8      
```

The data used in the scripts are defined inside the `myo_tools_testing\testing_version1\settings\logs_eval_gains_2019_10_09_10_18__2020_06_26_30` file

Two techniques were used:

* Controller estimation from reference position as input and torque as output
* Estimation of the whole model W, Controller C plus Mechanical system G, assuming a PD controller. Then performing algebraic operations to extract C

To have better results in second technique, the estimation of G was faced at first to validate our assumptions

## Data Cleaning
In the scripts were G and C were estimated, data cleaning operations were performed before the actual estimation

An experiment was considered an outlier if it was charaterized by one of the following aspects 

* The last value of the torque signal to be outside a small range.
* Considered an outlier by the `rmoutliers` matlab function, so for each temporal istant, if the mean of all the torque value of that istant, is more than 3 standard deviations computed for the istant

After the cleaning operation, the signal were trimmed at start to remove the non relevant portion of the signal

## Estimation pipeline
To estimate the models, the `tfest` function was exploited
From each experiment contained in the cleaned set, a model was estimated
Then every experiment of the cleaned set was used to validate the computed models

## Mechanical Model Estimation - File `G_estimation` 

We wanted to validate our assumption where G has the following form `1/(J*s^2+d*s)` 
To estimate G 3 approaches were used, the results of these approaches were saved in the `resultStructures` folder

1) Estimation from torque to position 
...This was the first estimation approach tested - results saved in the `G_resultTable` variable

2) Estimation from torque to position and setting the end tail of the torque signal to zero 
...From the first approach we noticed that sometimes the posisition output from the estimated mechanical model diverged even if the torque was set costant at the end of the signal for a consistent amount of time. The reason for this behaviour is that the mechanical model contains an integrator, so even if the torque signal is set to costant, the error accumulates making the position diverge - results saved in the `G_resultTable_zero` variable

3) Estimation from torque to velocity 
...Given the considerations discusses at point 2, we tried to estimate the mechanical model from torque and velocity, excluding the estimation of the pole in 0, the one responsible of the integration in the control architecture - results saved in the `G_resultTable_vel` variable

### G estimation Results

## Controller Estimation

### C estimation Results

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
