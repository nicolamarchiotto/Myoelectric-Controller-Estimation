# RPC_module_A Project
This repo contains the code for the controller estimation of a study regarding 
Myoelectric Control Architectures to Drive Upper Limb Exoskeletons.

The entry point of the project is the project_myo,m file.

In this file the following actions are performed:

* Data import from the logs files

* Selection and filtering for choosing the data of the decode and control architecture used in the experiments
    The cases should be the following
 
    CTRLArchitecture    Decoder     
    
    COMP                NONE (MULTI8)
    FORCE_INT           MULTI2/MULTI8
    IMPEDANCE           MULTI2/MULTI8
    VELOCITY/POS        MULTI2/MULTI8

    Ask Davide, the informations about CTRLArchitecture and Decoder are not right

# For Each of the follwing architecures
    C                   M = C*G
    1 poles             3 poles
    1 zeros, 1 poles    1 zeros, 3 poles
    1 zeros, 2 poles    1 zeros, 4 poles
    2 zeros, 2 poles    2 zeros, 4 poles 
    
    - Identification of the  brain controller C and overall plant C*G
        Chose randomly 80% of  the experiments, and for each estimate the plant and controller using the single experiment data
    - Use the remaining 20% of the experiments to validate the estimated controllers and plants and choose the best one
    - Save the best controller and plant for the given architecture

# Use the best model for each architectures and use all the data to validate them, the ones which give the best score will be
the best controller and plant for the given pair CTRLArchitecture/Decoder, C_best and P_best

* Test C_best and P_best on simulink architectures

Link to external resources: 
Davide Constanzi's PhD Thesis: https://iris.univr.it/handle/11562/1061781

N.B. This repository is not working by itself but it requires the myo_tools_testing
repository and  myo_interface log data, the first of which is present on GitLab, 
the access to this repository require access to be granted by the owner Davide Costanzi
myo_testing_tools repo: https://gitlab.com/davide.costanzi/myo_tools_testing
myo_interafce data: https://univr-my.sharepoint.com/:f:/g/personal/davide_costanzi_univr_it/EjxulHf-4XBNh5BEWGRDnK4B22zxYxDER_DAsW-j8O49UA?e=efy68k

To work the folder containig this repository must be placed into the myo_testing_tools folder,
the add to path with subfolders both myo_interface and myo_testing_tools folders