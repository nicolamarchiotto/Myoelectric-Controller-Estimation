# RPC Project
This repo contains the code for the controller estimation of a study regarding 
Myoelectric Control Architectures to Drive Upper Limb Exoskeletons.

## Purpose 
**THE GOAL OF THE PROJECT IS TRY TO ESTIMATE A GENERAL CONTROLLER WHICH ADAPTS WELL TO THE ARCHITECTURES USED IN THE STUDY**

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

The data used in the scripts are defined inside the file `myo_tools_testing\testing_version1\settings\logs_eval_gains_2019_10_09_10_18__2020_06_26_30`

The repository also contains code for the mechanical part G estimation

Full documentation here [.](/docs/Documentation.pdf)

## External resources and important guidelines
Davide Constanzi's PhD Thesis: https://iris.univr.it/handle/11562/1061781

N.B. This repository is not working by itself but it requires the myo_tools_testing
repository and  myo_interface log data, the first of which is present on GitLab, 
the access to this repository require access to be granted by the owner Davide Costanzi
myo_testing_tools repo: https://gitlab.com/davide.costanzi/myo_tools_testing
myo_interafce data: https://univr-my.sharepoint.com/:f:/g/personal/davide_costanzi_univr_it/EjxulHf-4XBNh5BEWGRDnK4B22zxYxDER_DAsW-j8O49UA?e=efy68k

To work the folder containig this repository must be placed into the myo_testing_tools folder,
the add to path with subfolders both myo_interface and myo_testing_tools folders
