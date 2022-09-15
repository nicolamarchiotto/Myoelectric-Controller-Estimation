%% C ESTIMATION SCRIPT
% DATA IMPORT
clc
clear
addpath(genpath('..\\myo_tools_testing'))
logs_eval_gains_2019_10_09_10_18__2020_06_26_30 %inside testing1\settings folder
batchMultProcessingTest1
batchMultElaborationTest1Long
actualFixedDataTableExpanded = fixedDataLongTableExpanded;

% logs_eval_2019_09_06 %inside settings folder
% batchProcessingTest1
% batchElaborationTest1
% actualFixedDataTableExpanded = fixedDataTableExpanded;

close all

pathToGitFolder = 'C:\\Users\\nicol\\Desktop\\rpcProject\\myo_tools_testing\\RPC_MN_project\\';

clearvars -except exoRefSplineCells actualFixedDataTableExpanded pathToGitFolder


%% DEFINITION OF G RESULT TABLE
C_sz = [1 9];
W_sz = [1 10];
C_freq_sz = [1 8];

C_varTypes = {'char', 'double', 'double', 'double', 'cell', 'double', 'cell', 'cell', 'cell'};
C_varNames = {'Architecture', '# zeros', '# poles', '# est exps', 'C model', 'C fit', 'C num', 'C den', 'C poles, (s+p)'};
C_freq_varTypes = {'char', 'double', 'double', 'cell', 'cell', 'cell', 'cell', 'cell'};

W_varTypes = {'char', 'double', 'double', 'double', 'cell', 'double', 'cell', 'cell', 'cell', 'cell'};
W_varNames = {'Architecture', '# zeros', '# poles', '# est exps', 'W model', 'W fit', 'C model', 'C num', 'C den', 'C poles, (s+p)'};
C_freq_varNames = {'Architecture', '# zeros', '# poles', 'C model', 'C num', 'C den', 'C zeros, (s+z)','C poles, (s+p)'};

C_resultTable = table('Size', C_sz, 'VariableTypes', C_varTypes, 'VariableNames', C_varNames);
W_resultTable = table('Size', W_sz, 'VariableTypes', W_varTypes, 'VariableNames', W_varNames);
C_freq_resultTable = table('Size', C_freq_sz, 'VariableTypes', C_freq_varTypes, 'VariableNames', C_freq_varNames);

clear C_varTypes C_varNames C_sz W_varTypes W_varNames W_sz C_freq_varTypes C_freq_varNames C_freq_sz

%% SAVE G RESULT TABLE - To store results of experiments up to now
save(pathToGitFolder + "resultStructures\\C_resultTable.mat", 'C_resultTable')
save(pathToGitFolder + "resultStructures\\W_resultTable.mat", 'W_resultTable')
save(pathToGitFolder + "resultStructures\\C_freq_resultTable.mat", 'C_freq_resultTable')

%% LOAD G RESULT TABLE - If already defined system
pathToGitFolder = 'C:\\Users\\nicol\\Desktop\\rpcProject\\myo_tools_testing\\RPC_MN_project\\';

load(pathToGitFolder + "resultStructures\\C_resultTable.mat", 'C_resultTable')
load(pathToGitFolder + "resultStructures\\W_resultTable.mat", 'W_resultTable')
load(pathToGitFolder + "resultStructures\\C_freq_resultTable.mat", 'C_freq_resultTable')

%
%
%
%% DATA FILTERING: Architecture - Decoder
%
%
%

% GRAVITY COMP - None -> 69 exps 
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "COMP" & actualFixedDataTableExpanded.decoder == "NONE";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.COMP_NONE;
trimStartIndex = 10;

%% FORCE - PLAIN_P -> 44 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FORCE" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.FORCE_PLAIN_P;
trimStartIndex = 30;

%% POS_V - PLAIN_P -> 36 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "POS_V" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.POS_V_PLAIN_P;
trimStartIndex = 30;

%% FIX_IMP - PLAIN_P -> 45 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FIX_IMP" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.FIX_IMP_PLAIN_P;
trimStartIndex = 10;

%% ADM - PLAIN_P -> decent data after cleaning - 43 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "ADM" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.ADM_PLAIN_P;
trimStartIndex = 10;

%% FORCE_INT - PLAIN_P -> 33 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FORCE_INT" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.FORCE_INT_PLAIN_P;
trimStartIndex = 10;

%% FORCE - MULTICH8 -> 47 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FORCE" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.FORCE_MULTICH8;
trimStartIndex = 10;

%% POS_V - MULTICH8 -> 42 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "POS_V" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.POS_V_MULTICH8;
trimStartIndex = 10;

%% FIX_IMP - MULTICH8 -> 42 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FIX_IMP" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.FIX_IMP_MULTICH8;
trimStartIndex = 10;

%% ADM - MULTICH8 -> 50 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "ADM" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.ADM_MULTICH8;
trimStartIndex = 10;

%% FORCE_INT - MULTICH8 -> 46 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FORCE_INT" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.FORCE_INT_MULTICH8;
trimStartIndex = 1;

%%
%
%
%
%
clearvars -except exoRefSplineCells actualFixedDataTableExpanded pathToGitFolder trimStartIndex architecture selected_indeces selection C_resultTable W_resultTable C_freq_resultTable 

C_est_divided_2