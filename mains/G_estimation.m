%% DATA IMPORT
clc
clear
addpath(genpath('..\\myo_tools_testing'))
logs_eval_gains_2019_10_09_10_18__2020_06_26_30 %inside testing1\settings folder
batchMultProcessingTest1
batchMultElaborationTest1Long
close all
clearvars -except fixedDataLongTableExpanded exoRefSplineCells

actualFixedDataTableExpanded = fixedDataLongTableExpanded;

%% DEFINITION OF RESULT TABLE
sz = [1 11];
varTypes = {'char', 'double', 'double', 'double', 'cell', 'double', 'cell', 'cell', 'cell', 'cell', 'cell'};
varNames = {'Type', '# zeros', '# poles', '# est exps', 'estimated', 'est fit', 'est num', 'est den', 'est poles, (s+p)', 'assumed', 'assumed poles (s+p)'};
resultTable = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);

clear varTypes varNames sz
%% SAVE RESULT TABLE - To store results of experiments up to now
save resultTable

%% LOAD RESULT TABLE - If already defined system
load resultTable.mat


%% DATA FILTERING: Architecture - Decoder

% GRAVITY COMP - None -> good data - 51 exps 
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "COMP" & actualFixedDataTableExpanded.decoder == "NONE";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.COMP;
torqueMaxAllowed = 5;
torqueMinAllowed = -5;
%% FORCE - PLAIN_P -> bad data after cleaning - 30 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FORCE" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.FORCE_PLAIN_P;
torqueMaxAllowed = 1;
torqueMinAllowed = -0.07;
%% POS_V - PLAIN_P -> decent data after cleaning - 31 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "POS_V" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.POS_V_PLAIN_P;
torqueMaxAllowed = 1.5;
torqueMinAllowed = -0.1;
%% FIX_IMP - PLAIN_P -> decent data after cleaning - 29 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FIX_IMP" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.FIX_IMP_PLAIN_P;
torqueMaxAllowed = 1;
torqueMinAllowed = -0.1;
%% ADM - PLAIN_P -> decent data after cleaning - 30 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "ADM" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.ADM_PLAIN_P;
torqueMaxAllowed = 0.3;
torqueMinAllowed = 0.015;
%% FORCE_INT - PLAIN_P -> what is torque doing? seems different gains or target angle - 21 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FORCE_INT" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.FORCE_INT_PLAIN_P;
torqueMaxAllowed = 5;
torqueMinAllowed = -5;
%% FORCE - MULTICH8 -> bad data  - 44 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FORCE" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.FORCE_MULTICH8;
torqueMaxAllowed = 0.5;
torqueMinAllowed = -0.2;
%% POS_V - MULTICH8 -> good data after cleaning - 38 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "POS_V" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.POS_V_MULTICH8;
torqueMaxAllowed = 1.5;
torqueMinAllowed = -0.2;
%% FIX_IMP - MULTICH8 -> good data after cleaning - 49 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FIX_IMP" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.FIX_IMP_MULTICH8;
torqueMaxAllowed = 1.5;
torqueMinAllowed = -0.2;
%% ADM - MULTICH8 -> bad data after cleaning - 51 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "ADM" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.ADM_MULTICH8;
torqueMaxAllowed = 0.3;
torqueMinAllowed = -0.05;
%% FORCE_INT - MULTICH8 -> what is torque doing? seems different gains or target angle - 27 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FORCE_INT" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.FORCE_INT_MULTICH8;
torqueMaxAllowed = 5;
torqueMinAllowed = -5;
%% IMPORTANT PARAMETERS

Ts = 0.01;

posEpsilon = 0.01;
torqueEndEpsilon = 0.2;

maxSignalLength = 140;

G_plot_signals = true;

% variable for removing or not torques which min or max values are outside
% the torqueMaxAllowed amd torqueMinAllowed declared in DATA FILTERING section
torqueRangeFiltering = false;

% params for discarding torques signal which decrease after the trimming operation
discardDecreaseStartingTorque = false;
discardDecreaseStartingTorqueIdx = 15;

% DATA CLEANING
close all
clc
allTrimmedTorque = {};
allTrimmedPos = {};

i = 1;

if G_plot_signals
    figure; 
    hAx1 = subplot(2,1,1);
    title('Torque output from experiments')
    hAx2 = subplot(2,1,2);
    title('Position of the experiments') 
    figure; 
    hAx3 = subplot(2,1,1);
    title('FILT: Torque output from experiments')
    hAx4 = subplot(2,1,2);
    title('FILT: Position of the experiments') 
    axes(hAx1)
    hold on
    axes(hAx2)
    hold on
    axes(hAx3)
    hold on
    axes(hAx4)
    hold on
end

for expId = selected_indeces  %for each experiment
    tagIdx = fixedDataLongTableExpanded(expId,:).tag_idx;
    
    %vector of positions of the experiment
    position = exoRefSplineCells{tagIdx}.thetaE(fixedDataLongTableExpanded(expId,:).start_idx:fixedDataLongTableExpanded(expId,:).end_idx);
    
    % vector of torques
    torque = exoRefSplineCells{tagIdx}.reference(fixedDataLongTableExpanded(expId,:).start_idx:fixedDataLongTableExpanded(expId,:).end_idx);

    trimIdx = 1;
    
    % Data cleaning, remove data of test before actual motion
    while (abs(position(1) - position(trimIdx+1)) < posEpsilon)   
       trimIdx = trimIdx + 1; 
    end
    
    trimmedPos = position(trimIdx:trimIdx+maxSignalLength);
    
    
    % Reference signal cut to correspond with trimmedPos length
    trimmedTorque = torque(trimIdx:trimIdx+maxSignalLength);
    
    % outliers torque elimination
    if abs(trimmedTorque(end)) > torqueEndEpsilon
        continue;
    end
    if torqueRangeFiltering && (max(trimmedTorque) > torqueMaxAllowed || min(trimmedTorque) < torqueMinAllowed)
        continue;
    end
    if discardDecreaseStartingTorque && trimmedTorque(1) > trimmedTorque(discardDecreaseStartingTorqueIdx)
        continue;
    end
    
    
    % Take the last value of the filtered data and lengthening the signals
    for k=1:1:3*maxSignalLength
       trimmedTorque = [trimmedTorque; trimmedTorque(end)];
       trimmedPos = [trimmedPos; trimmedPos(end)];
    end
    
    if G_plot_signals 
        plot(hAx1,torque);
        plot(hAx2,position);
        plot(hAx3,trimmedTorque);
        plot(hAx4,trimmedPos);
    end
    
    allTrimmedTorque{i} = trimmedTorque;
    allTrimmedPos{i} = trimmedPos;
    
    i = i + 1;
end
% length(allTrimmedTorque)
clear trimIdx expId i k;
length(selected_indeces)
length(allTrimmedTorque)
%% Model Estimation

clc;

% Number of poles and zeros for the estimated models
G_num_zeros = 0;
G_num_poles = 2;

% construct iddata structures
G_iddata = {};

for testIdx=1:1:length(allTrimmedTorque)
    G_iddata{testIdx} = iddata(allTrimmedPos{testIdx}, allTrimmedTorque{testIdx}, Ts); 
end

clear testIdx;

% For each esperiment estimate G, O(n)
G_sys_est = {};

% option to force estimated model to be stable
opt = tfestOptions('EnforceStability',true,'InitialCondition','estimate');

for expIdx=1:1:length(allTrimmedPos)
    G_est_iddata = iddata(allTrimmedPos{expIdx}, allTrimmedTorque{expIdx}, Ts); 
    G_sys_est{expIdx} = tfest(G_est_iddata, G_num_poles, G_num_zeros,opt); 
end

% initialization of variables for plot function
G_bestModelOutput = {};
clear expIdx;

% Find the best G testing on all experiments, O(n^2)
%
%
%% CAREFULL, each bestModelFinder is O(n^2)
%
%
[G_bestModel, G_bestModelFit, G_bestModelOutput] = oldBestModelFinder(G_sys_est, G_iddata);

% RESULTS
clc
fprintf('\nG: Best model fit from single experiment estimation: %.3f\n', G_bestModelFit);

s = tf('s');
J = 0.068;
D = 0.1; %da 0.1 a 0.001
G_assumed = zpk(zero(zpk(1/(J*s^2 + D*s))),pole(zpk(1/(J*s^2 + D*s))),1);
[G_Assumed_num, G_Assumed_den] = tfdata(G_assumed, 'v');
G_best = zpk(zero(G_bestModel),pole(G_bestModel),1);
[G_best_num, G_best_den] = tfdata(G_best, 'v');

%% PLOTS
clc
close all;
oldPlotFunction(false, false, true, allTrimmedTorque, allTrimmedPos, {}, {}, G_bestModelOutput)
