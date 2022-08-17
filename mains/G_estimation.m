%% DATA IMPORT
clc
clear
addpath(genpath('..\\myo_tools_testing'))
logs_eval_gains_2019_10_09_10_18__2020_06_26_30 %inside testing1\settings folder
batchMultProcessingTest1
batchMultElaborationTest1Long
close all
%%
clearvars -except fixedDataLongTableExpanded exoRefSplineCells

actualFixedDataTableExpanded = fixedDataLongTableExpanded;
clc
%% DEFINITION OF G RESULT TABLE
sz = [1 9];

varTypes = {'char', 'double', 'double', 'double', 'cell', 'double', 'cell', 'cell', 'cell'};
varNames = {'Type', '# zeros', '# poles', '# est exps', 'est model', 'est fit', 'est num', 'est den', 'est poles, (s+p)'};
G_resultTable = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);

clear varTypes varNames sz

%% SAVE G RESULT TABLE - To store results of experiments up to now
save G_resultTable

%% LOAD G RESULT TABLE - If already defined system
load G_resultTable.mat


%% DATA FILTERING: Architecture - Decoder

% GRAVITY COMP - None -> good data - 64 exps 
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "COMP" & actualFixedDataTableExpanded.decoder == "NONE";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.COMP_NONE;
torqueMaxAllowed = 5;
torqueMinAllowed = -5;

trimStartIndex = 100;

%% FORCE - PLAIN_P -> bad data after cleaning - 35 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FORCE" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.FORCE_PLAIN_P;
torqueMaxAllowed = 1;
torqueMinAllowed = -1;
trimStartIndex = 220;

%% POS_V - PLAIN_P -> decent data after cleaning - 22 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "POS_V" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.POS_V_PLAIN_P;
torqueMaxAllowed = 1.5;
torqueMinAllowed = -0.5;
trimStartIndex = 80;
%% FIX_IMP - PLAIN_P -> decent data after cleaning - 29 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FIX_IMP" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.FIX_IMP_PLAIN_P;
torqueMaxAllowed = 1.5;
torqueMinAllowed = -0.5;
trimStartIndex = 45;

%% ADM - PLAIN_P -> decent data after cleaning - 24 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "ADM" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.ADM_PLAIN_P;
torqueMaxAllowed = 0.5;
torqueMinAllowed = -0.075;
trimStartIndex = 50;

%% FORCE_INT - PLAIN_P -> 26 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FORCE_INT" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.FORCE_INT_PLAIN_P;
torqueMaxAllowed = 5;
torqueMinAllowed = -1;
trimStartIndex = 50;


%% FORCE - MULTICH8 -> bad data  - 23 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FORCE" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.FORCE_MULTICH8;
torqueMaxAllowed = 0.5;
torqueMinAllowed = -0.2;
trimStartIndex = 10;

%% POS_V - MULTICH8 -> good data after cleaning - 25 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "POS_V" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.POS_V_MULTICH8;
torqueMaxAllowed = 1.5;
torqueMinAllowed = -0.5;
trimStartIndex = 50;

%% FIX_IMP - MULTICH8 -> good data after cleaning - 28 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FIX_IMP" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.FIX_IMP_MULTICH8;
torqueMaxAllowed = 1.5;
torqueMinAllowed = -0.5;
trimStartIndex = 90;

%% ADM - MULTICH8 -> bad data after cleaning - 38 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "ADM" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.ADM_MULTICH8;
torqueMaxAllowed = 0.3;
torqueMinAllowed = -0.05;
trimStartIndex = 110;

%% FORCE_INT - MULTICH8 -> 46 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FORCE_INT" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
estimationCase = EstimationCasesEnum.FORCE_INT_MULTICH8;
torqueMaxAllowed = 5;
torqueMinAllowed = -5;
trimStartIndex = 60;

%% IMPORTANT PARAMETERS

Ts = 0.01;

torqueEndEpsilon = 0.2;

maxSignalLength = 400;

G_plot_signals = true;

% variable for removing or not torques which min or max values are outside
% the torqueMaxAllowed amd torqueMinAllowed declared in DATA FILTERING section
torqueRangeFiltering = true;

% params for discarding torques signal which decrease after the trimming operation
discardDecreaseStartingTorque = true;
discardDecreaseStartingTorqueIdx = 30;

% outliers data filtering based on standard deviation, 
stdOutlierRemoval = true;

% Align both pos at torque at zero
alignAtZero = false;

% DATA CLEANING
clc
close all

allTrimmedTorque = [];
allTrimmedPos = [];

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

% Find the torque index of the further max value of the entire set
maxIdx = 0;

% Find greater max torque value idx of the various signals 
for expId = selected_indeces 
    % vector of torques
    tagIdx = fixedDataLongTableExpanded(expId,:).tag_idx;

    torque = exoRefSplineCells{tagIdx}.reference(fixedDataLongTableExpanded(expId,:).start_idx:fixedDataLongTableExpanded(expId,:).end_idx);

    [maxVal, idx] = max(torque);
    if(idx > maxIdx)
        maxIdx = idx;
    end
end
clear tagIdx expId idx;

i = 1;
for expId = selected_indeces  %for each experiment
    % vector of torques
    tagIdx = fixedDataLongTableExpanded(expId,:).tag_idx;
    
    %vector of positions of the experiment
    position = exoRefSplineCells{tagIdx}.thetaE(fixedDataLongTableExpanded(expId,:).start_idx:fixedDataLongTableExpanded(expId,:).end_idx);
    
    % vector of torques
    torque = exoRefSplineCells{tagIdx}.reference(fixedDataLongTableExpanded(expId,:).start_idx:fixedDataLongTableExpanded(expId,:).end_idx);

    % Torque outliers elimination
    if abs(torque(end)) > torqueEndEpsilon
        continue;
    end
    if torqueRangeFiltering && (max(torque) > torqueMaxAllowed || min(torque) < torqueMinAllowed)
        continue;
    end
    if discardDecreaseStartingTorque && torque(1) > torque(discardDecreaseStartingTorqueIdx)
        continue;
    end
    
    [maxVal, idx] = max(torque);
    
    pointsToAddStart = maxIdx - idx; 
    
    firstValOfTorque = torque(1);
    firstValOfPos = position(1);
    
    lastValOfTorque = torque(end);
    lastValOfPos = position(end);
    
    % signal where the torque max value is at the same index for all exps
    augmTorque = [ones(pointsToAddStart,1)*firstValOfTorque; torque];
    augmPos = [ones(pointsToAddStart,1)*firstValOfPos; position];
  
    adjustedTorque = []; 
    adjustedPos = [];
    
    if(length(augmTorque) < maxSignalLength)
        pointsToAddAtEnd = maxSignalLength - length(augmTorque);
        adjustedTorque = [augmTorque; ones(pointsToAddAtEnd,1)*lastValOfTorque];
        adjustedPos = [augmPos; ones(pointsToAddAtEnd,1)*lastValOfPos]; 
    else
        adjustedTorque = augmTorque(1:maxSignalLength);
        adjustedPos = augmPos(1:maxSignalLength);
    end
    
    adjustedTorque = [adjustedTorque; ones(maxSignalLength/2,1)*adjustedTorque(end)];
    adjustedPos = [adjustedPos; ones(maxSignalLength/2,1)*adjustedPos(end)];
    
    allTrimmedTorque(i,:) = adjustedTorque;
    allTrimmedPos(i,:) = adjustedPos;
    
    i = i+1;
    if G_plot_signals 
        plot(hAx1,torque);
        plot(hAx2,position);
    end
end

% STANDARD DEVIATION OUTLIER REMOVAL
% Remove outliers of a vector where an outlier is defined as a point more 
% than three standard deviations from the mean of the data.
if stdOutlierRemoval
    [B,TF]=rmoutliers(allTrimmedTorque,'mean');
    idcs = flip(find(TF)');

    for expId = idcs 
        allTrimmedTorque(expId,:) = [];
        allTrimmedPos(expId,:) = [];
    end
end

% Align at zero both position and torque, methods tested for obtaining
% better starting conditions for responses.

if alignAtZero
    for expId=1:1:size(allTrimmedTorque,1)
        torqFirstVal=allTrimmedTorque(expId,1);
        posFirstVal=allTrimmedPos(expId,1);

        allTrimmedTorque(expId,:) = allTrimmedTorque(expId,:)-torqFirstVal;
        allTrimmedPos(expId,:) =  allTrimmedPos(expId,:)-posFirstVal;
    end
end


if alignAtZero
    for expId=1:1:size(allTrimmedTorque,1)
        torqFirstVal=allTrimmedTorque(expId,1);
        posFirstVal=allTrimmedPos(expId,1);

        allTrimmedTorque(expId,:) = allTrimmedTorque(expId,:)-torqFirstVal;
        allTrimmedPos(expId,:) =  allTrimmedPos(expId,:)-posFirstVal;
    end
end

% Trim at start 
supTorque = [];
supPos = [];
for expId=1:1:size(allTrimmedTorque,1)
    supTorque(expId,:) = allTrimmedTorque(expId, trimStartIndex:end);
    supPos(expId,:) =  allTrimmedPos(expId, trimStartIndex:end);
end
allTrimmedTorque = supTorque;
allTrimmedPos = supPos;
clear supTorque supPos;

if G_plot_signals 
    for testIdx=1:1:size(allTrimmedTorque,1)
        plot(hAx3, allTrimmedTorque(testIdx, :));
        plot(hAx4, allTrimmedPos(testIdx, :));
    end
end


size(allTrimmedPos,1)
%% Model Estimation

clc;

% Number of poles and zeros for the estimated models
G_num_zeros = 0;
G_num_poles = 2;

% For each esperiment estimate G, O(n) and construct iddata structures
G_iddata = {};
G_sys_est = {};

% option to force estimated model to be stable
opt = tfestOptions('EnforceStability',true,'InitialCondition','estimate');

for expIdx=1:1:size(allTrimmedTorque,1)
    G_est_iddata = iddata(allTrimmedPos(expIdx,:)', allTrimmedTorque(expIdx,:)', Ts); 
    G_iddata{expIdx} = G_est_iddata;
    G_sys_est{expIdx} = tfest(G_est_iddata, G_num_poles, G_num_zeros,opt); 
end

% initialization of variables for plot function
G_bestModelOutput = {};
clear expIdx G_est_iddata opt;

% Find the best G testing on all experiments, O(n^2)
%
%
% CAREFULL, each bestModelFinder is O(n^2)
%
%
[G_bestModel, G_bestModelFit, G_bestModelOutput] = bestModelFinder(G_sys_est, G_iddata);

% RESULTS
clc
fprintf('\nG: Best model fit from single experiment estimation: %.3f\n', G_bestModelFit);
%
% s = tf('s');
% J = 0.068;
% D = 0.1; %da 0.1 a 0.001
% G_assumed = zpk(zero(zpk(1/(J*s^2 + D*s))),pole(zpk(1/(J*s^2 + D*s))),1);
% [G_Assumed_num, G_Assumed_den] = tfdata(G_assumed, 'v');
G_best = zpk(zero(G_bestModel),pole(G_bestModel),1);
[G_best_num, G_best_den] = tfdata(G_best, 'v');

G_resultTable(estimationCase, :) = {char(estimationCase), G_num_zeros, G_num_poles, size(allTrimmedTorque,1), {G_best}, G_bestModelFit, {G_best_num}, {G_best_den}, {pole(G_best)}};

% PLOTS
close all
imageSavePath = "C:\\Users\\nicol\\Desktop\\rpcProject\\myo_tools_testing\\RPC_MN_project\\images\\";
saveImage = true;
plotFunction(false, false, true, allTrimmedTorque, allTrimmedPos, {}, {}, G_bestModelOutput, saveImage, imageSavePath, estimationCase)