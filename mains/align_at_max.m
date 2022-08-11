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

if G_plot_signals 
    for testIdx=1:1:size(allTrimmedTorque,1)
        plot(hAx3, allTrimmedTorque(testIdx, :));
        plot(hAx4, allTrimmedPos(testIdx, :));
    end
end
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
    G_est_iddata = iddata(allTrimmedTorque(expIdx,:)', allTrimmedPos(expIdx,:)', Ts); 
    G_iddata{expIdx} = G_est_iddata;
    G_sys_est{expIdx} = tfest(G_est_iddata, G_num_poles, G_num_zeros,opt); 
end

% initialization of variables for plot function
G_bestModelOutput = {};
clear expIdx G_sys_est opt;

% Find the best G testing on all experiments, O(n^2)
%
%
% CAREFULL, each bestModelFinder is O(n^2)
%
%
[G_bestModel, G_bestModelFit, G_bestModelOutput] = bestModelFinder(G_sys_est, G_iddata);

%% RESULTS
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
close all
plot_C = false;
plot_W = false;
plot_G = true;
C_bestModelOutput={}; 
W_bestModelOutput={};
plotFunction(plot_C, plot_W, plot_G, allTrimmedTorque, allTrimmedPos, C_bestModelOutput, W_bestModelOutput, G_bestModelOutput)
