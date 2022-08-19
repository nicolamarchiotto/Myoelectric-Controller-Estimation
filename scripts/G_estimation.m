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
G_resultTable_zero = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);
G_resultTable_vel = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);

clear varTypes varNames sz

%% SAVE G RESULT TABLE - To store results of experiments up to now
save G_resultTable
save G_resultTable_zero
save G_resultTable_vel

%% LOAD G RESULT TABLE - If already defined system
load('G_resultTable.mat')
load('G_resultTable_zero.mat')
load('G_resultTable_vel.mat')

%% DATA FILTERING: Architecture - Decoder

% GRAVITY COMP - None -> good data - 64 exps 
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "COMP" & actualFixedDataTableExpanded.decoder == "NONE";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.COMP_NONE;
torqueMaxAllowed = 5;
torqueMinAllowed = -5;

trimStartIndex = 100;

%% FORCE - PLAIN_P -> bad data after cleaning - 35 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FORCE" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.FORCE_PLAIN_P;
torqueMaxAllowed = 1;
torqueMinAllowed = -1;
trimStartIndex = 30;

%% POS_V - PLAIN_P -> decent data after cleaning - 22 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "POS_V" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.POS_V_PLAIN_P;
torqueMaxAllowed = 1.5;
torqueMinAllowed = -0.5;
trimStartIndex = 30;
%% FIX_IMP - PLAIN_P -> decent data after cleaning - 29 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FIX_IMP" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.FIX_IMP_PLAIN_P;
torqueMaxAllowed = 1.5;
torqueMinAllowed = -0.5;
trimStartIndex = 45;

%% ADM - PLAIN_P -> decent data after cleaning - 24 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "ADM" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.ADM_PLAIN_P;
torqueMaxAllowed = 0.5;
torqueMinAllowed = -0.075;
trimStartIndex = 50;

%% FORCE_INT - PLAIN_P -> 26 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FORCE_INT" & actualFixedDataTableExpanded.decoder == "PLAIN_P";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.FORCE_INT_PLAIN_P;
torqueMaxAllowed = 5;
torqueMinAllowed = -1;
trimStartIndex = 50;


%% FORCE - MULTICH8 -> bad data  - 23 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FORCE" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.FORCE_MULTICH8;
torqueMaxAllowed = 0.5;
torqueMinAllowed = -0.2;
trimStartIndex = 10;

%% POS_V - MULTICH8 -> good data after cleaning - 25 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "POS_V" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.POS_V_MULTICH8;
torqueMaxAllowed = 1.5;
torqueMinAllowed = -0.5;
trimStartIndex = 50;

%% FIX_IMP - MULTICH8 -> good data after cleaning - 28 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FIX_IMP" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.FIX_IMP_MULTICH8;
torqueMaxAllowed = 1.5;
torqueMinAllowed = -0.5;
trimStartIndex = 90;

%% ADM - MULTICH8 -> bad data after cleaning - 38 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "ADM" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.ADM_MULTICH8;
torqueMaxAllowed = 0.3;
torqueMinAllowed = -0.05;
trimStartIndex = 110;

%% FORCE_INT - MULTICH8 -> 46 exps
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "FORCE_INT" & actualFixedDataTableExpanded.decoder == "MULTICH8";
selected_indeces = find(selection)';
architecture = ArchitectureEnum.FORCE_INT_MULTICH8;
torqueMaxAllowed = 5;
torqueMinAllowed = -5;
trimStartIndex = 60;

%% IMPORTANT PARAMETERS

Ts = 0.01;

torqueEndEpsilon = 0.2;

maxSignalLength = 400;

% SET G ESTIMATION TECHNIQUE

% G_estCase=GEstCasesEnum.NO_TORQUE_EDIT;
G_estCase=GEstCasesEnum.TORQUE_END_AT_ZERO;
% G_estCase=GEstCasesEnum.TORQUE_VEL_ESTIMATION;

G_plot_signals = true;

% variable for removing or not torques which min or max values are outside
% the torqueMaxAllowed amd torqueMinAllowed declared in DATA FILTERING section
torqueRangeFiltering = false;

% params for discarding torques signal which decrease after the trimming operation
discardDecreaseStartingTorque = false;
discardDecreaseStartingTorqueIdx = 30;

% outliers data filtering based on standard deviation, 
stdOutlierRemoval = true;

% align the torque signal at max
alignAtMax = false;

% trim signals at start
trimAtStart = true;

% DATA CLEANING
clc
close all

allTrimmedTorque = [];
allTrimmedPos = [];
allTrimmedVel = [];

if G_plot_signals
    figure; 
    hAx1 = subplot(3,1,1);
    title('Torque output from experiments')
    hAx2 = subplot(3,1,2);
    title('Position of the experiments') 
    hAx3 = subplot(3,1,3);
    title('Velocity of the experiments') 
    
    figure; 
    hAx4 = subplot(3,1,1);
    title('FILT: Torque output from experiments')
    hAx5 = subplot(3,1,2);
    title('FILT: Position of the experiments')  
    hAx6 = subplot(3,1,3);
    title('FILT: Velocity of experiments')
    
    
    axes(hAx1)
    hold on
    axes(hAx2)
    hold on
    axes(hAx3)
    hold on
    axes(hAx4)
    hold on
    axes(hAx5)
    hold on
    axes(hAx6)
    hold on
    
end

% Find the torque index of the further max value of the entire set

% Find greater max torque value idx of the various signals 
if alignAtMax
    maxIdx = 0;
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
end

i = 1;
for expId = selected_indeces  %for each experiment
    tagIdx = fixedDataLongTableExpanded(expId,:).tag_idx;
    
    % vector of position
    position = exoRefSplineCells{tagIdx}.thetaE(fixedDataLongTableExpanded(expId,:).start_idx:fixedDataLongTableExpanded(expId,:).end_idx);
    
    % vector of velocity
    velocity = exoRefSplineCells{tagIdx}.dThetaE(fixedDataLongTableExpanded(expId,:).start_idx:fixedDataLongTableExpanded(expId,:).end_idx);
   
    % vector of torque
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
    
    
    augmTorque = [];
    augmPos = [];
    augmVel = [];

    firstValOfTorque = torque(1);
    firstValOfPos = position(1);

    lastValOfTorque = torque(end);
    lastValOfPos = position(end);
    
    % Align signal to the maximum value of torque
    if alignAtMax
        [maxVal, idx] = max(torque);

        pointsToAddStart = maxIdx - idx; 

        % signal where the torque max value is at the same index for all exps
        augmTorque = [ones(pointsToAddStart,1)*firstValOfTorque; torque];
        augmPos = [ones(pointsToAddStart,1)*firstValOfPos; position];
        augmVel = [zeros(pointsToAddStart,1); velocity];
    else
        augmTorque = torque;
        augmPos = position;
        augmVel = velocity;
    end
    
    adjustedTorque = []; 
    adjustedPos = [];
    adjustedVel = [];
    
    if(length(augmTorque) < maxSignalLength)
        pointsToAddAtEnd = maxSignalLength - length(augmTorque);
        adjustedPos = [augmPos; ones(pointsToAddAtEnd,1)*lastValOfPos]; 
        adjustedVel = [augmVel; zeros(pointsToAddAtEnd,1)];
        
        if G_estCase==GEstCasesEnum.NO_TORQUE_EDIT
            adjustedTorque = [augmTorque; ones(pointsToAddAtEnd,1)*lastValOfTorque];
        end
        
        if G_estCase==GEstCasesEnum.TORQUE_END_AT_ZERO || G_estCase==GEstCasesEnum.TORQUE_VEL_ESTIMATION
            adjustedTorque = [augmTorque; zeros(pointsToAddAtEnd,1)];
        end
    else
        adjustedTorque = augmTorque(1:maxSignalLength);
        adjustedPos = augmPos(1:maxSignalLength);
        adjustedVel = augmVel(1:maxSignalLength);
    end
    
    if G_estCase==GEstCasesEnum.TORQUE_END_AT_ZERO || G_estCase==GEstCasesEnum.TORQUE_VEL_ESTIMATION
       adjustedTorque = [adjustedTorque; zeros(maxSignalLength,1)];
    else
        adjustedTorque = [adjustedTorque; ones(maxSignalLength,1)*adjustedTorque(end)];
    end
    adjustedPos = [adjustedPos; ones(maxSignalLength,1)*adjustedPos(end)];    
    adjustedVel = [adjustedVel; ones(maxSignalLength,1)*adjustedVel(end)];    
    
    
    allTrimmedTorque(i,:) = adjustedTorque';
    allTrimmedPos(i,:) = adjustedPos';
    allTrimmedVel(i,:) = adjustedVel';
    
    i = i+1;
    if G_plot_signals
        plot(hAx1,torque);
        plot(hAx2,position);
        plot(hAx3,velocity);
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


% Trim at start, only of align at max is disabled 
if trimAtStart
    supTorque = [];
    supPos = [];
    supVel = [];

    for expId=1:1:size(allTrimmedTorque,1)
        supTorque(expId,:) = allTrimmedTorque(expId, trimStartIndex:end);
        supPos(expId,:) =  allTrimmedPos(expId, trimStartIndex:end);
        supVel(expId,:) =  allTrimmedVel(expId, trimStartIndex:end);
    end
    allTrimmedTorque = supTorque;
    allTrimmedPos = supPos;
    allTrimmedVel = supVel;

    clear supTorque supPos;
end
 
if G_plot_signals 
    for testIdx=1:1:size(allTrimmedTorque,1)
        plot(hAx4, allTrimmedTorque(testIdx, :));
        plot(hAx5, allTrimmedPos(testIdx, :));
        plot(hAx6, allTrimmedVel(testIdx, :));
    end
end
% clear hAx1 hAx2 hAx3 hAx4 hAx5 firstValOfPos firstValOfTorque expId augmPos augmTorque 
% clear adjustedPos adjustedVel adjustedTorque i maxIdx pointsToAddStart pointsToAddAtEnd
% clear torque position

size(allTrimmedPos,1)
%% Model Estimation

clc;

% Number of poles and zeros for the estimated models
G_num_zeros = 0;

if G_estCase == GEstCasesEnum.TORQUE_VEL_ESTIMATION
    G_num_poles = 1;
else
    G_num_poles = 2;
end

% For each esperiment estimate G, O(n) and construct iddata structures
G_iddata = {};
G_sys_est = {};

% option to force estimated model to be stable
opt = tfestOptions('EnforceStability',true,'InitialCondition','estimate');

for expIdx=1:1:size(allTrimmedTorque,1)
    if G_estCase == GEstCasesEnum.TORQUE_VEL_ESTIMATION
        G_est_iddata = iddata(allTrimmedVel(expIdx,:)', allTrimmedTorque(expIdx,:)', Ts); 
    else
        G_est_iddata = iddata(allTrimmedPos(expIdx,:)', allTrimmedTorque(expIdx,:)', Ts); 
    end
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


G_best = zpk(zero(G_bestModel),pole(G_bestModel),1);
[G_best_num, G_best_den] = tfdata(G_best, 'v');


% switch G_estCase
%     case GEstCasesEnum.NO_TORQUE_EDIT
%         G_resultTable(architecture, :) = {char(architecture), G_num_zeros, G_num_poles, size(allTrimmedTorque,1), {G_best}, G_bestModelFit, {G_best_num}, {G_best_den}, {pole(G_best)}};
%     case GEstCasesEnum.TORQUE_END_AT_ZERO
%         G_resultTable_zero(architecture, :) = {char(architecture), G_num_zeros, G_num_poles, size(allTrimmedTorque,1), {G_best}, G_bestModelFit, {G_best_num}, {G_best_den}, {pole(G_best)}};
%     case GEstCasesEnum.TORQUE_VEL_ESTIMATION
%         G_resultTable_vel(architecture, :) = {char(architecture), G_num_zeros, G_num_poles, size(allTrimmedTorque,1), {G_best}, G_bestModelFit, {G_best_num}, {G_best_den}, {pole(G_best)}};
% end

% PLOTS
close all
imageSavePath = "C:\\Users\\nicol\\Desktop\\rpcProject\\myo_tools_testing\\RPC_MN_project\\images\\G_estimation\\";
saveImage = false;

switch G_estCase
    case GEstCasesEnum.NO_TORQUE_EDIT
        imageSavePath = imageSavePath + "noTorqueEdit//";
    case GEstCasesEnum.TORQUE_END_AT_ZERO
        imageSavePath = imageSavePath + "torqueEndAtZero//";
    case GEstCasesEnum.TORQUE_VEL_ESTIMATION
        imageSavePath = imageSavePath + "torqueVelEstimation//";
end

G_plotFunction(allTrimmedTorque, allTrimmedPos, G_bestModelOutput, saveImage, imageSavePath, architecture )