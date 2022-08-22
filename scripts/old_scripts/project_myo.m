%% Data import
clc
addpath(genpath('..\\myo_tools_testing'))
logs_eval_2019_10_10 %inside testing1\settings folder
% logs_eval_2019_10_09
batchProcessingTest1
batchElaborationTest1
close all

% BatchKinematicAnalysis

% Requires:
%  - calling batchElaborationTest1 or batchAdjustElaborationTest1Long 
elabDir = ".\\elaborated_kinematics\\";
save_plots = 0;

time_limit_comp = 4.0;
time_limit = 10.0;
time_eps = 0.1;
normalize = 1;
reverse_dir = 1;
plot_mean = 1;
plot_median = 1;
%%
if fixed_data_long_format % both 1&2
    actualFixedDataTableExpanded = fixedDataLongTableExpanded_adjusted;
else
    actualFixedDataTableExpanded = fixedDataTableExpanded;
end

%% DATA FILTERING: F=forward direction of motion
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "COMP" & actualFixedDataTableExpanded.decoder == "NONE";
% plotKinematicAnalysis(logDir,elabDir,'COMP_F_POS1_',exoRefSplineCells,Ts_emg,time_limit_comp,actualFixedDataTableExpanded,selection,normalize,reverse_dir,plot_mean,plot_median,save_plots)
selected_indexes = find(selection)';

%% CONTROLLER/MODEL IDENTIFICATION
% For each experiment (single arm motion), estimate the controller C and
% model C*G

Ts = 0.01;

allTrimmedPosRef = {};
allTrimmedPos = {};

allTrimmedTorque = {};
allTrimmedPosError = {};


all_W_iddata = {};
all_W_est = {};

all_C_iddata = {};
all_C_est = {};

i = 0;

maxSignalLength = 140;

for testIdx=selected_indexes  %for each arm motion = 1 experiment, estimate controller and model
    i = i + 1;
    tagIdx = fixedDataTableExpanded(testIdx,:).tag_idx;
    %reference value
    referencePosVal = fixedDataTableExpanded(testIdx,:).target_angle;
    %vector of positions of the experiment
    position = exoRefSplineCells{tagIdx}.thetaE(fixedDataTableExpanded(testIdx,:).start_idx:fixedDataTableExpanded(testIdx,:).end_idx);
    
    % Reference torque
    torque = exoRefSplineCells{tagIdx}.reference(fixedDataTableExpanded(testIdx,:).start_idx:fixedDataTableExpanded(testIdx,:).end_idx);
    
    impIdx = 1;

    % Data cleaning, remove data of test before actual motion
    % idea: cut everything until we drop the value over 0.01, i.e. the movement has begun
    while (abs(position(1) - position(impIdx+1))<0.005)   
       impIdx = impIdx + 1; 
    end
    
    trimmedPos = position(impIdx:impIdx+maxSignalLength);
    
    % Vector of reference values
    trimmedRefPos = ones(length(trimmedPos),1)*referencePosVal;

    % Vector of position error
    trimmedPosErr = trimmedRefPos-trimmedPos;
    
    % Reference signal cut to correspond with cutPos length
    trimmedTorque = torque(impIdx:impIdx+maxSignalLength);

    % Structure needed by System Id Toolbox iddata(y,u,Ts)
    % For model estimation, tfest(iddata, num_poles, num_zeros)
    
    % CONTROLLER ESTIMATION C
    C_iddata = iddata(trimmedTorque, trimmedPosErr, Ts); 
    C_est = tfest(C_iddata,1,1); %1 pole(high freq hopefully), 1 zero
    
    % Save estimated controller
    all_C_iddata{i} = C_iddata;
    all_C_est{i} = C_est;
    allTrimmedTorque{i} = trimmedTorque;
    allTrimmedPosError{i} = trimmedPosErr;
    
    % MODEL ESTIMATION W = CG/1-CG
    
    W_iddata = iddata(trimmedPos, trimmedRefPos, Ts);
    W_est = tfest(W_iddata, 3, 1);
    
    % Save estimated model
    all_W_iddata{i} = W_iddata;
    all_W_est{i} = W_est;
    allTrimmedPos{i} = trimmedPos;
    allTrimmedPosRef{i} = trimmedRefPos;
end

%% PLOT: pose with removed delay
figure;
hold on;
for i = 1:length(allTrimmedPos)
    plot(allTrimmedPos{i});
end
%% PLOT: torque with removed delay
figure;
hold on;
for i = 1:length(allTrimmedTorque)
    plot(allTrimmedTorque{i});
end

%% TESTING: Search best estimatred model C*G, carefull O(n^2)
bestModel = 1;
bestModelFit = 0;
outliers = [];
for i = 1:length(all_W_est)
    % Compare the estimated model output of the single test with the real output from data
    modelFit = 0;
    for j = 1:length(all_W_iddata)    
        [y,fit] = compare(all_W_iddata{j}, all_W_est{i}); 
        modelFit = modelFit + fit;
    end
    modelFit = modelFit/length(all_W_iddata);
    if modelFit > bestModelFit
        bestModelFit = modelFit;
        bestModel = i;
    end
    if modelFit < 70
        outliers = [outliers i];
    end
    fprintf('Model %d mean fit: %.2f\n', i, modelFit);
end
fprintf('\nBest model %d with %.2f mean fit\n', bestModel, bestModelFit);
%legend

pos = allTrimmedPos{bestModel};
figure;
plot(0:length(pos)-1, pos);
figure;
[y,fit] = compare(all_W_iddata{bestModel}, all_W_est{bestModel});
y1 = cell2mat(get(y).OutputData);
plot(0:length(y1)-1,y1);

%% search best controller
[bestCont, bestContFit] = bestModelFinder(all_C_est, all_C_iddata);

%% PLOT: See best model on all positions
figure(10);
subplot(1,2,1);
title('Real position');
subplot(1,2,2);
title('Estimated position');
for i=1:length(all_W_iddata)
    if ismember(i, outliers)
        continue;
    end
    pos = allTrimmedPos{i};
    subplot(1,2,1);
    hold on;
    plot(0:length(pos)-1, pos);
    [y,fit] = compare(all_W_iddata{i}, all_W_est{bestModel});
    y1 = cell2mat(get(y).OutputData);
    subplot(1,2,2);
    hold on;
    plot(0:length(y1)-1,y1);
end

%% PLOT: See best controller on all torques
figure(30);
subplot(1,2,1);
title('Real torque');
subplot(1,2,2);
title('Estimated torque');
for i=1:length(allTrimmedTorque)
    torque = all_C_y{i};
    subplot(1,2,1);
    hold on;
    plot(0:length(torque)-1, torque);
    [y,fit] = compare(allTrimmedTorque{i}, all_C_est{bestCont});
    y1 = cell2mat(get(y).OutputData);
    subplot(1,2,2);
    hold on;
    plot(0:length(y1)-1,y1);
end



%% GRAV COMP: Remove outliers from myo signals, done for case of gravity compensation architecture, NOT appliable for other architectures
newSignals = {};
sys11 = {};
sys21 = {};
sys22 = {};
newSigData = {};
j = 1;
figure(111);
hold on;
title('Myo signals that end in 0');
for i=1:length(allTrimmedTorque)
    %almost end in zero
    if abs(allTrimmedTorque{i}(end))>0.002
        continue;
    end
    newSignals{j} = allTrimmedTorque{i};
    plot(newSignals{j});
    realPos = allTrimmedPos{i};
    trimmedPosErr = trimmedRefPos-realPos;
    C_iddata = iddata(newSignals, trimmedPosErr, Ts);%for controller estimation
    sys11{j} = tfest(C_iddata,1,1);
    sys21{j} = tfest(C_iddata,2,1);
    sys22{j} = tfest(C_iddata,2,2);
    newSigData{j} = C_iddata;
    j = j+1;
end

%% GRAV COMP: result
[best11idx, best11fit] = bestModelFinder(sys11, newSigData);
[best21idx, best21fit] = bestModelFinder(sys21, newSigData);
[best22idx, best22fit] = bestModelFinder(sys22, newSigData);

%     signal = models{bestModelIdx};
%     figure;
%     plot(0:length(signal)-1, signal);
%     figure;
%     [y,fit] = compare(data{bestCont}, models{bestCont});
%     y1 = cell2mat(get(y).OutputData);
%     plot(0:length(y1)-1,y1);

%% GRAV COMP: best controller model
bestCtrlModel = sys22{best22idx};
save bestController.mat bestCtrlModel
% forw = bestCtrlModel*G;
% H = minreal(forw/(1+forw), 0.1);
%% GRAV COMP
load bestController.mat
s = tf('s');
j = 0.068;
d = 0.01; %da 0.1 a 0.001
beta_pos = 4;
beta_force = 3;
beta_force_int = 4;
Kp = 140;
Kd = 2;
Ki = 0;