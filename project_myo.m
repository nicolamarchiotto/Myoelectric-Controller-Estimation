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
allData = {};
allSys = {};
allRef = {};
allPos = {};
allContData = {};
allContSys = {};
allSignal = {};
i = 0;

maxSignalLength = 140;

for testIdx=selected_indexes  %for each arm motion (Forward or Backward) estimate the controller
    i = i + 1;
    tagIdx = fixedDataTableExpanded(testIdx,:).tag_idx;
    target = fixedDataTableExpanded(testIdx,:).target_angle;
    test_position_v = exoRefSplineCells{tagIdx}.thetaE(fixedDataTableExpanded(testIdx,:).start_idx:fixedDataTableExpanded(testIdx,:).end_idx);
    
    % Reference voltage
    signal = exoRefSplineCells{tagIdx}.reference(fixedDataTableExpanded(testIdx,:).start_idx:fixedDataTableExpanded(testIdx,:).end_idx);
    
    impIdx = 1;

    % Data cleaning, remove data of test before actual motion
    while (abs(test_position_v(1) - test_position_v(impIdx+1))<0.005) 
        %idea: cut everything until we drop the value over 0.01, i.e. the
        %movement has begun
       impIdx = impIdx + 1; 
    end
    cutPos = test_position_v(impIdx:impIdx+maxSignalLength);
    
    % Vector of reference values
    ref = ones(length(cutPos),1)*target;
    allPos{i} = cutPos;
    allRef{i} = ref;

    % Vector of error
    err = ref-cutPos;
    
    % Reference signal cut to correspond with cutPos length
    cutSignal = signal(impIdx:impIdx+maxSignalLength);

    % Structure needed by System Id Toolbox iddata, output=cutSignal,
    % input=err, brain controller estimation 
    cont_data = iddata(cutSignal, err, Ts); 

    %CONTROLLER ESTIMATION C
    cont_sys = tfest(cont_data,1,1); %1 pole(high freq hopefully), 1 zero
    
    % Save estimated controller for the single test
    allContData{i} = cont_data;
    allContSys{i} = cont_sys;
    allSignal{i} = cutSignal;
    
    %MODEL ESTIMATION C*G
    full_data = iddata(cutPos, ref, Ts);
    full_sys = tfest(full_data, 3, 1);
    
    % Save estimated model for the single test
    allData{i} = full_data;
    allSys{i} = full_sys;
end

%% PLOT: see signals with removed delay, plotted pos of imported data
figure;
hold on;
for i = 1:length(allPos)
    plot(allPos{i});
end
%% TESTING: Search best estimatred model C*G, carefull O(n^2)
bestModel = 1;
bestModelFit = 0;
outliers = [];
for i = 1:length(allSys)
    % Compare the estimated model output of the single test with the real output from data
    modelFit = 0;
    for j = 1:length(allData)    
        [y,fit] = compare(allData{j}, allSys{i}); 
        modelFit = modelFit + fit;
    end
    modelFit = modelFit/length(allData);
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

pos = allPos{bestModel};
figure;
plot(0:length(pos)-1, pos);
figure;
[y,fit] = compare(allData{bestModel}, allSys{bestModel});
y1 = cell2mat(get(y).OutputData);
plot(0:length(y1)-1,y1);

%% search best controller
[bestCont, bestContFit] = bestModelFinder(allContSys, allContData);

%% PLOT: See best model on all positions
figure(10);
subplot(1,2,1);
title('Real position');
subplot(1,2,2);
title('Estimated position');
for i=1:length(allData)
    if ismember(i, outliers)
        continue;
    end
    pos = allPos{i};
    subplot(1,2,1);
    hold on;
    plot(0:length(pos)-1, pos);
    [y,fit] = compare(allData{i}, allSys{bestModel});
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
for i=1:length(allContData)
    signal = allSignal{i};
    subplot(1,2,1);
    hold on;
    plot(0:length(signal)-1, signal);
    [y,fit] = compare(allContData{i}, allContSys{bestCont});
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
for i=1:length(allSignal)
    %almost end in zero
    if abs(allSignal{i}(end))>0.002
        continue;
    end
    newSignals{j} = allSignal{i};
    plot(newSignals{j});
    realPos = allPos{i};
    err = ref-realPos;
    cont_data = iddata(newSignals, err, Ts);%for controller estimation
    sys11{j} = tfest(cont_data,1,1);
    sys21{j} = tfest(cont_data,2,1);
    sys22{j} = tfest(cont_data,2,2);
    newSigData{j} = cont_data;
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