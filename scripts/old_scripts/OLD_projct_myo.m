%%
addpath(genpath('../myo_tools_testing'))
logs_eval_2019_10_10 %inside settings folder
batchProcessingTest1
batchElaborationTest1
close all

%% batchKinematicAnalysis
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

% Use table & exo data 

%%
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "COMP" & actualFixedDataTableExpanded.decoder == "NONE";
% plotKinematicAnalysis(logDir,elabDir,'COMP_F_POS1_',exoRefSplineCells,Ts_emg,time_limit_comp,actualFixedDataTableExpanded,selection,normalize,reverse_dir,plot_mean,plot_median,save_plots)
selected_indexes = find(selection)';

%% collect allData and models
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
for testIdx=selected_indexes
    i = i + 1;
    tagIdx = fixedDataTableExpanded(testIdx,:).tag_idx;
    target = fixedDataTableExpanded(testIdx,:).target_angle;
    test_position_v = exoRefSplineCells{tagIdx}.thetaE(fixedDataTableExpanded(testIdx,:).start_idx:fixedDataTableExpanded(testIdx,:).end_idx);
    signal = exoRefSplineCells{tagIdx}.reference(fixedDataTableExpanded(testIdx,:).start_idx:fixedDataTableExpanded(testIdx,:).end_idx);
    
    impIdx = 1;
    while (abs(test_position_v(1) - test_position_v(impIdx+1))<0.005) 
        %idea: cut everything until we drop the value over 0.01, i.e. the
        %movement has begun
       impIdx = impIdx + 1; 
    end
    cutPos = test_position_v(impIdx:impIdx+maxSignalLength);
    ref = ones(length(cutPos),1)*target;
    allPos{i} = cutPos;
    allRef{i} = ref;
    err = ref-cutPos;
    
    cutSignal = signal(impIdx:impIdx+maxSignalLength);
    cont_data = iddata(cutSignal, err, Ts);%for controller estimation
    cont_sys = tfest(cont_data,1,1); %1 pole(high freq hopefully), 1 zero
    allContData{i} = cont_data;
    allContSys{i} = cont_sys;
    allSignal{i} = cutSignal;
    
    full_data = iddata(cutPos, ref, Ts);
    full_sys = tfest(full_data, 3, 1);
    allData{i} = full_data;
    allSys{i} = full_sys;
end
%% see signals with removed delay
figure;
hold on;
for i = 1:length(allPos)
    plot(allPos{i});
end
%% search best model
bestModel = 1;
bestModelFit = 0;
outliers = [];
for i = 1:length(allSys)
    %for each estimated model
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

%% see best model on all positions
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


%% see best controller on all torques
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

%% another approach, calculate mean signal and estimate on that
medPos = median(cell2mat(allPos),2);
%plot(medPos)
% move every signal up or down to the median value
newPos = {};

for i=1:length(allPos)
   diff = mean(allPos{i}-medPos);
   newPos{i} = allPos{i} - diff;
   figure(100);
   hold on;
   plot(newPos{i});
   figure(200);
   hold on;
   plot(allPos{i});
end
figure(100);
plot(medPos, 'LineWidth', 3);
figure(200);
plot(medPos, 'LineWidth', 3);

%% median model
full_data = iddata(medPos, ref, Ts);
medSys = tfest(full_data, 3, 1);
totalFit = 0;
for i=1:length(allData)
    figure(80);
    title('Real position');
    hold on;
    plot(0:length(allPos{i})-1, allPos{i});
    [y,fit] = compare(allData{i}, medSys);
    totalFit = totalFit + fit;
    y1 = cell2mat(get(y).OutputData);
    figure(90);
    title('Estimated position');
    hold on;
    plot(0:length(y1)-1,y1);
end
totalFit/length(allData)

%%
s = tf('s');
j = 0.068;
d = 0.01;
G = 1/(j*s^2 + d*s);
C = getC_from_G_and_W(G, allSys{bestModel});
newC = zpk(minreal(C,0.5));

%% remove outliers from myo signals
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
    if abs(allSignal{i}(end))>0.2
        continue;
    end
    newSignals{j} = allSignal{i};
    plot(newSignals{j});
    realPos = allPos{i};
    err = ref-realPos;
    cont_data = iddata(cutSignal, err, Ts);%for controller estimation
    sys11{j} = tfest(cont_data,1,1);
    sys21{j} = tfest(cont_data,2,1);
    sys22{j} = tfest(cont_data,2,2);
    newSigData{j} = cont_data;
    j = j+1;
end

%% result
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

%% best controller model
bestCtrlModel = sys22{best22idx};
save bestController.mat bestCtrlModel
% forw = bestCtrlModel*G;
% H = minreal(forw/(1+forw), 0.1);
%%
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