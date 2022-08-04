%% Data import
clc
clear
addpath(genpath('..\\myo_tools_testing'))
logs_eval_gains_2019_10_09_10_18__2020_06_26_30 %inside testing1\settings folder
batchMultProcessingTest1
batchMultElaborationTest1Long
close all

%% BatchKinematicAnalysis

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
    actualFixedDataTableExpanded = fixedDataLongTableExpanded;
else
    actualFixedDataTableExpanded = fixedDataTableExpanded;
end

%% DATA FILTERING: F=forward direction of motion
selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "COMP" & actualFixedDataTableExpanded.decoder == "NONE";
% plotKinematicAnalysis(logDir,elabDir,'COMP_F_POS1_',exoRefSplineCells,Ts_emg,time_limit_comp,actualFixedDataTableExpanded,selection,normalize,reverse_dir,plot_mean,plot_median,save_plots)
selected_indeces = find(selection)';

%% selection with given starting angle and targetAngle
startAng = -1.323240000000000;
targetAng = -2.094400000000000;

selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "COMP" & actualFixedDataTableExpanded.decoder == "NONE" & abs(actualFixedDataTableExpanded.start_angle - startAng) < 0.1 & actualFixedDataTableExpanded.target_angle==targetAng;
% plotKinematicAnalysis(logDir,elabDir,'COMP_F_POS1_',exoRefSplineCells,Ts_emg,time_limit_comp,actualFixedDataTableExpanded,selection,normalize,reverse_dir,plot_mean,plot_median,save_plots)
selected_indeces = find(selection)';

%% IMPORTANT PARAMETERS

Ts = 0.01;

posEpsilon = 0.005;
torqueEpsilon = 0.2;

maxSignalLength = 140;

C_plot_signals = true;
W_plot_signals = true;
G_plot_signals = true;
%% DATA CLEANING

allTrimmedTorque = {};
allTrimmedPosError = {};


allTrimmedPos = {};
allTrimmedRefPos={};

i = 1;

if C_plot_signals
    figure; 
    hAx1 = subplot(2,1,1);
    title('Position error of the experiments') 
    hAx2 = subplot(2,1,2);
    title('Torque output from experiments')
    figure; 
    hAx3 = subplot(2,1,1);
    title('FILT: Position error of the experiments') 
    hAx4 = subplot(2,1,2);
    title('FILT: Torque output from experiments')
    axes(hAx1)
    hold on
    axes(hAx2)
    hold on
    axes(hAx3)
    hold on
    axes(hAx4)
    hold on
end

if W_plot_signals
    figure; 
    hAx5 = subplot(2,1,1);
    title('Ref position of the experiments') 
    hAx6 = subplot(2,1,2);
    title('Position output from experiments')
    figure; 
    hAx7 = subplot(2,1,1);
    title('FILT: Ref position error of the experiments') 
    hAx8 = subplot(2,1,2);
    title('FILT: Position output from experiments')
    axes(hAx5)
    hold on
    axes(hAx6)
    hold on
    axes(hAx7)
    hold on
    axes(hAx8)
    hold on
end

for expId = selected_indeces  %for each experiment
    tagIdx = fixedDataLongTableExpanded(expId,:).tag_idx;
    
    %reference position value
    referencePosVal = fixedDataLongTableExpanded(expId,:).target_angle;
    
    %vector of positions of the experiment
    position = exoRefSplineCells{tagIdx}.thetaE(fixedDataLongTableExpanded(expId,:).start_idx:fixedDataLongTableExpanded(expId,:).end_idx);
    
    % vector of reference torques
    torque = exoRefSplineCells{tagIdx}.reference(fixedDataLongTableExpanded(expId,:).start_idx:fixedDataLongTableExpanded(expId,:).end_idx);

    trimIdx = 1;

    % Data cleaning, remove data of test before actual motion
    while (abs(position(1) - position(trimIdx+1)) < posEpsilon)   
       trimIdx = trimIdx + 1; 
    end
    
    trimmedPos = position(trimIdx:trimIdx+maxSignalLength);
    
    % Vector of reference values
    trimmedRefPos = ones(length(trimmedPos),1)*referencePosVal;

    % Vector of position error
    trimmedPosErr = trimmedRefPos-trimmedPos;
    
    % Reference signal cut to correspond with trimmedPos length
    trimmedTorque = torque(trimIdx:trimIdx+maxSignalLength);
    
    % outliers torque elimination
    if abs(trimmedTorque(end)) > torqueEpsilon
        continue;
    end
    
    % Take the last value of the filtered data and lengthening the signals
    for k=1:1:2*maxSignalLength
       trimmedTorque = [trimmedTorque; trimmedTorque(end)];
       trimmedPosErr = [trimmedPosErr; trimmedPosErr(end)];
       trimmedPos = [trimmedPos; trimmedPos(end)];
       trimmedRefPos = [trimmedRefPos; trimmedRefPos(end)];
    end
    
    if C_plot_signals 
        plot(hAx1,ones(length(position),1)*referencePosVal-position);
        plot(hAx2,torque);
        plot(hAx3,trimmedPosErr);
        plot(hAx4,trimmedTorque);
    end
    
    if W_plot_signals
        plot(hAx5,ones(length(position),1)*referencePosVal);
        plot(hAx6,position);
        plot(hAx7,trimmedRefPos);
        plot(hAx8,trimmedPos);
    end
    
    allTrimmedTorque{i} = trimmedTorque;
    allTrimmedPosError{i} = trimmedPosErr;
    allTrimmedPos{i} = trimmedPos;
    allTrimmedRefPos{i} = trimmedRefPos;
    
    i = i + 1;
end

clear expId i k;

%% Construct structures for estimation and testing

clc;

% Number of poles and zeros for the estimated models
C_num_zeros = 2;
C_num_poles = 2;

G_num_zeros = 0;
G_num_poles = 2;

W_num_zeros = C_num_zeros;
W_num_poles = C_num_poles + 2;

% construct iddata structures
C_iddata = {};
W_iddata = {};
G_iddata = {};

for testIdx=1:1:length(allTrimmedTorque)
    C_iddata{testIdx} = iddata(allTrimmedTorque{testIdx}, allTrimmedRefPos{testIdx}, Ts); 
    W_iddata{testIdx} = iddata(allTrimmedPos{testIdx}, allTrimmedRefPos{testIdx}, Ts); 
    G_iddata{testIdx} = iddata(allTrimmedPos{testIdx}, allTrimmedTorque{testIdx}, Ts); 

end

clear testIdx;

% For each esperiment estimate C and W and G, 3*O(n)

C_sys_est = {};
W_sys_est = {};
G_sys_est = {};
for expIdx=1:1:length(allTrimmedRefPos)
    C_est_iddata = iddata(allTrimmedTorque{expIdx}, allTrimmedRefPos{expIdx}, Ts); 
    C_sys_est{expIdx} = tfest(C_est_iddata, C_num_poles, C_num_zeros); 
    
    W_est_iddata = iddata(allTrimmedPos{expIdx}, allTrimmedRefPos{expIdx}, Ts); 
    W_sys_est{expIdx} = tfest(W_est_iddata, W_num_poles, W_num_zeros); 
    
    G_est_iddata = iddata(allTrimmedPos{expIdx}, allTrimmedTorque{expIdx}, Ts); 
    G_sys_est{expIdx} = tfest(G_est_iddata, G_num_poles, G_num_zeros); 
end

%initialization of variables for plot function
C_bestModelOutput = [];
W_bestModelOutput = [];
G_bestModelOutput = [];
clear expIdx;
%% Find the best C testing on all experiments, O(n^2)
%
%
% CAREFULL, each bestModelFinder is O(n^2)
%
%

[C_bestModel, C_bestModelFit, C_bestModelOutput] = bestModelFinder(C_sys_est, C_iddata);

% Find the best W testing on all experiments, O(n^2)
%
%
% CAREFULL, each bestModelFinder is O(n^2)
%
%
[W_bestModel, W_bestModelFit, W_bestModelOutput] = bestModelFinder(W_sys_est, W_iddata);

% Find the best G testing on all experiments, O(n^2)
%
%
% CAREFULL, each bestModelFinder is O(n^2)
%
%
[G_bestModel, G_bestModelFit, G_bestModelOutput] = bestModelFinder(G_sys_est, G_iddata);

%% RESULTS
clc 
fprintf('\nC: Best model fit from single experiment estimation: %.3f\n', C_bestModelFit);

fprintf('\nW: Best model fit from single experiment estimation: %.3f\n', W_bestModelFit);

fprintf('\nG: Best model fit from single experiment estimation: %.3f\n', G_bestModelFit);

%% Simulink, estimated controller from position error and torque
save bestController.mat C_bestModel
load bestController.mat
s = tf('s');
J = 0.068;
D = 0.1; %da 0.1 a 0.001
beta_pos = 4;
beta_force = 3;
beta_force_int = 4;
Kp = 140;
Kd = 2;
Ki = 0;

%G = 1/(J*s^2 + D*s);
G = zpk(G_bestModel);
% Estimated controller from posErr and torque
C_best =  zpk(C_bestModel);

% Extracted controller from estimated W
C_from_Wbest = zpk(getC_from_G_and_W(G, W_bestModel));
    
%% PLOTS

close all

% Controller C
if C_plot_signals
    for j=1:1:ceil(length(C_iddata)/10)
        figure(j)
        sgtitle('Torque of best estimated C'); 
        for i=1:1:10
            testIdx=(j-1)*10+i;
            if(length(allTrimmedTorque)<testIdx)
                continue
            end
            subplot(2,5,i)
            titleStr=['Test ',num2str(testIdx)]; 
            hold on;
            title(titleStr) 
            plot(0:length(allTrimmedTorque{testIdx})-1,allTrimmedTorque{testIdx});
            plot(0:length(C_bestModelOutput{testIdx})-1,C_bestModelOutput{testIdx});
            legend('Location','southoutside')
            legend('testing', 'estimated model output')
        end
    end
end

k=j;
clear j;

%% Whole model W
k = 0;
if W_plot_signals
     for j=1:1:ceil(length(W_iddata)/10)
        figure(j+k)
        sgtitle('Position best estimated W'); 
        for i=1:1:10
            testIdx=(j-1)*10+i;
            if(length(allTrimmedPos)<testIdx)
                continue
            end
            subplot(2,5,i)
            titleStr=['Test ',num2str(testIdx)]; 
            hold on;
            title(titleStr) 
            plot(0:length(allTrimmedPos{testIdx})-1,allTrimmedPos{testIdx});
            plot(0:length(W_bestModelOutput{testIdx})-1,W_bestModelOutput{testIdx});
            legend('Location','southoutside')
            legend('testing', 'estimated model output')
        end
     end
end
clear i j k;

%% Motor G
k = 0;
h = 0;
if G_plot_signals
     for j=1:1:ceil(length(G_iddata)/10)
        figure(j+k+h)
        sgtitle('Position best estimated G'); 
        for i=1:1:10
            testIdx=(j-1)*10+i;
            if(length(allTrimmedPos)<testIdx)
                continue
            end
            subplot(2,5,i)
            titleStr=['Test ',num2str(testIdx)]; 
            hold on;
            title(titleStr) 
            plot(0:length(allTrimmedPos{testIdx})-1,allTrimmedPos{testIdx});
            plot(0:length(G_bestModelOutput{testIdx})-1,G_bestModelOutput{testIdx});
            legend('Location','southoutside')
            legend('testing', 'estimated model output')
        end
     end
end
clear i j k h;
%%
clc
plot_C = true;
plot_W = true;
plot_G = true;

plotFunction(plot_C, plot_W, plot_G, allTrimmedTorque, allTrimmedPos, C_bestModelOutput, W_bestModelOutput, G_bestModelOutput)
