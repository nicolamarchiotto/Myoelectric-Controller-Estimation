%% Controller estimation script with concatenated experiments
%
% Script which confronts the controller estimation using different methods
%
% Preliminary operations:
% 
% Position data cleaning:
% Remove experiment data before actual motion, remove position data until a
% variation of posEpsilon is registered, adjust lenght of other signals
% accordingly, the maximum signal lenght allowed is 140 samples
%
% Torque outliers removal:
% Since all the architectures are based on gravity compensation, is
% reasonable to assume that all the torques signal must be end near the 0
% value. All experiments for which the last torque value is beyond a
% certain threshold torqueEpsilon is consider an outlier and discarded
%
% Divide the experiments in 2 sets, 80% are the estimation set, 20% are the
% testing set
%
% Perform the following estimation methods:
%
% 1) Controller estimation with single experiments 
% Take the estimation set, for each experiment estimate the controller, test each
% controller with the experiment from the testing set and normalize the fitting score, take
% the best fitting score
%
% 2) Controller estimation with concatenation of experiments
% Concat all the esperiments of the estimation set, with this new data 
% estimate the controller, test the estimated controller with the testing
% set experiments of the testing set, normalize the fitting score
%
% 3) Same as method 2 but flip the estimation data to concatenate 1 every 2 experiment,
% this method was tried to partially remove the discontinuities in the
% concatenated position signal
%
% Compare  the results of the methods 
%
%% CONCLUSIONS
%
% The method 3 gave bad results, thus it is not taken in cosideration
%
% Is not always clear which method between 1 and 2 is the best, these
% methods are strongly influenced by the estimatinon and testing experiments which
% are randomly choosen
%
% Using logs_eval_2019_10_10 only 5 experiments respected the filtering
% criterias, 5 were used for the estimation and 2 for testing,
%
% A BIGGER SET OF CLEANED DATA COULD IMPROVE BY A LOT THE CONTROLLER
% ESTIMATION
%
% The best fit was recorded with method 2: bestControllerFit = 83.092
% The bestControllerFit with Simone Cremasco procedure was 83.7869
%
%% DATA IMPORT

clc
clear
addpath(genpath('..\\myo_tools_testing'))
%inside testing1\settings folder
% logs_eval_2019_10_09 % DOES NOT WORK, problems in processingFncTest1 
logs_eval_2019_10_10 % works, bad torques data 
% logs_eval_2019_10_18 % works, only one testing sample
batchProcessingTest1
batchElaborationTest1
close all

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

%% DATA FILTERING: Select experimets from logs data

selection = actualFixedDataTableExpanded.direction == "F" & actualFixedDataTableExpanded.target == "POS1" & actualFixedDataTableExpanded.controller == "COMP" & actualFixedDataTableExpanded.decoder == "NONE";
% plotKinematicAnalysis(logDir,elabDir,'COMP_F_POS1_',exoRefSplineCells,Ts_emg,time_limit_comp,actualFixedDataTableExpanded,selection,normalize,reverse_dir,plot_mean,plot_median,save_plots)
selected_indeces = find(selection)';

%% IMPORTANT PARAMETERS

Ts = 0.01;

posEpsilon=0.005;
torqueEpsilon=0.2;

maxSignalLength = 140;

C_plot_signals=true;
W_plot_signals=false;

% variable which if set to true enables method 3
flip_1_every_2=false;

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
    tagIdx = fixedDataTableExpanded(expId,:).tag_idx;
    
    %reference position value
    referencePosVal = fixedDataTableExpanded(expId,:).target_angle;
    
    %vector of positions of the experiment
    position = exoRefSplineCells{tagIdx}.thetaE(fixedDataTableExpanded(expId,:).start_idx:fixedDataTableExpanded(expId,:).end_idx);
    
    % vector of reference torques
    torque = exoRefSplineCells{tagIdx}.reference(fixedDataTableExpanded(expId,:).start_idx:fixedDataTableExpanded(expId,:).end_idx);

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
    if abs(trimmedTorque(end))>torqueEpsilon
        continue;
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

clear expId i;

%% Construct estimation and testing set

%control estimation
estimationPosError={};
estimationTorque={};

testingPosError={};
testingTorque={};

%whole model estimation 
estimationRefPos={};
estimationPos={};

testingRefPos={};
testingPos={};

% Choosing randomly the testing experiments
numOfTest=ceil(length(allTrimmedTorque)*0.2);

% vector containig the testing indeces, randperm guarantees distinct
% indeces
testing_indeces=randperm(length(allTrimmedTorque),numOfTest);

%support variables used to not modify "all" data structures
deprecatedTrimmedPosError=allTrimmedPosError;
deprecatedTrimmedTorque=allTrimmedTorque;

deprecatedTrimmedRefPos=allTrimmedRefPos;
deprecatedTrimmedPos=allTrimmedPos;

i=1;
for idx = testing_indeces  %for each experiment    
    
    % controller
    testingPosError{i} = allTrimmedPosError{idx};
    testingTorque{i} = allTrimmedTorque{idx};
    
    deprecatedTrimmedPosError{idx}=[];
    deprecatedTrimmedTorque{idx}=[];
    
    % whole model
    testingRefPos{i} = allTrimmedRefPos{idx};
    testingPos{i} = allTrimmedPos{idx};
    
    deprecatedTrimmedRefPos{idx}=[];
    deprecatedTrimmedPos{idx}=[];
    
    
    i=i+1;
end

clear i idx;

i=1;

% For method 3 set flip_1_every_2 to true in IMPORTANT PARAMETERS paragraph
for idx = 1:1:length(deprecatedTrimmedPosError)
   if(isempty(deprecatedTrimmedPosError{idx}))
       continue;
   end
   if flip_1_every_2
        % method 3
       if mod(i,2)==0
           estimationPosError{i}=flip(deprecatedTrimmedPosError{idx});
           estimationTorque{i}=flip(deprecatedTrimmedTorque{idx});
           
           estimationRefPos{i} = flip(deprecatedTrimmedRefPos{idx});
           estimationPos{i} = flip(deprecatedTrimmedPos{idx}); 
       else 
           estimationPosError{i}=deprecatedTrimmedPosError{idx};
           estimationTorque{i}=deprecatedTrimmedTorque{idx};   
           
           estimationRefPos{i} = deprecatedTrimmedRefPos{idx};
           estimationPos{i} = deprecatedTrimmedPos{idx}; 
       end
   else
       % controller
       estimationPosError{i} = deprecatedTrimmedPosError{idx};
       estimationTorque{i} = deprecatedTrimmedTorque{idx}; 
       
       % whole model
       estimationRefPos{i} = deprecatedTrimmedRefPos{idx};
       estimationPos{i} = deprecatedTrimmedPos{idx}; 

   end
   i=i+1;
end

clear i idx;
clear deprecatedTrimmedPosError deprecatedTrimmedTorque;
clear deprecatedTrimmedRefPos deprecatedTrimmedPos; 

% Number of poles and zeros for the controller
C_num_zeros = 2;
C_num_poles = 2;

W_num_zeros = C_num_zeros;
W_num_poles = C_num_poles + 2;

% construct iddata structure for testing
C_testing_iddata={};
W_testing_iddata={};

for testIdx=1:1:length(testingTorque)
    C_testing_iddata{testIdx} = iddata(testingTorque{testIdx}, testingPosError{testIdx}, Ts); 
    W_testing_iddata{testIdx} = iddata(testingPos{testIdx}, testingRefPos{testIdx}, Ts); 
end

clear testIdx;

% Method 1) Model estimation using the single experiments
C_sys_est={};
W_sys_est={};
for expIdx=1:1:length(estimationPosError)
    C_est_iddata = iddata(estimationTorque{expIdx}, estimationPosError{expIdx}, Ts); 
    C_sys_est{expIdx} = tfest(C_est_iddata, C_num_poles, C_num_zeros); 
    
    W_est_iddata = iddata(estimationRefPos{expIdx}, estimationPos{expIdx}, Ts); 
    W_sys_est{expIdx} = tfest(W_est_iddata, W_num_poles, W_num_zeros); 
end

clear expIdx;

[C_bestModelSingle, C_bestModelFitSingle, C_bestModelOutputSingle] = bestModelFinder(C_sys_est, C_testing_iddata);

[W_bestModelSingle, W_bestModelFitSingle, W_bestModelOutputSingle] = bestModelFinder(W_sys_est, W_testing_iddata);

% Method 2) Controller estimation using the 80% of experiments
estPosError80 = [];
estTorque80 = [];
C_sys_est_80 = {};

estPos80 = [];
estRefPos80 = [];
W_sys_est_80 = {};

for expIdx=1:1:length(estimationPosError)
    estPosError80 = [estPosError80; estimationPosError{expIdx}];
    estTorque80 = [estTorque80; estimationTorque{expIdx}];
    
    estPos80 = [estPos80; estimationPos{expIdx}];
    estRefPos80 = [estRefPos80; estimationRefPos{expIdx}];
end

if C_plot_signals
    figure; 
    hAx1 = subplot(2,1,1);
    title('CONCATENATED Estimation Position error') 
    hAx2 = subplot(2,1,2);
    title('CONCATENATED Estimation Torque')
    axes(hAx1)
    hold on
    axes(hAx2)
    hold on
    plot(hAx1,estPosError80);
    plot(hAx2,estTorque80);
end

if W_plot_signals
    figure; 
    hAx1 = subplot(2,1,1);
    title('CONCATENATED Estimation Ref Position') 
    hAx2 = subplot(2,1,2);
    title('CONCATENATED Estimation Position')
    axes(hAx1)
    hold on
    axes(hAx2)
    hold on
    plot(hAx1,estRefPos80);
    plot(hAx2,estPos80);
end

clear expIdx;


C_sys_est_80{1} = tfest(iddata(estTorque80, estPosError80, Ts), C_num_poles, C_num_zeros); 
W_sys_est_80{1} = tfest(iddata(estPos80, estRefPos80, Ts), W_num_poles, W_num_zeros); 

[C_bestModel80, C_bestModelFit80, C_bestModelOutput80] = bestModelFinder(C_sys_est_80, C_testing_iddata);
[W_bestModel80, W_bestModelFit80, W_bestModelOutput80] = bestModelFinder(W_sys_est_80, W_testing_iddata);

% RESULTS
clc
fprintf('\nC: Best model fit from single experiment estimation: %.3f\n', C_bestModelFitSingle);
fprintf('C: Model fit from using 80 experiments for estimation: %.3f\n\n', C_bestModelFit80);

fprintf('\nW: Best model fit from single experiment estimation: %.3f\n', W_bestModelFitSingle);
fprintf('W: Model fit from using 80 experiments for estimation: %.3f\n\n', W_bestModelFit80);

% if(C_bestModelFitSingle>C_bestModelFit80)
%     C_bestModelSingle
% else
%     C_bestModel80
% end

close all

%controller
if C_plot_signals
    for i=1:1:length(C_testing_iddata)
        figure
        titleStr=['Test ',num2str(i)]; 
        if(C_bestModelFitSingle>C_bestModelFit80)
            titleStr=['Single is better - ', titleStr];
        else
            titleStr=['80% is better - ', titleStr];
        end
        if flip_1_every_2
            titleStr=['flip1every2 - ', titleStr]; 
        end
        titleStr=['C - ', titleStr]; 

        sgtitle(titleStr); 
        subplot(2,1,1)
        hold on;
        title('Position of best W estimated with a single experiment') 
        plot(0:length(testingTorque{i})-1,testingTorque{i});
        plot(0:length(C_bestModelOutputSingle{i})-1,C_bestModelOutputSingle{i});
        legend('testing', 'estimated model output')
        subplot(2,1,2)
        hold on;
        title('Position of W estimated with 80% of experiments') 
        plot(0:length(testingTorque{i})-1,testingTorque{i});
        plot(0:length(C_bestModelOutput80{i})-1,C_bestModelOutput80{i});
        legend('testing', 'estimated model output')
    end
end

clear i;

%whole model
if W_plot_signals
    for i=1:1:length(W_testing_iddata)
        figure
        titleStr=['Test ',num2str(i)]; 
        if(W_bestModelFitSingle>W_bestModelFit80)
            titleStr=['Single is better - ', titleStr];
        else
            titleStr=['80% is better - ', titleStr];
        end
        if flip_1_every_2
            titleStr=['flip1every2 - ', titleStr]; 
        end
        titleStr=['W - ', titleStr]; 

        sgtitle(titleStr); 
        subplot(2,1,1)
        hold on;
        title('Position best W estimated with a single experiment') 
        plot(0:length(testingPos{i})-1,testingPos{i});
        plot(0:length(W_bestModelOutputSingle{i})-1,W_bestModelOutputSingle{i});
        legend('testing', 'estimated model output')
        subplot(2,1,2)
        hold on;
        title('Position of W estimated with 80% of experiments') 
        plot(0:length(testingPos{i})-1,testingPos{i});
        plot(0:length(W_bestModelOutput80{i})-1,W_bestModelOutput80{i});
        legend('testing', 'estimated model output')
    end
end
clear i;

%% Extracting controllers from estimated W and compare them with the estimated C
clc;
s = tf('s');
j = 0.068;
d = 0.01;
G = 1/(j*s^2 + d*s);

C_single = invC(G, W_bestModelSingle);
newC_Single = zpk(minreal(C_single,0.5))
C_BM_Single = zpk(C_bestModelSingle)

C_80 = invC(G, W_bestModel80);
newC_80 = zpk(minreal(C_80,0.5))
C_BM_80 = zpk(C_bestModel80)
