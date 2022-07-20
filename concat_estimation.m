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

plot_signals=true;

% variable which if set to true enables method 3
flip_1_every_2=false;

%% DATA CLEANING

allTrimmedTorque = {};
allTrimmedPosError = {};

i = 1;

if plot_signals
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
    if plot_signals
        plot(hAx1,ones(length(position),1)*referencePosVal-position);
        plot(hAx2,torque);
        plot(hAx3,trimmedPosErr);
        plot(hAx4,trimmedTorque);
    end
    
    allTrimmedTorque{i} = trimmedTorque;
    allTrimmedPosError{i} = trimmedPosErr;
    i = i + 1;
end

clear expId i;

%% Construct estimation and testing set

estimationPosError={};
estimationTorque={};

testingPosError={};
testingTorque={};

numOfTest=ceil(length(allTrimmedTorque)*0.2);

% vector containig the testing indeces, randperm guarantees distinct
% indeces
testing_indeces=randperm(length(allTrimmedTorque),numOfTest);

%support variables used to not modify "all" data structures
deprecatedTrimmedPosError=allTrimmedPosError;
deprecatedTrimmedTorque=allTrimmedTorque;

i=1;
for idx = testing_indeces  %for each experiment    
    testingPosError{i} = allTrimmedPosError{idx};
    testingTorque{i} = allTrimmedTorque{idx};
    
    i=i+1;
    
    deprecatedTrimmedPosError{idx}=[];
    deprecatedTrimmedTorque{idx}=[];
end

clear i idx;

i=1;

% For method 3 set flip_1_every_2 to true in IMPORTANT PARAMETERS paragraph
for idx = 1:1:length(deprecatedTrimmedPosError)
   if(isempty(deprecatedTrimmedPosError{idx}))
       continue;
   end
   if flip_1_every_2
       if mod(i,2)==0
           estimationPosError{i}=flip(deprecatedTrimmedPosError{idx});
           estimationTorque{i}=flip(deprecatedTrimmedTorque{idx});
       else 
           estimationPosError{i}=deprecatedTrimmedPosError{idx};
           estimationTorque{i}=deprecatedTrimmedTorque{idx};   
       end
   else
       estimationPosError{i}=deprecatedTrimmedPosError{idx};
       estimationTorque{i}=deprecatedTrimmedTorque{idx};    
   end
   i=i+1;
end

clear i idx;
clear deprecatedTrimmedPosError deprecatedTrimmedTorque;

% Number of poles and zeros for the controller
num_zeros = 2;
num_poles = 2;

% construct iddata structure for testing
testing_iddata={};
allTestingTorque=[];
for testIdx=1:1:length(testingTorque)
    testing_iddata{testIdx} = iddata(testingTorque{testIdx}, testingPosError{testIdx}, Ts); 
    allTestingTorque=[allTestingTorque;testingTorque{testIdx}]; 
end

clear testIdx;

% Method 1) Controller estimation using the single experiments
C_sys_est={};
for expIdx=1:1:length(estimationPosError)
    est_iddata = iddata(estimationTorque{expIdx}, estimationPosError{expIdx}, Ts); 
    C_sys_est{expIdx} = tfest(est_iddata, num_poles, num_zeros); 
end

clear expIdx;

[bestModelSingle, bestModelFitSingle, bestModelOutputSingle] = bestModelFinder(C_sys_est, testing_iddata);

% Method 2) Controller estimation using the 80% of experiments
estPosError80=[];
estTorque80=[];
C_sys_est_80={};
for expIdx=1:1:length(estimationPosError)
    estPosError80 = [estPosError80; estimationPosError{expIdx}];
    estTorque80 = [estTorque80; estimationTorque{expIdx}];
end

if plot_signals
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

clear expIdx;


C_sys_est_80{1} = tfest(iddata(estTorque80, estPosError80, Ts), num_poles, num_zeros); 

[bestModel80, bestModelFit80, bestModelOutput80] = bestModelFinder(C_sys_est_80, testing_iddata);

% RESULTS
clc
fprintf('\nBest model fit from single experiment estimation: %.3f\n', bestModelFitSingle);
fprintf('Model fit from using 80 experiments for estimation: %.3f\n\n', bestModelFit80);

% if(bestModelFitSingle>bestModelFit80)
%     bestModelSingle
% else
%     bestModel80
% end

close all

for i=1:1:length(testing_iddata)
    figure
    titleStr=['Test ',num2str(i)]; 
    if(bestModelFitSingle>bestModelFit80)
        titleStr=['Single is better - ', titleStr];
    else
        titleStr=['80% is better - ', titleStr];
    end
    if flip_1_every_2
        titleStr=['flip1every2 - ', titleStr]; 
    end
   
    sgtitle(titleStr); 
    subplot(2,1,1)
    hold on;
    title('Torque of best C estimated with a single experiment') 
    plot(0:length(testingTorque{i})-1,testingTorque{i});
    plot(0:length(bestModelOutputSingle{i})-1,bestModelOutputSingle{i});
    legend('testing', 'estimated model output')
    subplot(2,1,2)
    hold on;
    title('Torque of C estimated with 80% of experiments') 
    plot(0:length(testingTorque{i})-1,testingTorque{i});
    plot(0:length(bestModelOutput80{i})-1,bestModelOutput80{i});
    legend('testing', 'estimated model output')
end
clear i;
