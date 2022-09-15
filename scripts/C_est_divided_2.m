
%% DATA FILTERING

Ts = 0.01;

torqueEndEpsilon = 0.2;

maxSignalLength = 400;

C_plot_signals = false;
W_plot_signals = false;

% Outliers data filtering based on standard deviation
% Remove outliers of a vector where an outlier is defined as a point more 
% than three standard deviations away from the mean of the data.
stdOutlierRemoval = true;

% Trim signals at start
trimAtStart = true;

clc
% close all

allTrimmedPosErr = [];
allTrimmedTorque = [];

allTrimmedPosRef = [];
allTrimmedPos = [];

if C_plot_signals
    figure; 
    hAx1 = subplot(2,1,1);   
    title('Position error of the experiments') 
    xlabel('time [cs]');
    ylabel('position [rad]');
    
    hAx2 = subplot(2,1,2);
    title('Torque output of experiments')
    xlabel('time [cs]');
    ylabel('position [rad]');
    
    figure; 
    hAx3 = subplot(2,1,1);    
    title('FILT: Position error of the experiments')
    xlabel('time [cs]');
    ylabel('position [rad]');
    
    hAx4 = subplot(2,1,2);
    title('FILT: Torque output of experiments')
    xlabel('time [cs]');
    ylabel('position [rad]');
   
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
    title('Reference position of the experiments') 
    xlabel('time [cs]');
    ylabel('position [rad]');
    
    hAx6 = subplot(2,1,2);
    title('Position output of experiments')
    xlabel('time [cs]');
    ylabel('position [rad]');
    
    figure; 
    hAx7 = subplot(2,1,1);    
    title('FILT: Reference position of the experiments')
    xlabel('time [cs]');
    ylabel('position [rad]');
    
    hAx8 = subplot(2,1,2);
    title('FILT: Position output of experiments')
    xlabel('time [cs]');
    ylabel('position [rad]');
    
    axes(hAx5)
    hold on
    axes(hAx6)
    hold on
    axes(hAx7)
    hold on
    axes(hAx8)
    hold on
end


i = 1;
for expId = selected_indeces  %for each experiment
    tagIdx = actualFixedDataTableExpanded(expId,:).tag_idx;
    
    % vector of position
    position = exoRefSplineCells{tagIdx}.thetaE(actualFixedDataTableExpanded(expId,:).start_idx:actualFixedDataTableExpanded(expId,:).end_idx);
    
    % vector of torque
    torque = exoRefSplineCells{tagIdx}.reference(actualFixedDataTableExpanded(expId,:).start_idx:actualFixedDataTableExpanded(expId,:).end_idx);
    
    %reference position value
    referencePosVal = actualFixedDataTableExpanded(expId,:).target_angle;
    
    % Torque outliers elimination
    if abs(torque(end)) > torqueEndEpsilon
        continue;
    end    

    firstValOfTorque = torque(1);
    firstValOfPos = position(1);

    lastValOfTorque = torque(end);
    lastValOfPos = position(end);
    
    adjustedTorque = []; 
    adjustedPos = [];
    
    if(length(torque) < maxSignalLength)
        pointsToAddAtEnd = maxSignalLength - length(torque);
        adjustedPos = [position; ones(pointsToAddAtEnd,1)*lastValOfPos]; 
        adjustedTorque = [torque; zeros(pointsToAddAtEnd,1)];
    else
        adjustedTorque = torque(1:maxSignalLength);
        adjustedPos = position(1:maxSignalLength);
    end
    
    adjustedTorque = [adjustedTorque; zeros(maxSignalLength,1)];
    adjustedPos = [adjustedPos; ones(maxSignalLength,1)*adjustedPos(end)];
    adjustedPosRef = ones(length(adjustedPos),1)*referencePosVal;
    adjustedPosErr = adjustedPosRef-adjustedPos;

    allTrimmedPosErr(i,:) = adjustedPosErr';
    allTrimmedTorque(i,:) = adjustedTorque';

    allTrimmedPosRef(i,:) = adjustedPosRef';
    allTrimmedPos(i,:) = adjustedPos';

    i = i+1;
    
    if C_plot_signals
        plot(hAx1, (ones(length(position),1)*referencePosVal)-position);
        plot(hAx2, torque);
    end
    
    if W_plot_signals
        plot(hAx5, ones(length(position),1)*referencePosVal);
        plot(hAx6, position);
    end
end

% STANDARD DEVIATION OUTLIER REMOVAL
% Remove outliers of a vector where an outlier is defined as a point more 
% than three standard deviations from the mean of the data.

if stdOutlierRemoval
    [B,TF]=rmoutliers(allTrimmedTorque,'mean');
    idcs = flip(find(TF)');

    for expId = idcs 
        allTrimmedPosErr(expId,:) = [];
        allTrimmedTorque(expId,:) = [];
        
        allTrimmedPosRef(expId,:) = [];
        allTrimmedPos(expId,:) = [];    
    end
end


% Trim at start, only of align at max is disabled 
if trimAtStart
    supPosErr = [];
    supTorque = [];
    
    
    supPosRef = [];
    supPos = [];

    for expId=1:1:size(allTrimmedTorque,1)
        supPosErr(expId,:) = allTrimmedPosErr(expId, trimStartIndex:end);
        supTorque(expId,:) =  allTrimmedTorque(expId, trimStartIndex:end);
        
        supPosRef(expId,:) = allTrimmedPosRef(expId, trimStartIndex:end);
        supPos(expId,:) =  allTrimmedPos(expId, trimStartIndex:end);
    end
    allTrimmedPosErr = supPosErr;
    allTrimmedTorque = supTorque;
    
    allTrimmedPosRef = supPosRef;
    allTrimmedPos = supPos;

    clear supTorque supPos supPosErr supPosRef;
end
 
if C_plot_signals || W_plot_signals
    for testIdx=1:1:size(allTrimmedTorque,1)
        if C_plot_signals
            plot(hAx3, allTrimmedPosErr(testIdx, :));
            plot(hAx4, allTrimmedTorque(testIdx, :));
        end
        if W_plot_signals
            plot(hAx7, allTrimmedPosRef(testIdx, :));
            plot(hAx8, allTrimmedPos(testIdx, :));
        end
    end
end
clear hAx1 hAx2 hAx3 hAx4 hAx5 hAx6 hAx7 hAx8 firstValOfPos firstValOfTorque expId 
clear adjustedPos adjustedTorque adjustedPosRef adjustedPosErr i pointsToAddStart pointsToAddAtEnd
clear torque position referencePosVal lastValOfPos lastValOfTorque

% size(allTrimmedPos,1)
% architecture
%  close all
%% CREATING IIDATA AND MODEL ESTIMATION

clc;

% Number of poles and zeros for the estimated models

% Assuming a PD for the controller estimation
C_num_zeros = 1;
C_num_poles = 1;

% The mechanical model G is assumed 1/(J*s^2+D*s)


switch architecture
        %GRAV COMP
    case ArchitectureEnum.COMP_NONE
        W_num_zeros = 1;
        W_num_poles = 3;
        
        % FORCE
    case ArchitectureEnum.FORCE_PLAIN_P 
        W_num_zeros = 1;
        W_num_poles = 3;
    case ArchitectureEnum.FORCE_MULTICH8
        W_num_zeros = 1;
        W_num_poles = 3;
    
        %POS_V
    case ArchitectureEnum.POS_V_PLAIN_P 
        W_num_zeros = 2;
        W_num_poles = 5;
    case ArchitectureEnum.POS_V_MULTICH8
        W_num_zeros = 2;
        W_num_poles = 5;
    
        % IMP
    case ArchitectureEnum.FIX_IMP_PLAIN_P
        W_num_zeros = 1;
        W_num_poles = 4;
    case ArchitectureEnum.FIX_IMP_MULTICH8
        W_num_zeros = 1;
        W_num_poles = 4;
        
        
        % ADM
    case ArchitectureEnum.ADM_PLAIN_P 
        W_num_zeros = 2;
        W_num_poles = 6;
    case ArchitectureEnum.ADM_MULTICH8
        W_num_zeros = 2;
        W_num_poles = 6;
    
        % FORCE INT
    case ArchitectureEnum.FORCE_INT_PLAIN_P 
        W_num_zeros = 1;
        W_num_poles = 4;
    case ArchitectureEnum.FORCE_INT_MULTICH8
        W_num_zeros = 1;
        W_num_poles = 4;
end

% For each esperiment estimate the models and construct iddata structures, O(n)
C_iddata = {};
C_sys_est = {};

W_iddata = {};
W_sys_est = {};

% option to force estimated model to be stable
opt = tfestOptions('EnforceStability',true,'InitialCondition','estimate');

for expIdx=1:1:size(allTrimmedTorque,1)

    C_est_iddata = iddata(allTrimmedTorque(expIdx,:)', allTrimmedPosErr(expIdx,:)', Ts);
    C_iddata{expIdx} = C_est_iddata;
    C_sys_est{expIdx} = tfest(C_est_iddata, C_num_poles, C_num_zeros, opt);

    W_est_iddata = iddata(allTrimmedPos(expIdx,:)', allTrimmedPosRef(expIdx,:)', Ts); 
    W_iddata{expIdx} = W_est_iddata;
    W_sys_est{expIdx} = tfest(W_est_iddata, W_num_poles, W_num_zeros, opt);    
end

% initialization of variables for plot function
C_bestModelOutput = {};
W_bestModelOutput = {};

clear expIdx C_est_iddata W_est_iddata opt;

%% FIND BEST MODEL, C AND W
%
%
% CAREFULL, each bestModelFinder is O(n^2)
%
%
[C_bestModel, C_bestModelFit, C_bestModelOutput] = bestModelFinder(C_sys_est, C_iddata);

[W_bestModel, W_bestModelFit, W_bestModelOutput] = bestModelFinder(W_sys_est, W_iddata);

% RESULTS
clc
fprintf('\nC best model fit: %.3f\n', C_bestModelFit);

fprintf('\nW best model fit: %.3f\n', W_bestModelFit);
%
%
%% GET C AND W FROM RESULT STRUCTURES
%
%
% clc
% C_bestModel = C_resultTable{architecture,5}{1};
% W_bestModel = W_resultTable{architecture,5}{1};
% 

%% SIMULINK PARAMS AND MECHANICAL MODEL
s = tf('s');
J = 0.068;
D = 0.085; %da 0.1 a 0.001

% saturation position parameter for simulink
posLim = 3;

% beta parameters for architectures
beta_pos = 4;
beta_force = 3;
beta_force_int = 2;
beta_adm = 0.2;
beta_imp = 4;

% position ctrl parameters
Kp_pos = 140;
Kd_pos = 2;
Ki_pos = 0;
% pole to make the pos ctrl fisible
posPole = 1000;

% admittance ctrl parameters
J_adm = 0.05;
D_adm = 0.25;

% impedance ctrl parameters
K_imp = 8;

G_assumed = 1/(J*s^2 + D*s);
%%
%
%
% Estimated controller from posErr and torque
%
%
clc
C_best =  zpk(C_bestModel);
[C_best_num, C_best_den] = tfdata(C_best, 'v');

%
%
% Extracted controller from estimated W using estimated G
%
%

zpk(W_bestModel);
C_from_Wbest = zpk(getC_from_G_and_W(G_assumed, W_bestModel, architecture, beta_pos, beta_force, beta_force_int, beta_adm, beta_imp, Kp_pos, Kd_pos, posPole, J_adm, D_adm, K_imp));
[C_from_Wbest_num, C_from_Wbest_den] = tfdata(C_from_Wbest, 'v');

%
%
%% FOURIER TRANSFORM ESTIMATION 
%
%
C_freq_num_zeros = 1;
C_freq_num_poles = 1;
C_from_Wbest_no_minreal = tf(getC_from_G_and_W_no_minreal(G_assumed, W_bestModel, architecture, beta_pos, beta_force, beta_force_int, beta_adm, beta_imp, Kp_pos, Kd_pos, posPole, J_adm, D_adm, K_imp));

N = 10;
t = 0:Ts:N-1;

f0=1;
f1=10;

sweep = exp(-t).*sin(pi*(f0*t + ((f1 - f0)*t.^2)/2*Ts));
[y,l] = lsim(C_from_Wbest_no_minreal, sweep, t);

% figure
% plot(t, sweep, t, y)
% legend('sweep', 'y')

% Use fft and ifft to transform existing iddata objects to and from the time and frequency domains.
% https://it.mathworks.com/help/ident/ref/iddata.html
F_iddata = fft(iddata(y, sweep', Ts));
% figure 
% plot(F_iddata)
C_freq_est = tfest(F_iddata, C_freq_num_zeros, C_freq_num_poles);
[C_freq_num, C_freq_den] = tfdata(zpk(C_freq_est), 'v');

%% RESULT SAVING

C_resultTable(architecture, :) = {char(architecture), C_num_zeros, C_num_poles, size(allTrimmedTorque,1), {C_best}, C_bestModelFit, {C_best_num}, {C_best_den}, {pole(C_best)}};

W_resultTable(architecture, :) = {char(architecture), W_num_zeros, W_num_poles, size(allTrimmedTorque,1), {zpk(W_bestModel)}, W_bestModelFit, {C_from_Wbest}, {C_from_Wbest_num}, {C_from_Wbest_den}, {pole(C_from_Wbest)}};

C_freq_resultTable(architecture, :) = {char(architecture), C_freq_num_zeros, C_freq_num_poles, {zpk(C_freq_est)}, {C_freq_num}, {C_freq_den}, {zero(C_freq_est)}, {pole(C_freq_est)}};
%%
% Validation plots
close all
imageSavePath = pathToGitFolder + "images\\";
saveImage = true;
plotC = true;
plotW = true;
betaPath = "allBetas\\";
CW_plotFunction(plotC, plotW, allTrimmedTorque, allTrimmedPos, C_bestModelOutput, W_bestModelOutput, saveImage, imageSavePath, architecture, betaPath )
close all

%% Simulink responses plot
clc
close all
sim_imageSavePath = pathToGitFolder + "images\\simulink_responses\\";
sim_saveImage = true;   
sim_betaPath = "allBetas\\";
stepAmp = allTrimmedPosRef(1,1)-allTrimmedPos(1,1);

out = sim('myo_model_architectures');

simulink_plot_function(sim_saveImage, sim_imageSavePath, architecture, sim_betaPath, out, allTrimmedPosRef, allTrimmedPos)
