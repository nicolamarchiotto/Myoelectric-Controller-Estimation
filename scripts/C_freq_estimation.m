clear
%%
clc
pathToGitFolder = 'C:\\Users\\nicol\\Desktop\\rpcProject\\myo_tools_testing\\RPC_MN_project\\';

load(pathToGitFolder + "resultStructures\\C_resultTable.mat", 'C_resultTable')
load(pathToGitFolder + "resultStructures\\W_resultTable.mat", 'W_resultTable')
load(pathToGitFolder + "resultStructures\\W2_resultTable.mat", 'W2_resultTable')
load(pathToGitFolder + "resultStructures\\C_freq_resultTable.mat", 'C_freq_resultTable')
%% Assumed mechanical model
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
clc
C_num_zeros = 1;
C_num_poles = 1;
Ts = 0.01;
t_arch = ArchitectureEnum.COMP_NONE;


C_best = C_resultTable{t_arch,5}{1};
C_from_Wbest = W_resultTable{t_arch,7}{1};
Wbest_no_minreal = W_resultTable{t_arch,5}{1};
C_from_Wbest_no_minreal = zpk(getC_from_G_and_W_no_minreal(G_assumed, Wbest_no_minreal, t_arch, beta_pos, beta_force, beta_force_int, beta_adm, beta_imp, Kp_pos, Kd_pos, posPole, J_adm, D_adm, K_imp));

N = 10;
t = 0:Ts:N-1;

f0=1;
f1=10;

sweep = exp(-t).*sin(pi*(f0*t + ((f1 - f0)*t.^2)/2*Ts));
[y,l] = lsim(C_from_Wbest_no_minreal, sweep, t);

% figure
% plot(t, sweep, t, y)
% legend('sweep', 'y')

F_iddata = fft(iddata(y, sweep', Ts));

C_freq_est = tfest(F_iddata, 1, 1);
zpk(C_freq_est)
[C_freq_num, C_freq_den] = tfdata(zpk(C_freq_est), 'v');


C_freq_resultTable(t_arch, :) = {char(t_arch), C_num_zeros, C_num_poles, {zpk(C_freq_est)}, {C_freq_num}, {C_freq_den}, {zero(C_freq_est)}, {pole(C_freq_est)}};

%% Simulink responses
close all
sim_imageSavePath = pathToGitFolder + "images\\simulink_responses\\";
sim_saveImage = true;   
sim_betaPath = "allBetas\\";


simulink_plot_function(sim_saveImage, sim_imageSavePath, t_arch, sim_betaPath, out)
