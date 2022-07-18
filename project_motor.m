%% ec45
f = load('ec45_sweep.csv');
t = f(:,1);
current = f(:,2);
ref = f(:,3);
theta = f(:,4);
dtheta = f(:,5);
ddtheta = f(:,6);

kt = 1.2685;
tau = current*kt;

ts = mean(diff(t));

half_sampling_freq = 1/(2.*ts);
cut_off_freq = 20; % [Hz]

[B,A] = butter(8,cut_off_freq/half_sampling_freq);

currentFiltered = filter(B,A, current);
tauFiltered = filter(B,A, tau);
thetaFiltered = filter(B,A, theta);
dthetaFiltered = filter(B,A, dtheta);
ddthetaFiltered = filter(B,A, ddtheta);

%expected controller
kd = 2;
kp = 50;
err = ref-thetaFiltered;

%%
kd = 2;
kp = 50;
jm = 0.02496;
dm = 0.03;
s = tf('s');
G = 1/(s*(jm*s+dm));
Tau = 1/(2*pi*100);
C = kp+kd*s/(Tau*s+1);
forw = C*G;%*kt;
H = forw/(1+forw);
%H = feedback(forw,1);
ts = 1e-4;


% 80.12 s + 2003
%  --------------------
%  s^2 + 81.32 s + 2003

Hd = c2d(H,ts,'tustin');
Gd = c2d(G,ts,'tustin');
Cd = c2d(C,ts,'tustin');

%% simulink model.slx
simRef = out.simRef.Data;
simTheta = out.simTheta.Data;
simErr = out.simErr.Data;
simCurr = out.simCurr.Data;

%%
systemIdentification('ec45');
%export tf_ec45, tf_ec45_filtered

%%
checkC = zpk(minreal(H/(G-H*G),-0.9));
model = tf_sim;
newC = zpk(minreal(model/(G-model*G),0.1));
%  1914.8 (s+24.33)
%  ----------------
%     (s+927.8)

%https://www.researchgate.net/post/How_I_could_know_the_dominant_poles_of_a_dynamic_control_system_if_I_have_the_transfer_function

%%
z = tf('z',ts); 
estCz = (1 + 0.9988*z^-1)/(1 - z^-1); %from simBj
estC = d2c(estCz, 'tustin');
estC = zpk(minreal(estC));

%%
z = tf('z',ts); 
estCz = (1 + 1.479* z^-1 + 0.8223 *z^-2)/(1 - 1.998 *z^-1 + z^-2);  %from bj22220
estC = d2c(estCz, 'tustin');
estC = zpk(minreal(estC));
%estC = (s^2 + 1625*s + 2.609e06)/(s^2 + 0.2224*s + 1185);











%%
%tf_ec45
%  34.7
%  ---------
%  s + 6.255
est_jm = 1/34.7 %0.0288, from script: 0.0290, from datasheet 0.024963
est_dm = 6.255/34.7 %0.1803, from script: 0.0818, from datasheet 0.03

%tf_ec45_filtered %almost identical to non filtered version


%%
%export bj_cut from system identification
%contEst = d2c(bj_cut, 'tustin');
%estC = (s^3 + 2.422e11* s^2 + 9.153e12 *s + 3.235e14)/(s^3 + 39.02 *s^2 + 1382 *s + 1650);

%% maxon dcx22L
f = load('maxon_sweep.csv');
t = f(:,1);
current = f(:,6);
ref = f(:,2);
theta = f(:,3);
dtheta = f(:,4);
ddtheta = f(:,5);

kt = 1.5038;
tau = current*kt;

ts = mean(diff(t));

half_sampling_freq = 1/(2.*ts);
cut_off_freq = 20; % [Hz]

[B,A] = butter(8,cut_off_freq/half_sampling_freq);

tauFiltered = filter(B,A, tau);
dthetaFiltered = filter(B,A, dtheta);

%%
kd = 0.3;
kp = 15;
jm = 0.0104;
dm = 0.0068;
s = tf('s');
G = 1/(s*(jm*s+dm));
C = kp+kd*s;%/(ts*s+1);
forw = C*G;%*kt;
H = forw/(1+forw);
% 80.12 s + 2003
%  --------------------
%  s^2 + 81.32 s + 2003

Hd = c2d(H,ts,'tustin');
Gd = c2d(G,ts,'tustin');
Cd = c2d(C,ts,'tustin');

%%
systemIdentification('maxon');
%export tf_maxon, tf_maxon_filtered

%%
tf_maxon
%  155.5
%  ---------
%  s + 5.297
jm = 1/155.5 %0.0064, from script: 0.0064, from datasheet 0.0104
dm = 5.297/155.5 %0.0341, from script: 0.0068 , from datasheet 0.0068

tf_maxon_filtered %almost identical to non filtered version