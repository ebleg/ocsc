% Robust Control
%% Assignment II: Multi-variable mixed-sensitivity control.
% Load the plant model.

clear; clc; close all;

load assignment_data_sc42145
sys = ss(A,B,C,D);
G_MIMO = tf(sys);
G = G_MIMO(:,1:2); % Check if they want us to do this

s = tf('s');

% Question 1 
% Compute the RGA for $$\omega=0$$ and $$\omega = 0.4\times2\times \pi$. What 
% can you conclude?

RGA = @(G, omega) evalfr(G, omega).*pinv(evalfr(G, omega)).';

disp('RGA (omega = 0):');
disp(RGA(G, 0))

disp('RGA (omega = 0.8*pi):');
disp(RGA(G, 0.8*pi))

%% 
% *Observations*
%% 
% * Diagonal entries are negative.
% * Diagonal entries are a lot smaller than the anti-diagonal entries.
% * Usually the positive pairings closest to one are desirable for control, 
% as such it seems we want to pair both channels with each other, so u1 to y2 
% and u2 to y1.
% * There are no large (> 10) entries.
% Question 2 
% Compute the zeros and poles of the system.

G = minreal(balreal(G), [], false);

disp('Zeros of the plant')
G_zeros = tzero(G);
disp(G_zeros)

disp('Poles of the plant')
G_poles = pole(G);
disp(G_poles)

% figure; hold on;
% plot(G_poles, 'x');
% plot(G_zeros, 'o');
% title('\textbf{Pole-zero map of the plant}', 'interpreter', 'latex');
% xlabel('Re', 'interpreter', 'latex'); ylabel('Im', 'interpreter', 'latex');
% grid;

%% 
% *Observations*
%% 
% * No RHP poles
% * No RHP zeros (?)
% * Two complex pole pairs
% * One real pole
% * Three complex zero pairs
% * One real zero
% Question 3 
% Design an appropriate $$W_{p_{11}}$$
% 
% For this we used the equation 3.81 from Skogestadt

A = 1e-4;
omega_B = 0.4*2*pi;
M = 1.8;

Wp11 = (s/M + omega_B)/(s + omega_B*A);
Wp = [Wp11, 0; 0 0.2];

% Question 4/5
% Build the generalized plant

% Controller sensitivity weight
Wu = [0.01 0; 0 tf([5e-3 7e-4 5e-5], [1 14e-4 1e-6])];

% Generalized plant
systemnames = 'G Wp Wu';
inputvar = '[w(2); [beta; tau_e]]';
input_to_G = '[beta; tau_e]';
input_to_Wu = '[beta; tau_e]';
input_to_Wp = '[w-G]';
outputvar = '[Wp; Wu; w-G]'; % Output generalized plant
sysoutname = 'P'; cleanupsysic = 'yes';
sysic; P = minreal(P, 1e-6, false);
% Question 6/7

[K1, ~, gam, ~] = hinfsyn(P, 2, 2);
K1 = minreal(K1, 1e-6, false);
K2 = minreal(mixsyn(G, Wp, Wu, 0), 1e-6, false);

% Question 8
% Test the closed loop.

% CL1 = feedback(G*K1, eye(2), -1);
% y1 = step(CL1);
% plot(y1(:,1,1)); hold on;

L = G*K2;
S = minreal(inv(eye(2) + L), [], false);
T = minreal(L*S, [], false);

% Step response
% figure
% tile = tiledlayout(1, 2, 'Padding', 'compact', 'Tilespacing', 'compact');
% nexttile
% [y1, t] = step(T(1,1));
% plot(t, y1);
% title('Step response');
% xlim([min(t) max(t)]);
% nexttile
% KS = tf(minreal(K2*S, [], false));
% y2 = step(KS(1,1), t);
% plot(t, y2);
% xlim([min(t) max(t)]);
% title('Controller effort');
% title(tile, '\textbf{Time simulation: step response}', 'interpreter', 'latex')

% Disturbance rejection
% figure
% tile = tiledlayout(1, 2, 'Padding', 'compact', 'Tilespacing', 'compact');
% nexttile
% GdS = tf(minreal(S*G_MIMO(:, 3), [], false));
% [y3, t] = step(GdS(1,1));
% plot(t, y3);
% title('Disturbance rejection');
% xlim([min(t) max(t)]);
% nexttile
% KGdS = tf(minreal(K2*S*G_MIMO(:, 3), [], false));
% y4 = step(KS(1,1), t);
% plot(t, y4);
% xlim([min(t) max(t)]);
% title('Controller effort');
% title(tile, '\textbf{Time simulation: disturbance rejection}', 'interpreter', 'latex')

% figure
% [S_mag, ~, w_out] = bode(S(1,1));
% [Wp_mag, ~] = bode(1/Wp11, w_out);
% semilogx(w_out, squeeze(20*log10(S_mag)), 'displayname', '$S_{11}$'); hold on;
% semilogx(w_out, squeeze(20*log10(Wp_mag)), 'displayname', '$W_{p}$');
% xlim([min(w_out) max(w_out)]);
% title('\textbf{Sensitivity}', 'interpreter', 'latex');
% xlabel('\textbf{Frequency (rad/s)}', 'interpreter', 'latex');
% ylabel('\textbf{Magnitude (dB)}', 'interpreter', 'latex');
% grid
% legend;

%% 
% Nyquist plot
% [re, im] = gen_nyquist(L);
% figure
% tile = tiledlayout(1, 2, 'Padding', 'compact', 'Tilespacing', 'compact');
% nexttile
% plot(re, im); hold on;
% plot(re, -im); hold on;
% title('Global');
% xlabel('Re'); ylabel('Im');
% nexttile
% plot(re, im); hold on;
% plot(re, -im); hold on;
% title('Focus on origin');
% xlabel('Re'); ylabel('Im');
% title(tile, 'Generalized Nyquist')


% Question 9
% Time domain simulations for both reference tracking and disturbance rejection.

Gd = G_MIMO(1,3);
% Part 2.1
% Design new controller and error weights. The idea is that the controller for 
% the torque $$\tau_e$$ reacts on the high-frequency disturbances and the controller 
% for the pitch angle is sensitive to low-frequency () disturbances)
% 
% *Scale the plant*

G2 = G_MIMO(1,1:2);

% SCALE_BETA = 20;     % degrees
% SCALE_TORQUE = 5e3;  % Nm
% SCALE_OUT = 1/90;   % rad/s

% G2_scaled = SCALE_OUT*G2*diag([SCALE_BETA, SCALE_TORQUE]);
% G2_scaled = G2;

%% 
% *Scale the disturbance dynamics*

SCALE_DIST = 2.5;
% Gd_scaled = SCALE_OUT*Gd*SCALE_DIST;
Gd_scaled = Gd;

%% 
% *Build and tune the weights*
% omega_LP = 2*pi/1000; % Tune roll-off
% omega_HP = 0.1*0.12*2*pi;
% omega_B2 = 0.1*omega_LP;
% 
% Wu211 = tf(makeweight(1, [omega_LP, 2/sqrt(2)], 1e10));
% Wu222 = tf(makeweight(1e10, [omega_LP, 2/sqrt(2)], 1));
% 
% % Wp2 = (s/1.8 + omega_B2)/(s + omega_B2*1e4);
% Wp2 = Wp11;%/Gd_scaled;
% 
% Wu2 = [Wu211 0; 0 Wu222];

% ------------------------- FROM GORRRRGGG --------------------------------

load Assignment_Data_SC42145
V_meas = Wind_Data.Data;
freq_meas = fft(V_meas);
sampling_mean = mean(diff(Wind_Data.Time));

n = numel(V_meas);
omega_rng = 1/sampling_mean*(0:(n/2))/n;
freq_meas = abs(freq_meas/n);

figure
plot(omega_rng, freq_meas(1:n/2+1)) 

omg1 = 2*pi/1000;
omg2 = 0.12*2*pi;
Wu211 = 1/100*(1 + s/omg1)/(1 + s/(100*omg1));
Wu222 = (1 + s/(100*omg1))/(1 + s/(omg1));

% figure;
% bodemag(Wu211); hold on;
% bodemag(Wu222);

Wp2 = Wp11;
Wu2 = [Wu211 0; 0 Wu222];

%% 
% *Connect the systems*
% systemnames = 'G2 Gd Wp2 Wu2';
% inputvar = '[V; beta; tau_e]';
% 
% input_to_G2 = '[beta; tau_e]';
% input_to_Gd = '[V]';
% 
% input_to_Wp2 = '[-G2 - Gd]';
% input_to_Wu2 = '[beta; tau_e]';
% 
% outputvar = '[Wp2; Wu2; -G2 - Gd]';
% sysoutname = 'P2'; cleanupsysic = 'yes';
% sysic;

P2 = [Wp2*Gd -Wp2*G2; zeros(2,1) -Wu2; Gd -G2];

P2 = balreal(minreal(P2, [], false));
%% 
% *Compute the* $$H_\infty$$*-controller* 

[K3, ~, gam2, ~] = hinfsyn(P2, 1, 2);
disp(gam2)
K3 = balreal(minreal(K3, 1e-6, false));

L2 = G2*K3;
S2 = inv(eye(size(L2, 1)) + L2);
T2 = eye(size(S2, 1)) - S2;

figure
bodemag(K3*S2*Gd); hold on;
bodemag([1/Wu211; 1/Wu222])

figure
% bodemag(S2*Gd); hold on;
bodemag(S2*Gd); hold on;
bodemag(1/Wp2)
% bodemag(1/Wp11)