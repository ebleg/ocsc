clear; clc; close all;

load assignment_data_sc42145

s = tf('s');

%% Transfer functions
sys = ss(A,B,C,D);
G_MIMO = tf(sys);
G = G_MIMO(1, 1);

% figure('Name', 'Plant frequency response')
% margin(G)

% pzmap(G)

% Kp = -5;
% Ti = 4;
% Td = -10;
% Tf = 50;
% K = pid(Kp, Kp/Ti, Kp*Td, Tf);
% figure('Name', 'Controller FRF')
% bode(K)


%% Low pass idea
% w_cut = 0.5;
% lowpass = tf([0 0 w_cut^2], [1 2*w_cut*0.7 w_cut^2]);
% bode(lowpass); hold on;
% bode(G); 
% bode(K);

% L = minreal(lowpass*K*G);

% bode(L); hold on;
% bode(G);
% figure('Name', 'Loop gain FRF');
% margin(L)

% G_new = lowpass*G;



% figure
% margin(T)

%% Requirements
OS = 1;
zeta_target = -log(OS/100)/sqrt(pi^2 + log(OS/100)^2);
PM_target = rad2deg(atan(2*zeta_target/...
                    sqrt(-2*zeta_target^2 + sqrt(1 + 4*zeta_target^4))));

BW_factor = 4/zeta_target*sqrt((1-2*zeta_target^2) ...
            + sqrt(4*zeta_target^4 - 4*zeta_target^2 + 2));
        
load('C_best_2.mat');
[num, den] = tfdata(C_best_2);
num = num{1}; den = den{1};
w_lp = sqrt(den(3));
zeta_lp = den(2)/2/w_lp;
Kp = num(3)/w_lp^2;
Ti = w_lp^2*Kp/num(4);

% Check
K = -w_lp^2/(s^2 + 2*zeta_lp*w_lp*s + w_lp^2) ... % Low pass filter
    *Kp*(1 + 1/Ti/s); % PI controller


L = minreal(K*G);
T = minreal(L/(1 + L), 1e-7);
S = minreal(1/(1 + L), 1e-7);

stepinfo(T, 'SettlingTimeThreshold', 0.01)

%% Plots
close all;

% Step response
[y_step, t_step] = step(T);
plot_step = plot(t_step, y_step, 'Color', '#27ae60');
set(plot_step, 'LineWidth', 1.7);
ax = gca;
xlim([0 max(t_step)])
xlabel('Time [s]', 'interpreter', 'latex', 'fontsize', 12);
ylabel('Amplitude', 'interpreter', 'latex', 'fontsize', 12);
title('\textbf{Step response}', 'interpreter', 'latex', 'fontsize', 13);
format_axes(ax);

% Controller effort
figure
[y_input] = step(K/(1+K*G), t_step);
plot_step = plot(t_step, y_input, 'Color', '#27ae60');
set(plot_step, 'LineWidth', 1.7);
ax = gca;
xlim([0 max(t_step)])
xlabel('Time [s]', 'interpreter', 'latex', 'fontsize', 12);
ylabel('Amplitude', 'interpreter', 'latex', 'fontsize', 12);
title('\textbf{Controller effort}', 'interpreter', 'latex', 'fontsize', 13);
format_axes(ax);

%% Sensitivity and complementary sensitivity
% sens_plot = bodeplot(S);
% opt = getoptions(sens_plot);
% opt.PhaseVisible = 'off';
% setoptions(sens_plot, opt);
% 
% ax = gca;
% xlabel('$\omega$', 'interpreter', 'latex', 'fontsize', 12);
% ylabel('$|S|$', 'interpreter', 'latex', 'fontsize', 12);
% title('\textbf{Sensitivity}', 'interpreter', 'latex', 'fontsize', 13);
% format_axes(ax);

%% Complementary sensitivity
comp_plot = bodeplot(T);
opt = getoptions(comp_plot);
opt.PhaseVisible = 'off';
setoptions(comp_plot, opt);

ax = gca;
xlabel('$\omega$', 'interpreter', 'latex', 'fontsize', 12);
ylabel('$|S|$', 'interpreter', 'latex', 'fontsize', 12);
title('\textbf{Complementary sensitivity}', 'interpreter', 'latex', 'fontsize', 13);
format_axes(ax);

%% Pole-zero plot
plot(pole(G), 'x');
hold on;
plot(zero(G), 'o');

%%
format_axes(gca);

%% Plant bode plot
h = bodeplot(G);
opt = getoptions(h);
opt.XLabel.Interpreter = 'latex';
opt.YLabel.Interpreter = 'latex';
opt.Title.Interpreter = 'latex';
setoptions(h, opt)
title('\textbf{Plant Bode diagram}', 'Fontsize', 13);

%% Plant Nyquist
nyquistplot(-G)


function format_axes(ax)
    set(ax, 'TickLabelInterpreter', 'latex', 'TickDir', 'both');
    set(ax, 'XGrid', 'on', 'YGrid', 'on');
    set(ax.XAxis, 'LineWidth', 1);
    set(ax.YAxis, 'LineWidth', 1);
    set(ax, 'Box', 'off');
end

