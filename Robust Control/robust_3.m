%% ROBUST ASSIGNMENT 3

clear; clc; close all;

load assignment_data_sc42145
load mixsyn_controller.mat;

run plot_settings

sys = ss(A,B,C,D);
G_MIMO = tf(sys);
G = G_MIMO(:,1:2);

s = tf('s');

%% 2. Build generalized plant
A = 1e-4;
omega_B = 0.4*2*pi;
M = 1.8;

% Performance weight
Wp11 = (s/M + omega_B)/(s + omega_B*A);
Wp = blkdiag(Wp11, 0.2);

% Input weight
Wu = blkdiag(0.01, tf([5e-3 7e-4 5e-5], [1 14e-4 1e-6]));

% Input uncertainty weight
Wi11 = (1/16/pi*s + 0.3)/(1/64/pi*s + 1);
Wi = blkdiag(Wi11, Wi11);

% Output uncertainty weight
Wo11 = (0.05*s + 0.2)/(0.01*s + 1);
Wo = blkdiag(Wo11, Wo11);

P = [zeros(2, 6) Wi; ... % delta y i
     Wo*G zeros(2, 4) Wo*G; ... % delta y o
     -Wp*G -Wp Wp -Wp*G; ... % z1
     zeros(2,6) Wu; ... % z2
     -G -eye(2) eye(2) -G]; % v

N = minreal(lft(P, K2), [], false);
Delta = ultidyn('Delta', [4 4], 'Bound', 1);

Delta_i = ultidyn('Delta_i', [2 2], 'Bound', 1);
Delta_o = ultidyn('Delta_o', [2 2], 'Bound', 1);

%% 4. NS, RS, RP analysis
% Analyse nominal stability
L = K2*G;
nominal_stability = max(real(eig(N))) < 0;
if nominal_stability
    disp('Nominal stability OK');
else
    disp('Nominal stability NOK');
end

% Nominal performance
N22 = N(5:8, 5:6);
[bnds_N22, ~] = mussv(N22, [2 4], 's');
nominal_performance = nominal_stability && all(norm(bnds_N22(:, 1), inf) < 1);
if nominal_performance
    disp('Nominal performance OK');
else
    disp('Nominal performance NOK');
end

% Robust stability
M = N(1:4, 1:4);
[bnds_M, ~] = mussv(M, [4 0], 's');

robust_stability = nominal_stability && all(norm(bnds_M(:, 1), inf) < 1);
if robust_stability
    disp('Robust stability OK');
else
    disp('Robust stability NOK');
end

% Robust performance
[bnds_N, ~] = mussv(N, [4 0; 2 4], 's');
robust_performance = nominal_stability && all(norm(bnds_N(:, 1), inf) < 1);
if robust_performance
    disp('Robust performance OK');
else
    disp('Robust performance NOK');
end

%% Visualisation
figure;
specialbode(Wi11, 'PhaseVisible', 'off');
exportgraphics(gcf, 'C:\Users\emiel\Dropbox (Personal)\Apps\Overleaf\Robust and Multivariable Control Design (Dropbox)\Media_3\wi11.eps');

figure;
specialbode(Wo11, 'PhaseVisible', 'off');
exportgraphics(gcf, 'C:\Users\emiel\Dropbox (Personal)\Apps\Overleaf\Robust and Multivariable Control Design (Dropbox)\Media_3\wo11.eps');

figure;
bodemag(Wi*Delta_i);
set(gcf, 'Position', get(gcf, 'Position').*[1 1 1.3 1.5])
exportgraphics(gcf, 'C:\Users\emiel\Dropbox (Personal)\Apps\Overleaf\Robust and Multivariable Control Design (Dropbox)\Media_3\wi_delta.eps');

figure
bodemag(Wi*Delta_o);
set(gcf, 'Position', get(gcf, 'Position').*[1 1 1.3 1.5])
exportgraphics(gcf, 'C:\Users\emiel\Dropbox (Personal)\Apps\Overleaf\Robust and Multivariable Control Design (Dropbox)\Media_3\wi_delta.eps');

figure
tmp = lft(Delta, P);
sigma((eye(2) + Wo*Delta_o*eye(size(G)))*G*(eye(2) + Wi*Delta_i*eye(size(G))));
set(gcf, 'Position', get(gcf, 'Position').*[1 1 1 1.3]); grid;
exportgraphics(gcf, 'C:\Users\emiel\Dropbox (Personal)\Apps\Overleaf\Robust and Multivariable Control Design (Dropbox)\Media_3\sigma.eps');

% Nyquist plot
% [re, im] = gen_nyquist(L);
% figure
% tile = tiledlayout(1, 2, 'Padding', 'compact', 'Tilespacing', 'compact');
% nexttile
% plot(re, im); hold on;
% plot(re, -im); hold on;
% title('Global');
% set(gca, 'XAxisLocation', 'origin');
% set(gca, 'YAxisLocation', 'origin');
% xlabel('Re'); ylabel('Im');
% nexttile
% plot(re, im); hold on;
% plot(re, -im); hold on;
% title('Focus on origin');
% xlabel('Re'); ylabel('Im');
% xlim([-20 10]);  ylim([-15 15])
% set(gca, 'XAxisLocation', 'origin');
% set(gca, 'YAxisLocation', 'origin');
% title(tile, '\textbf{Generalized Nyquist}', 'interpreter', 'latex')

