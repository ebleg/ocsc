%% FILTERING & IDENTIFICATION
% Feres Hassan & Emiel Legrand

clear; clc; close all;

student_number = min(4446100, 4362152);

exciteSystem2 = @(u, fs) exciteSystem(2512634, u, fs);

%% Part 0: Experiment
h = 1/30;
t_max = 20;
tin = 0:h:t_max;
% u = zeros(size(tin))';
% u(tin > 2) = 1000;

u = square(tin/2).'*1e3;
% y = exciteSystem2(u, h);
% save('system_output_inv', 'y');

%% Part 1: Data Preprocessing
load 'system_output_inv';

% Spike removal
y_filt = hampel(y, 4);
figure(1); hold on;
stairs(tin, y, 'DisplayName', 'Original');
stairs(tin, y_filt, 'DisplayName', 'Spikes removed');
legend('Location', 'best'); grid; 
xlabel('Time (s)'); ylabel('System output');
title('Spike removal using Hampel');

% Detrending
y_filt = y_filt - nanmean(y_filt);
u_filt = u - nanmean(u);

% y_filt = movmean(y_filt, 10);
figure(2)
yyaxis right; set(gca, 'YColor', 'black')
stairs(tin, u);
ylabel('Input'); axis('padded');
% xlim([2 4]);
title('Input delay');
xlabel('Time(s)'); grid;
yyaxis left; set(gca, 'YColor', 'black')
stairs(tin, y_filt); ylabel('System output');
hold on;

% Determine the sampling frequency
figure(3)
L = numel(y_filt);
P2 = abs(fft(y_filt))/L; P1 = P2(1:L/2+1);
plot(1/h*(0:(L/2))/L, P1)
xlabel('Frequency (Hz)');
ylabel('Amplitude'); grid;
title('FFT of system output');

% Check the persistency of excitation of the input signal
k = 60;
rank(hankel(u(1:k), u(k:end)))

% System order
n = 6;

% figure(4)
% open '../../../singular_values_500s.fig'
% title('Singular values for large number of samples'); grid;
% line([n n], get(gca, 'YLim'), 'Color', 'black');

%% System identification
K = zeros(n, 1);
[A, B, C, D, x0, sv] = subspaceID(u_filt, y_filt, 50, n, 'po-moesp');
% [pem.Abar, pem.Bbar, pem.C, pem.D, pem.K, pem.x0] = ...
%                    pem(A - K*C, B - K*D, C, D, K, x0, u_filt, y_filt, 300);

figure(4)
semilogy(sv);
title('Singular values'); grid;

figure(5)
y_ss = simsystem(A, B, C, D, x0, u_filt); hold on;
stairs(tin, y_filt);

% y_pem = simsystem(pem.Abar, [pem.Bbar pem.K], pem.C, ...
%                             [pem.D zeros(1, 1)], pem.x0, [u_filt y_filt]);

stairs(tin, y_ss);
% stairs(tin, y_pem);

disp(vaf(y_filt, y_ss));

%% Validation


