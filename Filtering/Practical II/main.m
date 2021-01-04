%% FILTERING & IDENTIFICATION
% Feres Hassan & Emiel Legrand

clear; clc; close all;

student_number = min(4446100, 4362152);

exciteSystem2 = @(u, fs) exciteSystem(4439260, u, fs);


%% Part 1: Data Preprocessing
h = 1/11;
t_max = 300;
tin = 0:h:t_max;
u = zeros(size(tin))';
u(tin > 2) = 1000;

% u = square(tin/2)';%.'*1e3;
y = exciteSystem2(u, h);

% Find sampling time
% Y = fft(y_filt)

% Spike removal
y_filt = hampel(y, 10);

% Detrending
y_filt = y_filt - nanmean(y_filt);
u_filt = u - nanmean(u);
u_filt = u;

% y_filt = movmean(y_filt, 10);
figure(1)
plot(tin, y_filt, '.-'); hold on;
% plot(tin(1:5:end), y_filt(1:5:end), 'LineWidth', 2);

figure(2)
L = numel(y_filt);
P2 = abs(fft(y_filt))/L; P1 = P2(1:L/2+1);
plot(1/h*(0:(L/2))/L, P1)

k = 60;
rank(hankel(u(1:k), u(k:end)))

n = 3;

%% System identification
K = zeros(n, 1);
[A, B, C, D, x0, sv] = subspaceID(u_filt, y_filt, 30, n, 'po-moesp');
% [Abar, Bbar, C, D, K, x0] = pem(A - K*C, ...
%                                 B - K*D, ...
%                                 C, D, K, x0, y, u, 500);
figure(3)
semilogy(sv);

figure(1)
y_ss = simsystem(A, B, C, D, x0, u_filt); hold on;
% y_pem = simsystem(Abar, [Bbar K], C, [D zeros(1, 1)], x0, [u y]);
plot(tin, y_ss);
% plot(tin, y_pem);

disp(vaf(y_filt, y_ss));

