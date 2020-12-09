%% -----------------------------------------------------
%  Filtering and Identification assignment 1
%  -----------------------------------------------------
%
%    Feres Hassan (4362152) 
%       & Emile-Bernard Legrand (4446100)
%
%    December 10, 2020
% -----------------------------------------------------
clear; clc; close all

%% Plot settings

%% Assignment 1
load('calibration.mat')

% Alternatively, without loops

% Question 1a: Error per microphone
y_calib = y; % This will be overwritten later by new data
TOA_true = mean(y_calib, 2);
mic_error = y_calib - TOA_true;

% Question 1b: Microphone bias
mic_bias = mean(mic_error);

% Question 1c: Microphone variance
mic_variance = var(mic_error - mic_bias);

% Question 1d: Visualization with histogram
[N, l] = hist(mic_error(:,1),20);
Wb = l(2)-l(1); % Bin width
Ny = length(mic_error(:,1)); % Nr of samples
fig1 = bar(l, N/(Ny*Wb));

% figure
% hist_microphone = 1;
% histogram(mic_error(:,hist_microphone), 100, 'Normalization', 'pdf');
% title(sprintf('Error of microphone %d', hist_microphone))

% figure
% bar(mic_bias)

%% Assignment 2: Nonlinear least-squares
load('experiment.mat')

% Question 2a
y_exp = y;
N_exp = size(y_exp, 1);
y_bias_correct = y_exp - repmat(mic_bias, N_exp, 1);

theta_hat = nan(3,N_exp+1);
diagP = nan(3, N_exp);
theta_hat(:,1) = [0.1 0.6 0]';

% Question 2b
for k = 2:(N_exp+1)
    [theta_hat(:,k), diagP(:,k-1), F] = nls(y_bias_correct(k-1,:)', ...
        sqrt(mic_variance), ...
        theta_hat(:,k-1), ...
        100, mic_locations);
end
theta_hat(:,1) = []; % Remove the initial guess to obtain the correct length

% Question 2c: Visualize the results
figure
fig2 = plotresults(theta_hat(1:2,:), diagP(1:2,:)', mic_locations');

%% Assignment 3: Kalman filtering using "position measurements"
A = eye(2);
C = eye(2);
Q3 = eye(2)*2e-7;

kf_pos = nan(2, N_exp+1);
kf_P = nan(2, 2, N_exp+1);

kf_pos(:,1) = [0.1 0.6]';
kf_P(:,:,1) = eye(2);

for k = 2:(N_exp+1)
    K = A*kf_P(:,:,k-1)*C'/(C*kf_P(:,:,k-1)*C' + diag(diagP(1:2,k-1))*0.5);
    kf_pos(:,k) = (A - K*C)*kf_pos(:,k-1) + K*theta_hat(1:2,k-1);
    kf_P(:,:,k) = A*kf_P(:,:,k-1)*A' + Q3 - K*C*kf_P(:,:,k-1)*A';
end

figure
kf_pos(:,1) = [];
kf_P(:,:,1) = [];
fig4 = plotresults(kf_pos, dimdiag(kf_P, 3)', mic_locations');

%% Assignment 4: Extended Kalman filter using TOA measurements
% States = [x, y, tau]
F = eye(3);
delta_tau = mean(diff(TOA_true));
var_tau = var(diff(TOA_true));
Q4 = blkdiag(eye(2)*3e-7, var_tau);
R4 = 1.8*diag(mic_variance);

% State 'dynamics'
f_ext = @(x) [x(1), x(2), x(3) + delta_tau]';
% f_ext_2 = @(xp, xpp) [xp + (xp - xpp); ...
%                       xp(3) + delta_tau]; % Cool idea, but the state must
%                       be extended to do this

% Measurement model (same as f from Assignment 2)
h_ext = @(x) f(x, mic_locations); % For readability

ekf_x = nan(3, N_exp+1);
ekf_P = nan(3, 3, N_exp+1);

ekf_x(:,1) = [0.1 0.6 0]';
ekf_P(:,:,1) = eye(3);

for k = 2:(N_exp+1)
    % Time update
    ekf_x(:,k) = f_ext(ekf_x(:,k-1));
    ekf_P(:,:,k) = F*ekf_P(:,:,k-1)*F' + Q4;
    
    % Measurement update
    H = Jacobian(ekf_x(:,k), mic_locations);
    K = ekf_P(:,:,k)*H'/(H*ekf_P(:,:,k)*H' + R4);
    ekf_x(:,k) = ekf_x(:,k) + K*(y_bias_correct(k-1,:)' - h_ext(ekf_x(:,k)));
    ekf_P(:,:,k) = ekf_P(:,:,k) ...
               - ekf_P(:,:,k)*H'*((H*ekf_P(:,:,k)*H' + R4)\H*ekf_P(:,:,k));
end

figure
ekf_x(:,1) = [];
ekf_P(:,:,1) = [];
fig5 = plotresults(ekf_x(1:2,:), ...
                         dimdiag(ekf_P(1:2,1:2,:), 3)', mic_locations');

%% Functions
function [theta, diagP, F] = nls(yk, stds, th_hat0, maxiter, mic_locations)
    % NLS algorithm
    converged = false;
    max_iter_reached = false;

    theta = th_hat0;
    count = 1;

    L = diag(stds);
    W = inv(L*L');  % Assume uncorrelated errors

    while ~converged && ~max_iter_reached
        % Solve a linear least-squares model
        F = Jacobian(theta, mic_locations);
        yhat = f(theta, mic_locations);
        
        % Unweighted solution
%         d_theta = F\(yk - yhat);

        % Weighted LSQ solution (BLUE)
        d_theta = (F'*W*F)\F'*W*(yk - yhat);

        % Check convergence
        if norm(d_theta) < 1e-6
            converged = true;
        elseif count > maxiter
            max_iter_reached = true;
%             warning('Max iterations reached')
%             disp(theta)
        end
        
        count = count + 1;
        theta = theta + d_theta;
    end
    
    % Variance of the BLUE
    P = inv(F'*W*F);
    diagP = diag(P);
end

function dF = Jacobian(theta, mic_locations)
    c = 343; % speed of sound in [m/s]

    norm_dist = vecnorm([theta(1) theta(2)] - mic_locations, 2, 2);

    dF = [(theta(1) - mic_locations(:,1))./norm_dist/c, ...
        (theta(2) - mic_locations(:,2))./norm_dist/c, ...
        ones(size(mic_locations, 1), 1)];
end

function ftheta = f(theta, mic_locations)
    c = 343; % speed of sound in [m/s]

    % Normalize for cm (????)
    ftheta = (theta(3) ...
        + vecnorm([theta(1) theta(2)] - mic_locations, 2, 2)/c);
    % Slightly verbose notation to avoid issues with theta being a column
    % vector
end

function res = dimdiag(A, dim)
    A_size = size(A);
    N = A_size(dim); A_size(dim) = [];
    res = nan(A_size(1));
    A = shiftdim(A, dim);
        
    for i = 1:N
        res(:,i) = diag(squeeze(A(:,:,i)));
    end
end