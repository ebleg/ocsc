%% ------------------------------------------------------------------------
%  Filtering and Identification assignment 1
%  ------------------------------------------------------------------------
%
%    Fères Hassan (4362152) & Emiel Legrand (4446100)
%
%    December 10, 2020
% -------------------------------------------------------------------------

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
% figure
% hist_microphone = 1;
% histogram(mic_error(:,hist_microphone), 6, 'Normalization', 'pdf');
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
    [theta_hat(:,k), diagP(:,k-1)] = nls(y_bias_correct(k-1,:)', ...
        sqrt(mic_variance), ...
        theta_hat(:,k-1), ...
        100, mic_locations);
end
theta_hat(:,1) = []; % Remove the initial guess to obtain the correct length

% Question 2c: Visualize the results
figure
plotresults(theta_hat(1:2,:), diagP(1:2,:)', mic_locations');

%% Assignment 3


%% A4

%% Functions
function [theta, diagP] = nls(yk, stds, th_hat0, maxiter, mic_locations)
    % NLS algorithm
    converged = false;
    max_iter_reached = false;

    theta = th_hat0;
    count = 1;
    P = eye(numel(th_hat0))*1e-6;

%     W = stds*stds';
    W = diag(stds.^2);

    while ~converged && ~max_iter_reached
        % Solve a linear least-squares model
        F = Jacobian(theta, mic_locations);
        yhat = f(theta, mic_locations);
        d_theta = F\(yk - yhat);

        P = (F'*F)\F'*W*F/(F'*F);

        % Check convergence
        if norm(d_theta) < 1e-6
            converged = true;
        elseif count > maxiter
            max_iter_reached = true;
            %            warning('Max iterations reached')
            %            disp(theta)
        end
        

        count = count + 1;
        theta = theta + d_theta;
    end

%     theta = lsqnonlin(@(x) yk - f(x, mic_locations), th_hat0);

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
