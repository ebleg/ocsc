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

c = 343; % speed of sound [m/s]
k = 59; % number of pulses
m = 7; % number of microphones

% y = [y zeros(size(y,1),1)]; % Append column of zero's for speed

% % Calculate measurement errors and estimate biases
% for i = 1:k
%     y(i,8) = mean(y(i,1:7)); % Average TOA for each time step
%     for j = 1:m
%         error(i,j) = y(i,j) - y(i,8); % Measurement errors at each timestep
%         bias(j) = mean(error(:,j)); % Mean error of all timesteps to estimate bias
%     end
% end
% 
% % Estimate variance
% for j = 1:m
%     for i = 1:k
%         efb(i,j) = error(i,j)-bias(j); % Measurement noise
%     end
%     variance(j) = var(efb(:,j)); % Variance of measurement noise
% end
% 
% % Visualize measurement errors of microphone 1
% H = histogram(efb(:,5),6);

% -------------------------------------------------------------------------
% Alternatively, without loops

% Question 1a: Error per microphone
y_calib = y; % This will be overwritten later by new data
TOA_true = mean(y_calib, 2);
mic_error = y_calib - repmat(TOA_true, 1, m);

% Question 1b: Microphone bias
mic_bias = mean(mic_error);

% Question 1c: Microphone variance
mic_variance = var(mic_error);

% Question 1d: Visualization with histogram
figure
hist_microphone = 5;
histogram(mic_error(:,hist_microphone), 6, 'Normalization', 'pdf');
title(sprintf('Error of microphone %d', hist_microphone))

figure
bar(mic_bias)

%% Assignment 2: Nonlinear least-squares
load('experiment.mat')

% Question 2a
y_exp = y;
N_exp = size(y_exp, 1);
y_bias_correct = y_exp - repmat(mic_bias, N_exp, 1);

%%
maxiter = 100;
% for k = 1:size(y,1) 
%     % initial estimate of theta
%     th_hat0 = [10 60 0];
%     % non linear LS estimate 
%     [th_hat(k,:),diagP(k,:)] = nls(yk,stds,th_hat0,maxiter,mic_locations);
% end


% A3

%% A4

%% Functions
function [theta, diagP] = nls(yk, stds, th_hat0, maxiter, mic_locations)
    % NLS algorithm
    converged = false;
    max_iter_reached = false;
    
    theta = th_hat0;
    count = 1;
    
    while ~converged && ~max_iter_reached       
       % Solve a linear least-squares model
       theta_new = Jacobian(theta, mic_locations)\yk;       
       
       % Check convergence
       if norm(theta_new - theta) < 1e-6
           converged = true;
       elseif count > maxiter
           max_iter_reached = true;
       end
       
       count = count + 1;
       theta = theta_new;
       
    end 
end

function dF = Jacobian(theta, mic_locations)
    c = 343; % speed of sound in [m/s]
    
    norm_dist = vecnorm([theta(1) theta(2)] - mic_locations, 2, 2)*100/c;
    
    dF = [(theta(1) - mic_locations(:,1))./norm_dist, ...
          (theta(2) - mic_locations(:,2))./norm_dist, ...
          ones(size(mic_locations, 1), 1)];
end

function ftheta = f(theta, mic_locations)
    c = 343; % speed of sound in [m/s]
    
    % Normalize for cm (????)
    ftheta = theta(3) ...
             + vecnorm([theta(1) theta(2)] - mic_locations, 2, 2)*100/c;
    % Slightly verbose notation to avoid issues with theta being a column
    % vector
end
