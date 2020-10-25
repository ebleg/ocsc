% addpath('bfgs');
% 
% close all; clear; clc;
% 
% run parameters
% run defaultPlotSettings
% 
% u0 = repmat(15, 60, 1);
% x0 = zeros(8, 1);
% lb = repmat(15, 60, 1);
% ub = repmat(45, 60, 1);
% tic
% % % Time the cost function for reference
% % t1 = toc;
% % for i = 1:50
% %     x = J(u0, x0, par);
% % end
% % t2 = toc;
% % fprintf('Average cost function evaluation time: %fs\n', (t2-t1)/1000);
% %     
% %% Default SQP solution
% sqp_sol = struct();
% sqp_sol.settings = optimoptions(@fmincon, 'Algorithm', 'sqp', ...
%                     'MaxFunctionEvaluations', 10000, 'Display', 'off');
%                 
% [sqp_sol.u, sqp_sol.fval, sqp_sol.exitflag] = ...
%     fmincon(@(u) J(u, x0, par), u0, [], [], [], [], lb, ub, [], sqp_sol.settings);
%                                        
% fprintf('SQP solution --> total cost %f\n', sqp_sol.fval);
% 
% %% Interior-point solution
% ip_sol = struct();
% ip_sol.settings = optimoptions(@fmincon, 'Algorithm', 'interior-point', ...
%                          'MaxFunctionEvaluations', 10000, 'Display', 'off');
%                      
% [ip_sol.u, ip_sol.fval, ip_sol.exitflag] = ...
%     fmincon(@(u) J(u, x0, par), u0, [], [], [], [], lb, ub, [], ip_sol.settings);
%                                        
% fprintf('Interior-point solution --> total cost %f\n', ip_sol.fval);
% 
% %% Simulated annealing
% sa_sol = struct();
% sa_sol.settings = optimoptions(@simulannealbnd, ...
%           'MaxFunctionEvaluations', 30000, 'Display', 'off');
% 
% [sa_sol.u, sa_sol.fval, sa_sol.exitflag] = ...
%     simulannealbnd(@(u) J(u, x0, par), u0, lb, ub);
% 
% fprintf('Simulated annealing solution --> total cost %f\n', sa_sol.fval);
% 
% %% Integer Genetic Algorithm 
% lb = ones(1, 60);      % Lower bound representing the first discrete values in the set
% ub = repmat(7, 60, 1); % Number of discrete values in set
% 
% ga_sol = struct();
% ga_sol.settings = optimoptions(@ga, ...
%                     'PopulationSize', 150, ...
%                     'MaxGenerations', 600, ...
%                     'EliteCount', 10, ...
%                     'FunctionTolerance', 1e-8, ...
%                     'PlotFcn', @gaplotbestf);
%                                 
% %rng(0, 'twister');
% %[ga_sol.u, ga_sol.fval, ga_sol.exitflag] = ga(@(u)J(u, x0, par), ...  
% %    60, [], [], [], [], lb, ub,[],[], ga_sol.settings);
% 
% rng(0, 'twister');
% [ga_sol.u, ga_sol.fval, ga_sol.exitflag] = ga(@(u)J_disc(u, x0, par), ...    
%     60, [], [], [], [], lb, ub, [],1:60 , ga_sol.settings);
% 
% ga_sol.u_map = mapvariables(ga_sol.u);
% 
% fprintf('Integer Genetic Algorithm solution --> total cost %f\n', ga_sol.fval):

%% Plot results
close all; figure;
plot_traffic_simulation(gcf, sa_sol.u, x0, 'Simulated annealing', par);
plot_traffic_simulation(gcf, sqp_sol.u, x0, 'SQP', par);
plot_traffic_simulation(gcf, ip_sol.u, x0, 'Interior-point', par);
plot_traffic_simulation(gcf, ga_sol.u, x0, 'Genetic Algorithm', par);
