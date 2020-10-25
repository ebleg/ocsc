addpath('bfgs');

close all; clear; clc;

run parameters
run defaultPlotSettings

u0 = repmat(15, 60, 1);
x0 = zeros(8, 1);
lb = repmat(1, 60, 1);% Lower bound representing the first discrete values in the set
ub = repmat(7, 60, 1);% Number of discrete values in set

%% Integer Genetic Algorithm 
ga_sol = struct();
ga_sol.settings = optimoptions(@ga, ...
                    'PopulationSize', 150, ...
                    'MaxGenerations', 300, ...
                    'EliteCount', 10, ...
                    'FunctionTolerance', 1e-8, ...
                    'PlotFcn', @gaplotbestf);
                                
%rng(0, 'twister');
%[ga_sol.u, ga_sol.fval, ga_sol.exitflag] = ga(@(u)J(u, x0, par), ...         % nog niet werkend.
%    60, [], [], [], [], lb, ub,[],[], ga_sol.settings);


rng(0, 'twister');
[ga_sol.u, ga_sol.fval, ga_sol.exitflag] = ga(@(u)J_disc(u, x0, par), ...         % nog niet werkend.
    60, [], [], [], [], lb, ub, [],1:60 , ga_sol.settings);

ga_sol.u_map = mapvariables(ga_sol.u);
