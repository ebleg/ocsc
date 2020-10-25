addpath('bfgs');

close all; clear; clc;

run parameters
run defaultPlotSettings

discrete_set = (15:5:45);

u0 = repmat(15, 60, 1);
x0 = zeros(8, 1);
lb = repmat(1, 60, 1);% Lower bound representing the first discrete values in the set
ub = repmat(7, 60, 1);% Number of discrete values in set

%% integer optimization

surr_sol = struct()

options.MinSurrogatePoints = 50;

surr_sol.setting = optimoptions('surrogateopt','CheckpointFile','C:\TEMP\checkfile.mat','PlotFcn',[]);
rng default % For reproducibility
objfun = @(u)J_disc(u, x0, par);
[surr_sol.u, surr_sol.Fval] = surrogateopt(objfun,lb,ub,1:60,surr_sol.setting);