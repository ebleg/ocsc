addpath('custom bfgs');

close all; clear; clc;

run parameters

u0 = repmat(30, 60, 1);
x0 = zeros(8, 1);
lb = repmat(15, 60, 1);
ub = repmat(45, 60, 1);

%% Time the cost function for reference
t1 = toc;
for i = 1:50
    x = J(u0, x0, par);
end
t2 = toc;
fprintf('Average cost function evaluation time: %fs\n', (t2-t1)/1000);

%% Default fmincon solution
tic
[u_opt, fval, exitflag, output] = fmincon(@(u) J(u, x0, par), u0, [], [], [], [], lb, ub);
toc

%% Custom bfgs / line search solution
% u_opt = bfgs(@(u) J(u, x0, par), u0);
