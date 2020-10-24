addpath('bfgs');

close all; clear; clc;

run parameters

u0_1 = repmat(15, 60, 1);
u0_2 = repmat(45, 60, 1);
x0 = zeros(8, 1);
lb = repmat(15, 60, 1);
ub = repmat(45, 60, 1);

%% Time the cost function for reference
% t1 = toc;
% for i = 1:50
%     x = J(u0, x0, par);
% end
% t2 = toc;
% fprintf('Average cost function evaluation time: %fs\n', (t2-t1)/1000); 
    
%% Default fmincon solution
settings = optimoptions('fmincon', 'Algorithm', 'sqp');
[u1, fval, exitflag, output] = fmincon(@(u) J(u, x0, par), ...
                                           u0_1, [], [], [], [], lb, ub, [], settings);

