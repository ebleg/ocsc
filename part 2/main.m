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
settings = optimoptions('fmincon', 'Algorithm', 'interior-point');
[u1, fval, exitflag, output] = fmincon(@(u) J(u, x0, par), ...
                                           u0_1, [], [], [], [], lb, ub, [], settings);

[u2, fval, exitflag, output] = fmincon(@(u) J(u, x0, par), ...
                                           u0_2, [], [], [], [], lb, ub, [], settings);

%% Plot the solutions
close all;
figure
[~, x_0] = J(repmat(30, 60, 1), x0, par);
[~, x_1] = J(u1, x0, par);
[~, x_2] = J(u2, x0, par);
t = 1:60;
subplot(421); hold on;
plot(t, x_0(1,:));
plot(t, x_1(1,:));
plot(t, x_2(1,:));
subplot(422); hold on;
plot(t, x_0(2,:));
plot(t, x_1(2,:));
plot(t, x_2(2,:));
subplot(423); hold on;
plot(t, x_0(3,:));
plot(t, x_1(3,:));
plot(t, x_2(3,:));
subplot(424); hold on;
plot(t, x_0(4,:));
plot(t, x_1(4,:));
plot(t, x_2(4,:));

subplot(425); hold on;
plot(t, x_0(5,:));
plot(t, x_1(5,:));
plot(t, x_2(5,:));
subplot(426); hold on;
plot(t, x_0(6,:));
plot(t, x_1(6,:));
plot(t, x_2(6,:));
subplot(427); hold on;
plot(t, x_0(7,:));
plot(t, x_1(7,:));
plot(t, x_2(7,:));
subplot(428); hold on;
plot(t, x_0(8,:));
plot(t, x_1(8,:));
plot(t, x_2(8,:));


