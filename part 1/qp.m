%% Determine custom parameters
clear; clc; close all;

id_fer = '4362152';
id_emiel = '4446100';

E3 = str2double(id_emiel(end)) + str2double(id_fer(end));
E2 = str2double(id_emiel(end-1)) + str2double(id_fer(end-1));
E1 = str2double(id_emiel(end-2)) + str2double(id_fer(end-2));

%% First question
data = readtable('measurements.csv');
delta_t = 3600;

q_dot_occ = table2array(data(:, 'q_dot_occ'));
q_dot_ac = table2array(data(:, 'q_dot_ac'));
q_dot_vent = table2array(data(:, 'q_dot_vent'));
T_amb = table2array(data(:, 'T_amb'));
T_b = table2array(data(:, 'T_b'));
q_dot_solar = table2array(data(:, 'q_dot_solar'));
cost = table2array(data(:, 'Phi'))/3600;

phi = [q_dot_solar, ...
       q_dot_occ + q_dot_ac - q_dot_vent, ...
       T_amb - T_b]*delta_t;

y = diff(T_b);
phi_org = phi;
phi(end, :) = []; % Remove to ensure compatibility with diff array

H = 2*(phi'*phi);
c = -2*phi'*y;

qpoptions = optimoptions('quadprog', 'Display', 'none', 'Algorithm', 'interior-point-convex', 'MaxIter', 2000);
[a, error, flag1, ~] = quadprog(H, c, [], [], [], [], [], [], [], qpoptions);

plot(T_b(1:(end-1)) + phi*a); hold on;
plot(T_b);

%% Second question
% Parameters
q_ac_max = 100;
N = 2160;
T_min = 15;  
T_max = 28;
T0 = 22.43;
Tref = 22;
alpha = 0.1 + E2/10;

% Cost function
H = [eye(N) zeros(N); zeros(N, 2*N)]*2*alpha;
c = [repmat(-2*alpha*Tref, N, 1); cost(1:N)*delta_t];

% Constraints
Aeq = [eye(N) zeros(N)] ...
    - [[zeros(1, N); [eye(N-1)*(1 - delta_t*a(3)),  zeros(N-1, 1)]] eye(N)*a(2)*delta_t];
beq = (a(1)*q_dot_solar(1:N) + a(2)*(q_dot_occ(1:N) - q_dot_vent(1:N)) + a(3)*T_amb(1:N))*delta_t;
beq(1) = beq(1) + (1 - a(3))*T0;

filt = [not(q_dot_occ(1:N) > 0); false(N, 1)];
ub = [repmat(T_max, N, 1); repmat(q_ac_max, N, 1)];
lb = [repmat(T_min, N, 1); zeros(N, 1)];
ub(filt) = inf;
lb(filt) = -inf;

[x, fval, flag2, ~] = quadprog(H, c, [], [], Aeq, beq, lb, ub, [], qpoptions);

%% Visualization
t = (1:N)/24;
close all;
subplot(411)
plot(t, x(1:N)); hold on;
plot(t, repmat(T_min, N, 1), 'k', 'LineWidth', 1.5);
plot(t, repmat(T_max, N, 1), 'k', 'LineWidth', 1.5);
title('Inside temperature');
xlabel('t [days]');
ylabel('T [$^\circ$C]');

subplot(412)
plot(t, T_amb);
title('Ambient temperature');
xlabel('t [days]');
ylabel('T [$^\circ$C]');

subplot(413)
plot(t, x((N+1):(2*N))); hold on;
title('AC Power');
xlabel('t [days]');
ylabel('$\dot{q}_\mathrm{ac}$ [W]');

cumul_cost = zeros(1,N);
q = x((N+1):(2*N));
Tb = x(1:N);
for i = 1:N
   cumul_cost(i) = cost(1:i)'*q(1:i)*delta_t + alpha*(Tb(1:i) - Tref)'*(Tb(1:i) - Tref);
end
subplot(414)
plot(t, cumul_cost);
title('Cumulative cost')
xlabel('t [days]');
ylabel('Cost [EUR]');
