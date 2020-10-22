%% ----------------------------------------------------------------------------
%  LINEAR AND QUADRATIC PROGRAMMING ASSIGNMENT
%  ----------------------------------------------------------------------------
% 
%    Fères Hassan (4362152) & Emiel Legrand (4446100)
% 
%    October 7, 2020
%
% ------------------------------------------------------------------------------

%% Clear workspace and plot settings
clear; clc; close all;

% Default plot settings
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultAxesXGrid','on');
set(groot, 'defaultAxesYGrid','on');
set(groot, 'defaultAxesZGrid','on');
set(groot, 'defaultAxesXMinorGrid','on');
set(groot, 'defaultAxesYMinorGrid','on');
set(groot, 'defaultAxesZMinorGrid','on');
set(groot, 'defaultAxesTitleFontSizeMultiplier', 1.5);
set(groot, 'defaultAxesXMinorGridMode','manual');
set(groot, 'defaultAxesYMinorGridMode','manual');
set(groot, 'defaultAxesYMinorGridMode','manual');

%% Determine custom parameters
id_fer = '4362152';
id_emiel = '4446100';

E3 = str2double(id_emiel(end)) + str2double(id_fer(end));
E2 = str2double(id_emiel(end-1)) + str2double(id_fer(end-1));
E1 = str2double(id_emiel(end-2)) + str2double(id_fer(end-2));

%% Options for optimization routines
lpoptions = optimoptions('linprog', 'Algorithm', 'interior-point', 'Display', 'none');
qpoptions = optimoptions('quadprog', 'Display', 'none', 'Algorithm', 'interior-point-convex', 'MaxIter', 2000);

%% Question 1a - b
fprintf('Question 1a - Installation plan for maximal power\n')

total_budget = 24000 + 300*E1;

% LP problem setup
f = -[4 2.5];
A = [3000 1500; 1 1];
b = [total_budget; 12];

% Perform optimization with linprog
[x, fval, flag, output] = linprog(f, A, b, [], [], [0 0], [], lpoptions);

if flag == 1 % Optimization succesful
    fprintf('Optimal solution found\nX = %f Y = %f\n', x(1), x(2)) 
    fprintf('Check closest integer solutions\n')
    % Check the closest integer solution (heuristic method to find the optimum)
    [X_int, Y_int] = meshgrid([floor(x(1)), ceil(x(1))], [floor(x(2)), ceil(x(2))]);
    X_int = X_int(:); Y_int = Y_int(:);
    feasible = all(A*[X_int Y_int]' <= b);
    for i = 1:numel(X_int)
    fprintf('X = %d  Y = %d => f = %f  feasible: %d\n', ...
             X_int(i)', Y_int(i)', f*[X_int(i) Y_int(i)]', feasible(i));
    end
    [~, minidx] = min(f*[X_int(feasible) Y_int(feasible)]');
    fprintf('OPTIMAL SOLUTION : X = %d  Y = %d\n\n', X_int(minidx), Y_int(minidx));

    % Plot the analysis
    [X1, X2] = meshgrid(0:12, 0:12);
    Y = zeros(size(X1));

    % Objective function
    for i = 1:numel(X1)
        Y(i) = f*[X1(i) X2(i)]';
    end
    
    figure; hold on; grid; grid minor;
    contourf(X1, X2, Y, 30, 'LineWidth', 0.5, 'HandleVisibility', 'off')
    xrange = 0:10;
    plot(xrange, 12-xrange, 'k', 'Linewidth', 1.8, 'HandleVisibility', 'off')
    plot(xrange, total_budget/1500-2*xrange, 'k', 'Linewidth', 1.8, 'HandleVisibility', 'off')
    [Xdots, Ydots] = meshgrid(1:10, 1:10);
    scatter(Xdots(:), Ydots(:), 'k', 'filled', 'HandleVisibility', 'off')
    scatter(x(1), x(2), 'b', 'filled', 'DisplayName', 'LP solution')
    scatter(X_int(minidx), Y_int(minidx), 'r', 'filled', 'DisplayName', 'Integer solution')
    xlabel('X'); ylabel('Y'); title('Analysis of integer solutions')
    colorbar; legend; xlim([0 10]); ylim([1 10]);
else
   fprintf('Optimization routine failed (flag %d)\n%s\n', flag, output.message)
end

%clearvars -except A b f E1 E2 E3 lpoptions qpoptions

%% Question 1c - Find the number of durable years
total_inst_budget = 24000 + 300*E1;

% Array for the maintenance costs per machine
X1 = [200 200 200 300 300 400 500 600 700 800];
X2 = [1:5 5 5 5 5 5];
Xcost = X1 + X2*E2;

Y1 = [50 50 100 150 150 200 250 300 350 400];
Ycost = Y1 + X2*E3;
clear X1 X2 Y1

% Solve LP problem for every year
maintenance_budget = 4000 + 100*E1;
x_arr = nan(length(Xcost),length(A));

% Loop over all possible durable years and save the x and fval values
for yrs = 1:10  
    A = [1 1; 3000 + sum(Xcost(1:yrs)) 1500 + sum(Ycost(1:yrs))];
    b = [12; total_inst_budget + yrs*maintenance_budget];
    
    % Perform the optimization
    [x_arr(yrs,:), fval, flag, ~] = linprog(f, A, b, [], [], [0 0], [], lpoptions);
    x_arr_int = round(x_arr(yrs,:));   

    if flag ~= 1
        error('Infeasible problem')
    end
end

% Compute the number of feasible years
power_max = -f*x_arr'; % Find theoretical maximum power
power_int = -f*round(x_arr)';  % Find the maximum power for integer values
filt = power_max >= power_int;
filt = filt & (power_int == max(power_int(filt))); 
fprintf('Question 1c - Number of feasible years\n');
fprintf('FEASIBLE YEARS = %d\n\n', find(filt, 1, 'last'));

% Visualization
figure; hold on; grid;
xrange = 1:10;
plot(xrange,power_int, 'g', 'Linewidth', 1.3);
plot(xrange,power_max, 'k', 'Linewidth', 1.3);
legend('Integer max power','Theoretical max power');
xlabel('Durable years'); ylabel('Power [kW]'); 
title('Analysis of max power','FontWeight','Normal','FontSize',15);
hold off;

clearvars -except E1 E2 E3 qpoptions

%% Question 3
fprintf('Question 3 - Parameter identification\n');
data = readtable('measurements_physical.csv', 'PreserveVariableNames', true);
delta_t = 3600;

% Construct arrays from the table data
q_dot_occ = table2array(data(:, 1));
q_dot_ac = table2array(data(:, 2));
q_dot_vent = table2array(data(:, 3));
q_dot_solar = table2array(data(:, 4));
T_amb = table2array(data(:, 5'));
T_b = table2array(data(:, 6));
cost = table2array(data(:, 7))/3600; % Make sure the units match
N = numel(T_b);

phi = [q_dot_solar, ...
       q_dot_occ + q_dot_ac - q_dot_vent, ...
       T_amb - T_b]*delta_t;

y = diff(T_b);
phi_org = phi;
phi(end, :) = []; % Remove to ensure compatibility with diff array

H = 2*(phi'*phi);
c = -2*phi'*y;

[a, error, flag1, ~] = quadprog(H, c, [], [], [], [], [], [], [], qpoptions);

t = (1:N)/24;
figure;
plot(t(1:(end-1)), T_b(1:(end-1)) + phi*a - T_b(1:(end-1)), 'LineWidth', 1.25)
title('Prediction error')
xlabel('t [days]'); ylabel('error [$^\circ$C]');
fprintf('MODEL COEFFICIENTS:\n\t a1 = %e, a2 = %e a3 = %e\n', a(1), a(2), a(3));
if abs(1 - delta_t*a(3)) < 1
   fprintf('The resulting model is stable\n');
else
    fprintf('The resulting model is unstable\n');
end

%% Question 4
fprintf('\nQuestion 4 - Optimal control problem\n')

% Parameters
q_ac_max = 100;
T_min = 15;  
% T_max = 28;
T_max = 30; % Relaxed constraint because of infeasibility
T0 = 22.43;
Tref = 22;
alpha = 0.1 + E2/10; % For readability

% Cost function
H = [eye(N) zeros(N); zeros(N, 2*N)]*2*alpha;
c = [repmat(-2*alpha*Tref, N, 1); cost(1:N)*delta_t];

% Establish the constraints
Aeq = [eye(N) zeros(N)] ...
    - [[zeros(1, N); [eye(N-1)*(1 - delta_t*a(3)),  zeros(N-1, 1)]] eye(N)*a(2)*delta_t];
beq = (a(1)*q_dot_solar(1:N) + a(2)*(q_dot_occ(1:N) - q_dot_vent(1:N)) + a(3)*T_amb(1:N))*delta_t;
beq(1) = beq(1) + (1 - a(3))*T0; % Take the initial condition into account

filt = [not(q_dot_occ(1:N) > 0); false(N, 1)];
ub = [repmat(T_max, N, 1); repmat(q_ac_max, N, 1)]; % Upper bounds
lb = [repmat(T_min, N, 1); zeros(N, 1)]; % Lower bounds
ub(filt) = inf;
lb(filt) = -inf;

[x, fval, flag2, ~] = quadprog(H, c, [], [], Aeq, beq, lb, ub, [], qpoptions);

figure;
% Visualization
subplot(411) % Plot the inside temperature
plot(t, x(1:N)); hold on;
plot(t, repmat(T_min, N, 1), 'k', 'LineWidth', 1.5);
plot(t, repmat(T_max, N, 1), 'k', 'LineWidth', 1.5);
title('Inside temperature');
xlabel('t [days]');
ylabel('T [$^\circ$C]');

subplot(412)
plot(t, T_amb); % Plot the ambient temperature
title('Ambient temperature');
xlabel('t [days]');
ylabel('T [$^\circ$C]');

subplot(413) % Plot q_dot_ac
plot(t, x((N+1):(2*N))); hold on;
title('AC Power');
xlabel('t [days]');
ylabel('$\dot{q}_\mathrm{ac}$ [W]');

cumul_cost = zeros(1,N); % Plot cumulative cost
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

fprintf('TOTAL COST: %f EUR\n', cumul_cost(end));