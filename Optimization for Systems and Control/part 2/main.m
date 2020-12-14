%% ------------------------------------------------------------------------
%  NONLINEAR PROGRAMMING ASSIGNMENT
%  ------------------------------------------------------------------------
% 
%    Fères Hassan (4362152) & Emiel Legrand (4446100)
% 
%    October 27, 2020
% -------------------------------------------------------------------------

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

% Group specific parameters
id_fer = '4362152';
id_emiel = '4446100';

par.E3 = str2double(id_emiel(end)) + str2double(id_fer(end));
par.E2 = str2double(id_emiel(end-1)) + str2double(id_fer(end-1));
par.E1 = str2double(id_emiel(end-2)) + str2double(id_fer(end-2));
clear id_fer id_emiel;

% Parameters for the first link
par.ud = struct();
par.ud.c = 60;
par.ud.N_lane = 3;
par.ud.v_free = 50/3.6;
par.ud.l_veh = 7;
par.ud.l = 1000;
par.ud.beta = [0.33 0.34 0.33]';
par.ud.mu = [1600 1800 1500]'/3600;
par.ud.Cp = par.ud.N_lane*par.ud.l/par.ud.l_veh;

% Parameters for the second link
par.o1d = struct();
par.o1d.c = 60;
par.o1d.N_lane = 3;
par.o1d.v_free = 60/3.6;
par.o1d.l_veh = 7;
par.o1d.l = 1000;
par.o1d.beta = [0.33 0.34 0.33]';
par.o1d.mu = [1600 1800 1500]'/3600;
par.o1d.Cp = par.o1d.N_lane*par.o1d.l/par.o1d.l_veh;

%% Optimization parameters
u0 = repmat(45, 60, 1);
x0 = zeros(8, 1);
lb = repmat(15, 60, 1);
ub = repmat(45, 60, 1);

<<<<<<< HEAD
% %% Time the cost function for reference
% t1 = toc;
% for i = 1:50
%     x = J(u0, x0, par);
% end
% t2 = toc;
% fprintf('Average cost function evaluation time: %fs\n', (t2-t1)/50);
% %%

%% Time the cost function for reference
t1 = toc;
for i = 1:50
    x = J(u0, x0, par);
end
t2 = toc;
fprintf('Average cost function evaluation time: %fs\n', (t2-t1)/50);

%% Default SQP solution
tic
sqp_sol = struct();
sqp_sol.settings = optimoptions(@fmincon, 'Algorithm', 'sqp', ...
                    'MaxFunctionEvaluations', 10000, 'Display', 'off');
                
[sqp_sol.u, sqp_sol.fval, sqp_sol.exitflag] = ...
    fmincon(@(u) J(u, x0, par), u0, [], [], [], [], lb, ub, [], sqp_sol.settings);
                                       
fprintf('SQP solution --> total cost %f\n', sqp_sol.fval);
toc
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

%% Simulated annealing
sa_sol = struct();
sa_sol.settings = optimoptions(@simulannealbnd, ...
          'MaxFunctionEvaluations', 30000, 'Display', 'off');

[sa_sol.u, sa_sol.fval, sa_sol.exitflag] = ...
    simulannealbnd(@(u) J(u, x0, par), u0, lb, ub);

fprintf('Simulated annealing solution --> total cost %f\n', sa_sol.fval);

%% Integer Genetic Algorithm 
lb = ones(1, 60);      % Lower bound representing the first discrete value in the set
ub = repmat(7, 60, 1); % Number of discrete values in set
 
ga_sol = struct();
ga_sol.settings = optimoptions(@ga, 'PopulationSize', 150, 'MaxGenerations', 600, ...
         'EliteCount', 10, 'FunctionTolerance', 1e-8, 'PlotFcn', @gaplotbestf);                                 

rng(0, 'twister');
[ga_sol.u, ga_sol.fval, ga_sol.exitflag] = ga(@(u)J_disc(u, x0, par), ...    
     60, [], [], [], [], lb, ub, [],1:60 , ga_sol.settings);
 
ga_sol.u_map = mapvariables(ga_sol.u);
 
 fprintf('Integer genetic algorithm solution --> total cost %f\n', sa_sol.fval);
 
%% Plot green light time

t = 1:60;
%sqp_sol_45 = load('nlp_results\sqp_sol_45.mat');
%sqp_sol_15 = load('nlp_results\sqp_sol_15.mat');
%sa_sol_45 = load('nlp_results\sa_sol_45.mat');
%sa_sol_15 = load('nlp_results\sa_sol_15.mat');
%ip_sol_45 = load('nlp_results\ip_sol_45.mat');
%ip_sol_15 = load('nlp_results\ip_sol_15.mat');
%ga_sol = load('nlp_results\ga_sol.mat');

% sqp_sol_45 = sqp_sol_45.sqp_sol;
% sqp_sol_15 = sqp_sol_15.sqp_sol;
% sa_sol_45 = sa_sol_45.sa_sol;
% sa_sol_15 = sa_sol_15.sa_sol;
% ip_sol_45 = ip_sol_45.ip_sol;
% ip_sol_15 = ip_sol_15.ip_sol;

figure;
hold on;

tile = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on;
title('Initialized at 15s')
plot(t, sa_sol_15.u, 'DisplayName', 'Simulated annealing');
plot(t, sqp_sol_15.u, 'DisplayName', 'Sequential Quadratic Programming');
plot(t, ip_sol_15.u, 'DisplayName', 'Interior-point');

nexttile; hold on;
title('Initialized at 45s')
plot(t, sa_sol_45.u, 'DisplayName', 'Simulated annealing');
plot(t, sqp_sol_45.u, 'DisplayName', 'Sequential Quadratic Programming');
plot(t, ip_sol_45.u, 'DisplayName', 'Interior-point');
plot(t, ga_sol.u_map, 'DisplayName', 'Genetic algorithm');
lgd = legend(gca);

lgd.Orientation = 'horizontal';
lgd.Layout.Tile = 'south';

title(tile, '\textbf{Solution comparison be different algorithms}', 'interpreter', 'latex');
ylabel(tile, 'Green light time [s]', 'interpreter', 'latex');
xlabel(tile, 'Time [min]', 'interpreter', 'latex');


%% Plot results
close all; figure;
% plot_traffic_simulation(gcf, sa_sol_15.u, x0, 'Simulated annealing (15s)', par);
% plot_traffic_simulation(gcf, sqp_sol_15.u, x0, 'SQP (15s)', par);
% plot_traffic_simulation(gcf, ip_sol_15.u, x0, 'Interior-point (15s)', par);
plot_traffic_simulation(gcf, sa_sol_15.u, x0, 'Simulated annealing (15s)', par);
plot_traffic_simulation(gcf, sqp_sol_15.u, x0, 'SQP (15s)', par);
plot_traffic_simulation(gcf, ip_sol_15.u, x0, 'Interior-point (15s)', par);
plot_traffic_simulation(gcf, repmat(30, 60, 1), x0, 'No control', par);
plot_traffic_simulation(gcf, ga_sol.u_map, x0, 'Genetic Algorithm', par);


%% plot system inputs 
t=1:60;
cp_ud_plot=zeros(3,60); cp_o1d_plot=zeros(3,60); alpha_ud_plot=zeros(1,60); alpha_o1d_plot=zeros(1,60);
for k=1:60
    cp_ud_plot(:,k)=Cp_ud_output(k, par);
    cp_o1d_plot(:,k)=Cp_o1d_output(k, par);
    alpha_ud_plot(1,k)=3600*alpha_ud_enter(k, par);
    alpha_o1d_plot(1,k)=3600*alpha_o1d_enter(k, par);
end

close all; figure;
tile = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on;
title('Capacity of downstream links','interpreter', 'latex')
plot(t, cp_ud_plot(1,:),'DisplayName','link $d,o_1$','lineWidth', .9);
plot(t, cp_ud_plot(2,:),'DisplayName','link $d,o_2$','lineWidth', .9);
plot(t, cp_ud_plot(3,:),'DisplayName','link $d,o_3$','lineWidth', .9);
plot(t, cp_o1d_plot(3,:),'DisplayName','link $d,u$','lineWidth', .9);
ylabel('Number of cars [-]', 'interpreter', 'latex');
legend('Location','northwest','Orientation','vertical')
ylim([5 58]);
grid off

nexttile; hold on;
title('Traffic flow of upstream link','interpreter', 'latex')
plot(t, alpha_ud_plot(1,:),'DisplayName','link $u,d$','lineWidth', .9);
plot(t, alpha_o1d_plot(1,:),'DisplayName','link $o_1 d$','lineWidth', .9);
ylabel('Vehicles / hour', 'interpreter', 'latex');
legend('Location','northwest','Orientation','vertical')
ylim([1780 2360])
grid off

title(tile, '\textbf{System intputs}', 'interpreter', 'latex');
xlabel(tile, 'Time [min]', 'interpreter', 'latex')

%% functions --------------------------------------------------------------

%% System inputs
function [out] = alpha_ud_enter(k, par)
% Rate of cars entering link ud
    if k <= 20
        out = (1800 + 10*par.E1)/3600;
    elseif k > 20 && k <= 40
        out = (2100 + 10*par.E2)/3600;
    else
        out = (2300 + 10*par.E3)/3600;
    end
end

function [out] = alpha_o1d_enter(k, par)
% Cars entering link o1d (constant in time)
    out = (2000 + 10*par.E1)/3600;
end

function [out] = Cp_ud_output(k, par)
% Capacity of output links for link ud
% CL = lane turning left , CS = straight, CR = turning right
    if k <= 20
        CL = 40 + par.E1;
    elseif k > 20 && k <= 35
        CL = 40 + par.E1 - 2*(k-20);
    elseif k > 35 && k <= 45
        CL = 10 + par.E1;
    else
        CL = 10 + par.E1 + 2*(k-45); 
    end
    
    CS = CL - par.E2;
    
    if k <= 30
        CR = 30 - par.E3;
    else
        CR = 30 + par.E3;
    end
    
    % Make sure the outputs are positive
    out = max([[CL CS CR]', zeros(3, 1)], [], 2);
end

function [out] = Cp_o1d_output(k, par)
% Capacity of the output links for link o1, d
    
    % Two capacities are identical to the ones computed 
    % for link ud
    tmp = Cp_ud_output(k, par);
    
    % Capacity for right link (ud)
    if k <= 30
        CR = 40 - par.E3;
    else
        CR = 40 + par.E3;
    end
    
    % Change the order (LEFT - STRAIGHT - RIGHT)
    % (o2 o3 ud)
    out = [tmp(2), tmp(3), CR]';
end


%% Cost functions 
function [x_new] = f(x, k, u, par)
% Overall state transition function, including both links ud and o1d.
    N = numel(x)/2;
    x_ud = x(1:N); % States for link ud
    x_o1d = x((N+1):end); % States for link 01d
    
    % Update states for link ud
    x_ud_new = state_transition(x_ud, k, u, ...
                                @(k) alpha_ud_enter(k, par), ...
                                @(k) Cp_ud_output(k, par), ...
                                par.ud);
    
    % Update states for link o1d
    x_o1d_new = state_transition(x_o1d, k, par.o1d.c - u, ...
                                 @(k) alpha_o1d_enter(k, par), ...
                                 @(k) Cp_o1d_output(k, par), ...
                                 par.o1d);
                             
    x_new = [x_ud_new; x_o1d_new];
end


function [x_new] = state_transition(x, k, u, alpha_enter_fcn, Co_fcn, par)
% alpha_enter_fcn = function for the entering cars
% C_fcn = function for the link capacity for a given
% Co_fcn = function (returning a vector) for the output link capacities
% x = [queue_left, queue_straight, queue_right, total # cars]

    % the parameter struct is link specific!

    % Compute tau and gamma
    tau = floor((par.Cp - [1 1 1 0]*x)*par.l_veh...
                 /par.N_lane/par.v_free/par.c);
    gamma = rem((par.Cp - [1 1 1 0]*x)*par.l_veh*par.N_lane/par.v_free, ...
                par.c);
            
    % Compute the cars arriving at the queues
    alpha_arrive = par.beta*[(par.c - gamma)/par.c, gamma/par.c]...
                           *[alpha_enter_fcn(k-tau); alpha_enter_fcn(k-tau-1)];
    
    % Update the state vector
    x_new = x + par.c*[eye(3); zeros(1,3)]*alpha_arrive ...
              + [0 0 0 par.c]'*alpha_enter_fcn(k) ...
              - par.c*[eye(3); ones(1,3)]...
                     *min([1/par.c*diag(par.mu)*[u u par.c]', ...
                           1/par.c*[eye(3) zeros(3,1)]*x + alpha_arrive, ...
                           Co_fcn(k)/par.c], [], 2);

end

function [y] = g(x, u, par)
% Output % number of vehicles on link (u,d) and link (o1,d)
    y = par.ud.c*[0 0 0 1 0 0 0 1]*x; % Both time steps are equal
end

function [cost, x_rec] = J(u, x0, par)
% Cost function
    T = numel(u); % the length of u determines the time horizon
    
    if nargout > 1
        x_rec = nan(numel(x0), numel(u));
    end
    
    x = x0; cost = 0;
    
    for k = 1:T
        x = f(x, k, u(k), par);
        cost = cost + g(x, u(k), par);
        
        if nargout > 1 % Don't do this if it's unnecessary to avoid overhead
            x_rec(:, k) = x;
        end
    end
end



%% Integer programming functions
function [cost, x_rec] = J_disc(u, x0, par)

% map the discrete variables
u = mapvariables(u);

[cost, x_rec] = J(u, x0, par);
end 


function u = mapvariables (u_disc)

% Discrete time set for green light time
discrete_set= (15:5:45);

% Map u from integer values used by GA to the discrete values required
u = discrete_set(u_disc);
end



%% Plot functions
function [ax] = plot_traffic_simulation(fig, u, x0, label, par)
    [~, x] = J(u, x0, par);
    
    t = 1:60;

    axes = flip(findall(fig,'type','axes'));
    
    titles = {'Link $ud$ - Queue for $o_1$', 'Link $o_1d$ - Queue for $o_2$', ...
              'Link $ud$ - Queue for $o_2$' 'Link $o_1d$ - Queue for $o_3$', ...
              '$n_{ud}$' '$n_{o_1d}$'};

    if isempty(axes) % No plots yet
        tiles = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
%       title('\textbf{State evolution over time}')
        for i = 1:6
           axes(i) = nexttile(tiles);
           hold(axes(i), 'on');
        end
        title(tiles, '\textbf{State evolution over time}', 'interpreter', 'latex');
        xlabel(tiles, 'Time [min]', 'interpreter', 'latex');
        ylabel(tiles, 'Number of cars [-]', 'interpreter', 'latex');
    end

    for state = 1:4
        for link = 1:2
            if state ~= 3
                lin_idx = min([(state-1), 2])*2 + link;
                plot(axes(lin_idx), t, x(state + 4*(link-1), :), ...
                    'DisplayName', label, 'LineWidth', 1)
                title(axes(lin_idx), titles{lin_idx});
            end
        end
    end
    
%     % Need Matlab R2020b
%     lgd = legend(axes(1));
%     lgd.Orientation = 'horizontal';
%     lgd.Layout.Tile = 'north';
    
end


