function [x_new] = state_transition2(x, u, par)
    % x = [enter; queue; total_cars];
    % enter = [enter_k enter_k+1 enter_k+2 ...]
    % queue = [queue_1 queue_2 queue_3]
    % total_cars = total number of cars in the link
    % u = [alpha_enter, C_ud, g_udo, C_udo1, C_udo2, C_udo3]
    
    N = numel(x);
    
    % Construct vectors for the three queues
    queue = x((end-3):(end-1));
    
    % Compute N_cars, tau and gamma
    waiting_total = sum(queue);
    tau = floor(((u(2) - waiting_total)*par.l_veh)/par.N_lane/par.v_free/par.c); % Eq 7
    gamma = rem((u(2) - waiting_total)*par.l_veh/par.N_lane/par.v_free, par.c); % Eq 8

    M = N - 4; % Maximum number of delayed states
    
    if M ~= 0 % Safety check; time delays are possible
        % Vectors for the cars entering
        enter = x(1:M);
        enter_new = [u(1); enter(1:end-1)]; % Newly arrived cars at the front
        % i.e. enter_new [arr_k, arr_k-1, arr_k-2, ...], zero delay means
        % the first entry

        % Get the cars arriving to the queue (Eq 6)
        arriving_total = (par.c - gamma)/par.c*enter_new(tau+1) ... % tau + 1 because for zero delay we want the first entry
                         + gamma/par.c*enter(tau+2);        
    end
    
    arriving = arriving_total*par.beta; % Vector of arriving cars per queue
    
    % Find the number of cars leaving per queue
    leaving = min([par.mu.*[u(3) u(3) par.c]'/par.c, ...
                   queue/par.c + arriving, ...
                   u((end-2):end)/par.c], [], 2); % Eq 3
    
    % Update the queue
    queue_new = queue + (arriving - leaving)*par.c;
    
    % Update the total number of cars
    total_cars_new = x(end) + (arriving_total - sum(leaving))*par.c;
    
    % Stack the new state vector
    x_new = [enter_new; queue_new; total_cars_new];
end

