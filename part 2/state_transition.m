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
                     *min([1/par.c*diag(par.mu)*[u u par.c]' ...
                           1/par.c*[eye(3) zeros(3,1)]*x + alpha_arrive, ...
                           Co_fcn(k)], [], 2);

end

