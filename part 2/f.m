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