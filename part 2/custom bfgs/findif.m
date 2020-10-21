function [grad] = findif(fcn, x0, feval, u, D)
    % Use finite differences to compute the gradient of a function
    % feval can be used if there is already a function evaluation available
    % at point x0, as is usually the case. 
    % In some cases, a directional derivative is already known from a line
    % search, which can be used to reduce the number of function
    % evaluations with 1. u represents the direction of the line search and
    % D the derivative. It is assumed that feval will also be present in 
    % that case. 
    
    h = 1e-8;
    
    % Check number of input arguments
    switch nargin
        case 2
            use_feval = false; dir_deriv = false;
            princ_idx = inf;
        case 3
            use_feval = true; dir_deriv = false;
            princ_idx = inf;
        case 5
            use_feval = true; dir_deriv = true;
            [~, princ_idx] = max(u);
        otherwise
            error('Incorrect number of arguments supplied');
    end
    
    grad = nan(size(x0));
    
    for i=1:numel(x0)
        if i ~= princ_idx % Will always happen at princ_idx = inf
            dist = zeros(size(x0));
            dist(i) = h; % Create disturbance vector

            if use_feval
                grad(i) = (fcn(x0 + dist) - feval)./h;
            else
                grad(i) = (fcn(x0 + dist) - fcn(x0))./h;
            end
        end
    end
    
    if dir_deriv
        filt = 1:numel(x0) ~= princ_idx; % All the elements but the one left out
        grad(princ_idx) =  (D - grad(filt)'*u(filt))/u(princ_idx); % Compute the remaining element based on the directional derivative
    end
    
end

