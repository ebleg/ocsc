function [xopt, count, xtraj] = bfgs_gradproj(fcn, x0, lb, ub)
% Initial guess for B
    N = numel(x0);
    B = eye(N);
    
    if size(x0, 2) ~= 1
        x0 = x0';
    end
    
    f0 = fcn(x0);
    g0 = findif(fcn, x0, f0);
    
    % Convergence criterion
    eps = 1e-8;
    count = 1; converged = false;
    maxiter = 200;
    
    xtraj = nan(N, maxiter);

    % Obtain search direction
    while ~converged && (count <= maxiter)
        p = -B*(g0);
        xtraj(:, count) = x0;
        
        % Check for active constraints
        check_lb = (x0 - lb) < 1e-8;
        check_ub = (ub - x0) < 1e-8;
        
        % Active constraint
        if any(check_lb) || any(check_ub)
            active_idx = [find(check_lb) find(check_ub)];
            p(active_idx) = 0; % Project on box
        end
        
        % Find value for the upper bound of alpha
        alpha_max = inf;
        for i = 1:numel(x0)
           if p(i) > 0
               tmp = (ub(i) - x0(i))/p(i);
           elseif p(i) < 0
               tmp = (lb(i) - x0(i))/p(i);
           else
               tmp = inf;
           end
           
           if tmp < alpha_max
               alpha_max = tmp;
           end
        end

        a = lineSearch(@(alph) fcn(x0 + alph*p),[0, alpha_max]);
        x1 = x0 + a*p;
        
        % check violations
        if any(x1 < lb) || any(x1 > ub) || any(isnan(x1))
            error('constraint violation');
        end
        
        % Gradient at new point (f1 already calculated in line search)
        g1 = findif(fcn, x1);
        
        % Vectors for hessian update
        s = x1 - x0;
        y = g1 - g0;
        
        B = B + ((s'*y + y'*B*y)*(s*s'))/(s'*y)^2 - (B*y*s' + s*y'*B)/(s'*y);
        
        if (norm(g1) < eps) || (norm(s) < 1e-10)
            converged = true;
        end
        
        x0 = x1;
        g0 = g1;
        count = count + 1;
    end    
    xopt = x1;
end