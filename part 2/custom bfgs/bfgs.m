function [xopt, count, xtraj, fval] = bfgs(fcn, x0)
% Initial guess for B
    N = numel(x0);
    B = eye(N);
    
    if size(x0, 2) ~= 1
        x0 = x0';
    end
    
    f0 = fcn(x0);
    g0 = findif(fcn, x0, f0);
    
    % Convergence criterion
    eps = 1e-5;
    count = 1; converged = false;
    maxiter = 200;
    
    xtraj = nan(N, maxiter);

    % Obtain search direction
    while ~converged && (count <= maxiter)
        p = -B*(g0);
        xtraj(:, count) = x0;

        [a, f1, df1] = lineSearch2(@(alph) fcn(x0 + alph*p), 1e100, f0, p'*g0);
        x1 = x0 + a*p;
        
        g1 = findif(fcn, x1, f1, p, df1);
        
        s = x1 - x0;
        y = g1 - g0;
        
        B = B + ((s'*y + y'*B*y)*(s*s'))/(s'*y)^2 - (B*y*s' + s*y'*B)/(s'*y);
        
        if (norm(g1) < eps) || (norm(s) < 1e-10)
            converged = true;
        end
        
        x0 = x1;
        f0 = f1;
        g0 = g1;
        count = count + 1;
        
    end    
    xopt = x1;
    fval = f1;
end