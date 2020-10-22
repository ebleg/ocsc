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

