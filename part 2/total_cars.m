function [q] = total_cars(x, par)
    % Determine the total number of cars in a link
    N = numel(x);
    m = N - 3;
    
    coeff = [repmat(par.c, 1, m) 1 1 1];
    q = dot(coeff, x);
end

