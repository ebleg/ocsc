function [xopt] = lineSearch(fcn, bounds)
    % Line search using Fibonacci method
    N_max = 20;
    fib = fastFibonacci(N_max+1);
    
    % Initialize
    a = min(bounds);
    b = max(bounds);
    L = b - a;
    
    for k = 2:N_max
        L2 = fib(N_max - k + 1)/fib(N_max + 1)*L;
        x1 = a + L2;
        x2 = b - L2;
        f1 = fcn(x1);
        f2 = fcn(x2);

        if f1 > f2
            a = x1;
        else
            b = x2;
        end
    end
    xopt = mean([a b]);
end

function row = fastFibonacci(n)
    row = ones(1, n);
    
    for i = 3:n
        row(i) = row(i-1) + row(i-2);
    end
end