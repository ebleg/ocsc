clear; clc;

run parameters

lb = repmat(15, 60, 1);
ub = repmat(45, 60, 1);
x0 = zeros(8, 1);

best_optimum = inf;
optima = [];
best_sol = zeros(60, 1);
best_init = zeros(60, 1);
sqp_sol.settings = optimoptions(@fmincon, 'Algorithm', 'sqp', ...
                    'MaxFunctionEvaluations', 10000, 'Display', 'off');
for i = 1:2000
    u0 = 15 + (45-15)*rand(60, 1);
    
    [sol, fval, flag] = ...
        fmincon(@(u) J(u, x0, par), u0, [], [], [], [], lb, ub, [], sqp_sol.settings);
    
    if (flag > 0) && (fval < best_optimum)
       best_optimum = fval;
       best_sol = sol;
       best_init = u0;
       optima = [optima best_optimum];
       scatter(1:numel(optima), optima);
       pause(0.1);
       fprintf('New best solution: %f\n', best_optimum)
    end
end
