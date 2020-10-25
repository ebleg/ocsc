function [cost, x_rec] = J_disc(u, x0, par)

% map the discrete variables
u = mapvariables(u);

[cost, x_rec] = J(u, x0, par);
end 