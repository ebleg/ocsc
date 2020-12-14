function [y] = g(x, u, par)
% Output % number of vehicles on link (u,d) and link (o1,d)
    y = par.ud.c*[0 0 0 1 0 0 0 1]*x; % Both time steps are equal
end