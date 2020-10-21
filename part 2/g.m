function [y] = g(x, u, par)
% Output
    y = par.ud.c*[0 0 0 1 0 0 0 1]*x; % Both time steps are equal
end