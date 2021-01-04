function [yhat, xhat] = simsystem(A, B, C, D, x0, u)
    [yhat, ~, xhat] = lsim(ss(A, B, C, D, -1), u, [], x0);
end