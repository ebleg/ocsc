function [out] = alpha_o1d_enter(k, par)
% Cars entering link o1d (constant in time)
    out = (2000 + 10*par.E1)/3600;
end

