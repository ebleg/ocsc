function [out] = alpha_ud_enter(k, par)
% Rate of cars entering link ud
    if k <= 20
        out = (1800 + 10*par.E1)/3600;
    elseif k > 20 && k <= 40
        out = (2100 + 10*par.E2)/3600;
    else
        out = (2300 + 10*par.E3)/3600;
    end
end