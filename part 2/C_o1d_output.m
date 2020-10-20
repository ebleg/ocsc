function [out] = C_o1d_output(k, par)
% Capacity of the output links for link o1, d
    
    % Two capacities are identical to the ones computed 
    % for link ud
    tmp = C_ud_output(k, par);
    
    % Capacity for right link (ud)
    if k <= 30
        CR = 40 - par.E3;
    else
        CR = 40 + par.E3;
    end
    
    % Change the order (LEFT - STRAIGHT - RIGHT)
    % (o2 o3 ud)
    out = [tmp(2), tmp(3), CR]';
end

