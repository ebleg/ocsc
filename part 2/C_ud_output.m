function [out] = C_ud_output(k, par)
% Capacity of output links for link ud
% CL = lane turning left , CS = straight, CR = turning right
    if k <= 20
        CL = 40 + par.E1;
    elseif k > 20 && k <= 35
        CL = 40 + par.E1 - 2*(k-20);
    elseif k > 35 && k <= 45
        CL = 10 + par.E1;
    else
        CL = 10 + par.E1 + 2*(k-45); 
    end
    
    CS = CL - par.E2;
    
    if k <= 30
        CR = 30 - par.E3;
    else
        CR = 30 + par.E3;
    end
    
    % Make sure the outputs are positive
    out = max([[CL CS CR]', zeros(3, 1)], [], 2);
end

