function [astar, phi_astar, dphi_astar] = lineSearch2(fnc, amax, phi_at0, dphi_at0)
    a0 = 0; 
    
    % Struct with constant values to reduce function parameters of zoom
    meta = struct();
    meta.phi = fnc;
    
    meta.phi_deriv = @(a, phival) findif(meta.phi, a, phival);
    
    % Constants for Wolfe condition (from Nocedal & Wright)
    meta.C1 = 1e-4;
    meta.C2 = 0.9;

    meta.phi_at0 = phi_at0;
    meta.dphi_at0 = dphi_at0;
    
    % Because a0 = 0
    phi0 = meta.phi_at0;
    dphi0 = meta.dphi_at0;      
    
    a1 = 1;
    
    if (a1 < 0)
        error('negative value')
    end
    
    phi1 = meta.phi(a1);
    while phi1 == inf
        a1 = 0.5*a1;
        phi1 = meta.phi(a1);
    end
    dphi1 = meta.phi_deriv(a1, phi1);

    max_iter = 10;
    
    for i = 1:max_iter
        % First Wolfe condition
        if (phi1 > (meta.phi_at0 + meta.C1*a1*meta.dphi_at0)) || ...
                ((phi1 >= phi0) && i > 1)
            [astar, phi_astar, dphi_astar] = zoom(meta, a0, phi0, dphi0, a1, phi1, dphi1); break;
        end
        
        % Second Wolfe condition
        if (abs(dphi1) <= -meta.C2*meta.dphi_at0)
           astar = a1;
           phi_astar = phi1;
           dphi_astar = dphi1;
           break; 
        end
        
        % Third condition
        if (dphi1 >= 0)
           [astar, phi_astar, dphi_astar] = zoom(meta, a1, phi1, dphi1, a0, phi0, dphi0);
           break;
        end
        
        a0 = a1;
        phi0 = phi1;
        dphi0 = dphi1;
        
        a2 = min(a1*2, amax); % From Python implementation
        phi2 = meta.phi(a2);
        while phi2 == inf
            a2 = a1 + (a2 - a1)*0.5;
            phi2 = meta.phi(a2);
        end
        a1 = a2;
        phi1 = phi2;
        dphi1 = meta.phi_deriv(a1, phi1);       
    end
    
    if i == max_iter
       astar = a1;
       phi_astar = phi1;
       dphi_astar = dphi1;
    end
end

function [astar, phi_astar, dphi_astar] = zoom(meta, alow, philow, dphilow, ahigh, phihigh, dphihigh)
    % Progressively reduce interval size if necessary
    max_iter = 10;
    for j = 1:max_iter
        anew = cubicintp(alow, philow, dphilow, ahigh, phihigh, dphihigh);
        phinew = meta.phi(anew);
        dphinew = meta.phi_deriv(anew, phinew);

        if (phinew > (meta.phi_at0 + meta.C1*anew*meta.dphi_at0)) || (phinew >= philow)
            ahigh = anew;
            phihigh = phinew;
            dphihigh = dphinew;
        else
            if abs(dphinew) <= -meta.C2*meta.dphi_at0
               astar = anew;
               phi_astar = phinew;
               dphi_astar = dphinew;
               break;
            end
            if (dphinew*(ahigh - alow)) >= 0
                ahigh = alow;
                phihigh = philow;
                dphihigh = dphilow;
            end
            alow = anew;
            philow = phinew;
            dphilow = dphinew;
        end
        
        if j == max_iter
            astar = ahigh;
            phi_astar = meta.phi(astar);  % Inefficient, but almost never happens
            dphi_astar = meta.phi_deriv(astar, phi_astar);
        end
    end
end

function a = cubicintp(a0, phi0, dphi0, a1, phi1, dphi1)
% Cubic interpolation
    d1 = dphi0 + dphi1 - 3*(phi0 - phi1)/(a0 - a1);
    d2 = sign(a1 - a0)*sqrt(d1^2 - dphi0*dphi1);
    a = a1 - (a1 - a0)*((dphi1 + d2 - d1)/(dphi1 - dphi0 + 2*d2));
end