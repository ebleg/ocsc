function u = mapvariables (u_disc)

% Discrete time set for green light time
discrete_set= (15:5:45);
u = zeros(60,1);
% Map u from integer values used by GA to the discrete values required

u = discrete_set(u_disc);
end

