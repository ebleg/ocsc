function x = mapvariables (x);

% Discrete time set for green light time
Discrete_set= (15:5:45);

% Map ... from integer values used by GA to the discrete values required

x = Discrete_set(x);
end
