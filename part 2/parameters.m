%% Parameter file
% Parameters for the first link
% sadfasfd
par_ud = struct();
par_ud.c = 60;
par_ud.N_lane = 3;
par_ud.v_free = 50/3.6;
par_ud.l_veh = 7;
par_ud.l = 1000;
par_ud.beta = [0.33 0.34 0.33];
par_ud.mu = [1600 1800 1500]/3600;

% Parameters for the second link
par_o1d = struct();
par_01d.c = 60;
par_01d.N_lane = 3;
par_01d.v_free = 60/3.6;
par_01d.l_veh = 7;
par_01d.l = 1000;
par_01d.beta = [0.33 0.34 0.33];
par_o1d.mu = [1600 1800 1500]/3600;
