%% Parameter file
% Parameters for the first link
% sadfasfd
par.ud = struct();
par.ud.c = 60;
par.ud.N_lane = 3;
par.ud.v_free = 50/3.6;
par.ud.l_veh = 7;
par.ud.l = 1000;
par.ud.beta = [0.33 0.34 0.33]';
par.ud.mu = [1600 1800 1500]'/3600;

% Parameters for the second link
par.o1d = struct();
par.o1d.c = 60;
par.o1d.N_lane = 3;
par.o1d.v_free = 60/3.6;
par.o1d.l_veh = 7;
par.o1d.l = 1000;
par.o1d.beta = [0.33 0.34 0.33]';
par.o1d.mu = [1600 1800 1500]'/3600;
