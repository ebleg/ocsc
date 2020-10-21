%% Parameter file
% Group specific parameters
id_fer = '4362152';
id_emiel = '4446100';

par.E3 = str2double(id_emiel(end)) + str2double(id_fer(end));
par.E2 = str2double(id_emiel(end-1)) + str2double(id_fer(end-1));
par.E1 = str2double(id_emiel(end-2)) + str2double(id_fer(end-2));

clear id_fer id_emiel;

% Parameters for the first link
par.ud = struct();
par.ud.c = 60;
par.ud.N_lane = 3;
par.ud.v_free = 50/3.6;
par.ud.l_veh = 7;
par.ud.l = 1000;
par.ud.beta = [0.33 0.34 0.33]';
par.ud.mu = [1600 1800 1500]'/3600;
par.ud.Cp = par.ud.N_lane*par.ud.l/par.ud.l_veh;

% Parameters for the second link
par.o1d = struct();
par.o1d.c = 60;
par.o1d.N_lane = 3;
par.o1d.v_free = 60/3.6;
par.o1d.l_veh = 7;
par.o1d.l = 1000;
par.o1d.beta = [0.33 0.34 0.33]';
par.o1d.mu = [1600 1800 1500]'/3600;
par.o1d.Cp = par.o1d.N_lane*par.o1d.l/par.o1d.l_veh;

