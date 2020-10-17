close all; clear; clc;

id_fer = '4362152';
id_emiel = '4446100';

E3 = str2double(id_emiel(end)) + str2double(id_fer(end));
E2 = str2double(id_emiel(end-1)) + str2double(id_fer(end-1));
E1 = str2double(id_emiel(end-2)) + str2double(id_fer(end-2));

run parameters

C_ud = par.ud.N_lane*par.ud.l/par.ud.l_veh;
g = 22;
alpha_enter = (1800 + 10*E1)/3600;
C_udo1 = 40 + E1; 
C_udo2 = C_udo1 - E2; 
C_udo3 = 30 - E3;

u = [alpha_enter C_ud g C_udo1 C_udo2 C_udo3]';

tau_max = floor(C_ud*par.ud.l_veh/par.ud.N_lane/par.ud.v_free/par.ud.c) + 2;

T = 60;
x = zeros(tau_max+4, T);

for k = 2:T
    x(:, k) = state_transition2(x(:, k-1), u, par.ud);
    bar(x(:, k));
    title(sprintf('k = %d', k))
    pause(0.5);
end