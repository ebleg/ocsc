clear;

data = readtable('measurements.csv', 'PreserveVariableNames', true);
data_phys = readtable('measurements_physical.csv', 'PreserveVariableNames', true);
delta_t = 3600;

plot(table2array(data(:, 'T_b'))); hold on;
plot(table2array(data_phys(:, 'T_b')));

% q_dot_solar = table2array(data(:, 'q_dot_solar'));
% q_dot_occ = table2array(data(:, 'q_dot_occ'));
% q_dot_ac = table2array(data(:, 'q_dot_ac'));
% q_dot_vent = table2array(data(:, 'q_dot_vent'));
% T_amb = table2array(data(:, 'T_amb'));
% T_b = table2array(data(:, 'T_b'));
% cost = table2array(data(:, 'Phi'))/3600;