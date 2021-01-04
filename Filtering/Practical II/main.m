%% FILTERING & IDENTIFICATION
% Feres Hassan & Emiel Legrand

clear; clc; close all;

student_number = min(4446100, 4362152);

exciteSystem2 = @(u, fs) exciteSystem(student_number, u, fs);

figure

%% Part 1: Data Preprocessing
h = 0.03;
tin = 0:h:3;
u = zeros(size(tin));
u(tin > 0.5) = 1000;

y = exciteSystem2(u, h);

% Find sampling time
% Y = fft(y_filt)

% Spike removal
y_filt = y;
spike_filt = abs(y_filt - mean(y_filt)) > 3*std(y_filt);
y_filt(spike_filt) = nan;
spike_filt = spike_filt & abs(y_filt - nanmean(y_filt)) > 3*nanstd(y_filt);
y_filt(spike_filt) = nan;
spike_filt = abs(y_filt - nanmean(y_filt)) > 3*nanstd(y_filt);
y_filt(spike_filt) = nan;
% Detrending
y_filt = y_filt - nanmean(y_filt);

disp(nanmean(y_filt(tin > 0.5)));
disp(nanmean(y_filt(tin < 0.5)));



plot(tin, y_filt); hold on;
plot(tin, y_filt)
% yyaxis right; set(gca, 'YColor', 'black'); yyaxis
% plot(tin, u);

function moving_median(x, w)
    % x: data vector
    % w: window size
    for i = 

end

