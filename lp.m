%% LP 
clear; clc; close all;

%% Determine custom parameters
id_fer = '4362152';
id_emiel = '4446100';

E3 = str2double(id_emiel(end)) + str2double(id_fer(end));
E2 = str2double(id_emiel(end-1)) + str2double(id_fer(end-1));
E1 = str2double(id_emiel(end-2)) + str2double(id_fer(end-2));

total_budget = 24000 + 300*E1;

%% First problem
f = -[4 2.5];
A = [3000 1500; 1 1];
b = [total_budget; 12];

lpoptions = optimoptions('linprog', 'Algorithm', 'interior-point', 'Display', 'none');
[x, fval, flag, ~] = linprog(f, A, b, [], [], [0 0], [], lpoptions);
xint = round(x);

[X1, X2] = meshgrid(0:12, 0:12);
Y = zeros(size(X1));

for i = 1:numel(X1)
    Y(i) = f*[X1(i) X2(i)]';
end

% Plot optimization domain
figure; hold on; grid; grid minor;
contourf(X1, X2, Y, 30, 'LineWidth', 0.5)
xrange = 0:10;
plot(xrange, 12-xrange, 'k', 'Linewidth', 1.8)
plot(xrange, total_budget/1500-2*xrange, 'k', 'Linewidth', 1.8)
[Xdots, Ydots] = meshgrid(1:10, 1:10);
scatter(Xdots(:), Ydots(:), 'k', 'filled')
xlabel('X'); ylabel('Y'); title('Analysis of integer solutions')
colorbar;

%% Second question

% Array for the maintenance costs per machine
X1 = [200 200 200 300 300 400 500 600 700 800];
X2 = [1:5 5 5 5 5 5];
Xcost = X1 + X2*E2;

Y1 = [50 50 100 150 150 200 250 300 350 400];
Ycost = Y1 + X2*E3;
clear X1 X2 Y1

% Solve LP problem for every year
maintenance_budget = 4000 + 100*E1; 
remaining = total_budget - [3000 1500]*xint;

xopt_arr = nan(10, 2);
A0 = A; b0 = b;

for yrs = 1:10
    A = [A; sum(Xcost(1:yrs)) sum(Ycost(1:yrs))];
    b = [b; remaining + yrs*maintenance_budget];
%     A = [A0; [Xcost(yrs) Ycost(yrs)]];
%     b = [b0; maintenance_budget];
    
    [xopt_arr(yrs,:) , fval, flag, ~] = linprog(f, A, b, [], [], [0 0], [], lpoptions);
    if flag ~= 1
        error('Infeasible problem')
    end
end

xopt_arr

    
