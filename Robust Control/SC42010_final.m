%% Robust and Multivariable Control
%Computer session 1: MIMO Part
%Carmen Petsch en Britt Krabbenborg
clear all
close all
clc

% %% Computer session 1: SISO
% Data=importdata('assignment_data_sc42145.mat');
% A=Data.A;
% B=Data.B;
% C=Data.C;
% D=Data.D;
% % E=Data.E;
% 
% SYS=ss(A,B,C,D);
% s=tf('s');
% 
% %G=C*(s*I-A)^-1*B+D
% transfer=tf(SYS);
% G1=transfer(1,1);
% pole1=pole(G1);
% zero1=zero(G1);
% 
% S1=1/(1+G1);
% figure
% margin(G1)
% figure
% step(G1)
% %% Design a controller with IMC
% Ga=(-s+zero1(3))*(-s+zero1(4))/((s+zero1(3))*(s+zero1(4))); %part with RHP zeros
% Gm=minreal(G1/Ga);  %needs to be minimum phase, minreal cancels pole/zero
% 
% tau_c=1/0.2;
% n=1;
% 
% f=1/(tau_c*s+1)^n;
% K=inv(Gm)*1/(inv(f)-Ga);
% 
% SYS1_OL=G1*K;           %open loop system
% SYS1_CL=G1*K/(1+G1*K);  %closed loop system
% SYS1_OL=minreal(SYS1_OL);
% SYS1_CL=minreal(SYS1_CL);
% 
% pzmap(SYS1_OL)
% 
% %% Plot figures
% S1=1/(1+G1*K);          %sensitivity function
% margin(S1)
% %title('Sensitivity function controlled system') %plot S controlled
% 
% figure
% margin(SYS1_OL)
% %title('Bode plot controlled system')      %plot Bode with controller
% 
% figure
% margin(SYS1_CL)
% %title('Complementary Sensitivity function controlled system')
% 
% figure
% step(SYS1_CL)
% %title('Step response controlled system')
% 
% figure
% pzmap(SYS1_CL)
% 
% 
% %% Control the non-disturbed plant, reference tracking
% t=0:0.1:300;
% u=zeros(size(t));               %reference is zero
% subplot(2,1,1)
% 
% x0=5*ones(1,21);
% x0_u=5*ones(1,15);
% x0_2=ones(1,5);
% 
% SS_sys1=ss(minreal(SYS1_CL));
% SS_sys2=ss(minreal(G3/(1+G3)));
% 
% subplot(4,1,1)
% lsim(SS_sys1,u,t,x0)                  %omega_r goes to zero
% legend('\omega_r sys1[rad/s]')
% 
% subplot(4,1,2)
% lsim(SS_sys2,u,t,x0_2)
% legend('z [m]')
% title('Output Reference tracking')
% 
% subplot(4,1,3)
% SS_ks=ss(minreal(K*S1));
% lsim(SS_ks,u,t,x0_u)                  %input beta
% legend('\beta[rad/s]')
% title('Input Reference tracking')
% 
% subplot(4,1,4)
% plot(t,zeros(size(t)))
% legend('\tau [Nm]')
% 
% %% control under a step disturbance, disturbance rejection
% GV = transfer(1,3);     %transferfunction from disturbance
% S2=1/(1+G2);
% figure
% t=0:0.1:300;
% d=ones(size(t));            
% subplot(4,1,1)
% lsim(S1*GV,d,t)
% legend('\omega_r [rad/s]')
% title('Output disturbance rejection')
% 
% subplot(4,1,2)
% lsim(S2*GV,d,t)
% legend('z [m]')
% title('Output disturbance rejection')
% 
% subplot(4,1,3)
% lsim(-K*GV*S1,d,t)                %input beta
% legend('\beta[rad/s]')
% title('Input disturbance rejection')
% 
% subplot(4,1,4)
% plot(t,zeros(length(t)),'b')
% legend('\tau [Nm]')
% title('Input disturbance rejection')
%% MIMO part Question 1
% omega = 0 
Data=importdata('assignment_data_sc42145.mat');
A=Data.A;
B=Data.B;
C=Data.C;
D=Data.D;
% E=Data.E;

SYS=ss(A,B,C,D);
transfer=tf(SYS);
s=tf('s');
system              = ss(A,B,C,D); 
G_mimo              = tf(system);

G1              = G_mimo(1,1); 
G2              = G_mimo(1,2); 
G3              = G_mimo(2,1); 
G4              = G_mimo(2,2); 
G               = [G1 G2; G3 G4];

lambda          = 1/(1-(tf(G2)*tf(G3)/(tf(G1)*tf(G4)))); 

RGA             = [lambda 1-lambda; 1-lambda lambda]; 
omega           = 0.4*2*pi; 
s1               = j*omega;
RGA1            = evalfr(RGA,s1); 
s2               = 0; 
RGA2            = evalfr(RGA,s2);


%% question 2 
% Gdata           =[Data(1,1) Data(1,2); Data(2,1) Data(2,2)]; 
p_mimo          = pole(G_mimo); 
z_mimo          = tzero(G_mimo);
G_red           = minreal(G_mimo); %no pole/zero cancellation

figure(2)
pzplot(G_red);
p_mimo_red      = pole(G_red);
z_mimo_red      = tzero(G_red);

%% question 3 

wB1             = 0.4*2*pi; % desired bandwith 
At              = 1E-4; % desired attenuation inside bandwith
M               = 1.8; % desired bound on hinfnorm(S)
x=10;
n=1;
Wp              = [(s/M^(1/n)+wB1)^n/(s+wB1*At^(1/n))^n 0; 0 0.2];
Wu              = [0.01 0;0 (5E-3*s^2+7E-4*s+5E-5)/(s^2+14E-4*s+10E-6)];

Wt              = []; 

Wp11            = Wp(1,1); 

%% question 5 

P5 = [Wp -Wp*G; zeros(2) Wu; eye(2) -G];  % generalized plant
 
%% question 8 mixed sensitivity 
%mixed sensitivity generalized controller
P_output        = 2; 
P_input         = 2; 

close all
systemnames = 'G Wp Wu'; % Define systems
inputvar ='[w(2); u(2)]'; % Input generalized plant
input_to_G = '[u]';
input_to_Wu = '[u]';
input_to_Wp = '[w-G]';
outputvar = '[Wp; Wu; -G+w]'; % Output generalized plant
sysoutname ='P8';
sysic;

[K8,CL8,GAM8,INFO8] = hinfsyn(P8,P_output,P_input); %hinf design

CL8_min             = minreal(CL8);
K8_min              = minreal(K8); 
P8_min              = minreal(P8); 
P8_ss               = minreal(ss(P8_min)); 

GK8                 = minreal(G*K8);
KG8                 =(minreal(K8)*minreal(G));
L8                  = minreal(eye(size(GK8))+GK8);
detL                = L8(1,1)*L8(2,2)-L8(1,2)*L8(2,1);
figure(1)
nyquist(L8)
figure(2)
nyquist(detL)
figure(3)
nyquist(detL)
axis([-50 50 -50 50])
      
S9=minreal(inv(eye(2)+GK8)); %sensitivity function
GV2=transfer(:,3);
T9=minreal(inv(eye(2)+GK8))*GK8;    %complementary sensitivity function

%% question 9 Time simulations
%Control the non-disturbed plant, reference tracking
figure
t=0:0.1:200;
u=[zeros(size(t));
   zeros(size(t))];          %reference is zero

T9_ss=ss(minreal(T9));
x0=5*ones(1,8);
subplot(4,1,1)
[yo,t1,x]=lsim(T9_ss,u,t,x0);
plot(t,yo(:,1))
title('Reference tracking') %omega_r should go to zero
legend('\omega_r[rad/s]')
title('Output reference tracking')

subplot(4,1,2)
plot(t,yo(:,2))             %z should go to zero               
legend('z[m]')
title('Output reference tracking')
 
T10_ss=ss(minreal(inv(eye(2)+KG8)*K8));
x0_u=ones(1,9);
[yi,t1,x]=lsim(T10_ss,u,t,x0_u);

subplot(4,1,3)
plot(t,yi(:,1))
legend('\beta[rad]')
title('Input reference tracking')

subplot(4,1,4)
plot(t,yi(:,2))
legend('\tau_e[rad]')
title('Input reference tracking')

%% Disturbance rejection
% output disturbance rejection
figure;
t=0:0.1:1000;
u = ones(size(t));
subplot(4,1,1)
[yo,~,x]=lsim(S9*GV2,u,t);
plot(t,yo(:,1))
title('Output disturbance rejection')
legend('\omega_r[rad/s]')

subplot(4,1,2)
plot(t,yo(:,2))
title('Output disturbance rejection')
legend('z[m]')

%  Input disturbance rejection
[yi,t1,x]=lsim(minreal(inv(eye(2)+KG8)*-K8*GV2),u,t);
subplot(4,1,3)
plot(t,yi(:,1))
title('Input disturbance rejection')
legend('\beta[rad]')
subplot(4,1,4)
plot(t,yi(:,2))
title('Input disturbance rejection')
legend('\tau_e[Nm]')

%% Question 10 improved weights
% plot S, 1/Wp, upperbound, omega en beta 
figure
subplot(2,2,1)
bode(S9(1,1))
legend('Sensitivity (1,1)')
hold on
bode(1/Wp11)
legend('Upper bound Wp(1,1)') 
hold on
bode(Wp11*S9(1,1))
legend('Sensitivity (1,1)','Upper bound(1,1)','Wp*S')
title('Upper bound on sensitivity function')

subplot(2,2,2)
bode(S9(1,2))
legend('Sensitivity (1,2)')
title('No upper bound on sensitivity function')

subplot(2,2,3)
bode(S9(2,1))
legend('Sensitivity (2,1)')
title('No upper bound on sensitivity function')

subplot(2,2,4)
bode(S9(2,2))
legend('Sensitivity (2,2)')
hold on
bode(1/Wp(2,2))
hold on
bode(Wp(2,2)*S9(2,2))
legend('Sensitivity (2,2)','Upper bound(2,2)','Wp*S')
title('Upper bound on sensitivity function')
%%  plot K*S, 1/Wu, upperbound 
figure
bound=minreal(K8*S9);
subplot(2,2,1)
bode(bound(1,1))
legend('K*S')
hold on
bode(1/Wu(1,1))
hold on
bode(Wu(1,1)*bound(1,1))
legend('K*S','Upper bound','Wu*K*S')
title('Upper bound on K*S')

subplot(2,2,2)
bode(bound(1,2))
legend('K*S')
title('No upper bound on K*S')

subplot(2,2,3)
bode(bound(2,1))
legend('K*S')
title('No upper bound on K*S')

subplot(2,2,4)
bode(bound(2,2))
legend('K*S')
hold on
bode(1/Wu(2,2))
hold on
bode(Wu(2,2)*bound(2,2))
legend('K*S','Upper bound','Wu*K*S')
title('Upper bound on K*S')

%% infinity norm system
N=[minreal(Wp*S9);minreal(Wu*K8*S9)];
Hinf=hinfnorm(N)
HinfWp=hinfnorm(minreal(Wp*S9))

%% Comparison SISO and MIMO
% sensitivity functions both systems
G2=transfer(2,2);
figure
bode(S1)
hold on
SYS2_OL=minreal(G1*K8(1,1)+G2*K8(2,1));
Gcomp=[G1 G2;G3 G4];
S2=inv(eye(2)+Gcomp*K8);
bode(S2(1,1))
hold on
bode(1/Wp11)
title('Sensitivy function system1 and system 2')
legend('SISO controller','MIMO controller','Upper bound S2')

%% Computer session 2
% importing variables 
Data                = importdata('DataFWT.mat');

A                   = Data.A;
B                   = Data.B;
C                   = Data.C; 
D                   = Data.D; 

s                   = tf('s');

system              = ss(A,B,C,D); 
G_mimo              = tf(system); 

G1                  = G_mimo(1,1); 
G2                  = G_mimo(1,2); 
G3                  = G_mimo(2,1); 
G4                  = G_mimo(2,2); 
G                   = [G1 G2; G3 G4];


wB                  = 0.4*2*pi;         % desired bandwith 
At                  = 1E-4;             % desired attenuation inside bandwith
M                   = 1.8;              % desired bound on hinfnorm(S)
Wp                  = [(s/M+wB)/(s+wB*At) 0; 0 0.2];
Wu                  = [0.01 0;0 (5E-3*s^2+7E-4*s+5E-5)/(s^2+14E-4*s+10E-6)];

%% Adjusted weights for DK iterations
wB                  = 0.4;              % desired bandwith 
At                  = 1E-4;             % desired attenuation inside bandwith
M                   = 1.8;              % desired bound on hinfnorm(S)
Wp                  = [(s/M+wB)/(s+wB*At) 0; 0 0.2];
Wu                  = [0.01 0;0 (5E-3*s^2+7E-4*s+5E-5)/(s^2+14E-4*s+10E-6)];

%% Generalized plant

close all
systemnames = 'G Wp Wu';                % Define systems
inputvar ='[w(2); u(2)]';               % Input generalized plant
input_to_G = '[u]';
input_to_Wu = '[u]';
input_to_Wp = '[w-G]';
outputvar = '[Wp; Wu; -G+w]';           % Output generalized plant
sysoutname ='P';
sysic;

[K,CL,GAM2,INFO2]   = hinfsyn(P,2,2);   % Hinf design
K                   = minreal(K);
N1                  = lft(P,K);

G                   = minreal(G);
L                   = minreal(G*K);
IL                  = minreal(eye(size(L))+L);
detIL               = IL(1,1)*IL(2,2) - IL(1,2)*IL(2,1);

figure(5)
nyquist(detIL) 
figure(6)
nyquist(detIL)
axis([-50 50 -50 50])

%% question 3 

W_i1        = tf([1/(16*pi)  0.3],[1/(64*pi) 1]); 
W_i2        = tf([1/(16*pi)  0.3],[1/(64*pi) 1]); 
W_o1        = tf([0.05  0.2],[0.01  1]); 
W_o2        = tf([0.05  0.2],[0.01  1]); 

Wi          = append(W_i1,W_i2); 
Wo          = blkdiag(W_o1,W_o2); 

Delta_o     = [ultidyn('D_o1',[1 1],'Bound',1), 0 ; 0,  ultidyn('D_o2',[1 1],'Bound',1)];
Delta_i     = [ultidyn('D_i1',[1 1],'Bound',1), 0 ; 0,  ultidyn('D_i2',[1 1],'Bound',1)];

Delta       = append(Delta_i, Delta_o);
Gp          = (eye(2) + Wo*Delta_o*eye(size(G)))*G*(eye(2) + Wi*Delta_i*eye(size(G)));

figure(1)
sigma(Gp)
grid on 

figure(2)
bodemag(Wi,Wo)
legend('W_{i}','W_{o}','Interpreter','tex')
grid on 

H           = ultidyn('H',[2,2],'Bound',1); 
HWi         = Wi*H; 
HWo         = Wo*H; 

figure(3) 
bodemag(HWi,Wi,'r'); 
legend('Uncertain Wi','Wi','Interpreter','tex')

figure(4) 
bodemag(HWo,Wo,'r'); 
legend('Uncertain Wo','Wo','Interpreter','tex')

%% question 4
systemnames = 'G Wp Wu Wi Wo' ;                     % Define systems
inputvar ='[udeltai(2); udeltao(2); w(2); u(2)]';   % Input generalized plant
input_to_G='[u+udeltai]';
input_to_Wu='[u]';
input_to_Wp='[w-G-udeltao]';
input_to_Wi='[u]';
input_to_Wo='[G]';
outputvar = '[Wi; Wo; Wp; Wu; w-G-udeltao]';        %output generalized plant
sysoutname = 'P4'; 
cleanupsysic= 'yes'; 
sysic;

P11                 = P4(1:6,1:6);
P12                 = P4(1:6,7:8);
P21                 = P4(8:10,1:6);
P22                 = P4(7:end,7:end);

N                   = lft(P4,K);
N11                 = N(1:4,1:4);
N12                 = N(1:4,5:6);
N21                 = N(5:8,1:4);
N22                 = N(5:8,5:6);

%% question 5

F                   = minreal(lft(Delta,N));
omega               = logspace(-3,3,100);
Nf                  = frd(N,omega);  


% I+GK is denominator in the entire N matrix, so check gen. Nyquist of
% det(I+GK) must not encircle origin

G                   = minreal(G);
L                   = minreal(G*K);
IL                  = minreal(eye(size(L))+L);
detIL               = IL(1,1)*IL(2,2) - IL(1,2)*IL(2,1);

figure(5)
nyquist(detIL) 
figure(6)
nyquist(detIL)
axis([-50 50 -50 50])


% Nominal stability check if N is internally stable 
% so all eigenvalues are negative  
if max(real(eig(N))) < 0
    disp('We have Nominal Stability')
else
    disp('We have no Nominal Stability')
end


% Nominal performance check if mu of N22 < 1
blk                 = [2 4];
[mubnds,muinfo]     = mussv(Nf(5:8,5:6),blk,'c');
muNP                = mubnds(:,1); 
[muNPinf, muNPw]    = norm(muNP,inf);

if muNPinf < 1 && max(real(eig(N))) < 0
    disp('We have Nominal Performance')
else
    disp('We have no Nominal Performance')
end

% Robust stability check if mu of N11 < 1

blk                 = [1 1; 1 1; 1 1; 1 1];
[mubnds,muinfo]     = mussv(Nf(1:4,1:4), blk, 'c');
muRS                = mubnds(:,1); 
[muRSinf, muRSw]    = norm(muRS,inf);

if muRSinf < 1 && max(real(eig(N))) < 0
    disp('We have Robust Stability')
else
    disp('We have no Robust Stability')
end

% Robust performance check if 

blk                 = [1 1; 1 1; 1 1; 1 1; 2 4];
[mubnds,muinfo]     = mussv(Nf(1:8,1:6), blk, 'c');
muRP                = mubnds(:,1); 
[muRPinf, muRPw]    = norm(muRP,inf);

if muRPinf < 1 && max(real(eig(N))) < 0
    disp('We have Robust Performance')
else
    disp('We have no Robust Performance')
end

%% question 6

hinfnorm(N1)
muNPinf

%% Computer session 3 
 
%% DK iterations automatic 

%check if robust performance with 4 inputs
%We now have robust performance with value of 0.9818
%so the generalized plant is not wrong

P4                      = [zeros(2) zeros(2) zeros(2) Wi;  Wo*G zeros(2) zeros(2) Wo*G; -Wp*G -Wp Wp -Wp*G;   zeros(2) zeros(2) zeros(2) Wu;    -G -eye(2) eye(2) -G]; 

omega                   = logspace(-3,3,130);
Punc                    = lft(Delta,P4);
opt                     = dkitopt('FrequencyVector', omega);
[K_DK,clp,bnd,dkinfo]    = dksyn(Punc,2,2,opt);  

N_DK                    = minreal(lft(P4,K_DK));
Nf_DK                   = frd(N_DK,omega);

blk                 = [1 1; 1 1; 1 1; 1 1; 2 4];
[mubnds,muinfo]     = mussv(Nf_DK, blk, 'c');
muRP                = mubnds(:,1); 
[muRPinf, muRPw]    = norm(muRP,inf);

if muRPinf < 1 && max(real(eig(N_DK))) < 0
    disp('We have Robust Performance')
else
    disp('We have no Robust Performance')
end

muRPinf_ok         = muRPinf

%% DK iterations manual

P4                  = [zeros(2) zeros(2) zeros(2) Wi;  Wo*G zeros(2) zeros(2) Wo*G; -Wp*G -Wp Wp -Wp*G;   zeros(2) zeros(2) zeros(2) Wu;    -G -eye(2) eye(2) -G]; 


blk                 = [1 1; 1 1; 1 1;1 1; 2 4];
omega               = logspace(-3,3,121); 
K3 = K; 

mu_path=[];


for i = 1:1:10
    Nf                  = frd(lft(P4,K3,2,2),omega); 
    [mubnds,muinfo]     = mussv(Nf,blk,'c');
    muRP                = mubnds(:,1);
    [muRPinf, muRPw]    = norm(muRP, inf); 
    [VDelta,VSigma,VLmi]= mussvextract(muinfo); 
    D                   = VSigma.DLeft; 
    dd1                 = fitmagfrd((D(1,1)/D(5,5)),6);
    dd2                 = fitmagfrd((D(2,2)/D(5,5)),6);
    dd3                 = fitmagfrd((D(3,3)/D(5,5)),6);     
    dd4                 = fitmagfrd((D(4,4)/D(5,5)),6);
    DscaleLeft          = minreal(append(dd1, dd2, dd3, dd4, tf(eye(6))));
    DscaleRight         = minreal(append(1/(dd1), 1/(dd2), 1/(dd3), 1/(dd4),tf(eye(4))));
    [K3,CL3,GAM3,INFO3] = hinfsyn(minreal(DscaleLeft*P4*DscaleRight),2,2);
   
    figure(7)
    bodemag(mubnds(1,1),omega)
    hold on
    mu_path(i)          = muRPinf;
    muRPinf
    if muRPinf < 1
        break 
    end 
end 

title('Structured singular value for each iteration n=1:10')
legend('1','2','3','4','5','6','7','8','9','10') 

N_DK = minreal(lft(P4,K3));
Nf_DK = frd(N_DK,omega);

% Nominal stability check if N is internally stable 
% so all eigenvalues are negative  
if max(real(eig(N))) < 0
    disp('We have Nominal Stability')
else
    disp('We have no Nominal Stability')
end


% Nominal performance check if mu of N22 < 1

blk                 = [2 4];
[mubnds,muinfo]     = mussv(Nf_DK(5:8,5:6),blk,'c');
muNP                = mubnds(:,1); 
[muNPinf, muNPw]    = norm(muNP,inf);

if muNPinf < 1 && max(real(eig(N_DK))) < 0
    disp('We have Nominal Performance')
else
    disp('We have no Nominal Performance')
end

% Robust stability check if mu of N11 < 1

blk                 = [1 1; 1 1; 1 1; 1 1];
[mubnds,muinfo]     = mussv(Nf_DK(1:4,1:4), blk, 'c');
muRS                = mubnds(:,1); 
[muRSinf, muRSw]    = norm(muRS,inf);

if muRSinf < 1 && max(real(eig(N_DK))) < 0
    disp('We have Robust Stability')
else
    disp('We have no Robust Stability')
end

% Robust performance check if 

if muRPinf < 1 && max(real(eig(N_DK))) < 0
    disp('We have Robust Performance')
else
    disp('We have no Robust Performance')
end


% generalized nyquist of new controller

L3                  = minreal(G*K3);
IL3                 = minreal(eye(size(L3))+L3);
detIL3              = IL3(1,1)*IL3(2,2) - IL3(1,2)*IL3(2,1);

figure(8)
nyquist(detIL3) 
figure(9)
nyquist(detIL3)
axis([-50 50 -50 50])

%% simulations and freq plots with robust controller

GK3             = minreal(G*K3); 
S3              = minreal(inv(eye(2)+GK3));

KdKG            = minreal(K3)*minreal(G);
KGinf            = minreal(K)*minreal(G); 

GK               =G*K;
S               = minreal(inv(eye(2)+GK));
w               = logspace(-4,4,1000)';
GV2              =G_mimo(:,3);

figure(10)
bode(-K*G,w);
hold on
bode(-K3*G,w);
grid on
legend('H\infty','DK-iteration')
title('MIMO Bode comparison')


figure(11)
sigmaplot(inv(Wp(1,1)),'g');
hold on
sigmaplot(S(1,1));
sigmaplot(S3(1,1));
grid on
legend('1/|Wp|','H\infty','DK-iteration');
title('MIMO Sensitivity comparison')

% Disturbance rejection
% output disturbance rejection
figure(12)
t=0:0.1:1000;
u = ones(size(t));
subplot(4,1,1)
[yo,t1,x]=lsim(S3*GV2,u,t);
[yo1,t1,x]=lsim(S*GV2,u,t);
plot(t,yo(:,1))
hold on
plot(t,yo1(:,1))
title('Output disturbance rejection')
legend('\omega_r[rad/s] DK iteration','\omega_r[rad/s] H_{inf} controller')

subplot(4,1,2)
plot(t,yo(:,2))
hold on 
plot(t,yo1(:,2))
title('Output disturbance rejection')
legend('z[m] DK iteration','z[m] H_{inf} controller')

%  Input disturbance rejection
[yi,t1,x]=lsim(minreal(inv(eye(2)+KdKG)*-K3*GV2),u,t);
[yi1,t1,x]=lsim(minreal(inv(eye(2)+KGinf)*-K*GV2),u,t);
subplot(4,1,3)
plot(t,yi(:,1))
hold on 
plot(t,yi1(:,1))
title('Input disturbance rejection')
legend('\beta[rad] DK iteration','\beta[rad] H_{inf} controller')
subplot(4,1,4)
plot(t,yi(:,2))
hold on 
plot(t,yi1(:,2))
title('Input disturbance rejection')
legend('\tau_e[Nm] DK iteration','\tau_e[Nm] H_{inf} controller')
