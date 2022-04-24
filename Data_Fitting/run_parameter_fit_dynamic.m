clear all;
close all;
clc;
format long;
global P lP
global F

% numerical config
tfinal = 50*365; age_max = 100*365; P.age_max = age_max;
dt = 20; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t; P.tfinal = tfinal;

%% load data - interpolated function
s =  load('Filipe_paper/F_Filipe.mat','F');
F = s.F;

%% optimization - fit to data 
Malaria_parameters_baseline;
options = optimset('Display','iter','TolX',10^-5,'MaxIter',20);
% [phi_s phi_r rho_s=psi_s rho_r=psi_r]
% x0 = [2.432473210664639   1.277554702429119   3.186642715992263   1.030298116795388]; % SIAP paper fitting
% x0 = [3.741225742715598   0.100038517193148   5.609668653723179   2.102940450587619]; % Tfinal = 10 years fitting
x0 = [3.924837426334165   0.106405511374953   6.153593975708337   2.515512993061040]; % Tfinal = 25 years fitting
lb = [0, 0.1, 0, 0.1];
ub = [8, 10, 8, 10];

[x,fval] = fmincon(@fun_Filipe_dynamic,x0,[],[], [], [], lb, ub, [], options);
% x =    0.064773667811742   3.506525885314898   0.145462635652681   9.293057095109972;
keyboard
% ----> update Malaria_parameters_baseline.m with the fitted results <-----

%% plot sigmoids with the populational average (in legend)
tfinal = 50*365; age_max = 100*365; P.age_max = age_max;
dt = 20; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t; P.tfinal = tfinal;

x = [3.792409077935755   0.100082759162754   7.001613781444139   2.799435902292193];
Malaria_parameters_baseline;
P.phi_s_2 = x(1);
P.phi_r_2 = x(2); 
P.rho_s_2 = x(3);
P.rho_r_2 = x(4); 
P.psi_s_2 = x(3);
P.psi_r_2 = x(4);
Malaria_parameters_transform;
[SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0] = age_structured_Malaria_IC('init');
[SH, EH, DH, AH, SM, EM, IM, ~, ~, Ctot] = age_structured_Malaria(P.da,P.na,P.tfinal,SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0);
PH = SH(:,end)+EH(:,end)+DH(:,end)+AH(:,end);
NH = trapz(PH)*P.da;
Ctot_pp = Ctot(:,end)./PH;
figure_setups; hold on
cc = linspace(0,max(Ctot_pp),100);
phi_curve = sigmoid_prob(cc, 'phi');
rho_curve = sigmoid_prob(cc, 'rho');
psi_curve = sigmoid_prob(cc, 'psi');
plot(cc,phi_curve,'-')
plot(cc,rho_curve,'--')
% plot(cc,psi_curve,'-.')
% population average sigmoids
phi_ave = trapz(sigmoid_prob(Ctot_pp, 'phi').*PH/NH)*P.da;
rho_ave = trapz(sigmoid_prob(Ctot_pp, 'rho').*PH/NH)*P.da;
psi_ave = trapz(sigmoid_prob(Ctot_pp, 'psi').*PH/NH)*P.da;
legend('$\phi(\tilde{C}_{H})$','$\rho(\tilde{C}_{H})=\psi(\tilde{C}_{H})$','Location','e')
axis([0 max(cc) 0 1])
xlabel('$\tilde{C}_{H}$')
ylabel('Probability')
%title('Calibrated linking functions')
[phi_ave rho_ave psi_ave]
keyboard

%% plotting heatmap (age, EIR, immunity level)
tfinal = 50*365; age_max = 100*365; P.age_max = age_max;
dt = 20; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t;

Malaria_parameters_baseline;
P.phi_s_2 = x(1);
P.phi_r_2 = x(2); 
P.rho_s_2 = x(3);
P.rho_r_2 = x(4); 
P.psi_s_2 = x(3);
P.psi_r_2 = x(4);
Malaria_parameters_transform;
immunity_feedback = 1;
if immunity_feedback == 0 
    % population average sigmoids: f0 = f1 = average
    P.phi_f_0 = 0.915792480087329; % value at zero
    P.phi_f_1 = 0.915792480087329; % value at L (function saturates to this value)
    
    P.rho_f_0 = 0.114825053290306; % value at zero
    P.rho_f_1 = 0.114825053290306; % value at L (function saturates to this value)  
    
    P.psi_f_0 = 0.114825053290306; % value at zero
    P.psi_f_1 = 0.114825053290306; % value at L (function saturates to this value)   
    var_list = [0.01:0.05:1.0].^2; 
else
    var_list = [0.01:0.05:5].^2;
end
x = P.a/365;
y = zeros(1,length(var_list));
z = zeros(na,length(var_list));
for jj = 1:length(var_list)
    P.betaM = var_list(jj);
    Malaria_parameters_transform;
    [SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0] = age_structured_Malaria_IC('init');
    [SH, EH, DH, AH, SM, EM, IM, ~, ~, Ctot] = age_structured_Malaria(P.da,P.na,P.tfinal,SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0);
    EIR = fit_EIR(SH,EH,DH,AH,SM, EM, IM);
    PH = SH+EH+DH+AH;
    y(1,jj) = EIR(end); % aEIR
    z(:,jj) = Ctot(:,end)./PH(:,end); % final Ctot at EE    
end
% figure_setups;
% plot(var_list,y,'-o');
% keyboard
figure_setups;
imagesc(x,y,z')
xlim([0 20])
xlabel('Age (years)')
ylabel('aEIR')
% title(['Immunity levels, feedback = ',num2str(immunity_feedback)]);
title('Immunity level per person');
set(gca,'YDir','normal');
colormap jet
colorbar('Ticks',0:2:15);

keyboard

function EIR = fit_EIR(SH,EH,DH,AH,SM,EM,IM)
global P
NH = trapz(SH+EH+DH+AH)*P.da;
NM = SM+EM+IM;
[bH,~] = biting_rate(NH,NM);
IM_frac = IM./NM;
EIR = bH.*IM_frac*365; % annual EIR
end