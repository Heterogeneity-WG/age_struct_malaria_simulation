2%% generate plots for review 1 comment 2 (not sure)
%% vac = constant, fixed number
close all
clear all
% clc
format long
global P

tic

%% numerical config
age_max = 100*365; % max ages in days
P.age_max = age_max;
dt = 5; % time/age step size in days, default = 5;
da = dt;
a = (0:da:age_max)';
na = length(a);

P.dt = dt;
P.a = a;
P.na = na;
P.da = da;

Malaria_parameters_baseline;
Malaria_parameters_baseline_Siaya; % SA based on Nanoro climate profile
Malaria_parameters_transform;
Malaria_parameters_transform_vac;

%% initial condition 'EE' - numerical EE
P.v0 = 0; 
Malaria_parameters_transform_vac;
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');

%% time evolution
P.v0 = 0; Malaria_parameters_transform_vac; % resetting vaccination distribution
tfinal = 365*3;
[t,SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(da,na,0,tfinal,...
    SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
PH = SH+EH+DH+AH+VH+UH;
PH_final = PH(:,end);
NH = trapz(PH,1)*da;
NM = SM+IM+EM;
[bH,bM] = biting_rate(PH,NM); % bH(age, time) matrix
EIR = bH.*IM./NM*365; % EIR(age, time) matrix
EIR_tot = trapz(EIR.*PH,1)*P.da./NH;
EIR_final = EIR_tot(end);
FOI_AH = bM.*trapz(P.betaA*AH)*da./NH; 
FOI_DH = bM.*trapz(P.betaD*DH)*da./NH; 
%%
rho = sigmoid_prob(Ctot(:,end)./PH(:,end), 'rho');
figure_setups; hold on
plot(P.a/365,1./(P.h*rho),'DisplayName','$E\to D$')
xlabel('age')
ylabel('days') % 1/(E->D rate)

psi = sigmoid_prob(Ctot(:,end)./PH(:,end), 'psi');
lamH = bH(:,end).*P.betaM.*IM(end)./NM(end); 
plot(P.a/365,1./(psi.*lamH),'DisplayName','$A\to D$')
legend
xlim([5 10])
