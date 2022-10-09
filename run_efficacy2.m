clearvars
% close all
clc

format long
global P

%% numerical config
age_max = 100*365; % max ages in days
P.age_max = age_max;
dt = 10; % time/age step size in days, default = 5;
da = dt;
a = (0:da:age_max)';
na = length(a);

P.dt = dt; 
P.a = a;
P.na = na;
P.da = da;

% model parameters
Malaria_parameters_baseline;
Malaria_parameters_transform; 
Malaria_parameters_transform_vac;

%% initial condition 'EE' - numerical EE
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE');
NN_S = trapz(SH0)*P.da;

P.v0 = 15; tfinal_vacc = 5*365; 
t = (0:dt:tfinal_vacc)'; 
Malaria_parameters_transform_vac; 
[SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH,SHv, EHv, DHv, AHv, VHv, UHv,SHc, EHc, DHc, AHc, Ctotv, Ctotc] = ...
    age_structured_Malaria_eff2(da, na, tfinal_vacc, ...
    SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0, 0,0,0,0,0,0,0,0,0,0);
% vacc = trapz(trapz(P.v.*SH,1)*P.da)*P.dt; %
% PHv = SHv+EHv+DHv+AHv+VHv+UHv;
% NHv = trapz(PHv,1)*da;
% NHv(end)
% PHc = SHc+EHc+DHc+AHc;
% NHc = trapz(PHc,1)*da;
% NHc(end)
PH = SH+EH+DH+AH+VH+UH;
NH = trapz(PH,1)*da;
NM = SM+EM+IM;
[bH,~] = biting_rate(NH,NM);
lamH = FOI_H(bH,IM,NM);
PHv = SHv+EHv+DHv+AHv+VHv+UHv;
rhov = sigmoid_prob(Ctotv./PHv, 'rho'); % prob. of severely infected, EH -> DH
psiv = sigmoid_prob(Ctotv./PHv, 'psi'); % prob. AH -> DH
rhov(PHv==0)=0; psiv(PHv==0)=0;
temp1 = rhov.*P.h.*EHv+psiv.*lamH.*AHv;
Incidence_vacc = trapz(temp1,1)*P.da;
PHc = SHc+EHc+DHc+AHc;
rhoc = sigmoid_prob(Ctotc./PHc, 'rho'); % prob. of severely infected, EH -> DH
psic = sigmoid_prob(Ctotc./PHc, 'psi'); % prob. AH -> DH
rhoc(PHc==0)=0; psic(PHc==0)=0;
temp2 = rhoc.*P.h.*EHc+psic.*lamH.*AHc;
Incidence_control = trapz(temp2,1)*P.da;
eff = (Incidence_control'-Incidence_vacc')./Incidence_control';
figure_setups;
plot(t/365,eff)
xlabel('Year')
ylabel('Efficacy')
title('with immunity')


