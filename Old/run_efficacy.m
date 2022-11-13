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
[SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH,SHv, EHv, DHv, AHv, VHv, UHv,SHc, EHc, DHc, AHc] = age_structured_Malaria_eff(da, na, tfinal_vacc, SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0, 0,0,0,0,0,0,0,0,0,0);
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
rho = sigmoid_prob(Ctot./PH, 'rho'); % prob. of severely infected, EH -> DH
psi = sigmoid_prob(Ctot./PH, 'psi'); % prob. AH -> DH
temp1 = rho.*P.h.*EHv+psi.*lamH.*AHv;
Incidence_vacc = trapz(temp1,1)*P.da;
temp2 = rho.*P.h.*EHc+psi.*lamH.*AHc;
Incidence_control = trapz(temp2,1)*P.da;
eff = (Incidence_control'-Incidence_vacc')./Incidence_control';
figure_setups;
plot(t/365,eff)
xlabel('Year')
ylabel('Efficacy')
title('w/o immunity')
keyboard
%% turn off vaccine
% P.v0 = 0;
Malaria_parameters_transform_vac; 
tfinal_count = 36*30;
t2 = (tfinal_vacc:dt:tfinal_vacc+tfinal_count)'; 
t = [t;t2];
[SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH,SHv, EHv, DHv, AHv, VHv, UHv,SHc, EHc, DHc, AHc] = age_structured_Malaria_eff(da, na, tfinal_count, SH(:,end), EH(:,end), DH(:,end), AH(:,end), VH(:,end), UH(:,end), SM(end), EM(end), IM(end), Cm(:,end), Cac(:,end), Cv(:,end), Ctot(:,end), MH(:,end),SHv(:,end), EHv(:,end), DHv(:,end), AHv(:,end), VHv(:,end), UHv(:,end),SHc(:,end), EHc(:,end), DHc(:,end), AHc(:,end));
% PHv = SHv+EHv+DHv+AHv+VHv+UHv;
% NHv = trapz(PHv,1)*da;
% NHv(end)
% PHc = SHc+EHc+DHc+AHc;
% NHc = trapz(PHc,1)*da;
% NHc(end)
%
PH = SH+EH+DH+AH+VH+UH;
NH = trapz(PH,1)*da;
NM = SM+EM+IM;
[bH,~] = biting_rate(NH,NM);
lamH = FOI_H(bH,IM,NM);
rho = sigmoid_prob(Ctot./PH, 'rho'); % prob. of severely infected, EH -> DH
psi = sigmoid_prob(Ctot./PH, 'psi'); % prob. AH -> DH
temp1 = rho.*P.h.*EHv+psi.*lamH.*AHv;
Incidence_vacc2 = trapz(temp1,1)*P.da;
temp2 = rho.*P.h.*EHc+psi.*lamH.*AHc;
Incidence_control2 = trapz(temp2,1)*P.da;

eff2 = (Incidence_control2'-Incidence_vacc2')./Incidence_control2';
eff = [eff;eff2];
figure_setups;
plot(t/365,eff)
xlabel('Year')
ylabel('Efficacy')

