%% generate plots for review 1 comment 5
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
% force psi = 0 all the time
% P.psif0 = 0;  P.psif1 = 0; 
Malaria_parameters_baseline_Nanoro; % SA based on Nanoro climate profile
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
NH = trapz(PH,1)*da;
NM = SM+IM+EM;
%% plot year vs. EH AH DH proportions
figure_setups;
EH_time = trapz(EH,1)*P.da;
AH_time = trapz(AH,1)*P.da;
DH_time = trapz(DH,1)*P.da;
infected = EH_time+AH_time+DH_time;
p = area(t/365, [EH_time'./infected',AH_time'./infected',DH_time'./infected']);
colour_mat3 = [0.9290 0.6940 0.1250]; % EH
colour_mat2 = [0.8500 0.3250 0.0980]; % AH
colour_mat7 = [0.6350 0.0780 0.1840]; % DH
p(1).FaceColor = colour_mat3;
p(2).FaceColor = colour_mat2;
p(3).FaceColor = colour_mat7;
xlabel('Time (years)')
ylabel('Prop. of infected pop')
legend({'$E_H$ prop','$A_H$ prop','$D_H$ prop'},'Location','e')
ylim([0 1])
%% compare contribution to DH, two routes
% % EH to DH 
% rho = sigmoid_prob(Ctot./PH, 'rho');
% EHDH= trapz(rho.*P.h.*EH,1)*P.da; % EH -> DH
% psi = sigmoid_prob(Ctot./PH, 'psi');
% NM = SM+EM+IM;
% [bH,~] = biting_rate(PH,NM);
% lamH = FOI_H(bH,IM,NM);
% AHDH = trapz(psi.*lamH.*AH,1)*P.da; % AH -> DH
% figure_setups;
% p=area(t/365,[EHDH'./(EHDH'+AHDH'),AHDH'./(EHDH'+AHDH')]);
% p(1).FaceColor = colour_mat3;
% p(2).FaceColor = colour_mat2;
% xlabel('Time (years)')
% ylabel('Frac. of incidence rates')
% legend({'$E_H \rightarrow D_H$','$A_H \rightarrow D_H$'},'Location','e')
% ylim([0 1])
%% compare contribution to DH, two routes, by age
% % EH to DH 
% rho = sigmoid_prob(Ctot./PH, 'rho');
% temp = rho.*P.h.*EH; % EH -> DH
% EHDH2 = temp(:,end);
% psi = sigmoid_prob(Ctot./PH, 'psi');
% NM = SM+EM+IM;
% [bH,~] = biting_rate(PH,NM);
% lamH = FOI_H(bH,IM,NM);
% temp2 = psi.*lamH.*AH; % AH -> DH
% AHDH2 = temp2(:,end);
% figure_setups;
% p=area(P.a/365,[EHDH2./(EHDH2+AHDH2),AHDH2./(EHDH2+AHDH2)]);
% p(1).FaceColor = colour_mat3;
% p(2).FaceColor = colour_mat2;
% xlabel('Age (years)')
% xlim([0 30])
% ylim([0 1])
% ylabel('Frac. of incidence rates')
% legend({'$E_H \rightarrow D_H$','$A_H \rightarrow D_H$'},'Location','e')
%% plot DH/PH vs age at tfinal
[bH,~] = biting_rate(PH,NM); % bH(age, time) matrix
EIR = bH.*IM./NM*365; % EIR(age, time) matrix
EIR_tot = trapz(EIR.*PH,1)*P.da./NH;
[~,ind_max] = max(EIR_tot);
[~,ind_min] = min(EIR_tot);
Infect = EH+AH+DH;
figure_setups; hold on
plot(P.a/365, DH(:,ind_max)./Infect(:,ind_max),'-','DisplayName','Peak season (max EIR $\approx$ 123)')
plot(P.a/365, DH(:,ind_min)./Infect(:,ind_min),':','DisplayName','Low season (min EIR $\approx$ 21)')
xlim([0 30])
xlabel('Age (years)')
ylabel('Proportion')
title('$D_H/(E_H+A_H+D_H$)')
legend