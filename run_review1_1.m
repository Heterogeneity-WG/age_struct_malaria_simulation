%% generate plots for review 1 comment 1
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
% Malaria_parameters_baseline_Siaya; % SA based on Nanoro climate profile
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
P.zeta = P.zeta_fun(a);
FOI_AH = bM.*trapz(P.betaA*AH)*da./NH; 
FOI_DH = bM.*trapz(P.betaD*DH)*da./NH; 
FOI_AHw = bM.*trapz(P.betaA*P.zeta.*AH)*da./NH; 
FOI_DHw = bM.*trapz(P.betaD*P.zeta.*DH)*da./NH; 

%% plot1 age vs. AH DH prevalences
figure_setups;
colour_mat2 = [0.8500 0.3250 0.0980];
colour_mat3 = [0.9290 0.6940 0.1250];
p = area(a/365,[AH(:,end)./(DH(:,end)+AH(:,end)),DH(:,end)./(DH(:,end)+AH(:,end))]);
p(1).FaceColor = colour_mat3;
p(2).FaceColor = colour_mat2;
legend({'$A_H/(A_H+D_H)$','$D_H/(A_H+D_H)$'},'Location','se')
title('Prevalence');
xlabel('Age (years)');
ylabel('Proportion')
grid on
axis([0 P.age_max/365 0 1]);
xlim([0 30])
%% plot2 time vs. AH DH FOIs
figure_setups;
colour_mat3 = [0.9290 0.6940 0.1250]; % EH
colour_mat2 = [0.8500 0.3250 0.0980]; % AH
colour_mat7 = [0.6350 0.0780 0.1840]; % DH
p = area(t/365,[FOI_AH'./(FOI_AH'+FOI_DH'),FOI_DH'./(FOI_AH'+FOI_DH')]);
p(1).FaceColor = colour_mat3;
p(2).FaceColor = colour_mat2;
legend({'FOI $A_H$','FOI $D_H$'},'Location','se')
title('FOIs (onward inf.)');
xlabel('Time (years)');
ylabel('Proportion')
ylim([0 1])
%% plot3 weighted by surface area: age vs. AH DH FOIs
% figure_setups;
% yyaxis left
% plot(a/365, P.zeta);
% xlabel('Age (years)');
% ylabel('weight')
% yyaxis right
% plot(a/365, DH(:,end),'-', a/365, AH(:,end),':');
% title('weight based on body surface area')
% xlim([0 30])
% legend('weight','DH','AH')
% figure_setups;
% colour_mat3 = [0.9290 0.6940 0.1250]; 
% colour_mat2 = [0.8500 0.3250 0.0980]; 
% colour_mat7 = [0.6350 0.0780 0.1840]; 
% p = area(t/365,[FOI_AHw'./(FOI_AHw'+FOI_DHw'),FOI_DHw'./(FOI_AHw'+FOI_DHw')]);
% p(1).FaceColor = colour_mat3;
% p(2).FaceColor = colour_mat2;
% legend({'FOI $A_H$','FOI $D_H$'},'Location','se')
% title('weighted FOIs (onward inf.)');
% xlabel('Time (years)');
% ylabel('Proportion')
% ylim([0 1])
%% plot2 time vs. AH DH FOIs
figure_setups;
FOI_AH_age = bM(:,end).*P.betaA*AH(:,end)./NH(end); 
FOI_DH_age = bM(:,end).*P.betaD*DH(:,end)./NH(end); 
FOI_AH_agew = FOI_AH_age.*P.zeta;
FOI_DH_agew = FOI_DH_age.*P.zeta;
p = area(a/365,[FOI_AH_agew./(FOI_AH_agew+FOI_DH_agew), FOI_DH_agew./(FOI_AH_agew+FOI_DH_agew)]);
p(1).FaceColor = colour_mat3;
p(2).FaceColor = colour_mat2;
legend({'FOI $A_H$','FOI $D_H$'},'Location','se')
title('FOIs (onward inf.) by age, weighted');
xlabel('Age (years)');
ylabel('Proportion')
xlim([0 30])
ylim([0 1])
%%
figure_setups;
FOI_AH_age = bM(:,end).*P.betaA*AH(:,end)./NH(end); 
FOI_DH_age = bM(:,end).*P.betaD*DH(:,end)./NH(end); 
p = area(a/365,[FOI_AH_age./(FOI_AH_age+FOI_DH_age), FOI_DH_age./(FOI_AH_age+FOI_DH_age)]);
p(1).FaceColor = colour_mat3;
p(2).FaceColor = colour_mat2;
legend({'FOI $A_H$','FOI $D_H$'},'Location','se')
title('FOIs (onward inf.) by age, unweighted');
xlabel('Age (years)');
ylabel('Proportion')
xlim([0 30])
ylim([0 1])
