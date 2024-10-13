%% generate plots for review 1 comment 4
close all;
clear all;
clc;
format long;
global P

% numerical config
tfinal = 10*365; age_max = 100*365; P.age_max = age_max;
dt = 20; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t; P.tfinal = tfinal;

Malaria_parameters_baseline;
Malaria_parameters_transform;
Malaria_parameters_transform_vac;

[SH, EH, DH, AH, ~, ~, SM, EM, IM, ~, ~, ~, Ctot, ~] = age_structured_Malaria_IC_vac('EE_reset');
PH = SH+EH+DH+AH;
NH = trapz(PH)*P.da;
Ctot_pp = Ctot./PH;
var_list = [0.01:0.01:1].^2;
EIR_plot = [1 10 50 80 120];
betaM_plot = NaN(size(EIR_plot));
%% 
xx = P.a/365;
EIR_list = zeros(1,length(var_list));
EH_plot = zeros(na,length(EIR_plot));
AH_plot = zeros(na,length(EIR_plot));
DH_plot = zeros(na,length(EIR_plot));
NewEH_plot = zeros(na,length(EIR_plot));
NewEHDH_plot = zeros(na,length(EIR_plot));
for jj = 1:length(var_list)
    P.betaM = var_list(jj);
    Malaria_parameters_transform;
    [SH, EH, DH, AH, ~, ~, SM, EM, IM, ~, ~, ~, Ctot, ~] = age_structured_Malaria_IC_vac('EE_reset');
    PH = SH+EH+DH+AH;
    NM = SM+EM+IM;
    [bH,~] = biting_rate(PH,NM);
    EIR = bH.*IM./NM*365; % EIR matrix
    EIR_tot = trapz(EIR.*PH)/trapz(PH); % EIR sum over age, at final time
    EIR_list(1,jj) = EIR_tot; % aEIR
end

for iEIR = 1:length(EIR_plot)
    EIR_target = EIR_plot(iEIR);
    [~,ind] = min(abs(EIR_target-EIR_list));
    betaM_plot(iEIR) = var_list(ind);
end

for ibetaM = 1:length(betaM_plot)
    P.betaM = betaM_plot(ibetaM);
    P.v0 = 0; Malaria_parameters_transform_vac;
    [SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');
    tfinal = 365*10;
    [t,SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(da,na,0,tfinal,...
    SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
    PH = SH+EH+DH+AH+VH+UH;
    NM = SM+IM+EM;
    [bH,~] = biting_rate(PH(:,end),NM(:,end));
    lamH = FOI_H(bH,IM(:,end),NM(:,end));
    EH_plot(:,ibetaM) = EH(:,end)./PH(:,end);
    AH_plot(:,ibetaM) = AH(:,end)./PH(:,end);
    DH_plot(:,ibetaM) = DH(:,end)./PH(:,end);    
    NewEH_plot(:,ibetaM) = lamH.*SH(:,end); % SH -> EH

    rho = sigmoid_prob(Ctot(:,end)./PH(:,end), 'rho');
    NewEHDH_plot(:,ibetaM) = rho.*P.h.*EH(:,end); % EH -> DH

end

% baseline
P.betaM = 0.35;
Malaria_parameters_transform;
[SH, EH, DH, AH, ~, ~, SM, EM, IM, ~, ~, ~, Ctot, ~] = age_structured_Malaria_IC_vac('EE_reset');
PH = SH+EH+DH+AH;
NM = SM+EM+IM;
[bH,~] = biting_rate(PH(:,end),NM(:,end));
lamH = FOI_H(bH,IM(:,end),NM(:,end));
EIR = bH.*IM./NM*365; % EIR matrix
EIR_baseline = trapz(EIR.*PH)/trapz(PH); % EIR sum over age, at final time
DH_plot_baseline = DH(:,end)./PH(:,end); 
AH_plot_baseline = AH(:,end)./PH(:,end); 
EH_plot_baseline = EH(:,end)./PH(:,end); 
NewEH_plot_baseline = lamH.*SH(:,end);
rho = sigmoid_prob(Ctot(:,end)./PH(:,end), 'rho');
NewEHDH_plot_baseline = rho.*P.h.*EH(:,end); % EH -> DH
%%
% DH/PH
figure_setups_34;
hold on
plot(xx,DH_plot,'LineWidth',5)
xlim([0 30])
xlabel('Age (years)');
ylabel('Proportion')
title('$D_H/P_H$')
strs = "EIR = " + string(EIR_plot);
legend(strs)
legend('AutoUpdate','on','location','e')
plot(xx,DH_plot_baseline,'m:','LineWidth',5,'DisplayName',['EIR = ',num2str(EIR_baseline,3),'(baseline)'])
saveas(gcf,'plot_baseline_DH.eps','epsc');
%% (AH + DH)/PH
% % figure_setups_34;
% hold on
% plot(xx,AH_plot+DH_plot)
% xlim([0 30])
% xlabel('Age (years)');
% ylabel('Proportion')
% title('$(A_H+D_H)/P_H$')
% strs = "EIR = " + string(EIR_plot);
% legend(strs)
% legend('AutoUpdate','on','location','e')
% plot(xx,DH_plot_baseline+AH_plot_baseline,'m','LineWidth',5,'DisplayName',['EIR = ',num2str(EIR_baseline,3),'(baseline)'])
%% (EH + AH + DH)/PH
figure_setups_34;
hold on
plot(xx,AH_plot+DH_plot+EH_plot)
xlim([0 30])
xlabel('Age (years)');
ylabel('Proportion')
title('$(E_H+A_H+D_H)/P_H$')
strs = "EIR = " + string(EIR_plot);
legend(strs)
legend('AutoUpdate','on','location','e')
plot(xx,DH_plot_baseline+AH_plot_baseline+EH_plot_baseline,'m:','LineWidth',5,'DisplayName',['EIR = ',num2str(EIR_baseline,3),'(baseline)'])
saveas(gcf,'plot_baseline_DH_AH_EH.eps','epsc');
%% New infection: SH -> EH  
% figure_setups_34;
% hold on
% plot(xx,NewEH_plot)
% xlim([0 30])
% xlabel('Age (years)');
% ylabel('Rate')
% title('New infections: $S_H \rightarrow E_H~(\Lambda_H S_H)$')
% strs = "EIR = " + string(EIR_plot);
% legend(strs)
% legend('AutoUpdate','on')
% plot(xx,NewEH_plot_baseline,'m','LineWidth',5,'DisplayName',['EIR = ',num2str(EIR_baseline,3),'(baseline)'])
%% Incidence: EH -> DH 
figure_setups_34;
hold on
plot(xx,NewEHDH_plot)
xlim([0 30])
xlabel('Age (years)');
ylabel('Rate')
title('New uncomplicated incidence: $E_H \rightarrow D_H$')
strs = "EIR = " + string(EIR_plot);
legend(strs)
legend('AutoUpdate','on')
plot(xx,NewEHDH_plot_baseline,'m:','LineWidth',5,'DisplayName',['EIR = ',num2str(EIR_baseline,3),'(baseline)'])
saveas(gcf,'plot_baseline_new_EHDH.eps','epsc');



