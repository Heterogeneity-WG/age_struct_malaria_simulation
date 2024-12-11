%% generate plots Figure 10, and for review 1 comment 4
close all;
clear all;
clc;
format long;
global P

flag_save = 1; % flag for saving the results or not (Note: it will overwrite previous results in the folder)

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
EIR_plot = [30 80 100];
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
Malaria_parameters_baseline;
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
figure_setups; 
hold on
plot(xx,DH_plot)
xlabel('Age (years)','fontsize',45);
ylabel('Proportion','fontsize',45);
title('$D_H/P_H$','fontsize',45);
strs = "EIR = " + string(EIR_plot);
legend(strs)
legend('AutoUpdate','on','location','ne','fontsize',45);
plot(xx,DH_plot_baseline,'m:','DisplayName',['EIR = ',num2str(EIR_baseline,3)])
xlim([0 30])
ylim([0, 0.32])
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
set(gca,'fontsize', 45) 
ax.GridAlpha = 1;  % Make grid lines transparent..
ax.GridColor = [224, 224, 224]/255; % change grid color
if flag_save; saveas(gcf,'fig_R2_A.eps','epsc'); end
%% (EH + AH + DH)/PH
figure_setups;
hold on
plot(xx,AH_plot+DH_plot+EH_plot)
xlim([0 30])
ylim([0 1.05])
xlabel('Age (years)','fontsize',45);
ylabel('Proportion','fontsize',45);
title('Total infection','fontsize',45);
strs = "EIR = " + string(EIR_plot);
legend(strs)
legend('AutoUpdate','on','location','se','fontsize',45);
plot(xx,DH_plot_baseline+AH_plot_baseline+EH_plot_baseline,'m:','DisplayName',['EIR = ',num2str(EIR_baseline,3)])
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
set(gca,'fontsize', 45) 
ax.GridAlpha = 1;  % Make grid lines transparent..
ax.GridColor = [224, 224, 224]/255; % change grid color
if flag_save; saveas(gcf,'fig_R2_B.eps','epsc'); end
%% Incidence: EH -> DH 
figure_setups;
hold on
plot(xx,NewEHDH_plot)
xlim([0 30])
ylim([0,4.5])
xlabel('Age (years)','fontsize',45);
ylabel('Rate','fontsize',45);
title('$E_H \rightarrow D_H$','fontsize',45);
strs = "EIR = " + string(EIR_plot);
legend(strs)
legend('AutoUpdate','on','location','ne','fontsize',45);
plot(xx,NewEHDH_plot_baseline,'m:','DisplayName',['EIR = ',num2str(EIR_baseline,3)])
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
set(gca,'fontsize', 45) 
ax.GridAlpha = 1;  % Make grid lines transparent..
ax.GridColor = [224, 224, 224]/255; % change grid color
if flag_save; saveas(gcf,'fig_R2_C.eps','epsc'); end



