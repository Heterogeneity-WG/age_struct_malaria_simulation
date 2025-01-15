%% Figure 5, first panel for various seasonality curves

close all
clear all
clc
format long
global P

% Directory for outputting SA results (sample, data matrices, plots)
% if data is available, the script will load results in the folder;
% otherwise, it will generate new results (could time and storage consuming)
direc = 'Data_SA/Results_local_SA_time/';
flag_save = 1; % flag for saving the results or not (Note: it will overwrite previous results in the folder)

%% numerical config
age_max = 100*365; % max ages in days
P.age_max = age_max;
dt = 1; % time/age step size in days, default = 1;
da = dt;
a = (0:da:age_max)';
na = length(a);

P.dt = dt;
P.a = a;
P.na = na;
P.da = da;

Malaria_parameters_baseline;
Malaria_parameters_baseline_Nanoro; % SA based on Nanoro climate profile
Malaria_parameters_transform;
Malaria_parameters_transform_vac;

%% initial condition 'EE' - numerical EE
P.v0 = 0; Malaria_parameters_transform_vac;
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

%% calculate quantities for plotting
[bH,~] = biting_rate(PH,NM); % bH(age, time) matrix
EIR = bH.*IM./NM*365; % EIR(age, time) matrix
EIR_tot = trapz(EIR.*PH,1)*P.da./NH;
EIR_final = EIR_tot(end);

%% plotting
% seasonlity plots
figure_setups_4;
yyaxis left
hold on 
plot(t/365,EIR_tot,'DisplayName','EIR');
death_rate_tot = trapz(P.muD.*DH,1)*da;
plot(t/365,death_rate_tot,'DisplayName','death incidence rate'); 
ylim([0 140])
ylabel('EIR/death rate')

yyaxis right
hold on
seasonal_profile = P.ss_S;
% plot(t/365, NM,'DisplayName','NM')
% plot(t/365, IM,'DisplayName','IM')
% plot(t/365, seasonal_profile(t),'DisplayName','seasonal profile')
infected_tot = trapz(AH+DH,1)*P.da./NH;
infected_AH =  trapz(AH,1)*P.da./NH;
infected_DH =  trapz(DH,1)*P.da./NH;
plot(t/365,infected_tot,'DisplayName','$A_H+D_H$ prop.')
plot(t/365,infected_AH,'DisplayName','$A_H$ prop.')
plot(t/365,infected_DH,'DisplayName','$D_H$ prop.')
xlabel('years')
xticks([0:3])
ylabel('Proportion')
ylim([0 0.8])

legend('Location','se')
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
% this sets an '(A)' right at the top left of the axes
text(ax,-12,temp_ax(end)+2,'(A)','Units','characters');
% if flag_save; saveas(gcf,[direc,'SA_seasonal_curves_A.eps'],'epsc'); end
if flag_save; saveas(gcf,[direc,'SA_seasonal_curves_A.svg']); end
