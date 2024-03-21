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
flag_save = 0; % flag for saving the results or not (Note: it will overwrite previous results in the folder)

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
Malaria_parameters_baseline_Nanoro; % SA based on Nanoro climate profile
Malaria_parameters_transform;
profile_Nanoro = P.gM_fun;
Malaria_parameters_baseline_Siaya; % SA based on Nanoro climate profile
Malaria_parameters_transform;
profile_Siaya = P.gM_fun;
%% plotting
figure_setups; hold on
t = linspace(0,365,100);
plot(t/365,profile_Nanoro(t),'-', 'DisplayName','Nanoro','Linewidth',5);
plot(t/365,profile_Siaya(t),'--','DisplayName','Siaya','Linewidth',5);
xlabel('Month','FontSize',30)
ylabel('Mosquito recruitment rate ($g_M$)','FontSize',35)
xlim([0 1]);
xticks([2 5 8 11]./12)
set(gca,'fontsize',30);
xticklabels({'Mar.','Jun.','Sep.','Dec.'});
legend('Location','ne','FontSize',35)
grid off
if flag_save; saveas(gcf,[direc,'Seasonal_curves_comp.eps'],'epsc'); end
%%
P.v0 = 1.2*10^4;
Malaria_parameters_transform_vac;
figure_setups; hold on
plot(P.a/30,P.v,'-','Linewidth',5);
xlabel('Ages (months)','FontSize',35)
ylabel('Vaccination rate ($\nu$)','FontSize',35)
xlim([0 24])
legend off
grid off
if flag_save; saveas(gcf,[direc,'Vacc_dist.eps'],'epsc'); end

