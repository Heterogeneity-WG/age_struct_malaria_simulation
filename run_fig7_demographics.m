%% vac = constant, fixed number
% close all
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
Malaria_parameters_transform;
direc = 'Data_SA/Results_local_SA_time/';
flag_save = 0;
%% plot calibrated fertility & death
figure_setups; hold on
grid off
plot(a/365, P.gH,':', 'DisplayName','Birth rate ($g_H$)','Linewidth',5)
plot(a/365, P.muD,'-', 'DisplayName','Malaria-induced mortality rate ($\mu_D$)','Linewidth',5);
plot(a/365, P.muH,'--','DisplayName','Non-malaria mortality rate ($\mu_H$)','Linewidth',5);
xlabel('Age (years)','FontSize',35)
ylabel('Rates','FontSize',35)
legend('Location','nw','FontSize',28)
ylim([-0.1 4]*10^(-4))
xlim([0 90])
if flag_save; saveas(gcf,[direc,'Demo_curves.eps'],'epsc'); end




