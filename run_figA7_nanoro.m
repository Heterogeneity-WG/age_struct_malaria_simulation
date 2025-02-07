clear all
% close all
% clc
global P

if ~exist('Results/season_booster','dir')
    disp('Create the folder Results/season_booster');
    disp('Download appropriate data or run "run_vac_booster_season.m" to generate necessary results.');
    return
end

% Data_verylow3 = load('Results/season_booster/season_booster_3_verylow.mat');
Data_low1 = load('Results/season_booster/season_booster_1_low.mat'); % max vax rate 2*10^4
Data_low3 = load('Results/season_booster/season_booster_3_low.mat');
Data_low6 = load('Results/season_booster/season_booster_6_low.mat'); 
% Data_high6 = load('Results/season_booster/season_booster_6_high.mat'); % max vax rate 12*10^4
% Data_low9 = load('Results/season_booster/season_booster_9_low.mat');
% Data_high3 = load('Results/season_booster/season_booster_3_high.mat');
% Data_high6 = load('Results/season_booster/season_booster_6_high.mat');
% load('Results/season_booster/season_booster_EIR.mat','t','EIR_tot')
% load('Results/season_booster/season_booster_mosquitoes.mat','NM')

t0_list = Data_low1.t0_list;
mat_color = [0 0.4470 0.7410;0.9290 0.6940 0.1250; 0.6350 0.0780 0.1840];
plot_year_ind = 3; % study vaccination impact in year two/ten of the program
disp(['Results shown for year ', num2str(plot_year_ind)]);

%% Plotting of cases prevented in the target age range
figure_setups_3;
hold on

plot(t0_list,Data_low1.cases_per_booster_target(:,plot_year_ind),'k-','DisplayName','1-month');
plot(t0_list,Data_low3.cases_per_booster_target(:,plot_year_ind),'k-.','DisplayName','3-month');
plot(t0_list,Data_low6.cases_per_booster_target(:,plot_year_ind),'k:','DisplayName','6-month');

plot(t0_list,Data_low1.cases_per_booster_target_constant(:,plot_year_ind),'r','DisplayName','year-long');
% plot(t0_list,Data_low3.cases_per_booster_target_constant(:,plot_year_ind),'m','DisplayName','year-long');
% plot(t0_list,Data_low6.cases_per_booster_target_constant(:,plot_year_ind),'b','DisplayName','year-long');

xlim([0 12]);
ylim([1.5 2.1]);
xlabel('start of protection');
ylabel('cases prevented per vac ');
legend('Location','northwest');
xticks([2 5 8 11]);
xticklabels({'Mar.','Jun.','Sep.','Dec.'});
title('2 $\sim$ 10 years cohort')
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
% this sets an 'a)' right at the top left of the axes
text(ax,-10,temp_ax(end)+5,'(A)','Units','characters');
save_string = strcat('figA7_nanoro','A','.svg');
saveas(gcf,save_string);

%% Plotting of deaths prevented in the target age range
figure_setups_3;
hold on
plot(t0_list,Data_low1.death_per_booster_target(:,plot_year_ind),'k-','DisplayName','1-month');
plot(t0_list,Data_low3.death_per_booster_target(:,plot_year_ind),'k-.','DisplayName','3-month');
plot(t0_list,Data_low6.death_per_booster_target(:,plot_year_ind),'k:','DisplayName','6-month');

plot(t0_list,Data_low1.death_per_booster_target_constant(:,plot_year_ind),'r','DisplayName','year-long vac');

xlim([0 12]);
ylim([3.2 4.9].*10^(-4));
xlabel('start of protection');
ylabel('deaths prevented per vac');
% legend('Location','southwest');
%yticks([0.3 0.4 0.5 0.6 0.7 0.8]);
xticks([2 5 8 11]);
xticklabels({'Mar.','Jun.','Sep.','Dec.'});
title('2 $\sim$ 10 years cohort')
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
% this sets an 'a)' right at the top left of the axes
text(ax,-10,temp_ax(end)+5,'(B)','Units','characters');
save_string = strcat('figA7_nanoro','B','.svg');
saveas(gcf,save_string);

% max((Data_low1.death_per_booster_target(:,plot_year_ind)-Data_low1.death_per_booster_target_constant(:,plot_year_ind))./Data_low1.death_per_booster_target_constant(:,plot_year_ind))
% max((Data_low3.death_per_booster_target(:,plot_year_ind)-Data_low3.death_per_booster_target_constant(:,plot_year_ind))./Data_low3.death_per_booster_target_constant(:,plot_year_ind))
% max((Data_low6.death_per_booster_target(:,plot_year_ind)-Data_low6.death_per_booster_target_constant(:,plot_year_ind))./Data_low6.death_per_booster_target_constant(:,plot_year_ind))

%% Plotting of cases prevented for the full population
figure_setups_3;
hold on
plot(t0_list,Data_low1.cases_per_booster_full(:,plot_year_ind),'k-','DisplayName','1-month');
plot(t0_list,Data_low3.cases_per_booster_full(:,plot_year_ind),'k-.','DisplayName','3-month');
plot(t0_list,Data_low6.cases_per_booster_full(:,plot_year_ind),'k:','DisplayName','6-month');

plot(t0_list,Data_low1.cases_per_booster_full_constant(:,plot_year_ind),'r','DisplayName','year-long vac');

xlim([0 12]);
ylim([1.5 2.1]);
xlabel('start of protection');
ylabel('cases prevented per vac ');
% legend('Location','southwest');
% yticks([0.3 0.4 0.5 0.6 0.7 0.8]);
xticks([2 5 8 11]);
xticklabels({'Mar.','Jun.','Sep.','Dec.'});
title('Full Population');
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
% this sets an 'a)' right at the top left of the axes
text(ax,-10,temp_ax(end)+5,'(C)','Units','characters');
save_string = strcat('figA7_nanoro','C','.svg');
saveas(gcf,save_string);

%% Plotting of deaths prevented in full population
figure_setups_3;
hold on
plot(t0_list,Data_low1.death_per_booster_full(:,plot_year_ind),'k-','DisplayName','1-month');
plot(t0_list,Data_low3.death_per_booster_full(:,plot_year_ind),'k-.','DisplayName','3-month');
plot(t0_list,Data_low6.death_per_booster_full(:,plot_year_ind),'k:','DisplayName','6-month');

plot(t0_list,Data_low1.death_per_booster_full_constant(:,plot_year_ind),'r','DisplayName','year-long vac');

xlim([0 12]);
ylim([3.2 4.9].*10^(-4));
xlabel('start of protection');
ylabel('deaths prevented per vac');
% legend('Location','southwest');
%yticks([0.3 0.4 0.5 0.6 0.7 0.8]);
xticks([2 5 8 11]);
xticklabels({'Mar.','Jun.','Sep.','Dec.'});
title('Full Population');
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
% this sets an 'a)' right at the top left of the axes
text(ax,-10,temp_ax(end)+5,'(D)','Units','characters');
save_string = strcat('figA7_nanoro','D','.svg');
saveas(gcf,save_string);

% max((Data_low1.death_per_booster_full(:,plot_year_ind)-Data_low1.death_per_booster_full_constant(:,plot_year_ind))./Data_low1.death_per_booster_full_constant(:,plot_year_ind))
% max((Data_low3.death_per_booster_full(:,plot_year_ind)-Data_low3.death_per_booster_full_constant(:,plot_year_ind))./Data_low3.death_per_booster_full_constant(:,plot_year_ind))
% max((Data_low6.death_per_booster_full(:,plot_year_ind)-Data_low6.death_per_booster_full_constant(:,plot_year_ind))./Data_low6.death_per_booster_full_constant(:,plot_year_ind))