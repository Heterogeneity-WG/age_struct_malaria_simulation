clear all
close all
clc

Data_low1 = load('Results/season_vacc/season_vacc_1_low.mat');
Data_low3 = load('Results/season_vacc/season_vacc_3_low.mat');
Data_low6 = load('Results/season_vacc/season_vacc_6_low.mat');
Data_high3 = load('Results/season_vacc/season_vacc_3_high.mat');
Data_high6 = load('Results/season_vacc/season_vacc_6_high.mat');
load('Results/season_vacc/season_vacc_EIR.mat','t','EIR_tot')
load('Results/season_vacc/season_vacc_mosquitoes.mat','NM')
%%
t0_list = Data_low1.t0_list;
mat_color = [0 0.4470 0.7410;0.9290 0.6940 0.1250; 0.6350 0.0780 0.1840];
plot_year_ind = 10;

[Opt1,ind1] = max(Data_low1.cases_per_vacc_target(:,plot_year_ind));
[Opt2,ind2] = max(Data_low3.cases_per_vacc_target(:,plot_year_ind));
[Opt3,ind3]= max(Data_low6.cases_per_vacc_target(:,plot_year_ind));
t1 = t0_list(ind1); t11 = t1+1;
t2 = t0_list(ind2); t22 = t2+3;
t3 = t0_list(ind3); t33 = t3+6;
%%
% h = figure_setups; 
% hold on
% ind = 1:10:length(t);
% yyaxis left
% p1 = plot(t/30,NM,'-^','MarkerIndices',ind,'DisplayName','Mosquitoes');
% ylim([0 2*10^8])
% ylabel('\# of mosquitoes')
% yyaxis right
% p2 = plot(t/30,EIR_tot,'-o','MarkerIndices',ind,'DisplayName','EIR');
% ylabel('EIR')
% sz = 100;
% scatter([t1,t2,t3],[80,70,60],sz,'k','filled','Marker','<');scatter([t11,t22,t33],[80,70,60],sz,'k','filled','Marker','>')
% f1 = plot([t1,t11],[80, 80],'k-','DisplayName','1-month');
% f2 = plot([t2,t22],[70, 70],'k--','DisplayName','3-month');
% f3 = plot([t3,t33],[60, 60],'k-.','DisplayName','6-month');
% xlim([0 12])
% ll = legendUnq(h);
% xlabel('Month')
% legend(ll,'Location','southeast')
% saveas(gcf,'Results/season_vacc/vac_optimal.eps','epsc')
% %%
% figure_setups;
% hold on
% plot(t0_list,Data_low1.cases_per_vacc_target(:,plot_year_ind),'k-','DisplayName','1-month');
% plot(t0_list,Data_low3.cases_per_vacc_target(:,plot_year_ind),'k--','DisplayName','3-month');
% plot(t0_list,Data_low6.cases_per_vacc_target(:,plot_year_ind),'k-.','DisplayName','6-month');
% plot(t0_list,Data_low1.cases_per_vacc_target_constant(:,plot_year_ind),'k:','DisplayName','year-long vac');
% xlim([0 12])
% ylim([0.3 0.8])
% xlabel('Starting month of vaccination')
% ylabel('Cases prevented per year per vac ')
% legend('Location','northeast')
% saveas(gcf,'Results/season_vacc/vac_cases_pp.eps','epsc')
% %%
% h = figure_setups;
% yyaxis left
% hold on
% ind = 1:4:length(t0_list);
% plot(t0_list,Data_low3.cases_per_vacc_target(:,plot_year_ind),':^','DisplayName','low vac count ($1.2\times 10^4$/year)','MarkerIndices',ind)
% plot(t0_list,Data_high3.cases_per_vacc_target(:,plot_year_ind),'--o','DisplayName','high vac count ($6\times 10^4$/year)','MarkerIndices',min(ind+2,length(t0_list)))
% ylim([0.35 0.8])
% ylabel('Cases prevented per year per vac ')
% yyaxis right
% plot(t0_list,Data_low3.cases_pp_py_target(:,plot_year_ind),':^','DisplayName','low vac count ($1.2\times 10^4$/year)','MarkerIndices',ind)
% plot(t0_list,Data_high3.cases_pp_py_target(:,plot_year_ind),'--o','DisplayName','high vac count ($6\times 10^4$/year)','MarkerIndices',min(ind+2,length(t0_list)))
% ylim([2.35 2.6])
% ylabel('Cases per person per year')
% ll = legendUnq(h);
% xlabel('Starting month of vaccination')
% legend(ll,'Location','northeast')
% saveas(gcf,'Results/season_vacc/vac_cases_counts.eps','epsc')

%% plotting (target pop)
plot_year_ind = [1,2,10];
mat_color = [0 0.4470 0.7410;0.9290 0.6940 0.1250; 0.6350 0.0780 0.1840];
figure_setups;  % plot cases prevented per vacc
h = subplot(1,2,1);
f1 = plot(t0_list,Data_low1.cases_per_vacc_target_constant(:,plot_year_ind),'-','DisplayName','year-long (target pop)');
hold on; grid on
f2 = plot(t0_list,Data_low1.cases_per_vacc_target(:,plot_year_ind),'--','DisplayName',' seasonal (target pop)');
f3 = plot(t0_list,Data_low1.cases_per_vacc_target_baseline(:,plot_year_ind),':','DisplayName','no vacc (target pop)');
xlabel('vacc start month')
ll = legendUnq(h);
legend(ll,'Location','best')
ylabel('cases prev. per vacc')
xlim([0 12])
ylim([0 0.8])
set(f1, {'color'}, num2cell(mat_color,2));
set(f2, {'color'}, num2cell(mat_color,2));
set(f3, {'color'}, num2cell(mat_color,2));

h = subplot(1,2,2); % plot death prevented per vacc
f1 = plot(t0_list,Data_low1.death_per_vacc_target_constant(:,plot_year_ind),'-','DisplayName','year-long (target pop)');
hold on; grid on
f2 = plot(t0_list,Data_low1.death_per_vacc_target(:,plot_year_ind),'--','DisplayName','seasonal (target pop)');
f3 = plot(t0_list,Data_low1.death_per_vacc_target_baseline(:,plot_year_ind),':','DisplayName',' no vacc (target pop)');
xlabel('vacc start month')
ll = legendUnq(h);
legend(ll,'Location','best')
ylabel('death prev. per vac')
xlim([0 12])
ylim([0 1.5/1000])
set(f1, {'color'}, num2cell(mat_color,2));
set(f2, {'color'}, num2cell(mat_color,2));
set(f3, {'color'}, num2cell(mat_color,2));

% h = subplot(1,3,3);
% f1 = plot(t0_list,Data_low1.cases_pp_py_target_constant(:,plot_year_ind),'-','DisplayName','year-long (target pop)');
% hold on; grid on
% f2 = plot(t0_list,Data_low1.cases_pp_py_target(:,plot_year_ind),'--','DisplayName','seasonal (target pop)');
% f3 = plot(t0_list,Data_low1.cases_pp_py_target_baseline(:,plot_year_ind),':','DisplayName',' no vacc (target pop)');
% xlabel('vacc start month')
% ll = legendUnq(h);
% legend(ll,'Location','best')
% ylabel('cases per person per year')
% xlim([0 12])
% ylim([2.4 2.65])
% set(f1, {'color'}, num2cell(mat_color,2));
% set(f2, {'color'}, num2cell(mat_color,2));
% set(f3, {'color'}, num2cell(mat_color,2));