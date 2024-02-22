clear all
close all
clc
global P

if ~exist('Results/season_vacc','dir')
    disp('Create the folder Results/season_vacc');
    disp('Download appropriate data or run "run_vacc_season.m" to generate necessary results.');
    return
end

Data_low1 = load('Results/season_vacc/season_vacc_1_low.mat');
Data_low3 = load('Results/season_vacc/season_vacc_3_low.mat');
Data_low6 = load('Results/season_vacc/season_vacc_6_low.mat');
Data_high3 = load('Results/season_vacc/season_vacc_3_high.mat');
Data_high6 = load('Results/season_vacc/season_vacc_6_high.mat');
%load('Results/season_vacc/season_vacc_EIR.mat','t','EIR_tot')
%load('Results/season_vacc/season_vacc_mosquitoes.mat','NM')

%% Generate mosquito and EIR data with no vaccination for comparison purposes
% numerical config
tfinal = 1*365; age_max = 100*365; P.age_max = age_max;
dt = 5; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t; P.tfinal = tfinal;

Malaria_parameters_baseline;
Malaria_parameters_transform;
Malaria_parameters_transform_vac;
P.v0 = 0;
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');
[t,SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(da,na,0,tfinal,...
    SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
PH = SH+EH+DH+AH+VH+UH;
PH_final = PH(:,end);
NH = trapz(PH,1)*da;
NM = SM+IM+EM;
[bH,~] = biting_rate(PH,NM); % bH(age, time) matrix
EIR = bH.*IM./NM*365; % EIR(age, time) matrix
EIR_tot = trapz(EIR.*PH,1)*P.da./NH;
EIR_final = EIR_tot(end);

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
h = figure_setups_3;
hold on
ind = 1:10:length(t);
yyaxis left
p1 = plot(t/30,NM,'-^','MarkerIndices',ind,'DisplayName','Mosquitoes');
ylim([0 3*10^8])
yticks([0 1*10^8 2*10^8 3*10^8]);
ylabel('Mosquitoes')
yyaxis right
p2 = plot(t/30,EIR_tot,'-o','MarkerIndices',ind,'DisplayName','EIR');
ylabel('EIR');
ylim([0 120])
yticks([0 20 40 60 80 100 120]);
xticks([0 2 4 6 8 10 12]);
sz = 100;
scatter([t1,t2,t3],[80,70,60],sz,'k','filled','Marker','<');scatter([t11,t22,t33],[80,70,60],sz,'k','filled','Marker','>')
f1 = plot([t1,t11],[80, 80],'k-','DisplayName','1-month');
f2 = plot([t2,t22],[70, 70],'k--','DisplayName','3-month');
f3 = plot([t3,t33],[60, 60],'k:','DisplayName','6-month');
xlim([0 12])
ll = legendUnq(h);
xlabel('Month')
xticks([2 5 8 11])
xticklabels({'Feb.','May','Aug.','Nov.'});
legend(ll,'Location','southeast');
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
% this sets an 'a)' right at the top left of the axes
text(ax,0,temp_ax(end)+3,'(B)','Units','characters');
save_string = strcat('fig6_','B','.svg');
saveas(gcf,save_string);
%%
figure_setups_3;
hold on
plot(t0_list,Data_low1.cases_per_vacc_target(:,plot_year_ind),'k-','DisplayName','1-month');
plot(t0_list,Data_low3.cases_per_vacc_target(:,plot_year_ind),'k-.','DisplayName','3-month');
plot(t0_list,Data_low6.cases_per_vacc_target(:,plot_year_ind),'k:','DisplayName','6-month');
plot(t0_list,Data_low1.cases_per_vacc_target_constant(:,plot_year_ind),'r','DisplayName','year-long vac');
xlim([0 12])
ylim([0.35 0.7])
xlabel('Vac. start month');
ylabel('Cases prevented per vac ')
legend('Location','southwest');
yticks([0.3 0.4 0.5 0.6 0.7 0.8]);
xticks([2 5 8 11])
xticklabels({'Feb.','May','Aug.','Nov.'});
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
% this sets an 'a)' right at the top left of the axes
text(ax,0,temp_ax(end)+3,'(A)','Units','characters');
save_string = strcat('fig6_','A','.svg');
saveas(gcf,save_string);
%%
h = figure_setups_3;
yyaxis left
hold on
ind = 1:4:length(t0_list);
plot(t0_list,Data_low3.cases_per_vacc_target(:,plot_year_ind),':^','DisplayName','low vac count (3 months)','MarkerIndices',ind)
plot(t0_list,Data_high3.cases_per_vacc_target(:,plot_year_ind),'--o','DisplayName','high vac count (3 months)','MarkerIndices',min(ind+2,length(t0_list)))
ylim([0.35 0.8])
ylabel('Cases prevented per vac');
yyaxis right
plot(t0_list,Data_low3.cases_pp_py_target(:,plot_year_ind),':^','MarkerIndices',ind)
plot(t0_list,Data_high3.cases_pp_py_target(:,plot_year_ind),'--o','MarkerIndices',min(ind+2,length(t0_list)))
ylim([2.35 2.6]);
xlim([0 12]);
xticks([2 5 8 11])
xticklabels({'Feb.','May','Aug.','Nov.'});
ylabel('Cases per person');
ll = legendUnq(h);
xlabel('Starting month of vaccination');
legend(ll,'Location','southeast');
%xticklabels(month)
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
% this sets an 'a)' right at the top left of the axes
text(ax,0,temp_ax(end)+3,'(C)','Units','characters');
save_string = strcat('fig6_','C','.svg');
saveas(gcf,save_string);

%% plotting (target pop)
% plot_year_ind = [1,2,10];
% mat_color = [0 0.4470 0.7410;0.9290 0.6940 0.1250; 0.6350 0.0780 0.1840];
% figure_setups;  % plot cases prevented per vacc
% h = subplot(1,2,1);
% f1 = plot(t0_list,Data_low1.cases_per_vacc_target_constant(:,plot_year_ind),'-','DisplayName','year-long (target pop)');
% hold on; grid on
% f2 = plot(t0_list,Data_low1.cases_per_vacc_target(:,plot_year_ind),'--','DisplayName',' seasonal (target pop)');
% f3 = plot(t0_list,Data_low1.cases_per_vacc_target_baseline(:,plot_year_ind),':','DisplayName','no vacc (target pop)');
% xlabel('vacc start month')
% ll = legendUnq(h);
% legend(ll,'Location','best')
% ylabel('cases prev. per vacc')
% xlim([0 12])
% ylim([0 0.8])
% set(f1, {'color'}, num2cell(mat_color,2));
% set(f2, {'color'}, num2cell(mat_color,2));
% set(f3, {'color'}, num2cell(mat_color,2));
%
% h = subplot(1,2,2); % plot death prevented per vacc
% f1 = plot(t0_list,Data_low1.death_per_vacc_target_constant(:,plot_year_ind),'-','DisplayName','year-long (target pop)');
% hold on; grid on
% f2 = plot(t0_list,Data_low1.death_per_vacc_target(:,plot_year_ind),'--','DisplayName','seasonal (target pop)');
% f3 = plot(t0_list,Data_low1.death_per_vacc_target_baseline(:,plot_year_ind),':','DisplayName',' no vacc (target pop)');
% xlabel('vacc start month')
% ll = legendUnq(h);
% legend(ll,'Location','best')
% ylabel('death prev. per vac')
% xlim([0 12])
% ylim([0 1.5/1000])
% set(f1, {'color'}, num2cell(mat_color,2));
% set(f2, {'color'}, num2cell(mat_color,2));
% set(f3, {'color'}, num2cell(mat_color,2));

%%
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