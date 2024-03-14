%% Distributions of Length Infectious and Pre Infectious
clear all; close all; clc

%% Treated MT data

treat = readtable('Treated.csv');
lenDH = treat.lenDH;
lenAH = treat.lenAH;
lenh = treat.lenh;

%% Plot Treated
ft1 = figure_setups_3;
%set(ft,'Position',[100 100 1400 500])
%subplot(1,3,1)
hT_DH = histogram(lenDH);
%xlabel({'days infectious','(with symptoms)'})
xlabel({'Days infectious (w/ symp.)'});
ylabel('Frequency');
%title('Infectious with symptoms')
axis square
ax = gca;
ymax = ax.YLim(2);
xmax = ax.XLim(2)*.35;
text(xmax,ymax*.92,sprintf('mean = %d days',round(nanmean(lenDH))),'fontsize',40,'Interpreter','latex');
text(xmax,ymax*.86,sprintf('n = %d ',sum(hT_DH.Values)),'fontsize',40);
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
% this sets an 'a)' right at the top left of the axes
text(ax,-12,temp_ax(end)+3,'(A)','Units','characters');
save_string = strcat('fig7_','A','.svg');
saveas(gcf,save_string);

ft2 = figure_setups_3;
hT_AH = histogram(lenAH);
%xlabel({'days infectious','(without symptoms)'})
xlabel({'Days infectious (w/o symp.)'})
ylabel('Frequency')
axis square
ax = gca;
ymax = ax.YLim(2);
xmax = ax.XLim(2)*.35;
text(xmax,ymax*.92,sprintf('mean = %d days',round(nanmean(lenAH))),'fontsize',40)
text(xmax,ymax*.86,sprintf('n = %d ',sum(hT_AH.Values)),'fontsize',40)
%title('Infectious without symptoms')
th = title({'Treated during parasitemia (and before 30 days)'},'fontsize',40);
set(th,'Position',[th.Position(1) th.Position(2)*1.025])
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
% this sets an 'a)' right at the top left of the axes
text(ax,-12,temp_ax(end)+3,'(B)','Units','characters');
save_string = strcat('fig7_','B','.svg');
saveas(gcf,save_string);

ft3 = figure_setups_3;
hT_h = histogram(lenh);
%xlabel({'days prior to','infectiousness'})
xlabel({'Days before infectious'})
ylabel('Frequency')
axis square
ax = gca;
ymax = ax.YLim(2);
xmax = ax.XLim(1)+(ax.XLim(2)-ax.XLim(1))*.35;
text(xmax,ymax*.92,sprintf('mean =  %d days',round(nanmean(lenh))),'fontsize',40)
text(xmax,ymax*.86,sprintf('n = %d ',sum(hT_h.Values)),'fontsize',40)
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
% this sets an 'a)' right at the top left of the axes
text(ax,-12,temp_ax(end)+3,'(C)','Units','characters');
%title('Prepatent + time to infectious')
save_string = strcat('fig7_','C','.svg');
saveas(gcf,save_string);

%sgtitle('Treated during parasitemia (and before 30 days)','fontsize',36)

%% Untreated MT Data

untreat = readtable('Untreated.csv');
lenDH = untreat.lenDH;
lenAH = untreat.lenAH;
lenh = untreat.lenh;

%% Plot Untreated

ft4 = figure_setups_3;
%set(ft,'Position',[100 100 1400 500])
%subplot(1,3,1)
hT_DH = histogram(lenDH);
%xlabel({'days infectious','(with symptoms)'})
xlabel({'Days infectious (w/ symp.)'});
ylabel('Frequency');
%title('Infectious with symptoms')
axis square
ax = gca;
ymax = ax.YLim(2);
xmax = ax.XLim(2)*.35;
text(xmax,ymax*.92,sprintf('mean = %d days',round(nanmean(lenDH))),'fontsize',40)
text(xmax,ymax*.86,sprintf('n = %d ',sum(hT_DH.Values)),'fontsize',40)
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
% this sets an 'a)' right at the top left of the axes
text(ax,-12,temp_ax(end)+3,'(D)','Units','characters');
save_string = strcat('fig7_','D','.svg');
saveas(gcf,save_string);

ft5 = figure_setups_3;
hT_AH = histogram(lenAH);
%xlabel({'days infectious','(without symptoms)'})
xlabel({'Days infectious (w/o symp.)'})
ylabel('Frequency')
axis square
ax = gca;
ymax = ax.YLim(2);
xmax = ax.XLim(2)*.35;
text(xmax,ymax*.92,sprintf('mean = %d days',round(nanmean(lenAH))),'fontsize',40)
text(xmax,ymax*.86,sprintf('n = %d ',sum(hT_AH.Values)),'fontsize',40)
%title('Infectious without symptoms')
th = title('Not treated during parasitemia','fontsize',40);
set(th,'Position',[th.Position(1) th.Position(2)*1.025])
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
% this sets an 'a)' right at the top left of the axes
text(ax,-12,temp_ax(end)+3,'(E)','Units','characters');
save_string = strcat('fig7_','E','.svg');
saveas(gcf,save_string);

ft6 = figure_setups_3;
hT_h = histogram(lenh);
%xlabel({'days prior to','infectiousness'})
xlabel({'Days before infectious'})
ylabel('Frequency')
axis square
ax = gca;
ymax = ax.YLim(2);
xmax = ax.XLim(1)+(ax.XLim(2)-ax.XLim(1))*.35;
text(xmax,ymax*.92,sprintf('mean =  %d days',round(nanmean(lenh))),'fontsize',40)
text(xmax,ymax*.86,sprintf('n = %d ',sum(hT_h.Values)),'fontsize',40)
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
% this sets an 'a)' right at the top left of the axes
text(ax,-12,temp_ax(end)+3,'(F)','Units','characters');
%title('Prepatent + time to infectious')
save_string = strcat('fig7_','F','.svg');
saveas(gcf,save_string);

%sgtitle('Not treated during parasitemia','fontsize',36)
