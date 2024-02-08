%% Distributions of Length Infectious and Pre Infectious
clear all; close all; clc

%% Treated MT data

treat = readtable('Treated.csv');
lenDH = treat.lenDH;
lenAH = treat.lenAH;
lenh = treat.lenh;

%% Plot Treated
ft = figure_setups();
set(ft,'Position',[100 100 1400 500])
subplot(1,3,1)
hT_DH = histogram(lenDH);
xlabel({'days infectious','(with symptoms)'})
ylabel('frequency')
%title('Infectious with symptoms')
axis square
ax = gca;
ymax = ax.YLim(2);
xmax = ax.XLim(2)*.5;
text(xmax,ymax*.92,sprintf('mean = %d days',round(nanmean(lenDH))),'fontsize',20)
text(xmax,ymax*.86,sprintf('n = %d ',sum(hT_DH.Values)),'fontsize',20)

subplot(1,3,2)
hT_AH = histogram(lenAH);
xlabel({'days infectious','(without symptoms)'})
ylabel('frequency')
axis square
ax = gca;
ymax = ax.YLim(2);
xmax = ax.XLim(2)*.5;
text(xmax,ymax*.92,sprintf('mean = %d days',round(nanmean(lenAH))),'fontsize',20)
text(xmax,ymax*.86,sprintf('n = %d ',sum(hT_AH.Values)),'fontsize',20)
%title('Infectious without symptoms')
th = title({'Treated during parasitemia (and before 30 days)'},'fontsize',36);
set(th,'Position',[th.Position(1) th.Position(2)*1.05])

subplot(1,3,3)
hT_h = histogram(lenh);
xlabel({'days prior to','infectiousness'})
ylabel('frequency')
axis square
ax = gca;
ymax = ax.YLim(2);
xmax = ax.XLim(1)+(ax.XLim(2)-ax.XLim(1))*.5;
text(xmax,ymax*.92,sprintf('mean =  %d days',round(nanmean(lenh))),'fontsize',20)
text(xmax,ymax*.86,sprintf('n = %d ',sum(hT_h.Values)),'fontsize',20)
%title('Prepatent + time to infectious')

%sgtitle('Treated during parasitemia (and before 30 days)','fontsize',36)

%% Untreated MT Data

untreat = readtable('Untreated.csv');
lenDH = untreat.lenDH;
lenAH = untreat.lenAH;
lenh = untreat.lenh;

%% Plot Untreated

ft = figure_setups();
set(ft,'Position',[100 100 1400 500])
subplot(1,3,1)
hT_DH = histogram(lenDH);
xlabel({'days infectious','(with symptoms)'})
ylabel('frequency')
%title('Infectious with symptoms')
axis square
ax = gca;
ymax = ax.YLim(2);
xmax = ax.XLim(2)*.5;
text(xmax,ymax*.92,sprintf('mean = %d days',round(nanmean(lenDH))),'fontsize',20)
text(xmax,ymax*.86,sprintf('n = %d ',sum(hT_DH.Values)),'fontsize',20)

subplot(1,3,2)
hT_AH = histogram(lenAH);
xlabel({'days infectious','(without symptoms)'})
ylabel('frequency')
axis square
ax = gca;
ymax = ax.YLim(2);
xmax = ax.XLim(2)*.5;
text(xmax,ymax*.92,sprintf('mean = %d days',round(nanmean(lenAH))),'fontsize',20)
text(xmax,ymax*.86,sprintf('n = %d ',sum(hT_AH.Values)),'fontsize',20)
%title('Infectious without symptoms')
th = title('Not treated during parasitemia','fontsize',36);
set(th,'Position',[th.Position(1) th.Position(2)*1.05])

subplot(1,3,3)
hT_h = histogram(lenh);
xlabel({'days prior to','infectiousness'})
ylabel('frequency')
axis square
ax = gca;
ymax = ax.YLim(2);
xmax = ax.XLim(1)+(ax.XLim(2)-ax.XLim(1))*.5;
text(xmax,ymax*.92,sprintf('mean =  %d days',round(nanmean(lenh))),'fontsize',20)
text(xmax,ymax*.86,sprintf('n = %d ',sum(hT_h.Values)),'fontsize',20)
%title('Prepatent + time to infectious')

%sgtitle('Not treated during parasitemia','fontsize',36)
