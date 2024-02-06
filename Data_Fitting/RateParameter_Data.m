%% Distributions of Length Infectious and Pre Infectious
clear all; close all; clc

%% Treated MT data

treat = readtable('Treated.csv');
lenDH = treat.lenDH;
lenAH = treat.lenAH;
lenh = treat.lenh;

%% Plot Treated
ft = figure_setups();
set(ft,'Position',[100 100 1400 400])
subplot(1,3,1)
hT_DH = histogram(lenDH);
xlabel('day')
ylabel('frequency')
title('Infectious with symptoms')
ymax = max(hT_DH.Values);
xmax = hT_DH.BinEdges(end)*.5;
text(xmax,ymax*.9,sprintf('mean = %d days',round(nanmean(lenDH))),'fontsize',20)
text(xmax,ymax*.8,sprintf('n = %d ',sum(hT_DH.Values)),'fontsize',20)

subplot(1,3,2)
hT_AH = histogram(lenAH);
xlabel('day')
ylabel('frequency')
ymax = max(hT_AH.Values);
xmax = hT_AH.BinEdges(end)*.5;
text(xmax,ymax*.9,sprintf('mean = %d days',round(nanmean(lenAH))),'fontsize',20)
text(xmax,ymax*.8,sprintf('n = %d ',sum(hT_AH.Values)),'fontsize',20)
title('Infectious without symptoms')

subplot(1,3,3)
hT_h = histogram(lenh);
xlabel('days')
ylabel('frequency')
ymax = max(hT_h.Values);
xmax = hT_h.BinEdges(end)*.65;
text(xmax,ymax*.9,sprintf('mean =  %d days',round(nanmean(lenh))),'fontsize',20)
text(xmax,ymax*.8,sprintf('n = %d ',sum(hT_h.Values)),'fontsize',20)
title('Prepatent + time to infectious')

sgtitle('Treated during parasitemia (and before 30 days)','fontsize',36)

%% Untreated MT Data

untreat = readtable('Untreated.csv');
lenDH = untreat.lenDH;
lenAH = untreat.lenAH;
lenh = untreat.lenh;

%% Plot Untreated

ft = figure_setups();
set(ft,'Position',[100 100 1400 400])
subplot(1,3,1)
hT_DH = histogram(lenDH);
xlabel('day')
ylabel('frequency')
title('Infectious with symptoms')
ymax = max(hT_DH.Values);
xmax = hT_DH.BinEdges(end)*.5;
text(xmax,ymax*.9,sprintf('mean = %d days',round(nanmean(lenDH))),'fontsize',20)
text(xmax,ymax*.8,sprintf('n = %d ',sum(hT_DH.Values)),'fontsize',20)

subplot(1,3,2)
hT_AH = histogram(lenAH);
xlabel('day')
ylabel('frequency')
ymax = max(hT_AH.Values);
xmax = hT_AH.BinEdges(end)*.5;
text(xmax,ymax*.9,sprintf('mean = %d days',round(nanmean(lenAH))),'fontsize',20)
text(xmax,ymax*.8,sprintf('n = %d ',sum(hT_AH.Values)),'fontsize',20)
title('Infectious without symptoms')

subplot(1,3,3)
hT_h = histogram(lenh);
xlabel('days')
ylabel('frequency')
ymax = max(hT_h.Values);
xmax = hT_h.BinEdges(end)*.65;
text(xmax,ymax*.9,sprintf('mean =  %d days',round(nanmean(lenh))),'fontsize',20)
text(xmax,ymax*.8,sprintf('n = %d ',sum(hT_h.Values)),'fontsize',20)
title('Prepatent + time to infectious')

sgtitle('Not treated during parasitemia','fontsize',36)
