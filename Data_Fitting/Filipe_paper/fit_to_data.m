% return an interpolated surface function as data for fitting
% F(age(years), EIR) = susceptibility (~ rho) from Filipe's paper, Fig 7

close all
clear all
clc
x1 = load('EIR_1.mat'); 
x2 = load('EIR_10.mat');
x3 = load('EIR_20.mat');
x4 = load('EIR_50.mat');
x5 = load('EIR_100.mat');
Data_1 = [x1.EIR_1(:,1), 1*ones(size(x1.EIR_1(:,1))), x1.EIR_1(:,2); 150, 1, x1.EIR_1(end,2)];
Data_2 = [x2.EIR_10(:,1), 10*ones(size(x2.EIR_10(:,1))), x2.EIR_10(:,2); 150, 10, x2.EIR_10(end,2)];
Data_3 = [x3.EIR_20(:,1), 20*ones(size(x3.EIR_20(:,1))), x3.EIR_20(:,2); 150, 20, x3.EIR_20(end,2)];
Data_4 = [x4.EIR_50(:,1), 50*ones(size(x4.EIR_50(:,1))), x4.EIR_50(:,2); 150, 50, x4.EIR_50(end,2)];
Data_5 = [x5.EIR_100(:,1), 100*ones(size(x5.EIR_100(:,1))), x5.EIR_100(:,2); 150, 100, x5.EIR_100(end,2)];
Data_mat = [Data_1;Data_2;Data_3;Data_4;Data_5];
x = Data_mat(:,1); y = Data_mat(:,2); v = Data_mat(:,3);
F = scatteredInterpolant(x,y,v,'linear','nearest');
xq = [0.1:0.2:0.9,1:0.1:10,11:1:100];
yq = 1:1:150;
[xx,yy] = meshgrid([0.1:0.2:0.9,1:0.1:10,11:1:100],1:1:150);
vq = F(xx,yy); vq_sm = vq;
for ind_y = 1:length(yq)
    vq_sm(:,ind_y) = smooth(vq(:,ind_y),'lowess');
end
figure_setups; hold on;
% plot3(x,y,v,'.')
contourf(xx,yy,vq)
% surf(xq,yq,vq)
axis([0 50 0 100])
figure_setups; hold on;
contourf(xx,yy,vq_sm)
axis([0 50 0 100])
% save('Filipe_paper/F_Filipe.mat','F'); 
%% reproduce Filipe's curves in paper (PLoS Computational Biology 2007)
figure_setups; hold on
xx = 0.3:0.2:50; yy = [1 10 20 50 100];
for iEIR = 1:5
    zz = smooth(F(xx,yy(iEIR)*ones(size(xx))),3);
    plot(xx,zz,'DisplayName',['EIR=',num2str(yy(iEIR))])
end
legend('Location','e')
xlabel('age (years)')
ylabel('susceptibility')
axis([0 50 0 1.05 ])




