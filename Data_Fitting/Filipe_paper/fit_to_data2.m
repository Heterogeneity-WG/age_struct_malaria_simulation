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
xq = [0.1:0.2:0.9,1:0.1:10,11:1:100]'; 
yq = [1:1:150]';
Data_1p = [xq, 1*ones(size(xq)), interp1(Data_1(:,1),Data_1(:,3),xq,'linear','extrap')];
Data_2p = [xq, 10*ones(size(xq)), interp1(Data_2(:,1),Data_2(:,3),xq,'linear','extrap')];
Data_3p = [xq, 20*ones(size(xq)), interp1(Data_3(:,1),Data_3(:,3),xq,'linear','extrap')];
Data_4p = [xq, 50*ones(size(xq)), interp1(Data_4(:,1),Data_4(:,3),xq,'linear','extrap')];
Data_5p = [xq, 100*ones(size(xq)), interp1(Data_5(:,1),Data_5(:,3),xq,'linear','extrap')];
Data_matp = [Data_1p;Data_2p;Data_3p;Data_4p;Data_5p];
xp = Data_matp(:,1); yp = Data_matp(:,2); vp = Data_matp(:,3);
F = scatteredInterpolant(xp,yp,vp,'linear','nearest');
% save('Data_Fitting/Filipe_paper/F_Filipe.mat','F'); 

xx = [0.1:0.05:20];
yy = [0:0.05:120];
[XX,YY] = meshgrid(xx,yy);
vq = F(XX,YY); 
figure_setups; hold on;
% plot3(x,y,v,'.')
imagesc(xx,yy,vq)
set(gca,'YDir','normal');
colormap jet
colorbar('Ticks',0:0.1:1);
axis([0 20 0 120])
% contourf(XX,YY,vq)
% surf(xq,yq,vq)
% axis([0 50 0 100])

%% reproduce Filipe's curves in paper (PLoS Computational Biology 2007)
% figure_setups; hold on
% xx_f = [0.3:0.2:0.9,1:0.1:10,11:1:100]';  yy_f = [1 10 20 50 100];
% for iEIR = 1:5
%     zz_f = F(xx_f,yy_f(iEIR)*ones(size(xx_f)));
%     plot(xx_f,zz_f,'DisplayName',['EIR=',num2str(yy_f(iEIR))])
% end
% legend('Location','e')
% xlabel('age (years)')
% ylabel('susceptibility')
% axis([0 50 0 1.05 ])




