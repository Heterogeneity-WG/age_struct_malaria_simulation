clearvars
% close all
clc

global P

% s =  load('Penny_ve_Siaya.mat','Data001'); Data = s.Data001;  Data(1,1) = 0; % point data
% s =  load('Penny_ve_Siaya_model.mat','Penny_ve_Siaya_model'); Data = s.Penny_ve_Siaya_model;  Data(1,1) = 0; % curve data
s =  load('Penny_ve_Siaya_fig2.mat','Penny_ve_Siaya_fig2'); Data = s.Penny_ve_Siaya_fig2; Data(:,1) = Data(:,1)/12;  % Penny point data in figure 2
Data(1,:)=[]; % remove the first data point

%% numerical config
age_max = 100*365; % max ages in days
P.age_max = age_max;
dt = 10; % time/age step size in days, default = 5;
da = dt;
a = (0:da:age_max)';
na = length(a);

P.dt = dt; 
P.a = a;
P.na = na;
P.da = da;

Malaria_parameters_baseline;
Malaria_parameters_transform; 
P.v0s = 0; P.v0c = 0;
Malaria_parameters_transform_vac;

% initial condition
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE');

%% fitting
% small bounds
% lb = [0.6, 1/(5*365)]; % P.etas; P.w
% ub = [1, 1/(0.6*365)]; % x = [0.657384317216716   0.004566197564364]; % 0.642807494823668   0.004565996907337

% lb = [0.6, 1/(5*365)]; % P.etas; P.w
% ub = [1, 1/(0.2*365)]; % x=[0.886936900157632   0.008643665536998] % 0.817146715222409   0.006612625446636 

% lb = [0.1, 1/(5*365)]; % P.etas; P.w
% ub = [1, 1/(0.6*365)];  % x = [0.657378448201810   0.004566148188910] 0.714385219513241   0.004566145122966  % 0.642800169482050   0.004565996823606

% larger bounds
lb = [0.1, 1/(5*365)]; % P.etas; P.w
ub = [1, 1/(0.2*365)]; % x = [0.885804714036094   0.008625159145255]; % 0.817207222138530   0.006613248607787 


% x0 = [P.etas, P.w];
x0 = [0.5, 1/(2*365)];
options = optimset('Display','iter','TolX',10^-8,'MaxIter',80);
[x,fval] = fmincon(@(x) fun_efficacy(x,Data,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0),x0,[],[], [], [], lb, ub, [], options);

%%

figure_setups; hold on
[~, xdata,ydata,yrun0] = fun_efficacy(x0,Data,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0);
[~, ~,~,yrun1] = fun_efficacy(x,Data,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0);
scatter(xdata/365,ydata,'filled')
plot(xdata/365,yrun0,'*-')
plot(xdata/365,yrun1,'*-')
legend('data - Penny','fitting - initial','fitting - final')
axis([0 1550/365 -0.3 0.8])
title('Fit to point data of efficacy')


% figure_setups; hold on
% [~, xdata,ydata,yrun0] = fun_efficacy(x0,Data,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0);
% [~, ~,~,yrun1] = fun_efficacy(x,Data,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0);
% scatter(xdata,ydata,'filled')
% plot(xdata,yrun0,xdata,yrun1)
% legend('data - Penny','fitting - initial','fitting - final')
% axis([0 1550 -0.3 0.8])
% title('Fit to prediction curve, $\eta=0.5590, w=1/(0.3243*365)$')
% % estimates using large bounds, including the 1st data point
