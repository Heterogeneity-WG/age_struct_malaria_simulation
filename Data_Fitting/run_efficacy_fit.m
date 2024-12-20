clearvars
close all
clc

global P

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
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');

%% fitting
% smaller bounds
% lb = [0.6, 1/(1*365)]; % P.etas; P.w
% ub = [1, 1/(0.5*365)]; 

% larger bounds
lb = [0.1, 1/(5*365)]; % P.etas; P.w
ub = [1, 1/(0.1*365)]; 

% x0 = [0.725, 1/(0.95*365)]; % estimates from Penny et al, Table 3 (half-life = 7.9 month)
x0 = (lb+ub)/2;
options = optimset('Display','iter','TolX',10^-8,'MaxIter',80);
[x,fval] = fmincon(@(x) fun_efficacy(x,Data,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0),x0,[],[], [], [], lb, ub, [], options);

%%
figure_setups; hold on
[~, xdata,ydata,yrun0,~,~] = fun_efficacy(x0,Data,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0);
[~, ~,~,yrun1,t,eff] = fun_efficacy(x,Data,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0);
scatter(xdata/365,ydata,100,'filled','DisplayName','data')
plot(xdata/365,yrun0,'*-','LineWidth',3,'DisplayName','initial - Penny estimate')
plot(t/365,eff,'-','LineWidth',3,'DisplayName','fitted')
axis([0 3 -0.2 1])
ylabel('Efficacy against clinical malaria')
xlabel('Time since vaccination (years)')
legend
%%
disp(['initial efficacy = ', num2str(x(1),2)])
disp(['ave immunity period = ', num2str(1/x(2)/365,2),' years'])

% 0.7137    0.0051  - large bounds