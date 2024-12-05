clearvars
% close all
% clc

global P
% data primary only
% s =  load('Penny_ve_Siaya_fig2.mat','Penny_ve_Siaya_fig2'); Data = s.Penny_ve_Siaya_fig2; Data(:,1) = Data(:,1)/12; % Penny point data in figure 2
% data primary + booster
s =  load('Penny_paper/Penny_ve_Siaya_booster.mat','Penny_ve_Siaya_booster'); Data = s.Penny_ve_Siaya_booster;  % remove the first data point
Data(Data(:,1)>32/12,:) = []; % remove the data points beyond 32 months
Data(Data(:,1)<5/12,:) = []; % remove the data points before 6 months, time since vacc (completion of 3rd dose, 9month old, 1st dose = 3months, 6 months apart)
%% numerical config 
age_max = 100*365; % max ages in days
P.age_max = age_max;
dt = 5; % time/age step size in days, default = 5;
da = dt;
a = (0:da:age_max)';
na = length(a);

P.dt = dt; 
P.a = a;
P.na = na;
P.da = da;

Malaria_parameters_baseline;
Malaria_parameters_transform; 
P.v0s = 0; P.v0c = 0; P.vb0 = 0;
Malaria_parameters_transform_vac;

% initial condition
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');

%% fitting
% larger bounds
lb = [0.1, 1/(5*365)]; % P.etab; P.wb
ub = [1, 1/(0.1*365)]; 
% x0 = [0.72, 1/(0.53*365)];
x0 = [0.8, 1/(2.1*365)];
% x0 = [0.71, 1/(0.95*365)]; % primary estimates from Penny et al, Table 3 (half-life = 7.9 month)
% % x0 = [0.392, 1/(0.95*365)]; % booster estimates from Penny et al, Table 3 (half-life = 7.9 month, same as primary doses)

% lb = [0.1, 1/(5*365),0.1, 1/(5*365)]; % P.etab; P.wb0; P.wbs; P.wbr
% ub = [1, 1/(0.1*365),1, 1/(0.1*365)]; 
% x0 = [0.71, 1/(0.95*365),0.71, 1/(0.95*365)]; 

% x0 = (lb+ub)/2;

options = optimset('Display','iter','TolX',10^-8,'MaxIter',0);
[x,fval] = fmincon(@(x) fun_efficacy_booster(x,Data,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0),x0,[],[], [], [], lb, ub, [], options);
%%
figure_setups_B; hold on
[~, xdata,ydata,yrun0,t0,eff0] = fun_efficacy_booster(x0,Data,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0);
[~, ~,~,yrun1,t,eff] = fun_efficacy_booster(x,Data,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0);
scatter(xdata/365,ydata,100,'filled','DisplayName','Data')
% plot(xdata/365,yrun0,'*-','LineWidth',3,'DisplayName','initial - Penny estimate')
% plot(t0/365,eff0,'--','LineWidth',3,'DisplayName','initial')
plot(t/365,eff,'-','LineWidth',3,'DisplayName','Fitted')
axis([0 3 -0.2 1])
ylabel('Efficacy')
xlabel('Time since vaccination (years)')
% title('No booster')
% title('With booster')
legend
%%
disp(['initial efficacy = ', num2str(x(1),2)])
disp(['immunity period = ', num2str(1/x(2)/365,2),' years'])
% disp(['initial efficacy = ', num2str(x(3),2)])
% disp(['immunity period  = ', num2str(1/x(4)/365,2),' years'])
%%
% figure_setups;
% plot(P.a/30, 1./P.wb/365)
% title('immunity period')
% xlim([0 38])
% xlabel('month')
% ylabel('years')


