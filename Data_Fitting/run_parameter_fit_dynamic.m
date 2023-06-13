close all;
clearvars;
clc;
format long;
global P
global F

% numerical config
tfinal = 10*365; age_max = 100*365; P.age_max = age_max;
dt = 20; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t; P.tfinal = tfinal;

%% load data - interpolated function
s =  load('Filipe_paper/F_Filipe.mat','F');
F = s.F;

%% optimization - fit to data 
Malaria_parameters_baseline;
Malaria_parameters_transform;
Malaria_parameters_transform_vac;
options = optimset('Display','iter','TolX',10^-5,'MaxIter',30);
% [phi_s phi_r rho_s psi_r rho_s psi_r]
lb = [0, 0.1, 0, 0.1, 0, 0.1];
ub = [5, 4, 5, 4, 5, 4];
x0 = (lb+ub)/2;
% x0 = [1.748400494250446   4.089465852051163   2.781182708408349   3.349185468908294   1.267935962166972   2.767371953595199]; % Tfinal = 10 years IC = [4, 5, 4, 5];
% x0 = [0.978511191108162   0.101270457374745   1.957064739448506   3.354032211122635   0.001254524505761   1.119700445524332];

[x,fval] = fmincon(@fun_Filipe_dynamic,x0,[],[], [], [], lb, ub, [], options);
keyboard
% ----> update Malaria_parameters_baseline.m with the fitted results <-----

%% plot sigmoids with the populational average (in legend)
tfinal = 10*365; age_max = 100*365; P.age_max = age_max;
dt = 20; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t; P.tfinal = tfinal;

x = [2.190590606590128   2.078437497848872   2.921061093971156   1.218264166958853   0.614949157082559   2.462790900030959];
%[0.978511191108162   0.101270457374745   1.957064739448506   3.354032211122635   0.001254524505761   1.119700445524332];
Malaria_parameters_baseline;
P.phis2 = x(1);
P.phir2 = x(2); 
P.rhos2 = x(3);
P.rhor2 = x(4); 
P.psis2 = x(5);
P.psir2 = x(6);
Malaria_parameters_transform;
Malaria_parameters_transform_vac;
[SH, EH, DH, AH, ~, ~, SM, EM, IM, ~, ~, ~, Ctot, ~] = age_structured_Malaria_IC_vac('EE_reset');
PH = SH+EH+DH+AH;
NH = trapz(PH)*P.da;
Ctot_pp = Ctot./PH;
figure_setups; hold on
cc = linspace(0,max(Ctot_pp),100);
phi_curve = sigmoid_prob(cc, 'phi');
rho_curve = sigmoid_prob(cc, 'rho');
psi_curve = sigmoid_prob(cc, 'psi');
plot(cc,phi_curve,'-')
plot(cc,rho_curve,'--')
plot(cc,psi_curve,'-.')
% population average sigmoids
phi_ave = trapz(sigmoid_prob(Ctot_pp, 'phi').*PH/NH)*P.da;
rho_ave = trapz(sigmoid_prob(Ctot_pp, 'rho').*PH/NH)*P.da;
psi_ave = trapz(sigmoid_prob(Ctot_pp, 'psi').*PH/NH)*P.da;
legend('$\phi(\tilde{C}_{H})$','$\rho(\tilde{C}_{H})$','$\psi(\tilde{C}_{H})$','Location','e')
axis([0 max(cc) 0 1])
xlabel('$\tilde{C}_{H}$')
ylabel('Probability')
% title('Calibrated linking functions')
[phi_ave rho_ave psi_ave]


%% plotting heatmap (age, EIR, immunity level)
tfinal = 20*365; age_max = 100*365; P.age_max = age_max;
dt = 20; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t; P.tfinal = tfinal;

Malaria_parameters_baseline;
Malaria_parameters_transform;
Malaria_parameters_transform_vac;
P.phis2 = x(1);
P.phir2 = x(2); 
P.rhos2 = x(3);
P.rhor2 = x(4); 
P.psis2 = x(5);
P.psir2 = x(6);
Malaria_parameters_transform;
immunity_feedback = 1;
if immunity_feedback == 0 
    % population average sigmoids: f0 = f1 = average
    P.phif0 = 0.886990461669485; % value at zero
    P.phif1 = 0.886990461669485; % value at L (function saturates to this value)
    
    P.rhof0 = 0.112790581535007; % value at zero
    P.rhof1 = 0.112790581535007; % value at L (function saturates to this value)  
    
    P.psif0 = 0.056984636045019; % value at zero
    P.psif1 = 0.056984636045019; % value at L (function saturates to this value)   
    var_list = [0.01:0.05:1.0].^2; 
else
    var_list = [0.01:0.1:1].^2;
end
xx = P.a/365;
yy = zeros(1,length(var_list));
zz = zeros(na,length(var_list));
for jj = 1:length(var_list)
    P.betaM = var_list(jj);
    Malaria_parameters_transform;
    [SH, EH, DH, AH, ~, ~, SM, EM, IM, ~, ~, ~, Ctot, ~] = age_structured_Malaria_IC_vac('EE_reset');
    PH = SH+EH+DH+AH;
    NM = SM+EM+IM;
    [bH,~] = biting_rate(PH,NM);
    EIR = bH.*IM./NM*365; % EIR matrix
    EIR_tot = trapz(EIR.*PH)/trapz(PH); % EIR sum over age, at final time
    yy(1,jj) = EIR_tot; % aEIR
    zz(:,jj) = Ctot./PH; % final Ctot at EE    
end
%%
% figure_setups;
% plot(var_list,y,'-o');
% keyboard
figure_setups;
imagesc(xx,yy,zz')
xlim([0 20])
xlabel('Age (years)')
ylabel('aEIR')
% title(['Immunity levels, feedback = ',num2str(immunity_feedback)]);
title('Immunity level per person');
set(gca,'YDir','normal');
colormap jet
colorbar('Ticks',0:2:15);
zz = sigmoid_prob(zz, 'rho');
%% plotting heatmap (age, EIR, rho) -> to compare with Filipe's curves
% figure_setups; hold on
% zz = sigmoid_prob(zz, 'rho');
% imagesc(xx,yy,zz')
% axis([0 20 0 120])
% xlabel('Age (years)')
% ylabel('aEIR')
% % title(['Immunity levels, feedback = ',num2str(immunity_feedback)]);
% title('Susceptibility ($\rho$)');
% set(gca,'YDir','normal');
% colormap jet
% colorbar('Ticks',0:0.1:1);
% EIR_list = [1 10 20 50 100];
% for iEIR = 1:length(EIR_list)
%     plot(xx,EIR_list(iEIR)*ones(size(xx)),'m--')
% end
%% slices of heatmap to compare with Filipe's curves
EIR_list = [1 10 20 50 100];
figure_setups; hold on
xx_data = 0.3:0.2:50; 
xx_fit = xx;
for iEIR = 4
    zz_data = F(xx_data,EIR_list(iEIR)*ones(size(xx_data)))';
    plot(xx_data,zz_data,'--')
    [~,ind] = min(abs(EIR_list(iEIR)-yy));
    zz_fit = zz(:,ind);
    plot(xx,zz_fit,'-.')
end
legend('Location','e')
xlabel('age (years)')
ylabel('susceptibility')
axis([0 50 0 1.05 ])

