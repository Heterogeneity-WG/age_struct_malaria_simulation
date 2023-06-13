% calibration of sigmoids 
% with multiple initial sampling of x0 for fmincon -> nrun
% and fixed beta_list values
% Note: set tfinal = 10 years for IC = EE_reset to cut down CPU time

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
load('Filipe_paper/F_Filipe.mat','F');
Malaria_parameters_baseline;
Malaria_parameters_transform;
Malaria_parameters_transform_vac;
np = 6; % # of sigmoid parameters to be fitted
nrun = 1; % # of LHC samples to generate

%% optimization - fit to data 
lb = [0, 0.1, 0, 0.1, 0, 0.1];
ub = [5, 5, 5, 5, 5, 5];

% x_lhsraw = lhsdesign(nrun,np); 
% x_lhs = NaN(size(x_lhsraw));
% for k = 1:np
%     x_lhs(:,k)=(x_lhsraw(:,k).*(ub(k)-lb(k)))+lb(k);
% end

% [phi_s phi_r rho_s rho_s psi_r psi_r phif0 phif1 rhof0 rhof1 psif0 psif1]
x0 = (lb+ub)/2;

options = optimset('Display','iter','TolX',10^-5,'MaxIter',30);
Xmat = NaN(nrun,np+1);
tic
for irun = 1:nrun
    disp(['run #',num2str(irun)])
    Malaria_parameters_baseline;
    % x0 = x_lhs(irun,:);
    [x,fval] = fmincon(@fun_Filipe_dynamic2,x0,[],[], [], [], lb, ub, [], options);
    Xmat(irun,1) = fval;
    Xmat(irun,2:end) = x;
    % save('Data_Fitting/Filipe_paper/Xmat_new.mat','Xmat','nrun','x_lhs')
    save('Data_Fitting/Filipe_paper/Xmat_new.mat','Xmat','nrun','x0')
end
toc
%% plot sigmoids and calculate the populational average 
% load('Data_Fitting/Filipe_paper/Xmat_new.mat','Xmat','nrun','x_lhs')
load('Data_Fitting/Filipe_paper/Xmat_new.mat','Xmat','nrun','x0')
for irun = 1:nrun
    x = Xmat(irun,2:end);
    Malaria_parameters_baseline;
    P.phis2 = x(1);
    P.phir2 = x(2);
    P.rhos2 = x(3);
    P.rhor2 = x(4);
    P.psis2 = x(5);
    P.psir2 = x(6);
    Malaria_parameters_transform;
    Malaria_parameters_transform_vac;
    [SH, EH, DH, AH, ~, ~, ~, ~, ~, ~, ~, ~, Ctot, ~] = age_structured_Malaria_IC_vac('EE_reset');
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
    title('Calibrated linking functions')
    [phi_ave rho_ave psi_ave]
    % fun_Filipe_dynamic2(x);
end

% keyboard
% ----> update Malaria_parameters_baseline.m with the fitted results <-----

%% generate data for heatmap (age, EIR, immunity level)
tfinal = 10*365; age_max = 100*365; P.age_max = age_max;
dt = 20; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt;    P.da = da; P.t = t; P.tfinal = tfinal;

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
    var_list = [0.01:0.01:1].^2;
end
xx = P.a/365;
yy = zeros(1,length(var_list));
yys = zeros(na,length(var_list));
xxs = zeros(na,length(var_list));
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
%% plotting heatmap (age, EIR, immunity) 
figure_setups;
imagesc(xx,yy,zz')
xlim([0 20])
ylim([0 max(yy)])
xlabel('Age (years)')
ylabel('aEIR')
title('Immunity level per person');
set(gca,'YDir','normal');
colormap jet
colorbar('Ticks',0:2:15);
%% plotting heatmap (age, EIR, rho)
figure_setups; hold on; grid off
zz_rho = sigmoid_prob(zz, 'rho');
imagesc(xx,yy,zz_rho')
xlabel('Age (years)')
ylabel('aEIR')
title('Susceptibility ($\rho$)');
set(gca,'YDir','normal');
colormap jet
colorbar('Ticks',0:0.2:1);
xlim([0 20])
caxis([0 max(max(zz_rho))])
%% slices of heatmap to compare with Filipe's curves
EIR_list = [1 10 20 50 100];
figure_setups; hold on
xx_data = 0.3:0.2:50; 
xx_fit = xx;
for iEIR = 1:5
    zz_data = F(xx_data,EIR_list(iEIR)*ones(size(xx_data)))';
    plot(xx_data,zz_data,'-')
end
for iEIR = 1:5
    [~,ind] = min(abs(EIR_list(iEIR)-yy));
    zz_fit = zz_rho(:,ind);    
    plot(xx,zz_fit,'-.')       
end

legend('Location','e')
xlabel('age (years)')
ylabel('susceptibility')
axis([0 50 0 1.05 ])
