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
Malaria_parameters_baseline;
Malaria_parameters_transform;
Malaria_parameters_transform_vac;
np = 6; % # of sigmoid parameters to be fitted
nrun = 60; % # of LHC samples to generate
%% optimization - fit to data 
% lb = [0, 0.1, 0, 0.1, 0, 0.1];
% ub = [5, 5, 5, 5, 5, 5];
% x_lhsraw = lhsdesign(nrun,np); 
% x_lhs = NaN(size(x_lhsraw));
% for k = 1:np
%     x_lhs(:,k)=(x_lhsraw(:,k).*(ub(k)-lb(k)))+lb(k);
% end
% 
% % [phi_s phi_r rho_s rho_s psi_r psi_r]
% x0 = [2, 1, 2, 1, 2, 1];
% 
% options = optimset('Display','iter','TolX',10^-3,'MaxIter',30);
% Xmat = NaN(nrun,np+1);
% for irun = 1:nrun
%     disp(['run #',num2str(irun)])
%     x0 = x_lhs(irun,:);
%     Malaria_parameters_baseline;
%     [x,fval] = fmincon(@fun_Filipe_dynamic,x0,[],[], [], [], lb, ub, [], options);
%     Xmat(irun,1) = fval;
%     Xmat(irun,2:end) = x;
% %     save('Data_Fitting/Filipe_paper/Xmat3.mat','Xmat','irun','x_lhs')
% end
% keyboard
%%
load('Data_Fitting/Filipe_paper/Xmat2.mat','Xmat','irun','x_lhs')
for irun = 55
    x = Xmat(irun,2:end);
    Malaria_parameters_baseline;
    P.phis2 = x(1);
    P.phir2 = x(2);
    P.rhos2 = x(3);
    P.rhor2 = x(4);
    P.psis2 = x(5);
    P.psir2 = x(6);
    Malaria_parameters_transform;
    [SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0] = age_structured_Malaria_IC('init');
    [SH, EH, DH, AH, SM, EM, IM, ~, ~, Ctot] = age_structured_Malaria(P.da,P.na,P.tfinal,SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0);
    PH = SH(:,end)+EH(:,end)+DH(:,end)+AH(:,end);
    NH = trapz(PH)*P.da;
    Ctot_pp = Ctot(:,end)./PH;
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
    %title('Calibrated linking functions')
    [phi_ave rho_ave psi_ave]
    fun_Filipe_dynamic(x);
end

keyboard
% ----> update Malaria_parameters_baseline.m with the fitted results <-----

%% plot sigmoids with the populational average (in legend)
tfinal = 10*365; age_max = 100*365; P.age_max = age_max;
dt = 20; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t; P.tfinal = tfinal;

% x = Xmat(21,2:end);
% x = Xmat(55,2:end);
% x = Xmat(39,2:end);
% x = [4.988291588990569   2.493769945649229   3.682573000812414   2.520124591290541   4.143564601288833   1.363574026818044]; % Xmat2 - 55
% x = [3.158023763585565   2.472630554404047   3.905495296989682   2.488459672461482   0.972329999517317   2.126640707617171]; % Xmat3 - 39
Malaria_parameters_baseline;
P.phis2 = x(1);
P.phir2 = x(2); 
P.rhos2 = x(3);
P.rhor2 = x(4); 
P.psis2 = x(5);
P.psir2 = x(6);
Malaria_parameters_transform;
[SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0] = age_structured_Malaria_IC('init');
[SH, EH, DH, AH, SM, EM, IM, ~, ~, Ctot] = age_structured_Malaria(P.da,P.na,P.tfinal,SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0);
PH = SH(:,end)+EH(:,end)+DH(:,end)+AH(:,end);
NH = trapz(PH)*P.da;
Ctot_pp = Ctot(:,end)./PH;
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
%title('Calibrated linking functions')
[phi_ave rho_ave psi_ave]


%% generate data for heatmap (age, EIR, immunity level)
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
    var_list = [0.01:0.01:1].^2;
end
xx = P.a/365;
yy = zeros(1,length(var_list));
zz = zeros(na,length(var_list));
for jj = 1:length(var_list)
    P.betaM = var_list(jj);
    Malaria_parameters_transform;
    [SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0] = age_structured_Malaria_IC('init');
    [SH, EH, DH, AH, SM, EM, IM, ~, ~, Ctot] = age_structured_Malaria(P.da,P.na,P.tfinal,SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0);
    EIR = fit_EIR(SH,EH,DH,AH,SM, EM, IM);   
    PH = SH+EH+DH+AH;    
    yy(1,jj) = EIR(end); % aEIR
    zz(:,jj) = Ctot(:,end)./PH(:,end); % final Ctot at EE    
end
%% plotting heatmap (age, EIR, immunity) 
% figure_setups;
% plot(var_list,y,'-o');
% keyboard
figure_setups;
imagesc(xx,yy,zz')
% axis([0 20 0 120])
xlim([0 20])
xlabel('Age (years)')
ylabel('aEIR')
% title(['Immunity levels, feedback = ',num2str(immunity_feedback)]);
title('Immunity level per person');
set(gca,'YDir','normal');
colormap jet
colorbar('Ticks',0:2:15);
zz_rho = sigmoid_prob(zz, 'rho');
%% plotting heatmap (age, EIR, rho) -> to compare with Filipe's curves
figure_setups; hold on
zz_rho = sigmoid_prob(zz, 'rho');
zz_psi = sigmoid_prob(zz, 'psi');
imagesc(xx,yy,zz_rho')
% axis([0 20 0 120])
xlim([0 20])
xlabel('Age (years)')
ylabel('aEIR')
% title(['Immunity levels, feedback = ',num2str(immunity_feedback)]);
title('Susceptibility ($\rho$)');
set(gca,'YDir','normal');
colormap jet
colorbar('Ticks',0:0.1:1);
% EIR_list = [1 10 20 50 100];
% for iEIR = 1:length(EIR_list)
%     plot(xx,EIR_list(iEIR)*ones(size(xx)),'m--')
% end
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
    ind
end
legend('Location','e')
xlabel('age (years)')
ylabel('susceptibility')
axis([0 50 0 1.05 ])


function EIR = fit_EIR(SH,EH,DH,AH,SM,EM,IM)
global P
NH = trapz(SH+EH+DH+AH)*P.da;
NM = SM+EM+IM;
[bH,~] = biting_rate(NH,NM);
IM_frac = IM./NM;
EIR = bH.*IM_frac*365; % annual EIR
end