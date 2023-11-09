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
options = optimset('Display','iter','TolX',10^-6,'MaxIter',30);
% [phi_s phi_r rho_s psi_r rho_s psi_r]
lb = [0, 0.1, 0, 0.1, 0, 0.1];
ub = [5, 5, 5, 5, 5, 5];
x0 = (lb+ub)/2;
tic
[x,fval] = fmincon(@fun_Filipe_dynamic,x0,[],[], [], [], lb, ub, [], options);
toc
keyboard
% ----> update Malaria_parameters_baseline.m with the fitted results <-----

%% plot sigmoids with the populational average (in legend)
tfinal = 10*365; age_max = 100*365; P.age_max = age_max;
dt = 20; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t; P.tfinal = tfinal;

% uniform EIR sampling
x = [2.567957971786876   2.487540758554113   3.649596968324358   1.395449806257184   2.332526365071812   2.150211932758257];

% refine EIR sampling linspace(0,0.1,50)
% x = [4.054028031572193   3.112804946380192   1.795855001626887   1.489513759827821   3.576364022306257   2.176976575527686];

% refine EIR sampling linspace(0,0.05,50)
% x = [2.447889755984070   2.544920188222090   1.841524526659086   2.359471176773747   2.390735795064254   2.419704159102402];

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
title('Calibrated linking functions')
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
    var_list = [0.01:0.01:1].^2;
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
%% calculate two slices
beta_slice1 = 0.02; % EIR ~ 32
P.betaM = beta_slice1;
Malaria_parameters_transform;
[SH, EH, DH, AH, ~, ~, SM, EM, IM, Cm1, Cac1, ~, Ctot1, ~] = age_structured_Malaria_IC_vac('EE_reset');
PH = SH+EH+DH+AH; PH1 = PH;
NM = SM+EM+IM;
[bH,~] = biting_rate(PH,NM);
EIR = bH.*IM./NM*365; % EIR matrix
EIR_tot1 = trapz(EIR.*PH)/trapz(PH); % EIR sum over age, at final time
beta_slice2 = 0.25; % EIR ~ 79
P.betaM = beta_slice2;
Malaria_parameters_transform;
[SH, EH, DH, AH, ~, ~, SM, EM, IM, Cm2, Cac2, ~, Ctot2, ~] = age_structured_Malaria_IC_vac('EE_reset');
PH = SH+EH+DH+AH; PH2 = PH;
NM = SM+EM+IM;
[bH,~] = biting_rate(PH,NM);
EIR = bH.*IM./NM*365; % EIR matrix
EIR_tot2 = trapz(EIR.*PH)/trapz(PH); % EIR sum over age, at final time
%% plotting heatmap (age, EIR, immunity) 
figure_setups; hold on; grid off
set(gcf,'Position',[353   307   552   446])
imagesc(xx,yy,zz')
xlim([0 20])
ylim([0 max(yy)])
xlabel('Age (years)')
ylabel('aEIR')
title('Immunity level per person');
set(gca,'YDir','normal');
colormap jet
colorbar('Ticks',0:2:15);
% plot([0 20],[EIR_tot1 EIR_tot1],'k--')
% plot([0 20],[EIR_tot2 EIR_tot2],'k--')
%% plotting slices of immunity heatmap
figure_setups;
set(gcf,'Position',[353   307   552   446])
plot(P.a/365,Cac2(:,end)./PH2(:,end),'-.',P.a/365,Cm2(:,end)./PH2(:,end),'--', P.a/365,Ctot2(:,end)./PH2(:,end))
xlim([0 10])
ylim([0 10])
legend off
title('Per-person immunity dist.');
xlabel('Age (years)')
ylabel('Immunity level')
figure_setups;
set(gcf,'Position',[353   307   552   446])
plot(P.a/365,Cac1(:,end)./PH1(:,end),'-.',P.a/365,Cm1(:,end)./PH1(:,end),'--', P.a/365,Ctot1(:,end)./PH1(:,end))
xlim([0 10])
ylim([0 10])
title('Per-person immunity dist.');
legend('$\widetilde{C}_e$ (Exposure-acquired immunity)','$\widetilde{C}_{m}$ (Maternal immunity)','$\widetilde{C}_{H}$ (Total immunity)','Location','North');
xlabel('Age (years)')
ylabel('Immunity level')
%% plotting heatmap (age, EIR, rho)
figure_setups; hold on; grid off
set(gcf,'Position',[353   307   552   446])
zz_rho = sigmoid_prob(zz, 'rho');
imagesc(xx,yy,zz_rho')
xlabel('Age (years)')
ylabel('aEIR')
title('Susceptibility ($\rho$)');
set(gca,'YDir','normal');
colormap jet
colorbar('Ticks',0:0.2:1);
xlim([0 20])
ylim([0 max(yy)])
caxis([0 max(max(zz_rho))])
%% slices of heatmap to compare with Filipe's curves
zz_rho = sigmoid_prob(zz, 'rho');
EIR_list = [1 10 20 50 100];
figure_setups; hold on
xx_data = 0.3:0.2:50; 
xx_fit = xx;
for iEIR = 1:5
    zz_data = F(xx_data,EIR_list(iEIR)*ones(size(xx_data)))';
    plot(xx_data,zz_data,'-','DisplayName',['EIR=',num2str(EIR_list(iEIR))])
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

