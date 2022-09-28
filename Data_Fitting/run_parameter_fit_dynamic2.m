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
for irun = 1:nrun
    disp(['run #',num2str(irun)])
    Malaria_parameters_baseline;
    [x,fval] = fmincon(@fun_Filipe_dynamic2,x0,[],[], [], [], lb, ub, [], options);
    Xmat(irun,1) = fval;
    Xmat(irun,2:end) = x;
%     save('Data_Fitting/Filipe_paper/Xmat_fine.mat','Xmat','irun','x_lhs')
end
keyboard
%%
load('Data_Fitting/Filipe_paper/Xmat.mat','Xmat','irun','x_lhs')
% for irun = 1:60
%     x = Xmat(irun,2:end);
%     Malaria_parameters_baseline;
%     P.phis2 = x(1);
%     P.phir2 = x(2);
%     P.rhos2 = x(3);
%     P.rhor2 = x(4);
%     P.psis2 = x(5);
%     P.psir2 = x(6);
%     Malaria_parameters_transform;
%     [SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0] = age_structured_Malaria_IC('init');
%     [SH, EH, DH, AH, SM, EM, IM, ~, ~, Ctot] = age_structured_Malaria(P.da,P.na,P.tfinal,SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0);
%     PH = SH(:,end)+EH(:,end)+DH(:,end)+AH(:,end);
%     NH = trapz(PH)*P.da;
%     Ctot_pp = Ctot(:,end)./PH;
%     figure_setups; hold on
%     cc = linspace(0,max(Ctot_pp),100);
%     phi_curve = sigmoid_prob(cc, 'phi');
%     rho_curve = sigmoid_prob(cc, 'rho');
%     psi_curve = sigmoid_prob(cc, 'psi');
%     plot(cc,phi_curve,'-')
%     plot(cc,rho_curve,'--')
%     plot(cc,psi_curve,'-.')
%     % population average sigmoids
%     phi_ave = trapz(sigmoid_prob(Ctot_pp, 'phi').*PH/NH)*P.da;
%     rho_ave = trapz(sigmoid_prob(Ctot_pp, 'rho').*PH/NH)*P.da;
%     psi_ave = trapz(sigmoid_prob(Ctot_pp, 'psi').*PH/NH)*P.da;
%     legend('$\phi(\tilde{C}_{H})$','$\rho(\tilde{C}_{H})$','$\psi(\tilde{C}_{H})$','Location','e')
%     axis([0 max(cc) 0 1])
%     xlabel('$\tilde{C}_{H}$')
%     ylabel('Probability')
%     %title('Calibrated linking functions')
%     [phi_ave rho_ave psi_ave]
%     fun_Filipe_dynamic2(x);
% end

% keyboard
% ----> update Malaria_parameters_baseline.m with the fitted results <-----

%% plot sigmoids with the populational average (in legend)
tfinal = 10*365; age_max = 100*365; P.age_max = age_max;
dt = 20; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t; P.tfinal = tfinal;
x = Xmat(55,2:end);
% [~,ind1] = min(abs(P.a-0.3*365)); 
% [~,ind0] = min(abs(P.a-3*365)); 
% [~,ind2] = min(abs(P.a-20*365)); 
% ind_a = [round(linspace(ind1,ind0,15)'); round(linspace(ind0+1,ind2,nsamp-15)')];
% x = [1.174354913102552   3.865058149466797   1.551219543369452   0.787630781652653   4.279593577031902   3.994787277041905];
% fine EIR mesh
% x = [3.266020807401758   3.248827198899857   1.594701633752857   0.837465044857498   4.664942886501461   4.300057715509009];
% [~,ind1] = min(abs(P.a-0.3*365)); 
% [~,ind0] = min(abs(P.a-10*365)); 
% [~,ind2] = min(abs(P.a-20*365)); 
% ind_a = [round(linspace(ind1,ind0,20)'); round(linspace(ind0+1,ind2,nsamp-20)')];
% x = [3.291597293958823   4.243807990790567   1.646220399576347   0.774696249293831   2.119792106093347   4.427487320234270];
% [~,ind1] = min(abs(P.a-0.3*365)); 
% [~,ind0] = min(abs(P.a-20*365)); 
% [~,ind2] = min(abs(P.a-50*365)); 
% ind_a = [round(linspace(ind1,ind0,20)'); round(linspace(ind0+1,ind2,nsamp-20)')];
% x = [1.457776066263742   4.097149381426517   1.314712272979270   0.796127082004893   2.874959267261557   4.173818846106313];
Malaria_parameters_baseline;
P.phis2 = x(1);
P.phir2 = x(2); 
P.rhos2 = x(3);
P.rhor2 = x(4); 
P.psis2 = x(5);
P.psir2 = x(6);
Malaria_parameters_transform;
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('init');
[SH, EH, DH, AH, ~, ~, SM, EM, IM, ~, ~, ~, Ctot, ~] = age_structured_Malaria_vac(P.da,P.na,P.tfinal,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
PH = SH(:,end)+EH(:,end)+DH(:,end)+AH(:,end);
NH = trapz(PH)*P.da;
% EIR = fit_EIR(SH,EH,DH,AH,SM, EM, IM);   
% keyboard
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
P.a = a; P.na = na; P.nt = nt; P.dt = dt;    P.da = da; P.t = t; P.tfinal = tfinal;

Malaria_parameters_baseline;
Malaria_parameters_transform;
% Malaria_parameters_transform_vac;
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
    [SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('init');
    [SH, EH, DH, AH, ~, ~, SM, EM, IM, ~, ~, ~, Ctot, ~] = age_structured_Malaria_vac(P.da,P.na,P.tfinal,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
    EIR = fit_EIR(SH,EH,DH,AH,SM, EM, IM);   
    PH = SH+EH+DH+AH;    
    yy(1,jj) = EIR(end); % aEIR
    xxs(:,jj) = xx;
    yys(:,jj) = EIR(end)*ones(size(yys(:,jj)));
    zz(:,jj) = Ctot(:,end)./PH(:,end); % final Ctot at EE    
end
figure_setups;
plot(var_list,yy)
%% sorting
% [yy,ind] = sort(yy);
% zz = zz(:,ind);
%% plotting heatmap (age, EIR, immunity) 
figure_setups; imagesc(xx,yy,zz')
% plot(var_list,y,'-o');
figure_setups; grid off; [XX,YY] = meshgrid(xx,yy); surf(XX,YY,zz'); shading interp; view(2)
figure_setups; grid off; scatter3(xxs(:),yys(:),zz(:))
axis([0 20 0 190])
xlabel('Age (years)')
ylabel('aEIR')
title('Immunity level per person');
set(gca,'YDir','normal');
colormap jet
caxis([0 12])
colorbar
zz_rho = sigmoid_prob(zz, 'rho');
%% plotting heatmap (age, EIR, rho) -> to compare with Filipe's curves
figure_setups; hold on; grid off
zz_rho = sigmoid_prob(zz, 'rho');
% zz_psi = sigmoid_prob(zz, 'psi');
% imagesc(xx,yy,zz_rho')
[XX,YY] = meshgrid(xx,yy);
surf(xx,yy,zz_rho')
shading interp
axis([0 20 0 190])
view(2)
xlabel('Age (years)')
ylabel('aEIR')
% title(['Immunity levels, feedback = ',num2str(immunity_feedback)]);
title('Susceptibility ($\rho$)');
set(gca,'YDir','normal');
colormap jet
colorbar('Ticks',0:0.2:1);
caxis([0 max(max(zz_rho))])
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