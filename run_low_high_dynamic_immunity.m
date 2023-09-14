% close all
clc
format long
global P
tic
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
% model parameters
Malaria_parameters_baseline;
Malaria_parameters_transform;
Malaria_parameters_transform_vac;

lowerB = P.rA_lower*0.1;
upperB = P.rA_upper*2;
base_value = P.rA;
param = lowerB:0.001:upperB;
%param = [0.01:0.01:sqrt(0.55)].^2; % used for betaM
% baseline: P.rD = 1/33.5, P.betaM = 0.25; P.betaD = 0.35; P.M_lower = 0.16; P.betaD_upper = 0.53;
tfinal = 10*365;
nt = length(0:dt:tfinal);
% create matrix to store results/output
EIR_timeseries = NaN(length(param),nt);
AH_matrix = NaN(length(param),na); % averages over one year, final year of simulation
DH_matrix = NaN(length(param),na);

for ii = 1:length(param)
    progressbar(ii,length(param));
    P.rA = param(ii); % parameter to be varied
    %% calculate initial condition 'EE' - numerical EE
    [SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');
    NN_S = trapz(SH0)*P.da;
    %% time evolution
    P.v0 = 0;
    P.z = 0; % z=0 sterile, z=1 blood-stage
    Malaria_parameters_transform_vac; % resetting vaccination distribution
    [t,SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(da,na,0,tfinal,...
        SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
    PH = SH+EH+DH+AH+VH+UH;
    PH_final = PH(:,end); % total human at age a, t = n
    NH = trapz(PH,1)*da;
    NM = SM+IM+EM;
    vacc_sterile = trapz(P.v*(1-P.z).*SH,1)*P.da*365*P.NN/1000;
    vacc_blood = trapz(P.v*P.z.*SH,1)*P.da*365*P.NN/1000;
    [bH,~] = biting_rate(PH,NM); % NB bH varies by age and time
    EIR = bH.*IM./NM*365; % EIR matrix
    EIR_tot = trapz(EIR.*PH,1)*P.da./NH;

    AH_tilde = AH./PH;
    DH_tilde = DH./PH;
    % store output
    AH_matrix(ii,:) = mean(AH_tilde(:,end-floor(365/dt):end),2);
    DH_matrix(ii,:) = mean(DH_tilde(:,end-floor(365/dt):end),2);
    EIR_timeseries(ii,:) = EIR_tot;
end
toc;
%% Plotting of results
figure; 
imagesc(param,a/365,AH_matrix');
clim([0 1]);
colorbar;
title('$\int^{T}_{T-1} A_H(\alpha,t)dt$');
xlabel('$r_A$');
ylabel('Age ($\alpha$)');
ylim([0 30]);
%xticks([0 0.04 0.08 0.12]);
%xticks([min(param) median(param) max(param)]);
set(gca,'YDir','normal')
colormap jetwhite;
hold on;
xline([P.rA_lower base_value P.rA_upper],'--r','LineWidth',2.5); % indicator for the baseline

figure; 
imagesc(param,a/365,DH_matrix');
clim([0 1]);
colorbar;
title('$\int^{T}_{T-1} D_H(\alpha,t)dt$');
xlabel('$r_A$');
ylabel('Age ($\alpha$)');
ylim([0 30]);
%xticks([0 0.04 0.08 0.12]);
set(gca,'YDir','normal')
colormap jetwhite;
hold on;
xline([P.rA_lower base_value P.rA_upper],'--r','LineWidth',3); % indicator for the baseline

figure;
years = (0:dt:tfinal)/365;
imagesc(param,years,EIR_timeseries');
clim([0 max(max(EIR_timeseries))]);
colorbar;
ylim([max(years)-1 max(years)]);
title('$EIR(t)$');
xlabel('$r_A$');
%xticks([0 0.04 0.08 0.12]);
ylabel('Time (years)');
%xticks([min(param) median(param) max(param)]);
set(gca,'YDir','normal')
colormap jetwhite;
hold on;
xline([P.rA_lower base_value P.rA_upper],'--r','LineWidth',3); % indicator for the baseline

%% Population proportions & EIR versus time
% figure_setups;
% figure(1);
% subplot(1,2,1), plot(t/365,EIR(end,:));
% hold on;
% plot(t/365,EIR(floor(end/10),:));
% plot(t/365,EIR(floor(end/100),:));
% legend('EIR age 100', 'EIR age 10','EIR age 1');
% xlabel('Time (years)');
% title('aEIR dynamics');
% figure_setups;
% plot(t/365,(PH(:,end)')*EIR/NH(end));
% title('Average EIR (pop. weighted)');
% figure(1);
% subplot(1,2,2), plot(t/365,trapz(SH,1)*da./NH,'-','Color',colour_mat1); hold on;
% plot(t/365,trapz(EH,1)*da./NH,'-','Color',colour_mat3);
% plot(t/365,trapz(AH,1)*da./NH,'-','Color',colour_mat2);
% plot(t/365,trapz(DH,1)*da./NH,'-','Color',colour_mat7);
% plot(t/365,trapz(VH,1)*da./NH,'-','Color',colour_mat6);
% plot(t/365,trapz(UH,1)*da./NH,'-','Color',colour_mat4);
% plot(t/365,(trapz(SH,1)+trapz(EH,1)+trapz(AH,1)+trapz(DH,1)+trapz(VH,1)+trapz(UH,1))*da./NH,'-.k');
% legend('$\tilde{S}_H$','$\tilde{E}_H$','$\tilde{A}_H$', '$\tilde{D}_H$', '$\tilde{V}_H$','$\tilde{U}_H$','$N_H$');
% title('Disease dynamics');
% xlabel('Time (years)');
% grid on
% axis([0 max(t)/365 0 1.1]);

%% Disease class dynamics
% figure_setups;
% 
% tempNorm = SH./PH;
% subplot(2,2,1), plot(t/365,tempNorm(end,:));
% hold on;
% plot(t/365,tempNorm(floor(end/10),:));
% plot(t/365,tempNorm(floor(end/100),:));
% ylim([0 1]);
% title('$\tilde{S}_{H}(\alpha,t)$');
% grid on;
% 
% tempNorm = EH./PH;
% subplot(2,2,2), plot(t/365,tempNorm(end,:));
% hold on;
% plot(t/365,tempNorm(floor(end/10),:));
% plot(t/365,tempNorm(floor(end/100),:));
% title('$\tilde{E}_{H}(\alpha,t)$');
% legend('Age 100', 'Age 10','Age 1');
% ylim([0 1]);
% grid on;
% 
% tempNorm = AH./PH;
% subplot(2,2,3), plot(t/365,tempNorm(end,:));
% hold on;
% plot(t/365,tempNorm(floor(end/10),:));
% plot(t/365,tempNorm(floor(end/100),:));
% title('$\tilde{A}_{H}(\alpha,t)$');
% ylim([0 1]);
% grid on;
% xlabel('Time (years)');
% 
% tempNorm = DH./PH;
% subplot(2,2,4), plot(t/365,tempNorm(end,:));
% hold on;
% plot(t/365,tempNorm(floor(end/10),:));
% plot(t/365,tempNorm(floor(end/100),:));
% title('$\tilde{D}_{H}(\alpha,t)$');
% ylim([0 1]);
% grid on;
% xlabel('Time (years)');

%% Immunity dynamics
% figure_setups;
% nt = length(t);
% subplot(2,2,1), plot(a/365,Ctot(:,floor(nt/40))./PH(:,floor(nt/40)));
% hold on;
% subplot(2,2,1), plot(a/365,Ctot(:,floor(3*nt/40))./PH(:,floor(3*nt/40)));
% subplot(2,2,1), plot(a/365,Ctot(:,floor(nt/20))./PH(:,floor(nt/20)));
% subplot(2,2,1), plot(a/365,Ctot(:,floor(nt/10))./PH(:,floor(nt/10)));
% %subplot(2,2,1), plot(a/365,Ctot(:,end)./PH(:,end),'-.');
% xlabel('age')
% legend(['t = ',num2str(t(end)/(40*365))],['t = ',num2str(t(end)/(20*365))],...
%     ['t = ',num2str(3*t(end)/(40*365))],['t = ',num2str(t(end)/(10*365))],'Location','SouthEast');
% title('$C_{total}(t)/P_H(t)$');
% grid on
% subplot(2,2,2), plot(t/365,(trapz(Ctot,1)*da)./NH);
% title('$\int C_{total}(\alpha,t)d\alpha / N_H(t)$');
% xlabel('time');
% grid on
% 
% subplot(2,2,3), imagesc(t/365,a/365,Ctot./PH);
% set(gca,'YDir','normal');
% colorbar;
% ylabel('age');
% xlabel('Time (years)');
% title('$C_{total}(\alpha,t)/P_H(\alpha,t)$');
% 
% Ctot_norm = Ctot./PH;
% subplot(2,2,4), plot(t/365,Ctot_norm(floor(na/100),:));
% hold on;
% subplot(2,2,4), plot(t/365,Ctot_norm(floor(na/20),:));
% subplot(2,2,4), plot(t/365,Ctot_norm(floor(na/10),:));
% subplot(2,2,4), plot(t/365,Ctot_norm(na,:),'-.');
% xlabel('Time(years)')
% legend('age 1','age 5','age 10', 'age 100'); % assume age_max=100
% title('$\tilde{C}_{total}(t,\cdot)$');

%% DALY calculation
% [DALY,YLL,YLD] = DALY_cal(SH, EH, DH, AH, VH, UH, SM, EM, IM, Ctot);
% figure(4);
% subplot(2,2,1), plot(t/365,DALY,t/365,YLL,t/365,YLD);
% grid on
% xlabel('Year')
% legend('DALY','YLL (death)','YLD (disability)')

%% Mosquito population size
% figure_setups;
% plot(t/365,SM,'b-'); hold on;
% plot(t/365,EM,'-','Color',colour_r1);
% plot(t/365,IM,'r-.');
% plot(t/365,SM+EM+IM,'-.')
% legend('$S_M$','$E_M$','$I_M$','$N_M$');
% title('mosquito population size by stages')
% grid on
% xlim([0 tfinal/365])

%% Mosquito population proportion
% figure_setups;
% plot(t/365,SM./NM,'b-'); hold on;
% plot(t/365,EM./NM,'-','Color',colour_r1);
% plot(t/365,IM./NM,'r-.');
% plot(t/365,(SM+EM+IM)./NM,'-.')
% legend('$S_M$','$E_M$','$I_M$','$N_M$');
% title('mosquito population prop by stages')
% grid on
% xlim([0 tfinal/365])

