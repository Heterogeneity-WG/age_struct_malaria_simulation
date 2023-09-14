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

lowerB = 0.01;
upperB = 0.25;
base_value = P.muM;
param = lowerB:0.005:upperB;
%param = [0.01:0.01:sqrt(0.55)].^2; % used for betaM
% baseline: P.rD = 1/33.5, P.betaM = 0.25; P.betaD = 0.35; P.M_lower = 0.16; P.betaD_upper = 0.53;
tfinal = 10*365;
nt = length(0:dt:tfinal);
% create matrix to store results/output
EIR_timeseries = NaN(length(param),nt);
AH_matrix = NaN(length(param),na); % averages over one year, final year of simulation
DH_matrix = NaN(length(param),na); % averages over one year, final year of simulation

AH_time = NaN(length(param),nt); % store population weighted AH timeseries
DH_time = NaN(length(param),nt); % store population weighted DH timeseries

for ii = 1:length(param)
    progressbar(ii,length(param));
    P.muM = param(ii); % parameter to be varied
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
    %keyboard;
    AH_time(ii,:) = da*trapz(PH.*AH_tilde)./NH;
    DH_time(ii,:) = da*trapz(PH.*DH_tilde)./NH;;
end
toc;
%% Plotting of results
figure_setups; % average of A_H for each age group over one period of the limit cycle
imagesc(param,a/365,AH_matrix');
clim([0 1]);
colorbar;
title('$\int^{T}_{T-1} A_H(\alpha,t)dt$');
xlabel('$\mu_M$');
ylabel('Age ($\alpha$)');
ylim([0 30]);
%xticks([0 0.04 0.08 0.12]);
%xticks([min(param) median(param) max(param)]);
set(gca,'YDir','normal')
colormap jetwhite;
hold on;
xline([P.muM_lower base_value P.muM_upper],'--r','LineWidth',2.5); % indicator for the baseline

figure; % average of D_H for each age group over one period of the limit cycle
imagesc(param,a/365,DH_matrix');
clim([0 1]);
colorbar;
title('$\int^{T}_{T-1} D_H(\alpha,t)dt$');
xlabel('$\mu_M$');
ylabel('Age ($\alpha$)');
ylim([0 30]);
%xticks([0 0.04 0.08 0.12]);
set(gca,'YDir','normal')
colormap jetwhite;
hold on;
xline([P.muM_lower base_value P.muM_upper],'--r','LineWidth',3); % indicator for the baseline

figure; % EIR dynamics over one period of the limit cycle
years = (0:dt:tfinal)/365;
imagesc(param,years,EIR_timeseries');
clim([0 max(max(EIR_timeseries))]);
colorbar;
ylim([max(years)-1 max(years)]);
title('$EIR(t)$');
xlabel('$\beta_M$');
%xticks([0 0.04 0.08 0.12]);
ylabel('Time (years)');
%xticks([min(param) median(param) max(param)]);
set(gca,'YDir','normal')
colormap jetwhite;
hold on;
xline([P.betaM_lower base_value P.betaM_upper],'--r','LineWidth',3); % indicator for the baseline

%% More plotting, dynamics over time as the given parameter varies
M = length(param);
figure;
plot(years,AH_time(floor(M/4),:));
hold on;
plot(years,AH_time(floor(M/2),:));
plot(years,AH_time(floor(3*M/4),:));
plot(years,AH_time(floor(M),:));
xlabel('Time (years)');
legend(strcat('$\mu_M = $', num2str(param(floor(M/4)),1)),...
    strcat('$\mu_M = $', num2str(param(floor(M/2)),2)),...
    strcat('$\mu_M = $', num2str(param(floor(3*M/4)),2)),...
    strcat('$\mu_M = $', num2str(param(floor(M)),2)));
ylim([0 1]);
%title(['$\int_0^A P_H(\alpha,t)\tilde{A}_H(\alpha,t)d\alpha$']);
title(['$\bar{A}_H(t)$']);
xticks([0 2.5 5 7.5 10]);
yticks([0 0.25 0.5 0.75 1]);
grid on;

figure;
plot(years,DH_time(floor(M/4),:));
hold on;
plot(years,DH_time(floor(M/2),:));
plot(years,DH_time(floor(3*M/4),:));
plot(years,DH_time(floor(M),:));
xlabel('Time (years)');
legend(strcat('$\mu_M = $', num2str(param(floor(M/4)),1)),...
    strcat('$\mu_M = $', num2str(param(floor(M/2)),2)),...
    strcat('$\mu_M = $', num2str(param(floor(3*M/4)),2)),...
    strcat('$\mu_M = $', num2str(param(floor(M)),2)));
ylim([0 1]);
title(['$\bar{D}_H(t)$']);
yticks([0 0.25 0.5 0.75 1]);
xticks([0 2.5 5 7.5 10]);
grid on;

%% AH and DH proportion population averages versus parameter
figure;
title('Seasonal Averages');
plot(param, mean(AH_time(:,end-floor(365/dt):end),2)' );
hold on;
plot(param, mean(DH_time(:,end-floor(365/dt):end),2)' );
ylim([0 1]);
yticks([0 0.25 0.5 0.75 1]);
yyaxis right
plot(param, mean(EIR_timeseries(:,end-floor(365/dt):end),2)' );
legend('$\bar{A}_H$','$\bar{D}_H$','EIR');
xlabel('$\mu_M$');
grid on;








