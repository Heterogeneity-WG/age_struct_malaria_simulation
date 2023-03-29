clearvars
% close all
clc

format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2
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
% using fitted value
P.etas = 0.659553156743644;
P.w = 0.004566152172269; 
P.v0s = 0; P.v0c = 0;
Malaria_parameters_transform;
Malaria_parameters_transform_vac;
t0_list=[1,5,8,11]; % the start time of the vaccination programv (1=Jan, 2=Feb, etc..)
%% initial condition 'EE' - numerical EE
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');
for it = 1:length(t0_list)
    t0_vacc = t0_list(it)*30; % time to start vaccination
    [~,SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(P.da,P.na, 0, t0_vacc,...
        SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
    SH0 = SH(:,end); EH0 = EH(:,end); DH0 = DH(:,end); AH0 = AH(:,end); VH0 = VH(:,end); UH0 = UH(:,end); MH0 = MH(:,end);
    SM0 = SM(end); EM0 = EM(end); IM0 = IM(end);
    Cac0 = Cac(:,end); Cm0 = Cm(:,end); Cv0 = Cv(:,end); Ctot0 = Ctot(:,end);
    %% simulation using three-group model - vaccine on
    tfinal_vacc = 3*30; total_vacc = 2.4*10^5; % vacc for three months
    P.v0s = total_vacc/tfinal_vacc; P.v0c = P.v0s; % define constant vaccination rate
    Malaria_parameters_transform_vac;
    temp0 = zeros(size(SH0));
    [t, SHr, EHr, DHr, AHr, Cmr, Cacr, Ctotr, SHv, EHv, DHv, AHv, VHv, UHv, Cmv, Cacv, Ctotv, SHc, EHc, DHc, AHc, Cmc, Cacc, Ctotc, SM, EM, IM] = ...
        age_structured_Malaria_eff(da, na, t0_vacc, t0_vacc+tfinal_vacc, SH0, EH0, DH0, AH0, Cm0, Cac0, Ctot0, ...
        temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, SM0, EM0, IM0);
    %% simulation using three-group model - vaccine off
    tfinal_conti = 3*365; total_vacc = 0;
    P.v0s = total_vacc/tfinal_vacc; P.v0c = P.v0s; % define constant vaccination rate
    Malaria_parameters_transform_vac;
    [t2, SHr2, EHr2, DHr2, AHr2, Cmr2, Cacr2, Ctotr2, SHv2, EHv2, DHv2, AHv2, VHv2, UHv2, Cmv2, Cacv2, Ctotv2, SHc2, EHc2, DHc2, AHc2, Cmc2, Cacc2, Ctotc2, SM2, EM2, IM2] = ...
        age_structured_Malaria_eff(da, na, t(end), t(end)+tfinal_conti, SHr(:,end), EHr(:,end), DHr(:,end), AHr(:,end), Cmr(:,end), Cacr(:,end), Ctotr(:,end), ...
        SHv(:,end), EHv(:,end), DHv(:,end), AHv(:,end), VHv(:,end), UHv(:,end), Cmv(:,end), Cacv(:,end), Ctotv(:,end), ...
        SHc(:,end), EHc(:,end), DHc(:,end), AHc(:,end), Cmc(:,end), Cacc(:,end), Ctotc(:,end), SM(end), EM(end), IM(end));
    % combine results
    SHr = [SHr,SHr2(:,2:end)]; EHr = [EHr,EHr2(:,2:end)]; DHr = [DHr, DHr2(:,2:end)]; AHr = [AHr, AHr2(:,2:end)]; Cmr = [Cmr, Cmr2(:,2:end)]; Cacr = [Cacr, Cacr2(:,2:end)]; Ctotr = [Ctotr,Ctotr2(:,2:end)];
    SHc = [SHc,SHc2(:,2:end)]; EHc = [EHc,EHc2(:,2:end)]; DHc = [DHc, DHc2(:,2:end)]; AHc = [AHc, AHc2(:,2:end)]; Cmc = [Cmc, Cmc2(:,2:end)]; Cacc = [Cacc, Cacc2(:,2:end)]; Ctotc = [Ctotc,Ctotc2(:,2:end)];
    SHv = [SHv,SHv2(:,2:end)]; EHv = [EHv,EHv2(:,2:end)]; DHv = [DHv, DHv2(:,2:end)]; AHv = [AHv, AHv2(:,2:end)]; VHv = [VHv, VHv2(:,2:end)]; UHv = [UHv, UHv2(:,2:end)]; Cmv = [Cmv, Cmv2(:,2:end)]; Cacv = [Cacv, Cacv2(:,2:end)]; Ctotv = [Ctotv,Ctotv2(:,2:end)];
    SM = [SM, SM2(2:end)]; EM = [EM, EM2(2:end)]; IM = [IM, IM2(2:end)];
    t = [t;t2(2:end)];
    %% calculate incidences
    PHr = SHr+EHr+DHr+AHr;
    PHc = SHc+EHc+DHc+AHc;
    PHv = SHv+EHv+DHv+AHv+VHv+UHv;
    PH = PHr+PHc+PHv;
    NH = trapz(PH,1)*da;
    NM = SM+EM+IM;
    [bH,~] = biting_rate(NH,NM);
    lamH = FOI_H(bH,IM,NM);
    
    rhov = sigmoid_prob(Ctotv./PHv, 'rho'); % prob. of severely infected, EH -> DH
    psiv = sigmoid_prob(Ctotv./PHv, 'psi'); % prob. AH -> DH
    rhov(PHv==0)=0; psiv(PHv==0)=0;
    temp1 = rhov.*P.h.*EHv+psiv.*lamH.*AHv; % incidence of DH
    Incidence_vacc = trapz(temp1,1)*P.da;
    
    rhoc = sigmoid_prob(Ctotc./PHc, 'rho'); % prob. of severely infected, EH -> DH
    psic = sigmoid_prob(Ctotc./PHc, 'psi'); % prob. AH -> DH
    rhoc(PHc==0)=0; psic(PHc==0)=0;
    temp2 = rhoc.*P.h.*EHc+psic.*lamH.*AHc; % incidence of DH
    Incidence_control = trapz(temp2,1)*P.da;
    
    %% calculate efficacy - using instantaneous incidence
    eff = (Incidence_control'-Incidence_vacc')./Incidence_control';
    figure_setups;
    plot(t/365,eff)
    xlabel('Year')
    ylabel('Efficacy')
    % title(['$\eta_s=$',num2str(P.etas), ',   $w=$',num2str(1/P.w/365),'years'])]
    title(['vacc starts at month=',num2str(t(2)/365*12,2)])
    % axis([0 max(t/365) -0.2 1])
    axis([t(2)/365 max(t/365) -0.2 1])
    save(['Results/season_vacc_',num2str(t0_vacc/30,2),'.mat'],'t','eff')
end
%% comparison vacc starting time
figure_setups; hold on
for it = 1:length(t0_list)
    load(['Results/season_vacc_',num2str(t0_list(it)),'.mat'],'t','eff')
    plot(t(2:end)/365,eff(2:end),'DisplayName',['t0=',num2str(t0_list(it))])   
%     plot((t(2:end)-t0_list(it)*30)/365,eff(2:end),'DisplayName',['t0=',num2str(t0_list(it))])    
end
xlabel('Time since vacc')
ylabel('Efficacy')
legend;
ylim([-0.2 1])
%% calculate efficacy - using aggregated incidence (within three month period) 
% time_period = 30*3;
% Incidence_vacc_agg = NaN(size(t));
% Incidence_control_agg = NaN(size(t));
% for it = 1:length(t)
%     period_end = t(it);
%     period_beg = t(it)-time_period;
%     if period_beg<0; continue; end
%     [~,ind_beg] = min(abs(t-period_beg)); ind_beg = ind_beg+1;
%     [~,ind_end] = min(abs(t-period_end)); 
%     Incidence_vacc_agg(it) = trapz(Incidence_vacc(ind_beg:ind_end))*P.dt;
%     Incidence_control_agg(it) = trapz(Incidence_control(ind_beg:ind_end))*P.dt;
% end
% eff_agg = (Incidence_control_agg'-Incidence_vacc_agg')./Incidence_control_agg';
% s =  load('Penny_ve_Siaya_fig2.mat','Penny_ve_Siaya_fig2'); Data = s.Penny_ve_Siaya_fig2; % Penny point data in figure 2
% figure_setups; hold on
% grid on
% scatter(Data(:,1)/12, Data(:,2),'filled')
% plot(t/365,eff,t/365,eff_agg)
% xlabel('Year')
% ylabel('Efficacy')
% title(['$\eta_s=$',num2str(P.etas), ',   $w=$',num2str(1/P.w/365),'years'])
% axis([t(2)/365 max(t/365) -0.2 1])
% legend('Penny data fig2','efficacy (instant cases)','efficacy (three-monthly cases)')

%% plots for debugging
% SH = SHr+SHv+SHc;
% EH = EHr+EHv+EHc;
% AH = AHr+AHv+AHc;
% DH = DHr+DHv+DHc;
% VH = VHv; UH = UHv;
% PH_final = PH(:,end);
%% Population proportions versus time
% figure_setups;
% plot(t/365,trapz(SH,1)*da./NH,'-','Color',colour_mat1); hold on;
% plot(t/365,trapz(EH,1)*da./NH,'--','Color',colour_mat3);
% plot(t/365,trapz(AH,1)*da./NH,'-.','Color',colour_mat2);
% plot(t/365,trapz(DH,1)*da./NH,'-','Color',colour_mat7);
% plot(t/365,trapz(VH,1)*da./NH,'-','Color',colour_mat6);
% plot(t/365,trapz(UH,1)*da./NH,'-','Color',colour_mat4);
% plot(t/365,(trapz(SH,1)+trapz(EH,1)+trapz(AH,1)+trapz(DH,1)+trapz(VH,1)+trapz(UH,1))*da./NH,'-.k');
% legend('SH-age','EH-age','AH-age', 'DH-age', 'VH-age','UH-age','$N_H$');
% title('Population proportions vs time');
% xlabel('years');
% grid on
% axis([0 max(t)/365 0 1.1]);
%% Age proportions at tfinal
% figure_setups;
% plot(a/365,SH(:,end)./PH_final,'-','Color',colour_mat1); hold on;
% plot(a/365,EH(:,end)./PH_final,'--','Color',colour_mat6);
% plot(a/365,DH(:,end)./PH_final,'-.','Color',colour_mat2);
% plot(a/365,AH(:,end)./PH_final,':','Color',colour_mat3);
% plot(a/365,VH(:,end)./PH_final,':','Color',colour_mat5);
% plot(a/365,UH(:,end)./PH_final,':','Color',colour_mat4);
% plot(a/365,(AH(:,end)+DH(:,end))./PH_final,'r-.');
% plot(a/365,PH_final./PH_final,'-k');
% legend('SH/PH','EH/PH','DH/PH', 'AH/PH', 'VH/PH','UH/PH','(AH+DH)/PH');
% title(['Final Age Dist. Proportions']); 
% xlabel('age (years)');
% grid on
% axis([0 P.age_max/365 0 1.1]);
% xlim([0 15])
%%
% figure_setups;
% plot(a/365,SHr(:,end)./PHr(:,end),'-','Color',colour_mat1); hold on;
% plot(a/365,EHr(:,end)./PHr(:,end),'--','Color',colour_mat6);
% plot(a/365,DHr(:,end)./PHr(:,end),'-.','Color',colour_mat2);
% plot(a/365,AHr(:,end)./PHr(:,end),':','Color',colour_mat3);
% plot(a/365,(AHr(:,end)+DHr(:,end))./PHr(:,end),'r-.');
% plot(a/365,PHr(:,end)./PHr(:,end),'-k');
% legend('SH/PH','EH/PH','DH/PH', 'AH/PH', '(AH+DH)/PH');
% title(['Final Age Dist. Proportions']); 
% % title(['Final Age Dist. Proportions ~~ feedback =',num2str(immunity_feedback)]); 
% xlabel('age (years)');
% grid on
% % axis([0 P.age_max/365 0 1.1]);
% xlim([0 15])
%% Immunity breakdown
% figure_setups;
% subplot(1,2,1); hold on
% plot(a/365,Cacv(:,end)./PHv(:,end),'-.r');
% hold on;
% plot(a/365,Cmv(:,end)./PHv(:,end),'-.b');
% xlabel('age (years)')
% legend('Acquired','Maternal','Location','SouthEast');
% title('vaccinated group');
% axis([0 5 0 11]);
% grid on
% subplot(1,2,2); hold on
% plot(a/365,Cacc(:,end)./PHc(:,end),'-.r');
% hold on;
% plot(a/365,Cmc(:,end)./PHc(:,end),'-.b');
% xlabel('age (years)')
% legend('Acquired','Maternal','Location','SouthEast');
% title('control group');
% axis([0 5 0 11]);
% grid on
% figure_setups;
% plot(a/365,Cacr(:,end)./PHr(:,end),'-.r');
% hold on;
% plot(a/365,Cmr(:,end)./PHr(:,end),'-.b');
% xlabel('age (years)')
% legend('Acquired','Maternal','Location','SouthEast');
% title('rest group');
% axis([0 5 0 11]);
% grid on
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