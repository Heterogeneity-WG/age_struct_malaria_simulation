close all
clc
format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

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

%% initial condition 'EE' - numerical EE
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');
NN_S = trapz(SH0)*P.da;
        
%% time evolution - initial run
tfinal = 6*365; t = (0:dt:tfinal)'; nt = length(t);
P.nt = nt;  P.t = t;
[~,SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(da,na,0,tfinal,...
    SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
PH = SH+EH+DH+AH+VH+UH;
PH_final = PH(:,end); % total human at age a, t = n
NH = trapz(PH,1)*da;
NM = SM+IM+EM;
vacc_sterile = trapz(P.v*(1-P.z).*SH,1)*P.da*365*P.NN/1000;
vacc_blood = trapz(P.v*P.z.*SH,1)*P.da*365*P.NN/1000;

%% DALY calculation
% [DALY,YLL,YLD] = DALY_cal(SH, EH, DH, AH, VH, UH, SM, EM, IM, Ctot); 
% figure_setups;
% plot(t/365,DALY,t/365,YLL,t/365,YLD)
% grid on
% xlabel('Year')
% legend('DALY','YLL (death)','YLD (disability)')
%% vac control
P.v0 = 0;
P.z = 1;
Malaria_parameters_transform_vac; 
tfinal_conti = 0*365; t2 = (tfinal:dt:tfinal+tfinal_conti)'; nt = length(t2);
P.nt = nt;  P.t = t2;
SH0 = SH(:,end); EH0 = EH(:,end); DH0 = DH(:,end); AH0 = AH(:,end); VH0 = VH(:,end); UH0 = UH(:,end); MH0 = MH(:,end);
SM0 = SM(end); EM0 = EM(end); IM0 = IM(end); 
Cac0 = Cac(:,end); Cm0 = Cm(:,end); Cv0 = Cv(:,end); Ctot0 = Ctot(:,end);
[~,SH2, EH2, DH2, AH2, VH2, UH2, SM2, EM2, IM2, Cm2, Cac2, Cv2, Ctot2, MH2] = age_structured_Malaria_vac(da,na,0,tfinal_conti,...
    SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
PH2 = SH2+EH2+DH2+AH2+VH2+UH2;
NH2 = trapz(PH2,1)*da;
NM2 = SM2+EM2+IM2;
vacc_sterile2 = trapz(P.v*(1-P.z).*SH2,1)*P.da;
vacc_blood2 = trapz(P.v*P.z.*SH2,1)*P.da;

t = [t;t2(2:end)];
SH = [SH,SH2(:,2:end)]; EH = [EH,EH2(:,2:end)]; DH = [DH, DH2(:,2:end)]; AH = [AH, AH2(:,2:end)]; VH = [VH, VH2(:,2:end)]; UH = [UH, UH2(:,2:end)]; MH = [MH, MH2(:,2:end)];
SM = [SM, SM2(2:end)]; EM = [EM, EM2(2:end)]; IM = [IM, IM2(2:end)]; NM = [NM, NM2(2:end)];
Cm = [Cm, Cm2(:,2:end)]; Cac = [Cac, Cac2(:,2:end)]; Cv = [Cv, Cv2(:,2:end)]; Ctot = [Ctot,Ctot2(:,2:end)];
PH = SH+EH+DH+AH+VH+UH;
PH_final = PH(:,end); % total human at age a, t = n
NH = [NH, NH2(2:end)]; 
vacc_sterile = [vacc_sterile, vacc_sterile2(2:end)];
vacc_blood = [vacc_blood, vacc_blood2(2:end)];
%% vaccine efficacy - reduction in incidence with or w/o vaccine (on population level)
% % initial condition = at the end of initial run
% %---- vaccinated group -----
% lgroup = 'DH'; tfinal_vacc = 10; tfinal_count = 3650;  % -> counting incidence for one year
% P.v0 = 100; % turn on vaccination
% Malaria_parameters_transform_vac; 
% [~,ind1] = min(abs(P.a-7*30));
% [~,ind2] = min(abs(P.a-19*30));
% [~,vacc,~,SHv, EHv, DHv, AHv, VHv, UHv, SMv, EMv, IMv, Cmv, Cacv, Cvv, Ctotv, MHv] = incidence_cal(da,na,tfinal_vacc,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0,ind1,ind2,lgroup);
% P.v0 = 0; % turn off vaccination
% Malaria_parameters_transform_vac; 
% Incidence_vacc = incidence_cal_time(da,na,tfinal_count,SHv, EHv, DHv, AHv, VHv, UHv, SMv, EMv, IMv, Cmv, Cacv, Cvv, Ctotv, MHv,ind1,ind2,lgroup);
% %---- control group -----
% P.v0 = 0;
% Malaria_parameters_transform_vac; 
% [~,~,~,SHc, EHc, DHc, AHc, VHc, UHc, SMc, EMc, IMc, Cmc, Cacc, Cvc, Ctotc, MHc] = incidence_cal(da,na,tfinal_vacc,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0,ind1,ind2,lgroup);
% Incidence_control = incidence_cal_time(da,na,tfinal_count,SHc, EHc, DHc, AHc, VHc, UHc, SMc, EMc, IMc, Cmc, Cacc, Cvc, Ctotc, MHc,ind1,ind2,lgroup);
% eff = (Incidence_control-Incid ence_vacc)./vacc
% figure_setups;
% plot(0:P.dt:tfinal_count,eff)
% figure_setups;
% plot(0:P.dt:tfinal_count,Incidence_control,0:P.dt:tfinal_count,Incidence_vacc)
% keyboard
%% output data to .mat file for analysis
% SH_EE = SH(:,end); EH_EE = EH(:,end); AH_EE = AH(:,end); DH_EE = DH(:,end); VH_EE = VH(:,end); UH_EE = UH(:,end); PH_EE = PH(:,end); v = P.v;
% Cm_EE = Cm(:,end); Cac_EE = Cac(:,end); Ctot_EE = Ctot(:,end);
% save(['Results/Vaccine/v0_',num2str(P.v0*100),'.mat'],'t','a','v','vacc_sterile','vacc_blood','SH_EE','EH_EE','AH_EE','DH_EE','VH_EE','UH_EE','PH_EE',...
%     'Cm_EE','Cac_EE','Ctot_EE');
%% EIR
[bH,~] = biting_rate(NH,NM);
EIR = bH.*IM./NM*365;
EIR_EE = EIR(end);

%% seasonlity plots
figure_setups;
yyaxis left
plot(t/365,EIR,'b-');
yyaxis right
plot(t/365, P.ss_S(t))
xlabel('Year')
legend('EIR','seasonal profile')
lamH = FOI_H(bH,IM,NM);
rho = sigmoid_prob(Ctot./PH, 'rho'); % prob. of severely infected, EH -> DH
psi = sigmoid_prob(Ctot./PH, 'psi'); % prob. AH -> DH
temp1 = psi.*lamH.*AH; % AH -> DH
temp2 = rho.*P.h.*EH;% EH -> DH, number of new cases - rate
[~,ind1] = min(abs(P.a-5*30));
[~,ind2] = min(abs(P.a-17*30));
cases_rate1 = trapz(temp1(ind1:ind2,:),1)*P.da;
cases_rate2 = trapz(temp2(ind1:ind2,:),1)*P.da;
%% calculate incidence per person per year
cases = cases_rate1+cases_rate2;
pop = trapz(PH(ind1:ind2,:),1)*P.da;
ind11 = length(t);
[~,ind22] = min(abs(t-(t(end)-365)));
cases_pp_py = trapz(cases(ind22:ind11))*P.dt/mean(pop(ind22:ind11));% should be around ~3 for Nanoro seasonlity study
%trapz(cases./pop)*P.dt/t(end)*365
%trapz(cases)*P.dt/(pop(end))/t(end)*365
%%
% figure_setups; hold on;
% plot(t/365,cumsum(cases_rate1)*dt);
% plot(t/365,cumsum(cases_rate2)*dt);
% title('Cumulative new (cohort) cases');
% legend('AH to DH','EH to DH');

%%
figure_setups;
bar(t'/365,(cases_rate1+cases_rate2).*30); % corresponding to Figure S1 in White et al. (2015)
xlabel('year')
ylabel('total cases (per 30days)')
title('for entire people cohort')
figure_setups;
plot(t/365,cases./pop*365)
ylim([0, 8])
xlabel('year')
ylabel('Incidence pp per year')
title(['Incidence pp per year = ', num2str(cases_pp_py)]);

%%
% figure_setups;
% plot(t/365,NM./NH,t/365,bH) % mosquito/human population ratio
% legend('$N_M/N_H$ ratio','biting')

%% vaccine #
% figure_setups; hold on
% subplot(1,2,1)
% yyaxis left
% plot(t/365,cumsum(vacc_blood)*dt/NN_S,t/365,cumsum(vacc_sterile)*dt/NN_S)
% legend('Blood-stage','Sterile')
% ylabel('fraction')
% axis([0 max(t)/365 0 1]);
% xlabel('years')
% title('cumulative vaccinated')
% yyaxis right
% plot(t/365,cumsum(vacc_blood)*dt,t/365,cumsum(vacc_sterile)*dt)
% ylabel('count')
% axis([0 max(t)/365 0 1*NN_S]);
% subplot(1,2,2)
% plot(P.a/30,P.v)
% xlim([0 25])
% xlabel('month')
% title('$\nu(\alpha)$ daily per-capita vacc. rate')
%% Population size versus time
% figure_setups;
% plot(t/365,trapz(SH,1)*da,'-','Color',colour_mat1); hold on;
% plot(t/365,trapz(EH,1)*da,'--','Color',colour_mat3);
% plot(t/365,trapz(AH,1)*da,'-.','Color',colour_mat2);
% plot(t/365,trapz(DH,1)*da,'-','Color',colour_mat7);
% plot(t/365,trapz(VH,1)*da,'-','Color',colour_mat6);
% plot(t/365,trapz(UH,1)*da,'-','Color',colour_mat4);
% plot(t/365,NH,'-.k')
% plot(t/365,trapz(MH,1)*da,'-.r'); % diagnostic
% legend('$S_H$','$E_H$','$A_H$', '$D_H$', '$V_H$','$U_H$','$N_H$','$M_H$ ($\mu_D$)', 'Location','e');
% title(['Population size vs time']); 
% grid on; grid minor
% axis([0 max(t)/365 0 max(NH)+0.1]);
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
%% Age profiles at tfinal
% figure_setups;
% plot(a/365,SH(:,end),'-','Color',colour_mat1); hold on;
% plot(a/365,EH(:,end),'--','Color',colour_mat3);
% plot(a/365,DH(:,end),'-.','Color',colour_mat2);
% plot(a/365,AH(:,end),':','Color',colour_mat7);
% plot(a/365,VH(:,end),':','Color',colour_mat6);
% plot(a/365,UH(:,end),':','Color',colour_mat4); 
% plot(a/365,PH_final,'-k');
% plot(a/365,MH(:,end),'-r');
% hold off;
% title(num2str(t(end)/365));
% xlim([0 15])
% legend('$S_H$','$E_H$','$D_H$','$A_H$','$V_H$','$U_H$','$P_H$','$M_H$ ($\mu_D$)');
% title(['Final Age Distribution']);
% xlabel('age (years)');
% grid on
% axis([0 age_max/365 0 max(PH_final)]);
%% Age proportions at tfinal
% figure_setups;
% plot(a/365,SH(:,end)./PH_final,'-','Color',colour_mat1); hold on;
% plot(a/365,EH(:,end)./PH_final,'--','Color',colour_mat6);
% plot(a/365,DH(:,end)./PH_final,'-.','Color',colour_mat2);
% plot(a/365,AH(:,end)./PH_final,':','Color',colour_mat3);
% plot(a/365,VH(:,end)./PH_final,':','Color',colour_mat5);
% plot(a/365,UH(:,end)./PH_final,':','Color',colour_mat4);
% plot(a/365,MH(:,end)./PH_final,'r-');
% plot(a/365,(AH(:,end)+DH(:,end))./PH_final,'r-.');
% plot(a/365,PH_final./PH_final,'-k');
% legend('SH/PH','EH/PH','DH/PH', 'AH/PH', 'VH/PH','UH/PH','$M_H$ ($\mu_D$)','(AH+DH)/PH');
% title(['Final Age Dist. Proportions']); 
% % title(['Final Age Dist. Proportions ~~ feedback =',num2str(immunity_feedback)]); 
% xlabel('age (years)');
% grid on
% axis([0 P.age_max/365 0 1.1]);
% xlim([0 15])
%% Age proportions in time - movie
% figure_setups;
% for iplot = 1:10:length(t)
%     plot(a/365,SH(:,iplot)./PH(:,iplot),'-','Color',colour_mat1); hold on;
%     plot(a/365,EH(:,iplot)./PH(:,iplot),'--','Color',colour_mat6);
%     plot(a/365,DH(:,iplot)./PH(:,iplot),'-.','Color',colour_mat2);
%     plot(a/365,AH(:,iplot)./PH(:,iplot),':','Color',colour_mat3);
%     plot(a/365,VH(:,iplot)./PH(:,iplot),':','Color',colour_mat5);
%     plot(a/365,UH(:,iplot)./PH(:,iplot),':','Color',colour_mat4);
%     plot(a/365,(AH(:,iplot)+DH(:,iplot))./PH(:,iplot),'r-.');
%     plot(a/365,PH(:,iplot)./PH(:,iplot),'-k');hold off
%     legend('SH/PH','EH/PH','DH/PH', 'AH/PH', 'VH/PH','UH/PH','(AH+DH)/PH');
%     title(['time = ', num2str(round(t(iplot)/365)), ' years']);
%     xlabel('age (years)');
%     grid on
%     axis([0 P.age_max/365 0 1.1]);
%     xlim([0 15])
%     pause(0.01)
% end
%% Immunity related figures
% figure_setups;
% subplot(2,2,1), plot(a/365,Ctot(:,floor(nt/4))./PH(:,floor(nt/4)));
% hold on;
% subplot(2,2,1), plot(a/365,Ctot(:,floor(nt/2))./PH(:,floor(nt/2)));
% subplot(2,2,1), plot(a/365,Ctot(:,floor(3*nt/4))./PH(:,floor(3*nt/4)));
% subplot(2,2,1), plot(a/365,Ctot(:,end)./PH(:,end),'-.');
% xlabel('age')
% legend(['t = ',num2str(t(end)/(4*365))],['t = ',num2str(t(end)/(2*365))],...
%     ['t = ',num2str(3*t(end)/(4*365))],['t = ',num2str(t(end)/365)],'Location','NorthWest');
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
% xlabel('time');
% title('$C_{total}(\alpha,t)/P_H(\alpha,t)$');
% 
% subplot(2,2,4), plot(a/365,Ctot(:,floor(nt/4)));
% hold on;
% subplot(2,2,4), plot(a/365,Ctot(:,floor(nt/2)));
% subplot(2,2,4), plot(a/365,Ctot(:,floor(3*nt/4)));
% subplot(2,2,4), plot(a/365,Ctot(:,end),'-.');
% xlabel('age')
% legend(['t = ',num2str(t(end)/(4*365))],['t = ',num2str(t(end)/(2*365))],...
%     ['t = ',num2str(3*t(end)/(4*365))],['t = ',num2str(t(end)/365)]);
% title('$C_{total}(t)$');
%% 
% figure_setups;
% subplot(2,2,1), plot(a/365,Cv(:,floor(nt/4))./PH(:,floor(nt/4)));
% hold on;
% subplot(2,2,1), plot(a/365,Cv(:,floor(nt/2))./PH(:,floor(nt/2)));
% subplot(2,2,1), plot(a/365,Cv(:,floor(3*nt/4))./PH(:,floor(3*nt/4)));
% subplot(2,2,1), plot(a/365,Cv(:,end)./PH(:,end),'-.');
% xlabel('age')
% legend(['t = ',num2str(t(end)/(4*365))],['t = ',num2str(t(end)/(2*365))],...
%     ['t = ',num2str(3*t(end)/(4*365))],['t = ',num2str(t(end)/365)],'Location','NorthWest');
% title('$C_{v}(t)/P_H(t)$');
% grid on
% subplot(2,2,2), plot(t/365,(trapz(Cv,1)*da)./NH);
% title('$\int C_{v}(\alpha,t)d\alpha / N_H(t)$');
% xlabel('time');
% grid on
% 
% subplot(2,2,3), imagesc(t/365,a/365,Cv./PH);
% set(gca,'YDir','normal');
% colorbar;
% ylabel('age');
% xlabel('time');
% title('$C_{v}(\alpha,t)/P_H(\alpha,t)$');
% 
% subplot(2,2,4), plot(a/365,Cv(:,floor(nt/4)));
% hold on;
% subplot(2,2,4), plot(a/365,Cv(:,floor(nt/2)));
% subplot(2,2,4), plot(a/365,Cv(:,floor(3*nt/4)));
% subplot(2,2,4), plot(a/365,Cv(:,end),'-.');
% xlabel('age')
% legend(['t = ',num2str(t(end)/(4*365))],['t = ',num2str(t(end)/(2*365))],...
%     ['t = ',num2str(3*t(end)/(4*365))],['t = ',num2str(t(end)/365)]);
% title('$C_{v}(t)$');
%% Immunity breakdown
% figure_setups;
% plot(a/365,Cac(:,end),'-.r');
% hold on;
% plot(a/365,Cm(:,end),'-.b');
% plot(a/365,Cv(:,end),'-.c');
% plot(a/365,Ctot(:,end),'-.k');
% xlabel('age (years)')
% legend('Acquired','Maternal','Vaccine-derived','Total','Location','SouthEast');
% title(['Immun dist.']);
% axis([0 age_max/365 0 max(Ctot(:,end))*1.1]);
% grid on
%%
% figure_setups;
% plot(a/365,Cac(:,end)./PH_final,'-.');
% hold on;
% plot(a/365,Cm(:,end)./PH_final,'--');
% plot(a/365,Cv(:,end)./PH_final,'--');
% plot(a/365,Ctot(:,end)./PH_final,'-');
% xlabel('age (years)')
% ylabel('immunity level')
% legend('Acquired (pp)','Maternal (pp)','Vaccine-derived (pp)','Total (pp)','Location','SouthEast');
% % title(['Per-person Immun dist.~~ feedback =',num2str(immunity_feedback)]);
% title('Per-person Immun distribution');
% axis([0 age_max/365 0 max(Ctot(:,end)./PH_final)*1.1]);
% xlim([0 10])
% ylim([0 7])
% grid on
%% plot sigmoids
% figure_setups; hold on;
% plot(a/365,sigmoid_prob(Ctot(:,end)./PH_final, 'rho'),'-');
% % plot(a/365,sigmoid_prob(Ctot(:,end)./PH_final, 'psi'));
% % plot(a/365,sigmoid_prob(Ctot(:,end)./PH_final, 'phi'));
% grid on
% legend('rho (suscept.)');
% axis([0 age_max/365 0 1]);
% title(['EIR = ',num2str(EIR_EE)])
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

