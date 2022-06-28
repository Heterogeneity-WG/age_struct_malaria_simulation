clear all
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
dt = 20; % time/age step size in days, default = 5;
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

immunity_feedback = 1;
if immunity_feedback == 0
    P.phi_f_0 = 0.570320665853183; % value at zero
    P.phi_f_1 = 0.570320665853183; % value at L (function saturates to this value)
    
    P.rho_f_0 = 0.088575583518581; % value at zero
    P.rho_f_1 = 0.088575583518581; % value at L (function saturates to this value)  
    
    P.psi_f_0 = 0.409302219871934; % value at zero
    P.psi_f_1 = 0.409302219871934; % value at L (function saturates to this value)   
end

% [SH, EH, DH, AH, VH, SM, EM, IM, Cm, Cac, Cv, Ctot] = age_structured_Malaria_IC_vac('EE');
% PH = SH+EH+DH+AH+VH;
% PH_final = PH(:,end); % total human at age a, t = n
% NH = trapz(PH,1)*da;
%% initial condition 'init' 'EE'
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('init');
%% time evolution - initial run
tfinal = 500*365; t = (0:dt:tfinal)'; nt = length(t);
P.nt = nt;  P.t = t;
[SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(da,na,tfinal,...
    SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
PH = SH+EH+DH+AH+VH+UH;
PH_final = PH(:,end); % total human at age a, t = n
NH = trapz(PH,1)*da;
NM = SM+IM+EM;
vacc_sterile = trapz(P.v*(1-P.z).*SH,1)*P.da*365*P.NN/1000;
vacc_blood = trapz(P.v*P.z.*SH,1)*P.da*365*P.NN/1000;
%% time evolution - continuous run
% non-vac control
% cont_level = 0.5;
% P.betaA = (1-cont_level)*P.betaA;
% P.betaD = (1-cont_level)*P.betaD;
% P.betaM = (1-cont_level)*P.betaM;
% Malaria_parameters_transform;
%% vac control
% P.v0 = 0.8;
% Malaria_parameters_transform_vac; 
% tfinal_conti = 10*365; t2 = (tfinal:dt:tfinal+tfinal_conti)'; nt = length(t2);
% P.nt = nt;  P.t = t2;
% SH0 = SH(:,end); EH0 = EH(:,end); DH0 = DH(:,end); AH0 = AH(:,end); VH0 = VH(:,end); UH0 = UH(:,end); MH0 = MH(:,end);
% SM0 = SM(end); EM0 = EM(end); IM0 = IM(end); 
% Cac0 = Cac(:,end); Cm0 = Cm(:,end); Cv0 = Cv(:,end); Ctot0 = Ctot(:,end);
% [SH2, EH2, DH2, AH2, VH2, UH2, SM2, EM2, IM2, Cm2, Cac2, Cv2, Ctot2, MH2] = age_structured_Malaria_vac(da,na,tfinal_conti,...
%     SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
% PH2 = SH2+EH2+DH2+AH2+VH2+UH2;
% NH2 = trapz(PH2,1)*da;
% vacc_sterile2 = trapz(P.v*(1-P.z).*SH2,1)*P.da*365*P.NN/1000;
% vacc_blood2 = trapz(P.v*P.z.*SH2,1)*P.da*365*P.NN/1000;
% 
% t = [t;t2];
% SH = [SH,SH2]; EH = [EH,EH2]; DH = [DH, DH2]; AH = [AH, AH2]; VH = [VH, VH2]; UH = [UH, UH2]; MH = [MH, MH2];
% SM = [SM, SM2]; EM = [EM, EM2]; IM = [IM, IM2];
% Cm = [Cm, Cm2]; Cac = [Cac, Cac2]; Cv = [Cv, Cv2]; Ctot = [Ctot,Ctot2];
% PH = SH+EH+DH+AH+VH+UH;
% PH_final = PH(:,end); % total human at age a, t = n
% NH = [NH, NH2];
% vacc_sterile = [vacc_sterile, vacc_sterile2];
% vacc_blood = [vacc_blood, vacc_blood2];
%% output data to .mat file for analysis
% SH_EE = SH(:,end); EH_EE = EH(:,end); AH_EE = AH(:,end); DH_EE = DH(:,end); VH_EE = VH(:,end); UH_EE = UH(:,end); PH_EE = PH(:,end); v = P.v;
% Cm_EE = Cm(:,end); Cac_EE = Cac(:,end); Ctot_EE = Ctot(:,end);
% save(['Results/Vaccine/v0_',num2str(P.v0*100),'.mat'],'t','a','v','vacc_sterile','vacc_blood','SH_EE','EH_EE','AH_EE','DH_EE','VH_EE','UH_EE','PH_EE',...
%     'Cm_EE','Cac_EE','Ctot_EE');
%% EIR
% [bh,bm] = biting_rate(NH,NM);
% EIR = bh.*IM./NM*365;
% EIR_EE = EIR(end)
% tic
% R0 = R0_cal();
% keyboard
% toc

% figure_setups;
% plot(t,EIR,'b-'); hold on;
%% vaccine #
% figure_setups; hold on
% grid on
% plot(t/365,vacc_blood,t/365,vacc_sterile)
% legend('Blood-stage','Sterile')
% title('vaccination counts')
% xlim([0 max(t)/365]);
% xlabel('years')
%% Population size versus time
figure_setups;
plot(t/365,trapz(SH,1)*da,'-','Color',colour_mat1); hold on;
plot(t/365,trapz(EH,1)*da,'--','Color',colour_mat3);
plot(t/365,trapz(AH,1)*da,'-.','Color',colour_mat2);
plot(t/365,trapz(DH,1)*da,'-','Color',colour_mat7);
plot(t/365,trapz(VH,1)*da,'-','Color',colour_mat6);
plot(t/365,trapz(UH,1)*da,'-','Color',colour_mat4);
plot(t/365,NH,'-.k')
plot(t/365,trapz(MH,1)*da,'-.r'); % diagnostic
legend('$S_H$','$E_H$','$A_H$', '$D_H$', '$V_H$','$U_H$','$N_H$','$M_H$ ($\mu_D$)', 'Location','e');
title(['Population size vs time', '~~feedback = ',num2str(immunity_feedback)]); 
grid on; grid minor
axis([0 max(t)/365 0 max(NH)+0.1]);
%% Population proportions versus time
figure_setups;
plot(t/365,trapz(SH,1)*da./NH,'-','Color',colour_mat1); hold on;
plot(t/365,trapz(EH,1)*da./NH,'--','Color',colour_mat3);
plot(t/365,trapz(AH,1)*da./NH,'-.','Color',colour_mat2);
plot(t/365,trapz(DH,1)*da./NH,'-','Color',colour_mat7);
plot(t/365,trapz(VH,1)*da./NH,'-','Color',colour_mat6);
plot(t/365,trapz(UH,1)*da./NH,'-','Color',colour_mat4);
plot(t/365,(trapz(SH,1)+trapz(EH,1)+trapz(AH,1)+trapz(DH,1)+trapz(VH,1)+trapz(UH,1))*da./NH,'-.k');
legend('SH-age','EH-age','AH-age', 'DH-age', 'VH-age','UH-age','$N_H$');
title('Population proportions vs time');
xlabel('years');
grid on
axis([0 max(t)/365 0 1.1]);
%% Age profiles at tfinal
% figure_setups;
% plot(a/365,SH2(:,end),'-','Color',colour_mat1); hold on;
% plot(a/365,EH2(:,end),'--','Color',colour_mat3);
% plot(a/365,DH2(:,end),'-.','Color',colour_mat2);
% plot(a/365,AH2(:,end),':','Color',colour_mat7);
% plot(a/365,VH2(:,end),':','Color',colour_mat6);
% plot(a/365,UH2(:,end),':','Color',colour_mat4); hold off;
% title(num2str(t(end)/365));
% xlim([0 15])
% plot(a/365,PH_final,'-k');
% plot(a/365,MH(:,end),'-r');
% legend('$D_H$','$A_H$','$V_H$','$U_H$');
% legend('$S_H$','$E_H$','$A_H$','$D_H$','$V_H$','$U_H$','$P_H$','$M_H$ ($\mu_D$)');
% title(['Final Age Dist.,~~ feedback =',num2str(immunity_feedback)]);
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
% plot(a/365,(AH(:,end)+DH(:,end))./PH_final,'r-.');
% plot(a/365,PH_final./PH_final,'-k');
% legend('SH/PH','EH/PH','DH/PH', 'AH/PH', 'VH/PH','UH/PH','(AH+DH)/PH');
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
% legend(['t = ',num2str(tfinal/(4*365))],['t = ',num2str(tfinal/(2*365))],...
%     ['t = ',num2str(3*tfinal/(4*365))],['t = ',num2str(tfinal/365)],'Location','NorthWest');
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
% legend(['t = ',num2str(tfinal/(4*365))],['t = ',num2str(tfinal/(2*365))],...
%     ['t = ',num2str(3*tfinal/(4*365))],['t = ',num2str(tfinal/365)]);
% title('$C_{total}(t)$');
%% Immunity breakdown
% figure_setups;
% plot(a/365,Cac(:,end),'-.r');
% hold on;
% plot(a/365,Cm(:,end),'-.b');
% plot(a/365,Cv(:,end),'-.c');
% plot(a/365,Ctot(:,end),'-.k');
% xlabel('age (years)')
% legend('Acquired','Maternal','Vaccine-derived','Total','Location','SouthEast');
% title(['Immun dist.~~ feedback =',num2str(immunity_feedback)]);
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
% toc
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
% toc
