clearvars
% close all
clc

format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7

%% numerical config
age_max = 100*365; % max ages in days
P.age_max = age_max;
dt = 10; % time/age step size in days, default = 5;
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
P.v0s = 0; P.v0c = 0;
Malaria_parameters_transform_vac;

%% initial condition 'EE' - numerical EE
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE');

tfinal_vacc = 5*365;  
P.v0s = 2000/(2*365); P.v0c = P.v0s; % define constant vaccination rate
Malaria_parameters_transform_vac;
t = (0:dt:tfinal_vacc)'; 
[SHr, EHr, DHr, AHr, Cmr, Cacr, Ctotr, SHv, EHv, DHv, AHv, VHv, UHv, Cmv, Cacv, Ctotv, SHc, EHc, DHc, AHc, Cmc, Cacc, Ctotc, SM, EM, IM] = ...
    age_structured_Malaria_eff3(da, na, tfinal_vacc, SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0);
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

eff = (Incidence_control'-Incidence_vacc')./Incidence_control';
figure_setups;
plot(t/365,eff)
xlabel('Year')
ylabel('Efficacy')

%% plots for debugging
SH = SHr+SHv+SHc;
EH = EHr+EHv+EHc;
AH = AHr+AHv+AHc;
DH = DHr+DHv+DHc;
VH = VHv; UH = UHv;
PH_final = PH(:,end);
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
