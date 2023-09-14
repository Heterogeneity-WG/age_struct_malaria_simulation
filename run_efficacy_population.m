% close all
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

%% time evolution
P.v0 = 0; Malaria_parameters_transform_vac; % resetting vaccination distribution
tfinal = 3*365;
[t,SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(da,na,0,tfinal,...
    SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
PH = SH+EH+DH+AH+VH+UH;
PH_final = PH(:,end); % total human at age a, t = n
NH = trapz(PH,1)*da;
NM = SM+IM+EM;

%% vaccine efficacy - reduction in incidence with or w/o vaccine (on population level)
% initial condition = at the end of initial run
% ---- vaccinated group -----
lgroup = 'DH'; tfinal_vacc = 365; tfinal_count = 3650;  % -> counting incidence for one year
P.v0 = 100; % turn on vaccination
Malaria_parameters_transform_vac;
[~,ind1] = min(abs(P.a-7*30));
[~,ind2] = min(abs(P.a-19*30));
[~,vacc,~,SHv, EHv, DHv, AHv, VHv, UHv, SMv, EMv, IMv, Cmv, Cacv, Cvv, Ctotv, MHv] = incidence_cal(da,na,0,tfinal_vacc,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0,ind1,ind2,lgroup);
P.v0 = 0; % turn off vaccination
Malaria_parameters_transform_vac;
Incidence_vacc = incidence_cal_time(da,na,0,tfinal_count,SHv, EHv, DHv, AHv, VHv, UHv, SMv, EMv, IMv, Cmv, Cacv, Cvv, Ctotv, MHv,ind1,ind2,lgroup);
%---- control group -----
P.v0 = 0;
Malaria_parameters_transform_vac;
[~,~,~,SHc, EHc, DHc, AHc, VHc, UHc, SMc, EMc, IMc, Cmc, Cacc, Cvc, Ctotc, MHc] = incidence_cal(da,na,0,tfinal_vacc,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0,ind1,ind2,lgroup);
Incidence_control = incidence_cal_time(da,na,0,tfinal_count,SHc, EHc, DHc, AHc, VHc, UHc, SMc, EMc, IMc, Cmc, Cacc, Cvc, Ctotc, MHc,ind1,ind2,lgroup);
eff = (Incidence_control-Incidence_vacc)./vacc;

%% vaccine efficacy plots
figure_setups;
subplot(1,2, 1), plot((0:P.dt:tfinal_count),eff), axis tight;
grid on;
xlabel('Time');
title('Incidence reduction per vaccine');

subplot(1,2, 2), plot((0:P.dt:tfinal_count),Incidence_control,(0:P.dt:tfinal_count),Incidence_vacc), axis tight;
grid on;
xlabel('Time');
legend('Control','With vaccination');
title('Incidence with/without vaccination');

%% output data to .mat file for analysis
% SH_EE = SH(:,end); EH_EE = EH(:,end); AH_EE = AH(:,end); DH_EE = DH(:,end); VH_EE = VH(:,end); UH_EE = UH(:,end); PH_EE = PH(:,end); v = P.v;
% Cm_EE = Cm(:,end); Cac_EE = Cac(:,end); Ctot_EE = Ctot(:,end);
% save(['Results/Vaccine/v0_',num2str(P.v0*100),'.mat'],'t','a','v','vacc_sterile','vacc_blood','SH_EE','EH_EE','AH_EE','DH_EE','VH_EE','UH_EE','PH_EE',...
%     'Cm_EE','Cac_EE','Ctot_EE');
