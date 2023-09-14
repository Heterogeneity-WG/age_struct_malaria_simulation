% close all
clear all
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

Malaria_parameters_baseline;
Malaria_parameters_transform;
Malaria_parameters_transform_vac;

%% initial condition 'EE' - numerical EE
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');

%% time evolution
P.v0 = 0; Malaria_parameters_transform_vac; % resetting vaccination distribution
tfinal = 365*5;
[t,SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(da,na,0,tfinal,...
    SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
PH = SH+EH+DH+AH+VH+UH;
PH_final = PH(:,end);
NH = trapz(PH,1)*da;
NM = SM+IM+EM;
vacc_sterile = trapz(P.v*(1-P.z).*SH,1)*P.da;

%% ---- simulation for vac control on and off
% P.v0 = 10; Malaria_parameters_transform_vac;
% vacc_fun = P.v;
% tfinal_conti = 365*20;
% SH0 = SH(:,end); EH0 = EH(:,end); DH0 = DH(:,end); AH0 = AH(:,end); VH0 = VH(:,end); UH0 = UH(:,end); MH0 = MH(:,end);
% SM0 = SM(end); EM0 = EM(end); IM0 = IM(end);
% Cac0 = Cac(:,end); Cm0 = Cm(:,end); Cv0 = Cv(:,end); Ctot0 = Ctot(:,end);
% [t2,SH2, EH2, DH2, AH2, VH2, UH2, SM2, EM2, IM2, Cm2, Cac2, Cv2, Ctot2, MH2] = age_structured_Malaria_vac(da,na,t(end),t(end)+tfinal_conti,...
%     SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
% PH2 = SH2+EH2+DH2+AH2+VH2+UH2;
% NH2 = trapz(PH2,1)*da;
% NM2 = SM2+EM2+IM2;
% vacc_sterile2 = trapz(P.v*(1-P.z).*SH2,1)*P.da;
% 
% t = [t;t2(2:end)];
% SH = [SH, SH2(:,2:end)]; EH = [EH,EH2(:,2:end)]; DH = [DH, DH2(:,2:end)]; AH = [AH, AH2(:,2:end)]; VH = [VH, VH2(:,2:end)]; UH = [UH, UH2(:,2:end)]; MH = [MH, MH2(:,2:end)];
% SM = [SM, SM2(2:end)]; EM = [EM, EM2(2:end)]; IM = [IM, IM2(2:end)]; NM = [NM, NM2(2:end)];
% Cm = [Cm, Cm2(:,2:end)]; Cac = [Cac, Cac2(:,2:end)]; Cv = [Cv, Cv2(:,2:end)]; Ctot = [Ctot,Ctot2(:,2:end)];
% PH = SH+EH+DH+AH+VH+UH;
% PH_final = PH(:,end); % total human at age a, t = n
% NH = [NH, NH2(2:end)];
% vacc_sterile = [vacc_sterile, vacc_sterile2(2:end)];

% --- continued run - vacc off
% P.v0 = 0;
% Malaria_parameters_transform_vac;
% tfinal_conti2 = 365-t(end);
% SH0 = SH(:,end); EH0 = EH(:,end); DH0 = DH(:,end); AH0 = AH(:,end); VH0 = VH(:,end); UH0 = UH(:,end); MH0 = MH(:,end);
% SM0 = SM(end); EM0 = EM(end); IM0 = IM(end);
% Cac0 = Cac(:,end); Cm0 = Cm(:,end); Cv0 = Cv(:,end); Ctot0 = Ctot(:,end);
% [t2,SH2, EH2, DH2, AH2, VH2, UH2, SM2, EM2, IM2, Cm2, Cac2, Cv2, Ctot2, MH2] = age_structured_Malaria_vac(da,na,t(end),t(end)+tfinal_conti2,...
%     SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
% PH2 = SH2+EH2+DH2+AH2+VH2+UH2;
% NH2 = trapz(PH2,1)*da;
% NM2 = SM2+EM2+IM2;
% vacc_sterile2 = trapz(P.v*(1-P.z).*SH2,1)*P.da;
% 
% t = [t;t2(2:end)];
% SH = [SH,SH2(:,2:end)]; EH = [EH,EH2(:,2:end)]; DH = [DH, DH2(:,2:end)]; AH = [AH, AH2(:,2:end)]; VH = [VH, VH2(:,2:end)]; UH = [UH, UH2(:,2:end)]; MH = [MH, MH2(:,2:end)];
% SM = [SM, SM2(2:end)]; EM = [EM, EM2(2:end)]; IM = [IM, IM2(2:end)]; NM = [NM, NM2(2:end)];
% Cm = [Cm, Cm2(:,2:end)]; Cac = [Cac, Cac2(:,2:end)]; Cv = [Cv, Cv2(:,2:end)]; Ctot = [Ctot,Ctot2(:,2:end)];
% PH = SH+EH+DH+AH+VH+UH;
% PH_final = PH(:,end); % total human at age a, t = n
% NH = [NH, NH2(2:end)];
% vacc_sterile = [vacc_sterile, vacc_sterile2(2:end)];

%%
[bH,~] = biting_rate(PH,NM); % bH(age, time) matrix
EIR = bH.*IM./NM*365; % EIR(age, time) matrix
EIR_tot = trapz(EIR.*PH,1)*P.da./NH;
EIR_final = EIR_tot(end);

%% plotting
% plot_EIR_time;
% plot_human_prop_time; % human proportion in time (total and by age group) 
% plot_human_popsize_time; % human pop size in time (total and by age group) 
plot_human_pop_age_tfinal; % age distributions (in prop and size) at final time (and movie)
% plot_immunity; 
% plot_DALY;
% plot_sigmoids;
% plot_incidences;
% plot_seasonality;
% plot_mosquitoes;
% plot_vacc_counts_time;



