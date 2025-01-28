%% vac = constant, fixed number, with booster
close all
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
% Malaria_parameters_baseline_Nanoro;
Malaria_parameters_transform;
Malaria_parameters_transform_vac;

%% initial condition 'EE' - numerical EE
P.v0 = 0; 
Malaria_parameters_transform_vac;
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');

%% time evolution
P.v0 = 0; P.vb0 = 0; Malaria_parameters_transform_vac; % resetting vaccination distribution
tfinal = 365*3;
temp0 = zeros(size(SH0));
[t,SHv, EHv, DHv, AHv, VHv, UHv, Cmv, Cacv, Ctotv, SHb, EHb, DHb, AHb, VHb, UHb, Cmb, Cacb, Ctotb, SM, EM, IM, vb_temp] = ...
    age_structured_Malaria_booster(da,na,0,tfinal,...
    SH0, EH0, DH0, AH0, VH0, UH0, Cm0, Cac0, Ctot0, ...
    temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, SM0, EM0, IM0);
PH = SHv+EHv+DHv+AHv+VHv+UHv+SHb+EHb+DHb+AHb+VHb+UHb;
% PH_final = PH(:,end);
NH = trapz(PH,1)*da;
NM = SM+IM+EM;
vacc_sterile = repmat(P.v,1,length(t)); % record vac rate
vacc_booster = trapz(vb_temp,1)*P.da;

%% ---- simulation for vac control on and off
% vaccinate at least 1.2*10^5 children per year in the selected areas
% target population total = 2.56*10^5;  max = 2.2*10^5/365 to avoid negative SH
P.v0 = 2.2*10^5/365; P.vb0 = P.v0; Malaria_parameters_transform_vac; 
tfinal_conti = 365*5;
SH0 = SHv(:,end); EH0 = EHv(:,end); DH0 = DHv(:,end); AH0 = AHv(:,end); VH0 = VHv(:,end); UH0 = UHv(:,end); 
Cac0 = Cacv(:,end); Cm0 = Cmv(:,end); Ctot0 = Ctotv(:,end);
SHb0 = SHb(:,end); EHb0 = EHb(:,end); DHb0 = DHb(:,end); AHb0 = AHb(:,end); VHb0 = VHb(:,end); UHb0 = UHb(:,end); 
Cacb0 = Cacb(:,end); Cmb0 = Cmb(:,end); Ctotb0 = Ctotb(:,end);
SM0 = SM(end); EM0 = EM(end); IM0 = IM(end);
[t2,SHv2, EHv2, DHv2, AHv2, VHv2, UHv2, Cmv2, Cacv2, Ctotv2, SHb2, EHb2, DHb2, AHb2, VHb2, UHb2, Cmb2, Cacb2, Ctotb2, SM2, EM2, IM2, vb_temp2] = ...
    age_structured_Malaria_booster(da,na,t(end),t(end)+tfinal_conti,...
    SH0, EH0, DH0, AH0, VH0, UH0, Cm0, Cac0, Ctot0, ...
    SHb0, EHb0, DHb0, AHb0, VHb0, UHb0, Cmb0, Cacb0, Ctotb0, SM0, EM0, IM0);
PH2 = SHv2+EHv2+DHv2+AHv2+VHv2+UHv2+SHb2+EHb2+DHb2+AHb2+VHb2+UHb2;
NH2 = trapz(PH2,1)*da;
NM2 = SM2+EM2+IM2;
vacc_sterile2 = repmat(P.v,1,length(t2));
vacc_booster2 = trapz(vb_temp2,1)*P.da;

t = [t;t2(2:end)];
SHv = [SHv, SHv2(:,2:end)]; EHv = [EHv,EHv2(:,2:end)]; DHv = [DHv, DHv2(:,2:end)]; AHv = [AHv, AHv2(:,2:end)]; VHv = [VHv, VHv2(:,2:end)]; UHv = [UHv, UHv2(:,2:end)];
Cmv = [Cmv, Cmv2(:,2:end)]; Cacv = [Cacv, Cacv2(:,2:end)]; Ctotv = [Ctotv,Ctotv2(:,2:end)];
SHb = [SHb, SHb2(:,2:end)]; EHb = [EHb,EHb2(:,2:end)]; DHb = [DHb, DHb2(:,2:end)]; AHb = [AHb, AHb2(:,2:end)]; VHb = [VHb, VHb2(:,2:end)]; UHb = [UHb, UHb2(:,2:end)];
Cmb = [Cmb, Cmb2(:,2:end)]; Cacb = [Cacb, Cacb2(:,2:end)]; Ctotb = [Ctotb,Ctotb2(:,2:end)];
SM = [SM, SM2(2:end)]; EM = [EM, EM2(2:end)]; IM = [IM, IM2(2:end)]; NM = [NM, NM2(2:end)];
PH = SHv+EHv+DHv+AHv+VHv+UHv+SHb+EHb+DHb+AHb+VHb+UHb;
PH_final = PH(:,end); % total human at age a, t = n
NH = [NH, NH2(2:end)];
vacc_sterile = [vacc_sterile, vacc_sterile2(:,2:end)];

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
% vacc_sterile2 = repmat(P.v,1,length(t2));
% 
% t = [t;t2(2:end)];
% SH = [SH,SH2(:,2:end)]; EH = [EH,EH2(:,2:end)]; DH = [DH, DH2(:,2:end)]; AH = [AH, AH2(:,2:end)]; VH = [VH, VH2(:,2:end)]; UH = [UH, UH2(:,2:end)]; MH = [MH, MH2(:,2:end)];
% SM = [SM, SM2(2:end)]; EM = [EM, EM2(2:end)]; IM = [IM, IM2(2:end)]; NM = [NM, NM2(2:end)];
% Cm = [Cm, Cm2(:,2:end)]; Cac = [Cac, Cac2(:,2:end)]; Cv = [Cv, Cv2(:,2:end)]; Ctot = [Ctot,Ctot2(:,2:end)];
% PH = SH+EH+DH+AH+VH+UH;
% PH_final = PH(:,end); % total human at age a, t = n
% NH = [NH, NH2(2:end)];
% vacc_sterile = [vacc_sterile, vacc_sterile2(:,2:end)];

%%
[bH,~] = biting_rate(PH,NM); % bH(age, time) matrix
EIR = bH.*IM./NM*365; % EIR(age, time) matrix
EIR_tot = trapz(EIR.*PH,1)*P.da./NH;
EIR_final = EIR_tot(end);
SH = SHv+SHb; EH = EHv+EHb; AH = AHv+AHb; DH = DHv+DHb; VH = VHv+VHb; UH = UHv+UHb;
%% plotting
% plot_EIR_time;
% plot_human_prop_time; % human proportion in time (total and by age group) 
% plot_human_popsize_time; % human pop size in time (total and by age group) 
plot_human_pop_age_tfinal; % age distributions (in prop and size) at final time (and movie)
% plot_immunity; 
% plot_DALY;
% plot_sigmoids;
% plot_FOI;
% plot_incidences;
% plot_seasonality;
% plot_mosquitoes;
% plot_vacc_counts_time;



