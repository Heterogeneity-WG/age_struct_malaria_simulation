function [err, xdata,ydata,yrun] = fun_efficacy(x,Data,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0)

global P 

P.etas = x(1);
P.w = x(2);

dt = P.dt;
da = P.da;
na = P.na;

tfinal_vacc = 2*365; total_vacc = 2000;
P.v0s = total_vacc/tfinal_vacc; P.v0c = P.v0s; % define constant vaccination rate
Malaria_parameters_transform_vac;
t = (0:dt:tfinal_vacc)';
t_index = length(t);
temp0 = zeros(size(SH0));
[SHr, EHr, DHr, AHr, Cmr, Cacr, Ctotr, SHv, EHv, DHv, AHv, VHv, UHv, Cmv, Cacv, Ctotv, SHc, EHc, DHc, AHc, Cmc, Cacc, Ctotc, SM, EM, IM] = ...
    age_structured_Malaria_eff(da, na, tfinal_vacc, SH0, EH0, DH0, AH0, Cm0, Cac0, Ctot0, ...
    temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, SM0, EM0, IM0);

%% simulation using three-group model - vaccine off
tfinal_conti = 5*365; total_vacc = 0;
P.v0s = total_vacc/tfinal_vacc; P.v0c = P.v0s; % define constant vaccination rate
Malaria_parameters_transform_vac;
t2 = (tfinal_vacc:dt:tfinal_vacc+tfinal_conti)';
[SHr2, EHr2, DHr2, AHr2, Cmr2, Cacr2, Ctotr2, SHv2, EHv2, DHv2, AHv2, VHv2, UHv2, Cmv2, Cacv2, Ctotv2, SHc2, EHc2, DHc2, AHc2, Cmc2, Cacc2, Ctotc2, SM2, EM2, IM2] = ...
    age_structured_Malaria_eff(da, na, tfinal_conti, SHr(:,end), EHr(:,end), DHr(:,end), AHr(:,end), Cmr(:,end), Cacr(:,end), Ctotr(:,end), ...
    SHv(:,end), EHv(:,end), DHv(:,end), AHv(:,end), VHv(:,end), UHv(:,end), Cmv(:,end), Cacv(:,end), Ctotv(:,end), ...
    SHc(:,end), EHc(:,end), DHc(:,end), AHc(:,end), Cmc(:,end), Cacc(:,end), Ctotc(:,end), SM(end), EM(end), IM(end));
% combine results
SHr = [SHr,SHr2]; EHr = [EHr,EHr2]; DHr = [DHr, DHr2]; AHr = [AHr, AHr2]; Cmr = [Cmr, Cmr2]; Cacr = [Cacr, Cacr2]; Ctotr = [Ctotr,Ctotr2];
SHc = [SHc,SHc2]; EHc = [EHc,EHc2]; DHc = [DHc, DHc2]; AHc = [AHc, AHc2]; Cmc = [Cmc, Cmc2]; Cacc = [Cacc, Cacc2]; Ctotc = [Ctotc,Ctotc2];
SHv = [SHv,SHv2]; EHv = [EHv,EHv2]; DHv = [DHv, DHv2]; AHv = [AHv, AHv2]; VHv = [VHv, VHv2]; UHv = [UHv, UHv2]; Cmv = [Cmv, Cmv2]; Cacv = [Cacv, Cacv2]; Ctotv = [Ctotv,Ctotv2];
SM = [SM, SM2]; EM = [EM, EM2]; IM = [IM, IM2];
t = [t;t2];

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

%% calculate efficacy
eff = (Incidence_control'-Incidence_vacc')./Incidence_control';

t(t_index) = []; % remove duplication
eff(t_index) = []; % remove duplication

eff(1)= []; % remove day 1 
t(1) = [];

ind = Data(:,1)>0;

xdata = Data(ind,1)*365; % time in years -> days
ydata = Data(ind,2);

yrun = interp1(t,eff,xdata); 

err = norm(ydata-yrun);

end