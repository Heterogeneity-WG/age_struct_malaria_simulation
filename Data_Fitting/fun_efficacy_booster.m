function [err, xdata,ydata,yrun, t,eff] = fun_efficacy_booster(x,Data,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0)

global P 

P.etab = x(1);
P.wb = x(2);

dt = P.dt;
da = P.da;
na = P.na;

%% simulation using three-group model - vaccine on
tfinal_vacc = 2*365; total_vacc = 2000;
P.v0s = total_vacc/tfinal_vacc; P.v0c = P.v0s; % define constant vaccination rate
P.vb0 = 1; % vaccinate all the eligible people
Malaria_parameters_transform_vac;
t = (0:dt:tfinal_vacc)';
temp0 = zeros(size(SH0));
[~,SHr, EHr, DHr, AHr, Cmr, Cacr, Ctotr, SHv, EHv, DHv, AHv, VHv, UHv, Cmv, Cacv, Ctotv,...
    SHb, EHb, DHb, AHb, VHb, UHb, Cmb, Cacb, Ctotb, SHc, EHc, DHc, AHc, Cmc, Cacc, Ctotc, SM, EM, IM] = ...
    age_structured_Malaria_eff_booster(da, na, 0, tfinal_vacc, SH0, EH0, DH0, AH0, Cm0, Cac0, Ctot0, ...% rest
    temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0,...%vacc 
    temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0, temp0,...%boost 
    temp0, temp0, temp0, temp0, temp0, temp0, temp0,...%control
    SM0, EM0, IM0);

%% simulation using three-group model - vaccine off the primary dose (keep vacc booster)
tfinal_conti = 5*365; total_vacc = 0;
P.v0s = total_vacc/tfinal_vacc; P.v0c = P.v0s; % define constant vaccination rate
Malaria_parameters_transform_vac;
t2 = (tfinal_vacc:dt:tfinal_vacc+tfinal_conti)';
[~,SHr2, EHr2, DHr2, AHr2, Cmr2, Cacr2, Ctotr2, SHv2, EHv2, DHv2, AHv2, VHv2, UHv2, Cmv2, Cacv2, Ctotv2, ...
    SHb2, EHb2, DHb2, AHb2, VHb2, UHb2, Cmb2, Cacb2, Ctotb2, SHc2, EHc2, DHc2, AHc2, Cmc2, Cacc2, Ctotc2, SM2, EM2, IM2] = ...
    age_structured_Malaria_eff_booster(da, na, 0, tfinal_conti, SHr(:,end), EHr(:,end), DHr(:,end), AHr(:,end), Cmr(:,end), Cacr(:,end), Ctotr(:,end), ...
    SHv(:,end), EHv(:,end), DHv(:,end), AHv(:,end), VHv(:,end), UHv(:,end), Cmv(:,end), Cacv(:,end), Ctotv(:,end), ...
    SHb(:,end), EHb(:,end), DHb(:,end), AHb(:,end), VHb(:,end), UHb(:,end), Cmb(:,end), Cacb(:,end), Ctotb(:,end), ...
    SHc(:,end), EHc(:,end), DHc(:,end), AHc(:,end), Cmc(:,end), Cacc(:,end), Ctotc(:,end), SM(end), EM(end), IM(end));
% combine results
SHr = [SHr,SHr2(:,2:end)]; EHr = [EHr,EHr2(:,2:end)]; DHr = [DHr, DHr2(:,2:end)]; AHr = [AHr, AHr2(:,2:end)]; Cmr = [Cmr, Cmr2(:,2:end)]; Cacr = [Cacr, Cacr2(:,2:end)]; Ctotr = [Ctotr,Ctotr2(:,2:end)];
SHc = [SHc,SHc2(:,2:end)]; EHc = [EHc,EHc2(:,2:end)]; DHc = [DHc, DHc2(:,2:end)]; AHc = [AHc, AHc2(:,2:end)]; Cmc = [Cmc, Cmc2(:,2:end)]; Cacc = [Cacc, Cacc2(:,2:end)]; Ctotc = [Ctotc,Ctotc2(:,2:end)];
SHv = [SHv,SHv2(:,2:end)]; EHv = [EHv,EHv2(:,2:end)]; DHv = [DHv, DHv2(:,2:end)]; AHv = [AHv, AHv2(:,2:end)]; VHv = [VHv, VHv2(:,2:end)]; UHv = [UHv, UHv2(:,2:end)]; Cmv = [Cmv, Cmv2(:,2:end)]; Cacv = [Cacv, Cacv2(:,2:end)]; Ctotv = [Ctotv,Ctotv2(:,2:end)];
SHb = [SHb,SHb2(:,2:end)]; EHb = [EHb,EHb2(:,2:end)]; DHb = [DHb, DHb2(:,2:end)]; AHb = [AHb, AHb2(:,2:end)]; VHb = [VHb, VHb2(:,2:end)]; UHb = [UHb, UHb2(:,2:end)]; Cmb = [Cmb, Cmb2(:,2:end)]; Cacb = [Cacb, Cacb2(:,2:end)]; Ctotb = [Ctotb,Ctotb2(:,2:end)];
SM = [SM, SM2(2:end)]; EM = [EM, EM2(2:end)]; IM = [IM, IM2(2:end)];
t = [t;t2(2:end)];
%% calculate incidences
PHr = SHr+EHr+DHr+AHr;
PHc = SHc+EHc+DHc+AHc;
PHv = SHv+EHv+DHv+AHv+VHv+UHv;
PHb = SHb+EHb+DHb+AHb+VHb+UHb;
PH = PHr+PHc+PHv+PHb;
NM = SM+EM+IM;
[bH,~] = biting_rate(PH,NM);
lamH = FOI_H(bH,IM,NM);
rhov = sigmoid_prob(Ctotv./PHv, 'rho'); % prob. of severely infected, EH -> DH
psiv = sigmoid_prob(Ctotv./PHv, 'psi'); % prob. AH -> DH
rhov(PHv==0)=0; psiv(PHv==0)=0;
temp1 = rhov.*P.h.*EHv+psiv.*lamH.*AHv; % incidence of DH in Vacc group

rhob = sigmoid_prob(Ctotb./PHb, 'rho'); % prob. of severely infected, EH -> DH
psib = sigmoid_prob(Ctotb./PHb, 'psi'); % prob. AH -> DH
rhob(PHb==0)=0; psib(PHb==0)=0;
temp2 = rhob.*P.h.*EHb+psib.*lamH.*AHb; % incidence of DH in Boosted group

Incidence_vacc = trapz(temp1+temp2,1)*P.da;
Incidence_vacc_pp = Incidence_vacc./(trapz(PHv+PHb,1)*P.da);

rhoc = sigmoid_prob(Ctotc./PHc, 'rho'); % prob. of severely infected, EH -> DH
psic = sigmoid_prob(Ctotc./PHc, 'psi'); % prob. AH -> DH
rhoc(PHc==0)=0; psic(PHc==0)=0;
temp3 = rhoc.*P.h.*EHc+psic.*lamH.*AHc; % incidence of DH
Incidence_control = trapz(temp3,1)*P.da; 
Incidence_control_pp = Incidence_control./(trapz(PHc,1)*P.da);
%% calculate efficacy
xdata = Data(:,1)*365; % time in years -> days
ydata = Data(:,2);

%% using aggregated incidence (within three month period) for residual calculation
% Incidence_vacc_3mon = NaN(size(ydata));
% Incidence_control_3mon = NaN(size(ydata));
% for i = 1:length(xdata)
%     period_end = xdata(i);
%     period_beg = xdata(i)-30*3;
%     [~,ind_beg] = min(abs(t-period_beg)); ind_beg = ind_beg+1;
%     [~,ind_end] = min(abs(t-period_end)); 
%     Incidence_vacc_3mon(i) = trapz(Incidence_vacc(ind_beg:ind_end))*P.dt;
%     Incidence_control_3mon(i) = trapz(Incidence_control(ind_beg:ind_end))*P.dt;
% end
% eff_3mon = (Incidence_control_3mon'-Incidence_vacc_3mon')./Incidence_control_3mon';
% ydata(isnan(eff_3mon))=[];
% eff_3mon(isnan(eff_3mon))=[];
% yrun = eff_3mon';

%% using instantaneous incidence for residual calculation
% eff = (Incidence_control'-Incidence_vacc')./Incidence_control';
eff = (Incidence_control_pp'-Incidence_vacc_pp')./Incidence_control_pp';
t(isnan(eff))=[]; % remove NaN if needed
eff(isnan(eff))=[];
yrun = interp1(t,eff,xdata,'pchip'); 

err = norm(ydata-yrun);

end