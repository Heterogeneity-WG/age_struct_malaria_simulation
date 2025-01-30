% repeat time evolution from tnow to tnow+tconti
% preserve the history of solutions
[t2, SHv2, EHv2, DHv2, AHv2, VHv2, UHv2, Cmv2, Cacv2, Ctotv2, SHb2, EHb2, DHb2, AHb2, VHb2, UHb2, Cmb2, Cacb2, Ctotb2, SM2, EM2, IM2, vb_temp2] = ...
    age_structured_Malaria_booster(da,na,tnow,tnow+tconti,...
    SH0, EH0, DH0, AH0, VH0, UH0, Cm0, Cac0, Ctot0, ...
    SHb0, EHb0, DHb0, AHb0, VHb0, UHb0, Cmb0, Cacb0, Ctotb0, SM0, EM0, IM0);
vacc_sterile2 = repmat(P.v,1,length(t2));
vacc_booster2 = vb_temp2;

t = [t;t2(2:end)];
SHv = [SHv, SHv2(:,2:end)]; EHv = [EHv, EHv2(:,2:end)]; DHv = [DHv, DHv2(:,2:end)]; AHv = [AHv, AHv2(:,2:end)]; VHv = [VHv, VHv2(:,2:end)]; UHv = [UHv, UHv2(:,2:end)];
Cmv = [Cmv, Cmv2(:,2:end)]; Cacv = [Cacv, Cacv2(:,2:end)]; Ctotv = [Ctotv,Ctotv2(:,2:end)];
SHb = [SHb, SHb2(:,2:end)]; EHb = [EHb, EHb2(:,2:end)]; DHb = [DHb, DHb2(:,2:end)]; AHb = [AHb, AHb2(:,2:end)]; VHb = [VHb, VHb2(:,2:end)]; UHb = [UHb, UHb2(:,2:end)];
Cmb = [Cmb, Cmb2(:,2:end)]; Cacb = [Cacb, Cacb2(:,2:end)]; Ctotb = [Ctotb, Ctotb2(:,2:end)];
SM  = [SM, SM2(2:end)]; EM = [EM, EM2(2:end)]; IM = [IM, IM2(2:end)]; 

vacc_sterile = [vacc_sterile, vacc_sterile2(:,2:end)];
vacc_booster = [vacc_booster, vacc_booster2(:,2:end)];

tnow = t(end);

SH0 = SHv(:,end); EH0 = EHv(:,end); DH0 = DHv(:,end); AH0 = AHv(:,end); VH0 = VHv(:,end); UH0 = UHv(:,end); 
Cac0 = Cacv(:,end); Cm0 = Cmv(:,end); Ctot0 = Ctotv(:,end);
SHb0 = SHb(:,end); EHb0 = EHb(:,end); DHb0 = DHb(:,end); AHb0 = AHb(:,end); VHb0 = VHb(:,end); UHb0 = UHb(:,end); 
Cacb0 = Cacb(:,end); Cmb0 = Cmb(:,end); Ctotb0 = Ctotb(:,end);
SM0 = SM(end); EM0 = EM(end); IM0 = IM(end);
