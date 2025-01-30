% repeat time evolution from tnow to tnow+tconti
% do not preserve the history of solutions
[t, SHv, EHv, DHv, AHv, VHv, UHv, Cmv, Cacv, Ctotv, SHb, EHb, DHb, AHb, VHb, UHb, Cmb, Cacb, Ctotb, SM, EM, IM, vb_temp] = ...
    age_structured_Malaria_booster(da,na,tnow,tnow+tconti,...
    SH0, EH0, DH0, AH0, VH0, UH0, Cm0, Cac0, Ctot0, ...
    SHb0, EHb0, DHb0, AHb0, VHb0, UHb0, Cmb0, Cacb0, Ctotb0, SM0, EM0, IM0);
vacc_sterile = repmat(P.v,1,length(t));
vacc_booster = vb_temp;
tnow = t(end);

SH0 = SHv(:,end); EH0 = EHv(:,end); DH0 = DHv(:,end); AH0 = AHv(:,end); VH0 = VHv(:,end); UH0 = UHv(:,end); 
Cac0 = Cacv(:,end); Cm0 = Cmv(:,end); Ctot0 = Ctotv(:,end);
SHb0 = SHb(:,end); EHb0 = EHb(:,end); DHb0 = DHb(:,end); AHb0 = AHb(:,end); VHb0 = VHb(:,end); UHb0 = UHb(:,end); 
Cacb0 = Cacb(:,end); Cmb0 = Cmb(:,end); Ctotb0 = Ctotb(:,end);
SM0 = SM(end); EM0 = EM(end); IM0 = IM(end);