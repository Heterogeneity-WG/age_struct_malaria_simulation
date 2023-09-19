% repeat time evolution from tnow to tnow+tfinal
[t,SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(da,na,tnow,tnow+tconti,...
    SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, zeros(size(SH0)));
vacc_sterile = trapz(P.v*(1-P.z).*SH,1)*P.da;

tnow = t(end);

SH0 = SH(:,end); EH0 = EH(:,end); DH0 = DH(:,end); AH0 = AH(:,end); VH0 = VH(:,end); UH0 = UH(:,end); MH0 = MH(:,end);
SM0 = SM(end); EM0 = EM(end); IM0 = IM(end);
Cac0 = Cac(:,end); Cm0 = Cm(:,end); Cv0 = Cv(:,end); Ctot0 = Ctot(:,end);