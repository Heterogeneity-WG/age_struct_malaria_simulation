% repeat time evolution from tnow to tnow+tfinal
[t2,SH2, EH2, DH2, AH2, VH2, UH2, SM2, EM2, IM2, Cm2, Cac2, Cv2, Ctot2, MH2] = age_structured_Malaria_vac(da,na,tnow,tnow+tconti,...
    SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
vacc_sterile2 = trapz(P.v*(1-P.z).*SH2,1)*P.da;

t = [t;t2(2:end)];
SH = [SH,SH2(:,2:end)]; EH = [EH,EH2(:,2:end)]; DH = [DH, DH2(:,2:end)]; AH = [AH, AH2(:,2:end)]; VH = [VH, VH2(:,2:end)]; UH = [UH, UH2(:,2:end)]; MH = [MH, MH2(:,2:end)];
SM = [SM, SM2(2:end)]; EM = [EM, EM2(2:end)]; IM = [IM, IM2(2:end)];
Cm = [Cm, Cm2(:,2:end)]; Cac = [Cac, Cac2(:,2:end)]; Cv = [Cv, Cv2(:,2:end)]; Ctot = [Ctot,Ctot2(:,2:end)];
vacc_sterile = [vacc_sterile, vacc_sterile2(2:end)];

tnow = t(end);

SH0 = SH(:,end); EH0 = EH(:,end); DH0 = DH(:,end); AH0 = AH(:,end); VH0 = VH(:,end); UH0 = UH(:,end); MH0 = MH(:,end);
SM0 = SM(end); EM0 = EM(end); IM0 = IM(end);
Cac0 = Cac(:,end); Cm0 = Cm(:,end); Cv0 = Cv(:,end); Ctot0 = Ctot(:,end);
