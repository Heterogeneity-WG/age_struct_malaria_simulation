function [y,vacc_sterile, vacc_blood,SH1, EH1, DH1, AH1, VH1, UH1, SM1, EM1, IM1, Cm1, Cac1, Cv1, Ctot1, MH1] =...
    incidence_cal(da,na,tinitial,tfinal_conti,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0,ind1,ind2,group)
global P

[~, SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(da,na,tinitial,tfinal_conti,...
    SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
PH = SH+EH+DH+AH+VH+UH;
%NH = trapz(PH,1)*da;
NM = SM+EM+IM;
[bH,~] = biting_rate(PH,NM);
lamH = FOI_H(bH,IM,NM);
vacc_sterile = trapz(trapz(P.v*(1-P.z).*SH,1)*P.da)*P.dt;
vacc_blood = trapz(trapz(P.v*P.z.*SH,1)*P.da)*P.dt;

switch group
    case 'SH' % new infection       
        y = trapz(trapz(lamH.*SH(ind1:ind2,:),1)*P.da,2)*P.dt;
    case 'DH' % new clinical infection
        rho = sigmoid_prob(Ctot./PH, 'rho'); % prob. of severely infected, EH -> DH
        psi = sigmoid_prob(Ctot./PH, 'psi'); % prob. AH -> DH
        temp = rho.*P.h.*EH+psi.*lamH.*AH;
        y = trapz(trapz(temp(ind1:ind2,:),1)*P.da,2)*P.dt;
end

SH1 = SH(:,end); EH1 = EH(:,end); DH1 = DH(:,end); AH1 = AH(:,end); VH1 = VH(:,end); UH1 = UH(:,end); MH1 = MH(:,end);
SM1 = SM(end); EM1 = EM(end); IM1 = IM(end); 
Cac1 = Cac(:,end); Cm1 = Cm(:,end); Cv1 = Cv(:,end); Ctot1 = Ctot(:,end);

end