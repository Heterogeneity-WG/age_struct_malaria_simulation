function [y,vacc_sterile, vacc_blood] = incidence_cal(da,na,tfinal_conti,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0,ind1,ind2,group)
global P

[SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(da,na,tfinal_conti,...
    SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
PH = SH+EH+DH+AH+VH+UH;
NH = trapz(PH,1)*da;
NM = SM+EM+IM;
[bH,~] = biting_rate(NH,NM);
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
end