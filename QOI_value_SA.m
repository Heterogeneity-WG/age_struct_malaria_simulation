function Q_val = QOI_value_SA(lQ,time_points)
global P

da = P.da;
a = P.a;

for iQ = 1:length(lQ)
    if strcmp(lQ{iQ}(1:2),'EE')
        flag_EE = 1;
        break
    end
end

if flag_EE
    [SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');
    [SH_solu, EH_solu, DH_solu, AH_solu, VH_solu, UH_solu, SM_solu, EM_solu, IM_solu, Cm_solu, Cac_solu, Cv_solu, Ctot_solu, MH_solu] = ...
        age_structured_Malaria_vac(P.da,P.na,P.tfinal,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);   
    SH = SH_solu(:,end); EH = EH_solu(:,end); DH = DH_solu(:,end); AH = AH_solu(:,end); MH = MH_solu(:,end);
    VH = VH_solu(:,end); UH = UH_solu(:,end); Cm = Cm_solu(:,end); Cac = Cac_solu(:,end); Cv = Cv_solu(:,end);
    Ctot = Ctot_solu(:,end); 
    PH = SH+EH+DH+AH+VH+UH;
    SM = SM_solu(:,end); EM = EM_solu(:,end); IM = IM_solu(:,end);
end

Q_val = NaN(length(time_points),length(lQ));
ind0210y = age_range_ind(a,2,10);
ind0924m = age_range_ind(a,9/12,24/12);

for iQ = 1:length(lQ)
    switch lQ{iQ}
        case 'R0'
            Q_val(:,iQ)  = R0_cal();
        case 'RHM'
            [~,Q_val(:,iQ),~]  = R0_cal();
        case 'RMH'
            [~,~,Q_val(:,iQ)]  = R0_cal();     
        case 'EE-D'
            Q_val(:,iQ) = trapz(DH)*da;
        case 'EE-D-02-10'
            Q_val(:,iQ) = trapz(DH(ind0210y))*da;
        case 'EE-D-09-24'
            Q_val(:,iQ) = trapz(DH(ind0924m))*da;
        case 'EE-DA'
            Q_val(:,iQ) = trapz(DH+AH)*da;    
        case 'EE-DA-02-10'            
            Q_val(:,iQ) = trapz(DH(ind0210y)+AH(ind0210y))*da;  
        case 'EE-DA-09-24'
            Q_val(:,iQ) = trapz(DH(ind0924m)+AH(ind0924m))*da;  
        case 'EE-D-frac'
            Q_val(:,iQ) = trapz(DH)/trapz(DH+AH);
        case 'EE-D-frac-02-10'
            Q_val(:,iQ) = trapz(DH(ind0210y))/trapz(DH(ind0210y)+AH(ind0210y));
        case 'EE-D-frac-09-24'
            Q_val(:,iQ) = trapz(DH(ind0924m))/trapz(DH(ind0924m)+AH(ind0924m));
        case 'EE-death' % Cumulative disease-induced mortality, diagnostic eqn MH
            Q_val(:,iQ) = trapz(MH)*da;
        case 'EE-death-02-10' % Cumulative disease-induced mortality, diagnostic eqn MH
            Q_val(:,iQ) = trapz(MH(ind0210y))*da;
        case 'EE-death-09-24' % Cumulative disease-induced mortality, diagnostic eqn MH
            Q_val(:,iQ) = trapz(MH(ind0924m))*da;
        case 'EE-EIR'
            NH = trapz(PH)*da;
            NM = SM+EM+IM;
            [bH,~] = biting_rate(NH,NM);
            IM_frac = IM./NM;
            Q_val(:,iQ) = bH.*IM_frac*365; % annual EIR
        case 'EE-Ctot-pp'
            Q_val(:,iQ) = trapz(Ctot)/trapz(PH);
        case 'EE-Ctot-pp-02-10'
            Q_val(:,iQ) = trapz(Ctot(ind0210y))/trapz(PH(ind0210y));
        case 'EE-Ctot-pp-09-24'
            Q_val(:,iQ) = trapz(Ctot(ind0924m))/trapz(PH(ind0924m));
        case 'year5-D'
            tfinal = 5*365;
            [SH_solu, EH_solu, DH_solu, AH_solu, VH_solu, UH_solu, SM_solu, EM_solu, IM_solu, Cm_solu, Cac_solu, Cv_solu, Ctot_solu, MH_solu] = ...
                age_structured_Malaria_vac(P.da,P.na,tfinal,SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH);   
            DH = DH_solu(:,end);
            Q_val(:,iQ) = trapz(DH)*da;
        otherwise
            keyboard
    end
end


end

function ind = age_range_ind(a,a_start,a_end)

[~,ind1] = min(abs(a-a_start*365)); % start from a_start years old
[~,ind2] = min(abs(a-a_end*365)); % end at a_end years old

ind = ind1:ind2;
end