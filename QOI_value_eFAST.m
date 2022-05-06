function Q_val = QOI_value_eFAST(lQ,time_points)
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
    %     [SH,EH,DH,AH,Cac,Cm,Ctot] = steady_state('EE','numerical');
    [SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0] = age_structured_Malaria_IC('init');
    [SH_solu, EH_solu, DH_solu, AH_solu, SM_solu, EM_solu, IM_solu, ~, ~, Ctot_solu] = age_structured_Malaria(P.da,P.na,P.tfinal,SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0);
    SH = SH_solu(:,end); EH = EH_solu(:,end); DH = DH_solu(:,end); AH = AH_solu(:,end); Ctot = Ctot_solu(:,end);
    SM = SM_solu(:,end); EM = EM_solu(:,end); IM = IM_solu(:,end);
end

Q_val = NaN(length(time_points),length(lQ));
ind0210 = age_range_ind(a,2,10);

for iQ = 1:length(lQ)
    switch lQ{iQ}
        case 'R0'
            Q_val(:,iQ)  = R0_cal();
        case 'RHM'
            [~,Q_val(:,iQ),~]  = R0_cal();
        case 'RMH'
            [~,~,Q_val(:,iQ)]  = R0_cal();
        case 'EE-infected'
            Q_val(:,iQ) = 1-da*trapz(SH);       
        case 'EE-D'
            Q_val(:,iQ) = trapz(DH)*da;
        case 'EE-D-02-10'
            Q_val(:,iQ) = trapz(DH(ind0210))*da;
        case 'EE-DA'
            Q_val(:,iQ) = trapz(DH+AH)*da;    
        case 'EE-DA-02-10'            
            Q_val(:,iQ) = trapz(DH(ind0210)+AH(ind0210))*da;  
        case 'EE-D-frac'
            Q_val(:,iQ) = trapz(DH)/trapz(DH+AH);
        case 'EE-D-frac-02-10'
            Q_val(:,iQ) = trapz(DH(ind0210))/trapz(DH(ind0210)+AH(ind0210));
        case 'EE-EDA'
            NH = trapz(SH+EH+DH+AH)*da;
            Q_val(:,iQ) = [trapz(EH)*da/NH; trapz(DH)*da/NH; trapz(AH)*da/NH];
        case 'EE-EIR'
            NH = trapz(SH+EH+DH+AH)*P.da;
            NM = SM+EM+IM;
            [bH,~] = biting_rate(NH,NM);
            IM_frac = IM./NM;
            Q_val(:,iQ) = bH.*IM_frac*365; % annual EIR
        case 'EE_Ctot_pp'
            Q_val(:,iQ) = Ctot./P.PH_stable;
        case 'EE_rho'
            Q_val(:,iQ) = sigmoid_prob(Ctot./P.PH_stable, 'rho');
        otherwise
            keyboard
    end
end


end

function ind = age_range_ind(a,a_start,a_end)

[~,ind1] = min(abs(a-a_start*365)); % start from 2 years old
[~,ind2] = min(abs(a-a_end*365)); % end at 10 years old

ind = ind1:ind2;
end