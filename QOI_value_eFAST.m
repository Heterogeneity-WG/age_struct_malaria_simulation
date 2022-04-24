function Q_val = QOI_value_eFAST(lQ)
global P 

da = P.da;
a = P.a;

if strcmp(lQ(1:2),'EE')
%     [SH,EH,DH,AH,Cac,Cm,Ctot] = steady_state('EE','numerical');
    [SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0] = age_structured_Malaria_IC('init');
    [SH_solu, EH_solu, DH_solu, AH_solu, SM_solu, EM_solu, IM_solu, ~, ~, Ctot_solu] = age_structured_Malaria(P.da,P.na,P.tfinal,SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0);
    SH = SH_solu(:,end); EH = EH_solu(:,end); DH = DH_solu(:,end); AH = AH_solu(:,end); Ctot = Ctot_solu(:,end);
    SM = SM_solu(:,end); EM = EM_solu(:,end); IM = IM_solu(:,end); 
 end

switch lQ
    case 'R0'    
        Q_val  = R0_cal();
    case 'RHM'
        [~,Q_val,~]  = R0_cal();
    case 'RMH'
        [~,~,Q_val]  = R0_cal();
    case 'EE-infected'
        Q_val = 1-da*trapz(SH);       
    case 'EE-DA'
        Q_val = 1-da*trapz(SH)-da*trapz(EH);
    case 'EE-D-frac'
        Q_val = trapz(DH)/trapz(DH+AH);
    case 'EE-D'
        Q_val = trapz(DH)*da;
    case 'EE-EDA'
        NH = trapz(SH+EH+DH+AH)*da; 
        Q_val = [trapz(EH)*da/NH; trapz(DH)*da/NH; trapz(AH)*da/NH];
    case 'EE-EIR'
        NH = trapz(SH+EH+DH+AH)*da; 
        NM = P.gM/P.muM;
        [bH,bM] = biting_rate(NH,NM);
        Lambda_M = bM*trapz(P.betaD*DH + P.betaA*AH)*da;
        IM_frac_EE = P.sigma/(P.sigma+P.muM)*(Lambda_M/(Lambda_M + P.muM));   
        Q_val = bH*IM_frac_EE*365; % annual EIR           
    case 'EE_Ctot_pp'
        Q_val = Ctot./P.PH_stable;
    case 'EE_rho'
        Q_val = sigmoid_prob(Ctot./P.PH_stable, 'rho');
    otherwise
        keyboard
end



end
