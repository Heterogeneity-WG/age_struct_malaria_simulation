function Q_val = QOI_value_eFAST(lQ)
global P 

da = P.da;
a = P.a;
if strcmp(lQ(1:2),'EE')
    if R0_cal()<1
        disp('EE DNE')
        Q_val = NaN;
        return
    end
    [SH,EH,DH,AH,Cac,Cm,Ctot] = steady_state('EE','numerical');
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
