function Q_val = QOI_value_SA_2(lQ,time_points,run_num,lmethod,direc)
% compress saving output, trapz on the age dimension to reduce storage
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
    if ~exist([direc, lmethod,'_',num2str(run_num),'.mat'],'file')
        if run_num == 1
            disp('Ready to generate new results? Could be time and storage consuming...If so, press Continue')
            keyboard
        end
        temp_v0 = P.v0;
        % find malaria EE
        P.v0 = 0;
        Malaria_parameters_transform_vac;
        [SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');
        P.v0 = temp_v0;
        Malaria_parameters_transform_vac;
        [~, SH_solu, EH_solu, DH_solu, AH_solu, VH_solu, UH_solu, SM_solu, EM_solu, IM_solu, ~, ~, ~, Ctot_solu, MH_solu] = ...
            age_structured_Malaria_vac(P.da,P.na, 0, P.tfinal,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
        SH = SH_solu(:,time_points); EH = EH_solu(:,time_points); DH = DH_solu(:,time_points); AH = AH_solu(:,time_points); MH = MH_solu(:,time_points);
        VH = VH_solu(:,time_points); UH = UH_solu(:,time_points);
        Ctot = Ctot_solu(:,time_points);
        SM = SM_solu(:,time_points); EM = EM_solu(:,time_points); IM = IM_solu(:,time_points);
        NM = SM+EM+IM;
        PH = SH+EH+DH+AH+VH+UH;
        AH_tot = trapz(AH,1)*da;
        DH_tot = trapz(DH,1)*da;
        muDH_tot = trapz(P.muD.*DH,1)*da;
        NH = trapz(PH,1)*P.da;
        [bH,~] =  biting_rate(PH,NM);
        EIR_age = bH.*IM./NM*365;
        EIR_tot = trapz(EIR_age.*PH,1)*P.da./NH;
        save([direc,lmethod,'_',num2str(run_num),'.mat'],'AH_tot','DH_tot','muDH_tot','EIR_tot');
    else
        % load quantities - already integrated over age dimension
        load([direc,lmethod,'_',num2str(run_num),'.mat'],'AH_tot','DH_tot','muDH_tot','EIR_tot');
    end   
end

Q_val = NaN(length(time_points),length(lQ));

for iQ = 1:length(lQ)
    switch lQ{iQ}
        case 'EE-D'
            Q_val(:,iQ) = DH_tot;
        case 'EE-A'
            Q_val(:,iQ) = AH_tot;
        case 'EE-DA'
            Q_val(:,iQ) = AH_tot+DH_tot;
        case 'EE-death-rate' % disease-induced mortality rate, diagnostic eqn MH
            Q_val(:,iQ) = muDH_tot;
        case 'EE-EIR'
            Q_val(:,iQ) = EIR_tot; % annual EIR
        otherwise
            keyboard
    end
end


end