function Q_val = QOI_value_SA(lQ,time_points,run_num,lmethod)
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
    if ~exist(['Results/vaccine_no/',lmethod,'_',num2str(run_num),'.mat'],'file')
        [SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');
        [~, SH_solu, EH_solu, DH_solu, AH_solu, VH_solu, UH_solu, SM_solu, EM_solu, IM_solu, ~, ~, ~, Ctot_solu, MH_solu] = ...
            age_structured_Malaria_vac(P.da,P.na, 0, P.tfinal,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
        SH = SH_solu(:,time_points); EH = EH_solu(:,time_points); DH = DH_solu(:,time_points); AH = AH_solu(:,time_points); MH = MH_solu(:,time_points);
        VH = VH_solu(:,time_points); UH = UH_solu(:,time_points);
        % Cm = Cm_solu(:,time_points); Cac = Cac_solu(:,time_points); Cv = Cv_solu(:,time_points);
        Ctot = Ctot_solu(:,time_points);
        SM = SM_solu(:,time_points); EM = EM_solu(:,time_points); IM = IM_solu(:,time_points);
        save(['Results/vaccine_no/',lmethod,'_',num2str(run_num),'.mat'],'SH','EH','DH','AH','MH','VH','UH','Ctot','SM','EM','IM');
    else
        load(['Results/vaccine_no/',lmethod,'_',num2str(run_num),'.mat'],'SH','EH','DH','AH','MH','VH','UH','Ctot','SM','EM','IM');
    end
    PH = SH+EH+DH+AH+VH+UH;
    % figure_setups;
    % plot(t/365,trapz(SH_solu,1)*da); hold on;
    % plot(t/365,trapz(EH_solu,1)*da);
    % plot(t/365,trapz(AH_solu,1)*da);
    % plot(t/365,trapz(DH_solu,1)*da);
    % plot(t/365,trapz(VH_solu,1)*da);
    % plot(t/365,trapz(UH_solu,1)*da);
    % plot(t/365,trapz(MH_solu,1)*da,'-.r'); % diagnostic
    % title(['Population size vs time']);
    % grid on; grid minor
    % keyboard
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
            Q_val(:,iQ) = trapz(DH,1)*da;
        case 'EE-D-02-10'
            Q_val(:,iQ) = trapz(DH(ind0210y,:),1)*da;
        case 'EE-D-09-24'
            Q_val(:,iQ) = trapz(DH(ind0924m,:),1)*da;
        case 'EE-DA'
            Q_val(:,iQ) = trapz(DH+AH,1)*da;
        case 'EE-DA-02-10'
            Q_val(:,iQ) = trapz(DH(ind0210y,:)+AH(ind0210y,:),1)*da;
        case 'EE-DA-09-24'
            Q_val(:,iQ) = trapz(DH(ind0924m,:)+AH(ind0924m,:),1)*da;
        case 'EE-D-frac'
            Q_val(:,iQ) = trapz(DH,1)/trapz(DH+AH,1);
        case 'EE-D-frac-02-10'
            Q_val(:,iQ) = trapz(DH(ind0210y,:),1)./trapz(DH(ind0210y,:)+AH(ind0210y,:),1);
        case 'EE-D-frac-09-24'
            Q_val(:,iQ) = trapz(DH(ind0924m,:),1)./trapz(DH(ind0924m,:)+AH(ind0924m,:),1);
        case 'EE-death' % Cumulative disease-induced mortality, diagnostic eqn MH
            Q_val(:,iQ) = trapz(MH,1)*da;
        case 'EE-death-02-10' % Cumulative disease-induced mortality, diagnostic eqn MH
            Q_val(:,iQ) = trapz(MH(ind0210y,:),1)*da;
        case 'EE-death-09-24' % Cumulative disease-induced mortality, diagnostic eqn MH
            Q_val(:,iQ) = trapz(MH(ind0924m,:),1)*da;
        case 'EE-EIR'
            NH = trapz(PH,1)*da;
            NM = SM+EM+IM;
            [bH,~] = biting_rate(NH,NM);
            IM_frac = IM./NM;
            Q_val(:,iQ) = bH.*IM_frac*365; % annual EIR
        case 'EE-Ctot-pp'
            Q_val(:,iQ) = trapz(Ctot,1)./trapz(PH,1);
        case 'EE-Ctot-pp-02-10'
            Q_val(:,iQ) = trapz(Ctot(ind0210y,:),1)./trapz(PH(ind0210y,:),1);
        case 'EE-Ctot-pp-09-24'
            Q_val(:,iQ) = trapz(Ctot(ind0924m,:),1)./trapz(PH(ind0924m,:),1);
        case 'DALY' % instantaneous (time series) DALY based on daily cases
            [DALY,~,~] = DALY_cal(SH, EH, DH, AH, VH, UH, SM, EM, IM, Ctot);
            Q_val(:,iQ) = DALY;
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