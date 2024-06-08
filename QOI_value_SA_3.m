function Q_val = QOI_value_SA_3(lQ,time_points,run_num,lmethod,direc)
% compress saving output, trapz on the age dimension to reduce storage
% for cluster runs, saving stats for age groups
global P

da = P.da;
a = P.a;

for iQ = 1:length(lQ)
    if strcmp(lQ{iQ}(1:2),'EE')
        flag_EE = 1;
        break
    end
end

flag_run = 0;

ind0210y = age_range_ind(a,2,10);

if flag_EE
    if ~exist([direc, lmethod,'_',num2str(run_num),'.mat'],'file')
        if run_num == 1
            disp('Ready to generate new results? Could be time and storage consuming...If so, press Continue')
            %keyboard
        end
        flag_run = 1;
        done = NaN;
        % save dummy file for placeholder
        save([direc,lmethod,'_',num2str(run_num),'.mat'],'done');
    else % file exists, can move on
        flag_run = 0;
    end

    if flag_run == 1
        disp(run_num)
        temp_v0 = P.v0;
        % find malaria EE
        P.v0 = 0;
        Malaria_parameters_transform_vac;
        [SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');
        P.v0 = temp_v0;
        Malaria_parameters_transform_vac;
        [~, SH_solu, EH_solu, DH_solu, AH_solu, VH_solu, UH_solu, SM_solu, EM_solu, IM_solu, ~, ~, ~, Ctot_solu, MH_solu] = ...
            age_structured_Malaria_vac(P.da,P.na, 0, P.tfinal,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
        SH = SH_solu(:,time_points); EH = EH_solu(:,time_points); 
        DH = DH_solu(:,time_points); AH = AH_solu(:,time_points); MH = MH_solu(:,time_points);
        VH = VH_solu(:,time_points); UH = UH_solu(:,time_points);
        % Ctot = Ctot_solu(:,time_points);
        % SM = SM_solu(:,time_points); EM = EM_solu(:,time_points); IM = IM_solu(:,time_points);
        % NM = SM+EM+IM;
        % PH = SH+EH+DH+AH+VH+UH;
        death_2_10 = trapz(MH(ind0210y,:),1)*da;
        death_10plus = trapz(MH(ind0210y(end)+1:end,:),1)*da;        
        DA_2_10 = trapz(DH(ind0210y,:)+AH(ind0210y,:),1)*da;
        DA_10plus(:,iQ) = trapz(DH(ind0210y(end)+1:end,:)+AH(ind0210y(end),:),1)*da;
        save([direc,lmethod,'_',num2str(run_num),'.mat'],'death_2_10','death_10plus','DA_2_10','DA_10plus');
    else % if file already exists
        if ~isempty(who('-file', [direc,lmethod,'_',num2str(run_num),'.mat'], 'death_2_10'))
            % load quantities - already integrated over age dimension
            load([direc,lmethod,'_',num2str(run_num),'.mat'],'death_2_10','death_10plus','DA_2_10','DA_10plus');
        else
            % file started but interrupted
            Q_val = NaN;
            return
        end
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
        case 'EE-DA-02-10'
            Q_val(:,iQ) = DA_2_10;
        case 'EE-DA-10+'
            Q_val(:,iQ) = DA_10plus;
        case 'EE-death-rate' % disease-induced mortality rate, diagnostic eqn MH
            Q_val(:,iQ) = muDH_tot;
        case 'EE-death-02-10' % Cumulative disease-induced mortality, diagnostic eqn MH
            Q_val(:,iQ) = death_2_10;
        case 'EE-death-10+' % Cumulative disease-induced mortality, diagnostic eqn MH
            Q_val(:,iQ) = death_10plus;    
        case 'EE-EIR'
            Q_val(:,iQ) = EIR_tot; % annual EIR
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