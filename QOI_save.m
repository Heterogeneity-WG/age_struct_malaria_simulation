function Q_val = QOI_save(lQ,soln)
global P

da = P.da;
a = P.a;

ind0210y = age_range_ind(a,2,10);
ind0924m = age_range_ind(a,9/12,24/12);

SH=soln.SH(:,end); EH=soln.EH(:,end); DH=soln.DH(:,end); AH=soln.AH(:,end); VH=soln.VH(:,end); UH=soln.UH(:,end);
SM=soln.SM(:,end); EM=soln.EM(:,end); IM=soln.IM(:,end);
Cm=soln.Cm(:,end); Cac=soln.Cac(:,end); Cv=soln.Cv(:,end); Ctot=soln.Ctot(:,end);
MH=soln.MH(:,end);

for iQ = 1:length(lQ)
    switch lQ{iQ}    
        case 'EE-D'
            Q_val.EE_D = trapz(DH)*da;
        case 'EE-D-02-10'
            Q_val.EE_D_02_10 = trapz(DH(ind0210y))*da;
        case 'EE-D-09-24'
            Q_val.EE_D_09_24 = trapz(DH(ind0924m))*da;
        case 'EE-DA'
            Q_val.EE_DA = trapz(DH+AH)*da;    
        case 'EE-DA-02-10'            
            Q_val.DA_02_10 = trapz(DH(ind0210y)+AH(ind0210y))*da;  
        case 'EE-DA-09-24'
            Q_val.EE_DA_09_24 = trapz(DH(ind0924m)+AH(ind0924m))*da;  
        case 'EE-D-frac'
            Q_val.EE_D_frac = trapz(DH)/trapz(DH+AH);
        case 'EE-D-frac-02-10'
            Q_val.D_frac_02_10 = trapz(DH(ind0210y))/trapz(DH(ind0210y)+AH(ind0210y));
        case 'EE-D-frac-09-24'
            Q_val.D_frac_09_24 = trapz(DH(ind0924m))/trapz(DH(ind0924m)+AH(ind0924m));
        case 'EE-death' % Cumulative disease-induced mortality, diagnostic eqn MH
            Q_val.EE_death = trapz(MH)*da;
        case 'EE-death-02-10' % Cumulative disease-induced mortality, diagnostic eqn MH
            Q_val.EE_death_02_10 = trapz(MH(ind0210y))*da;
        case 'EE-death-09-24' % Cumulative disease-induced mortality, diagnostic eqn MH
            Q_val.EE_death_09_24 = trapz(MH(ind0924m))*da;
        case 'EE-EIR'
            NH = trapz(PH)*da;
            NM = SM+EM+IM;
            [bH,~] = biting_rate(NH,NM);
            IM_frac = IM./NM;
            Q_val.EE_EIR = bH.*IM_frac*365; % annual EIR
        case 'EE-Ctot-pp'
            Q_val.EE_Ctot_pp = trapz(Ctot)/trapz(PH);
        case 'EE-Ctot-pp-02-10'
            Q_val.EE_Ctot_pp_02_10 = trapz(Ctot(ind0210y))/trapz(PH(ind0210y));
        case 'EE-Ctot-pp-09-24'
            Q_val.EE_Ctot_pp_09_24 = trapz(Ctot(ind0924m))/trapz(PH(ind0924m));
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