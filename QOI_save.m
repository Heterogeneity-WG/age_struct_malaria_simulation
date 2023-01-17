function Q_val = QOI_save(lQ,soln)
global P

da = P.da;
a = P.a;
h =  P.h;

ind0210y = age_range_ind(a,2,10);
ind0924m = age_range_ind(a,9/12,24/12);

% Unpack full solution
SH=soln.SH; EH=soln.EH; DH=soln.DH; AH=soln.AH(:,end); VH=soln.VH; UH=soln.UH;
SM=soln.SM; EM=soln.EM; IM=soln.IM;
Cm=soln.Cm; Cac=soln.Cac; Cv=soln.Cv; Ctot=soln.Ctot;
MH=soln.MH;

% Unpack final time solution
SHend=SH(:,end); EHend=EH(:,end); DHend=DH(:,end); AHend=AH(:,end); VHend=VH(:,end); UHend=UH(:,end);
SMend=SM(:,end); EMend=EM(:,end); IMend=IM(:,end);
Cmend=Cm(:,end); Cacend=Cac(:,end); Cvend=Cv(:,end); Ctotend=Ctot(:,end);
MHend=MH(:,end);

PH = SH + EH + DH + AH + VH + UH;
NH = trapz(PH,1)*da;
NM = SM + EM + IM;

for iQ = 1:length(lQ)
    switch lQ{iQ}    
        case 'EE-D'
            Q_val.EE_D = trapz(DHend)*da;
        case 'EE-D-02-10'
            Q_val.EE_D_02_10 = trapz(DHend(ind0210y))*da;
        case 'EE-D-09-24'
            Q_val.EE_D_09_24 = trapz(DHend(ind0924m))*da;
        case 'EE-DA'
            Q_val.EE_DA = trapz(DHend+AHend)*da;    
        case 'EE-DA-02-10'            
            Q_val.DA_02_10 = trapz(DHend(ind0210y)+AHend(ind0210y))*da;  
        case 'EE-DA-09-24'
            Q_val.EE_DA_09_24 = trapz(DHend(ind0924m)+AHend(ind0924m))*da;  
        case 'EE-D-frac'
            Q_val.EE_D_frac = trapz(DHend)/trapz(DHend+AHend);
        case 'EE-D-frac-02-10'
            Q_val.EE_D_frac_02_10 = trapz(DHend(ind0210y))/trapz(DHend(ind0210y)+AHend(ind0210y));
        case 'EE-D-frac-09-24'
            Q_val.EE_D_frac_09_24 = trapz(DHend(ind0924m))/trapz(DHend(ind0924m)+AHend(ind0924m));
        case 'EE-A-frac'
            Q_val.EE_A_frac = trapz(AHend)/trapz(DHend+AHend);
        case 'EE-A-frac-02-10'
            Q_val.EE_A_frac_02_10 = trapz(AHend(ind0210y))/trapz(DHend(ind0210y)+AHend(ind0210y));
        case 'EE-A-frac-09-24'
            Q_val.EE_A_frac_09_24 = trapz(AHend(ind0924m))/trapz(DHend(ind0924m)+AHend(ind0924m));
        case 'EE-death' % Cumulative disease-induced mortality, diagnostic eqn MH
            Q_val.EE_death = trapz(MHend)*da;
        case 'EE-death-02-10' % Cumulative disease-induced mortality, diagnostic eqn MH
            Q_val.EE_death_02_10 = trapz(MHend(ind0210y))*da;
        case 'EE-death-09-24' % Cumulative disease-induced mortality, diagnostic eqn MH
            Q_val.EE_death_09_24 = trapz(MHend(ind0924m))*da;
        case 'EE-EIR'
            NH = trapz(PH)*da;
            NM = SMend+EMend+IMend;
            [bH,~] = biting_rate(NH,NM);
            IM_frac = IMend./NM;
            Q_val.EE_EIR = bH.*IM_frac*365; % annual EIR
        case 'EE-Ctot-pp'
            Q_val.EE_Ctot_pp = trapz(Ctotend)/trapz(PH);
        case 'EE-Ctot-pp-02-10'
            Q_val.EE_Ctot_pp_02_10 = trapz(Ctotend(ind0210y))/trapz(PH(ind0210y));
        case 'EE-Ctot-pp-09-24'
            Q_val.EE_Ctot_pp_09_24 = trapz(Ctotend(ind0924m))/trapz(PH(ind0924m));
        case 'DH_Incidence' % Incidence of epidsodes (like cases)
            [bH,~] = biting_rate(NH,NM);
            lamH = FOI_H(bH,IM,NM);
            
            rho = sigmoid_prob(Ctot./PH, 'rho'); % prob. of severely infected, EH -> DH
            psi = sigmoid_prob(Ctot./PH, 'psi'); % prob. AH -> DH
            rho(PH==0)=0; psi(PH==0)=0;
            temp1 = rho.*h.*EH+psi.*lamH.*AH; % incidence of DH
            Q_val.DH_Incidence = trapz(temp1,1)*da;
        case 'Vaccines'
            vacc_sterile = trapz(P.v*(1-P.z).*SH,1)*da*365*P.NN/1000;
            vacc_blood = trapz(P.v*P.z.*SH,1)*da*365*P.NN/1000;

            Q_val.Vaccines_sterile = vacc_sterile;
            Q_val.Vaccines_blood = vacc_blood;
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