function Q_val = QOI_save(lQ,soln)
global P

da = P.da;
a = P.a;
h =  P.h;

ind0210y = age_range_ind(a,2,10);
ind0924m = age_range_ind(a,9/12,24/12);

% Unpack full solution
SH=soln.SH; EH=soln.EH; DH=soln.DH; AH=soln.AH; VH=soln.VH; UH=soln.UH;
SM=soln.SM; EM=soln.EM; IM=soln.IM;
Cm=soln.Cm; Cac=soln.Cac; Cv=soln.Cv; Ctot=soln.Ctot;
MH=soln.MH;

PH = SH + EH + DH + AH + VH + UH;
NH = trapz(PH,1)*da;
NM = SM + EM + IM;

% extract the solution at the last time step
SHend=SH(:,end); EHend=EH(:,end); DHend=DH(:,end); AHend=AH(:,end); VHend=VH(:,end); UHend=UH(:,end);
SMend=SM(:,end); EMend=EM(:,end); IMend=IM(:,end);
Cmend=Cm(:,end); Cacend=Cac(:,end); Cvend=Cv(:,end); Ctotend=Ctot(:,end);
MHend=MH(:,end);

PHend = SHend+EHend+DHend+AHend+VHend+UHend;

for iQ = 1:length(lQ)
    switch lQ{iQ}    
        case {'EE-D','year5-D'}
            Q_val = trapz(DHend)*da;
        case {'year5-D-process'}
            Q_val = trapz(DH,1)*da;
        case {'year5-A-process'}
            Q_val = trapz(AH,1)*da;
        case 'EE-D-02-10'
            Q_val = trapz(DHend(ind0210y))*da;
        case 'EE-D-09-24'
            Q_val = trapz(DHend(ind0924m))*da;
        case 'EE-DA'
            Q_val= trapz(DHend+AHend)*da;    
        case 'EE-DA-02-10'            
            Q_val = trapz(DHend(ind0210y)+AHend(ind0210y))*da;  
        case 'EE-DA-09-24'
            Q_val = trapz(DHend(ind0924m)+AHend(ind0924m))*da;  
        case {'EE-D-frac','year5-D-frac'}
            Q_val = trapz(DHend)/trapz(DHend+AHend);
        case {'year5-D-frac-process'}
            Q_val = trapz(DH,1)./trapz(DH+AH,1);
        case 'EE-D-frac-02-10'
            Q_val = trapz(DHend(ind0210y))/trapz(DHend(ind0210y)+AHend(ind0210y));
        case 'EE-D-frac-09-24'
            Q_val = trapz(DHend(ind0924m))/trapz(DHend(ind0924m)+AHend(ind0924m));
        case 'EE-A-frac'
            Q_val = trapz(AHend)/trapz(DHend+AHend);
        case 'EE-A-frac-02-10'
            Q_val = trapz(AHend(ind0210y))/trapz(DHend(ind0210y)+AHend(ind0210y));
        case 'EE-A-frac-09-24'
            Q_val = trapz(AHend(ind0924m))/trapz(DHend(ind0924m)+AHend(ind0924m));
        case 'EE-death' % Cumulative disease-induced mortality, diagnostic eqn MH
            Q_val = trapz(MHend)*da;
        case 'EE-death-02-10' % Cumulative disease-induced mortality, diagnostic eqn MH
            Q_val = trapz(MHend(ind0210y))*da;
        case 'EE-death-09-24' % Cumulative disease-induced mortality, diagnostic eqn MH
            Q_val = trapz(MHend(ind0924m))*da;
        case 'EE-EIR'
            NHend = trapz(PHend)*da;
            NMend = SMend+EMend+IMend;
            [bH,~] = biting_rate(NHend,NMend);
            IM_frac = IMend./NMend;
            Q_val = bH.*IM_frac*365; % annual EIR
        case 'EE-Ctot-pp'
            Q_val = trapz(Ctotend)/trapz(PHend);
        case 'EE-Ctot-pp-02-10'
            Q_val = trapz(Ctotend(ind0210y))/trapz(PHend(ind0210y));
        case 'EE-Ctot-pp-09-24'
            Q_val = trapz(Ctotend(ind0924m))/trapz(PHend(ind0924m));
        case 'DH_Incidence' % Incidence of epidsodes (like cases)
            [bH,~] = biting_rate(NH,NM);
            lamH = FOI_H(bH,IM,NM);
            
            rho = sigmoid_prob(Ctot./PH, 'rho'); % prob. of severely infected, EH -> DH
            psi = sigmoid_prob(Ctot./PH, 'psi'); % prob. AH -> DH
            temp1 = rho.*h.*EH+psi.*lamH.*AH; % incidence of DH
            Q_val = trapz(temp1,1)*da;
        case 'Vaccines'
            vacc_sterile = trapz(P.v*(1-P.z).*SH,1)*da*365*P.NN/1000;
            vacc_blood = trapz(P.v*P.z.*SH,1)*da*365*P.NN/1000;
            keyboard
            Q_val.Vaccines_sterile = vacc_sterile;
            Q_val.Vaccines_blood = vacc_blood;
        case 'DALY' % instantaneous (time series) DALY based on daily cases
            [DALY,~,~] = DALY_cal(SH, EH, DH, AH, VH, UH, SM, EM, IM, Ctot); 
            
            Q_val = DALY; 
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