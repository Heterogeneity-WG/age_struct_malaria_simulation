function [DALY,YLL,YLD] = DALY_cal(SH, EH, DH, AH, VH, UH, SM, EM, IM, Ctot)

global P

Malaria_parameters_DALY;

% use time series solutions
PH = SH+EH+DH+AH+VH+UH;
NH = trapz(PH,1)*P.da;
NM = SM+EM+IM;

[bH,~] = biting_rate(NH,NM);
lamH = FOI_H(bH,IM,NM);
rho = sigmoid_prob(Ctot./PH, 'rho'); % prob. of severely infected, EH -> DH
psi = sigmoid_prob(Ctot./PH, 'psi'); % prob. AH -> DH
DH_cases = rho.*P.h.*EH+psi.*lamH.*AH; % new DH counts in time (age, time)
MH_cases = P.muD.*DH;
MH_cases_cutoff = MH_cases.*(P.a<P.LifeEx); % exclude people who died at the standard age or older
age_of_death_byage = P.a.*MH_cases_cutoff;
life_expectancy = P.LifeEx*MH_cases_cutoff;

YLL_byage = (life_expectancy-age_of_death_byage)/365; % time series 
YLL = trapz(YLL_byage,1)*P.da;

YLD_byage = P.DW.*DH_cases*P.L/365;
YLD = trapz(YLD_byage,1)*P.da; 

DALY = YLL+YLD;
% keyboard
% plot(P.a, YLD_byage(:,1),P.a, YLL_byage(:,1))
% legend('YLD','YLL')
end