%% objective function for model calibration
function err = fun_Filipe_dynamic(x)
global F
global P
% F(age(years), EIR) 
% =  immunity (~ Ctot) from RB's paper, Fig 5
% (~ rho) from Filipe's paper

P.phis2 = x(1);
P.phir2 = x(2);
P.rhos2 = x(3);
P.rhor2 = x(4); 
P.psis2 = x(5);
P.psir2 = x(6);

Malaria_parameters_transform;
nsamp = 20;
[~,ind1] = min(abs(P.a-0.5*365)); % start from 0.5 years old
[~,ind2] = min(abs(P.a-10*365)); % end at 10 years old
ind_a = round(linspace(ind1,ind2,nsamp)');

[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('init');
P.betaM = 0.025; % low EIR region ~ 25
Malaria_parameters_transform;
[SH, EH, DH, AH, ~, ~, SM, EM, IM, ~, ~, ~, Ctot, ~] = age_structured_Malaria_vac(P.da,P.na,P.tfinal,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
EIR = fit_EIR(SH,EH,DH,AH,SM,EM,IM);
if EIR(end)<0.5; keyboard; end
PH = SH+EH+DH+AH;  % no vaccine
x1 = P.a(ind_a)/365;
y1 = EIR(end); % aEIR
[X1,Y1] = ndgrid(x1,y1);
Z1 = Ctot(ind_a,end)./PH(ind_a,end); % final Ctot at EE

P.betaM = 0.05; % med EIR region ~ 50
Malaria_parameters_transform;
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('init');
[SH, EH, DH, AH, ~, ~, SM, EM, IM, ~, ~, ~, Ctot, ~] = age_structured_Malaria_vac(P.da,P.na,P.tfinal,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
EIR = fit_EIR(SH,EH,DH,AH,SM,EM,IM);
if EIR(end)<1; keyboard; end
PH = SH+EH+DH+AH;
x2 = P.a(ind_a)/365;
y2 = EIR(end); % aEIR
[X2,Y2] = ndgrid(x2,y2);
Z2 = Ctot(ind_a,end)./PH(ind_a,end); % final Ctot at EE

P.betaM = 0.5; % high EIR region ~ 90
Malaria_parameters_transform;
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('init');
[SH, EH, DH, AH, ~, ~, SM, EM, IM, ~, ~, ~, Ctot, ~] = age_structured_Malaria_vac(P.da,P.na,P.tfinal,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
EIR = fit_EIR(SH,EH,DH,AH,SM,EM,IM);
if EIR(end)<1; keyboard; end
PH = SH+EH+DH+AH;
x3 = P.a(ind_a)/365;
y3 = EIR(end); % aEIR
[X3,Y3] = ndgrid(x3,y3);
Z3 = Ctot(ind_a,end)./PH(ind_a,end); % final Ctot at EE

Z1_data = F(X1,Y1); % rho from data
Z2_data = F(X2,Y2); % rho from data
Z3_data = F(X3,Y3); % rho from data

Z1_samp = sigmoid_prob(Z1, 'rho'); % rho from samples
Z2_samp = sigmoid_prob(Z2, 'rho'); % rho from samples
Z3_samp = sigmoid_prob(Z3, 'rho'); % rho from samples

res = abs([(Z1_data(:)-Z1_samp(:))./Z1_data(:); (Z2_data(:)-Z2_samp(:))./Z2_data(:); (Z3_data(:)-Z3_samp(:))./Z3_data(:)]);
w = ones(size(res));%./(res+eps);
err = sum(w.*(res.^2));

% [y1, y2, y3]
end

function EIR = fit_EIR(SH,EH,DH,AH,SM,EM,IM)
global P
NH = trapz(SH+EH+DH+AH)*P.da;
NM = SM+EM+IM;
[bH,~] = biting_rate(NH,NM);
IM_frac = IM./NM;
EIR = bH.*IM_frac*365; % annual EIR
end

% fun = @(x) f(x,10,SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0);
% P.betaM = fminbnd(fun,0,1);  
% P.betaM_low = P.betaM;

% function err = f(x,EIR_target,SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0)
% global P
% 
% P.betaM = x;
% Malaria_parameters_transform;
% [SH, EH, DH, AH, SM, EM, IM, ~, ~, ~] = age_structured_Malaria(P.da,P.na,P.tfinal,SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0);
% EIR = fit_EIR(SH,EH,DH,AH,SM, EM, IM);
% err = abs(EIR(end)-EIR_target);
% end