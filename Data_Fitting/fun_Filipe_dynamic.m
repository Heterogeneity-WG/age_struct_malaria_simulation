%% objective function for model calibration
function err = fun_Filipe_dynamic(x)
global F
global P
% F(age(years), EIR) 
% =  immunity (~ Ctot) from RB's paper, Fig 5
% (~ rho) from Filipe's paper

P.phi_s_2 = x(1);
P.phi_r_2 = x(2);
P.rho_s_2 = x(3);
P.rho_r_2 = x(4); 
P.psi_s_2 = x(3);
P.psi_r_2 = x(4);

Malaria_parameters_transform;
nsamp = 10;
[~,ind1] = min(abs(P.a-0.5*365)); % start from 0.5 years old
[~,ind2] = min(abs(P.a-10*365)); % end at 10 years old
ind_a = round(linspace(ind1,ind2,nsamp)');

P.betaM = 0.1; % low EIR region ~ 2
Malaria_parameters_transform;
[SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0] = age_structured_Malaria_IC('init');
[SH, EH, DH, AH, SM, EM, IM, ~, ~, Ctot] = age_structured_Malaria(P.da,P.na,P.tfinal,SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0);
EIR = fit_EIR(SH,EH,DH,AH,SM, EM, IM);
if EIR(end)<1; keyboard; end
PH = SH+EH+DH+AH;
x1 = P.a(ind_a)/365;
y1 = EIR(end); % aEIR
[X1,Y1] = ndgrid(x1,y1);
Z1 = Ctot(ind_a,end)./PH(ind_a,end); % final Ctot at EE

P.betaM = 2; % middle EIR region ~ 30
Malaria_parameters_transform;
[SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0] = age_structured_Malaria_IC('init');
[SH, EH, DH, AH, SM, EM, IM, ~, ~, Ctot] = age_structured_Malaria(P.da,P.na,P.tfinal,SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0);
EIR = fit_EIR(SH,EH,DH,AH,SM, EM, IM);
if EIR(end)<1; keyboard; end
PH = SH+EH+DH+AH;
x2 = P.a(ind_a)/365;
y2 = EIR(end); % aEIR
[X2,Y2] = ndgrid(x2,y2);
Z2 = Ctot(ind_a,end)./PH(ind_a,end); % final Ctot at EE

P.betaM = 5; % high EIR region ~ 100
Malaria_parameters_transform;
[SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0] = age_structured_Malaria_IC('init');
[SH, EH, DH, AH, SM, EM, IM, ~, ~, Ctot] = age_structured_Malaria(P.da,P.na,P.tfinal,SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0);
EIR = fit_EIR(SH,EH,DH,AH,SM, EM, IM);
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

res = abs([Z1_data(:)-Z1_samp(:); Z2_data(:)-Z2_samp(:);Z3_data(:)-Z3_samp(:)]);
w = ones(size(res));%./(res+eps);
err = sum(w.*(res.^2));

end

function EIR = fit_EIR(SH,EH,DH,AH,SM,EM,IM)
global P
NH = trapz(SH+EH+DH+AH)*P.da;
NM = SM+EM+IM;
[bH,~] = biting_rate(NH,NM);
IM_frac = IM./NM;
EIR = bH.*IM_frac*365; % annual EIR
end