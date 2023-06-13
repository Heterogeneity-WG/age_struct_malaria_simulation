%% objective function for model calibration
function err = fun_Filipe_dynamic2(x)
global F
global P
% F(age(years), EIR) 
% =  immunity (~ Ctot) from RB's paper, Fig 5
% (~ rho) from Filipe's paper

% [phi_s phi_r rho_s rho_s psi_r psi_r phif0 phif1 rhof0 rhof1 psif0 psif1]
P.phis2 = x(1);
P.phir2 = x(2);
P.rhos2 = x(3);
P.rhor2 = x(4); 
P.psis2 = x(5);
P.psir2 = x(6);

Malaria_parameters_transform;
nsamp = 30;
[~,ind1] = min(abs(P.a-0.3*365)); 
[~,ind0] = min(abs(P.a-10*365)); 
[~,ind2] = min(abs(P.a-20*365)); 
ind_a = [round(linspace(ind1,ind0,20)'); round(linspace(ind0+1,ind2,nsamp-20)')];

res_list = [];
EIR_list = [1 10:20:90 100];
for iEIR = 1:length(EIR_list)
    EIR = EIR_list(iEIR);
    x = bisection(@(beta) Cost_func(beta,EIR), 0, 1, 10, 10^-1, 10^-3, 0);
    P.betaM = x; 
    Malaria_parameters_transform;
    [SH, EH, DH, AH, ~, ~, SM, EM, IM, ~, ~, ~, Ctot, ~] = age_structured_Malaria_IC_vac('EE_reset');
    PH = SH+EH+DH+AH;
    NM = SM+EM+IM;
    [bH,~] = biting_rate(PH,NM);
    EIR = bH.*IM./NM*365; % EIR matrix
    EIR_tot = trapz(EIR.*PH)/trapz(PH); % EIR sum over age, at final time
    EIR_ss = EIR_tot;
    x = P.a(ind_a)/365;
    y = EIR_ss; % aEIR
    Z = Ctot(ind_a,end)./PH(ind_a,end); % final Ctot at EE
    Z_samp = sigmoid_prob(Z, 'rho'); % rho from samples
    [X,Y] = ndgrid(x,y);
    Z_data = F(X,Y); % rho from data
    res_list = [res_list; abs(Z_data(:)-Z_samp(:))];
    EIR_list(iEIR) = EIR_ss;
end

w = ones(size(res_list));%./(res+eps);
err = sum(w.*(res_list.^2));
% disp(round([EIR_list'; err]',3))
end

function y = Cost_func(beta,EIR_target)
global P
P.betaM = beta;
Malaria_parameters_transform;
[SH, EH, DH, AH, ~, ~, SM, EM, IM, ~, ~, ~, ~, ~] = age_structured_Malaria_IC_vac('EE_reset');
PH = SH+EH+DH+AH;
NM = SM+EM+IM;
[bH,~] = biting_rate(PH,NM);
EIR = bH.*IM./NM*365; % EIR matrix
EIR_ss = trapz(EIR.*PH)/trapz(PH); % EIR sum over age, at final time
y = EIR_ss - EIR_target;
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