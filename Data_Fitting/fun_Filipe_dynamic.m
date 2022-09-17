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
[~,ind1] = min(abs(P.a-0.3*365)); % start from 0.5 years old
[~,ind0] = min(abs(P.a-3*365)); % start from 0.5 years old
[~,ind2] = min(abs(P.a-10*365)); % end at 10 years old
ind_a = [round(linspace(ind1,ind0,15)'); round(linspace(ind0+1,ind2,nsamp-15)')];
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('init');

betaM_list = [0.001:0.002:0.03, 0.04:0.01:0.3, 0.5, 1];
res_list = [];
EIR_list = [];
% low EIR region ~ 25
for ibeta = 1:length(betaM_list)
    P.betaM = betaM_list(ibeta); 
    Malaria_parameters_transform;
    [SH, EH, DH, AH, ~, ~, SM, EM, IM, ~, ~, ~, Ctot, ~] = age_structured_Malaria_vac(P.da,P.na,P.tfinal,SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
    EIR = fit_EIR(SH,EH,DH,AH,SM,EM,IM);
    EIR_ss = EIR(end);
    if EIR_ss<0.5; keyboard; end
    PH = SH+EH+DH+AH;  % no vaccine
    x = P.a(ind_a)/365;
    y = EIR(end); % aEIR
    Z = Ctot(ind_a,end)./PH(ind_a,end); % final Ctot at EE
    Z_samp = sigmoid_prob(Z, 'rho'); % rho from samples
    [X,Y] = ndgrid(x,y);
    Z_data = F(X,Y); % rho from data
    res_list = [res_list; abs(Z_data(:)-Z_samp(:))];
    EIR_list = [EIR_list; EIR_ss];
end

w = ones(size(res_list));%./(res+eps);
err = sum(w.*(res_list.^2));
disp(round([EIR_list', err]',3))
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