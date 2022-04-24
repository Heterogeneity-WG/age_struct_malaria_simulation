function [R0,RHM,RMH] = R0_cal_muD()
global P

alphamax = P.age_max; %% Inf
muH_int = P.muH_int_fun;
muD_int = P.muD_int_fun;
[~,~,~,~,~,~,~,~, CH] = steady_state_vac('DFE','handle'); % steady state for CH - total immunity (function handle)
CH = @(a) CH(a)./P.PH_stable_fun(a);
sigmoid_rho = sigmoid_prob_fun('rho'); % return a function handle
sigmoid_phi = sigmoid_prob_fun('phi'); % return a function handle
rho = @(a) sigmoid_rho(CH(a));
phi = @(a) sigmoid_phi(CH(a));
%% For D part ----- integral (alpha, a, s)
D_alpha_a_x = @(alpha,a,x) exp(-muH_int(alpha)-muD_int(alpha)).*exp(muD_int(a)).*rho(a).*P.h.*exp(-(P.rD.*alpha+(P.h-P.rD).*a+(P.p_hat-P.h).*x));
amax = @(alpha) alpha;
xmax = @(alpha,a) a;
D_int = integral3(D_alpha_a_x,0,alphamax,0,amax,0,xmax);
%%  For A part --- integral on (alpha, a, x) + integral on (alpha, a, x, s)...
A1_alpha_a_x = @(alpha,a,x) exp(-muH_int(alpha)).*(1-rho(a)).*P.h.*exp(-(P.rA.*alpha+(P.h-P.rA).*a+(P.p_hat-P.h).*x));
A2_alpha_a_x_s = @(alpha,a,x,s) exp(-muH_int(alpha)-muD_int(a)).*exp(muD_int(x)).*P.rD.*(1-phi(a)).*rho(x).*P.h...
    .*exp(-P.rA.*(alpha-a)-P.rD.*(a-x)-P.h.*(x-s)-P.p_hat.*s);
amax = @(alpha) alpha;
xmax = @(alpha,a) a;
smax = @(alpha,a,x) x;
A1_int = integral3(A1_alpha_a_x,0,alphamax,0,amax,0,xmax);
A2_int = integralN(A2_alpha_a_x_s,0,alphamax,0,amax,0,xmax,0,smax);
A_int = A1_int+A2_int;
%% calculate R0
[bH,bM] = biting_rate(1,P.gM/P.muM);  % assume NH=1; NH(end) for numerical simulation is > 1
K = 1/integral(@(a) exp(-muH_int(a)-P.p_hat.*a), 0, alphamax);
RHM = bH*K*(P.betaD*D_int+P.betaA*A_int);
RMH = bM*P.betaM*P.sigma/(P.sigma+P.muM)./P.muM;
R0 = RHM*RMH;
R0 = sqrt(R0);

end



