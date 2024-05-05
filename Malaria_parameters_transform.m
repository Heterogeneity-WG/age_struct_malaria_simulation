function Malaria_parameters_transform
global P

a = P.a;

%% age-dependent biting

P.zeta_fun = @(a) 1-0.85*exp(-a/8/365); %ones(size(a));
P.zeta = P.zeta_fun(a);

%% seasonality function
if sum(contains(fieldnames(P),'ss_S0')) && sum(contains(fieldnames(P),'ss_v')) % check if seasonality is turned on
    P.ss_S = @(t) P.ss_S0*(P.ss_c+P.ss_v*(1-P.ss_c)*((1+cos(2*pi*((t+P.ss_t0)/365-P.ss_u1)))./2).^P.ss_k1+...
        (1-P.ss_v)*(1-P.ss_c)*((1+cos(2*pi*((t+P.ss_t0)/365-P.ss_u2)))./2).^P.ss_k2);
    P.gM_fun = @(t) 0.5*P.NN*P.ss_S(t);
else
    P.gM_fun = @(t) 0.5*P.NN;
end
% P.ss_S = @(t) P.ss_S0*(P.ss_c+P.ss_v*(1-P.ss_c)*((1+cos(2*pi*((t+P.ss_t0)/365-P.ss_u1)))./2).^P.ss_k1+...
%     (1-P.ss_v)*(1-P.ss_c)*((1+cos(2*pi*((t+P.ss_t0)/365-P.ss_u2)))./2).^P.ss_k2);
% P.gM_fun = @(t) 0.5*P.NN*P.ss_S(t); % recruitment rate of mosquitoes;

%%
P.c2 = P.c1; % weight for maternal immunity
P.c3 = P.c1; % weight for vaccine-derived immunity

P.rho = sigmoid_prob(zeros(size(a)), 'rho');
P.phi = sigmoid_prob(zeros(size(a)), 'phi');
P.psi = sigmoid_prob(zeros(size(a)), 'psi');

%% natural mortality functions
muH =  P.b0 + P.b1*exp(-P.b2*a/365) + P.b3*exp(P.b4*a/365); % natural human mortality rate
muH_fun = @(age) P.b0 + P.b1*exp(-P.b2*age/365) + P.b3*exp(P.b4*age/365);
muH = muH/365;
muH_int_fun = @(age) (age./365).*P.b0 + (P.b1./P.b2).*(1-exp(-P.b2.*age./365)) + (P.b3./P.b4).*(-1+exp(P.b4.*age./365));

%% malaria mortality functions
if P.disease_mortality
    a74 = 74*365;
    temp_muD = P.b0D + P.b1D*exp(-P.b2D*a74/365) + P.b3D*exp(P.b4D*a74/365);
    muD =  (P.b0D + P.b1D*exp(-P.b2D*a/365) + P.b3D*exp(P.b4D*a/365)).*(a/365<=74)...
        + temp_muD.*(exp(-P.c0D*(a-74*365)./365)).*(a/365>74);
    muD_fun = @(age) (P.b0D + P.b1D*exp(-P.b2D*age/365) + P.b3D*exp(P.b4D*age/365)).*(age/365<=74)...
        + temp_muD.*(exp(-P.c0D*(age/365-74))).*(age/365>74);
    muD = muD/365;
    %muD_int_fun = @(age) (age./365).*P.b0D + (P.b1D./P.b2D).*(1-exp(-P.b2D.*age./365)) + (P.b3D./P.b4D).*(-1+exp(P.b4D.*age./365));
    % not correct, commented out for now
else
    muD = 0*ones(size(a));
    muD_fun = @(age) 0*ones(size(age));
end

%% fertility rate
gH_fun = @(age) (2.*P.cc.*normpdf((age./365-P.zz)./P.ww).*normcdf(P.alpha.*(age./365-P.zz)./P.ww)./P.ww)./365/2;
gH =  gH_fun(a); % human fertility rate

% approximation of theta at DFE - needed for analytical purpose: DFE, R0, bifurcation
% pi_fun = @(x) P.w+P.e.*v_fun(x);
% pi_int_a = intf(pi_fun,P.a);
% pi_int_fun = @(x) interp1(P.a,pi_int_a,x);
% exp_pi_int_a = intf(@(x) P.w.*exp(pi_int_fun(x)),P.a);
% exp_pi_int = @(x) interp1(P.a,exp_pi_int_a,x);
% theta_fun = @(x) exp(-pi_int_fun(x)).*(1+exp_pi_int(x));
% theta0 = theta_fun(P.a);

P.muH = muH;
P.muH_fun = muH_fun;
P.muD = muD;
P.muD_fun = muD_fun;
P.gH = gH;
% P.theta_fun = theta_fun;
% P.theta = theta0;
P.gH_fun = gH_fun;
P.muH_int_fun = muH_int_fun;
P.muH_int = muH_int_fun(a);
%P.muD_int_fun = muD_int_fun;
%P.muD_int = muD_int_fun(a);

% Update the fertility and stable age dist. if balanced option is selected
if P.balance_fertility == 1
    balance_fertility;
end

if P.balance_mortality == 1
    balance_mortality;
end

find_stable_age; % when mu_D>0, this gives the age distribution when no disease presents in the population, only demographic effects

%% plot calibrated fertility & death
% figure_setups; hold on
% grid on
% plot(a/365, P.muH, a/365, P.muD, a/365, P.gH)
% xlabel('Age (years)')
% ylabel('Daily rates')
% legend('$\mu_H(\alpha)$','$\mu_D(\alpha)$','Fertility rate $g_H(\alpha)$')
% keyboard
% %% plot stable age distribution PH
% figure_setups;
% plot(a/365,P.PH_stable,'-k');
% %title(['Stable age distribution']);
% legend('$P_H(\alpha)$')
% xlabel('Age (years)');
% ylabel('Population density')
% grid on
% axis([0 age_max/365 0 max(P.PH_stable)]);
%
% %% flip
% figure_setups;
% plot(P.PH_stable,a/365,'-k');
% title(['Age Distribution']);
% ylabel('age (years)');
% title('Kenya')
% grid on
% axis([0 max(P.PH_stable) 0 age_max/365 ]);
% % keyboard
end

% function intfx = intf(f,xs)
% ns = length(xs);
% fx = f(xs);
% e = ones(ns,ns+1)/2;
% e(2:end,1:end-1)=1;e(1,end)=0;
% A = spdiags(e,-ns:0,ns,ns);
% intfx= A*fx;
% end

