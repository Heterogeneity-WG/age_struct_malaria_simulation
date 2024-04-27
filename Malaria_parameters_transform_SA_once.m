% common calculations/quantities that aren't impacted by the POIs
% calculate here once for efficiency

P.zeta_fun = @(a) 1-0.85*exp(-a/8/365);
P.zeta = P.zeta_fun(a);
if sum(contains(fieldnames(P),'ss_S0')) % check if seasonality is turned on
    P.ss_S = @(t) P.ss_S0*(P.ss_c+P.ss_v*(1-P.ss_c)*((1+cos(2*pi*((t+P.ss_t0)/365-P.ss_u1)))./2).^P.ss_k1+...
        (1-P.ss_v)*(1-P.ss_c)*((1+cos(2*pi*((t+P.ss_t0)/365-P.ss_u2)))./2).^P.ss_k2);
    P.gM_fun = @(t) 0.5*P.NN*P.ss_S(t);
else
    P.gM_fun = @(t) 0.5*P.NN;
end
gH_fun = @(age) (2.*P.cc.*normpdf((age./365-P.zz)./P.ww).*normcdf(P.alpha.*(age./365-P.zz)./P.ww)./P.ww)./365/2;
gH =  gH_fun(a); % human fertility rate
muH =  P.b0 + P.b1*exp(-P.b2*a/365) + P.b3*exp(P.b4*a/365); % natural human mortality rate
muH_fun = @(age) P.b0 + P.b1*exp(-P.b2*age/365) + P.b3*exp(P.b4*age/365);
muH = muH/365;
muH_int_fun = @(age) (age./365).*P.b0 + (P.b1./P.b2).*(1-exp(-P.b2.*age./365)) + (P.b3./P.b4).*(-1+exp(P.b4.*age./365));
P.muH_int = muH_int_fun(a);

a74 = 74*365;
temp_muD = P.b0D + P.b1D*exp(-P.b2D*a74/365) + P.b3D*exp(P.b4D*a74/365);
muD =  (P.b0D + P.b1D*exp(-P.b2D*a/365) + P.b3D*exp(P.b4D*a/365)).*(a/365<=74)...
    + temp_muD.*(exp(-P.c0D*(a-74*365)./365)).*(a/365>74);
muD_fun = @(age) (P.b0D + P.b1D*exp(-P.b2D*age/365) + P.b3D*exp(P.b4D*age/365)).*(age/365<=74)...
    + temp_muD.*(exp(-P.c0D*(age/365-74))).*(age/365>74);
muD = muD/365;
P.muH = muH;
P.muH_fun = muH_fun;
P.muH_int_fun = muH_int_fun;

P.muD = muD;
P.muD_fun = muD_fun;
P.gH = gH;
P.gH_fun = gH_fun;
find_stable_age; % depends on gH and muH