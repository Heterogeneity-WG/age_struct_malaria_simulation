%% seasonality parameters  - White et al 2015 supp Table S6
% for Siaya (Kenya), EIR has modes at May and November
P.ss_c = 0.31;
P.ss_v =  0.393; % v = 1 or 0 -> one peak; (0,1) -> two peaks
P.ss_k1 = 4.08;
P.ss_k2 = 3.66;
P.ss_u1 = 0.003;
P.ss_u2 = 0.456;
P.ss_S0 = 2.6; % magnitude of seasonlity profile, aim for 3.15 (3.09 calibrated) incidence rate
P.ss_t0 = 100; % time shift to incoporate delay 