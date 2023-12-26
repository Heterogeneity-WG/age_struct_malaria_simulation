%% seasonality parameters  - White et al 2015 supp Table S6
% % for Nanoro (Burkina Faso), EIR peaks at August
P.ss_c = 0.02;
P.ss_v =  0.55; % v = 1 or 0 -> one peak; (0,1) -> two peaks
P.ss_k1 = 6.73;
P.ss_k2 = 1.68;
P.ss_u1 = 0.656;
P.ss_u2 = 0.841;
P.ss_S0 = 3.19; % magnitude of seasonlity profile, aim for 2.69 incidence rate
P.ss_t0 = 140; %  for correct peaking, time shift to incoporate delay