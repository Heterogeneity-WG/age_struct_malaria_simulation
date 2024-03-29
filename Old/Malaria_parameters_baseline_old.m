global P

P.disease_mortality = 1; % turn on mu_D (disease induced mortality)
P.balance_fertility = 0; % balanced fertility or not
P.balance_mortality = 0; % balanced mortality or not
P.verbose = 1; % turn on the warning messages. Error messages from the check routines will display regardless

%% seasonality parameters  - White et al 2015 supp Table S6
% for Nanoro (Burkina Faso), EIR peaks at August
P.ss_c = 0.02;
P.ss_v =  0.55; % v = 1 or 0 -> one peak; (0,1) -> two peaks
P.ss_k1 = 6.73;
P.ss_k2 = 1.68;
P.ss_u1 = 0.656;
P.ss_u2 = 0.841;
P.ss_S0 = 0.68; % magnitude of seasonlity profile, aim for 2.69 incidence rate
P.ss_t0 = 100; % time shift to incoporate delay 

% for Siaya (Kenya), EIR has modes at May and November
% P.ss_c = 0.31;
% P.ss_v =  0.393; % v = 1 or 0 -> one peak; (0,1) -> two peaks
% P.ss_k1 = 4.08;
% P.ss_k2 = 3.66;
% P.ss_u1 = 0.003;
% P.ss_u2 = 0.456;
% P.ss_S0 = 0.32; % magnitude of seasonlity profile, aim for 3.15 incidence rate
% P.ss_t0 = 50; % time shift to incoporate delay 

% turn off seasonlity
P.ss_c = 1; P.ss_S0 = 1;

%% system configuration
P.lMsystem = 'ss'; % 'full' or 'ss'  full mosquito system or keep at quasi-SS
P.lMHfix = 'off'; % 'off' (default) or 'on' turn on/off the assumption on fixed mosquito-human ratio; off -> constant mosquito population; on -> exponentially grow with NH

%% dummy varaible for eFast SA
P.dummy = 1; P.dummy_lower = 0.65; P.dummy_upper = 2.1; % dummy parameter for global SA eFAST - values from original code 

%% vaccine related parameters
% RTS,S in Kenya --> Homa bay, Kisumu, Migori, Siaya, Busia, Bungoma, Vihiga, and Kakamega counties from wiki 2019 census
P.NN = 1131950+1155574+1116436+993183+893681+1670570+590013+1867579; % total Kenya population = 47,564,296; P.NN = 9,418,986; cohort 8500
P.v0 = 0; % default vaccination rate
P.v0s = P.v0; P.v0c = P.v0;
P.z = 0; P.z_lower = 0; P.z_upper = 1; % switch between sterilizing (1-z) and blood-stage (z)
P.dv = 5*365; % Half-life of vaccine-boosted immunity (Cv)
P.etas = 0.73; P.etas_lower = 0.4; P.etas_upper = 1; % Vaccine efficacy for the sterlizing immunity (VH) for children 
P.etab = 0.73; P.etab_lower = 0.4; P.etab_upper = 1; % Vaccine efficacy for the blood-stage immunity (Cv) for children 
P.w = 1/(0.66*365); P.w_lower = 1/(1.0*365); P.w_upper = 1/(0.31*365); % Waning rate for the sterlizing immunity (VH) for children 
%%
P.rD = 1/33.5; P.rD_lower = 1/51.5; P.rD_upper = 1/16; % recovery rate for DH (syptomatic, e.g. fever) P.rD_upper = 1/7;
P.rA = 1/85; P.rA_lower = 1/130; P.rA_upper = 1/40; % recovery rate for AH (clearance of parasite) 
P.h = 1/26; % incubation rate in human

%% immunity parameters/rates
P.dac = 5*365; P.dac_lower = 3.25*365; P.dac_upper = 10.5*365; % half life of acquired immunity   
P.dm = 0.25*365; % half life of maternal immunity
P.c1 = 1; % weight for acquired immunity
P.cX = 0.1; P.cX_lower = 0.065; P.cX_upper = 0.21; % free parameter - weight for boosting
P.m = 1; % fraction of new-born immunity relative to mother’s
P.uc = 10; % Duration in which immunity is not boosted
%% progression probabilities parameters, sigmoid parameters
% fitted values using Tfinal = 10 years
x = [2.567957971786876   2.487540758554113   3.649596968324358   1.395449806257184   2.332526365071812   2.150211932758257];
P.phis2 = x(1);
P.phir2 = x(2); 
P.rhos2 = x(3);
P.rhor2 = x(4); 
P.psis2 = x(5);
P.psir2 = x(6);

P.phif0 = 0.01; 
P.phif1 = 1;
P.rhof0 = 0.01; 
P.rhof1 = 1; 
P.psif0 = 0.01; 
P.psif1 = 1; 

P.phis2_lower = P.phis2*0.65; P.phis2_upper = P.phis2*2.1;
P.phir2_lower = P.phir2*0.65; P.phir2_upper = P.phir2*2.1;
P.rhos2_lower = P.rhos2*0.65; P.rhos2_upper = P.rhos2*2.1;
P.rhor2_lower = P.rhor2*0.65; P.rhor2_upper = P.rhor2*2.1;
P.psis2_lower = P.psis2*0.65; P.psis2_upper = P.psis2*2.1;
P.psir2_lower = P.psir2*0.65; P.psir2_upper = P.psir2*2.1;

%% mosquito related parameters/rates
P.bh = 5; % tolerated biting rate per human
P.bm = 0.6; % desired biting rate per mosquito
P.betaM = 0.25; P.betaM_lower = 0.16; P.betaM_upper = 0.53; % infectivity of mosquitoes P.betaM = 0.25;
P.betaD = 0.35; P.betaD_lower = 0.23; P.betaD_upper = 0.74; % infectivity of DH   P.betaD = 0.35;
P.betaA = 0.03; P.betaA_lower = 0.02; P.betaA_upper = 0.06; % infectivity of AH    P.betaA = 0.03;

P.muM = 1/14; P.muM_lower = 1/21.5; P.muM_upper = 1/6.7; % natural mortality rate of mosquitoes
P.gM = 0.5*P.NN; % recruitment rate of mosquitoes;
P.MHm = P.gM/P.muM/P.NN; % assume NH = 1; % mosquito/human ratio
P.sigma = 1/10; P.sigma_lower = 1/15.4; P.sigma_upper = 1/4.8; % incubation rate for mosquitoes
%% muH: non-malaria related mortality rate parameters
% use GHO life tables, nMx data
% P.b0 = 0.0024214446844162;
% P.b1 = 0.0887924178445357;
% P.b2 = 2.09862983723212;
% P.b3 = 6.87709371762464e-05;
% P.b4 = 0.0901695513967616;
% use IHME data
P.b0 = -0.00060391;
P.b1 = 0.056826;
P.b2 = 1.2853;
P.b3 = 0.00028799;
P.b4 = 0.072469;

%% mu_D: malaria related mortality rate parameters
% use IHME data
P.b0D = 1.1812*10^-5;
P.b1D = 0.055085;
P.b2D = 1.0282;
P.b3D = 0.00012962;
P.b4D = 0.064575;
P.c0D = 0.181055030942671;

%% fertility rate parameters
% use IHME data fit
P.cc = 3.2018;
P.zz = 18.213;
P.alpha = 5.7911;
P.ww = 12.603;
% use rawdata - Kenya
% P.cc = 4.024086261410830;
% P.zz = 17.963601264000353;
% P.alpha = 4.083610527018673;
% P.ww = 13.196127635937707;

%%
% Malaria_parameters_transform; % commented out for running SA

