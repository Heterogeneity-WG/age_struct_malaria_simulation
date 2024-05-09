global P

P.disease_mortality = 1; % turn on mu_D (disease induced mortality)
P.balance_mortality = 0; % balanced mortality or not
P.balance_fertility = 0; % balanced fertility or not
P.verbose = 1; % turn on the warning messages. Error messages from the check routines will display regardless

%% system configuration
P.lMsystem = 'full'; % 'full' or 'ss'  full mosquito system or keep at quasi-SS
P.lMHfix = 'off'; % 'off' (default) or 'on' turn on/off the assumption on fixed mosquito-human ratio; off -> constant mosquito population; on -> exponentially grow with NH

%% dummy varaible for eFast SA
P.dummy = 1; P.dummy_lower = 0.65; P.dummy_upper = 2.1; % dummy parameter for global SA eFAST - values from original code 

%% vaccine related parameters
% RTS,S in Kenya --> Homa bay, Kisumu, Migori, Siaya, Busia, Bungoma, Vihiga, and Kakamega counties from wiki 2019 census
P.NN = 1131950+1155574+1116436+993183+893681+1670570+590013+1867579; % total Kenya population = 47,564,296; P.NN = 9,418,986; 
% (5, 17) =216,893;  cohort = 8500
P.v0 = 0; 
% P.v0 = 1.2*10^5/365; 
P.v0_lower = 0; P.v0_upper = 1.2*10^5*2.1/365; % vaccination rate 
P.v0s = P.v0; P.v0c = P.v0;
P.z = 0; P.z_lower = 0; P.z_upper = 1; % switch between sterilizing (P.z = 0) and blood-stage (P.z = 1)
P.dv = 5*365; % averaged period of vaccine-boosted immunity (Cv)
P.etas = 0.82; P.etas_lower = 0.53; P.etas_upper = 1; % Vaccine efficacy for the sterlizing immunity (VH) for children P.etas = 0.73
P.etab = 0.82; P.etab_lower = 0.53; P.etab_upper = 1; % Vaccine efficacy for the blood-stage immunity (Cv) for children 
P.w = 1/(0.41*365); P.w_lower = 1/(0.63*365); P.w_upper = 1/(0.2*365); % Waning rate for the sterlizing immunity (VH) for children  
%%
P.rD = 1/31; P.rD_lower = 1/48; P.rD_upper = 1/15; % recovery rate for DH (syptomatic, e.g. fever) P.rD_upper = 1/7;
P.rA = 1/85; P.rA_lower = 1/130; P.rA_upper = 1/40; % recovery rate for AH (clearance of parasite) 
P.h = 1/26; % incubation rate in human

%% immunity parameters/rates
P.dac = 5*365; P.dac_lower = 3.25*365; P.dac_upper = 10.5*365; % averaged waning period of acquired immunity   
P.dm = 0.25*365; % half life of maternal immunity
P.c1 = 1; % weight for acquired immunity
%%
cX = 0.1; P.cX = cX;
P.cS = (1-2.5*cX)/2; P.cS_lower = P.cS*0.65; P.cS_upper = P.cS*2.1;
P.cE = cX; P.cE_lower = P.cE*0.65; P.cE_upper = P.cE*2.1;
P.cA = cX; P.cA_lower = P.cA*0.65; P.cA_upper = P.cA*2.1;
P.cD = 0.5*cX; P.cD_lower = P.cD*0.65; P.cD_upper = P.cD*2.1;
P.cU = P.cS; P.cU_lower = P.cU*0.65; P.cU_upper = P.cU*2.1;
P.cV = P.cS; % weight for vaccination ~~ SH
P.m = 1; P.m_lower = P.m*0.65; P.m_upper = P.m*2.1; % fraction of new-born immunity relative to mother's
P.uc = 10; P.uc_lower = P.uc*0.65; P.uc_upper = P.uc*2.1; % Duration in which immunity is not boosted
%% progression probabilities parameters, sigmoid parameters
% fitted values using Tfinal = 10 years
x_fit = [2.541054908661499   2.521263428837834   3.160535717325357...
    1.870698094604934   2.390238620104622   2.282195667153324];
P.phis2 = x_fit(1);
P.phir2 = x_fit(2); 
P.rhos2 = x_fit(3);
P.rhor2 = x_fit(4); 
P.psis2 = x_fit(5);
P.psir2 = x_fit(6);

% dynamic
P.phif0 = 0.01; 
P.phif1 = 1;
P.rhof0 = 0.01; 
P.rhof1 = 1; 
P.psif0 = 0.01; 
P.psif1 = 1; 

% fixed high immunity scenario
% P.phif0 = 0.908066527085272; % value at zero
% P.phif1 = 0.908066527085272; % value at L (function saturates to this value)
% P.rhof0 = 0.090858132882051; % value at zero
% P.rhof1 = 0.090858132882051; % value at L (function saturates to this value
% P.psif0 = 0.082503847123926; % value at zero
% P.psif1 = 0.082503847123926; % value at L (function saturates to this value)

% fixed low immunity scenario
% P.phif0 = 0.270581008833871; % value at zero
% P.phif1 = 0.270581008833871; % value at L (function saturates to this value)
% P.rhof0 = 0.932188361414638; % value at zero
% P.rhof1 = 0.932188361414638; % value at L (function saturates to this value)
% P.psif0 = 0.749277341090716; % value at zero
% P.psif1 = 0.749277341090716; % value at L (function saturates to this value)

P.phis2_lower = P.phis2*0.65; P.phis2_upper = P.phis2*2.1;
P.phir2_lower = P.phir2*0.65; P.phir2_upper = P.phir2*2.1;
P.rhos2_lower = P.rhos2*0.65; P.rhos2_upper = P.rhos2*2.1;
P.rhor2_lower = P.rhor2*0.65; P.rhor2_upper = P.rhor2*2.1;
P.psis2_lower = P.psis2*0.65; P.psis2_upper = P.psis2*2.1;
P.psir2_lower = P.psir2*0.65; P.psir2_upper = P.psir2*2.1;
%% mosquito related parameters/rates
P.bh = 5; % tolerated biting rate per human
P.bm = 0.6; % desired biting rate per mosquito
P.betaM = 0.35; P.betaM_lower = 0.23; P.betaM_upper = 0.74; % infectivity of mosquitoes P.betaM = 0.25;
P.betaD = 0.35; P.betaD_lower = 0.23; P.betaD_upper = 0.74; % infectivity of DH   P.betaD = 0.35;
P.betaA = 0.03; P.betaA_lower = 0.02; P.betaA_upper = 0.06; % infectivity of AH    P.betaA = 0.03;

P.muM = 1/14; P.muM_lower = 1/21.5; P.muM_upper = 1/6.7; % natural mortality rate of mosquitoes
% P.MHm = P.gM/P.muM/P.NN; % assume NH = 1; % mosquito/human ratio
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
% Malaria_parameters_transform_SA; % commented out for running SA

