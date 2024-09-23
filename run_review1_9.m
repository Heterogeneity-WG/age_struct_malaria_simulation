%% vac = constant, fixed number
% close all
clear all
% clc
format long
global P

tic

%% numerical config
age_max = 100*365; % max ages in days
P.age_max = age_max;
dt = 5; % time/age step size in days, default = 5;
da = dt;
a = (0:da:age_max)';
na = length(a);

P.dt = dt;
P.a = a;
P.na = na;
P.da = da;

Malaria_parameters_baseline;
Malaria_parameters_transform;
Malaria_parameters_transform_vac;

%% initial condition 'EE' - numerical EE
P.v0 = 0;
Malaria_parameters_transform_vac;
[SHEE, EHEE, DHEE, AHEE, VHEE, UHEE, SMEE, EMEE, IMEE, CmEE, CacEE, CvEE, CtotEE, ~] = age_structured_Malaria_IC_vac('EE_reset');
MHEE = zeros(size(SHEE));
%% seasonal vacc
t0 = 0; % starting month
vac_param_annual = 1.2*10^4;%6*10^4;
vac_period = 2; % seasonal vacc implementation period (months)
nyear  = 10; % implement vaccination strategy for nyear
%
tnow = 0;
SH0 = SHEE; EH0 = EHEE; DH0 = DHEE; AH0 = AHEE; VH0 = VHEE; UH0 = UHEE; SM0 = SMEE; EM0 = EMEE; IM0 = IMEE; Cm0 = CmEE; Cac0 = CacEE; Cv0 = CvEE; Ctot0 = CtotEE; MH0 = MHEE;
%% initial run (no vacc)
P.v0 = 0;
Malaria_parameters_transform_vac;
tconti = t0*30;
vac_season_time_evolution_init;
for iyear = 1:nyear
    % vac on - at the prescribed month and continue for three months
    vac_param = vac_param_annual/(vac_period*30);
    P.v0 = vac_param;
    Malaria_parameters_transform_vac;
    tconti = vac_period*30;
    vac_season_time_evolution_conti;
    % vac off - simulate the rest of the full year
    P.v0 = 0; Malaria_parameters_transform_vac;
    tconti = 365-vac_period*30;
    vac_season_time_evolution_conti;
end
%%
PH = SH+EH+DH+AH+VH+UH;
%% plotting
% plot_EIR_time;
% plot_human_prop_time; % human proportion in time (total and by age group)
% plot_human_popsize_time; % human pop size in time (total and by age group)
plot_human_pop_age_tfinal; % age distributions (in prop and size) at final time (and movie)
% plot_immunity;
% plot_DALY;
% plot_sigmoids;
% plot_FOI;
% plot_incidences;
% plot_seasonality;
% plot_mosquitoes;
% plot_vacc_counts_time;



