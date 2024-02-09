%% Code to generate Figure 2
% (summary of seasonal dynamics for low fixed, dynamic and
%   high fixed immunity scenarios)

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

for immunity_feedback = 1:3 % select the immune feedback scenario
    Malaria_parameters_baseline;
    Malaria_parameters_transform;
    Malaria_parameters_transform_vac;
    % sigmoid static values calculated with disease induced mortality, no
    % balancing and over a 5 year time horizon (see plot_sigmoids.m for averages)
    if immunity_feedback == 1
        % immunity values from betaM = 0.0196, fixed low immunity scenario
        P.phif0 = 0.270581008833871; % value at zero
        P.phif1 = 0.270581008833871; % value at L (function saturates to this value)
        P.rhof0 = 0.932188361414638; % value at zero
        P.rhof1 = 0.932188361414638; % value at L (function saturates to this value)
        P.psif0 = 0.749277341090716; % value at zero
        P.psif1 = 0.749277341090716; % value at L (function saturates to this value)
    elseif immunity_feedback == 3
        % immunity values from betaM = 0.75 (saturated), fixed high immunity scenario
        P.phif0 = 0.908066527085272; % value at zero
        P.phif1 = 0.908066527085272; % value at L (function saturates to this value)
        P.rhof0 = 0.090858132882051; % value at zero
        P.rhof1 = 0.090858132882051; % value at L (function saturates to this value
        P.psif0 = 0.082503847123926; % value at zero
        P.psif1 = 0.082503847123926; % value at L (function saturates to this value)
    elseif immunity_feedback == 2
        % don't change the sigmoids at all, dynamic immunity in force
    end

    %% initial condition 'EE' - numerical EE
    [SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');

    %% time evolution
    P.v0 = 0;
    Malaria_parameters_transform_vac; % resetting vaccination distribution
    tfinal = 365*5;
    [t,SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(da,na,0,tfinal,...
        SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
    PH = SH+EH+DH+AH+VH+UH;
    PH_final = PH(:,end);
    NH = trapz(PH,1)*da;
    NM = SM+IM+EM;
    %%
    [bH,~] = biting_rate(PH,NM); % bH(age, time) matrix
    EIR = bH.*IM./NM*365; % EIR(age, time) matrix
    EIR_tot = trapz(EIR.*PH,1)*P.da./NH;
    EIR_final = EIR_tot(end);

    %% plotting
    plot_incidences; 
    % plotting AH and DH over the target age ranges, plotting age cohort
    % dynamics

    %plot_EIR_time;
    % plot_human_prop_time; % human proportion in time (total and by age group)
    % plot_human_popsize_time; % human pop size in time (total and by age group)
    % plot_human_pop_age_tfinal; % age distributions (in prop and size) at final time (and movie)
    % plot_immunity;
    % plot_DALY;
    % plot_sigmoids;
    % plot_seasonality;
    % plot_mosquitoes;
    % plot_vacc_counts_time;
end



