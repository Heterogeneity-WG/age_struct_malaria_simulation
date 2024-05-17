%% Code to generate Figure 2
% (summary of seasonal dynamics for low fixed, dynamic and
%   high fixed immunity scenarios)

global P
tic
%% numerical config
age_max = 100*365; % max ages in days
P.age_max = age_max;
dt = 0.5; % time/age step size in days, default = 1;
da = dt;
a = (0:da:age_max)';
na = length(a);

P.dt = dt;
P.a = a;
P.na = na;
P.da = da;

immunity_feedback = 2;
for count = 1:2 % select the immune feedback scenario
    % 1 = fixed low, 2 = dynamic (standard), 3 = fixed high
    Malaria_parameters_baseline;
    Malaria_parameters_baseline_Nanoro;
    Malaria_parameters_transform;
    Malaria_parameters_transform_vac;

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

    plot_EIR_time;
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



