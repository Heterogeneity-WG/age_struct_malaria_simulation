%% Compare different/no seasonality profiles (gM, EIR, immunity)
clear all
close all
% clc

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

figure_setups;
figure_setups;
figure_setups;

for seasonality = 0:1 % select the seasonality scenario
    Malaria_parameters_baseline;
    Malaria_parameters_transform;
    Malaria_parameters_transform_vac;
    if seasonality == 1 % turn on the seasonality
        % EIR average ~ 65
        % Malaria_parameters_baseline_Nanoro;
        Malaria_parameters_baseline_Siaya;
        Malaria_parameters_transform;
        Malaria_parameters_transform_vac;
        label_legend = 'Mild seasonality (Siaya)';
        % label_legend = 'Strong seasonality (Nanoro)';
        label_line = '-';
    elseif  seasonality == 0 % turn off the seasonality
        Malaria_parameters_baseline_Siaya;
        P.ss_c = 1; P.ss_S0 = 1;
        % P.betaM = 0.1; % adjust EIR to match ~ 65
        Malaria_parameters_transform;
        Malaria_parameters_transform_vac;
        label_legend = 'No seasonality';
        % label_legend = 'Strong seasonality (Nanoro)';
        label_line = '--';
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
    % plot_EIR_time;
    % plot_incidences;
    % plot_human_prop_time; % human proportion in time (total and by age group)
    % plot_human_popsize_time; % human pop size in time (total and by age group)
    % plot_human_pop_age_tfinal; % age distributions (in prop and size) at final time (and movie)

    % plot_immunity;
    [~,ind_EIR] = max(EIR_tot);
    ['EIR=', num2str(mean(EIR_tot))]
    ['gM=', num2str(mean(P.gM_fun(t)),'%e')]

    figure(1); hold on
    plot(a/365,Ctot(:,ind_EIR)./PH_final,label_line,'DisplayName',label_legend);
    xlabel('Age (years)')
    ylabel('Immunity level')
    title('Per-person immunity distribution');
    axis([0 age_max/365 0 max(Ctot(:,ind_EIR)./PH_final)*1.1]);
    grid on

    figure(2); hold on
    plot(a/365,DH(:,ind_EIR)./(AH(:,ind_EIR)+DH(:,ind_EIR)),label_line,'DisplayName',label_legend);
    xlabel('Age (years)')
    ylabel('Proportion')
    title('Fraction of symptomatic infection');
    axis([0 age_max/365 0 1]);
    grid on

    figure(3); hold on
    yyaxis left
    hold on
    plot(t/365,EIR_tot,label_line,'DisplayName',['EIR(t)',label_legend]);
    ylabel('EIR')
    % ylim([0 120])
    yyaxis right
    hold on
    plot(t/365, P.gM_fun(t),label_line,'DisplayName',['$g_M(t)$',label_legend])
    ylabel('Mosquito recruitment rate')

    % plot_DALY;
    % plot_sigmoids;
    % plot_seasonality;
    % plot_mosquitoes;
    % plot_vacc_counts_time;
end
figure(1)
legend('Location','North');
xlim([0 10])
ylim([0 11])

figure(2)
legend('Location','North');
xlim([0 10])

figure(3)
legend('Location','East');
xlim([0 1])