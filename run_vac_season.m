% close all
clear all
clc
format long
global P

tic

%% numerical config
age_max = 100*365; % max ages in days
P.age_max = age_max;
dt = 1; % time/age step size in days, default = 5;
da = dt;
a = (0:da:age_max)';
na = length(a);

P.dt = dt;
P.a = a;
P.na = na;
P.da = da;

%% seasonal vac with constant vac rates
% model parameters
Malaria_parameters_baseline;
Malaria_parameters_baseline_Siaya; % choose seasonality profile here
Malaria_parameters_transform;
Malaria_parameters_transform_vac;

t0_list= (0:0.5:12)';
nyear  = 10; % implement vaccination strategy for nyear
vac_period = 3; % seasonal vacc implementation period (months)
vac_param_annual = (6*10^4); % annual vacc number % target population total = 2.56*10^5;
% pick 1.2*10^4 baseline acc count to avoid negative SH

% NB "low vacc" = 1.2*10^4, "high vacc = 6*10^4"

% allocation
cases_py_target = NaN(length(t0_list),nyear);
cases_pp_py_target = NaN(length(t0_list),nyear);
vac_py_target = NaN(length(t0_list),nyear);
death_py_target = NaN(length(t0_list),nyear);
cases_py_full = NaN(length(t0_list),nyear);
cases_pp_py_full = NaN(length(t0_list),nyear);
vac_py_full = NaN(length(t0_list),nyear);
death_py_full = NaN(length(t0_list),nyear);

cases_py_target_baseline = NaN(length(t0_list),nyear);
cases_pp_py_target_baseline = NaN(length(t0_list),nyear);
vac_py_target_baseline = NaN(length(t0_list),nyear);
death_py_target_baseline = NaN(length(t0_list),nyear);
cases_py_full_baseline = NaN(length(t0_list),nyear);
cases_pp_py_full_baseline = NaN(length(t0_list),nyear);
vac_py_full_baseline = NaN(length(t0_list),nyear);
death_py_full_baseline = NaN(length(t0_list),nyear);

cases_py_target_constant = NaN(length(t0_list),nyear);
cases_pp_py_target_constant = NaN(length(t0_list),nyear);
vac_py_target_constant = NaN(length(t0_list),nyear);
death_py_target_constant = NaN(length(t0_list),nyear);
cases_py_full_constant = NaN(length(t0_list),nyear);
cases_pp_py_full_constant = NaN(length(t0_list),nyear);
vac_py_full_constant = NaN(length(t0_list),nyear);
death_py_full_constant = NaN(length(t0_list),nyear);

P.v0 = 0; Malaria_parameters_transform_vac;
[SHEE, EHEE, DHEE, AHEE, VHEE, UHEE, SMEE, EMEE, IMEE, CmEE, CacEE, CvEE, CtotEE, ~] = age_structured_Malaria_IC_vac('EE_reset');
MHEE = zeros(size(SHEE));
tic

%% run vaccination program seasonally
for it = 1:length(t0_list)
    disp(['I am working on month ',num2str(t0_list(it))])
    tnow = 0;
    SH0 = SHEE; EH0 = EHEE; DH0 = DHEE; AH0 = AHEE; VH0 = VHEE; UH0 = UHEE; SM0 = SMEE; EM0 = EMEE; IM0 = IMEE; Cm0 = CmEE; Cac0 = CacEE; Cv0 = CvEE; Ctot0 = CtotEE; MH0 = MHEE;
    %% initial run (no vacc)
    P.v0 = 0;
    Malaria_parameters_transform_vac;
    tconti = t0_list(it)*30;
    vac_season_time_evolution_init;
    for iyear = 1:nyear
        % vac on - at the prescribed month and continue for three months
        vac_param = vac_param_annual/(vac_period*30);
        P.v0 = vac_param; 
        Malaria_parameters_transform_vac;
        tconti = vac_period*30;
        vac_season_time_evolution_init;
        % vac off - simulate the rest of the full year
        P.v0 = 0; Malaria_parameters_transform_vac;
        tconti = 365-vac_period*30;
        vac_season_time_evolution_conti;

        [cases_py_target(it,iyear), cases_pp_py_target(it,iyear), vac_py_target(it,iyear), death_py_target(it, iyear)] = ...
            vac_season_time_stats(t,'target', SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH, vacc_sterile);
        [cases_py_full(it,iyear), cases_pp_py_full(it,iyear), vac_py_full(it,iyear), death_py_full(it, iyear)] = ...
            vac_season_time_stats(t,'full', SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH, vacc_sterile);
    end
end

%% run baseline (no vacc)
for it = 1:length(t0_list)
    disp(['I am working on baseline - month ',num2str(t0_list(it))])
    tnow = 0;
    SH0 = SHEE; EH0 = EHEE; DH0 = DHEE; AH0 = AHEE; VH0 = VHEE; UH0 = UHEE; SM0 = SMEE; EM0 = EMEE; IM0 = IMEE; Cm0 = CmEE; Cac0 = CacEE; Cv0 = CvEE; Ctot0 = CtotEE; MH0 =  MHEE;
    % initial run (no vacc)
    P.v0 = 0; Malaria_parameters_transform_vac;
    tconti = t0_list(it)*30;
    vac_season_time_evolution_init;
    for iyear = 1:nyear
        P.v0 = 0; 
        Malaria_parameters_transform_vac;
        tconti = 365;
        vac_season_time_evolution_init;
        [cases_py_target_baseline(it,iyear), cases_pp_py_target_baseline(it,iyear), vac_py_target_baseline(it,iyear),...
            death_py_target_baseline(it,iyear)] = vac_season_time_stats(t,'target', SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH, vacc_sterile);
        [cases_py_full_baseline(it,iyear), cases_pp_py_full_baseline(it,iyear), vac_py_full_baseline(it,iyear),...
            death_py_full_baseline(it,iyear)] = vac_season_time_stats(t,'full', SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH, vacc_sterile);
    end
end

%% run constant vacc
for it = 1:length(t0_list)
    disp(['I am working on constant - month ',num2str(t0_list(it))])
    tnow = 0;
    SH0 = SHEE; EH0 = EHEE; DH0 = DHEE; AH0 = AHEE; VH0 = VHEE; UH0 = UHEE; SM0 = SMEE; EM0 = EMEE; IM0 = IMEE; Cm0 = CmEE; Cac0 = CacEE; Cv0 = CvEE; Ctot0 = CtotEE; MH0 =  MHEE;
    % initial run (no vacc)
    P.v0 = 0; Malaria_parameters_transform_vac;
    tconti = t0_list(it)*30;
    vac_season_time_evolution_init;
    for iyear = 1:nyear
        vac_param = vac_param_annual/365;
        P.v0 = vac_param; Malaria_parameters_transform_vac;
        tconti = 365;
        vac_season_time_evolution_init;
        [cases_py_target_constant(it,iyear), cases_pp_py_target_constant(it,iyear), vac_py_target_constant(it,iyear),...
            death_py_target_constant(it,iyear)] = vac_season_time_stats(t,'target', SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH, vacc_sterile);
        [cases_py_full_constant(it,iyear), cases_pp_py_full_constant(it,iyear), vac_py_full_constant(it,iyear),...
            death_py_full_constant(it,iyear)] = vac_season_time_stats(t,'full', SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH, vacc_sterile);
    end
end
toc
%% compare vacc seaonal vs baseline
cases_per_vacc_target = (cases_py_target_baseline-cases_py_target)./vac_py_target;
death_per_vacc_target = (death_py_target_baseline-death_py_target)./vac_py_target;
cases_per_vacc_target_baseline = zeros(length(t0_list),nyear);
death_per_vacc_target_baseline = zeros(length(t0_list),nyear);
cases_per_vacc_target_constant = (cases_py_target_baseline-cases_py_target_constant)./vac_py_target_constant;
death_per_vacc_target_constant = (death_py_target_baseline-death_py_target_constant)./vac_py_target_constant;

cases_per_vacc_full = (cases_py_full_baseline-cases_py_full)./vac_py_full;
death_per_vacc_full = (death_py_full_baseline-death_py_full)./vac_py_full;
cases_per_vacc_full_baseline = zeros(length(t0_list),nyear);
death_per_vacc_full_baseline = zeros(length(t0_list),nyear);
cases_per_vacc_full_constant = (cases_py_full_baseline-cases_py_full_constant)./vac_py_full_constant;
death_per_vacc_full_constant = (death_py_full_baseline-death_py_full_constant)./vac_py_full_constant;

%% saving data
direc = 'Results/season_vacc/';
if vac_param_annual == 6*10^4 % high vacc count
    filename = [direc,'season_vacc_',num2str(vac_period),'_high.mat'];
elseif vac_param_annual == 1.2*10^4 % low vacc count
    filename = [direc,'season_vacc_',num2str(vac_period),'_low.mat'];
elseif vac_param_annual == (1.2*10^4)/100 % very low vacc count
    filename = [direc,'season_vacc_',num2str(vac_period),'_verylow.mat'];
else
    keyboard;
end
save(filename,'cases_per_vacc_target','death_per_vacc_target','cases_per_vacc_target_constant','death_per_vacc_target_constant',...
    'cases_per_vacc_target_baseline','death_per_vacc_target_baseline',...
    'cases_per_vacc_full', 'death_per_vacc_full','cases_per_vacc_full_constant','death_per_vacc_full_constant',...
    'cases_per_vacc_full_baseline','death_per_vacc_full_baseline',...
    'cases_pp_py_target','cases_pp_py_full','cases_pp_py_target_baseline','cases_pp_py_full_baseline',...
    'cases_pp_py_target_constant','cases_pp_py_full_constant',...
    'vac_py_target','vac_py_full','vac_py_target_baseline','vac_py_full_baseline','vac_py_target_constant','vac_py_full_constant',...
    't0_list','nyear','vac_period','vac_param_annual','P');

%% plotting (entire pop)
% figure_setups;  % plot cases prevented per vacc
% h = subplot(1,3,1);
% f1 = plot(t0_list,cases_per_vacc_full_constant,'-','DisplayName','year-long (full pop)');
% hold on; grid on
% f2 = plot(t0_list,cases_per_vacc_full,'--','DisplayName',' seasonal (full pop)');
% f3 = plot(t0_list,cases_per_vacc_full_baseline,':','DisplayName','no vacc (full pop)');
% xlabel('vacc start month')
% ll = legendUnq(h);
% legend(ll,'Location','best')
% title('cases prev. per vacc')
% xlim([0 12])
% set(f1, {'color'}, num2cell(winter(nyear),2));
% set(f2, {'color'}, num2cell(winter(nyear),2));
% set(f3, {'color'}, num2cell(winter(nyear),2));
%
% h = subplot(1,3,2); % plot death prevented per vacc
% f1 = plot(t0_list,death_per_vacc_full_constant,'-','DisplayName','year-long (full pop)');
% hold on; grid on
% f2 = plot(t0_list,death_per_vacc_full,'--','DisplayName','seasonal (full pop)');
% f3 = plot(t0_list,death_per_vacc_full_baseline,':','DisplayName',' no vacc (full pop)');
% xlabel('vacc start month')
% ll = legendUnq(h);
% legend(ll,'Location','best')
% title('death prev. per vacc')
% xlim([0 12])
% set(f1, {'color'}, num2cell(winter(nyear),2));
% set(f2, {'color'}, num2cell(winter(nyear),2));
% set(f3, {'color'}, num2cell(winter(nyear),2));
%
% h = subplot(1,3,3);
% f1 = plot(t0_list,cases_pp_py_full_constant,'-','DisplayName','year-long (full pop)');
% hold on; grid on
% f2 = plot(t0_list,cases_pp_py_full,'--','DisplayName','seasonal (full pop)');
% f3 = plot(t0_list,cases_pp_py_full_baseline,':','DisplayName',' no vacc (full pop)');
% xlabel('vacc start month')
% ll = legendUnq(h);
% legend(ll,'Location','best')
% title('cases per person per year')
% xlim([0 12])
% set(f1, {'color'}, num2cell(winter(nyear),2));
% set(f2, {'color'}, num2cell(winter(nyear),2));
% set(f3, {'color'}, num2cell(winter(nyear),2));
%% raw stats
% figure_setups;
% subplot(1,3,1)
% h1 = plot(t0_list,cases_py_target,'--');
% hold on
% h2 = plot(t0_list,cases_py_target_constant,'-');
% h3 = plot(t0_list,cases_py_target_baseline,':');
% set(h1, {'color'}, num2cell(winter(nyear),2));
% set(h2, {'color'}, num2cell(winter(nyear),2));
% set(h3, {'color'}, num2cell(winter(nyear),2));
% title('cases per year')
%
% subplot(1,3,2)
% h1 = plot(t0_list,death_py_target,'--');
% hold on
% h2 = plot(t0_list,death_py_target_constant,'-');
% h3 = plot(t0_list,death_py_target_baseline,':');
% set(h1, {'color'}, num2cell(winter(nyear),2));
% set(h2, {'color'}, num2cell(winter(nyear),2));
% set(h3, {'color'}, num2cell(winter(nyear),2));
% title('death per year')
%
% subplot(1,3,3)
% h1 = plot(t0_list,vac_py_target,'--');
% hold on
% h2 = plot(t0_list,vac_py_target_constant,'-');
% h3 = plot(t0_list,vac_py_target_baseline,':');
% set(h1, {'color'}, num2cell(winter(nyear),2));
% set(h2, {'color'}, num2cell(winter(nyear),2));
% set(h3, {'color'}, num2cell(winter(nyear),2));
% title('vac per year')

function [cases_py, cases_pp_py, vac_py, death_py] = vac_season_time_stats(t,lage, SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH, vacc_sterile)
global P
% calculate QOIs
PH = SH+EH+DH+AH+VH+UH;
NM = SM+EM+IM;
[bH,~] = biting_rate(PH,NM); % NB bH varies by age and time
lamH = FOI_H(bH,IM,NM);
rho = sigmoid_prob(Ctot./PH, 'rho'); % prob. of severely infected, EH -> DH
psi = sigmoid_prob(Ctot./PH, 'psi'); % prob. AH -> DH
temp1 = psi.*lamH.*AH; % AH -> DH
temp2 = rho.*P.h.*EH;% EH -> DH, number of new cases - rate

if strcmp(lage,'target')
    [~,ind1] = min(abs(P.a-5*30)); %[5,17] in White paper % [7 19] previous setting from RTS,S trial?
    [~,ind2] = min(abs(P.a-17*30));
elseif strcmp(lage, 'full')
    ind1 = 1;
    ind2 = length(P.a);
end

cases_rate1 = trapz(temp1(ind1:ind2,:),1)*P.da;
cases_rate2 = trapz(temp2(ind1:ind2,:),1)*P.da;
cases = cases_rate1+cases_rate2;
cases_py = trapz(cases)*P.dt/(t(end)-t(1))*365; % total cases in cohort this year
vac_py = trapz(trapz(vacc_sterile,1)*P.da)*P.dt/(t(end)-t(1))*365; % total # vacc doeses this year
pop = trapz(PH(ind1:ind2,:),1)*P.da;
cases_pp_py = cases_py/mean(pop);
death_py = trapz(trapz(P.muD.*DH,1))*P.dt*P.da; % trapz(MH(:,end),1)*P.da;
end
