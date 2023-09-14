close all
clear all
clc
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

% model parameters
Malaria_parameters_baseline;
Malaria_parameters_transform;
Malaria_parameters_transform_vac;

vac_param = 20;

t0_list= [0:0.5:9];

cases_total = NaN(size(t0_list));
vacc_total = NaN(size(t0_list));
cases_pp_py = NaN(size(t0_list));

%% run vaccination program seasonally
tic
for it = 1:length(t0_list)
    disp(['I am working on month ',num2str(t0_list(it))])
    [SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');
    %% time evolution - initial run (no vacc)
    tfinal = t0_list(it)*30;
    P.v0 = 0; Malaria_parameters_transform_vac;
    [t,SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(da,na,0,tfinal,...
        SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
    vacc_sterile = trapz(P.v*(1-P.z).*SH,1)*P.da;
    
    %% vac control - turn on at the prescribed month and continue for three months
    P.v0 = vac_param; Malaria_parameters_transform_vac;
    tfinal_conti = 3*30;
    SH0 = SH(:,end); EH0 = EH(:,end); DH0 = DH(:,end); AH0 = AH(:,end); VH0 = VH(:,end); UH0 = UH(:,end); MH0 = MH(:,end);
    SM0 = SM(end); EM0 = EM(end); IM0 = IM(end);
    Cac0 = Cac(:,end); Cm0 = Cm(:,end); Cv0 = Cv(:,end); Ctot0 = Ctot(:,end);
    [t2,SH2, EH2, DH2, AH2, VH2, UH2, SM2, EM2, IM2, Cm2, Cac2, Cv2, Ctot2, MH2] = age_structured_Malaria_vac(da,na,tfinal,tfinal+tfinal_conti,...
        SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
    vacc_sterile2 = trapz(P.v*(1-P.z).*SH2,1)*P.da; 

    t = [t;t2(2:end)];
    SH = [SH,SH2(:,2:end)]; EH = [EH,EH2(:,2:end)]; DH = [DH, DH2(:,2:end)]; AH = [AH, AH2(:,2:end)]; VH = [VH, VH2(:,2:end)]; UH = [UH, UH2(:,2:end)]; MH = [MH, MH2(:,2:end)];
    SM = [SM, SM2(2:end)]; EM = [EM, EM2(2:end)]; IM = [IM, IM2(2:end)]; 
    Cm = [Cm, Cm2(:,2:end)]; Cac = [Cac, Cac2(:,2:end)]; Cv = [Cv, Cv2(:,2:end)]; Ctot = [Ctot,Ctot2(:,2:end)];
    vacc_sterile = [vacc_sterile, vacc_sterile2(2:end)];
    %% continued run - vacc off and simulate the rest of the full year
    P.v0 = 0; Malaria_parameters_transform_vac;
    tfinal_conti2 = 365-t(end);
    if tfinal_conti2<0; keyboard;end
    SH0 = SH(:,end); EH0 = EH(:,end); DH0 = DH(:,end); AH0 = AH(:,end); VH0 = VH(:,end); UH0 = UH(:,end); MH0 = MH(:,end);
    SM0 = SM(end); EM0 = EM(end); IM0 = IM(end);
    Cac0 = Cac(:,end); Cm0 = Cm(:,end); Cv0 = Cv(:,end); Ctot0 = Ctot(:,end);
    [t2,SH2, EH2, DH2, AH2, VH2, UH2, SM2, EM2, IM2, Cm2, Cac2, Cv2, Ctot2, MH2] = age_structured_Malaria_vac(da,na,t(end),t(end)+tfinal_conti2,...
        SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
    vacc_sterile2 = trapz(P.v*(1-P.z).*SH2,1)*P.da; 
    t = [t;t2(2:end)];
    SH = [SH,SH2(:,2:end)]; EH = [EH,EH2(:,2:end)]; DH = [DH, DH2(:,2:end)]; AH = [AH, AH2(:,2:end)]; VH = [VH, VH2(:,2:end)]; UH = [UH, UH2(:,2:end)]; MH = [MH, MH2(:,2:end)];
    SM = [SM, SM2(2:end)]; EM = [EM, EM2(2:end)]; IM = [IM, IM2(2:end)]; 
    Cm = [Cm, Cm2(:,2:end)]; Cac = [Cac, Cac2(:,2:end)]; Cv = [Cv, Cv2(:,2:end)]; Ctot = [Ctot,Ctot2(:,2:end)];
    vacc_sterile = [vacc_sterile, vacc_sterile2(2:end)];

    % calculate QOIs
    PH = SH+EH+DH+AH+VH+UH;
    NM = SM+EM+IM;
    [bH,~] = biting_rate(PH,NM); % NB bH varies by age and time
    lamH = FOI_H(bH,IM,NM);
    rho = sigmoid_prob(Ctot./PH, 'rho'); % prob. of severely infected, EH -> DH
    psi = sigmoid_prob(Ctot./PH, 'psi'); % prob. AH -> DH
    temp1 = psi.*lamH.*AH; % AH -> DH
    temp2 = rho.*P.h.*EH;% EH -> DH, number of new cases - rate
    [~,ind1] = min(abs(P.a-5*30)); %[5,17] in White paper % [7 19] previous setting from RTS,S trial?
    [~,ind2] = min(abs(P.a-17*30));
    cases_rate1 = trapz(temp1(ind1:ind2,:),1)*P.da;
    cases_rate2 = trapz(temp2(ind1:ind2,:),1)*P.da;
    cases = cases_rate1+cases_rate2;
    cases_total(it) = trapz(cases)*P.dt; % total cases in cohort this year    
    vacc_total(it) = trapz(vacc_sterile)*P.dt; % total # vacc doeses this year
    pop = trapz(PH(ind1:ind2,:),1)*P.da;
    cases_pp_py(it) = cases_total(it)/mean(pop);
end

%% run baseline (no vaccination) for one year
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');
tfinal = 365;
P.v0 = 0; Malaria_parameters_transform_vac;
[t,SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(da,na,0,tfinal,...
    SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
PH = SH+EH+DH+AH+VH+UH;
NM = SM+IM+EM;
[bH,~] = biting_rate(PH,NM); % NB bH varies by age and time
lamH = FOI_H(bH,IM,NM);
rho = sigmoid_prob(Ctot./PH, 'rho'); % prob. of severely infected, EH -> DH
psi = sigmoid_prob(Ctot./PH, 'psi'); % prob. AH -> DH
temp1 = psi.*lamH.*AH; % AH -> DH
temp2 = rho.*P.h.*EH;% EH -> DH, number of new cases - rate
[~,ind1] = min(abs(P.a-5*30)); %[5,17] in White paper % [7 19] previous setting from RTS,S trial?
[~,ind2] = min(abs(P.a-17*30));
cases_rate1 = trapz(temp1(ind1:ind2,:),1)*P.da;
cases_rate2 = trapz(temp2(ind1:ind2,:),1)*P.da;
cases = cases_rate1+cases_rate2;
cases_total_baseline = trapz(cases)*P.dt; % total cases in cohort this year
vacc_total_baseline = trapz(trapz(P.v*(1-P.z).*SH,1)*P.da)*P.dt; % total # vacc doeses this year
pop = trapz(PH(ind1:ind2,:),1)*P.da;
cases_pp_py_baseline = cases_total_baseline/mean(pop); % should be around ~3 for Nanoro seasonlity study

%% run vaccination year long
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');
tfinal = 365;
P.v0 = vac_param; Malaria_parameters_transform_vac;
[t,SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(da,na,0,tfinal,...
    SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
PH = SH+EH+DH+AH+VH+UH;
NM = SM+IM+EM;
[bH,~] = biting_rate(PH,NM); % NB bH varies by age and time
lamH = FOI_H(bH,IM,NM);
rho = sigmoid_prob(Ctot./PH, 'rho'); % prob. of severely infected, EH -> DH
psi = sigmoid_prob(Ctot./PH, 'psi'); % prob. AH -> DH
temp1 = psi.*lamH.*AH; % AH -> DH
temp2 = rho.*P.h.*EH;% EH -> DH, number of new cases - rate
[~,ind1] = min(abs(P.a-5*30)); %[5,17] in White paper % [7 19] previous setting from RTS,S trial?
[~,ind2] = min(abs(P.a-17*30));
cases_rate1 = trapz(temp1(ind1:ind2,:),1)*P.da;
cases_rate2 = trapz(temp2(ind1:ind2,:),1)*P.da;
cases = cases_rate1+cases_rate2;
cases_total_yearlong = trapz(cases)*P.dt; % total cases in cohort this year
vacc_total_yearlong= trapz(trapz(P.v*(1-P.z).*SH,1)*P.da)*P.dt; % total # vacc doeses this year
pop = trapz(PH(ind1:ind2,:),1)*P.da;
cases_pp_py_yearlong = cases_total_yearlong/mean(pop); % should be around ~3 for Nanoro seasonlity study

%% compare vacc seaonal vs baseline
cases_per_vacc = (cases_total_baseline-cases_total)./vacc_total;
cases_per_vacc_baseline = zeros(size(t0_list)); 
cases_per_vacc_yearlong = (cases_total_baseline-cases_total_yearlong)/vacc_total_yearlong; 

toc
%% plotting
figure_setups; 
subplot(1,2,1)
plot(t0_list,cases_per_vacc_yearlong*ones(size(t0_list)),'-','DisplayName','year-long');
hold on
plot(t0_list,cases_per_vacc,'--','DisplayName',' seasonal');
plot(t0_list,cases_per_vacc_baseline,':','DisplayName','no vacc');
xlabel('vacc start month')
legend;
title('cases prevented per vacc')
subplot(1,2,2)
plot(t0_list,cases_pp_py_yearlong*ones(size(t0_list)),'-','DisplayName','year-long');
hold on
plot(t0_list,cases_pp_py,'--','DisplayName','seasonal');
plot(t0_list,cases_pp_py_baseline*ones(size(t0_list)),':','DisplayName',' no vacc');
legend;
xlabel('vacc start month')
title('cases per person per year')