% close all
% clc
format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

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

t0_list= [1:0.5:9];
for it = 1:length(t0_list)
    disp(['I am working on month ',num2str(t0_list(it))])
    %% initial condition 'EE' - numerical EE
    [SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');
    NN_S = trapz(SH0)*P.da;
    %% time evolution - initial run (no vacc)
    tfinal = t0_list(it)*30;
    P.v0 = 0; P.z = 0; Malaria_parameters_transform_vac;
    [t,SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(da,na,0,tfinal,...
        SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
    PH = SH+EH+DH+AH+VH+UH;
    PH_final = PH(:,end); % total human at age a, t = n
    NH = trapz(PH,1)*da;
    NM = SM+IM+EM;
    vacc_sterile = trapz(P.v*(1-P.z).*SH,1)*P.da*365*P.NN/1000;
    vacc_blood = trapz(P.v*P.z.*SH,1)*P.da*365*P.NN/1000;
    
    %% vac control
    P.v0 = 20; 
    P.z = 0; % z=0 sterile, z=1 blood-stage
    Malaria_parameters_transform_vac;
    vacc_fun = P.v;
    tfinal_conti = 365;%3*30;
    SH0 = SH(:,end); EH0 = EH(:,end); DH0 = DH(:,end); AH0 = AH(:,end); VH0 = VH(:,end); UH0 = UH(:,end); MH0 = MH(:,end);
    SM0 = SM(end); EM0 = EM(end); IM0 = IM(end);
    Cac0 = Cac(:,end); Cm0 = Cm(:,end); Cv0 = Cv(:,end); Ctot0 = Ctot(:,end);
    [t2,SH2, EH2, DH2, AH2, VH2, UH2, SM2, EM2, IM2, Cm2, Cac2, Cv2, Ctot2, MH2] = age_structured_Malaria_vac(da,na,tfinal,tfinal+tfinal_conti,...
        SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
    PH2 = SH2+EH2+DH2+AH2+VH2+UH2;
    NH2 = trapz(PH2,1)*da;
    NM2 = SM2+EM2+IM2;
    vacc_sterile2 = trapz(P.v*(1-P.z).*SH2,1)*P.da;
    vacc_blood2 = trapz(P.v*P.z.*SH2,1)*P.da;
    
    t = [t;t2(2:end)];
    SH = [SH,SH2(:,2:end)]; EH = [EH,EH2(:,2:end)]; DH = [DH, DH2(:,2:end)]; AH = [AH, AH2(:,2:end)]; VH = [VH, VH2(:,2:end)]; UH = [UH, UH2(:,2:end)]; MH = [MH, MH2(:,2:end)];
    SM = [SM, SM2(2:end)]; EM = [EM, EM2(2:end)]; IM = [IM, IM2(2:end)]; NM = [NM, NM2(2:end)];
    Cm = [Cm, Cm2(:,2:end)]; Cac = [Cac, Cac2(:,2:end)]; Cv = [Cv, Cv2(:,2:end)]; Ctot = [Ctot,Ctot2(:,2:end)];
    PH = SH+EH+DH+AH+VH+UH;
    PH_final = PH(:,end); % total human at age a, t = n
    NH = [NH, NH2(2:end)];
    vacc_sterile = [vacc_sterile, vacc_sterile2(2:end)];
    vacc_blood = [vacc_blood, vacc_blood2(2:end)];
    %% continued run - vacc off
    P.v0 = 0;
    P.z = 0; % z=0 sterile, z=1 blood-stage
    Malaria_parameters_transform_vac;
    tfinal_conti2 = 365-t(end);
    SH0 = SH(:,end); EH0 = EH(:,end); DH0 = DH(:,end); AH0 = AH(:,end); VH0 = VH(:,end); UH0 = UH(:,end); MH0 = MH(:,end);
    SM0 = SM(end); EM0 = EM(end); IM0 = IM(end);
    Cac0 = Cac(:,end); Cm0 = Cm(:,end); Cv0 = Cv(:,end); Ctot0 = Ctot(:,end);
    [t2,SH2, EH2, DH2, AH2, VH2, UH2, SM2, EM2, IM2, Cm2, Cac2, Cv2, Ctot2, MH2] = age_structured_Malaria_vac(da,na,t(end),t(end)+tfinal_conti2,...
        SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
    PH2 = SH2+EH2+DH2+AH2+VH2+UH2;
    NH2 = trapz(PH2,1)*da;
    NM2 = SM2+EM2+IM2;
    vacc_sterile2 = trapz(P.v*(1-P.z).*SH2,1)*P.da;
    vacc_blood2 = trapz(P.v*P.z.*SH2,1)*P.da;
    
    t = [t;t2(2:end)];
    SH = [SH,SH2(:,2:end)]; EH = [EH,EH2(:,2:end)]; DH = [DH, DH2(:,2:end)]; AH = [AH, AH2(:,2:end)]; VH = [VH, VH2(:,2:end)]; UH = [UH, UH2(:,2:end)]; MH = [MH, MH2(:,2:end)];
    SM = [SM, SM2(2:end)]; EM = [EM, EM2(2:end)]; IM = [IM, IM2(2:end)]; NM = [NM, NM2(2:end)];
    Cm = [Cm, Cm2(:,2:end)]; Cac = [Cac, Cac2(:,2:end)]; Cv = [Cv, Cv2(:,2:end)]; Ctot = [Ctot,Ctot2(:,2:end)];
    PH = SH+EH+DH+AH+VH+UH;
    PH_final = PH(:,end); % total human at age a, t = n
    NH = [NH, NH2(2:end)];
    vacc_sterile = [vacc_sterile, vacc_sterile2(2:end)];
    vacc_blood = [vacc_blood, vacc_blood2(2:end)];
    
    % calculate QOIs
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
    cases_total = trapz(cases)*P.dt; % total cases in cohort this year    
    cum_vacc = cumsum(vacc_sterile)*dt;
    vacc_total = cum_vacc(end); % total # vacc doeses this year
    pop = trapz(PH(ind1:ind2,:),1)*P.da;
    cases_pp_py = cases_total/mean(pop);% should be around ~3 for Nanoro seasonlity study
    keyboard
    save(['Results/season_vacc_',num2str(t0_list(it)),'.mat'],'cases_total','vacc_total','cases_pp_py')
end

%% comparison vacc starting time
cases_baseline = 6.545272849130029e+05;
cases_pp_py_baseline = 2.925170112290933;
for it = 1:length(t0_list)
    load(['Results/season_vacc_',num2str(t0_list(it)),'.mat'],'cases_total','vacc_total','cases_pp_py')
    cases_per_vacc(it) = (cases_baseline-cases_total)/vacc_total; % # of cases prevented per vacc 
    cases_pp_py_list(it) = cases_pp_py;
end
load(['Results/season_vacc_yearlong.mat'],'cases_total','vacc_total','cases_pp_py')
cases_per_vacc_yearlong = (cases_baseline-cases_total)/vacc_total; 
cases_pp_py_yearlong = cases_pp_py;

figure_setups; 
yyaxis left
h1 = plot(t0_list,cases_per_vacc,'--','DisplayName','cases prevented per vacc (seasonal)');
hold on
h2 = plot(t0_list,cases_per_vacc_yearlong*ones(size(t0_list)),'-','DisplayName','cases prevented per vacc (yearlong)');
ylim([0 1])
yyaxis right
h3 = plot(t0_list,cases_pp_py_list,'--','DisplayName','cases pp py (seasonal)');
hold on
h4 = plot(t0_list,cases_pp_py_baseline*ones(size(t0_list)),':','DisplayName','cases pp py (no vacc)');
h5 = plot(t0_list,cases_pp_py_yearlong*ones(size(t0_list)),'-','DisplayName','cases pp py (yearlong)');
ylim([2 3])
legend;
xlabel('vacc start month')