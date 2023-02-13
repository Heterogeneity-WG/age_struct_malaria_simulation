clearvars
close all;  
clc; 
format long;
global P lP

%% numerical config
age_max = 100*365; % max ages in days
P.age_max = age_max;
dt = 20; % time/age step size in days, default = 5;
da = dt;
a = (0:da:age_max)';
na = length(a);

P.a = a;
P.na = na;
P.da = da;
P.dt = dt;

Malaria_parameters_baseline;
Malaria_parameters_transform;
Malaria_parameters_transform_vac;

lQ = {'year5-D-process'}; % metric to compare
tfinal = 5*365;
lP_list = {'v0'}; % etas - sterile, etab - blood stage

%% initial condition - numerical EE with reset, no vaccine
P.v0 = 0;
Malaria_parameters_transform;
Malaria_parameters_transform_vac;
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');
% soln = solution_pack(SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
% QOI_initial = QOI_save(lQ,soln);

%% baseline run
disp('baseline run ---');
P.z = 0; % RTS,S only
P.v0 = 15; 
Malaria_parameters_transform;
Malaria_parameters_transform_vac;

if strcmp(lQ{1}(1:5),'year5')
    [SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = ...
        age_structured_Malaria_vac(P.da,P.na,tfinal, SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
    soln = solution_pack(SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH);
end
% P_baseline = P.(lP);
Q_baseline = QOI_save(lQ,soln);

%% comparison runs
P.z = 1; % blood stage only
P.v0_lower = 1;
P.v0_upper = 100;
Malaria_parameters_transform;
Malaria_parameters_transform_vac;

tic
for iP = 1:length(lP_list)
    lP = lP_list{iP};   
    
    %% extended SA
    P_lower = P.([lP,'_lower']);
    P_upper = P.([lP,'_upper']);
    ngrid = 5;
    
    % allocation
    P_vals = linspace(P_lower,P_upper,ngrid)';
    Q_vals = NaN(length(P_vals),length(Q_baseline));
    
    for i=1:ngrid
        display(['I am working on simulation ', num2str(i)])
        P.(lP) = P_vals(i);
        Malaria_parameters_transform;
        Malaria_parameters_transform_vac;
        if strcmp(lQ{1}(1:5),'year5')
            [SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = ...
                age_structured_Malaria_vac(P.da,P.na,tfinal, SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
            soln = solution_pack(SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH);
        end
        Q_vals(i,:) = QOI_save(lQ,soln);
    end
    
end
toc

%% plotting
figure_setups; hold on
t = (0:dt:tfinal)';
plot(t/365,Q_vals);
plot(t/365,Q_baseline,'r','MarkerSize',20);
xlabel('time (years)')
ylabel(lQ)
title (['vary ', lP])
legendStrings = lP + " = " + string(P_vals);
legend(legendStrings)