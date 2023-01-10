clearvars
close all


%Load Parameters

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

%% Parameters to Vary

P.etas = 0.73;
P.etab = 0.73;
P.years=5;

param1_value=P.etas;
param2_value=P.etab;
param3_value=P.years;


Malaria_parameters_transform; 
Malaria_parameters_transform_vac;

%% Run Model

% initial condition 'EE' - numerical EE
[SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('EE_reset');

% time evolution - initial run
tfinal = P.years*365; t = (0:dt:tfinal)'; nt = length(t);
P.nt = nt;  P.t = t;
[SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(da,na,tfinal,...
    SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);

initial.SH=SH(:,1); initial.EH=EH(:,1); initial.DH=DH(:,1); initial.AH=AH(:,1); initial.VH=VH(:,1); initial.UH=UH(:,1);
initial.SM=SM(:,1); initial.EM=EM(:,1); initial.IM=IM(:,1);
initial.Cm=Cm(:,1); initial.Cac=Cac(:,1); initial.Cv=Cv(:,1); initial.Ctot=Ctot(:,1);
initial.MH=MH(:,1);

% lQ = 'EE-D-frac';  % R0 RHM RMH EE-EIR EE-EDA EE-infected EE-D-frac
lQ = {'EE-D','EE-DA'};%,'EE-D-frac','EE-EIR',...
%     'EE-D-02-10','EE-DA-02-10','EE-D-frac-02-10',...
%     'EE-D-09-24','EE-DA-09-24','EE-D-frac-09-24',...
%     'EE-death','EE-death-02-10','EE-death-09-24'};  


%% Output to Record

soln.SH=SH; soln.EH=EH; soln.DH=DH; soln.AH=AH; soln.VH=VH; soln.UH=UH;
soln.SM=SM; soln.EM=EM; soln.IM=IM;
soln.Cm=Cm; soln.Cac=Cac; soln.Cv=Cv; soln.Ctot=Ctot;
soln.MH=MH;

Qval=QOI_save(lQ,soln);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preparing output files
% Delete files for output data if they exist

% Change parameter values from decimals to strings --> 0.25=0p25
clearvars param1_value_print
if mod(param1_value,1)
    temp = num2str(param1_value);
    for i = 1:length(temp)
        if temp(i) == '.'
            param1_value_print(i) = 'p';
        else
            param1_value_print(i) = temp(i);
        end
    end
else
    param1_value_print = num2str(param1_value);
end

clearvars param2_value_print
if mod(param2_value,1)
    temp = num2str(param2_value);
    for i = 1:length(temp)
        if temp(i) == '.'
            param2_value_print(i) = 'p';
        else
            param2_value_print(i) = temp(i);
        end
    end
else
    param2_value_print = num2str(param2_value);
end


%Delete any previous copies
filename = ['Output/Output_etas_',param1_value_print,'_etab_',param2_value_print,'.mat'];
if exist(filename, 'file')==2
  delete(filename);
end

%Write the output to the file - where Output is everthing we want in the
%file
save([filename],'P','initial','Qval');


%% To save a figure

% str = sprintf('MC_Ebola NumDeaths PeakInf TimePeak: R_0 = %.2f, I1_0=%i, I2_0=%i',R0, I1_0, I2_0);
% set(gcf,'Name',str,'NumberTitle','off')
% saveas(gcf,strcat(str,'.fig'))
