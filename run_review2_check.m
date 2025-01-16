%% generate SA results for reviewer 2 comment (now figure A2)

close all
clearvars
clc
format long
global P

% Directory for outputting SA results (sample, data matrices, plots)
% if data is available, the script will load results in the folder;
% otherwise, it will generate new results (could time and storage consuming)
direc = 'Data_SA/Results_local_SA_vac/SA_4POI_PRCC/';
flag_save = 1; % flag for saving the results or not (Note: it will overwrite previous results in the folder)

% numerical config
tfinal = 3*365; % time for integration beyond EE (e.g. vaccination)
age_max = 100*365; % max ages in days
P.age_max = age_max;
dt = 1; % time/age step size in days, default = 5;
da = dt;
t = (0:dt:tfinal)';
nt = length(t);
a = (0:da:age_max)';
na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t; P.tfinal = tfinal;

% SA setting

% code will calculate PRCC results for all the quantities below but only
% plotting a subset of this list. To modify the plotting, turn on
% "Size_QOI_plot" in the for loop of last section
lQ = {'EE-death','EE-death-00-02','EE-DA','EE-DA-00-02'}; %'EE-death-02-10','EE-death-10+','EE-DA-02-10','EE-DA-10+'
Size_QOI = length(lQ); % length of the QOI.
time_points = length(t); % default # time_points = at tfinal, unless if wants to check QOI at particular time points
lP_list = {'v0','w','etas'};
% lP_list = {'v0','w','etas','muM','betaM'};
% lP_list = {'rA','rD','cS','cA','cU','psis2','psir2','dac','uc','betaD','betaA','v0','w','etas'};
% lP_list = {'rA','rD','cS','cA','cU','psis2','psir2','dac','uc','muM','betaM','betaD','betaA','v0','w','etas'};
lP_list{end+1} = 'dummy'; % add dummy to the POIs
Malaria_parameters_baseline;
Malaria_parameters_transform_SA;
Malaria_parameters_transform_SA_once;
pmin = NaN(length(lP_list),1); 
pmax = pmin; 
pmean = pmin;
for iP = 1:length(lP_list)
    pmin(iP,1) = P.([lP_list{iP},'_lower']);
    pmax(iP,1) = P.([lP_list{iP},'_upper']);
    pmean(iP,1) = P.([lP_list{iP}]);
    if strcmp(lP_list{iP},'w')
        % POI = w wanning rate, sample the 1/w = average period instead.
        pmin(iP,1) = 1/P.([lP_list{iP},'_upper']);
        pmax(iP,1) = 1/P.([lP_list{iP},'_lower']);
        pmean(iP,1) = 1/P.([lP_list{iP}]);
        index_w = iP; % index for the POI = w
    end
end
% PRCC config
NS = 1000; % number of samples, min = k+1, 100~1000
k = length(lP_list); % # of POIs
% Pre-allocation
Size_timepts = length(time_points); % # of time points to check QOI value;
Y = NaN(NS,Size_timepts,Size_QOI);  % For each model evaluation, the QOI has dimension [Size_timepts, Size_QOI]

%% Generate parameter samples, stored in matrix X
disp('load parameter samples...')
load([direc,'PRCC_sample_',num2str(NS),'_',num2str(k),'.mat'],'X')

%% check any NaN runs in Ymat
load([direc,'PRCC_result_Ymat_',num2str(NS),'_',num2str(k),'.mat'],'Y')
total = NS;
(sum(isnan(Y(:)))/length(lQ))/total % fraction of files that are interrupted
for ind = 1:NS % Loop through each parameter sample
    S = load([direc,'PRCC','_',num2str(ind),'.mat']);
    if length(fieldnames(S))==1
        disp(['okay to delete file ',num2str(ind),'?'])
        keyboard
        delete([direc,'PRCC','_',num2str(ind),'.mat']);
    end
end
