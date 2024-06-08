%% Figure 5, global SA time series using eFAST (selected eFAST curves)
% 'Data_SA/Results_local_SA_time/SA_22POI_eFAST/'

% adapted from
% Dr. Denise Kirschner - University of Michigan http://malthus.micro.med.umich.edu/lab/usadata/
% Output:
% SI[POI,timepoints,QOI] : first order sensitivity indices
% STI[POI,timepoints,QOI] : total effect sensitivity indices
% CV[POI,timepoints,QOI]

close all
clearvars
clc
format long
global P

% Directory for outputting SA results (sample, data matrices, plots)
% if data is available, the script will load results in the folder;
% otherwise, it will generate new results (could time and storage consuming)
% direc = 'D:/Results_local_SA/SA_22POI_eFAST/';
% direc = 'C:/Users/lbs042/Downloads/Results_local_SA/SA_22POI_eFAST/';
direc = 'Data_SA/Results_local_SA_no vac/SA_22POI_eFAST/dt1/';
flag_save = 0; % flag for saving the results or not (Note: it will overwrite previous results in the folder)

% numerical config
tfinal = 3*365; % time for integration beyond EE (e.g. vaccination)
age_max = 100*365; % max ages in days
P.age_max = age_max;
dt = 1; % time/age step size in days, default = 20 for SA (=5 for individual simulation);
da = dt;
t = (0:dt:tfinal)';
nt = length(t);
a = (0:da:age_max)';
na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t; P.tfinal = tfinal;

% SA setting

% NOTE: code will calculate results for all the quantities below but only
% plotting a subset of this list. To modify the plotting, turn on
% "Size_QOI_plot" in the for loop of last section
lQ = {'EE-death-02-10','EE-death-10+','EE-DA-02-10','EE-DA-10+'};
Size_QOI = length(lQ); % length of the QOI.
time_points = length(t); % default # time_points = at tfinal, unless if wants to check QOI at particular time points
% NOTE: code will calculate results for all the parameters below but only
% plotting a subset of this list. To modify the plotting, turn on
% "POI_index" in the for loop of last section
lP_list = {'cS','cE','cA','cD','cU','phis2','phir2','rhos2','rhor2','psis2','psir2','dac','uc','m',...
    'rA','rD','muM','sigma','betaM','betaD', 'betaA'};
lP_list{end+1} = 'dummy'; % add dummy to the POIs
Malaria_parameters_baseline;
Malaria_parameters_transform_SA;
Malaria_parameters_transform_SA_once;
pmin = NaN(length(lP_list),1); pmax = pmin; pmean = pmin;
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

% eFAST config
NR = 5;    % # of search curves - Resampling - keep the value
k = length(lP_list); % # of POIs + dummy parameter, keep it in the range 5~11
NS = 257;   % # of samples per search curve - keep 2^n+1, min = 65
wantedN = NS*k*NR; % wanted no. of sample points
MI = 4; %: maximum number of fourier coefficients that may be retained in calculating the partial variances without interferences between the assigned frequencies
% Computation of the frequency for the group of interest OMi and the # of sample points NS (here N=NS)
OMi = floor(((wantedN/NR)-1)/(2*MI)/k); % >=8
NS = 2*MI*OMi+1;
if(NS*NR < 65)
    fprintf(['Error: sample size must be >= 65 per factor.\n']);
    return;
end
% Pre-allocation
Size_timepts = length(time_points); % # of time points to check QOI value;
%% load sample
disp('load parameter samples...')
load([direc,'eFAST_samples_',num2str(NS),'_',num2str(k),'_',num2str(NR),'.mat'],'X')
%% check any NaN runs in Ymat
load([direc,'eFAST_result_Ymat_',num2str(NS),'_',num2str(k),'.mat'],'Y');
total = NS*NR*k;
(total-sum(isnan(Y(:)))/length(lQ))/total % fraction of files that are interrupted
keyboard
for i=1:k % Loop over POIs, including dummy
    for L=1:NR  % Loop over the NR search curves
        for run_num=1:NS % Loop through each parameter sample
            ind = NS*NR*(i-1)+NS*(L-1)+run_num;
            if mod(ind,round(NS*NR*k/100))==0
                disp([num2str(ind/(NS*NR*k)*100,3),' %']); % display progress
            end
            S = load([direc,'eFAST','_',num2str(ind),'.mat']);
            if length(fieldnames(S))==1
                disp(['okay to delete file ',num2str(ind),'?'])
                keyboard
                delete([direc,'eFAST','_',num2str(ind),'.mat']);
            end
        end
    end
end


