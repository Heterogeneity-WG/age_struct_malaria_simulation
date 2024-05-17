%% Figure A3 (supplement for Figure 5, selected PRCC curves)
% global SA time series at quasi-endmic state using using PRCC
% 'Data_SA/Results_local_SA_time/SA_22POI_PRCC/'

close all
clearvars
clc
format long
global P

% Directory for outputting SA results (sample, data matrices, plots)
% if data is available, the script will load results in the folder;
% otherwise, it will generate new results (could time and storage consuming)
direc = 'Data_SA/Results_local_SA_time/SA_22POI_PRCC/';
% direc = 'C:/Users/lbs042/Downloads/Results_local_SA/SA_22POI_PRCC/';
flag_save = 1; % flag for saving the results or not (Note: it will overwrite previous results in the folder)

% numerical config
tfinal = 3*365; % time for integration beyond EE (e.g. vaccination)
age_max = 100*365; % max ages in days
P.age_max = age_max;
dt = 1; % time/age step size in days, default = 1 for SA;
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
lQ = {'EE-A','EE-D','EE-DA','EE-EIR','EE-death-rate'};
Size_QOI = length(lQ); % length of the QOI.
time_points = 1:nt; % default # time_points = at tfinal, unless if wants to check QOI at particular time points
% NOTE: code will calculate results for all the parameters below but only
% plotting a subset of this list. To modify the plotting, turn on
% "POI_index" in the for loop of last section
lP_list = {'cS','cE','cA','cD','cU','phis2','phir2','rhos2','rhor2','psis2','psir2','dac','uc','m',...
    'rA','rD','muM','sigma','betaM','betaD', 'betaA'};
lP_list{end+1} = 'dummy'; % add dummy to the POIs
Malaria_parameters_baseline;
Malaria_parameters_baseline_Nanoro; % SA based on Nanoro climate profile
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
% PRCC config
NS = 1000; % number of samples, min = k+1, 100~1000
k = length(lP_list); % # of POIs
% Pre-allocation
Size_timepts = length(time_points); % # of time points to check QOI value;
Y = NaN(NS,Size_timepts,Size_QOI);  % For each model evaluation, the QOI has dimension [Size_timepts, Size_QOI]

%% Generate parameter samples, stored in matrix X
if ~exist([direc,'PRCC_sample_',num2str(NS),'_',num2str(k),'.mat'],'file')
    disp('generate parameter samples...')
    LHS_raw = lhsdesign(NS,k); % uniform random draw with LHS sampling in [0,1]
    X = parameterdist(LHS_raw,pmax,pmin,pmean,1,'triangular'); % 'unif' 'triangular'
    if flag_save; save([direc,'PRCC_sample_',num2str(NS),'_',num2str(k),'.mat'],'X'); end
    if flag_save; save([direc,'parameters_P_baseline.mat'],'P'); end
else
    disp('load parameter samples...')
    load([direc,'PRCC_sample_',num2str(NS),'_',num2str(k),'.mat'],'X')
end

%% model evaluations
tic
if exist([direc,'PRCC_result_Ymat_',num2str(NS),'_',num2str(k),'.mat'],'file')
    disp('load model evaluation results ...')
    load([direc,'PRCC_result_Ymat_',num2str(NS),'_',num2str(k),'.mat'],'Y')
else
    for run_num = 1:NS % Loop through each parameter sample
        if mod(run_num,NS/10)==0; disp([num2str(run_num/NS*100,3),' %']); end % display progress
        Malaria_parameters_baseline; % reset parameters to baseline
        for iP = 1:length(lP_list) % update parameter with sampled values
            P.(lP_list{iP}) = X(run_num,iP);
        end
        if ismember('w',lP_list)
            % POI = w wanning rate: X stores the sampled value 1/w = period, model uses w = rate.
            P.(lP_list{index_w}) = 1/X(run_num,index_w);
        end
        Malaria_parameters_transform_SA; % update dependent parameters
        Q_val = QOI_value_SA_2(lQ,time_points,run_num,'PRCC',direc); % calculate QOI values
        Y(run_num,:,:) = Q_val;
    end
    % Y(NS,Size_timepts,Size_QOI,length(pmin),NR)
    if flag_save; save([direc,'PRCC_result_Ymat_',num2str(NS),'_',num2str(k),'.mat'],'Y'); end
end
toc

%% PRCC on output matrix Y
PRCC = NaN(k,Size_timepts,Size_QOI); stat_p = PRCC;
for itime = 1:Size_timepts
    for iQOI = 1:Size_QOI
        [rho,p] = partialcorr([X Y(:,itime,iQOI)],'type','Spearman');
        PRCC(:,itime,iQOI) = rho(1:end-1,end); % correlations between parameters and QOI
        stat_p(:,itime,iQOI) = p(1:end-1,end); % associated p-value
    end
end

%% Sorting
if Size_timepts==1
    lP_order = {'dac','rD','cS','psir2','uc','muM','cA','rhos2','betaM','rA','psis2',...
        'cE','betaD','m','cD','betaA','rhor2','sigma','phis2','cU','phir2','v0','w','etas'};
else % ordering for time-series SA plots
    lP_order = {'muM','cS','betaM','psir2','uc','cA','dac','rhos2', 'psis2',...
        'cE','betaD','m','cD','betaA','rD','rhor2','sigma','phis2','rA','cU','phir2','v0','w','etas'};
end
[~,index] = ismember(lP_order,lP_list); index = index';
index(index==0)=[];
PRCC = PRCC([index;k],:,:); stat_p = stat_p(index,:,:,:);
lP_list = lP_list([index;k]);

%% PRCC plot
X = categorical(lP_list);
X = reordercats(X,lP_list);
palpha = 0.05; % alpha for t-test
QOI_plot = 1:length(lQ);
Size_QOI_plot = length(QOI_plot);
[lP_list_name,lQ,lQ_title] = SA_output_formatting(lP_list,lQ,1);

for iQOI = 1:Size_QOI_plot
    mylinestyles = ["-"; "--"; ":"; "-."];
    mycolors = lines(5);
    figure_setups_4; hold on;
    xlabel('years')
    title(['QOI = ', lQ_title{QOI_plot(iQOI)}])
    ylim([-1 1])
    POI_index = [1 3 7 15 19]; %1:k
    for iPOI = 1:length(POI_index)
        PRCC_vec = PRCC(POI_index(iPOI),:,QOI_plot(iQOI))';
        time_pts = t(time_points)/365;
        PRCC_fun = csape(time_pts,PRCC_vec,'not-a-knot');
        time_pts_fine = linspace(time_pts(1),time_pts(end),nt*3);
        PRCC_vec_fine = ppval(PRCC_fun,time_pts_fine);
        color_style = mod(POI_index(iPOI),5); if color_style==0; color_style=5; end
        line_style = ceil(POI_index(iPOI)/5);
        plot(time_pts_fine,PRCC_vec_fine,'color', mycolors(color_style,:),'LineStyle',mylinestyles(line_style),'DisplayName',lP_list_name{POI_index(iPOI)})
    end
    legend('Orientation','horizontal','NumColumns',10)
    plot(time_pts_fine,zeros(size(time_pts_fine)),'k-','DisplayName','')
    if flag_save; saveas(gcf,[direc,'PRCC_result_',num2str(NS),'_',num2str(k),'_',lQ{QOI_plot(iQOI)},'.eps'],'epsc'); end
end
