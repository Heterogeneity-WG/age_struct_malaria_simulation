%% Figure 3, global SA at quasi-endmic state using PRCC
% 'Data_SA/Results_local_SA_no vac/SA_22POI_PRCC/'

close all
clearvars
clc
format long
global P

% Directory for outputting SA results (sample, data matrices, plots)
% if data is available, the script will load results in the folder;
% otherwise, it will generate new results (could time and storage consuming)
direc = 'Data_SA/Results_local_SA_no vac/SA_22POI_PRCC/';
flag_save = 1; % flag for saving the results or not (Note: it will overwrite previous results in the folder)

% numerical config
tfinal = 3*365; % time for integration beyond EE (e.g. vaccination)
age_max = 100*365; % max ages in days
P.age_max = age_max;
dt = 5; % time/age step size in days, default = 5;
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
lQ = {'EE-death','EE-death-09-24','EE-death-02-10','EE-death-10+',...
    'EE-DA','EE-DA-09-24','EE-DA-02-10','EE-DA-10+'};
Size_QOI = length(lQ); % length of the QOI.
time_points = length(t); % default # time_points = at tfinal, unless if wants to check QOI at particular time points
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
    disp('generate parameter samples...');
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
        Q_val = QOI_value_SA(lQ,time_points,run_num,'PRCC',direc); % calculate QOI values
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
    lP_order = {'dac','rD','psir2','cS','uc','muM','betaM','cA','rhos2',...
        'psis2','rA','cE','betaD','rhor2','sigma','m','betaA','cD',...
        'phis2','phir2','cU','v0'};
else % ordering for time-series SA plots
    lP_order = {'muM','cS','betaM','psir2','uc','cA','dac','rhos2', 'psis2',...
        'cE','betaD','m','cD','betaA','rD','rhor2','sigma','phis2','rA','cU','phir2','v0','w','etas'};
end
[~,index] = ismember(lP_order,lP_list);
index = index';
index(index==0)=[];
PRCC = PRCC([index;k],:,:);
stat_p = stat_p(index,:,:,:);
lP_list = lP_list([index;k]);

%% PRCC plot
X = categorical(lP_list);
X = reordercats(X,lP_list);
palpha = 0.05; % alpha for t-test
QOI_plot = 1:length(lQ);
Size_QOI_plot = length(QOI_plot);
[lP_list_name,lQ,lQ_title] = SA_output_formatting(lP_list,lQ,1);

% order the xticks by abs(PRCC index) based on QOI = malaria deaths
%temp_sort = abs(PRCC(:,1,QOI_plot(1)));
%[~,Isort] = sort(temp_sort,'descend');

for iQOI = [1 5] % 1:Size_QOI_plot
    figure_setups;
    hold on;
    b = bar(X,PRCC(:,1,QOI_plot(iQOI))); % need to reorder PRCC also
    ylim([-1.1 1.1])
    xtips = b.XEndPoints;
    ytips = b.YEndPoints;
    ytips(PRCC(:,1,QOI_plot(iQOI))<0) = ytips(PRCC(:,1,QOI_plot(iQOI))<0)-0.25;
    labels = cell(1, length(lP_list_name));
    labels(stat_p(:,1,QOI_plot(iQOI))<palpha) = {'*'};
    text(xtips,ytips,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','fontsize',35)
    title(['QOI = ', lQ_title{QOI_plot(iQOI)}])
    ax = gca;
    ax.XAxis.FontSize = 26;
    xticklabels(lP_list_name);
    xtickangle(90)
    grid off
    ax=gca;
    % read out the position of the axis in the unit "characters"
    set(ax,'Units','characters'); temp_ax=get(ax,'Position');
    % this sets an 'a)' right at the top left of the axes
    if iQOI == 1
        text(ax,-12,temp_ax(end)+2,'(A)','Units','characters');
        save_string = strcat('fig3_PRCC_','A','.svg');
    else
        text(ax,-12,temp_ax(end)+2,'(B)','Units','characters');
        save_string = strcat('fig3_PRCC_','B','.svg');
    end
    saveas(gcf,save_string);
    if flag_save
        saveas(gcf,[direc,'PRCC_result_',num2str(NS),'_',num2str(k),'_',lQ{QOI_plot(iQOI)},'.eps'],'epsc');
    end
    %keyboard;
end


