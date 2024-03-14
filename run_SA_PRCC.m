%% global SA using PRCC

close all
clearvars
clc
format long
global P

% numerical config
tfinal = 3*365; % time for integration beyond EE (e.g. vaccination)
age_max = 100*365; % max ages in days
P.age_max = age_max;
dt = 20; % time/age step size in days, default = 5;
da = dt;
t = (0:dt:tfinal)';
nt = length(t);
a = (0:da:age_max)';
na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t; P.tfinal = tfinal; 

% SA setting 
lQ = {'EE-death','EE-death-09-24','EE-death-02-10','EE-death-10+',...
      'EE-DA','EE-DA-09-24','EE-DA-02-10','EE-DA-10+'};  
% lQ = {'EE-A','EE-D','EE-DA','EE-EIR','EE-death-rate'};
Size_QOI = length(lQ); % length of the QOI. Default = 1, unless it is an age distribution, or wants to test multiple QOIs at once
time_points = length(t); % default # time_points = at tfinal, unless if wants to check QOI at particular time points
% time_points = 1:nt;  % for time-series SA index
% lP_list = {'cS','cE','cA','cD','cU','phis2','phir2','rhos2','rhor2','psis2','psir2','dac','uc','m',...
%     'rA','rD','muM','sigma','betaM','betaD', 'betaA'};
lP_list = {'rA','rD','cS','cA','cU','psis2','psir2','dac','uc','muM','betaM','betaD','betaA','v0','w','etas'};
lP_list{end+1} = 'dummy'; % add dummy to the POIs
Malaria_parameters_baseline;
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
direc = 'D:/Results_local_SA_vac/SA_17POI_PRCC/';
% direc = ['C:/Users/lbs042/Downloads/Results_local_SA/SA_22POI_PRCC/'];
if ~exist([direc,'PRCC_sample_',num2str(NS),'_',num2str(k),'.mat'],'file')
    disp('generate parameter samples...')
    LHS_raw = lhsdesign(NS,k); % uniform random draw with LHS sampling in [0,1]
    X = parameterdist(LHS_raw,pmax,pmin,pmean,1,'triangular'); % 'unif' 'triangular'
    save([direc,'PRCC_sample_',num2str(NS),'_',num2str(k),'.mat'],'X')
else
    disp('load parameter samples...')
    load([direc,'PRCC_sample_',num2str(NS),'_',num2str(k),'.mat'],'X')
end
tic
%% model evaluations
% common calculations/quantities that don't be impacted by the POIs
% calculate here once for efficiency
Malaria_parameters_transform_SA;
P.zeta_fun = @(a) 1-0.85*exp(-a/8/365); %ones(size(a));
P.zeta = P.zeta_fun(a);
P.ss_S = @(t) P.ss_S0*(P.ss_c+P.ss_v*(1-P.ss_c)*((1+cos(2*pi*((t+P.ss_t0)/365-P.ss_u1)))./2).^P.ss_k1+...
    (1-P.ss_v)*(1-P.ss_c)*((1+cos(2*pi*((t+P.ss_t0)/365-P.ss_u2)))./2).^P.ss_k2);
P.gM_fun = @(t) 0.5*P.NN*P.ss_S(t);
gH_fun = @(age) (2.*P.cc.*normpdf((age./365-P.zz)./P.ww).*normcdf(P.alpha.*(age./365-P.zz)./P.ww)./P.ww)./365/2;
gH =  gH_fun(a); % human fertility rate
muH =  P.b0 + P.b1*exp(-P.b2*a/365) + P.b3*exp(P.b4*a/365); % natural human mortality rate
muH_fun = @(age) P.b0 + P.b1*exp(-P.b2*age/365) + P.b3*exp(P.b4*age/365);
muH = muH/365;
muH_int_fun = @(age) (age./365).*P.b0 + (P.b1./P.b2).*(1-exp(-P.b2.*age./365)) + (P.b3./P.b4).*(-1+exp(P.b4.*age./365));
P.muH_int = muH_int_fun(a);

a74 = 74*365;
temp_muD = P.b0D + P.b1D*exp(-P.b2D*a74/365) + P.b3D*exp(P.b4D*a74/365);
muD =  (P.b0D + P.b1D*exp(-P.b2D*a/365) + P.b3D*exp(P.b4D*a/365)).*(a/365<=74)...
    + temp_muD.*(exp(-P.c0D*(a-74*365)./365)).*(a/365>74); 
muD_fun = @(age) (P.b0D + P.b1D*exp(-P.b2D*age/365) + P.b3D*exp(P.b4D*age/365)).*(age/365<=74)...
    + temp_muD.*(exp(-P.c0D*(age/365-74))).*(age/365>74); 
muD = muD/365;
P.muH = muH;
P.muH_fun = muH_fun;
P.muH_int_fun = muH_int_fun;

P.muD = muD;
P.muD_fun = muD_fun;
P.gH = gH;
P.gH_fun = gH_fun;
find_stable_age; % depends on gH and muH
save([direc,'parameters_P_baseline.mat'],'P')

for run_num = 1:NS % Loop through each parameter sample
    if mod(run_num,NS/10)==0; disp([num2str(run_num/NS*100,3),' %']); end % display progress
    Malaria_parameters_baseline; % reset parameters to baseline
    for iP = 1:length(lP_list) % update parameter with sampled values
        P.(lP_list{iP}) = X(run_num,iP);
    end
    % POI = w wanning rate: X stores the sampled value 1/w = period, model uses w = rate.
    P.(lP_list{index_w}) = 1/X(run_num,index_w);
    Malaria_parameters_transform_SA; % update dependent parameters
    Q_val = QOI_value_SA(lQ,time_points,run_num,'PRCC',direc); % calculate QOI values
    Y(run_num,:,:) = Q_val;
end
% Y(NS,Size_timepts,Size_QOI,length(pmin),NR)
save([direc,'PRCC_result_Ymat_',num2str(NS),'_',num2str(k),'.mat'],'Y')
%% PRCC on output matrix Y
% load([direc,'PRCC_result_Ymat_',num2str(NS),'_',num2str(k),'.mat'],'Y')
PRCC = NaN(k,Size_timepts,Size_QOI); stat_p = PRCC;
for itime = 1:Size_timepts
    for iQOI = 1:Size_QOI
        [rho,p] = partialcorr([X Y(:,itime,iQOI)],'type','Spearman');
        PRCC(:,itime,iQOI) = rho(1:end-1,end); % correlations between parameters and QOI
        stat_p(:,itime,iQOI) = p(1:end-1,end); % associated p-value
    end
end
save([direc,'PRCC_result_',num2str(NS),'_',num2str(k),'.mat'],'PRCC','stat_p','lP_list','lQ')
toc

%% Plotting POIs vs. QOIs: check monotonic relationships
% for itime = 1:Size_timepts
%     figure_setups;
%     for j = 1:k
%         for iQOI = 1:Size_QOI
%             subplot(Size_QOI,k,(iQOI-1)*k+j)
%             plot(X(:,j),Y(:,itime,iQOI),'.')
%             xlabel(lP_list{j},'fontsize',14)
%             ylabel(lQ{iQOI},'fontsize',14)
%             set(gca,'fontsize',14)
%             xlim([pmin(j) pmax(j)])
%         end
%     end
%     sgtitle(['Check monotonic relationships (rows = QOIs, columns = POIs)'])
% %     sgtitle(['Check monotonic relationships at timepoint \#', num2str(itime),' (rows = QOIs, columns = parameters)'])
% end

%%
%% Sorting 
load(['Results/vaccine_no/PRCC_result_',num2str(NS),'_',num2str(k),'.mat'],'PRCC','stat_p','lP_list','lQ')
if Size_timepts==1
    lP_eFAST = {'dac','rD','cS','psir2','uc','muM','cA','rhos2','betaM','rA','psis2',...
    'cE','betaD','m','cD','betaA','rhor2','sigma','phis2','cU','phir2','v0','w','etas'};
else % ordering for time-series SA plots
    lP_eFAST = {'muM','cS','betaM','psir2','uc','cA','dac','rhos2', 'psis2',...
    'cE','betaD','m','cD','betaA','rD','rhor2','sigma','phis2','rA','cU','phir2','v0','w','etas'};
end
[~,index] = ismember(lP_eFAST,lP_list); index = index';
index(index==0)=[]; 
% [~,index] = sort(abs(PRCC(1:end-1,1,1)),'descend'); % sort using QOI #11, sort all the POIs except dummy
PRCC = PRCC([index;k],:,:); stat_p = stat_p(index,:,:,:);
lP_list = lP_list([index;k]);

%% PRCC plot 
X = categorical(lP_list);
X = reordercats(X,lP_list);
palpha = 0.05; % alpha for t-test
QOI_plot = 1:length(lQ); % [11:13] deaths; [2,6,9]: total infected
Size_QOI_plot = length(QOI_plot);
% lP_list_name = lP_list;
[lP_list_name,lQ,lQ_title] = SA_output_formatting(lP_list,lQ,1);
if Size_timepts == 1  % bar plot (for one time point)
    for iQOI = 1:Size_QOI_plot
        figure_setups; hold on
        % subplot(Size_timepts,Size_QOI_plot,(itime-1)*Size_QOI_plot+iQOI)
        b = bar(X,PRCC(:,1,QOI_plot(iQOI)));
        ylim([-1.1 1.1])
        xtips = b.XEndPoints;
        ytips = b.YEndPoints;
        ytips(PRCC(:,1,QOI_plot(iQOI))<0) = ytips(PRCC(:,1,QOI_plot(iQOI))<0)-0.15;
        labels = cell(1, length(lP_list_name));
        labels(stat_p(:,1,QOI_plot(iQOI))<palpha) = {'*'};
        %         for j = 1:k
        %             labels{j} = sprintf('p = %3.1e',stat_p(j,itime,QOI_plot(iQOI)));
        %         end
        text(xtips,ytips,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom','fontsize',14)
        title(['QOI = ', lQ_title{QOI_plot(iQOI)}])
        xticklabels(lP_list_name)
        grid off
        saveas(gcf,[direc,'PRCC_result_',num2str(NS),'_',num2str(k),'_',lQ{QOI_plot(iQOI)},'.eps'],'epsc')
    end
else
    for iQOI = 1:Size_QOI_plot  
        mylinestyles = ["-"; "--"; ":"; "-."];
        mycolors = lines(5);
        figure_setups_2; hold on; 
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
            % if iPOI<k
            %     labels = cell(size(PRCC_vec));
            %     % labels(stat_p(iPOI,:,QOI_plot(iQOI))'<palpha) = {'*'};
            %     % text(t(time_points)/365,PRCC_vec,labels,'HorizontalAlignment','center',...
            %     %     'VerticalAlignment','bottom','fontsize',14)
            % end
        end
        legend('Location','eastoutside')
        plot(time_pts_fine,zeros(size(time_pts_fine)),'k-','DisplayName','')
        saveas(gcf,[direc,'PRCC_result_',num2str(NS),'_',num2str(k),'_',lQ{QOI_plot(iQOI)},'.eps'],'epsc')
    end
end

