%% global SA using PRCC

% close all
clearvars
clc
format long
global P

% numerical config
tfinal = 20*365; % time for "EE" approximation
tfinal_vac = 10*365; % time to achieve "new EE" w/ vaccination
age_max = 100*365; % max ages in days
P.age_max = age_max;
dt = 20; % time/age step size in days, default = 5;
da = dt;
t = (0:dt:tfinal+tfinal_vac)';
nt = length(t);
a = (0:da:age_max)';
na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t; P.tfinal = tfinal; P.tfinal_vac = tfinal_vac;

% SA setting
lQ = {'EE-D','EE-DA','EE-D-frac','EE-EIR',...
    'EE-D-02-10','EE-DA-02-10','EE-D-frac-02-10',...
    'EE-D-09-24','EE-DA-09-24','EE-D-frac-09-24',...
    'EE-death','EE-death-02-10','EE-death-09-24'};  
% lQ = {'EE-death','EE-death-02-10','EE-death-09-24'};  
Size_QOI = length(lQ); % length of the QOI. Default = 1, unless it is an age distribution, or wants to test multiple QOIs at once
time_points = 1; % default time_points = 1, unless if wants to check QOI at particular time points
% lP_list = { 'dac','cX','w','etas'};
lP_list = {'rA','rD','muM','sigma','betaM','betaD', 'betaA','dac','cX','phis2','phir2','rhos2','rhor2','psis2','psir2','w','etas'};
% 'rA','muM','sigma','betaM','betaD', 'betaA','dac','cX','phis2','phir2','rhos2','rhor2','psis2','psir2','w'
lP_list{end+1} = 'dummy'; % add dummy to the POIs
Malaria_parameters_baseline;
pmin = NaN(length(lP_list),1); pmax = pmin; pmean = pmin;
for iP = 1:length(lP_list)
    pmin(iP,1) = P.([lP_list{iP},'_lower']);
    pmax(iP,1) = P.([lP_list{iP},'_upper']);
    pmean(iP,1) = P.([lP_list{iP}]);
end

% PRCC config
NS = 1000; % number of samples, min = k+1, 100~1000
k = length(lP_list); % # of POIs
% Pre-allocation
Size_timepts = length(time_points); % # of time points to check QOI value;
Y = NaN(NS,Size_timepts,Size_QOI);  % For each model evaluation, the QOI has dimension [Size_timepts, Size_QOI]

%% Generate parameter samples, stored in matrix X
LHS_raw = lhsdesign(NS,k); % uniform random draw with LHS sampling in [0,1]
X = parameterdist(LHS_raw,pmax,pmin,pmean,1,'triangular'); % 'unif' 'triangular'
tic
%% model evaluations
% common calculations/quantities that don't be impacted by the POIs
% calculate here once for efficiency
Malaria_parameters_baseline;
gH_fun = @(age) (2.*P.cc.*normpdf((age./365-P.zz)./P.ww).*normcdf(P.alpha.*(age./365-P.zz)./P.ww)./P.ww)./365/2;
gH =  gH_fun(a); % human fertility rate
muH =  P.b0 + P.b1*exp(-P.b2*a/365) + P.b3*exp(P.b4*a/365); % natural human mortality rate
muH_fun = @(age) P.b0 + P.b1*exp(-P.b2*age/365) + P.b3*exp(P.b4*age/365);
muH = muH/365;
muH_int_fun = @(age) (age./365).*P.b0 + (P.b1./P.b2).*(1-exp(-P.b2.*age./365)) + (P.b3./P.b4).*(-1+exp(P.b4.*age./365));
muD =  P.b0D + P.b1D*exp(-P.b2D*a/365) + P.b3D*exp(P.b4D*a/365); % natural human mortality rate
muD_fun = @(age) P.b0D + P.b1D*exp(-P.b2D*age/365) + P.b3D*exp(P.b4D*age/365);
muD = muD/365;
muD_int_fun = @(age) (age./365).*P.b0D + (P.b1D./P.b2D).*(1-exp(-P.b2D.*age./365)) + (P.b3D./P.b4D).*(-1+exp(P.b4D.*age./365));
v_fun = @(age) P.v0.*ones(size(age));
v = v_fun(a);
P.muH = muH;
P.muH_fun = muH_fun;
P.muD = muD;
P.muD_fun = muD_fun;
P.gH = gH;
P.v = v;
P.v_fun = v_fun;
P.gH_fun = gH_fun;
P.muH_int_fun = muH_int_fun;
P.muH_int = muH_int_fun(a);
P.muD_int_fun = muD_int_fun;
P.muD_int = muD_int_fun(a);
% Update the fertility and stable age dist. if balanced option is selected
if P.balance_fertility == 1
    balance_fertility;
end

if P.balance_mortality == 1
    balance_mortality;
end
find_stable_age; % depends on gH and muH

for run_num = 1:NS % Loop through each parameter sample
    if mod(run_num,NS/10)==0; disp([num2str(run_num/NS*100,3),' %']); end % display progress
    Malaria_parameters_baseline; % reset parameters to baseline
    for iP = 1:length(lP_list) % update parameter with sampled values
        P.(lP_list{iP}) = X(run_num,iP);
    end
    Malaria_parameters_transform_SA; % update dependent parameters
    Q_val = QOI_value_SA(lQ,time_points); % calculate QOI values
    Y(run_num,:,:) = Q_val;
end

%% PRCC on output matrix Y
PRCC = NaN(k,Size_timepts,Size_QOI); stat_p = PRCC;
for itime = 1:Size_timepts
    for iQOI = 1:Size_QOI
        [rho,p] = partialcorr([X Y(:,itime,iQOI)],'type','Spearman');
        PRCC(:,itime,iQOI) = rho(1:end-1,end); % correlations between parameters and QOI
        stat_p(:,itime,iQOI) = p(1:end-1,end); % associated p-value
    end
end
% save(['Results/vaccine_no/PRCC_result_',num2str(NS),'_',num2str(k),'.mat'],'PRCC','stat_p','lP_list','lQ')

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

%% Sorting (using order from eFAST)
load(['Results/vaccine_yes/PRCC_result_',num2str(NS),'_',num2str(k),'.mat'],'PRCC','stat_p','lP_list','lQ')
lP_eFAST = {'rD','cX','dac','psir2','muM','betaM', 'rhor2','betaD','rhos2','rA','sigma','psis2','phir2','phis2','betaA','w','etas'}; 
[~,index] = ismember(lP_eFAST,lP_list); index = index';
% [~,index] = sort(abs(PRCC(1:end-1,1,11)),'descend'); % sort using QOI #11, sort all the POIs except dummy
PRCC = PRCC([index;k],:,:); stat_p = stat_p(index,:,:,:);
lP_list = lP_list([index;k]);

% PRCC bar plot
palpha = 0.05; % alpha for t-test
QOI_plot = [2,6,9]; % [11:13] deaths; [2,6,9]: total infected
Size_QOI_plot = length(QOI_plot);
X = categorical(lP_list);
X = reordercats(X,lP_list);
for itime = 1:Size_timepts
    for iQOI = 1:Size_QOI_plot
        figure_setups; hold on
%         subplot(Size_timepts,Size_QOI_plot,(itime-1)*Size_QOI_plot+iQOI)
        b = bar(X,PRCC(:,itime,QOI_plot(iQOI)),'FaceAlpha',0.5);
        ylim([-1.1 1.1])
        xtips = b.XEndPoints;
        ytips = b.YEndPoints;
        labels = cell(1, length(lP_list));
        labels(stat_p(:,itime,QOI_plot(iQOI))<palpha) = {'*'};
%         for j = 1:k
%             labels{j} = sprintf('p = %3.1e',stat_p(j,itime,QOI_plot(iQOI)));
%         end
        text(xtips,ytips,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom','fontsize',14)
        title(['QOI=', lQ{QOI_plot(iQOI)}])
        saveas(gcf,['Results/vaccine_yes/PRCC_result_',num2str(NS),'_',num2str(k),'_',lQ{QOI_plot(iQOI)},'.png'])
    end
end

