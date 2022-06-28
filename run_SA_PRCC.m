%% global SA using PRCC

close all
clear all
clc
format long
global P

%% numerical config
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

%% SA setting
lQ = {'EE-D','EE-DA'};  % R0 RHM RMH EE-EIR EE-D-frac EE-D EE-DA EE-death
% 'EE-D','EE-D-frac','EE-DA'
Size_QOI = length(lQ); % length of the QOI. Default = 1, unless it is an age distribution, or wants to test multiple QOIs at once
time_points = 1; % default time_points = 1, unless if wants to check QOI at particular time points
lP_list = {'muM','betaM'};
% 'rA','muM','sigma','betaM','betaD', 'betaA','dac','cX','phis2','phir2','rhos2','rhor2','psis2','psir2','w'
Malaria_parameters_baseline;
pmin = NaN(length(lP_list),1); pmax = pmin; pmean = pmin;
for iP = 1:length(lP_list)
    pmin(iP,1) = P.([lP_list{iP},'_lower']);
    pmax(iP,1) = P.([lP_list{iP},'_upper']);
    pmean(iP,1) = P.([lP_list{iP}]);
end

%% PRCC config
NS = 100; % number of samples
k = length(lP_list); % # of POIs
%% Pre-allocation
Size_timepts = length(time_points); % # of time points to check QOI value;
Y = NaN(NS,Size_timepts,Size_QOI);  % For each model evaluation, the QOI has dimension [Size_timepts, Size_QOI]

%% Generate parameter samples, stored in matrix X
LHS_raw = lhsdesign(NS,k); % uniform random draw with LHS sampling in [0,1]
X = parameterdist(LHS_raw,pmax,pmin,pmean,1,'unif'); % 'unif' 'triangular'
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
    disp([num2str(run_num/NS*100,3),' %']) % display progress
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
toc

%% Plotting POIs vs. QOIs: check monotonic relationships
for itime = 1:Size_timepts
    figure_setups;
    for j = 1:k
        for iQOI = 1:Size_QOI
            subplot(Size_QOI,k,(iQOI-1)*k+j)
            plot(X(:,j),Y(:,itime,iQOI),'.')
            xlabel(lP_list{j},'fontsize',14)
            ylabel(lQ{iQOI},'fontsize',14)
            set(gca,'fontsize',14)
            xlim([pmin(j) pmax(j)])
        end
    end
    sgtitle(['Check monotonic relationships (rows = QOIs, columns = POIs)'])
%     sgtitle(['Check monotonic relationships at timepoint \#', num2str(itime),' (rows = QOIs, columns = parameters)'])
end

%% PRCC bar plot
figure_setups;
for itime = 1:Size_timepts
    for iQOI = 1:Size_QOI
        subplot(Size_timepts,Size_QOI,(itime-1)*Size_QOI+iQOI)
        b = bar(PRCC(:,itime,iQOI),'FaceAlpha',0.5);
        set(gca,'xticklabel',lP_list,'fontsize',14)
        ylim([-1.1 1.1])
        xtips = b.XEndPoints;
        ytips = b.YEndPoints;
        for j = 1:k
            labels{j} = sprintf('p = %3.1e',stat_p(j,itime,iQOI));
        end
        text(xtips,ytips,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom','fontsize',14)
        title(['QOI=', lQ{iQOI}])
    end
end

