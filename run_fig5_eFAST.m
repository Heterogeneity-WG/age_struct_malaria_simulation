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
direc = 'Data_SA/Results_local_SA_time/SA_22POI_eFAST/';
flag_save = 0; % flag for saving the results or not (Note: it will overwrite previous results in the folder)

% numerical config
tfinal = 3*365; % time for integration beyond EE (e.g. vaccination)
age_max = 100*365; % max ages in days
P.age_max = age_max;
dt = 20; % time/age step size in days, default = 20 for SA (=5 for individual simulation);
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
lQ = {'EE-EIR','EE-DA','EE-death-rate','EE-A','EE-D'};
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
Y(NS,Size_timepts,Size_QOI,length(pmin),NR)=0;  % For each model evaluation, the QOI has dimension [Size_timepts*Size_QOI]
X(NS,k,k,NR) = 0;

%% Generate parameter samples, stored in matrix X
if ~exist([direc,'eFAST_samples_',num2str(NS),'_',num2str(k),'_',num2str(NR),'.mat'],'file')
    disp('generate parameter samples...')
    X = efast_gensamples(X,OMi,MI,pmin,pmax,pmean,'triangular'); % triangular distribution for POIs
    % X = efast_gensamples(X,OMi,MI,pmin,pmax,pmean,'unif'); % uniform distribution for POIs
    if flag_save; save([direc,'eFAST_samples_',num2str(NS),'_',num2str(k),'_',num2str(NR),'.mat'],'X'); end
    if flag_save; save([direc,'parameters_P_baseline.mat'],'P'); end
else
    disp('load parameter samples...')
    load([direc,'eFAST_samples_',num2str(NS),'_',num2str(k),'_',num2str(NR),'.mat'],'X')
end

%% model evaluations
tic
if exist([direc,'eFAST_result_Ymat_',num2str(NS),'_',num2str(k),'.mat'],'file')
    disp('load model evaluation results ...')
    load([direc,'eFAST_result_Ymat_',num2str(NS),'_',num2str(k),'.mat'],'Y')
else
    for i=1:k % Loop over POIs, including dummy
        for L=1:NR  % Loop over the NR search curves
            for run_num=1:NS % Loop through each parameter sample
                ind = NS*NR*(i-1)+NS*(L-1)+run_num;
                if mod(ind,round(NS*NR*k/100))==0
                    disp([num2str(ind/(NS*NR*k)*100,3),' %']); % display progress
                end
                % disp([num2str([i, run_num, L]), '//parameter run NR']) % keeps track of [parameter run NR]
                Malaria_parameters_baseline; % reset parameters to baseline
                for iP = 1:length(lP_list) % update parameter with sampled values
                    P.(lP_list{iP}) = X(run_num,iP,i,L);
                end
                if ismember('w',lP_list)
                    % POI = w wanning rate: X stores the sampled value 1/w = period, model uses w = rate.
                    P.(lP_list{index_w}) = 1/X(run_num,index_w);
                end
                Malaria_parameters_transform_SA; % update dependent parameters
                Q_val = QOI_value_SA(lQ,time_points,ind,'eFAST',direc); % calculate QOI values
                Y(run_num,:,:,i,L) = Q_val;
            end
        end
    end
    % Y(NS,Size_timepts,Size_QOI,length(pmin),NR)
    if flag_save; save([direc,'eFAST_result_Ymat_',num2str(NS),'_',num2str(k),'.mat'],'Y'); end
end

toc

%% eFAST on output matrix Y
var = 1:length(lQ); % index of QOIs to analyze (among Size_QOI) (default = 1)
palpha = 0.05; % alpha for t-test
[Si,Sti,rangeSi,rangeSti] = efast_sd(Y,OMi,MI,time_points,var);
[CVsi, CVsti] = CVmethod(Si, rangeSi,Sti,rangeSti,var); % Coeff. of Variance;
% SI[POI,timepoints,QOI]
s_struct = efast_ttest(Si,rangeSi,Sti,rangeSti,time_points,lP_list,var,lQ,palpha); % T-test on Si and STi

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
s_struct.Si = s_struct.Si([index;k],:,:); s_struct.p_Si = s_struct.p_Si(index,:,:,:); s_struct.rangeSi = s_struct.rangeSi([index;k],:,:,:);
s_struct.Sti = s_struct.Sti([index;k],:,:); s_struct.p_Sti = s_struct.p_Sti(index,:,:,:); s_struct.rangeSti = s_struct.rangeSti([index;k],:,:,:);
lP_list = lP_list([index;k]);

%% eFAST plot
X = categorical(lP_list);
X = reordercats(X,lP_list);
palpha = 0.05; % alpha for t-test
QOI_plot = 1:length(lQ);
Size_QOI_plot = length(QOI_plot);
[lP_list_name,lQ,lQ_title] = SA_output_formatting(lP_list,lQ,1);

for iQOI = 1:Size_QOI_plot
    mylinestyles = ["-"; "--"; ":"; "-."];
    mycolors = lines(5);
    figure_setups_2; hold on;
    xlabel('years')
    title(['QOI = ', lQ_title{QOI_plot(iQOI)}])
    ylim([0 1.1])
    POI_index = [1 3 7 15 19];
    for iPOI = 1:length(POI_index)
        Sti_vec = s_struct.Sti(POI_index(iPOI),:,QOI_plot(iQOI))';
        time_pts = t(time_points)/365;
        Sti_fun = csape(time_pts,Sti_vec,'not-a-knot');
        time_pts_fine = linspace(time_pts(1),time_pts(end),nt*3);
        Sti_vec_fine = ppval(Sti_fun,time_pts_fine);
        color_style = mod(POI_index(iPOI),5); if color_style==0; color_style=5; end
        line_style = ceil(POI_index(iPOI)/5); if line_style==5; line_style=1; end
        if iPOI>=21
            mycolors(6,:) = [0.3010 0.7450 0.9330];
            color_style = 6;
        end
        plot(time_pts_fine,Sti_vec_fine,'color', mycolors(color_style,:),'LineStyle',mylinestyles(line_style),'DisplayName',lP_list_name{POI_index(iPOI)})
    end
    legend('Orientation','horizontal','NumColumns',10)
    % legend('Location','eastoutside')
    if flag_save; saveas(gcf,[direc,'/eFAST_result_',num2str(NS),'_',num2str(k),'_',lQ{iQOI},'.eps'],'epsc'); end
    % if flag_save; saveas(gcf,[direc,'/eFAST_result_',num2str(NS),'_',num2str(k),'_',lQ{iQOI},'_all.eps'],'epsc'); end
end

