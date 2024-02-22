%% Figure A.1 (supplement for Figure 4)
% global SA (with vaccination-related parameters) using eFAST
% 'Data_SA/Results_local_SA_vac/SA_17POI_eFAST/'

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
direc = 'Data_SA/Results_local_SA_vac/SA_17POI_eFAST/';
flag_save = 1; % flag for saving the results or not (Note: it will overwrite previous results in the folder)

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

% code will calculate results for all the quantities below but only
% plotting a subset of this list. To modify the plotting, turn on
% "Size_QOI_plot" in the for loop of last section
lQ = {'EE-death','EE-death-09-24','EE-death-02-10','EE-death-10+',...
    'EE-DA','EE-DA-09-24','EE-DA-02-10','EE-DA-10+'};
Size_QOI = length(lQ); % length of the QOI.
time_points = length(t); % default # time_points = at tfinal, unless if wants to check QOI at particular time points
lP_list = {'rA','rD','cS','cA','cU','psis2','psir2','dac','uc','muM','betaM','betaD','betaA','v0','w','etas'};
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

for iQ = [1 2 5 6] % 1:Size_QOI_plot
    figure_setups; hold on
    model_series = [s_struct.Si(:,1,iQ)';s_struct.Sti(:,1,iQ)']';
    b = bar(X,model_series,'FaceColor','flat');
    b(1).FaceColor = 0.75*[1 1 1]; % light gray
    b(1).EdgeColor = 0.75*[1 1 1];
    b(2).FaceColor = [0.2 0.2 0.2]; % dark gray
    b(2).EdgeColor = [0.2 0.2 0.2];
    % mark p-values for Si
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YData;
    labels1 = cell(1, length(lP_list));
    labels1(s_struct.p_Si(:,1,1,iQ)<palpha) = {'*'};
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
    % mark p-values for Sti
    xtips2 = b(2).XEndPoints;
    ytips2 = b(2).YData;
    labels2 = cell(1, length(lP_list));
    labels2(s_struct.p_Sti(:,1,1,iQ)<palpha) = {'*'};
    text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom')
    title(['QOI = ',lQ_title{iQ}])
    xticklabels(lP_list_name)
    ylim([0 0.6])
    grid off
    legend('First Order','Total Order')
    temp = k+0.5;
    plot([0.5 temp],[0.1, 0.1],':','Linewidth',1.5,'Color', 0.4*[1 1 1])
    ax=gca;
    % read out the position of the axis in the unit "characters"
    set(ax,'Units','characters'); temp_ax=get(ax,'Position');
    % this sets an 'a)' right at the top left of the axes
    if iQ == 1
        text(ax,-12,temp_ax(end)+3,'(C)','Units','characters');
        save_string = strcat('figA1_','C','.svg');
    elseif iQ == 2
        text(ax,-12,temp_ax(end)+3,'(B)','Units','characters');
        save_string = strcat('figA1_','B','.svg');
    elseif iQ == 5
        text(ax,-12,temp_ax(end)+3,'(F)','Units','characters');
        save_string = strcat('figA1_','F','.svg');
    else
        text(ax,-12,temp_ax(end)+3,'(E)','Units','characters');
        save_string = strcat('figA1_','E','.svg');
    end
    saveas(gcf,save_string);
    if flag_save; saveas(gcf,[direc,'eFAST_result_',num2str(NS),'_',num2str(k),'_',lQ{iQ},'.eps'],'epsc'); end
end


