%% global SA using eFAST, adapted from 
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

% numerical config
tfinal = 10*365; % time for "EE" approximation
tfinal_vac = 0*365; % time to achieve "new EE" w/ vaccination
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
% lQ = {'EE-D','EE-DA','EE-D-frac','EE-EIR',...
%     'EE-D-02-10','EE-DA-02-10','EE-D-frac-02-10',...
%     'EE-D-09-24','EE-DA-09-24','EE-D-frac-09-24',...
%     'EE-death','EE-death-02-10','EE-death-09-24',...
%     'DALY'};  
lQ = {'EE-death','EE-death-02-10','EE-death-09-24','EE-death-10+',...
      'EE-DA','EE-DA-02-10','EE-DA-09-24','EE-DA-10+'};  
Size_QOI = length(lQ); % length of the QOI. Default = 1, unless it is an age distribution, or wants to test multiple QOIs at once
time_points = length(t); % default time_points = 1, unless if wants to check QOI at particular time points
% time_points = 1:nt; 
lP_list = {'cS','cE','cA','cD','cU','phis2','phir2','rhos2','rhor2','psis2','psir2','dac','uc','m',...
    'rA','rD','muM','sigma','betaM','betaD', 'betaA'};
lP_list{end+1} = 'dummy'; % add dummy to the POIs
Malaria_parameters_baseline;
pmin = NaN(length(lP_list),1); pmax = pmin; pmean = pmin;
for iP = 1:length(lP_list)
    pmin(iP,1) = P.([lP_list{iP},'_lower']);
    pmax(iP,1) = P.([lP_list{iP},'_upper']);
    pmean(iP,1) = P.([lP_list{iP}]);
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
direc = 'D:/Results_local_SA/SA_22POI_eFAST/';
if ~exist([direc,'eFAST_samples_',num2str(NS),'_',num2str(k),'_',num2str(NR),'.mat'],'file')
    disp('generate parameter samples...')
    X = efast_gensamples(X,OMi,MI,pmin,pmax,pmean,'triangular'); % triangular distribution for POIs
    % X = efast_gensamples(X,OMi,MI,pmin,pmax,pmean,'unif'); % uniform distribution for POIs
    save([direc,'/eFAST_samples_',num2str(NS),'_',num2str(k),'_',num2str(NR),'.mat'],'X')
else
    disp('load parameter samples...')
    load([direc,'eFAST_samples_',num2str(NS),'_',num2str(k),'_',num2str(NR),'.mat'],'X')
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
P.gM0 = P.gM;
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

for i=1:k % Loop over POIs, including dummy
   for L=1:NR  % Loop over the NR search curves 
        for run_num=1:NS % Loop through each parameter sample
            ind = NS*NR*(i-1)+NS*(L-1)+run_num;
            disp([num2str(ind/(NS*NR*k)*100,3),' %']); % display progress
            % disp([num2str([i, run_num, L]), '//parameter run NR']) % keeps track of [parameter run NR]
            Malaria_parameters_baseline; % reset parameters to baseline
            for iP = 1:length(lP_list) % update parameter with sampled values
                P.(lP_list{iP}) = X(run_num,iP,i,L);
            end
            Malaria_parameters_transform_SA; % update dependent parameters
            Q_val = QOI_value_SA(lQ,time_points,ind,'eFAST'); % calculate QOI values
            Y(run_num,:,:,i,L) = Q_val; 
        end 
   end
end
        
%% eFAST on output matrix Y
var = 1:length(lQ); % index of QOIs to analyze (among Size_QOI) (default = 1)
palpha = 0.05; % alpha for t-test
[Si,Sti,rangeSi,rangeSti] = efast_sd(Y,OMi,MI,time_points,var);
[CVsi, CVsti] = CVmethod(Si, rangeSi,Sti,rangeSti,var); % Coeff. of Variance; See online Supplement A.5 for details
s_struct = efast_ttest(Si,rangeSi,Sti,rangeSti,time_points,lP_list,var,lQ,palpha); % T-test on Si and STi 
save([direc,'eFAST_result_',num2str(NS),'_',num2str(k),'.mat'],'s_struct','lP_list','lQ','palpha')
toc

%% Sorting
% load([direc,'eFAST_result_',num2str(NS),'_',num2str(k),'.mat'],'s_struct','lP_list','lQ','palpha')
lP_eFAST = {'rD','dac','uc','cS','psir2','betaM','muM','cA','rhos2','psis2','rA','cE',...
    'betaD','m','sigma','rhor2','betaA','cD','phir2','cU','phis2'}; 
[~,index] = ismember(lP_eFAST,lP_list); index = index';
% [~,index] = sort(abs(s_struct.Si(1:end-1,1,1)),'descend');
s_struct.Si = s_struct.Si([index;k],:,:); s_struct.p_Si = s_struct.p_Si(index,:,:,:); s_struct.rangeSi = s_struct.rangeSi([index;k],:,:,:);
s_struct.Sti = s_struct.Sti([index;k],:,:); s_struct.p_Sti = s_struct.p_Sti(index,:,:,:); s_struct.rangeSti = s_struct.rangeSti([index;k],:,:,:);
lP_list = lP_list([index;k]);

% Output to table
% fname = sprintf('Results/vaccine_no/SA_eFAST_Si.tex');
% latextable(squeeze(s_struct.Si)', 'Horiz', lP_list, 'Vert', lQ,...
%     'Hline', [0,1,NaN], 'Vline', [0,1,NaN],...
%     'name', fname, 'format', '%3.2f');

%% eFAST plot 
X = categorical(lP_list);
X = reordercats(X,lP_list);
QOI_plot = 1:length(lQ);
Size_QOI_plot = length(QOI_plot);
% lP_list_name = lP_list;
[lP_list_name,lQ] = SA_output_formatting(lP_list,lQ,1);
for iQ = 1:Size_QOI_plot
    figure_setups; hold on
    model_series = [s_struct.Si(:,1,iQ)';s_struct.Sti(:,1,iQ)']';
%     model_error = [std(s_struct.rangeSi(:,1,:,iQ),0,3)';std(s_struct.rangeSti(:,1,:,iQ),0,3)']';
    b = bar(X,model_series);
    % mark p-values for Si
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels1 = cell(1, length(lP_list));
    labels1(s_struct.p_Si(:,1,1,iQ)<palpha) = {'*'};
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
    % mark p-values for Sti
    xtips2 = b(2).XEndPoints;
    ytips2 = b(2).YEndPoints;
    labels2 = cell(1, length(lP_list));
    labels2(s_struct.p_Sti(:,1,1,iQ)<palpha) = {'*'};
    text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom')
    title(['QOI=',lQ{iQ}])
    xticklabels(lP_list_name)
    ylim([0 0.5])
    saveas(gcf,[direc,'/eFAST_result_',num2str(NS),'_',num2str(k),'_',lQ{iQ},'.eps'],'epsc')
end