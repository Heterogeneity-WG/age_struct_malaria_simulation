%% global SA using eFAST, adapted from 
% Dr. Denise Kirschner - University of Michigan http://malthus.micro.med.umich.edu/lab/usadata/
% Output:
% SI[POI,timepoints,QOI] : first order sensitivity indices
% STI[POI,timepoints,QOI] : total effect sensitivity indices
% CV[POI,timepoints,QOI]

close all
clear all
clc
format long
global P

use_X = 0; % load pre-generated parameter samples, stored in matrix X "eFAST_samples.mat"

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
lQ = {'EE-D','EE-D-02-10','EE-D-09-24'};  % R0 RHM RMH EE-EIR EE-D-frac EE-D EE-DA EE-death
% 'EE-D','EE-D-frac','EE-DA','EE-D-02-10','EE-D-frac-02-10','EE-DA-02-10','EE-EIR','EE-death'
Size_QOI = length(lQ); % length of the QOI. Default = 1, unless it is an age distribution, or wants to test multiple QOIs at once
time_points = 1; % default time_points = 1, unless if wants to check QOI at particular time points
lP_list = {'dac','cX','w','etas','v0'};
% 'rA','muM','sigma','betaM','betaD', 'betaA','dac','cX','phis2','phir2','rhos2','rhor2','psis2','psir2','w','etas'
lP_list{end+1} = 'dummy'; % add dummy to the POIs
Malaria_parameters_baseline;
pmin = NaN(length(lP_list),1); pmax = pmin; pmean = pmin;
for iP = 1:length(lP_list)
    pmin(iP,1) = P.([lP_list{iP},'_lower']);
    pmax(iP,1) = P.([lP_list{iP},'_upper']);
    pmean(iP,1) = P.([lP_list{iP}]);
%     sprintf('%5s  %2f  %2f', lP_list{iP}, pmin(iP,1)./pmean(iP,1), pmax(iP,1)./pmean(iP,1))
end

%% eFAST config
NR = 5;    % # of search curves - Resampling - keep the value
k = length(lP_list); % # of POIs + dummy parameter, keep it in the range 5~11
NS = 513;   % # of samples per search curve - keep 2^n+1
wantedN=NS*k*NR; % wanted no. of sample points
MI = 4; %: maximum number of fourier coefficients that may be retained in calculating the partial variances without interferences between the assigned frequencies
% Computation of the frequency for the group of interest OMi and the # of sample points NS (here N=NS)
OMi = floor(((wantedN/NR)-1)/(2*MI)/k); % >=8
NS = 2*MI*OMi+1;
if(NS*NR < 65)
    fprintf(['Error: sample size must be >= 65 per factor.\n']);
    return;
end
%% Pre-allocation 
Size_timepts = length(time_points); % # of time points to check QOI value; 
Y(NS,Size_timepts,Size_QOI,length(pmin),NR)=0;  % For each model evaluation, the QOI has dimension [Size_timepts*Size_QOI]
X(NS,k,k,NR) = 0;

%% Generate parameter samples, stored in matrix X
if use_X
    load(['Results/eFAST_samples_',num2str(NS),'_',num2str(k),'_',num2str(NR),'.mat'],'X') 
else
    X = efast_gensamples(X,OMi,MI,pmin,pmax,pmean,'triangular'); % triangular distribution for POIs
%     X = efast_gensamples(X,OMi,MI,pmin,pmax,pmean,'unif'); % uniform distribution for POIs
%     save(['Results/eFAST_samples_',num2str(NS),'_',num2str(k),'_',num2str(NR),'.mat'],'X')
end
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

for i=1:k % Loop over POIs, including dummy
   for L=1:NR  % Loop over the NR search curves 
        for run_num=1:NS % Loop through each parameter sample
            disp([num2str([i, run_num, L]), '//parameter run NR']) % keeps track of [parameter run NR]
            Malaria_parameters_baseline; % reset parameters to baseline
            for iP = 1:length(lP_list) % update parameter with sampled values
                P.(lP_list{iP}) = X(run_num,iP,i,L);
            end
            Malaria_parameters_transform_SA; % update dependent parameters
            Q_val = QOI_value_SA(lQ,time_points); % calculate QOI values
%             if Q_val>1.1; keyboard;end % sanity check
            Y(run_num,:,:,i,L) = Q_val; 
        end 
   end
end
        
%% eFAST on output matrix Y
var = 1:length(lQ); % index of QOIs to analyze (among Size_QOI) (default = 1)
palpha = 0.05; % alpha for t-test
[Si,Sti,rangeSi,rangeSti] = efast_sd(Y,OMi,MI,time_points,var);
[CVsi, CVsti] = CVmethod(Si, rangeSi,Sti,rangeSti,var); % Coeff. of Variance; See online Supplement A.5 for details
s_struct = efast_ttest(Si,rangeSi,Sti,rangeSti,1:length(time_points),lP_list,var,lQ,palpha); % T-test on Si and STi 
% save(['Results/eFAST_result_',num2str(NS),'_',num2str(k),'_',num2str(NR),'.mat'],'s_struct','lP_list','lQ','palpha')
toc

%% Plotting
% load(['Results/eFAST_result_',num2str(NS),'_',num2str(k),'_',num2str(NR),'.mat'],'s_struct','lP_list','lQ','palpha')
X = categorical(lP_list);
X = reordercats(X,lP_list);
for iQ = 1:length(lQ)
    figure_setups; hold on
    model_series = [s_struct.Si(:,:,iQ)';s_struct.Sti(:,:,iQ)']';
    model_error = [std(s_struct.rangeSi(:,:,:,iQ),0,3)';std(s_struct.rangeSti(:,:,:,iQ),0,3)']';
    b = bar(X,model_series);
    % mark p-values for Si
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels1 = cell(1, length(lP_list));
    labels1(s_struct.p_Si(:,:,:,iQ)<palpha) = {'*'};
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
    % mark p-values for Sti
    xtips2 = b(2).XEndPoints;
    ytips2 = b(2).YEndPoints;
    labels2 = cell(1, length(lP_list));
    labels2(s_struct.p_Sti(:,:,:,iQ)<palpha) = {'*'};
    text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom')
    title(['QOI=',lQ{iQ}])
%     saveas(gcf,['Results/',lQ{iQ},'_',num2str(NS),'.png'])
end