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

%% numerical config
tfinal = 100*365; % final time in days
age_max = 80*365; % max ages in days
P.age_max = age_max;
dt = 20; % time/age step size in days, default = 5;
da = dt;
t = (0:dt:tfinal)';
nt = length(t);
a = (0:da:age_max)';
na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t;

%% SA setting
lQ = 'EE-D';  % R0 RHM RMH EE-EIR EE-EDA EE-infected EE-D-frac
Size_QOI = 1; % length of the QOI. Default = 1, unless it is an age distribution, or wants to test multiple QOIs at once
time_points = 1; % default time_points = 1, unless if wants to check QOI at particular time points
lP_list = {'bh', 'bm', 'betaM', 'betaD', 'betaA', 'muM','sigma'};  % 'bh', 'bm', 'betaM', 'betaD', 'betaA', 'muM', 'MHr', 'sigma'
lP_list{end+1} = 'dummy'; % add dummy to the POIs
Malaria_parameters_baseline;
pmin = NaN(length(lP_list),1); pmax = pmin;
for iP = 1:length(lP_list)
    pmin(iP,1) = P.([lP_list{iP},'_lower']);
    pmax(iP,1) = P.([lP_list{iP},'_upper']);
end

%% eFAST config
NR = 3;    % # of search curves - Resampling
k = length(lP_list); % # of POIs + dummy parameter, keep it in the range 5~11
NS = 65;   % # of samples per search curve
wantedN=NS*k*NR; % wanted no. of sample points
MI = 4; %: maximum number of fourier coefficients that may be retained in calculating the partial variances without interferences between the assigned frequencies
% Computation of the frequency for the group of interest OMi and the # of sample points NS (here N=NS)
OMi = floor(((wantedN/NR)-1)/(2*MI)/k);
NS = 2*MI*OMi+1;
if(NS*NR < 65)
    fprintf(['Error: sample size must be >= 65 per factor.\n']);
    return;
end

%% Pre-allocation 
Size_timepts = length(time_points); % # of time points to check QOI value; 
Y(NS,Size_timepts,Size_QOI,length(pmin),NR)=0;  % For each model evaluation, the QOI has dimension [Size_timepts*Size_QOI]
X(NS,k,k,NR) = 0;
for i=1:k % Loop over POIs, including dummy
    OMci = SETFREQ(k,OMi/2/MI,i);  % selecting the set of frequencies. OMci(i), i=1:k-1, contains the set of frequencies to be used by the complementary group.
    for L=1:NR  % Loop over the NR search curves
        % Setting the vector of frequencies OM for the k parameters
        cj = 1;
        for j=1:k
            if(j==i)
                % For the parameter (factor) of interest
                OM(i) = OMi;
            else
                % For the complementary group.
                OM(j) = OMci(cj);
                cj = cj+1;
            end
        end
        % Setting the relation between the scalar variable S and the coordinates {X(1),X(2),...X(k)} of each sample point.
        FI = rand(1,k)*2*pi; % random phase shift
        S_VEC = pi*(2*(1:NS)-NS-1)/NS;
        OM_VEC = OM(1:k);
        FI_MAT = FI(ones(NS,1),1:k)';
        ANGLE = OM_VEC'*S_VEC+FI_MAT;
        
        X(:,:,i,L) = 0.5+asin(sin(ANGLE'))/pi; % between 0 and 1        
        % Transform distributions from standard uniform to general.
        X(:,:,i,L) = parameterdist(X(:,:,i,L),pmax,pmin,0,1,NS,'unif'); % this is what assigns 'our' values rather than 0:1 dist
        
        % model evaluations
        for run_num=1:NS
            disp([num2str([i, run_num, L]), '//parameter run NR']) % keeps track of [parameter run NR]
            % reset parameters to baseline
            Malaria_parameters_baseline;
            % update parameter with sampled values
            for iP = 1:length(lP_list)
                P.(lP_list{iP}) = X(run_num,iP,i,L);
            end
            % calculate QOI values
            Q_val = QOI_value_eFAST(lQ);
            if isnan(Q_val)
                % what to do to makeup the sample?
            end
            Y(run_num,:,:,i,L) = Q_val; 
        end 
    end 
end 

%% eFAST on output matrix Y
var = 1; % index of QOIs to analyze (among Size_QOI) (default = 1)
palpha = 0.05; % alpha for t-test
[Si,Sti,rangeSi,rangeSti] = efast_sd(Y,OMi,MI,time_points,var);
[CVsi, CVsti] = CVmethod(Si, rangeSi,Sti,rangeSti,var); % Coeff. of Variance; See online Supplement A.5 for details
s_HIV = efast_ttest(Si,rangeSi,Sti,rangeSti,1:length(time_points),lP_list,var,lQ,palpha); % T-test on Si and STi 
% save('eFAST_result.mat','s_HIV','lP_list','lQ','palpha')
keyboard

%% Plotting
clear all
close all
clc
load('eFAST_result.mat','s_HIV','lP_list','lQ','palpha')
figure_setups; hold on
X = categorical(lP_list);
X = reordercats(X,lP_list);
model_series = [s_HIV.Si';s_HIV.Sti']';
model_error = [std(s_HIV.rangeSi,0,3)';std(s_HIV.rangeSti,0,3)']';
b = bar(X,model_series);
% mark p-values for Si
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = cell(1, length(lP_list)); 
labels1(s_HIV.p_Si<palpha) = {'*'};
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
% mark p-values for Sti
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = cell(1, length(lP_list)); 
labels2(s_HIV.p_Sti<palpha) = {'*'};
text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom')
legend('first-order $S_i$','total-order $S_{Ti}$','Location','nw')
