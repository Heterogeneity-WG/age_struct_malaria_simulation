clearvars
% close all;  
% clc; 
format long;
global P lP
flag_disp = 1;

%% numerical config
tfinal = 30*365; % final time in days
age_max = 100*365; % max ages in days
P.age_max = age_max;
dt = 20; % time/age step size in days, default = 5;
da = dt;
t = (0:dt:tfinal)';
nt = length(t);
a = (0:da:age_max)';
na = length(a);

P.a = a;
P.na = na;
P.nt = nt;
P.dt = dt;
P.da = da;
P.t = t;
P.tfinal = tfinal;

%% SA setting
% lQ = 'EE-D-frac';  % R0 RHM RMH EE-EIR EE-EDA EE-infected EE-D-frac
lQ = {'EE-D'};%,'EE-D-frac','EE-EIR',...
%     'EE-D-02-10','EE-DA-02-10','EE-D-frac-02-10',...
%     'EE-D-09-24','EE-DA-09-24','EE-D-frac-09-24',...
%     'EE-death','EE-death-02-10','EE-death-09-24'};  
lP_list = {'etas'};%{'bh', 'bm', 'betaM', 'betaD', 'betaA', 'muM', 'sigma'}; 
% 'bh', 'bm', 'betaM', 'betaD', 'betaA', 'muM', 'MHr', 'sigma'
% 'rho0','psi0','phi0'
SA_index = 0;
Malaria_parameters_baseline;
Malaria_parameters_transform;
Malaria_parameters_transform_vac;
%%
tic
if flag_disp; disp(['QOI =  ', lQ]); end
for iP = 1:length(lP_list)
    lP = lP_list{iP};   
    if flag_disp; disp(['POI =  ', lP]); end
    
    %% baseline run
    if flag_disp; disp('baseline run ---'); end
    Malaria_parameters_baseline;
    Malaria_parameters_transform_SA;

    P_baseline = P.(lP);
    Q_baseline = QOI_value_SA(lQ,1);
    
%     if SA_index
%         % SI index at baseline
%         Qp_baseline = QOIp_value(lQ,lP);
%         if length(Qp_baseline)==1
%             Qp_rescale_baseline(iP,:) = Qp_baseline'.*P_baseline'./Q_baseline;
%         else
%             Qp_rescale_baseline(iP,:) = norm(Qp_baseline'.*P_baseline'./Q_baseline,Inf);
%         end
%     end
    
    %% extended SA
    P_lower = P.([lP,'_lower']);
    P_upper = P.([lP,'_upper']);
    ngrid = 5; % default = 11
    
    % allocation
    P_vals = linspace(P_lower,P_upper,ngrid)';%
    Q_vals = NaN(length(P_vals),length(lQ));
    
    for i=1:ngrid
        display(['I am working on simulation ', num2str(i)])
        P.(lP) = P_vals(i);
        Malaria_parameters_transform_SA;
        Q_vals(i,:) = QOI_value_SA(lQ,1)';
    end
    
    %% plotting
    figure_setups; hold on
    plot(P_vals,Q_vals);
    plot(P_baseline,Q_baseline,'r*','MarkerSize',20);
    xlabel(lP)
    legend(lQ)    
%     figure_setups; hold on
%     area(P_vals,Q_vals)
%     plot([P_baseline, P_baseline],[0,1],'k-')
%     legend('E','D','A')
%     xlabel(lP)
%     ylabel(lQ)
end
toc
if SA_index
    [~,index] = sort(abs(Qp_rescale_baseline),'descend');
    Qp_rescale_baseline = Qp_rescale_baseline(index);
    lP_list = lP_list(index);
    disp('========')
    disp(lQ)
    SI_index = Qp_rescale_baseline;
    SS = num2str(SI_index);
    for ip = 1:length(lP_list)
        display([SS(ip,:), '    ', lP_list{ip}])
    end
end
%% generate output
Qlen = 1; Plen  = length(lP_list);
fname = sprintf('Results/SA/sensitivity_table_h.tex');
%(SI_index', 'Horiz', lP_list, 'Vert', {lQ},'Hline', [0:Qlen,NaN], 'Vline', [0:Plen,NaN],'name', fname, 'format', '%.2g');
% fname = sprintf('Results/SA/sensitivity_table_v.tex');
% latextable(SI_index, 'Horiz', {lQ}, 'Vert', lP_list,...
%     'Hline', [0:Plen,NaN], 'Vline', [0:Qlen,NaN],...
%     'name', fname, 'format', '%.2g');