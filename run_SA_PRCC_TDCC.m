%% TTDC test to check the convergence of the results in terms of increasing sample size NS
% first run the PRCC for different NS, save the results 
% close all
clearvars
clc
% R[POIs, N]: matrix of ranks; (data replace by their ranks - from smallest to largest, tiedrank in MATLAB)
% PRCC[POI,timepoints,QOI]
NS_list = [100,200,500,1000];
k = 18; % k = 16 for vaccine_no; k = 18 for vaccine_yes
nQ = 13; % # of QOIs
rT_PRCC(length(NS_list)-1,nQ) = 0; 
p_PRCC(length(NS_list)-1,nQ) = 0;

for iNS = 1:length(NS_list)-1
    NS1 = NS_list(iNS);
    NS2 = NS_list(iNS+1);
    Data1 = load(['Results/vaccine_yes/PRCC_result_',num2str(NS1),'_',num2str(k),'.mat'],'PRCC');
    Data2 = load(['Results/vaccine_yes/PRCC_result_',num2str(NS2),'_',num2str(k),'.mat'],'PRCC');
    for iQ = 1:nQ 
        PRCC_1 = Data1.PRCC(:,1,iQ);
        PRCC_2 = Data2.PRCC(:,1,iQ); 
        R_PRCC = tiedrank([PRCC_1, PRCC_2]); 
        [rT_PRCC(iNS,iQ), p_PRCC(iNS,iQ)] = topdowncorr(R_PRCC);
    end
end

T_PRCC = cell(length(NS_list)-1,nQ); 

for iQ = 1:nQ
    for iNS = 1:length(NS_list)-1
        if p_PRCC(iNS,iQ)<0.01
            label_PRCC = '**';
        elseif p_PRCC(iNS,iQ)<0.05
            label_PRCC = '*';
        else
            label_PRCC = '';
        end
        T_PRCC{iNS,iQ} = [num2str(rT_PRCC(iNS,iQ)),label_PRCC];
    end
end
% output table, like Table D.7 in the review paper
T_PRCC

% test example from the review paper
% Si_1 = [0.0173, 0.035, 0.0482, 0.0455, 0.087, 0.0383, 0.1145, 0.2676,0.0456];
% Si_2 = [0.0325, 0.0127, 0.0286, 0.0293, 0.0709, 0.0173,0.0702,0.3271,0.0237];
% Sti_1 = [0.3845, 0.5007, 0.5535, 0.5158, 0.6066, 0.46, 0.7098, 0.8538, 0.5102];
% Sti_2 = [0.5347, 0.3488, 0.4809, 0.4766, 0.5164, 0.4239,0.5321, 0.8339, 0.4667];

