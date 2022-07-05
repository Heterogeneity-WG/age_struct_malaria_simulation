%% TTDC test to check the convergence of the results in terms of increasing sample size NS
% first run the eFast for different NS, save the results s_struct_NS.mat
close all
clearvars
clc
% R[POIs, N]: matrix of ranks; (data replace by their ranks - from smallest to largest, tiedrank in MATLAB)
% eFAST: SI[POI,timepoints,QOI]; SI[POI,timepoints,QOI] 
NS_list = [65, 129, 257, 513];
k = 18; % k = 16 for vaccine_no; k = 18 for vaccine_yes
NR = 5;
nQ = 13; % # of QOIs
rT_si(length(NS_list)-1,nQ) = 0; rT_sti = rT_si;
p_si(length(NS_list)-1,nQ) = 0; p_sti = p_si;

for iNS = 1:length(NS_list)-1
    NS1 = NS_list(iNS);
    NS2 = NS_list(iNS+1);
    Data1 = load(['Results/vaccine_yes/eFAST_result_',num2str(NS1),'_',num2str(k),'.mat'],'s_struct');
    Data2 = load(['Results/vaccine_yes/eFAST_result_',num2str(NS2),'_',num2str(k),'.mat'],'s_struct');
    for iQ = 1:nQ 
        Si_1 = Data1.s_struct.Si(:,1,iQ); Sti_1 = Data1.s_struct.Sti(:,1,iQ);
        Si_2 = Data2.s_struct.Si(:,1,iQ); Sti_2 = Data2.s_struct.Sti(:,1,iQ);
        Rsi = tiedrank([Si_1, Si_2]); Rsti = tiedrank([Sti_1, Sti_2]); 
        [rT_si(iNS,iQ), p_si(iNS,iQ)] = topdowncorr(Rsi);
        [rT_sti(iNS,iQ), p_sti(iNS,iQ)] = topdowncorr(Rsti);
    end
end

T_Si = cell(length(NS_list)-1,nQ); T_Sti = T_Si;

for iQ = 1:nQ
    for iNS = 1:length(NS_list)-1
        if p_si(iNS,iQ)<0.01
            label_si = '**';
        elseif p_si(iNS,iQ)<0.05
            label_si = '*';
        else
            label_si = '';
        end
        T_Si{iNS,iQ} = [num2str(rT_si(iNS,iQ)),label_si];
        if p_sti(iNS,iQ)<0.01
            label_sti = '**';
        elseif p_sti(iNS,iQ)<0.05
            label_sti = '*';
        else
            label_sti = '';
        end
        T_Sti{iNS,iQ} = [num2str(rT_sti(iNS,iQ)),label_sti];
    end
end
% output table, like Table D.7 in the review paper
T_Si
T_Sti

% test example from the review paper
% Si_1 = [0.0173, 0.035, 0.0482, 0.0455, 0.087, 0.0383, 0.1145, 0.2676,0.0456];
% Si_2 = [0.0325, 0.0127, 0.0286, 0.0293, 0.0709, 0.0173,0.0702,0.3271,0.0237];
% Sti_1 = [0.3845, 0.5007, 0.5535, 0.5158, 0.6066, 0.46, 0.7098, 0.8538, 0.5102];
% Sti_2 = [0.5347, 0.3488, 0.4809, 0.4766, 0.5164, 0.4239,0.5321, 0.8339, 0.4667];

