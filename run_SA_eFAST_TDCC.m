%% TTDC test to check the convergence of the results in terms of increasing sample size NS
% first run the eFast for different NS, save the results s_struct_NS.mat
close all
clear all
clc
% R[POIs, N]: matrix of ranks; (data replace by their ranks - from smallest to largest, tiedrank in MATLAB)
% SI[POI,timepoints,QOI] 
NS_list = [65, 129, 257];
k = 8;
NR = 5;
nQ = 2; % # of QOIs
rT_si(length(NS_list)-1,nQ) = 0; rT_sti = rT_si;
p_si(length(NS_list)-1,nQ) = 0; p_sti = p_si;

for iNS = 1:length(NS_list)-1
    NS1 = NS_list(iNS);
    NS2 = NS_list(iNS+1);
    Data1 = load(['Results/eFAST_result_',num2str(NS1),'_',num2str(k),'_',num2str(NR),'.mat'],'s_struct');
    Data2 = load(['Results/eFAST_result_',num2str(NS2),'_',num2str(k),'_',num2str(NR),'.mat'],'s_struct');
    for iQ = 1:nQ 
        Si_1 = Data1.s_struct.Si(:,iQ,4); Sti_1 = Data1.s_struct.Sti(:,iQ,4);
        Si_2 = Data2.s_struct.Si(:,iQ,4); Sti_2 = Data2.s_struct.Sti(:,iQ,4);
        Rsi = tiedrank([Si_1, Si_2]); Rsti = tiedrank([Sti_1, Sti_2]); 
        [rT_si(iNS,iQ), p_si(iNS,iQ)] = topdowncorr(Rsi);
        [rT_sti(iNS,iQ), p_sti(iNS,iQ)] = topdowncorr(Rsti);
    end
end

% output table, like Table D.7 in the review paper
[rT_si, double(p_si<0.01)+double(p_si<0.05), rT_sti, double(p_sti<0.01)+double(p_sti<0.05)] 

% test example from the review paper
% Si_1 = [0.0173, 0.035, 0.0482, 0.0455, 0.087, 0.0383, 0.1145, 0.2676,0.0456];
% Si_2 = [0.0325, 0.0127, 0.0286, 0.0293, 0.0709, 0.0173,0.0702,0.3271,0.0237];
% Sti_1 = [0.3845, 0.5007, 0.5535, 0.5158, 0.6066, 0.46, 0.7098, 0.8538, 0.5102];
% Sti_2 = [0.5347, 0.3488, 0.4809, 0.4766, 0.5164, 0.4239,0.5321, 0.8339, 0.4667];

