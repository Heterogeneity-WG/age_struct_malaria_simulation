function S = efast_ttest(Si,rangeSi,Sti,rangeSti,time_points,efast_var,y_var,y_var_label,alpha) % use rangeSi or ramgeSti from efast_sd
%Si=Si(:,time_points);
%Sti=Sti(:,time_points);
%rangeSi=rangeSi(:,time_points,:);
%rangeSti=rangeSti(:,time_points,:);
%Si_struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
function pairwise_ttest

input: matrix of sensitivity indices Si(# parameters,# search curves)
output: square matrix of p-values, comparing each parameter
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SAefast_struct=struct;

[k,t,NR,output]=size(rangeSi); %record k number of parameters, NR number of search curves

SAefast_struct.Si=Si;
SAefast_struct.rangeSi=rangeSi;
SAefast_struct.Sti=Sti;
SAefast_struct.rangeSti=rangeSti;

for iQOI=1:length(y_var)
    for itime=1:length(time_points)
        ['time = ',num2str(itime),'  ', y_var_label(iQOI)]
        %% Compare Si or STi of parameter j with the dummy
        for iPOI=1:k-1
            % Si
            [a,p_Si(iPOI,itime,iQOI)]=ttest2(squeeze(rangeSi(iPOI,itime,:,iQOI)),squeeze(rangeSi(k,itime,:,iQOI)),alpha,'right','unequal');
            % Sti
            [a,p_Sti(iPOI,itime,iQOI)]=ttest2(squeeze(rangeSti(iPOI,itime,:,iQOI)),squeeze(rangeSti(k,itime,:,iQOI)),alpha,'right','unequal');
        end % for i
        SAefast_struct.p_Si(:,:,itime,iQOI)=p_Si(:,itime,iQOI);
        SAefast_struct.p_Sti(:,:,itime,iQOI)=p_Sti(:,itime,iQOI);
    end % for t
end
S=SAefast_struct;
% efast_var  % POIs
% Si_out=squeeze(S.Si(:,:,y_var))
% p_Si_out=squeeze(S.p_Si(:,:,:,y_var))
% Sti_out=squeeze(S.Sti(:,:,y_var))
% p_Sti_out=squeeze(S.p_Si(:,:,:,y_var))
end