function s = efast_ttest(Si,rangeSi, Sti,rangeSti,time_points,efast_var,y_var,y_var_label,alpha) % use rangeSi or ramgeSti from efast_sd
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
g=1;
for u=y_var
    for s=1:length(time_points)
        z=['time = ',num2str(s)];
        [z y_var_label(u)]
        %% Compare Si or STi of parameter j with the dummy
        for i=1:k-1
            %% Si
            [rangeSi(i,s,:,u) rangeSi(k,s,:,u)];
            [squeeze(rangeSi(i,s,:,u)) squeeze(rangeSi(k,s,:,u))];
            [mean(squeeze(rangeSi(i,s,:,u))) mean(squeeze(rangeSi(k,s,:,u)))];
            [a,p_Si(i,s,u)]=ttest2(squeeze(rangeSi(i,s,:,u)),squeeze(rangeSi(k,s,:,u)),alpha,'right','unequal');
            %% Sti
            [rangeSti(i,s,:,u) rangeSti(k,s,:,u)];
            [squeeze(rangeSti(i,s,:,u)) squeeze(rangeSti(k,s,:,u))];
            [mean(squeeze(rangeSti(i,s,:,u))) mean(squeeze(rangeSti(k,s,:,u)))];
            [a,p_Sti(i,s,u)]=ttest2(squeeze(rangeSti(i,s,:,u)),squeeze(rangeSti(k,s,:,u)),alpha,'right','unequal');
        end % for i
        SAefast_struct.p_Si(:,:,s,u)=p_Si(:,s,u);
        SAefast_struct.p_Sti(:,:,s,u)=p_Sti(:,s,u);
    end % for t
end
s=SAefast_struct;
% output results
efast_var  % POIs
Si_out=squeeze(s.Si(:,:,y_var))'
p_Si_out=squeeze(s.p_Si(:,:,:,y_var))'
Sti_out=squeeze(s.Sti(:,:,y_var))'
p_Sti_out=squeeze(s.p_Si(:,:,:,y_var))'
end