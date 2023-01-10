clearvars
close all


%Load Parameters



%Parameters to Vary

param1_value=0.1;
param2_value=2;



%Run Model



%Output to Record

Output=ones(4,5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preparing output files
% Delete files for output data if they exist

% Change parameter values from decimals to strings --> 0.25=0p25
clearvars param1_value_print
if mod(param1_value,1)
    temp = num2str(param1_value);
    for i = 1:length(temp)
        if temp(i) == '.'
            param1_value_print(i) = 'p';
        else
            param1_value_print(i) = temp(i);
        end
    end
else
    param1_value_print = num2str(param1_value);
end

clearvars param2_value_print
if mod(param2_value,1)
    temp = num2str(param2_value);
    for i = 1:length(temp)
        if temp(i) == '.'
            param2_value_print(i) = 'p';
        else
            param2_value_print(i) = temp(i);
        end
    end
else
    param2_value_print = num2str(param2_value);
end


%Delete any previous copies
filename = ['Output_param1_',param1_value_print,'_param2_',param2_value_print,'.mat'];
if exist(filename, 'file')==2
  delete(filename);
end

%Write the output to the file - where Output is everthing we want in the
%file
save([filename],'Output');
%% To save a figure

% str = sprintf('MC_Ebola NumDeaths PeakInf TimePeak: R_0 = %.2f, I1_0=%i, I2_0=%i',R0, I1_0, I2_0);
% set(gcf,'Name',str,'NumberTitle','off')
% saveas(gcf,strcat(str,'.fig'))
