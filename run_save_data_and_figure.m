clearvars
close all


%Load Parameters

format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

tic

%% numerical config
age_max = 100*365; % max ages in days
P.age_max = age_max;
dt = 10; % time/age step size in days, default = 5;
da = dt;
a = (0:da:age_max)';
na = length(a);

P.dt = dt; 
P.a = a;
P.na = na;
P.da = da;

% model parameters
Malaria_parameters_baseline;
Malaria_parameters_transform; 
Malaria_parameters_transform_vac;

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
