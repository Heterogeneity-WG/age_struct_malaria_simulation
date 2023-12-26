clear all
% clc
%% Kenya malaria-related mortality data 
% downloaded from IHME
% @misc{GBD_data,
%   title = {{Global Burden of Disease Study 2019 (GBD 2019) Population Estimates 1950-2019}},
%   author = {{Institute for Health Metrics and Evaluation}},
%   howpublished = {\url{https://ghdx.healthdata.org/record/ihme-data/gbd-2019-population-estimates-1950-2019}},
%   note = {Accessed: 2022-03-29},
%   year = 2020
% }
% We use age-specific malaria-related deaths between ages x and x+n for both sexes in 2019
% Age Group: <1 year, 1-4 years, 5-9 years, 10-14 years,..., 70-74. We do
% not consider higher ages bc the curve trend changes
%nMDx: deaths from malaria/prevalence by age groups 0 - 74
nMDx = [0.03309 0.00437 0.00052 0.00048 0.00046 0.00057 0.00061 0.00085 0.00130...
    0.00187 0.00227 0.00381 0.00501 0.00826 0.00920 0.01357 0.00774857153383546...
    0.00482388846199508 0.00161635740909626 0.000238069451209632 1.77836885969402E-05];
muD0 = nMDx/365;
alpha0 = [0.5 2.5 7 12 17 22 27 32 37 42 47 52 57 62 67 72 77 82 87 92 97]*365;

%% estimate parameter values for the mortality rate function for ages 0 - 74: b0 + b1*exp(-b2*a) + b3*exp(b4*a)
modelfun = @(b,x) (b(1) + b(2)*exp(-b(3)*x./365) + b(4)*exp(b(5)*x./365))./365;

% initial guess for parameter values
b_0 = [0, 0.005, 0.0505, 0.001, 0.05];

options = statset('Display','final','TolFun',1e-12, 'MaxIter', 1e4, 'MaxFunEvals', 1e4);
md1 = fitnlm(alpha0(1:16), muD0(1:16), modelfun, b_0, 'Options', options);

% extract estimated parameter values
b_est = md1.Coefficients{:, 1};

% examine residuals
% figure(1)
% plotResiduals(md1, 'fitted')

% slice plot
% one can drag the vertical dashed line to see the effect of a change in alpha on mu
% plotSlice(md1)
%% estimate parameter values for the mortality rate function for ages 74+: b0 + b1*exp(-b2*74) + b3*exp(b4*74)
modelfun2 = @(c,x) modelfun(b_est,74*365).*(exp(-c*(x-74*365)./365));

% initial guess for parameter values
c_0 = 0.2;

options = statset('Display','final','TolFun',1e-12, 'MaxIter', 1e4, 'MaxFunEvals', 1e4);
md2 = fitnlm(alpha0(17:21), muD0(17:21), modelfun2, c_0, 'Options', options);

% extract estimated parameter values
c_est = md2.Coefficients{:, 1};

% examine residuals
% figure(1)
% plotResiduals(md1, 'fitted')

% slice plot
% one can drag the vertical dashed line to see the effect of a change in alpha on mu
% plotSlice(md1)
%% Check fitting for ages 0 - 74
figure_setups; hold on
scatter(alpha0(1:16)/365,muD0(1:16),'ko')
age = 0:20:74*365; % days
muD =  modelfun(b_est,age);
plot(age/365,muD)
muD =  modelfun(b_0,age);
plot(age/365,muD)
legend('data','estimate','initial')
title('Malaria mortality fit Kenya (0 - 74)')
xlabel('age')
ylabel('daily malaria-related death rate')
%% Check fitting for ages 74+
figure_setups; 
hold on
scatter(alpha0(17:21)/365,muD0(17:21),'ko')
age = 74*365:20:100*365; % days
muD =  modelfun2(c_est,age);
plot(age/365,muD)
muD =  modelfun2(c_0,age);
plot(age/365,muD)
legend('data','estimate','initial')
title('Malaria mortality fit Kenya (74+)')
xlabel('age')
ylabel('daily malaria-related death rate')