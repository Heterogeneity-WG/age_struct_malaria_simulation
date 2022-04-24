clear all
% clc
%% Kenya malaria-related mortality data 
% downloaded from IHME
% We use age-specific malaria-related deaths between ages x and x+n for both sexes in 2019
% Age Group: <1 year, 1-4 years, 5-9 years, 10-14 years,..., 70-74. We do
% not consider higher ages bc the curve trend changes
%nMDx: deaths from malaria/prevalence by age groups
nMDx = [0.03309 0.00437 0.00052 0.00048 0.00046 0.00057 0.00061 0.00085 0.00130 0.00187 0.00227 0.00381 0.00501 0.00826 0.00920 0.01357];
%%
muD0 = nMDx/365;%-log(1-nqx)/365;
%Q: why are we multiplying alpha by 365?
alpha0 = [0.5 2.5 7 12 17 22 27 32 37 42 47 52 57 62 67 72]*365;

%% estimate parameter values for the mortality rate function b0 + b1*exp(-b2*a) + b3*exp(b4*a)
modelfun = @(b,x) (b(1) + b(2)*exp(-b(3)*x./365) + b(4)*exp(b(5)*x./365))./365;

% initial guess for parameter values
b_0 = [0, 0.005, 0.0505, 0.001, 0.05];

options = statset('Display','final','TolFun',1e-12, 'MaxIter', 1e4, 'MaxFunEvals', 1e4);
md1 = fitnlm(alpha0, muD0, modelfun, b_0, 'Options', options)

% extract estimated parameter values
b_est = md1.Coefficients{:, 1};

% examine residuals
% figure(1)
% plotResiduals(md1, 'fitted')

% slice plot
% one can drag the vertical dashed line to see the effect of a change in alpha on mu
% plotSlice(md1)
%% 
figure_setups; hold on
scatter(alpha0/365,muD0,'ko')
age = 0:20:100*365; % days
muD =  modelfun(b_est,age);
plot(age/365,muD)
muD =  modelfun(b_0,age);
plot(age/365,muD)
legend('data','estimate','initial')
title('Malaria mortality fit Kenya')
xlabel('age')
ylabel('daily malaria-related death rate')