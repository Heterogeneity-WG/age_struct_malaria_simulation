clear all
% clc
%% 2019 Deaths KENYA all causes minus Malaria
% downloaded from GHO life tables
% nMx - age-specific death rate between ages x and x+n for both sexes in 2019
% Age Group: 1-4 years, 5-9 years, 10-14 years,..., 80-84 years and 85+ years 
%nMHx: Deaths all causes minus malaria / pop by age group
nMHx = [0.02958 0.00202 0.00060 0.00058 0.00129 0.00170 0.00223 0.00321 0.00479 0.00687 0.00935 0.01222...
0.01581 0.02197 0.03129 0.04764 0.07258 0.11234 0.16533 0.22615 0.32164];
%%
muH0 = nMHx/365;%-log(1-nqx)/365;
alpha0 = [0.5 2.5 7 12 17 22 27 32 37 42 47 52 57 62 67 72 77 82 87 92 97]*365;

%% estimate parameter values for the mortality rate function b0 + b1*exp(-b2*a) + b3*exp(b4*a)
modelfun = @(b,x) (b(1) + b(2)*exp(-b(3)*x./365) + b(4)*exp(b(5)*x./365))./365;

% initial guess for parameter values
b_0 = [0, 0.005, 0.0505, 0.001, 0.05];

options = statset('Display','final','TolFun',1e-12, 'MaxIter', 1e4, 'MaxFunEvals', 1e4);
md1 = fitnlm(alpha0, muH0, modelfun, b_0, 'Options', options)

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
scatter(alpha0/365,muH0,'ko')
age = 0:20:100*365; % days
muH =  modelfun(b_est,age);
plot(age/365,muH)
muH =  modelfun(b_0,age);
plot(age/365,muH)
legend('data','estimate','initial')
title('Non Malaria mortality fit Kenya')
xlabel('age')
ylabel('daily nonmalaria-related death rate')