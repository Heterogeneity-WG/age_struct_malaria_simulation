%% Plot Sigmoidal Probability Function
clear all
% close all
% clc
f_0 = 0.01; % value at zero
f_1 = 1; % value at L (function saturates to this value)
t_2 = 0.5 ; % threshold value (as a fraction of L)
s_2 = 0.1; % sigmoid steepness, smaller is steeper
x = 0:0.01:8; %input
%            
y_inc = f_0 + (f_1-f_0)./(1 + exp(-(x./L-t_2)./s_2));
y_dec = f_1 + (f_0-f_1)./(1 + exp(-(x./L-t_2)./s_2)); % switch f0 f1

figure_setups; hold on
plot(x,y_inc,'linewidth',2)
plot(x,y_dec,'linewidth',2)
plot(x,f_0*ones(size(x)),'k--')
plot(x,f_1*ones(size(x)),'k--')
text(0.85*x(end),f_0+0.05,'fmin')
text(0.85*x(end),f_1+0.05,'fmax')
axis([min(x) max(x) 0, 1])
xlabel('$C_H$, population immunity level')
ylabel('$\rho$, probability')

