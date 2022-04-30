%% Plot Sigmoidal Probability Function
clear all
close all
clc
x0 = [4   5   4   5];
f_0 = 0.01; % value at zero
f_1 = 1; % value at L (function saturates to this value)
x = 0:0.01:13; %input
s_2 = x0(1);
r_2 = x0(2);
y_inc1 = f_0 + (f_1-f_0)./(1 + exp(-(x-s_2)./r_2));
s_2 = x0(3);
r_2 = x0(4);
y_dec1 = f_1 + (f_0-f_1)./(1 + exp(-(x-s_2)./r_2)); % switch f0 f1
% s_2 = s_2*1.5;
% r_2 = r_2*1.5;
% y_inc2 = f_0 + (f_1-f_0)./(1 + exp(-(x-s_2)./r_2));
% y_dec2 = f_1 + (f_0-f_1)./(1 + exp(-(x-s_2)./r_2)); % switch f0 f1
% % 
% s_2 = s_2*0.5;
% r_2 = r_2*0.5;
% y_inc3 = f_0 + (f_1-f_0)./(1 + exp(-(x-s_2)./r_2));
% y_dec3 = f_1 + (f_0-f_1)./(1 + exp(-(x-s_2)./r_2)); % switch f0 f1

figure_setups; hold on
plot(x,y_inc1,'linewidth',2)
plot(x,y_dec1,'linewidth',2)
% plot(x,y_inc2,'linewidth',2)
% plot(x,y_dec2,'linewidth',2)
% plot(x,y_inc3,'linewidth',2)
% plot(x,y_dec3,'linewidth',2)
% plot(x,f_0*ones(size(x)),'k--')
% plot(x,f_1*ones(size(x)),'k--')
% text(0.85*x(end),f_0+0.05,'fmin')
% text(0.85*x(end),f_1+0.05,'fmax')
axis([min(x) max(x) 0, 1])
xlabel('$C_H$, population immunity level')
ylabel('$\rho$, probability')

