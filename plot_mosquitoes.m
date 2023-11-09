colour_r1 = [1, 0.6, 0.6]; % 80% red [255, 153, 153]
colour_r2 = [1, 0.4, 0.4]; % 70% red [255, 102, 102]
colour_r3 = [1, 0, 0]; % 50% red  red [255, 0, 0]
colour_r4 = [0.8, 0, 0]; % 40% red  red [204, 0, 0]
%% mosquito/human population ratio
figure_setups;
plot(t/365,NM./NH) 
legend('$N_M/N_H$ ratio')

%% Mosquito population size
figure_setups;
plot(t/365,SM,'b-'); hold on;
plot(t/365,EM,'-','Color',colour_r1);
plot(t/365,IM,'r-.');
plot(t/365,SM+EM+IM,'-.')
plot(t/365, P.gM_fun(t)./P.muM,'*')
legend('$S_M$','$E_M$','$I_M$','$N_M$');
title('mosquito population size by stages')
grid on
xlim([0 tfinal/365])
%% Mosquito population proportion
figure_setups;
plot(t/365,SM./NM,'b-'); hold on;
plot(t/365,EM./NM,'-','Color',colour_r1);
plot(t/365,IM./NM,'r-.');
plot(t/365,(SM+EM+IM)./NM,'-.')
legend('$S_M$','$E_M$','$I_M$','$N_M$');
title('mosquito population prop by stages')
grid on
xlim([0 tfinal/365])