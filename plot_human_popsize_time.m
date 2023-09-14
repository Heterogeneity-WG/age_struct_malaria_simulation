%% Population size versus time
colour_mat1 = [0 0.4470 0.7410];
colour_mat2 = [0.8500 0.3250 0.0980];
colour_mat3 = [0.9290 0.6940 0.1250];
colour_mat4 = [0.4940 0.1840 0.5560];
colour_mat5 = [0.4660 0.6740 0.1880];
colour_mat6 = [0.3010 0.7450 0.9330];
colour_mat7 = [0.6350 0.0780 0.1840];

figure_setups;
plot(t/365,trapz(SH,1)*da,'-','Color',colour_mat1); hold on;
plot(t/365,trapz(EH,1)*da,'--','Color',colour_mat3);
plot(t/365,trapz(AH,1)*da,'-.','Color',colour_mat2);
plot(t/365,trapz(DH,1)*da,'-','Color',colour_mat7);
plot(t/365,trapz(VH,1)*da,'-','Color',colour_mat6);
plot(t/365,trapz(UH,1)*da,'-','Color',colour_mat4);
plot(t/365,NH,'-.k')
plot(t/365,trapz(MH,1)*da,'-.r'); % diagnostic
legend('$S_H$','$E_H$','$A_H$', '$D_H$', '$V_H$','$U_H$','$N_H$','$M_H$ ($\mu_D$)', 'Location','e');
title(['Population size vs time']);
grid on; grid minor
axis([0 max(t)/365 0 max(NH)+0.1]);