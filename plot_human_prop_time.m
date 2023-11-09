[~,ind1] = min(abs(P.a-5*30));
[~,ind2] = min(abs(P.a-17*30));
%% proportion for the entire population
figure_setups_2;
colour_mat1 = [0 0.4470 0.7410];
colour_mat2 = [0.8500 0.3250 0.0980];
colour_mat3 = [0.9290 0.6940 0.1250];
colour_mat4 = [0.4940 0.1840 0.5560];
colour_mat5 = [0.4660 0.6740 0.1880];
colour_mat6 = [0.3010 0.7450 0.9330];
colour_mat7 = [0.6350 0.0780 0.1840];
plot(t/365,trapz(SH,1)*da./NH,'-','Color',colour_mat1); hold on;
plot(t/365,trapz(EH,1)*da./NH,'-','Color',colour_mat3);
plot(t/365,trapz(AH,1)*da./NH,'-','Color',colour_mat2);
plot(t/365,trapz(DH,1)*da./NH,'-','Color',colour_mat7);
plot(t/365,trapz(VH,1)*da./NH,'-','Color',colour_mat6);
plot(t/365,trapz(UH,1)*da./NH,'-','Color',colour_mat4);
% plot(t/365,(trapz(SH,1)+trapz(EH,1)+trapz(AH,1)+trapz(DH,1)+trapz(VH,1)+trapz(UH,1))*da./NH,'-.k');
% legend('$\tilde{S}_H$','$\tilde{E}_H$','$\tilde{A}_H$', '$\tilde{D}_H$', '$\tilde{V}_H$','$\tilde{U}_H$','$N_H$');
legend('$\tilde{S}_H$','$\tilde{E}_H$','$\tilde{A}_H$', '$\tilde{D}_H$', '$\tilde{V}_H$','$\tilde{U}_H$');
title('Disease dynamics');
xlabel('years');
ylabel('proportion')
grid on
% axis([0 max(t)/365 0 1.1]);
axis([0 max(t)/365 0 0.6]);
%% plotting different age groups 
% figure_setups;
% 
% tempNorm = SH./PH;
% subplot(2,2,1), plot(t/365,tempNorm(end,:));
% hold on;
% plot(t/365,tempNorm(floor(end/10),:));
% plot(t/365,tempNorm(floor(end/100),:));
% plot(t/365,trapz(SH(ind1:ind2,:),1)./trapz(PH(ind1:ind2,:),1));
% ylim([0 1]);
% title('$\tilde{S}_{H}(\alpha,t)$');
% grid on;
% 
% tempNorm = EH./PH;
% subplot(2,2,2), plot(t/365,tempNorm(end,:));
% hold on;
% plot(t/365,tempNorm(floor(end/10),:));
% plot(t/365,tempNorm(floor(end/100),:));
% plot(t/365,trapz(EH(ind1:ind2,:),1)./trapz(PH(ind1:ind2,:),1));
% title('$\tilde{E}_{H}(\alpha,t)$');
% legend('Age 100', 'Age 10','Age 1','Age target');
% ylim([0 1]);
% grid on;
% 
% tempNorm = AH./PH;
% subplot(2,2,3), plot(t/365,tempNorm(end,:));
% hold on;
% plot(t/365,tempNorm(floor(end/10),:));
% plot(t/365,tempNorm(floor(end/100),:));
% plot(t/365,trapz(AH(ind1:ind2,:),1)./trapz(PH(ind1:ind2,:),1));
% title('$\tilde{A}_{H}(\alpha,t)$');
% ylim([0 1]);
% grid on;
% xlabel('Time (years)');
% 
% tempNorm = DH./PH;
% subplot(2,2,4), plot(t/365,tempNorm(end,:));
% hold on;
% plot(t/365,tempNorm(floor(end/10),:));
% plot(t/365,tempNorm(floor(end/100),:));
% plot(t/365,trapz(DH(ind1:ind2,:),1)./trapz(PH(ind1:ind2,:),1));
% title('$\tilde{D}_{H}(\alpha,t)$');
% ylim([0 1]);
% grid on;
% xlabel('Time (years)');