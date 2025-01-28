%% Age profiles at tfinal
% figure_setups;
colour_mat1 = [0 0.4470 0.7410];
colour_mat2 = [0.8500 0.3250 0.0980];
colour_mat3 = [0.9290 0.6940 0.1250];
colour_mat4 = [0.4940 0.1840 0.5560];
colour_mat5 = [0.4660 0.6740 0.1880];
colour_mat6 = [0.3010 0.7450 0.9330];
colour_mat7 = [0.6350 0.0780 0.1840];
% plot(a/365,SH(:,end),'-','Color',colour_mat1); hold on;
% plot(a/365,EH(:,end),'--','Color',colour_mat3);
% plot(a/365,DH(:,end),'-.','Color',colour_mat2);
% plot(a/365,AH(:,end),':','Color',colour_mat7);
% plot(a/365,VH(:,end),':','Color',colour_mat6);
% plot(a/365,UH(:,end),':','Color',colour_mat4);
% plot(a/365,PH_final,'-k');
% % plot(a/365,MH(:,end),'-r');
% hold off;
% title(num2str(t(end)/365));
% xlim([0 15])
% legend('$S_H$','$E_H$','$D_H$','$A_H$','$V_H$','$U_H$','$P_H$');
% title(['Final Age Distribution']);
% xlabel('age (years)');
% grid on
% axis([0 age_max/365 0 max(PH_final)]);
%% Age proportions at tfinal
% figure_setups;
% plot(a/365,SH(:,end)./PH_final,'-','Color',colour_mat1); hold on;
% plot(a/365,EH(:,end)./PH_final,'--','Color',colour_mat6);
% plot(a/365,DH(:,end)./PH_final,'-.','Color',colour_mat2);
% plot(a/365,AH(:,end)./PH_final,':','Color',colour_mat3);
% % plot(a/365,VH(:,end)./PH_final,':','Color',colour_mat5);
% % plot(a/365,UH(:,end)./PH_final,':','Color',colour_mat4);
% % plot(a/365,MH(:,end)./PH_final,'r-');
% % plot(a/365,(AH(:,end)+DH(:,end))./PH_final,'r-.');
% % plot(a/365,PH_final./PH_final,'-k');
% % legend('SH/PH','EH/PH','DH/PH', 'AH/PH', 'VH/PH','UH/PH','$M_H$ ($\mu_D$)','(AH+DH)/PH');
% legend('SH/PH','EH/PH','DH/PH', 'AH/PH');
% title(['Final Age Dist. Proportions']);
% % title(['Final Age Dist. Proportions ~~ feedback =',num2str(immunity_feedback)]);
% xlabel('age (years)');
% grid on
% axis([0 P.age_max/365 0 1.1]);
% xlim([0 30])

%% Age proportions in time - movie
figure_setups;
for iplot = 1:12:length(t)
    plot(a/365,SH(:,iplot),'-','Color',colour_mat1); hold on;
    plot(a/365,EH(:,iplot),'--','Color',colour_mat6);
    plot(a/365,DH(:,iplot),'-.','Color',colour_mat2);
    plot(a/365,AH(:,iplot),':','Color',colour_mat3);
    plot(a/365,VH(:,iplot),':','Color',colour_mat5);
    plot(a/365,UH(:,iplot),':','Color',colour_mat4);
    plot(a/365,(AH(:,iplot)+DH(:,iplot)),'r-.');
    plot(a/365,PH(:,iplot),'-k');
    % plot(a/365,MH(:,end)./PH_final,'r-');
    hold off
    legend('SH/P H','EH/PH','DH/PH', 'AH/PH', 'VH/PH','UH/PH','(AH+DH)/PH');
    title(['time = ', num2str(t(iplot)/30), ' months']);
    xlabel('age (years)');
    grid on
    % axis([0 P.age_max/365 0 1.1]);
    xlim([0 5])
%     saveas(gcf,['Results/seasonality_',num2str(t(iplot)/30),'.png'])
    pause
end  
%% Age proportions in time - movie
figure_setups;
for iplot = 1:12:length(t)
    subplot(1,2,1)
    plot(a/365,SHv(:,iplot),'-','Color',colour_mat1); hold on;
    plot(a/365,EHv(:,iplot),'--','Color',colour_mat6);
    plot(a/365,DHv(:,iplot),'-.','Color',colour_mat2);
    plot(a/365,AHv(:,iplot),':','Color',colour_mat3);
    plot(a/365,VHv(:,iplot),':','Color',colour_mat5);
    plot(a/365,UHv(:,iplot),':','Color',colour_mat4);
    plot(a/365,(AHv(:,iplot)+DHv(:,iplot)),'r-.');
    plot(a/365,PH(:,iplot),'-k');
    % plot(a/365,MH(:,end)./PH_final,'r-');
    hold off
    legend('SH/PH','EH/PH','DH/PH', 'AH/PH', 'VH/PH','UH/PH','(AH+DH)/PH');
    title(['time = ', num2str(t(iplot)/30), ' months (primary)']);
    xlabel('age (years)');
    grid on
    % axis([0 P.age_max/365 0 1.1]);
    xlim([0 5])
%     saveas(gcf,['Results/seasonality_',num2str(t(iplot)/30),'.png'])

    subplot(1,2,2)
    plot(a/365,SHb(:,iplot)./PH(:,iplot),'-','Color',colour_mat1); hold on;
    plot(a/365,EHb(:,iplot)./PH(:,iplot),'--','Color',colour_mat6);
    plot(a/365,DHb(:,iplot)./PH(:,iplot),'-.','Color',colour_mat2);
    plot(a/365,AHb(:,iplot)./PH(:,iplot),':','Color',colour_mat3);
    plot(a/365,VHb(:,iplot)./PH(:,iplot),':','Color',colour_mat5);
    plot(a/365,UHb(:,iplot)./PH(:,iplot),':','Color',colour_mat4);
    plot(a/365,(AHb(:,iplot)+DHb(:,iplot))./PH(:,iplot),'r-.');
    plot(a/365,PH(:,iplot)./PH(:,iplot),'-k');
    % plot(a/365,MH(:,end)./PH_final,'r-');
    hold off
    legend('SH/PH','EH/PH','DH/PH', 'AH/PH', 'VH/PH','UH/PH','(AH+DH)/PH');
    title(['time = ', num2str(t(iplot)/30), ' months (vacc.)']);
    xlabel('age (years)');
    grid on
    axis([0 P.age_max/365 0 1.1]);
    xlim([0 5])
%     saveas(gcf,['Results/seasonality_',num2str(t(iplot)/30),'.png'])
    pause
end  