%% calculate cases per person per year
lamH = FOI_H(bH,IM,NM);
rho = sigmoid_prob(Ctot./PH, 'rho'); % prob. of severely infected, EH -> DH
psi = sigmoid_prob(Ctot./PH, 'psi'); % prob. AH -> DH
temp1 = psi.*lamH.*AH; % AH -> DH, number of new symptomatic cases - rate
temp2 = rho.*P.h.*EH;  % EH -> DH, number of new symptomatic cases - rate
[~,ind1] = min(abs(P.a-5*30));
[~,ind2] = min(abs(P.a-17*30));
% ind1 = 1;
% ind2 = length(P.a);
cases_rate1 = trapz(temp1(ind1:ind2,:),1)*P.da;
cases_rate2 = trapz(temp2(ind1:ind2,:),1)*P.da;
cases = cases_rate1+cases_rate2;
pop = trapz(PH(ind1:ind2,:),1)*P.da;
ind11 = length(t);
[~,ind22] = min(abs(t-(t(end)-365)));
cases_pp_py = trapz(cases(ind22:ind11))*P.dt/mean(pop(ind22:ind11));

%% plot cumulative cases
% figure_setups; hold on;
% plot(t/365,cumsum(cases_rate1)*dt);
% plot(t/365,cumsum(cases_rate2)*dt);
% title('Cumulative new (cohort) cases');
% legend('AH to DH','EH to DH');

%% Incidence per person per year
% corresponding to Figure S1 in White et al. (2015)

% figure_setups;
% subplot(1,2,1)
% bar(t'/365,(cases_rate1+cases_rate2).*30); 
% hold on
% ylim([0 10^5])
% plot(tfinal/365*ones(2),ylim,'m-')
% % plot((tfinal+tfinal_conti)/365*ones(2),ylim,'m-')
% xlabel('year')
% ylabel('total cases (per 30days)')
% 
% subplot(1,2,2)
% plot(t/365,cases./pop*365)
% ylim([0, 8])
% hold on
% plot(tfinal/365*ones(2),ylim,'m-')
% % plot((tfinal+tfinal_conti)/365*ones(2),ylim,'m-')
% xlabel('year')
% ylabel('Incidence pp per year')
% title(['Incidence pp per year = ', num2str(cases_pp_py)]);

%% Heatmaps of population disease burden over 3 seasons (years)
subfigure_strings1 = ["(N/A)","(A)","(B)"];

% figure_setups; 
% imagesc(t/365,a/365,AH./PH);
% clim([0 1]);
% colorbar;
% title('Fraction of asymptomatic infections');
% xlabel('Time (years)');
% ylabel('Age (years)');
% xlim([0 3]);
% ylim([0 20]);
% set(gca,'YDir','normal')
% colormap jetwhite;
% grid off;
% yticks([0 5 10 15 20]);
% xticks([0 1 2 3]);

figure_setups_3;
ax = imagesc(t/365,a/365,DH./(AH+DH));
clim([0 0.6]);
colorbar;
xlim([0 3]);
ylim([0 20]);
xlabel('Time (years)');
ylabel('Age (years)');
%title('Symptomatic infections');
set(gca,'YDir','normal')
colormap jet;
str_temp = subfigure_strings1(immunity_feedback);
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
% this sets an 'a)' right at the top left of the axes
text(ax,0,temp_ax(end)+3,str_temp,'Units','characters');
grid off;
yticks([0 5 10 15 20]);
xticks([0 1 2 3]);
save_string = strcat('fig2_',str_temp,'.svg');
saveas(gcf,save_string);
%%
subfigure_strings2 = ["(N/A)","(C)","(D)"];
str_temp = subfigure_strings2(immunity_feedback);

[~,age2] = min(abs(P.a-2*365));
[~,age10] = min(abs(P.a-10*365));
[~,age20] = min(abs(P.a-20*365));

figure_setups_3;
temp_plot = AH./(PH); 
plot(t/365,temp_plot(age2,:),'Color',[0.9290, 0.6940, 0.1250]);
hold on;
plot(t/365,temp_plot(age10,:),'-.','Color',[0.9290, 0.6940, 0.1250]);
plot(t/365,temp_plot(age20,:),':','Color',[0.9290, 0.6940, 0.1250]);

temp_plot = DH./(PH); 
plot(t/365,temp_plot(age2,:),'Color',[0.8500, 0.3250, 0.0980]);
plot(t/365,temp_plot(age10,:),'-.','Color',[0.8500, 0.3250, 0.0980]);
plot(t/365,temp_plot(age20,:),':','Color',[0.8500, 0.3250, 0.0980]);
hold off;

xlabel('Time (years)');
ylim([0 0.6]);
xlim([0 3]);
yticks([0 0.2 0.4 0.6 0.8 1]);
xticks([0 1 2 3]);
ylabel('Proportion infected');
grid on;
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
% this sets an 'a)' right at the top left of the axes
text(ax,0,temp_ax(end)+3,str_temp,'Units','characters');
% legend('Age 2 (asymptomatic)','Age 10 (asymptomatic)','Age 20 (asymptomatic)'...
%      ,'Age 2 (symptomatic)','Age 10 (symptomatic)','Age 20 (symptomatic)'...
%      ,'NumColumns', 2);
save_string = strcat('fig2_',str_temp,'.svg');
saveas(gcf,save_string);