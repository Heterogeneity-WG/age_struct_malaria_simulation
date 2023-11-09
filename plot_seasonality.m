% seasonlity plots
figure_setups_2;
yyaxis left
hold on 
plot(t/365,EIR_tot);
death_rate_tot = trapz(P.muD.*DH,1)*da;
plot(t/365,death_rate_tot); 
ylim([0 120])
yyaxis right
hold on
% seasonal_profile = P.ss_S;
% plot(t/365, seasonal_profile(t))
infected_tot = trapz(AH+DH,1)*P.da./NH;
infected_AH =  trapz(AH,1)*P.da./NH;
infected_DH =  trapz(DH,1)*P.da./NH;
plot(t/365,infected_tot)
plot(t/365,infected_AH)
plot(t/365,infected_DH)
% ylim([0.4 0.7])
legend('EIR','death incidence rate','$A_H+D_H$ prop.','$A_H$ prop.','$D_H$ prop.','Location','se')
xlabel('years')
% legend('EIR','death incidence rate','seasonal profile','$A_H+D_H$ prop.','Location','se')
% save('seasonality.mat','t','EIR_tot','death_rate_tot','seasonal_profile','infected_tot')
