%% vaccine #
figure_setups; hold on
subplot(1,2,1)
yyaxis left
plot(t/365,cumsum(vacc_blood)*dt/NN_S,t/365,cumsum(vacc_sterile)*dt/NN_S)
legend('Blood-stage','Sterile')
ylabel('fraction')
axis([0 max(t)/365 0 1]);
xlabel('years')
title('cumulative vaccinated')
yyaxis right
plot(t/365,cumsum(vacc_blood)*dt,t/365,cumsum(vacc_sterile)*dt)
ylabel('count')
axis([0 max(t)/365 0 1*NN_S]);
subplot(1,2,2)
vacc_fun = P.v;
plot(P.a/30,vacc_fun)
xlim([0 25])
xlabel('month')
title('$\nu(\alpha)$ daily per-capita vacc. rate')
