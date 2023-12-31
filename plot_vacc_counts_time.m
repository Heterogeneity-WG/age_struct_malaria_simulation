%% vaccine #
[~,ind1] = min(abs(P.a-5*30));
[~,ind2] = min(abs(P.a-17*30));
NN_target = trapz(PH_final(ind1:ind2))*P.da;

figure_setups; 
hold on
subplot(1,2,1)
yyaxis left
plot(t/365,cumsum(trapz(vacc_sterile,1)*P.da)*dt/P.NN)
ylabel('fraction (total pop)')
% axis([0 max(t)/365 0 0.1]);
xlabel('years')
title('cumu. vaccinated')
yyaxis right
plot(t/365,cumsum(trapz(vacc_sterile,1)*P.da)*dt/NN_target)
ylabel('fraction (target pop)')
axis([0 max(t)/365 0 1]);
subplot(1,2,2)
vacc_fun = P.v;
plot(P.a/30,vacc_fun)
xlim([0 25])
xlabel('month')
title('$\nu(\alpha)$ daily vacc. rate')
