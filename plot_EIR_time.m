pop = trapz(PH,1)*P.da;
figure_setups;
plot(t/365,EIR(end,:));
hold on;
plot(t/365,EIR(floor(end/10),:));
plot(t/365,EIR(floor(end/100),:));

plot(t/365,trapz(PH.*EIR,1)*P.da./pop);
legend('EIR age 100', 'EIR age 10','EIR age 1','Average EIR (pop. weighted)');
xlabel('Time (years)');
title('aEIR dynamics');

%%
figure_setups;
% plot(t/365,round(365/dt)*movmean(EIR(end,:),round(365/dt)));
% hold on;
% plot(t/365,EIR(floor(end/10),:));
% plot(t/365,EIR(floor(end/100),:));
% pop = trapz(PH,1)*P.da;
% 
pop_avg_EIR = trapz(PH.*EIR,1)*P.da./pop;
plot(t/365,movmean(pop_avg_EIR,[round(365/dt) 0]));
title('Moving average population EIR');
% 
% legend('EIR age 100', 'EIR age 10','EIR age 1','Average EIR (pop. weighted)');
% xlabel('Time (years)');
% title('aEIR dynamics');


