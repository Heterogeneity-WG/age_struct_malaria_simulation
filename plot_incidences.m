%% calculate cases per person per year
lamH = FOI_H(bH,IM,NM);
rho = sigmoid_prob(Ctot./PH, 'rho'); % prob. of severely infected, EH -> DH
psi = sigmoid_prob(Ctot./PH, 'psi'); % prob. AH -> DH
temp1 = psi.*lamH.*AH; % AH -> DH, number of new symptomatic cases - rate
temp2 = rho.*P.h.*EH;  % EH -> DH, number of new symptomatic cases - rate
[~,ind1] = min(abs(P.a-5*30));
[~,ind2] = min(abs(P.a-17*30));
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
figure_setups;
subplot(1,2,1)
bar(t'/365,(cases_rate1+cases_rate2).*30); 
hold on
ylim([0 10^5])
plot(tfinal/365*ones(2),ylim,'m-')
% plot((tfinal+tfinal_conti)/365*ones(2),ylim,'m-')
xlabel('year')
ylabel('total cases (per 30days)')

subplot(1,2,2)
plot(t/365,cases./pop*365)
ylim([0, 8])
hold on
plot(tfinal/365*ones(2),ylim,'m-')
% plot((tfinal+tfinal_conti)/365*ones(2),ylim,'m-')
xlabel('year')
ylabel('Incidence pp per year')
title(['Incidence pp per year = ', num2str(cases_pp_py)]);
