%% Immunity dynamics for Ctot
figure_setups;
nt = length(t);
subplot(2,2,1), plot(a/365,Ctot(:,floor(nt/40))./PH(:,floor(nt/40)));
hold on;
subplot(2,2,1), plot(a/365,Ctot(:,floor(3*nt/40))./PH(:,floor(3*nt/40)));
subplot(2,2,1), plot(a/365,Ctot(:,floor(nt/20))./PH(:,floor(nt/20)));
subplot(2,2,1), plot(a/365,Ctot(:,floor(nt/10))./PH(:,floor(nt/10)));
%subplot(2,2,1), plot(a/365,Ctot(:,end)./PH(:,end),'-.');
xlabel('age')
legend(['t = ',num2str(t(end)/(40*365))],['t = ',num2str(t(end)/(20*365))],...
    ['t = ',num2str(3*t(end)/(40*365))],['t = ',num2str(t(end)/(10*365))],'Location','SouthEast');
title('$C_{total}(t)/P_H(t)$');
grid on
subplot(2,2,2), plot(t/365,(trapz(Ctot,1)*da)./NH);
title('$\int C_{total}(\alpha,t)d\alpha / N_H(t)$');
xlabel('time');
grid on

subplot(2,2,3), imagesc(t/365,a/365,Ctot./PH);
set(gca,'YDir','normal');
colorbar;
ylabel('age');
xlabel('Time (years)');
title('$C_{total}(\alpha,t)/P_H(\alpha,t)$');

Ctot_norm = Ctot./PH;
subplot(2,2,4), plot(t/365,Ctot_norm(floor(na/100),:));
hold on;
subplot(2,2,4), plot(t/365,Ctot_norm(floor(na/20),:));
subplot(2,2,4), plot(t/365,Ctot_norm(floor(na/10),:));
subplot(2,2,4), plot(t/365,Ctot_norm(na,:),'-.');
xlabel('Time(years)')
legend('age 1','age 5','age 10', 'age 100'); % assume age_max=100
title('$\tilde{C}_{total}(t,\cdot)$');
 
%% immunity level per age group
ind0924m = age_range_ind(a,9/12,24/12);
ind0210y = age_range_ind(a,2,10);

figure_setups_2; hold on
plot(t/365,trapz(Cac(1:ind0924m(end)-1,:),1)./trapz(PH(1:ind0924m(end)-1,:),1));
plot(t/365,trapz(Cac(ind0924m:end,:),1)./trapz(PH(ind0924m:end,:),1));
plot(t/365,trapz(Cac(ind0210y:end,:),1)./trapz(PH(ind0210y:end,:),1));
plot(t/365,trapz(Cac(ind0210y(end)+1:end,:),1)./trapz(PH(ind0210y(end)+1:end,:),1));
plot(t/365,trapz(Cac,1)./trapz(PH,1),'k--');
% ylim([0 5])
legend('0-9', '9-24','2-10','10+','all')
xlabel('years')
ylabel('immunity level')
title('Per-person immunity level')

%% Averaged Immunity breakdown
figure_setups;
plot(a/365,Cac(:,end)./PH_final,'-.');
hold on;
plot(a/365,Cm(:,end)./PH_final,'--');
plot(a/365,Cv(:,end)./PH_final,'--');
plot(a/365,Ctot(:,end)./PH_final,'-');
xlabel('age (years)')
ylabel('immunity level')
legend('Acquired (pp)','Maternal (pp)','Vaccine-derived (pp)','Total (pp)','Location','SouthEast');
% title(['Per-person Immun dist.~~ feedback =',num2str(immunity_feedback)]);
title('Per-person Immune distribution');
axis([0 age_max/365 0 max(Ctot(:,end)./PH_final)*1.1]);
xlim([0 10])
ylim([0 7])
grid on

%% Immunity dynamics for Cv (no longer relevant)
% figure_setups;
% subplot(2,2,1), plot(a/365,Cv(:,floor(nt/4))./PH(:,floor(nt/4)));
% hold on;
% subplot(2,2,1), plot(a/365,Cv(:,floor(nt/2))./PH(:,floor(nt/2)));
% subplot(2,2,1), plot(a/365,Cv(:,floor(3*nt/4))./PH(:,floor(3*nt/4)));
% subplot(2,2,1), plot(a/365,Cv(:,end)./PH(:,end),'-.');
% xlabel('age')
% legend(['t = ',num2str(t(end)/(4*365))],['t = ',num2str(t(end)/(2*365))],...
%     ['t = ',num2str(3*t(end)/(4*365))],['t = ',num2str(t(end)/365)],'Location','NorthWest');
% title('$C_{v}(t)/P_H(t)$');
% grid on
% subplot(2,2,2), plot(t/365,(trapz(Cv,1)*da)./NH);
% title('$\int C_{v}(\alpha,t)d\alpha / N_H(t)$');
% xlabel('time');
% grid on
% 
% subplot(2,2,3), imagesc(t/365,a/365,Cv./PH);
% set(gca,'YDir','normal');
% colorbar;
% ylabel('age');
% xlabel('time');
% title('$C_{v}(\alpha,t)/P_H(\alpha,t)$');
% 
% subplot(2,2,4), plot(a/365,Cv(:,floor(nt/4)));
% hold on;
% subplot(2,2,4), plot(a/365,Cv(:,floor(nt/2)));
% subplot(2,2,4), plot(a/365,Cv(:,floor(3*nt/4)));
% subplot(2,2,4), plot(a/365,Cv(:,end),'-.');
% xlabel('age')
% legend(['t = ',num2str(t(end)/(4*365))],['t = ',num2str(t(end)/(2*365))],...
%     ['t = ',num2str(3*t(end)/(4*365))],['t = ',num2str(t(end)/365)]);
% title('$C_{v}(t)$');

function ind = age_range_ind(a,a_start,a_end)

[~,ind1] = min(abs(a-a_start*365)); % start from a_start years old
[~,ind2] = min(abs(a-a_end*365)); % end at a_end years old

ind = ind1:ind2;
end