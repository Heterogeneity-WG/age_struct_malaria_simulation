%% plot death incidence by age groups
ind0924m = age_range_ind(a,9/12,24/12);
ind0210y = age_range_ind(a,2,10);

figure_setups;
hold on
plot(t/365,trapz(P.muD(ind0924m:end).*DH(ind0924m:end,:),1)*da);
plot(t/365,trapz(P.muD(ind0210y:end).*DH(ind0210y:end,:),1)*da);
plot(t/365,trapz(P.muD(ind0210y(end)+1:end).*DH(ind0210y(end)+1:end,:),1)*da);
legend('9-24','2-10','10+')
title('death incidence by age groups')


function ind = age_range_ind(a,a_start,a_end)

[~,ind1] = min(abs(a-a_start*365)); % start from a_start years old
[~,ind2] = min(abs(a-a_end*365)); % end at a_end years old

ind = ind1:ind2;
end