figure_setups;
xx = 0:0.01:10;
rho = sigmoid_prob(xx, 'rho');
phi = sigmoid_prob(xx, 'phi');
psi = sigmoid_prob(xx, 'psi');
plot(xx,rho, xx, phi, xx, psi);
legend('rho','phi','psi')

%% averaged values of sigmoids
Ctot_pp = Ctot./PH;
PH_norm = PH./NH;

phi_ave_temp = trapz(sigmoid_prob(Ctot_pp, 'phi').*PH_norm ,1)*P.da;
phi_ave = mean(phi_ave_temp((end-round(365/dt)):end));

rho_ave_temp = trapz(sigmoid_prob(Ctot_pp, 'rho').*PH_norm ,1)*P.da;
rho_ave = mean(rho_ave_temp((end-round(365/dt)):end));

psi_ave_temp = trapz(sigmoid_prob(Ctot_pp, 'psi').*PH_norm ,1)*P.da;
psi_ave = mean(psi_ave_temp((end-round(365/dt)):end));
