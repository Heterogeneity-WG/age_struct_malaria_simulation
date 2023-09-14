figure_setups;
xx = 0:0.01:10;
rho = sigmoid_prob(xx, 'rho');
phi = sigmoid_prob(xx, 'phi');
psi = sigmoid_prob(xx, 'psi');
plot(xx,rho, xx, phi, xx, psi);
legend('rho','phi','psi')