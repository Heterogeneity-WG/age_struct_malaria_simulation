function Malaria_parameters_transform_SA
% reduced version of the Malaria_parameters_transform.m for running the eFAST. 
global P

a = P.a;

P.c2 = P.c1; % weight for maternal immunity
P.c3 = P.c1; % weight for vaccine-derived immunity
P.cV = P.cS; % weight for vaccination ~~ SH

P.rho = sigmoid_prob(zeros(size(a)), 'rho');
P.phi = sigmoid_prob(zeros(size(a)), 'phi');
P.psi = sigmoid_prob(zeros(size(a)), 'psi');

%% vaccination functions
Malaria_parameters_transform_vac;

end

