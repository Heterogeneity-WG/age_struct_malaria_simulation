function Malaria_parameters_transform_eFAST
% reduced version of the Malaria_parameters_transform.m for running the eFAST. 
global P

a = P.a;

P.c2 = P.c1; % weight for maternal immunity
P.c3 = P.c1; % weight for vaccine-derived immunity

P.cS = (1-2.5*P.cX)/2; % SH weight
P.cE = P.cX; % EH weight  ~~ AH
P.cA = P.cX; % AH weight
P.cD = 0.5*P.cX; % DH weight
P.cU = P.cS; % UH weight ~~ SH
P.cV = P.cS; % weight for vaccination ~~ SH

P.rho = sigmoid_prob(zeros(size(a)), 'rho');
P.phi = sigmoid_prob(zeros(size(a)), 'phi');
P.psi = sigmoid_prob(zeros(size(a)), 'psi');

%% vaccination functions
Malaria_parameters_transform_vac;

end

