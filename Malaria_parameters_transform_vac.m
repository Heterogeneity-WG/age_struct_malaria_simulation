function Malaria_parameters_transform_vac

global P

%% vaccination functions (baseline = RTS,S trial data)
% parameter based on the RTS,S trial report, children   (Penny et al. 2015, Table 1) doi.org/10.1016/S0140-6736(15)60721-8.
% first dose: age 5–17 months with a mode at 6-month old (Kenya brochure, Ministry of Health of Kenya)
% third dose: age 8-20 months with a mode at 9-month old (Kenya brochure, Ministry of Health of Kenya)

pd = makedist("Triangular","a",8*30,"b",9*30,"c",20*30); % age range for the completion of 3rd dose

v = pdf(pd,P.a)*P.v0;
P.v = v;

vs = pdf(pd,P.a)*P.v0s;
P.vs = vs; % for vaccine efficacy calculation only

vc = pdf(pd,P.a)*P.v0c;
P.vc = vc; % for vaccine efficacy calculation only

% booster dose: 18 months after 3rd dose [26, 38] months old
pd_booster = @(age) (age>26*30).*(age<38*30); % age range for the completion of 3rd dose
P.vb = pd_booster(P.a)*P.vb0; % indicator of age * fraction vb0

end