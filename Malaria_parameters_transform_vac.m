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
P.vs = vs; % vacc rate from Rest -> Vacc patch

vc = pdf(pd,P.a)*P.v0c;
P.vc = vc; % vacc rate from Rest -> Control patch

% booster dose: 18 months after 3rd dose [26, 38] months old
pd_booster = makedist("Triangular","a",26*30,"b",27*30,"c",38*30);
P.vb = pdf(pd_booster,P.a)*P.vb0; 
% pd_booster = @(age) (age>38*30).*(age<39*30); % debugging only
% P.vb = pd_booster(P.a)*P.vb0; 
end