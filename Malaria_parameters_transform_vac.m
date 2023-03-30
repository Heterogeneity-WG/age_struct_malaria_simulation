function Malaria_parameters_transform_vac

global P

%% vaccination functions (baseline = RTS,S trial data)
% parameter based on the RTS,S trial report, children (age 5–17 months) doi.org/10.1016/S0140-6736(15)60721-8.
% with a mode at 9month old (Kenya brochure)

pd = makedist("Triangular","a",5*30,"b",9*30,"c",17*30);
% pd = makedist("Uniform","lower",15*365,"upper",50*365);

v = pdf(pd,P.a)*P.v0;
P.v = v;

vs = pdf(pd,P.a)*P.v0s;
P.vs = vs; % for vaccine efficacy calculation only


vc = pdf(pd,P.a)*P.v0c;
P.vc = vc; % for vaccine efficacy calculation only

end