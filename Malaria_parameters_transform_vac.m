function Malaria_parameters_transform_vac

global P

%% vaccination functions (baseline = RTS,S trial data)
pd = makedist("Triangular","a",7*30,"b",9*30,"c",19*30);
v = pdf(pd,P.a)*P.v0;
P.v = v;

end