%% This ODE represents the HIV model in Section 4.2
function dydt=ODE_efast(t,y,X,run_num)

%% PARAMETERS -- baseline
Parameter_settings_EFAST;

% update parameter with sampled values
s=X(run_num,1);
muT=X(run_num,2);
r=X(run_num,3);
k1=X(run_num,4);
k2=X(run_num,5);
mub=X(run_num,6);
N=X(run_num,7);
muV=X(run_num,8);
dummy=X(run_num,9);

% [T] CD4+ uninfected: Tsource + Tprolif - Tinf
Tsource = s - muT*y(1);
Tprolif = r*y(1)*(1-(y(1)+y(2)+y(3))/Tmax);
Tinf = k1*y(1)*y(4);

% [T1] CD4+ latently infected: Tinf - T1death - T1inf
T1death = muT*y(2);
T1inf = k2*y(2);

% [T2] CD4+ actively infected: T1inf - T2death
T2death = mub*y(3);

% [V] Free infectious virus: Vrelease - Tinf - Vdeath
Vrelease = N*T2death;
Vdeath = muV*y(4);

dydt = [Tsource + Tprolif - Tinf;
        Tinf - T1death - T1inf;
        T1inf - T2death;
        Vrelease - Tinf - Vdeath];