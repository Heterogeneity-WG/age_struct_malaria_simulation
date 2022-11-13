function [SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH,SHv, EHv, DHv, AHv, VHv, UHv,SHc, EHc, DHc, AHc, Ctotv, Ctotc] = ...
    age_structured_Malaria_eff2(da, na, tfinal, ...
    SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0,SHv0, EHv0, DHv0, AHv0, VHv0, UHv0,SHc0, EHc0, DHc0, AHc0)
global P

dt = da;
t = (0:dt:tfinal)';
nt = length(t);

% allocation
% SH, EH, etc.. = cell averages
SH = NaN(na,nt); EH = NaN(na,nt); DH = NaN(na,nt); AH = NaN(na,nt); VH = NaN(na,nt); UH = NaN(na,nt); MH = NaN(na,nt);
SM = NaN(1,nt); EM = NaN(1,nt); IM = NaN(1,nt);
Cm = NaN(na,nt); Cac = NaN(na,nt); Cv = NaN(na,nt); Ctot = NaN(na,nt);
NH = NaN(1,nt); NM = NaN(nt);

SH(:,1) = SH0; EH(:,1) = EH0; DH(:,1) = DH0; AH(:,1) = AH0; VH(:,1) = VH0; UH(:,1) = UH0; MH(:,1) = MH0; 
SM(1) = SM0; EM(1) = EM0; IM(1) = IM0;
Cm(:,1) = Cm0; Cac(:,1) = Cac0; Cv(:,1) = Cv0; Ctot(:,1) = Ctot0;
NH(1) = trapz(SH(:,1)+EH(:,1)+DH(:,1)+AH(:,1)+VH(:,1)+UH(:,1))*da;
NM(1) = SM(1)+EM(1)+IM(1);

% diagnostic for vaccinated cohort
SHv = NaN(na,nt); EHv = NaN(na,nt); DHv = NaN(na,nt); AHv = NaN(na,nt); VHv = NaN(na,nt); UHv = NaN(na,nt);
Cmv = NaN(na,nt); Cacv = NaN(na,nt); Cvv = NaN(na,nt); Ctotv = NaN(na,nt);

SHv(:,1) = SHv0; EHv(:,1) = EHv0; DHv(:,1) = DHv0; AHv(:,1) = AHv0; VHv(:,1) = VHv0; UHv(:,1) = UHv0; 
Cmv(:,1) = Cm0; Cacv(:,1) = Cac0; Cvv(:,1) = Cv0; Ctotv(:,1) = Ctot0;

% diagnostic for control cohort
SHc = NaN(na,nt); EHc = NaN(na,nt); DHc = NaN(na,nt); AHc = NaN(na,nt); 
Cmc = NaN(na,nt); Cacc = NaN(na,nt); Cvc = NaN(na,nt); Ctotc = NaN(na,nt);

SHc(:,1) = SHc0; EHc(:,1) = EHc0; DHc(:,1) = DHc0; AHc(:,1) = AHc0; 
Cmc(:,1) = Cm0; Cacc(:,1) = Cac0; Cvc(:,1) = Cv0; Ctotc(:,1) = Ctot0;

% update progression probability based on immunity Ctot
PH0 = SH(:,1)+EH(:,1)+DH(:,1)+AH(:,1)+VH(:,1)+UH(:,1);
P.phi = sigmoid_prob(Ctot0./PH0, 'phi'); % prob. of DH -> RH
P.rho = sigmoid_prob(Ctot0./PH0, 'rho'); % prob. of severely infected, EH -> DH
P.psi = sigmoid_prob(Ctot0./PH0, 'psi'); % prob. AH -> DH
            
P.phiv = P.phi; P.rhov = P.rho; P.psiv = P.psi;
P.phic = P.phi; P.rhoc = P.rho; P.psic = P.psi;


%% time evolution
for n = 1:nt-1
    PH = SH(:,n)+EH(:,n)+DH(:,n)+AH(:,n)+VH(:,n)+UH(:,n); % total human at age a, t = n
    NH(n) = trapz(PH)*da; % total human population at t=n;
    NM(n) = SM(1,n)+EM(1,n)+IM(1,n);
    [bH,~] = biting_rate(NH(n),NM(n));
    lamH = FOI_H(bH,IM(1,n),NM(n));  % force of infection at t=n
    
    % human birth terms
    SH(1,n+1) = trapz(P.gH.*PH)*da;
    EH(1,n+1) = 0;
    DH(1,n+1) = 0;
    AH(1,n+1) = 0;
    VH(1,n+1) = 0;
    UH(1,n+1) = 0;
    MH(1,n+1) = 0; % diagnostic, disease-induced mortality
    SHv(1,n+1) = 0;
    EHv(1,n+1) = 0;
    DHv(1,n+1) = 0;
    AHv(1,n+1) = 0;
    VHv(1,n+1) = 0;
    UHv(1,n+1) = 0;
    SHc(1,n+1) = 0;
    EHc(1,n+1) = 0;
    DHc(1,n+1) = 0;
    AHc(1,n+1) = 0;
    
    %% human time evolution
    SH(2:end,n+1) = (SH(1:end-1,n)+dt*(P.phi(1:end-1)*P.rD.*DH(1:end-1,n)+P.rA*AH(1:end-1,n)))...
        ./(1+(lamH+P.v(2:end)+P.muH(2:end))*dt); 
    VH(2:end,n+1) = (VH(1:end-1,n)+dt*P.etas*(1-P.z)*P.v(2:end).*SH(2:end,n+1))...
        ./(1+(P.muH(2:end)+P.w)*dt);
    temp1 = (P.z*P.v(2:end)+(1-P.etas)*(1-P.z)*P.v(2:end)).*SH(2:end,n+1)+P.w*VH(2:end,n+1);
    UH(2:end,n+1) = (UH(1:end-1,n)+dt*temp1)./(1+(lamH+P.muH(2:end))*dt);
    EH(2:end,n+1) = (EH(1:end-1,n)+dt*lamH*(SH(2:end,n+1)+UH(2:end,n+1)))...
        ./(1+(P.h+P.muH(2:end))*dt);
    temp2 = (1-P.rho(1:end-1))*P.h.*EH(2:end,n+1)+(1-P.phi(1:end-1)).*P.rD.*DH(1:end-1,n);
    AH(2:end,n+1) = ((1-dt*P.rA)*AH(1:end-1,n)+dt*temp2)...
        ./(1+dt*(P.psi(1:end-1)*lamH+P.muH(2:end)));
    temp3 = P.rho(1:end-1)*P.h.*EH(2:end,n+1)+P.psi(1:end-1).*lamH.*AH(2:end,n+1);
    DH(2:end,n+1) = ((1-dt*P.rD)*DH(1:end-1,n)+dt*temp3)...
        ./(1+dt*(P.muH(2:end)+P.muD(2:end)));
    %% diagnostic, disease-induced mortality counts
    MH(2:end,n+1) = MH(1:end-1,n)+dt*P.muD(2:end).*DH(2:end,n+1);
    %% vaccinated cohort
    SHv(2:end,n+1) = (SHv(1:end-1,n)+dt*(P.phiv(1:end-1)*P.rD.*DHv(1:end-1,n)+P.rA*AHv(1:end-1,n)))...
        ./(1+(lamH+P.muH(2:end))*dt); 
    VHv(2:end,n+1) = (VHv(1:end-1,n)+dt*P.etas*(1-P.z)*P.v(2:end).*SH(2:end,n+1))...
        ./(1+(P.muH(2:end)+P.w)*dt);
    temp1 = (P.z*P.v(2:end)+(1-P.etas)*(1-P.z)*P.v(2:end)).*SH(2:end,n+1)+P.w*VHv(2:end,n+1);
    UHv(2:end,n+1) = (UHv(1:end-1,n)+dt*temp1)./(1+(lamH+P.muH(2:end))*dt);
    EHv(2:end,n+1) = (EHv(1:end-1,n)+dt*lamH*(SHv(2:end,n+1)+UHv(2:end,n+1)))...
        ./(1+(P.h+P.muH(2:end))*dt);
    temp2 = (1-P.rhov(1:end-1))*P.h.*EHv(2:end,n+1)+(1-P.phiv(1:end-1)).*P.rD.*DHv(1:end-1,n);
    AHv(2:end,n+1) = ((1-dt*P.rA)*AHv(1:end-1,n)+dt*temp2)...
        ./(1+dt*(P.psiv(1:end-1)*lamH+P.muH(2:end)));
    temp3 = P.rhov(1:end-1)*P.h.*EHv(2:end,n+1)+P.psiv(1:end-1).*lamH.*AHv(2:end,n+1);
    DHv(2:end,n+1) = ((1-dt*P.rD)*DHv(1:end-1,n)+dt*temp3)...
        ./(1+dt*(P.muH(2:end)+P.muD(2:end)));
    %% control group
    SHc(2:end,n+1) = (SHc(1:end-1,n)+dt*(P.phic(1:end-1)*P.rD.*DHc(1:end-1,n)+P.rA*AHc(1:end-1,n)+P.v(2:end).*SH(2:end,n+1)))...
        ./(1+(lamH+P.muH(2:end))*dt);      
    EHc(2:end,n+1) = (EHc(1:end-1,n)+dt*lamH*(SHc(2:end,n+1)))...
        ./(1+(P.h+P.muH(2:end))*dt);
    temp2 = (1-P.rhoc(1:end-1))*P.h.*EHc(2:end,n+1)+(1-P.phic(1:end-1)).*P.rD.*DHc(1:end-1,n);
    AHc(2:end,n+1) = ((1-dt*P.rA)*AHc(1:end-1,n)+dt*temp2)...
        ./(1+dt*(P.psic(1:end-1)*lamH+P.muH(2:end)));
    temp3 = P.rhoc(1:end-1)*P.h.*EHc(2:end,n+1)+P.psic(1:end-1).*lamH.*AHc(2:end,n+1);
    DHc(2:end,n+1) = ((1-dt*P.rD)*DHc(1:end-1,n)+dt*temp3)...
        ./(1+dt*(P.muH(2:end)+P.muD(2:end)));

    %%
    PHp1 = SH(:,n+1)+EH(:,n+1)+DH(:,n+1)+AH(:,n+1)+VH(:,n+1)+UH(:,n+1); % total human at age a, t = n+1
    NHp1 = trapz(PHp1)*da; % total human population at t=n+1;
    % mosquito time evolution
    [SM(1,n+1),EM(1,n+1),IM(1,n+1)] = mosquito_ODE(SM(1,n), EM(1,n), IM(1,n), DH(:,n), AH(:,n), NH(n), NHp1, NM(n));
    NM(n+1) = SM(1,n+1)+EM(1,n+1)+IM(1,n+1);
    
    % immunity gained at age = 0 
    Cm(1,n+1) = P.m*trapz(P.gH.*(P.c1*Cac(:,n)+P.c3*Cv(:,n)))*da;
    Cac(1,n+1) = 0;
    Cv(1,n+1) = P.cV*P.etab*P.z*P.v(1)*SH(1,n+1);
    
    Cmv(1,n+1) = 0; Cacv(1,n+1) = 0; Cvv(1,n+1) = 0;
    Cmc(1,n+1) = 0; Cacc(1,n+1) = 0; Cvc(1,n+1) = 0;
    
    % acquired immunity - use Qn+1
    [bHp1,~] = biting_rate(NHp1,NM(n+1));
    lamHp1 = FOI_H(bHp1,IM(1,n+1),NM(n+1));
    % Cm and Cac are both pooled immunity
    Bnp1 = f(lamHp1).*(P.cS*SH(2:end,n+1) + P.cE*EH(2:end,n+1) + P.cA*AH(2:end,n+1) + P.cD*DH(2:end,n+1) + P.cU*UH(2:end,n+1));
    Dnp1 = P.muH(2:end) + P.muD(2:end).*DH(2:end,n+1)./PHp1(2:end);
    Cac(2:end,n+1) = (Cac(1:end-1,n)+dt*Bnp1)./(1+dt*(1/P.dac+Dnp1));
    Cm(2:end,n+1) = Cm(1:end-1,n)./(1+dt*(1/P.dm+Dnp1));
    Cv(2:end,n+1) = (Cv(1:end-1,n)+P.cV*P.etab*P.z*P.v(2:end).*SH(2:end,n+1))...
        ./(1+dt*(1/P.dv+Dnp1));
    Ctot(:,n+1) = P.c1*Cac(:,n+1)+P.c2*Cm(:,n+1)+P.c3*Cv(:,n+1); % total immunity from acquired, maternal and vaccine-derived sources
    % 
    PHp1v = SHv(:,n+1)+EHv(:,n+1)+DHv(:,n+1)+AHv(:,n+1)+VHv(:,n+1)+UHv(:,n+1); 
    Bnp1v = f(lamHp1).*(P.cS*SHv(2:end,n+1) + P.cE*EHv(2:end,n+1) + P.cA*AHv(2:end,n+1) + P.cD*DHv(2:end,n+1) + P.cU*UHv(2:end,n+1));
    Dnp1v = P.muH(2:end) + P.muD(2:end).*DHv(2:end,n+1)./PHp1v(2:end);
    Cacv(2:end,n+1) = (Cacv(1:end-1,n)+dt*Bnp1v)./(1+dt*(1/P.dac+Dnp1v));
    Cmv(2:end,n+1) = Cmv(1:end-1,n)./(1+dt*(1/P.dm+Dnp1v));
    Cvv(2:end,n+1) = (Cvv(1:end-1,n)+P.cV*P.etab*P.z*P.v(2:end).*SH(2:end,n+1))...
        ./(1+dt*(1/P.dv+Dnp1v));
    Cmv(PHp1v==0,n+1) = 0; Cacv(PHp1v==0,n+1) = 0; Cvv(PHp1v==0,n+1) = 0;
    Ctotv(:,n+1) = P.c1*Cacv(:,n+1)+P.c2*Cmv(:,n+1)+P.c3*Cvv(:,n+1); 
    
    %
    PHp1c = SHc(:,n+1)+EHc(:,n+1)+DHc(:,n+1)+AHc(:,n+1);
    Bnp1c = f(lamHp1).*(P.cS*SHc(2:end,n+1) + P.cE*EHc(2:end,n+1) + P.cA*AHc(2:end,n+1) + P.cD*DHc(2:end,n+1));
    Dnp1c = P.muH(2:end) + P.muD(2:end).*DHc(2:end,n+1)./PHp1c(2:end);
    Cacc(2:end,n+1) = (Cacc(1:end-1,n)+dt*Bnp1c)./(1+dt*(1/P.dac+Dnp1c));
    Cmc(2:end,n+1) = Cmc(1:end-1,n)./(1+dt*(1/P.dm+Dnp1c));
    Cvc(2:end,n+1) = Cvc(1:end-1,n)./(1+dt*(1/P.dv+Dnp1c));
    Cmc(PHp1c==0,n+1) = 0; Cacc(PHp1c==0,n+1) = 0; Cvc(PHp1c==0,n+1) = 0;
    Ctotc(:,n+1) = P.c1*Cacc(:,n+1)+P.c2*Cmc(:,n+1)+P.c3*Cvc(:,n+1);     
    
    % update progression probability based on immunity Ctot
    P.phi = sigmoid_prob(Ctot(:,n+1)./PHp1, 'phi'); % prob. of DH -> RH
    P.rho = sigmoid_prob(Ctot(:,n+1)./PHp1, 'rho'); % prob. of severely infected, EH -> DH
    P.psi = sigmoid_prob(Ctot(:,n+1)./PHp1, 'psi'); % prob. AH -> DH
    
    P.phiv = sigmoid_prob(Ctotv(:,n+1)./PHp1v, 'phi'); % prob. of DH -> RH
    P.rhov = sigmoid_prob(Ctotv(:,n+1)./PHp1v, 'rho'); % prob. of severely infected, EH -> DH
    P.psiv = sigmoid_prob(Ctotv(:,n+1)./PHp1v, 'psi'); % prob. AH -> DH
    P.phiv(PHp1v==0) = 0; P.rhov(PHp1v==0) = 0; P.psiv(PHp1v==0) = 0;
    
    P.phic = sigmoid_prob(Ctotc(:,n+1)./PHp1c, 'phi'); % prob. of DH -> RH
    P.rhoc = sigmoid_prob(Ctotc(:,n+1)./PHp1c, 'rho'); % prob. of severely infected, EH -> DH
    P.psic = sigmoid_prob(Ctotc(:,n+1)./PHp1c, 'psi'); % prob. AH -> DH
    P.phic(PHp1c==0) = 0; P.rhoc(PHp1c==0) = 0; P.psic(PHp1c==0) = 0;

end

end
