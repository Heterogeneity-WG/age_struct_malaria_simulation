function [t, SHr, EHr, DHr, AHr, Cmr, Cacr, Ctotr, SHv, EHv, DHv, AHv, VHv, UHv, Cmv, Cacv, Ctotv, SHb, EHb, DHb, AHb, VHb, UHb, Cmb, Cacb, Ctotb, SHc, EHc, DHc, AHc, Cmc, Cacc, Ctotc, SM, EM, IM] = ...
    age_structured_Malaria_eff_booster(da, na, t0, tfinal, SHr0, EHr0, DHr0, AHr0, Cmr0, Cacr0, Ctotr0, ...
    SHv0, EHv0, DHv0, AHv0, VHv0, UHv0, Cmv0, Cacv0, Ctotv0, ...
    SHb0, EHb0, DHb0, AHb0, VHb0, UHb0, Cmb0, Cacb0, Ctotb0, ...
    SHc0, EHc0, DHc0, AHc0, Cmc0, Cacc0, Ctotc0, SM0, EM0, IM0)

global P

dt = da;
t = (t0:dt:tfinal)';
nt = length(t);

% allocation
SM = NaN(1,nt); EM = NaN(1,nt); IM = NaN(1,nt);
SM(1) = SM0; EM(1) = EM0; IM(1) = IM0;

NH = NaN(1,nt); NM = NaN(nt);
NM(1) = SM(1)+EM(1)+IM(1);

% Rest group
SHr = NaN(na,nt); EHr = NaN(na,nt); DHr = NaN(na,nt); AHr = NaN(na,nt);
Cmr = NaN(na,nt); Cacr = NaN(na,nt); Cvr = NaN(na,nt);  Ctotr = NaN(na,nt);
SHr(:,1) = SHr0; EHr(:,1) = EHr0; DHr(:,1) = DHr0; AHr(:,1) = AHr0; 
Cmr(:,1) = Cmr0; Cacr(:,1) = Cacr0; Cvr(:,1) = 0; Ctotr(:,1) = Ctotr0;

% vaccinated group
SHv = NaN(na,nt); EHv = NaN(na,nt); DHv = NaN(na,nt); AHv = NaN(na,nt); VHv = NaN(na,nt); UHv = NaN(na,nt);
Cmv = NaN(na,nt); Cacv = NaN(na,nt); Cvv = NaN(na,nt); Ctotv = NaN(na,nt);
SHv(:,1) = SHv0; EHv(:,1) = EHv0; DHv(:,1) = DHv0; AHv(:,1) = AHv0; VHv(:,1) = VHv0; UHv(:,1) = UHv0; 
Cmv(:,1) = Cmv0; Cacv(:,1) = Cacv0; Cvv(:,1) = 0; Ctotv(:,1) = Ctotv0;

% boosted group
SHb = NaN(na,nt); EHb = NaN(na,nt); DHb = NaN(na,nt); AHb = NaN(na,nt); VHb = NaN(na,nt); UHb = NaN(na,nt);
Cmb = NaN(na,nt); Cacb = NaN(na,nt); Cvb = NaN(na,nt); Ctotb = NaN(na,nt);
SHb(:,1) = SHb0; EHb(:,1) = EHb0; DHb(:,1) = DHb0; AHb(:,1) = AHb0; VHb(:,1) = VHb0; UHb(:,1) = UHb0; 
Cmb(:,1) = Cmb0; Cacb(:,1) = Cacb0; Cvb(:,1) = 0; Ctotb(:,1) = Ctotb0;

% control group
SHc = NaN(na,nt); EHc = NaN(na,nt); DHc = NaN(na,nt); AHc = NaN(na,nt); 
Cmc = NaN(na,nt); Cacc = NaN(na,nt); Cvc = NaN(na,nt); Ctotc = NaN(na,nt);
SHc(:,1) = SHc0; EHc(:,1) = EHc0; DHc(:,1) = DHc0; AHc(:,1) = AHc0; 
Cmc(:,1) = Cmc0; Cacc(:,1) = Cacc0; Cvc(:,1) = 0; Ctotc(:,1) = Ctotc0;

PHr = SHr0+EHr0+DHr0+AHr0;
PHc = SHc0+EHc0+DHc0+AHc0;
PHv = SHv0+EHv0+DHv0+AHv0+VHv0+UHv0;
PHb = SHb0+EHb0+DHb0+AHb0+VHb0+UHb0;

% update progression probability based on immunity Ctot
P.phir = sigmoid_prob(Ctotr(:,1)./PHr, 'phi'); % prob. of DH -> RH
P.rhor = sigmoid_prob(Ctotr(:,1)./PHr, 'rho'); % prob. of severely infected, EH -> DH
P.psir = sigmoid_prob(Ctotr(:,1)./PHr, 'psi'); % prob. AH -> DH

P.phiv = sigmoid_prob(Ctotv(:,1)./PHv, 'phi'); % prob. of DH -> RH
P.rhov = sigmoid_prob(Ctotv(:,1)./PHv, 'rho'); % prob. of severely infected, EH -> DH
P.psiv = sigmoid_prob(Ctotv(:,1)./PHv, 'psi'); % prob. AH -> DH
P.phiv(PHv==0) = P.phir(PHv==0); P.rhov(PHv==0) = P.rhor(PHv==0); P.psiv(PHv==0) = P.psir(PHv==0);

P.phib = sigmoid_prob(Ctotb(:,1)./PHb, 'phi'); % prob. of DH -> RH
P.rhob = sigmoid_prob(Ctotb(:,1)./PHb, 'rho'); % prob. of severely infected, EH -> DH
P.psib = sigmoid_prob(Ctotb(:,1)./PHb, 'psi'); % prob. AH -> DH
P.phib(PHb==0) = P.phir(PHb==0); P.rhob(PHb==0) = P.rhor(PHb==0); P.psib(PHb==0) = P.psir(PHb==0);

P.phic = sigmoid_prob(Ctotc(:,1)./PHc, 'phi'); % prob. of DH -> RH
P.rhoc = sigmoid_prob(Ctotc(:,1)./PHc, 'rho'); % prob. of severely infected, EH -> DH
P.psic = sigmoid_prob(Ctotc(:,1)./PHc, 'psi'); % prob. AH -> DH
P.phic(PHc==0) = P.phir(PHc==0); P.rhoc(PHc==0) = P.rhor(PHc==0); P.psic(PHc==0) = P.psir(PHc==0);

%% time evolution
for n = 1:nt-1
    PHr = SHr(:,n)+EHr(:,n)+DHr(:,n)+AHr(:,n); 
    PHc = SHc(:,n)+EHc(:,n)+DHc(:,n)+AHc(:,n); 
    PHv = SHv(:,n)+EHv(:,n)+DHv(:,n)+AHv(:,n)+VHv(:,n)+UHv(:,n); 
    PHb = SHb(:,n)+EHb(:,n)+DHb(:,n)+AHb(:,n)+VHb(:,n)+UHb(:,n); 

    PH = PHr+PHc+PHv+PHb;
    NH(n) = trapz(PH)*da;
    NM(n) = SM(1,n)+EM(1,n)+IM(1,n);
    [bH,~] = biting_rate(PH,NM(n));
    lamH = FOI_H(bH,IM(1,n),NM(n));
    
    % human birth terms
    SHr(1,n+1) = trapz(P.gH.*PH)*da;
    EHr(1,n+1) = 0;
    DHr(1,n+1) = 0;
    AHr(1,n+1) = 0;
    
    SHv(1,n+1) = 0;
    EHv(1,n+1) = 0;
    DHv(1,n+1) = 0;
    AHv(1,n+1) = 0;
    VHv(1,n+1) = 0;
    UHv(1,n+1) = 0;

    SHb(1,n+1) = 0;
    EHb(1,n+1) = 0;
    DHb(1,n+1) = 0;
    AHb(1,n+1) = 0;
    VHb(1,n+1) = 0;
    UHb(1,n+1) = 0;
    
    SHc(1,n+1) = 0;
    EHc(1,n+1) = 0;
    DHc(1,n+1) = 0;
    AHc(1,n+1) = 0;
    
    %% rest group
    SHr(2:end,n+1) = (SHr(1:end-1,n)+dt*(P.phir(1:end-1)*P.rD.*DHr(1:end-1,n)+P.rA*AHr(1:end-1,n))-(P.vc(2:end)+P.vs(2:end))*dt)...
        ./(1+(lamH(1:end-1)+P.muH(2:end))*dt); 
    if min(SHr(2:end,n+1))<0
        keyboard
    end
    EHr(2:end,n+1) = (EHr(1:end-1,n)+dt*lamH(1:end-1).*SHr(2:end,n+1))...
        ./(1+(P.h+P.muH(2:end))*dt);
    temp2 = (1-P.rhor(1:end-1))*P.h.*EHr(2:end,n+1)+(1-P.phir(1:end-1)).*P.rD.*DHr(1:end-1,n);
    AHr(2:end,n+1) = ((1-dt*P.rA)*AHr(1:end-1,n)+dt*temp2)...
        ./(1+dt*(P.psir(1:end-1).*lamH(1:end-1)+P.muH(2:end)));
    temp3 = P.rhor(1:end-1)*P.h.*EHr(2:end,n+1)+P.psir(1:end-1).*lamH(1:end-1).*AHr(2:end,n+1);
    DHr(2:end,n+1) = ((1-dt*P.rD)*DHr(1:end-1,n)+dt*temp3)...
        ./(1+dt*(P.muH(2:end)+P.muD(2:end)));  
    
    %% vaccinated group
    % the maximum vb to avoid SHV<0
    vb_temp = max(min(P.vb(2:end), SHv(1:end-1,n)/dt+(P.phiv(1:end-1)*P.rD.*DHv(1:end-1,n)+P.rA*AHv(1:end-1,n))-10^-5),0);
    SHv(2:end,n+1) = (SHv(1:end-1,n)+dt*(P.phiv(1:end-1)*P.rD.*DHv(1:end-1,n)+P.rA*AHv(1:end-1,n)-vb_temp))...
        ./(1+(lamH(1:end-1)+P.muH(2:end))*dt); 
    if min(SHv(2:end,n+1))<0
        keyboard
    end
    VHv(2:end,n+1) = (VHv(1:end-1,n)+dt*P.etas*P.vs(2:end))./(1+(P.muH(2:end)+P.w)*dt);
    temp1 = (1-P.etas)*P.vs(2:end)+P.w*VHv(2:end,n+1);
    UHv(2:end,n+1) = (UHv(1:end-1,n)+dt*temp1)./(1+(lamH(1:end-1)+P.muH(2:end))*dt);
    EHv(2:end,n+1) = (EHv(1:end-1,n)+dt*lamH(1:end-1).*(SHv(2:end,n+1)+UHv(2:end,n+1)))...
        ./(1+(P.h+P.muH(2:end))*dt);
    temp2 = (1-P.rhov(1:end-1))*P.h.*EHv(2:end,n+1)+(1-P.phiv(1:end-1)).*P.rD.*DHv(1:end-1,n);
    AHv(2:end,n+1) = ((1-dt*P.rA)*AHv(1:end-1,n)+dt*temp2)...
        ./(1+dt*(P.psiv(1:end-1).*lamH(1:end-1)+P.muH(2:end)));
    temp3 = P.rhov(1:end-1)*P.h.*EHv(2:end,n+1)+P.psiv(1:end-1).*lamH(1:end-1).*AHv(2:end,n+1);
    DHv(2:end,n+1) = ((1-dt*P.rD)*DHv(1:end-1,n)+dt*temp3)...
        ./(1+dt*(P.muH(2:end)+P.muD(2:end)));
    %% boosted group
    SHb(2:end,n+1) = (SHb(1:end-1,n)+dt*(P.phib(1:end-1)*P.rD.*DHb(1:end-1,n)+P.rA*AHb(1:end-1,n)))...
        ./(1+(lamH(1:end-1)+P.muH(2:end))*dt); 
    VHb(2:end,n+1) = (VHb(1:end-1,n)+dt*P.etab*vb_temp)./(1+(P.muH(2:end)+P.wb)*dt);
    temp1 = (1-P.etab)*vb_temp+P.wb*VHb(2:end,n+1);
    UHb(2:end,n+1) = (UHb(1:end-1,n)+dt*temp1)./(1+(lamH(1:end-1)+P.muH(2:end))*dt);
    EHb(2:end,n+1) = (EHb(1:end-1,n)+dt*lamH(1:end-1).*(SHb(2:end,n+1)+UHb(2:end,n+1)))...
        ./(1+(P.h+P.muH(2:end))*dt);
    temp2 = (1-P.rhob(1:end-1))*P.h.*EHb(2:end,n+1)+(1-P.phib(1:end-1)).*P.rD.*DHb(1:end-1,n);
    AHb(2:end,n+1) = ((1-dt*P.rA)*AHb(1:end-1,n)+dt*temp2)...
        ./(1+dt*(P.psib(1:end-1).*lamH(1:end-1)+P.muH(2:end)));
    temp3 = P.rhob(1:end-1)*P.h.*EHb(2:end,n+1)+P.psib(1:end-1).*lamH(1:end-1).*AHb(2:end,n+1);
    DHb(2:end,n+1) = ((1-dt*P.rD)*DHb(1:end-1,n)+dt*temp3)...
        ./(1+dt*(P.muH(2:end)+P.muD(2:end)));
    %% control group
    SHc(2:end,n+1) = (SHc(1:end-1,n)+dt*(P.phic(1:end-1)*P.rD.*DHc(1:end-1,n)+P.rA*AHc(1:end-1,n)+P.vc(2:end)))...
        ./(1+(lamH(1:end-1)+P.muH(2:end))*dt);      
    EHc(2:end,n+1) = (EHc(1:end-1,n)+dt*lamH(1:end-1).*SHc(2:end,n+1))./(1+(P.h+P.muH(2:end))*dt);
    temp2 = (1-P.rhoc(1:end-1))*P.h.*EHc(2:end,n+1)+(1-P.phic(1:end-1)).*P.rD.*DHc(1:end-1,n);
    AHc(2:end,n+1) = ((1-dt*P.rA)*AHc(1:end-1,n)+dt*temp2)...
        ./(1+dt*(P.psic(1:end-1).*lamH(1:end-1)+P.muH(2:end)));
    temp3 = P.rhoc(1:end-1)*P.h.*EHc(2:end,n+1)+P.psic(1:end-1).*lamH(1:end-1).*AHc(2:end,n+1);
    DHc(2:end,n+1) = ((1-dt*P.rD)*DHc(1:end-1,n)+dt*temp3)...
        ./(1+dt*(P.muH(2:end)+P.muD(2:end)));

    %%
    PHrp1 = SHr(:,n+1)+EHr(:,n+1)+DHr(:,n+1)+AHr(:,n+1); 
    PHcp1 = SHc(:,n+1)+EHc(:,n+1)+DHc(:,n+1)+AHc(:,n+1); 
    PHvp1 = SHv(:,n+1)+EHv(:,n+1)+DHv(:,n+1)+AHv(:,n+1)+VHv(:,n+1)+UHv(:,n+1); 
    PHbp1 = SHb(:,n+1)+EHb(:,n+1)+DHb(:,n+1)+AHb(:,n+1)+VHb(:,n+1)+UHb(:,n+1); 
    PHp1 = PHrp1+PHcp1+PHvp1+PHbp1;    
    NHp1 = trapz(PHp1)*da; % total human population at t=n+1;
    
    % mosquito time evolution
    P.gM = P.gM_fun(t(n+1)); % incorporate seasonlity
    [SM(1,n+1),EM(1,n+1),IM(1,n+1)] = mosquito_ODE(SM(1,n), EM(1,n), IM(1,n), DHr(:,n)+DHc(:,n)+DHv(:,n)+DHb(:,n), AHr(:,n)+AHc(:,n)+AHv(:,n)+AHb(:,n), PHp1, NHp1, NM(n));
    NM(n+1) = SM(1,n+1)+EM(1,n+1)+IM(1,n+1);
    
    % immunity gained at age = 0 
    Cacn = Cacr(:,n)+Cacv(:,n)+Cacc(:,n)+Cacb(:,n);
    Cvn = Cvr(:,n)+Cvv(:,n)+Cvc(:,n)+Cvb(:,n);
    Cmr(1,n+1) = P.m*trapz(P.gH.*(P.c1*Cacn+P.c3*Cvn))*da;
    Cacr(1,n+1) = 0;
    Cvr(1,n+1) = 0;
    
    Cmv(1,n+1) = 0; Cacv(1,n+1) = 0; Cvv(1,n+1) = 0;
    Cmb(1,n+1) = 0; Cacb(1,n+1) = 0; Cvb(1,n+1) = 0;
    Cmc(1,n+1) = 0; Cacc(1,n+1) = 0; Cvc(1,n+1) = 0;
    
    % acquired immunity - use Qn+1
    [bHp1,~] = biting_rate(PHp1,NM(n+1));
    lamHp1 = FOI_H(bHp1,IM(1,n+1),NM(n+1));
    
    % pooled immunity
    % rest group
    Brnp1 = f(lamHp1(2:end)).*(P.cS*SHr(2:end,n+1) + P.cE*EHr(2:end,n+1) + P.cA*AHr(2:end,n+1) + P.cD*DHr(2:end,n+1));
    Drnp1 = P.muH(2:end) + P.muD(2:end).*DHr(2:end,n+1)./PHrp1(2:end);
    temp = (P.vc(2:end)+P.vs(2:end))./PHrp1(2:end);
    Cacr(2:end,n+1) = (Cacr(1:end-1,n)+dt*Brnp1)./(1+dt*(1/P.dac+Drnp1+temp));
    Cmr(2:end,n+1) = Cmr(1:end-1,n)./(1+dt*(1/P.dm+Drnp1+temp));
    Cvr(2:end,n+1) = 0;
    Ctotr(:,n+1) = P.c1*Cacr(:,n+1)+P.c2*Cmr(:,n+1)+P.c3*Cvr(:,n+1); 
    % vacc group
    Bvnp1 = f(lamHp1(2:end)).*(P.cS*SHv(2:end,n+1) + P.cE*EHv(2:end,n+1) + P.cA*AHv(2:end,n+1) + P.cD*DHv(2:end,n+1) + P.cU*UHv(2:end,n+1));
    Dvnp1 = P.muH(2:end) + P.muD(2:end).*DHv(2:end,n+1)./PHvp1(2:end)+vb_temp./PHvp1(2:end);
    temp1 = P.vs(2:end).*Cacr(2:end,n+1)./PHrp1(2:end);
    Cacv(2:end,n+1) = (Cacv(1:end-1,n)+dt*(Bvnp1+temp1))./(1+dt*(1/P.dac+Dvnp1));
    temp2 = P.vs(2:end).*Cmr(2:end,n+1)./PHrp1(2:end);
    Cmv(2:end,n+1) = (Cmv(1:end-1,n)+temp2*dt)./(1+dt*(1/P.dm+Dvnp1));
    Cvv(2:end,n+1) = 0;
    Cmv(PHvp1==0,n+1) = 0; Cacv(PHvp1==0,n+1) = 0; Cvv(PHvp1==0,n+1) = 0;
    Ctotv(:,n+1) = P.c1*Cacv(:,n+1)+P.c2*Cmv(:,n+1)+P.c3*Cvv(:,n+1); 
    % boost
    Bbnp1 = f(lamHp1(2:end)).*(P.cS*SHb(2:end,n+1) + P.cE*EHb(2:end,n+1) + P.cA*AHb(2:end,n+1) + P.cD*DHb(2:end,n+1) + P.cU*UHb(2:end,n+1));
    temp1 = vb_temp.*Cacv(2:end,n+1)./PHvp1(2:end);
    Dbnp1 = P.muH(2:end) + P.muD(2:end).*DHb(2:end,n+1)./PHbp1(2:end);
    Cacb(2:end,n+1) = (Cacb(1:end-1,n)+dt*(Bbnp1+temp1))./(1+dt*(1/P.dac+Dbnp1));
    temp2 = vb_temp.*Cmv(2:end,n+1)./PHvp1(2:end);
    Cmb(2:end,n+1) = (Cmb(1:end-1,n)+temp2*dt)./(1+dt*(1/P.dm+Dbnp1));
    Cvb(2:end,n+1) = 0;
    Cmb(PHbp1==0,n+1) = 0; Cacb(PHbp1==0,n+1) = 0; Cvb(PHbp1==0,n+1) = 0;
    Ctotb(:,n+1) = P.c1*Cacb(:,n+1)+P.c2*Cmb(:,n+1)+P.c3*Cvb(:,n+1);
    % control
    Bcnp1 = f(lamHp1(2:end)).*(P.cS*SHc(2:end,n+1) + P.cE*EHc(2:end,n+1) + P.cA*AHc(2:end,n+1) + P.cD*DHc(2:end,n+1));
    Dcnp1 = P.muH(2:end) + P.muD(2:end).*DHc(2:end,n+1)./PHcp1(2:end);
    temp1 = P.vc(2:end).*Cacr(2:end,n+1)./PHrp1(2:end);
    Cacc(2:end,n+1) = (Cacc(1:end-1,n)+dt*(Bcnp1+temp1))./(1+dt*(1/P.dac+Dcnp1));
    temp2 = P.vc(2:end).*Cmr(2:end,n+1)./PHrp1(2:end);
    Cmc(2:end,n+1) = (Cmc(1:end-1,n)+temp2*dt)./(1+dt*(1/P.dm+Dcnp1));
    Cvc(2:end,n+1) = Cvc(1:end-1,n)./(1+dt*(1/P.dv+Dcnp1));
    Cmc(PHcp1==0,n+1) = 0; Cacc(PHcp1==0,n+1) = 0; Cvc(PHcp1==0,n+1) = 0;
    Ctotc(:,n+1) = P.c1*Cacc(:,n+1)+P.c2*Cmc(:,n+1)+P.c3*Cvc(:,n+1);     
    
    % update progression probability based on immunity Ctot
    P.phir = sigmoid_prob(Ctotr(:,n+1)./PHrp1, 'phi'); % prob. of DH -> RH
    P.rhor = sigmoid_prob(Ctotr(:,n+1)./PHrp1, 'rho'); % prob. of severely infected, EH -> DH
    P.psir = sigmoid_prob(Ctotr(:,n+1)./PHrp1, 'psi'); % prob. AH -> DH
    
    P.phiv = sigmoid_prob(Ctotv(:,n+1)./PHvp1, 'phi'); % prob. of DH -> RH
    P.rhov = sigmoid_prob(Ctotv(:,n+1)./PHvp1, 'rho'); % prob. of severely infected, EH -> DH
    P.psiv = sigmoid_prob(Ctotv(:,n+1)./PHvp1, 'psi'); % prob. AH -> DH
    P.phiv(PHvp1==0) = 0; P.rhov(PHvp1==0) = 0; P.psiv(PHvp1==0) = 0;
    
    P.phib = sigmoid_prob(Ctotb(:,n+1)./PHbp1, 'phi'); % prob. of DH -> RH
    P.rhob = sigmoid_prob(Ctotb(:,n+1)./PHbp1, 'rho'); % prob. of severely infected, EH -> DH
    P.psib = sigmoid_prob(Ctotb(:,n+1)./PHbp1, 'psi'); % prob. AH -> DH
    P.phib(PHbp1==0) = 0; P.rhob(PHbp1==0) = 0; P.psib(PHbp1==0) = 0;
    
    P.phic = sigmoid_prob(Ctotc(:,n+1)./PHcp1, 'phi'); % prob. of DH -> RH
    P.rhoc = sigmoid_prob(Ctotc(:,n+1)./PHcp1, 'rho'); % prob. of severely infected, EH -> DH
    P.psic = sigmoid_prob(Ctotc(:,n+1)./PHcp1, 'psi'); % prob. AH -> DH
    P.phic(PHcp1==0) = 0; P.rhoc(PHcp1==0) = 0; P.psic(PHcp1==0) = 0;

end

end
