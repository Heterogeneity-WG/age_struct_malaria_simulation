function [t, SH, EH, DH, AH, VH, UH, Cm, Cac, Ctot, SH2, EH2, DH2, AH2, VH2, UH2, ...
    Cm2, Cac2, Ctot2, SM, EM, IM, vb_temp] = ...
    age_structured_Malaria_booster(da, na, t0, tfinal, ...
    SH0, EH0, DH0, AH0, VH0, UH0, Cm0, Cac0, Ctot0, ...
    SH20, EH20, DH20, AH20, VH20, UH20, Cm20, Cac20, Ctot20, SM0, EM0, IM0)
% SH - never receive any vac; VH - primary dose
% SH2 - vaccinated; VH2 - booster dose
global P

dt = da;
t = (t0:dt:tfinal)';
nt = length(t);

% allocation
SM = NaN(1,nt); EM = NaN(1,nt); IM = NaN(1,nt);
SM(1) = SM0; EM(1) = EM0; IM(1) = IM0;

NH = NaN(1,nt); NM = NaN(nt);
NM(1) = SM(1)+EM(1)+IM(1);

% SH - never receive any vac; VH - primary dose
SH = NaN(na,nt); EH = NaN(na,nt); DH = NaN(na,nt); AH = NaN(na,nt); VH = NaN(na,nt); UH = NaN(na,nt);
Cm = NaN(na,nt); Cac = NaN(na,nt); Cv = NaN(na,nt); Ctot = NaN(na,nt);
SH(:,1) = SH0; EH(:,1) = EH0; DH(:,1) = DH0; AH(:,1) = AH0; VH(:,1) = VH0; UH(:,1) = UH0; 
Cm(:,1) = Cm0; Cac(:,1) = Cac0; Cv(:,1) = 0; Ctot(:,1) = Ctot0;

% SH2 - vaccinated; VH2 - booster dose
SH2 = NaN(na,nt); EH2 = NaN(na,nt); DH2 = NaN(na,nt); AH2 = NaN(na,nt); VH2 = NaN(na,nt); UH2 = NaN(na,nt);
Cm2 = NaN(na,nt); Cac2 = NaN(na,nt); Cv2 = NaN(na,nt); Ctot2 = NaN(na,nt);
vb_temp = NaN(na,nt);
SH2(:,1) = SH20; EH2(:,1) = EH20; DH2(:,1) = DH20; AH2(:,1) = AH20; VH2(:,1) = VH20; UH2(:,1) = UH20; 
Cm2(:,1) = Cm20; Cac2(:,1) = Cac20; Cv2(:,1) = 0; Ctot2(:,1) = Ctot20;
vb_temp(:,1) = zeros(size(SH20));

PH = SH0+EH0+DH0+AH0+VH0+UH0;
PH2 = SH20+EH20+DH20+AH20+VH20+UH20;

% update progression probability based on immunity Ctot
P.phi = sigmoid_prob(Ctot(:,1)./PH, 'phi'); % prob. of DH -> RH
P.rho = sigmoid_prob(Ctot(:,1)./PH, 'rho'); % prob. of severely infected, EH -> DH
P.psi = sigmoid_prob(Ctot(:,1)./PH, 'psi'); % prob. AH -> DH
if PH == 0; keyboard; end
P.phi(PH==0) = 0; P.rho(PH==0) = 0; P.psi(PH==0) = 0;

P.phi2 = sigmoid_prob(Ctot2(:,1)./PH2, 'phi'); % prob. of DH -> RH
P.rho2 = sigmoid_prob(Ctot2(:,1)./PH2, 'rho'); % prob. of severely infected, EH -> DH
P.psi2 = sigmoid_prob(Ctot2(:,1)./PH2, 'psi'); % prob. AH -> DH
P.phi2(PH2==0) = 0; P.rho2(PH2==0) = 0; P.psi2(PH2==0) = 0;

%% time evolution
for n = 1:nt-1
    PH = SH(:,n)+EH(:,n)+DH(:,n)+AH(:,n)+VH(:,n)+UH(:,n); 
    PH2 = SH2(:,n)+EH2(:,n)+DH2(:,n)+AH2(:,n)+VH2(:,n)+UH2(:,n); 

    PHtot = PH+PH2;
    NH(n) = trapz(PHtot)*da;
    NM(n) = SM(1,n)+EM(1,n)+IM(1,n);
    [bH,~] = biting_rate(PHtot,NM(n));
    lamH = FOI_H(bH,IM(1,n),NM(n));
    
    % human birth terms
    SH(1,n+1) = trapz(P.gH.*PHtot)*da;
    EH(1,n+1) = 0;
    DH(1,n+1) = 0;
    AH(1,n+1) = 0;
    VH(1,n+1) = 0;
    UH(1,n+1) = 0;

    SH2(1,n+1) = 0;
    EH2(1,n+1) = 0;
    DH2(1,n+1) = 0;
    AH2(1,n+1) = 0;
    VH2(1,n+1) = 0;
    UH2(1,n+1) = 0;
    vb_temp(1,n+1) = 0;

    %% primary group
    SH(2:end,n+1) = (SH(1:end-1,n)+dt*(P.phi(1:end-1)*P.rD.*DH(1:end-1,n)+P.rA*AH(1:end-1,n)-P.v(1:end-1)))...
        ./(1+(lamH(1:end-1)+P.muH(2:end))*dt);
    if min(SH(2:end,n+1))<0
        keyboard
    end
    VH(2:end,n+1) = (VH(1:end-1,n)+dt*P.etas*(1-P.z)*P.v(1:end-1))...
        ./(1+(P.muH(2:end)+P.w)*dt);
    temp1 = (P.z*P.v(1:end-1)+(1-P.etas)*(1-P.z)*P.v(1:end-1))+P.w*VH(2:end,n+1);
    UH(2:end,n+1) = (UH(1:end-1,n)+dt*temp1)./(1+(lamH(1:end-1)+P.muH(2:end))*dt);
    EH(2:end,n+1) = (EH(1:end-1,n)+dt*lamH(1:end-1).*(SH(2:end,n+1)))...
        ./(1+(P.h+P.muH(2:end))*dt);
    temp2 = (1-P.rho(1:end-1))*P.h.*EH(2:end,n+1)+(1-P.phi(1:end-1)).*P.rD.*DH(1:end-1,n);
    AH(2:end,n+1) = ((1-dt*P.rA)*AH(1:end-1,n)+dt*temp2)...
        ./(1+dt*(P.psi(1:end-1).*lamH(1:end-1)+P.muH(2:end)));
    temp3 = P.rho(1:end-1)*P.h.*EH(2:end,n+1)+P.psi(1:end-1).*lamH(1:end-1).*AH(2:end,n+1);
    DH(2:end,n+1) = ((1-dt*P.rD)*DH(1:end-1,n)+dt*temp3)...
        ./(1+dt*(P.muH(2:end)+P.muD(2:end)));
    % diagnostic, disease-induced mortality counts, accumulative (ODE)
    % MH(:,n+1) = MH(:,n)+dt*P.muD.*DH(:,n);
    
    %% boosted group
    % the maximum vb to avoid SHV<0
    vb_temp(2:end,n+1) = max(min(P.vb(2:end), SH2(1:end-1,n)/dt+(P.phi2(1:end-1)*P.rD.*DH2(1:end-1,n)+P.rA*AH2(1:end-1,n))-10^-5),0);
    SH2(2:end,n+1) = (SH2(1:end-1,n)+dt*(P.phi2(1:end-1)*P.rD.*DH2(1:end-1,n)+P.rA*AH2(1:end-1,n)-vb_temp(2:end,n+1)))...
        ./(1+(lamH(1:end-1)+P.muH(2:end))*dt); 
    if min(SH2(2:end,n+1))<0
        keyboard
    end
    VH2(2:end,n+1) = (VH2(1:end-1,n)+dt*P.etab*vb_temp(2:end,n+1))./(1+(P.muH(2:end)+P.wb)*dt);
    temp1 = (1-P.etab)*vb_temp(2:end,n+1)+P.wb*VH2(2:end,n+1);
    UH2(2:end,n+1) = (UH2(1:end-1,n)+dt*temp1)./(1+(lamH(1:end-1)+P.muH(2:end))*dt);
    
    
    
    EH2(2:end,n+1) = (EH2(1:end-1,n)+dt*lamH(1:end-1).*(SH2(2:end,n+1)+UH2(2:end,n+1)+UH(2:end,n+1)))...
        ./(1+(P.h+P.muH(2:end))*dt);
    temp2 = (1-P.rho2(1:end-1))*P.h.*EH2(2:end,n+1)+(1-P.phi2(1:end-1)).*P.rD.*DH2(1:end-1,n);
    AH2(2:end,n+1) = ((1-dt*P.rA)*AH2(1:end-1,n)+dt*temp2)...
        ./(1+dt*(P.psi2(1:end-1).*lamH(1:end-1)+P.muH(2:end)));
    temp3 = P.rho2(1:end-1)*P.h.*EH2(2:end,n+1)+P.psi2(1:end-1).*lamH(1:end-1).*AH2(2:end,n+1);
    
    DH2(2:end,n+1) = ((1-dt*P.rD)*DH2(1:end-1,n)+dt*temp3)...
        ./(1+dt*(P.muH(2:end)+P.muD(2:end)));   

    %%
    PHp1 = SH(:,n+1)+EH(:,n+1)+DH(:,n+1)+AH(:,n+1)+VH(:,n+1)+UH(:,n+1); 
    PH2p1 = SH2(:,n+1)+EH2(:,n+1)+DH2(:,n+1)+AH2(:,n+1)+VH2(:,n+1)+UH2(:,n+1); 
    PHtotp1 = PHp1+PH2p1;    
    NHtotp1 = trapz(PHtotp1)*da; % total human population at t=n+1;
    
    % mosquito time evolution
    P.gM = P.gM_fun(t(n+1)); % incorporate seasonlity
    [SM(1,n+1),EM(1,n+1),IM(1,n+1)] = mosquito_ODE(SM(1,n), EM(1,n), IM(1,n), DH(:,n)+DH2(:,n), AH(:,n)+AH2(:,n), PHtotp1, NHtotp1, NM(n));
    NM(n+1) = SM(1,n+1)+EM(1,n+1)+IM(1,n+1);
    
    % immunity gained at age = 0 
    Cacn = Cac(:,n)+Cac2(:,n);
    Cvn = Cv(:,n)+Cv2(:,n);
    Cm(1,n+1) = P.m*trapz(P.gH.*(P.c1*Cacn+P.c3*Cvn))*da;
    Cac(1,n+1) = 0;
    Cv(1,n+1) = 0;
    
    Cm2(1,n+1) = 0; Cac2(1,n+1) = 0; Cv2(1,n+1) = 0;
        
    % acquired immunity - use Qn+1
    [bHp1,~] = biting_rate(PHtotp1,NM(n+1));
    lamHp1 = FOI_H(bHp1,IM(1,n+1),NM(n+1));
    
    % pooled immunity
    % primary group
    Bnp1 = f(lamHp1(2:end)).*(P.cS*SH(2:end,n+1) + P.cE*EH(2:end,n+1) + P.cA*AH(2:end,n+1) + P.cD*DH(2:end,n+1) + P.cU*UH(2:end,n+1));
    Dnp1 = P.muH(2:end) + P.muD(2:end).*DH(2:end,n+1)./PHp1(2:end)+lamHp1(2:end).*UH(2:end,n+1)./PHp1(2:end);
    Cac(2:end,n+1) = (Cac(1:end-1,n)+dt*(Bnp1))./(1+dt*(1/P.dac+Dnp1));
    Cm(2:end,n+1) = (Cm(1:end-1,n))./(1+dt*(1/P.dm+Dnp1));
    Cv(2:end,n+1) = 0;
    Cm(PHp1==0,n+1) = 0; Cac(PHp1==0,n+1) = 0; Cv(PHp1==0,n+1) = 0;
    Ctot(:,n+1) = P.c1*Cac(:,n+1)+P.c2*Cm(:,n+1)+P.c3*Cv(:,n+1); 
    % boost group
    B2np1 = f(lamHp1(2:end)).*(P.cS*SH2(2:end,n+1) + P.cE*EH2(2:end,n+1) + P.cA*AH2(2:end,n+1) + P.cD*DH2(2:end,n+1) + P.cU*UH2(2:end,n+1));
    temp1 = lamHp1(2:end).*UH(2:end,n+1).*Cac(2:end,n+1)./PHp1(2:end);
    D2np1 = P.muH(2:end) + P.muD(2:end).*DH2(2:end,n+1)./PH2p1(2:end);
    Cac2(2:end,n+1) = (Cac2(1:end-1,n)+dt*(B2np1+temp1))./(1+dt*(1/P.dac+D2np1));
    temp2 = lamHp1(2:end).*UH(2:end,n+1).*Cm(2:end,n+1)./PHp1(2:end);
    Cm2(2:end,n+1) = (Cm2(1:end-1,n)+temp2*dt)./(1+dt*(1/P.dm+D2np1));
    Cv2(2:end,n+1) = 0;
    Cm2(PH2p1==0,n+1) = 0; Cac2(PH2p1==0,n+1) = 0; Cv2(PH2p1==0,n+1) = 0;
    Ctot2(:,n+1) = P.c1*Cac2(:,n+1)+P.c2*Cm2(:,n+1)+P.c3*Cv2(:,n+1);
        
    % update progression probability based on immunity Ctot
    P.phi = sigmoid_prob(Ctot(:,n+1)./PHp1, 'phi'); % prob. of DH -> RH
    P.rho = sigmoid_prob(Ctot(:,n+1)./PHp1, 'rho'); % prob. of severely infected, EH -> DH
    P.psi = sigmoid_prob(Ctot(:,n+1)./PHp1, 'psi'); % prob. AH -> DH
    P.phi(PHp1==0) = 0; P.rho(PHp1==0) = 0; P.psi(PHp1==0) = 0;
    
    P.phi2 = sigmoid_prob(Ctot2(:,n+1)./PH2p1, 'phi'); % prob. of DH -> RH
    P.rho2 = sigmoid_prob(Ctot2(:,n+1)./PH2p1, 'rho'); % prob. of severely infected, EH -> DH
    P.psi2 = sigmoid_prob(Ctot2(:,n+1)./PH2p1, 'psi'); % prob. AH -> DH
    P.phi2(PH2p1==0) = 0; P.rho2(PH2p1==0) = 0; P.psi2(PH2p1==0) = 0;
    
end

end
