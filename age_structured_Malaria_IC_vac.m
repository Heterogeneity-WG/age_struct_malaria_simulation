function  [SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_IC_vac(state)
% define the initial condition of the simulation

global P

na = P.na;
da = P.da;
NH = P.NN;
gM = P.gM_fun(0);
NM = gM/P.muM;
switch state
    case 'init' %
        SH = 0.4*P.PH_stable*NH; %0.9*NH/na/da*ones(na,1); % cell averages; %
        EH = 0.1*P.PH_stable*NH; %0.1*NH/na/da*ones(na,1); %
        DH = 0.25*P.PH_stable*NH;
        AH = 0.25*P.PH_stable*NH;
        VH = 0*NH/na/da*ones(na,1);
        UH = 0*NH/na/da*ones(na,1);
        % diagnostic compartment
        MH = 0*NH/na/da*ones(na,1);
        
        % for mosquitoes
        SM = NM; EM = 0; IM = 0;
        
        %         NH = trapz(SH+EH+DH+AH+VH+UH)*da;
        %         temp = P.lMsystem;
        %         P.lMsystem = 'ss';
        %         [SM,EM,IM] = mosquito_ODE(SM, EM, IM, DH, AH, NH, NH, NM);
        %         P.lMsystem = temp;
        
        Cm = 0*ones(na,1);
        Cac = 0*ones(na,1);
        Cv = 0*ones(na,1);
        Ctot = P.c1*Cac+P.c2*Cm+P.c3*Cv;
        
    case 'EE' % start from EE
        [SH, EH, DH, AH, VH, UH, SM, EM, IM, Cac, Cm, Cv, Ctot, MH] = steady_state_vac('EE','numerical');
    case 'EE_reset' % start from EE
        [SH, EH, DH, AH, VH, UH, SM, EM, IM, Cac, Cm, Cv, ~, ~] = steady_state_vac('EE','numerical');
        PH = SH + EH + DH + AH + VH + UH;
        SH0 = (SH./PH).*P.PH_stable*NH; %0.9*NH/na/da*ones(na,1); % cell averages; %
        EH0 = (EH./PH).*P.PH_stable*NH; %0.1*NH/na/da*ones(na,1); %
        DH0 = (DH./PH).*P.PH_stable*NH;
        AH0 = (AH./PH).*P.PH_stable*NH;
        VH0 = (VH./PH).*P.PH_stable*NH;
        UH0 = (UH./PH).*P.PH_stable*NH;
        % diagnostic compartment
        MH0 = 0*NH/na/da*ones(na,1);
        
        % for mosquitoes
        SM0 = SM; EM0 = EM; IM0 = IM;
        
        Cm0 = (Cm./PH).*P.PH_stable*NH;
        Cac0 = (Cac./PH).*P.PH_stable*NH;
        Cv0 = (Cv./PH).*P.PH_stable*NH;
        Ctot0 = P.c1*Cac0 + P.c2*Cm0 + P.c3*Cv0;
        
        % tfinal= 0.5*365; % for non-seasonal case 
        tfinal= 2*365; % for seasonal case
        [~,SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(P.da,P.na, 0, tfinal,...
            SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
        SH = SH(:,end); EH = EH(:,end); DH = DH(:,end); AH = AH(:,end);  UH = UH(:,end);  VH = VH(:,end); MH = 0*MH(:,end);
        SM = SM(end); EM = EM(end);  IM = IM(end);
        Cac = Cac(:,end); Cm = Cm(:,end); Cv = Cv(:,end); Ctot = Ctot(:,end);
    otherwise
        keyboard
end
PH = SH+EH+DH+AH+VH+UH;
P.phi = sigmoid_prob(Ctot./PH, 'phi'); % prob. of DH -> RH
P.rho = sigmoid_prob(Ctot./PH, 'rho'); % prob. of severely infected, EH -> DH
P.psi = sigmoid_prob(Ctot./PH, 'psi'); % prob. AH -> DH

end