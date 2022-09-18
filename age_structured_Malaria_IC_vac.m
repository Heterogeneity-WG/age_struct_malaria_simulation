function  [SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_IC_vac(state)
% define the initial condition of the simulation

global P

na = P.na;
da = P.da;
NH = P.NN;
NM = P.gM/P.muM;
switch state
    case 'init' %
        SH = 0.97*P.PH_stable*NH; %0.9*NH/na/da*ones(na,1); % cell averages; %
        EH = 0.01*P.PH_stable*NH; %0.1*NH/na/da*ones(na,1); % 
        DH = 0.01*P.PH_stable*NH;
        AH = 0.01*P.PH_stable*NH;
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
        
    otherwise
        keyboard
end
PH = SH+EH+DH+AH+VH+UH;
P.phi = sigmoid_prob(Ctot./PH, 'phi'); % prob. of DH -> RH
P.rho = sigmoid_prob(Ctot./PH, 'rho'); % prob. of severely infected, EH -> DH
P.psi = sigmoid_prob(Ctot./PH, 'psi'); % prob. AH -> DH

end