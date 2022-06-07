function  [SH, EH, DH, AH, SM, EM, IM, Cm, Cac, Ctot] = age_structured_Malaria_IC(state)
% define the initial condition of the simulation

global P

na = P.na;
da = P.da;
NM = P.gM/P.muM;
NH = 1;

switch state
    case 'init' %        
        SH = 0.99*P.PH_stable*NH; 
        EH = 0.01*P.PH_stable*NH;
        DH = 0.0*P.PH_stable*NH;
        AH = 0.0*P.PH_stable*NH; 
%         SH = NH/na/da*ones(na,1); % 
%         EH = 0*NH/na/da*ones(na,1); % 
%         DH = 0*NH/na/da*ones(na,1); % 
%         AH = 0*NH/na/da*ones(na,1); % 
        
        % for mosquitoes - assume at equilibrium
        NH = trapz(SH+EH+DH+AH)*da;
        [SM,EM,IM] = mosquito_ODE(DH,AH,NH,NM);
        
        Cm = 0*ones(na,1);
        Cac = 0*ones(na,1);
        Ctot = P.c1*Cac+P.c2*Cm;
        
    case 'EE' % start from EE    
        [SH, EH, DH, AH, Cac, Cm, Ctot] = steady_state('EE','numerical');
        NH = trapz(SH+EH+DH+AH)*da;
        [SM,EM,IM] = mosquito_ODE(DH,AH,NH,NM);        

    otherwise
        keyboard
end

PH = SH+EH+DH+AH;
P.phi = sigmoid_prob(Ctot./PH, 'phi'); % prob. of DH -> RH
P.rho = sigmoid_prob(Ctot./PH, 'rho'); % prob. of severely infected, EH -> DH
P.psi = sigmoid_prob(Ctot./PH, 'psi'); % prob. AH -> DH

end