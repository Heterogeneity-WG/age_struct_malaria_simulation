function [SMp1,EMp1,IMp1] = mosquito_ODE(SM, EM, IM, DH, AH, NH, NHp1, NM)
% evolve mosquito quantities

global P
dt = P.dt;

[~,bM] = biting_rate(NH,NM); % at t=n
lamM = FOI_M(bM,DH,AH,NH); % at t=n

if strcmp(P.lMHfix,'on') % if the assumption is on: mosquito/human ratio is fixed
    P.gM = P.MHm*NHp1*P.muM; % update gM at t = n+1;
elseif ~strcmp(P.lMHfix,'off')
    error('undefined label for lMHfix')  
end

switch P.lMsystem
    case 'ss' % keep at the steady state 
        SMp1 = P.gM/(lamM+P.muM);
        IMp1 = (P.gM/P.muM)*(P.sigma/(P.sigma+P.muM))*(lamM/(lamM+P.muM));
        EMp1 = (P.muM/P.sigma)*IMp1;
    case 'full' % evolve using the full system
        SMp1 = (SM+dt*P.gM)/(1+dt*(lamM+P.muM));
        EMp1 = (EM+dt*lamM*SMp1)/(1+dt*(P.sigma+P.muM));
        IMp1 = (IM+dt*P.sigma*EMp1)/(1+dt*P.muM);
    otherwise
        error('undefined label for lMsystem')
end

end