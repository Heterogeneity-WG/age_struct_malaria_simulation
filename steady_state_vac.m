function [S,E,D,A,V,U,SM, EM, IM,Cac,Cm,Cv,Ctot, M] = steady_state_vac(lstate,lreturn)
% lreturn = 'handle' return function handles for 'DFE'
% lreturn = 'numerical' return numerical values
% all variables are pop size
global P

gH = P.gH_fun; % feritlity function handle
muH_int = P.muH_int_fun; % death function handle M
switch lstate
    case 'DFE'
        if strcmp(lreturn,'handle')
            S = @(alpha) P.theta_fun(alpha).*ones(size(alpha)).*P.PH_stable_fun(alpha);
            E = @(alpha) 0*ones(size(alpha)).*P.PH_stable_fun(alpha);
            D = @(alpha) 0*ones(size(alpha)).*P.PH_stable_fun(alpha);
            A = @(alpha) 0*ones(size(alpha)).*P.PH_stable_fun(alpha);
            V = @(alpha) (1-P.theta_fun(alpha)).*ones(size(alpha)).*P.PH_stable_fun(alpha);
            U = @(alpha) (1-P.theta_fun(alpha)).*ones(size(alpha)).*P.PH_stable_fun(alpha);
            % Cac
            Cac_prop = @(alpha) 0*ones(size(alpha));
            % Cv exact expression obtained based on simple v(alpha)
            Cv_prop = @(alpha) P.cV*exp(-alpha/P.dv).*(P.vb0+P.vb0*P.dv.*P.theta_fun(alpha).*(exp(alpha/P.dv)-1));
            % Cm
            Cm0 = P.m*integral(@(alpha) gH(alpha).*exp(-muH_int(alpha)).*P.c3.*Cv_prop(alpha), 0, P.age_max);
            Cm_prop = @(alpha) Cm0.*exp(-alpha./P.dm);
            % CH
            Cac = @(alpha) Cac_prop(alpha).*P.PH_stable_fun(alpha);
            Cm = @(alpha) Cm_prop(alpha).*P.PH_stable_fun(alpha);
            Cv = @(alpha) Cv_prop(alpha).*P.PH_stable_fun(alpha);
            Ctot = @(alpha) P.c1*Cac(alpha)+P.c2*Cm(alpha)+P.c3*Cv(alpha);
        elseif strcmp(lreturn,'numerical')
            a = P.a;
            S = 0*ones(size(a)).*P.PH_stable; V=S; U=S;
            E = 0*ones(size(a)).*P.PH_stable;
            D = 0*ones(size(a)).*P.PH_stable;
            A = 0*ones(size(a)).*P.PH_stable;
            keyboard
            for i = 1:length(a)
                S(i) = trapz(P.a(1:i), P.v(1:i));
                int1 = trapz(a(1:i),exp(P.w*a(1:i)).*P.etas.*(1-P.z).*P.v(1:i).*S(1:i));
                V(i) = 0;
            end
            S = S.*P.PH_stable;
            % Cac
            Cac_prop = 0*ones(size(a));
            % Cv exact expression obtained based on simple v(alpha)
            Cv_prop = P.cV*exp(-a./P.dv).*(P.vb0+P.vb0*P.dv.*P.theta.*(exp(a./P.dv)-1));
            % Cm
            Cm0 = P.m*trapz(P.gH.*exp(-P.muH_int).*P.c3.*Cv_prop)*P.da;
            Cm = Cm0.*exp(-a./P.dm).*P.PH_stable;
            Cac = Cac_prop.*P.PH_stable;
            Cv = Cv_prop.*P.PH_stable;
            Ctot = P.c1*Cac+P.c2*Cm+P.c3*Cv;
        end
    case 'EE'
        if strcmp(lreturn,'numerical')
            [SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0] = age_structured_Malaria_IC_vac('init');
            tfinal= 30*365;  % run for a long time; numerical EE
            [~,SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = age_structured_Malaria_vac(P.da,P.na,0,tfinal,...
                SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0, MH0);
            S = SH(:,end); E = EH(:,end); D = DH(:,end); A = AH(:,end);  U = UH(:,end);  V = VH(:,end); M = MH(:,end); 
            SM = SM(end); EM = EM(end);  IM = IM(end); 
            Cac = Cac(:,end); Cm = Cm(:,end); Cv = Cv(:,end); Ctot = Ctot(:,end);
        elseif strcmp(lreturn,'fsolve')
            % use numerical simulation for an initial guess
            dt = 20; tfinal= 15*365;  % run for a few years to get closer to EE
            da = dt; a = (0:da:P.age_max)'; na = length(a);
            [SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0] = age_structured_Malaria_IC_vac('init');
            [SH, EH, DH, AH, VH, UH, ~, ~, ~, ~, Cac, Cv, ~] = age_structured_Malaria_vac(da,na,tfinal, SH0, EH0, DH0, AH0, VH0, UH0, SM0, EM0, IM0, Cm0, Cac0, Cv0, Ctot0);
            %% run solver on a coarser grid to speed up
            da_fine = P.da; da_coarse = 80; P.da = da_coarse;
            a_fine = P.a; a_coarse = (0:da_coarse:P.age_max)'; P.a = a_coarse;
            na_fine = P.na; na_coarse = length(a_coarse); P.na = na_coarse;
            Malaria_parameters_transform;
            SH0 = interp1(a_fine,SH(:,end),a_coarse);
            EH0 = interp1(a_fine,EH(:,end),a_coarse);
            DH0 = interp1(a_fine,DH(:,end),a_coarse);
            AH0 = interp1(a_fine,AH(:,end),a_coarse);
            VH0 = interp1(a_fine,VH(:,end),a_coarse);
            UH0 = interp1(a_fine,UH(:,end),a_coarse);
            Cac0 = interp1(a_fine,Cac(:,end),a_coarse);
            Cv0 = interp1(a_fine,Cv(:,end),a_coarse);
            PH0 = SH0 + EH0 + AH0 + DH0 + VH0 + UH0;
            x0 = [SH0./PH0; EH0./PH0; DH0./PH0; AH0./PH0; VH0./PH0; UH0./PH0; Cac0./PH0; Cv0./PH0];
            options = optimoptions('fsolve','Display','none','OptimalityTolerance', 1e-25);
            F_prop = @(x) human_model_der_prop_vac(x);
            [xsol,err,~,~,~] = fsolve(F_prop,x0,options);
            if max(max(abs(err)))>10^-5
                disp('not converged')
                keyboard
            end
            x_EE = reshape(xsol,[P.na,7]);
            S_coarse = x_EE(:,1);
            E_coarse = x_EE(:,2);
            D_coarse = x_EE(:,3);
            A_coarse = x_EE(:,4);
            V_coarse = x_EE(:,5);
            Cac_coarse = x_EE(:,6);
            Cv_coarse = x_EE(:,7);
            %% return the values on fine grid with interpolation
            % recover numerical config
            P.da = da_fine; P.a = a_fine; P.na = na_fine;
            Malaria_parameters_transform;
            S = interp1(a_coarse,S_coarse,a_fine).*P.PH_stable;
            E = interp1(a_coarse,E_coarse,a_fine).*P.PH_stable;
            D = interp1(a_coarse,D_coarse,a_fine).*P.PH_stable;
            A = interp1(a_coarse,A_coarse,a_fine).*P.PH_stable;
            V = interp1(a_coarse,V_coarse,a_fine).*P.PH_stable;
            Cac = interp1(a_coarse,Cac_coarse,a_fine).*P.PH_stable; % pooled
            Cv = interp1(a_coarse,Cv_coarse,a_fine).*P.PH_stable; % pooled
            Cm0 = P.m*trapz(P.gH.*(P.c1*Cac+P.c3*Cv))*P.da/P.PH_stable(1);
            Cm = Cm0.*exp(-P.a./P.dm).*P.PH_stable; % pooled
            Ctot = P.c1*Cac+P.c2*Cm+P.c3*Cv;
        end
        
    otherwise
        error('undefined steady state label...')
end
end
