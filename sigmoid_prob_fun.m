%% sigmoidal function for converting immunity to probabilties
% return function handle
function fun = sigmoid_prob_fun(lprob)
global P

sigmoid_inc = @(f_0, f_1, x, s_2, r_2) f_0 + (f_1-f_0)./(1 + exp(-(x-s_2)/r_2));
sigmoid_dec = @(f_0, f_1, x, s_2, r_2) f_1 + (f_0-f_1)./(1 + exp(-(x-s_2)/r_2));

% NB set f1 = f0 to get a constant function (equal to f0 everywhere)

switch lprob
    case 'phi'
        fun = @(x) sigmoid_inc(P.phif0, P.phif1, x, P.phis2, P.phir2);
        %         fun = @(x) P.phi0*ones(size(x));
    case 'rho'
        fun = @(x) sigmoid_dec(P.rhof0, P.rhof1, x, P.rhos2, P.rhor2);
        %         fun = @(x)  P.rho0*ones(size(x));
    case 'psi'
        fun = @(x) sigmoid_dec(P.psif0, P.psif1, x, P.psis2, P.psir2);
         %         fun = @(x) P.psi0*ones(size(x));
    otherwise
        error('not defined probability parameter')
end

end

% f_0  value min
% f_1  value max
% t_2  shift
% s_2  sigmoid steepness, smaller is steeper


