%% sigmoidal function for converting immunity to probabilties
% evaluate the sigmoid function - return discrete values
function y = sigmoid_prob(x, lprob)
global P

switch lprob
    case 'phi'
        y = sigmoid_inc(P.phif0, P.phif1, x, P.phis2, P.phir2);
        % y = P.phi0*ones(size(x));
    case 'rho'
        y = sigmoid_dec(P.rhof0, P.rhof1, x, P.rhos2, P.rhor2);
        % y = P.rho0*ones(size(x));
    case 'psi'
        y = sigmoid_dec(P.psif0, P.psif1, x, P.psis2, P.psir2);
        % y = P.psi0*ones(size(x));
    otherwise
        error('not defined probability parameter')

end

end

% NB set f1 = f0 to get a constant function (equal to f0 everywhere)
function y = sigmoid_inc(f_0, f_1, x, s_2, r_2)
y = f_0 + (f_1-f_0)./(1 + exp(-(x-s_2)/r_2));
end

function y = sigmoid_dec(f_0, f_1, x, s_2, r_2)
y = f_1 + (f_0-f_1)./(1 + exp(-(x-s_2)/r_2));
end
