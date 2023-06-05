function [bH,bM] = biting_rate(PH,NM)
global P
% happy mosquitoes:
% bM = P.bm;
% bH = P.bm*NM/NH;

% compromise model from Chitnis et al. (2006)
%b_tot = P.bm*P.bh*NH.*NM./(P.bm*NM+P.bh*NH);
%bH = b_tot./NH; % bites per human

% age-dependent biting rate (surface area factor)
temp = trapz(P.zeta.*PH,1)*P.da;
b_tot = P.bm*P.bh.*NM.*temp./(P.bm*NM+P.bh*temp);
bH = (b_tot./temp).*P.zeta; % bites per human

bM = b_tot./NM; % bites per mosquito

end