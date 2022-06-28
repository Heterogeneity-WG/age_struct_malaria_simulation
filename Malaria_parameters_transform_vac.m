function Malaria_parameters_transform_vac

global P

%% vaccination functions (baseline = RTS,S trial data)
age_range = 30; % # of days to finish vaccination [9 month, 9 month + age_range]
[~,vage_ind1] = min(abs(P.a-P.vage));
[~,vage_ind2] = min(abs(P.a-(P.vage+age_range)));
% NH = trapz(P.PH_stable)*P.da;
% [SH_EE,~,~,~,~,~,~] = steady_state('EE','numerical');
% if vage_ind1 == vage_ind2
%     P.v0 = P.vyear/365/(P.NN*SH_EE(vage_ind1)/NH); % per-capita daily vaccine rate
% else
%     P.v0 = P.vyear/365/(P.NN*trapz(SH_EE(vage_ind1:vage_ind2)*P.da)/NH); % per-capita daily vaccine rate
% end
v = zeros(size(P.a));
v(vage_ind1:vage_ind2) = P.v0;
v_fun = @(age) P.v0.*(age>=P.a(vage_ind1)).*(age<=P.a(vage_ind2));

P.v = v;
P.v_fun = v_fun;

end