function [lP_list,lQ,lQ_title] = SA_output_formatting(lP_list,lQ,flag_dollar)
if flag_dollar
    lP_list = cellfun(@(x) join(["$\",x,"$"],""),lP_list,'UniformOutput',false);
end

for ip = 1:length(lP_list)
    lP_list{ip} = strrep(lP_list{ip},'\dac','d_e');
    lP_list{ip} = strrep(lP_list{ip},'\rD','r_D');
    lP_list{ip} = strrep(lP_list{ip},'\rA','r_A');
    lP_list{ip} = strrep(lP_list{ip},'\muM','\mu_M');
    lP_list{ip} = strrep(lP_list{ip},'betaM','beta_M');
    lP_list{ip} = strrep(lP_list{ip},'betaA','beta_A');
    lP_list{ip} = strrep(lP_list{ip},'betaD','beta_D');
    lP_list{ip} = strrep(lP_list{ip},'\psir2','r_{\psi}');
    lP_list{ip} = strrep(lP_list{ip},'\psis2','s_{\psi}');
    lP_list{ip} = strrep(lP_list{ip},'\rhor2','r_{\rho}');
    lP_list{ip} = strrep(lP_list{ip},'\rhos2','s_{\rho}');
    lP_list{ip} = strrep(lP_list{ip},'\phir2','r_{\phi}');
    lP_list{ip} = strrep(lP_list{ip},'\phis2','s_{\phi}');
    lP_list{ip} = strrep(lP_list{ip},'\c','c');
    lP_list{ip} = strrep(lP_list{ip},'cD','c_D');
    lP_list{ip} = strrep(lP_list{ip},'cE','c_E');
    lP_list{ip} = strrep(lP_list{ip},'cA','c_A');
    lP_list{ip} = strrep(lP_list{ip},'cS','c_S');
    lP_list{ip} = strrep(lP_list{ip},'cU','c_U');
    if strcmp(lP_list{ip}, '$\m$');  lP_list{ip} = '$m_0$'; end
    lP_list{ip} = strrep(lP_list{ip},'\uc','\gamma');
    lP_list{ip} = strrep(lP_list{ip},'\dummy','dummy');
    lP_list{ip} = strrep(lP_list{ip},'\v0','\nu_0');
    lP_list{ip} = strrep(lP_list{ip},'\w','d_\nu'); % sample the period d_\nu (notation in paper), not the wanning rate w (notation in the code)
    lP_list{ip} = strrep(lP_list{ip},'\etas','\eta');
end

lQ_title = lQ;
for iq = 1:length(lQ)
    if strcmp(lQ{iq}, 'EE-death');  lQ_title{iq} = 'malaria death'; end
    if strcmp(lQ{iq}, 'EE-EIR');  lQ_title{iq} = 'EIR'; end
    if strcmp(lQ{iq}, 'EE-death-rate');  lQ_title{iq} = 'malaria death incidence'; end
    if strcmp(lQ{iq}, 'EE-death-05-17');  lQ_title{iq} = 'malaria death ($5 \sim 17$m)'; end
    if strcmp(lQ{iq}, 'EE-death-09-24');  lQ_title{iq} = 'malaria death ($9\sim 24$m)'; end
    if strcmp(lQ{iq}, 'EE-death-02-10');  lQ_title{iq} = 'malaria death ($2\sim 10$y)'; end
    if strcmp(lQ{iq}, 'EE-death-10+');  lQ_title{iq} = 'malaria death (10y+)'; end
    if strcmp(lQ{iq}, 'EE-D');  lQ_title{iq} = 'symptomatic infection'; end
    if strcmp(lQ{iq}, 'EE-A');  lQ_title{iq} = 'asymptomatic infection'; end
    if strcmp(lQ{iq}, 'EE-DA');  lQ_title{iq} = 'malaria prevalence'; end
    if strcmp(lQ{iq}, 'EE-DA-05-17');  lQ_title{iq} = 'malaria prevalence ($5\sim 17$m)'; end
    if strcmp(lQ{iq}, 'EE-DA-09-24');  lQ_title{iq} = 'malaria prevalence ($9\sim 24$m)'; end
    if strcmp(lQ{iq}, 'EE-DA-02-10');  lQ_title{iq} = 'malaria prevalence ($2\sim 10$y)'; end
    if strcmp(lQ{iq}, 'EE-DA-10+');  lQ_title{iq} = 'malaria prevalence (10y+)'; end
end

end
