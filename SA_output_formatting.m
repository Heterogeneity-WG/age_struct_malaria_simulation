function [lP_list,lQ] = SA_output_formatting(lP_list,lQ,flag_dollar)
if flag_dollar
    lP_list = cellfun(@(x) join(["$\",x,"$"],""),lP_list,'UniformOutput',false);
end

for ip = 1:length(lP_list)
    lP_list{ip} = strrep(lP_list{ip},'\dac','d_e');
    lP_list{ip} = strrep(lP_list{ip},'\rD','r_D');
    lP_list{ip} = strrep(lP_list{ip},'\rA','r_A');
    lP_list{ip} = strrep(lP_list{ip},'muM','mu_M');
    lP_list{ip} = strrep(lP_list{ip},'betaM','beta_M');
    lP_list{ip} = strrep(lP_list{ip},'betaA','beta_A');
    lP_list{ip} = strrep(lP_list{ip},'betaD','beta_D');
    lP_list{ip} = strrep(lP_list{ip},'psir2','psi(r_2)');
    lP_list{ip} = strrep(lP_list{ip},'psis2','psi(s_2)');
    lP_list{ip} = strrep(lP_list{ip},'rhor2','rho(r_2)');
    lP_list{ip} = strrep(lP_list{ip},'rhos2','rho(s_2)');
    lP_list{ip} = strrep(lP_list{ip},'phir2','phi(r_2)');
    lP_list{ip} = strrep(lP_list{ip},'phis2','phi(s_2)');
    lP_list{ip} = strrep(lP_list{ip},'\cX','c_X');
    lP_list{ip} = strrep(lP_list{ip},'\dummy','dummy');
end
lQ = strrep(lQ,'_',' ');
end
