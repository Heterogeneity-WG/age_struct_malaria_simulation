global P

flag_save = 1; % flag for saving the results or not (Note: it will overwrite previous results in the folder)
% numerical config
tfinal = 10*365; age_max = 100*365; P.age_max = age_max;
dt = 1; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t; P.tfinal = tfinal;

Malaria_parameters_baseline;
Malaria_parameters_transform;
Malaria_parameters_transform_vac;
[SH, EH, DH, AH, ~, ~, SM, EM, IM, ~, ~, ~, Ctot, ~] = age_structured_Malaria_IC_vac('EE_reset');
PH = SH+EH+DH+AH;
NH = trapz(PH)*P.da;
Ctot_pp = Ctot./PH;
cc = linspace(0,max(Ctot_pp),100);
phi_curve = sigmoid_prob(cc, 'phi');
rho_curve = sigmoid_prob(cc, 'rho');
psi_curve = sigmoid_prob(cc, 'psi');
% population average sigmoids
phi_ave = trapz(sigmoid_prob(Ctot_pp, 'phi').*PH/NH)*P.da;
rho_ave = trapz(sigmoid_prob(Ctot_pp, 'rho').*PH/NH)*P.da;
psi_ave = trapz(sigmoid_prob(Ctot_pp, 'psi').*PH/NH)*P.da;
%% plot link functions
figure_setups; 
%set(gcf,'Position',[353   307   552   446]);
hold on;
plot(cc,phi_curve,'-')
plot(cc,rho_curve,'--')
plot(cc,psi_curve,'-.')
legend('$\phi(\tilde{C}_{H})$','$\rho(\tilde{C}_{H})$','$\psi(\tilde{C}_{H})$','Location','e','fontsize',45)
axis([0 max(cc) 0 1])
xlabel('$\tilde{C}_{H}$','fontsize',45)
ylabel('Probability','fontsize',45)
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
% this sets an 'a)' right at the top left of the axes
text(ax,-12,temp_ax(end)+0.5,'(A)','Units','characters','fontsize',45);
xlim([0 10]);
xticks([0 5 10]);
ylim([0 1]);
yticks([0 0.5 1]);
set(gca,'fontsize', 45) 
ax.GridAlpha = 1;  % Make grid lines transparent..
ax.GridColor = [224, 224, 224]/255; % change grid color
title('Link functions','fontsize',45);
save_string = strcat('fig10_','A','.svg');
if flag_save; saveas(gcf,save_string); end
%% plotting heatmap (age, EIR, immunity level)
tfinal = 20*365; age_max = 100*365; P.age_max = age_max;
dt = 10; da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);
P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t; P.tfinal = tfinal;

Malaria_parameters_baseline;
P.ss_c = 1; P.ss_S0 = 1; % turn off seasonality
Malaria_parameters_transform;
Malaria_parameters_transform_vac;
Malaria_parameters_transform;

var_list = (0.01:0.01:1).^2;
xx = P.a/365;
yy = zeros(1,length(var_list));
zz = zeros(na,length(var_list));
for jj = 1:length(var_list)
    P.betaM = var_list(jj);
    Malaria_parameters_transform;
    [SH, EH, DH, AH, ~, ~, SM, EM, IM, ~, ~, ~, Ctot, ~] = age_structured_Malaria_IC_vac('EE_reset');
    PH = SH+EH+DH+AH;
    NM = SM+EM+IM;
    [bH,~] = biting_rate(PH,NM);
    EIR = bH.*IM./NM*365; % EIR matrix
    EIR_tot = trapz(EIR.*PH)/trapz(PH); % EIR sum over age, at final time
    yy(1,jj) = EIR_tot; % aEIR
    zz(:,jj) = Ctot./PH; % final Ctot at EE    
end

%% plotting heatmap (age, EIR, immunity) 
figure_setups; 
hold on; 
grid off;
%set(gcf,'Position',[353   307   552   446])
imagesc(xx,yy,zz')
xlim([0 20])
ylim([0 max(yy)])
xlabel('Age (years)','fontsize',45);
ylabel('aEIR','fontsize',45);
title('Immunity p.p.','fontsize',45);
set(gca,'YDir','normal');
colormap jet
colorbar('limits',[0 15],'Ticks',[0:5:15],'TickLabelInterpreter','latex');
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
set(gca,'fontsize', 45) 
% this sets an 'a)' right at the top left of the axes
text(ax,-12,temp_ax(end)+2.5,'(B)','Units','characters','fontsize',45);
save_string = strcat('fig10_','B','.svg');
if flag_save; saveas(gcf,save_string); end
%% plotting heatmap (age, EIR, rho)
figure_setups; hold on; grid off
% set(gcf,'Position',[353   307   552   446])
zz_rho = sigmoid_prob(zz, 'rho');
imagesc(xx,yy,zz_rho')
xlabel('Age (years)','fontsize',45);
ylabel('aEIR','fontsize',45);
title('Susceptibility','fontsize',45);
set(gca,'YDir','normal');
colormap jet
colorbar('limits',[0 1],'Ticks',0:0.5:1,'TickLabelInterpreter','latex');
xlim([0 20])
ylim([0 max(yy)])
clim([0 max(max(zz_rho))]);
ax=gca;
% read out the position of the axis in the unit "characters"
set(ax,'Units','characters'); temp_ax=get(ax,'Position');
% this sets an 'a)' right at the top left of the axes
text(ax,-12,temp_ax(end)+2.5,'(C)','Units','characters','fontsize',45);
set(gca,'fontsize', 45) 
save_string = strcat('fig10_','C','.svg');
if flag_save; saveas(gcf,save_string); end
