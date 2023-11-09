figure_setups; hold on
[bH,bM] = biting_rate(PH,NM);
lamH = FOI_H(bH,IM,NM);
lamM = FOI_M(bM,DH,AH,NH);
plot(t/365,trapz(lamH.*PH,1)./NH*P.da);
plot(t/365,lamM);
legend('$\lambda_H$','$\lambda_M$')

figure_setups;
plot(t/365, trapz(lamH.*SH,1)*P.da)
legend('$\lambda_H S_H$')