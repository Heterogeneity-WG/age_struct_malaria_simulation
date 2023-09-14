[DALY,YLL,YLD] = DALY_cal(SH, EH, DH, AH, VH, UH, SM, EM, IM, Ctot);
figure(4);
subplot(2,2,1), plot(t/365,DALY,t/365,YLL,t/365,YLD);
grid on
xlabel('Year')
legend('DALY','YLL (death)','YLD (disability)')