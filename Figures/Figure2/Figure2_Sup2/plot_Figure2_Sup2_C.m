clear all
close all
clc

[num,txt,raw] = xlsread('Figure2_Sup2_Data_C.xlsx');
m_PM = num(1:end,3);
s_mod_aux = num(1:end,1);
f_mod_aux = num(1:end,2);

[rho(1),pval(1)]=corr(s_mod_aux,m_PM,'type','spearman');
[rho(2),pval(2)]=corr(f_mod_aux,m_PM,'type','spearman');    

figure(1)
plot(s_mod_aux,m_PM,'ro','markersize',8);
hold on
plot(f_mod_aux,m_PM,'b.','markersize',20);

xlabel('pathway metric','fontsize',16)
ylabel('PM','fontsize',16)
legend('pathway sum','pathway length')

set(gca,'fontsize',14)
