clc
clear all
close all

[num,txt,raw] = xlsread('.\Figure2_Sup1_Data.xlsx');
PM = num(1:end,1);
ND_in = num(1:end,4);
ND_out = num(1:end,5);
ND_tot = num(1:end,2);
ND_tot_noRev = num(1:end,3);

%% Scatter Plot
 r = 4;
 
 figure(1)
 subplot(2,2,1)
 plot(PM,ND_in,'k.','markersize',7);
 xlabel('PM')
 ylabel('ND in')
[rho1,pval1] = corr(PM,ND_in,'type','spearman');
text(0.02,37,strrep(strcat('spearman''s rho =_',num2str(round(rho1,r))),'_',' '))
text(0.02,32,strrep(strcat('p-val =_',num2str(round(pval1,r))),'_',' '))

 subplot(2,2,2)
 plot(PM,ND_out,'k.','markersize',7);
 xlabel('PM')
 ylabel('ND out')
[rho2,pval2] = corr(PM,ND_out,'type','spearman');
text(0.02,28,strrep(strcat('spearman''s rho =_',num2str(round(rho2,r))),'_',' '))
text(0.02,24,strrep(strcat('p-val =_',num2str(round(pval2,r))),'_',' '))

subplot(2,2,3)
 plot(PM,ND_tot,'k.','markersize',7);
 xlabel('PM')
 ylabel('ND total')
[rho3,pval3] = corr(PM,ND_tot,'type','spearman');
text(0.02,74,strrep(strcat('spearman''s rho =_',num2str(round(rho3,r))),'_',' '))
text(0.02,64,strrep(strcat('p-val =_',num2str(round(pval3,r))),'_',' '))

 subplot(2,2,4)
 plot(PM,ND_tot_noRev,'k.','markersize',7);
 xlabel('PM')
 ylabel('ND total (no Rev)')
[rho4,pval4] = corr(PM,ND_tot_noRev,'type','spearman');
text(0.02,37,strrep(strcat('spearman''s rho =_',num2str(round(rho4,r))),'_',' '))
text(0.02,32,strrep(strcat('p-val =_',num2str(round(pval4,r))),'_',' '))

 