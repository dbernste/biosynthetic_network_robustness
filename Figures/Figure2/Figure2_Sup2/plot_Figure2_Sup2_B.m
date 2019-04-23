clc
clear all
close all

[num,txt,raw] = xlsread('.\Figure2_Sup2_Data_B.xlsx');

m_PM = num(1:68,5:14);

figure(1)
hold on
plot(10,m_PM(32,1),'color','b','marker','.','markersize',20)
for I = 2:size(m_PM,2)
    plot((I-1),m_PM(32,I),'color','b','marker','.','markersize',20)
end
ylabel('PM')
axis([1,10,0,1])
xnames ={'aux 1','aux 2','aux 3','aux 4','aux 5','aux 6','aux 7','aux 8','aux 9','WT'};
set(gca,'xtick',1:10,'xticklabels',[xnames],'XTickLabelRotation',90,'fontsize',14)
% Plot theoretical values
s = 0:8;
theory_PM_set = (1-0.5).^(1./s);
for I = 1:size(theory_PM_set,2)
    plot((I),theory_PM_set(I),'color','k','marker','x','markersize',12)
end
