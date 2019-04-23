clear all
close all
clc

load('Figure4_Sup5_Data')

diff_prot = m_PM_prot-m_PM_orig;
diff_carb = m_PM_carb-m_PM_orig;

%% Figures
figure(1)
clims = [0,1];
imagesc(diff_prot(2:3,:),clims)
colorbar
yticks([1,2])
yticklabels(strrep(mod_names_carb([2,3]),'_',' '))
xticks(1:97)
xticklabels(metabolite_names)
xtickangle(90)
set(gca,'fontsize',8)
colormap(flipud(hot))

figure(2)
clims = [0,1];
imagesc(diff_carb(2:3,:),clims)
colorbar
yticks([1,2])
yticklabels(strrep(mod_names_carb([2,3]),'_',' '))
xticks(1:97)
xticklabels(metabolite_names)
xtickangle(90)
set(gca,'fontsize',8)
colormap(flipud(hot))
