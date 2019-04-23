clear all
close all
clc

load('Figure4_Sup7_Data.mat')

%% Label
label = cell(size(PM_complementarity));
for I = 1:size(label,1)
    for J = 1:size(label,2)
        label{I,J} = 'All to All';
    end
end
ids = [1:16]; % other uncultivated
for I = 1:size(label,1)
    for J = 1:length(ids)
        label{I,ids(J)} = 'All to Fastidious';
    end
end
ids = [6:13]; % Mycoplasma
for I = 1:size(label,1)
    for J = 1:length(ids)
        label{I,ids(J)} = 'All to Mycoplasma';
    end
end
ids = [2:4]; % TM7
for I = 1:size(label,1)
    for J = 1:length(ids)
        label{I,ids(J)} = 'All to TM7';
    end
end

%% Plot
figure()
h = scatterhist(flipud(PM_complementarity(:)),flipud(Seed_complementarity(:)),'Group',flipud(label(:)),...
    'Kernel','off','PlotGroup','off','Location','SouthEast','Direction','out'...
    ,'Color','brmk','Marker','...','MarkerSize',[7,7,7],...
    'Style','bar','NBins',100);
xlabel('PM complementarity')
ylabel('Seed complementarity')
set(h(2),'YScale','log')
set(h(3),'YScale','log')
h(2).Children(1).FaceColor = [0.5 0.5 0.5];
h(3).Children(1).FaceColor = [0.5 0.5 0.5];
h(2).Children(1).EdgeColor = [0.5 0.5 0.5];
h(3).Children(1).EdgeColor = [0.5 0.5 0.5];
h(2).Children(1).FaceAlpha = 1;
h(3).Children(1).FaceAlpha = 1;

