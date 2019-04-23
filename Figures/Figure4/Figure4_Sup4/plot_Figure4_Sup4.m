clc
clear all
close all

load('./Figure4_Sup4_Data')

%% Plot Data

% convert names to labels
met_labels = {'Acetate';
              'L-Lactate';
              'D-Lactate';
              'Butyrate';
              'Isobutyrate';
              'Succinate';
              'Formate';
              'Propionate';
              'Isovalerate'};
mod_labels = cell(size(mod_names));
for I = 1:length(mod_labels)
    tmp = mod_names{I};
    tmp = strrep(tmp,'_',' ');
    mod_labels{I} = tmp;
end

% Jitter Plot
% sort data by median
med = median(m_PM);
[~,b] = sort(med,'descend');

jitter = 0.25;

% Analyze specific genera
genus = cell(size(mod_names));
species = cell(size(mod_names));
for I = 1:length(mod_names)
    tmp = mod_names{I};
    tmp = strsplit(tmp,'_');
    genus{I} = tmp{1};
    species{I} = tmp{2};
end

% ALL + Fusobacteria, Porphyromonas, Prevotella Butyrate
figure(1)
F1 = ismember(genus,'Fusobacterium');
F2 = ismember(genus,'Porphyromonas');
F3 = ismember(genus,'Prevotella');
hold on
plotted1 = m_PM(F1==1,b);
plotted1 = plotted1(:,7);
plotted2= m_PM(F2==1,b);
plotted2 = plotted2(:,7);
plotted3 = m_PM(F3==1,b);
plotted3 = plotted3(:,7);
plottedA = m_PM(:,b(1:7));
plottedB = m_PM(:,b(8:9));

plotted = plottedA;
for I = 1:size(plotted,2)
    for J = 1:size(plotted,1)
        rand_jitter = ((rand-0.5)*2)*jitter;
        plot(I+rand_jitter,plotted(J,I),'b.','markersize',5);
    end
    plot(I,median(plotted(:,I)),'r.','markersize',16);
end

plotted = plotted1;
for J = 1:size(plotted,1)
    rand_jitter = ((rand-0.5)*2)*jitter;
    plot(8+rand_jitter,plotted(J,1),'m.','markersize',5);
end
plot(8,median(plotted(:,1)),'r.','markersize',16);

plotted = plotted2;
for J = 1:size(plotted,1)
    rand_jitter = ((rand-0.5)*2)*jitter;
    plot(9+rand_jitter,plotted(J,1),'m.','markersize',5);
end
plot(9,median(plotted(:,1)),'r.','markersize',16);

plotted = plotted3;
for J = 1:size(plotted,1)
    rand_jitter = ((rand-0.5)*2)*jitter;
    plot(10+rand_jitter,plotted(J,1),'m.','markersize',5);
end
plot(10,median(plotted(:,1)),'r.','markersize',16);

plotted = plottedB;
for I = 1:size(plotted,2)
    for J = 1:size(plotted,1)
        rand_jitter = ((rand-0.5)*2)*jitter;
        plot(10+I+rand_jitter,plotted(J,I),'b.','markersize',5);
    end
    plot(10+I,median(plotted(:,I)),'r.','markersize',16);
end

met_labels_b = met_labels(b);
met_labels_2 = [met_labels_b(1:7);
                {'Butyrate F'};
                {'Butyrate Po'};
                {'Butyrate Pr'};
                met_labels_b(8:9)];

set(gca,'YLim',[0 1])
ylabel('PM distribution')
xlabel('organic acids')
xticks(1:size(b,2)+3)
xticklabels(met_labels_2);
xtickangle(90)
