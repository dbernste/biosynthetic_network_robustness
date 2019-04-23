clc
clear all
close all

[num,txt,raw] = xlsread('./Figure3_Sup1_Data.xlsx');
m_PM = num;
met_names = txt(1,2:end);
mod_names = txt(2:end,1);

%% Plot Data
mod_labels = cell(size(mod_names));
for I = 1:length(mod_names)
    tmp = strrep(mod_names{I},'_',' ');
    tmp = strrep(tmp,'iJO1366 AA Aux ','');
    mod_labels{I} = tmp;
end
met_labels = cell(size(met_names));
for I = 1:length(met_names)
    tmp = strrep(met_names{I},'_c',' ');
    tmp = strrep(tmp,'__',' ');
    tmp = strrep(tmp,'_',' ');
    met_labels{I} = tmp;
end

% cluster
[H2,T2,P2] = dendrogram(linkage(pdist(m_PM','euclidean'),'average'),0);
[H1,T1,P1] = dendrogram(linkage(pdist(m_PM,'euclidean'),'average'),0);
m_PM_c = m_PM(P1,P2);
mod_labels_c = mod_labels(P1);
met_labels_c = met_labels(P2);

%
figure(1)
imagesc(m_PM_c)
xticks(1:length(met_labels_c))
xticklabels(met_labels_c)
xtickangle(90)
yticks(1:length(mod_labels_c))
yticklabels(mod_labels_c)
colormap hot
colorbar
xlabel('metabolites')
ylabel('auxotrophs')
set(gca,'fontsize',8)

%% Load Experimental Data
data = zeros(4,2,48,48);
for I = 1:8
    [num,~,~] = xlsread('.\Wintermute_et_al_edited.xls',I+1);
    a = floor((I-1)/2)+1;
    b = mod(I-1,2)+1;
    data(a,b,:,:) = num;
end
[~,txt,~] = xlsread('.\Wintermute_et_al_edited.xls',10);
Aux = txt(2:end,1);

total_growth = squeeze(mean(data(4,1:2,:,:),2));

%% Compare Data
% align data
[a,b] = ismember(mod_labels,Aux);

% Distance metric
d = squareform(pdist(m_PM(2:end,:),'euclidean'));
growth = (squeeze(data(4,1,:,:))+squeeze(data(4,2,:,:)))./2; %mean of 2 replicates
growth_a = growth(b(2:end),b(2:end));
count = 0;
distance = zeros(size(m_PM,1)-1,size(m_PM,1)-1);
for i = 1:size(m_PM,1)-1
    for j = i+1:size(m_PM,1)-1 %all combinations of 2 not including self
        count = count + 1;
        % calculate distance
        distance(i,j) = norm((m_PM(i+1,:)-m_PM(j+1,:)),1);
    end
end
[rho_mantel, pval_mantel] = f_mantel(distance,(growth_a+growth_a')./2,triu(ones(size(distance)),1),'spearman',10000)

% cooperation = zeros(size(distance));
% for i = 1:size(m_PM,1)-1
%     for j = 1:size(m_PM,1)-1 %all combinations of 2 not including self
%         count = count + 1;
%         % calculate distance
%         num = 0;
%         den = 0;
%         for k = 1:size(m_PM,2)
%             num = num + max(m_PM(i,k),m_PM(j,k))-m_PM(j,k);
%             den = den + m_PM(i,k);
%         end
%         cooperation(i,j) = num/den;
%     end
% end
% [rho3_mantel, pval3_mantel] = f_mantel(cooperation,(growth_a+growth_a')./2,triu(ones(size(cooperation)),1)+tril(ones(size(cooperation)),-1),'spearman',10000)


%% Plot data
% Amino Acids [auxotroph index, corresponding metabolite index]
AA = [2,7;
      3,7;
      4,7;
      5,7;
      6,19;
      7,19;
      8,19;
      9,27;
      12,32;
      13,32;
      14,32;
      15,32;
      17,63;
      17,33;
      17,48;
      18,63;
      18,33;
      18,48;
      19,35;
      20,35;
      21,35;
      22,36;
      23,37;
      28,50;
      40,53;
      41,53;
      42,59;
      43,59;
      44,59;
      45,59;
      46,59;
      47,60];
AA1 = unique(AA(:,1));
AA2 = unique(AA(:,2));
AA2 = AA2([1:4,14,5:6,9,7:8,10:13]);
m_PM_AA = m_PM(AA1,AA2);
mod_labels_AA = mod_labels(AA1);
met_labels_AA = met_labels(AA2);
figure(2)
imagesc(m_PM_AA)
xticks(1:length(met_labels_AA))
xticklabels(met_labels_AA)
xtickangle(90)
yticks(1:length(mod_labels_AA))
yticklabels(mod_labels_AA)
colormap hot
colorbar
axis square
xlabel('metabolites')
ylabel('auxotrophs')
set(gca,'fontsize',8)

% plot PM for partner amino acid vs fold growth
growth_AA = growth(b(AA1),b(AA1));
Aux_AA = Aux(b(AA1));

figure(3)
imagesc(growth_AA)
xticks(1:length(Aux_AA))
xticklabels(Aux_AA)
xtickangle(90)
yticks(1:length(Aux_AA))
yticklabels(Aux_AA)
colormap hot
colorbar
axis square
xlabel('auxotrophs')
ylabel('auxotrophs')
set(gca,'fontsize',8)

%% average growth vs PM

% compare average PM with average growth
avg_PM = mean(m_PM(2:end,:),2);
avg_grow = mean((squeeze(data(4,1,:,:))+squeeze(data(4,2,:,:)))./2);
avg_grow_a = avg_grow(b(2:end))';
[rho_avg,pval_avg] = corr(avg_PM,avg_grow_a,'type','spearman')

% Trp
PM_trp = m_PM(42:46,59);
avg_g_trp = mean(growth(b(42:46),:),2);
[rho_trp,pval_trp] = corr(PM_trp,avg_g_trp,'type','Spearman')
figure(4)
plot(PM_trp,avg_g_trp,'.','markersize',20)
xlabel('trp PM')
ylabel('trp average growth')
axis square
set(gca,'fontsize',8)
