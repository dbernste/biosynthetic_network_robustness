clear all
close all
clc

load('Figure4_Sup6_Data.mat')
load('Supplementary_Table_4_Data.mat')

%% Load SparCC data
% change data source
data_source = 'Buccal Mucosa_Sparse';
[num,txt,~] = xlsread(strcat('.\Friedman_et_al_selected\',data_source,'_pairs.csv'));

C = num(1:end,1); 
OTUID1 = num(1:end,2);
OTUID2 = num(1:end,3);
LIN1 = txt(2:end,4);
LIN2 = txt(2:end,5);
% clean OTUIDs
OTU_Cat = [OTUID1;OTUID2];
OTU_LIN_Cat = [LIN1;LIN2];
OTU_IDS = unique(OTU_Cat);
OTU_LINS = cell(size(OTU_IDS));
Genus_LINS = cell(size(OTU_IDS));
for I = 1:length(OTU_IDS)
    [~,b] = ismember(OTU_IDS(I),OTU_Cat);
    OTU_LINS(I) = OTU_LIN_Cat(b);
    split = strsplit(OTU_LIN_Cat{b},';');
    Genus_LINS(I) = split(7);
end

%% Match Data
% Aggregate at the Genus Level
[all_genera,ia,ic] = unique(PM_genus);

% map genera
genera_map_PM = zeros(size(PM_genus));
genera_map_SparCC = zeros(size(OTU_IDS));
for I = 1:length(all_genera)
    % search PM for matched genera
    for J = 1:length(PM_genus)
        if contains(PM_genus{J},all_genera{I})==1
            genera_map_PM(J) = I;
        end
    end
    % search SparCC data for matched genera
    for J = 1:length(Genus_LINS)
        if contains(Genus_LINS{J},all_genera{I})==1
            genera_map_SparCC(J) = I;
        end
    end
end

%% Aggregate data

% PM distance
SUM1 = zeros(length(all_genera));
N1 = zeros(length(all_genera));
for I = 1:size(PM_distance,1)
    for J = 1:size(PM_distance,2)
        SUM1(genera_map_PM(I),genera_map_PM(J)) = SUM1(genera_map_PM(I),genera_map_PM(J)) + PM_distance(I,J);
        N1(genera_map_PM(I),genera_map_PM(J)) = N1(genera_map_PM(I),genera_map_PM(J)) + 1;    
    end
end
genera_PM_distance = SUM1./N1;

% PM_complementarity
SUM2 = zeros(length(all_genera));
N2 = zeros(length(all_genera));
for I = 1:size(PM_complementarity,1)
    for J = 1:size(PM_complementarity,2)
        SUM2(genera_map_PM(I),genera_map_PM(J)) = SUM2(genera_map_PM(I),genera_map_PM(J)) + PM_complementarity(I,J);
        N2(genera_map_PM(I),genera_map_PM(J)) = N2(genera_map_PM(I),genera_map_PM(J)) + 1;    
    end
end
genera_PM_complementarity = SUM2./N2;

% Seed distance
SUM1s = zeros(length(all_genera));
N1s = zeros(length(all_genera));
for I = 1:size(Seed_distance,1)
    for J = 1:size(Seed_distance,2)
        SUM1s(genera_map_PM(I),genera_map_PM(J)) = SUM1s(genera_map_PM(I),genera_map_PM(J)) + Seed_distance(I,J);
        N1s(genera_map_PM(I),genera_map_PM(J)) = N1s(genera_map_PM(I),genera_map_PM(J)) + 1;    
    end
end
genera_Seed_distance = SUM1s./N1s;

% Seed competition
SUM1s2 = zeros(length(all_genera));
N1s2 = zeros(length(all_genera));
for I = 1:size(Seed_distance,1)
    for J = 1:size(Seed_distance,2)
        SUM1s2(genera_map_PM(I),genera_map_PM(J)) = SUM1s2(genera_map_PM(I),genera_map_PM(J)) + Seed_competition(I,J);
        N1s2(genera_map_PM(I),genera_map_PM(J)) = N1s2(genera_map_PM(I),genera_map_PM(J)) + 1;    
    end
end
genera_Seed_competition = SUM1s2./N1s2;

% Seed complementarity
SUM1s3 = zeros(length(all_genera));
N1s3 = zeros(length(all_genera));
for I = 1:size(Seed_distance,1)
    for J = 1:size(Seed_distance,2)
        SUM1s3(genera_map_PM(I),genera_map_PM(J)) = SUM1s3(genera_map_PM(I),genera_map_PM(J)) + Seed_complementarity(I,J);
        N1s3(genera_map_PM(I),genera_map_PM(J)) = N1s3(genera_map_PM(I),genera_map_PM(J)) + 1;    
    end
end
genera_Seed_complementarity = SUM1s3./N1s3;

% rxn distance and jaccard
SUM1 = zeros(length(all_genera));
N1 = zeros(length(all_genera));
for I = 1:size(PM_distance,1)
    for J = 1:size(PM_distance,2)
        SUM1(genera_map_PM(I),genera_map_PM(J)) = SUM1(genera_map_PM(I),genera_map_PM(J)) + Rxn_jaccard(I,J);
        SUM2(genera_map_PM(I),genera_map_PM(J)) = SUM2(genera_map_PM(I),genera_map_PM(J)) + Rxn_distance(I,J);
        N1(genera_map_PM(I),genera_map_PM(J)) = N1(genera_map_PM(I),genera_map_PM(J)) + 1;    
    end
end
genera_Rxn_jaccard = SUM1./N1;
genera_Rxn_distance = SUM2./N1;

% SparCC
SUM3 = zeros(length(all_genera));
N3 = zeros(length(all_genera));
for I = 1:length(C)
    [~,b1] = ismember(OTUID1(I),OTU_IDS);
    [~,b2] = ismember(OTUID2(I),OTU_IDS);
    if genera_map_SparCC(b1) ~= 0 && genera_map_SparCC(b2) ~= 0 % genera can be mapped
        [a] = sort([genera_map_SparCC(b1),genera_map_SparCC(b2)],'ascend');
        SUM3(a(1),a(2)) = SUM3(a(1),a(2)) + C(I);
        N3(a(1),a(2)) = N3(a(1),a(2)) + 1;
    end
end
genera_SparCC = SUM3./N3;
genera_SparCC(isnan(genera_SparCC))=0;

%%
mapped_index = unique(genera_map_SparCC);
mapped_index(1) = [];

m_genera_PM_distance = genera_PM_distance(mapped_index,mapped_index);
m_genera_PM_complementarity = genera_PM_complementarity(mapped_index,mapped_index);
m_genera_Seed_distance = genera_Seed_distance(mapped_index,mapped_index);
m_genera_Seed_competition = genera_Seed_competition(mapped_index,mapped_index);
m_genera_Seed_complementarity = genera_Seed_complementarity(mapped_index,mapped_index);
m_genera_Rxn_jaccard = genera_Rxn_jaccard(mapped_index,mapped_index);
m_genera_Rxn_distance = genera_Rxn_distance(mapped_index,mapped_index);
m_genera_SparCC = genera_SparCC(mapped_index,mapped_index);
m_genera_SparCC = m_genera_SparCC + m_genera_SparCC'; % Add to transpose because original only has upper triangle. Some diagonal elements will get messed up but they are never used for correlation calculations so it doesn't matter.

%% Correlations
[n1,n2] = size(m_genera_PM_distance);
METRICSg = zeros(n1,n2,7);
METRICSg(:,:,1) = m_genera_PM_distance;
METRICSg(:,:,2) = m_genera_PM_complementarity;
METRICSg(:,:,3) = m_genera_Seed_distance;
METRICSg(:,:,4) = m_genera_Seed_competition;
METRICSg(:,:,5) = m_genera_Seed_complementarity;
METRICSg(:,:,6) = m_genera_Rxn_distance;
METRICSg(:,:,7) = m_genera_Rxn_jaccard;
CORR_CO = zeros(size(METRICSg,3),1);
PVAL_CO = zeros(size(METRICSg,3),1);
% Plot data
for I = 1:size(METRICSg,3)
    M1 = METRICSg(:,:,I);
    M2 = m_genera_SparCC;
    [CORR_CO(I),PVAL_CO(I)] = f_mantel(M1,M2,triu(ones(size(m_genera_PM_distance)),1)+tril(ones(size(m_genera_PM_distance)),-1),'spearman',10000);
end

%% Partial Correlations
PAR_CORR_CO = zeros(size(METRICSg,3)-1,1);
PAR_PVAL_CO = zeros(size(METRICSg,3)-1,1);
PAR_CORR_CO_I = zeros(size(METRICSg,3)-1,1);
PAR_PVAL_CO_I = zeros(size(METRICSg,3)-1,1);
for I = 2:size(METRICSg,3)
    M1 = METRICSg(:,:,1);
    M2 = m_genera_SparCC;
    M3 = METRICSg(:,:,I);
    [PAR_CORR_CO(I-1),PAR_PVAL_CO(I-1)] = f_mantel_partial(M1,M2,M3,triu(ones(size(m_genera_PM_distance)),1)+tril(ones(size(m_genera_PM_distance)),-1),'spearman',10000);    
    [PAR_CORR_CO_I(I-1),PAR_PVAL_CO_I(I-1)] = f_mantel_partial(M3,M2,M1,triu(ones(size(m_genera_PM_distance)),1)+tril(ones(size(m_genera_PM_distance)),-1),'spearman',10000);    
end

%% Organize Correlation Data
X = [[CORR_CO;PAR_CORR_CO;PAR_CORR_CO_I],[PVAL_CO;PAR_PVAL_CO;PAR_PVAL_CO_I]]
