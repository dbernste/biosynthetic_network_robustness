clear all
close all
clc

%% PM
load('calculate_Metrics_PM_Data.mat')
% PM distance
PM_distance = zeros(size(m_PM,2));
for i = 1:size(m_PM,2)
    for j = 1:size(m_PM,2)
        % calculate distance
        PM_distance(i,j) = norm((m_PM(:,i)-m_PM(:,j)),1);
    end
end
% PM complementarity
PM_complementarity = zeros(size(m_PM,2),size(m_PM,2));
for I1 = 1:size(m_PM,2)
    for I2 = 1:size(m_PM,2)
        num = 0;
        den = 0;
        for J = 1:size(m_PM,1)
            if GN(J) == 1 && GP(J) == 1
                num = num + max(m_PM(J,I1),m_PM(J,I2)) - m_PM(J,I2);
                den = den + m_PM(J,I1);
            end
        end
        PM_complementarity(I1,I2) = num/den;
    end
end

%% Seed
% Seed complementarity
Seed_complementarity = dlmread('calculate_Metrics_Seed_Complementarity.tsv','\t');
Seed_complementarity = Seed_complementarity';
fid = fopen('calculate_Metrics_ET_names.txt','r');
count = 1;
tline = fgetl(fid);
tline = tline(4:end-8);
Seed_mods{1,1} = tline;
while ischar(tline)
    tline = fgetl(fid);
    tline = tline(4:end-8);
    count = count + 1;
    Seed_mods{count,1} = tline;
end
fclose(fid);
Seed_mods(end) = [];
[~,sort_Seed] = ismember(PM_mods,Seed_mods);
Seed_mods_S = Seed_mods(sort_Seed);
Seed_complementarity = Seed_complementarity(sort_Seed,sort_Seed);

% Seed distance and competition
[num,~,~] = xlsread('calculate_Metrics_Seeds.csv');
seeds = num;
seeds_S = seeds(:,sort_Seed);
% Seed distance
Seed_distance = zeros(size(m_PM,2));
for i = 1:size(seeds_S,2)
    for j = 1:size(seeds_S,2)
        Seed_distance(i,j) = norm((seeds_S(:,i)-seeds_S(:,j)),1);
    end
end
% Seed competition
Seed_competition = zeros(size(seeds_S,2));
for i = 1:size(seeds_S,2)
    for j = 1:size(seeds_S,2)
        Seed_competition(i,j) = length(find(seeds_S(:,i)~=0))-length(intersect(find(seeds_S(:,i)~=0),find(seeds_S(:,j)~=0)));
    end
end

%% Rxn
load('calculate_Metrics_rxn_Data.mat')
rxn_mods = strrep(rxn_mods,'.mat','');
[~,sort_rxn] = ismember(PM_mods,rxn_mods);
rxn_mods_S = rxn_mods(sort_rxn);

% rxn jaccard
Rxn_jaccard = squareform(pdist(rxnmat,'jaccard'));
Rxn_jaccard = Rxn_jaccard(sort_rxn,sort_rxn);

% rxn distance
rxnmat_T = rxnmat';
Rxn_distance = zeros(size(rxnmat_T,2));
for i = 1:size(rxnmat_T,2)
    for j = 1:size(rxnmat_T,2)
        % calculate distance
        Rxn_distance(i,j) = norm((rxnmat_T(:,i)-rxnmat_T(:,j)),1);
    end
end
Rxn_distance = Rxn_distance(sort_rxn,sort_rxn);

%% Save
save('Figure4_Sup6_Data1','PM_distance','PM_complementarity','Seed_complementarity','Seed_distance','Seed_competition','Rxn_jaccard','Rxn_distance');
