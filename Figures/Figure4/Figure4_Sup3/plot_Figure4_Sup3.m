clc
clear all
close all

% load data
load('./Figure4_Sup3_Data')

%%
% Re-order metadata into original order
% oo = num1(:,1);
[~,b] = sort(oo,'ascend');
gSize = gSize(b);
gram_stain = gram_stain(b);
phylum = phylum(b);
class = class(b);

[~,b1] = sort(oo_m,'ascend');
[~,b2] = sort(oo_o,'ascend');
data = data(b1,b2);
metNames = metNames(b1);

%%
log_gSize = log(gSize);

% fit model to predict the PM for each metabolite using a quadrative fit on
% log genome size and phylogenetic parameters

% Base model
LL = zeros(size(data,1),1);
for I = 1:size(data,1)
    tbl = table(data(I,:)',log_gSize,'VariableNames',{'m_PM','log_gSize'});
    lm = fitlm(tbl,'m_PM~log_gSize+log_gSize^2');
    LL(I) = lm.LogLikelihood;
end

%% Gram Stain
LL_GSP = zeros(size(data,1),1);
COEF_GSP = zeros(size(data,1),1);
% reconstruct binary variable
bin = cell(length(gram_stain),1);
for P1 = 1:length(gram_stain)
    if strcmp(gram_stain{P1},'Positive')==1
        bin{P1} = 'yes';
    else
        bin{P1} = 'no';
    end
end
for I = 1:size(data,1)
    tbl = table(data(I,:)',log_gSize,'VariableNames',{'m_PM','log_gSize'});
    tbl.GS = nominal(bin);
    lm = fitlm(tbl,'m_PM~GS+log_gSize+log_gSize^2');
    LL_GSP(I) = lm.LogLikelihood;
    COEF_GSP(I) = lm.Coefficients{3,1};
end

LL_GSN = zeros(size(data,1),1);
COEF_GSN = zeros(size(data,1),1);
% reconstruct binary variable
bin = cell(length(gram_stain),1);
for P1 = 1:length(gram_stain)
    if strcmp(gram_stain{P1},'Negative')==1
        bin{P1} = 'yes';
    else
        bin{P1} = 'no';
    end
end
for I = 1:size(data,1)
    tbl = table(data(I,:)',log_gSize,'VariableNames',{'m_PM','log_gSize'});
    tbl.GS = nominal(bin);
    lm = fitlm(tbl,'m_PM~GS+log_gSize+log_gSize^2');
    LL_GSN(I) = lm.LogLikelihood;
    COEF_GSN(I) = lm.Coefficients{3,1};
end

%% Phylum
taxa = phylum;
u_taxa = unique(taxa);
u_Phy = u_taxa;
LL_Phy = zeros(size(data,1),length(u_taxa));
COEF_Phy = zeros(size(data,1),length(u_taxa));
for P = 1:length(u_taxa)
    P
    % reconstruct binary variable
    bin = cell(length(taxa),1);
    for P1 = 1:length(taxa)
        if strcmp(taxa{P1},u_taxa{P})==1
            bin{P1} = 'yes';
        else
            bin{P1} = 'no';
        end
    end
    for I = 1:size(data,1)
        tbl = table(data(I,:)',log_gSize,'VariableNames',{'m_PM','log_gSize'});
        tbl.P = nominal(bin);
        lm = fitlm(tbl,'m_PM~P+log_gSize+log_gSize^2');
        LL_Phy(I,P) = lm.LogLikelihood;
        COEF_Phy(I,P) = lm.Coefficients{3,1};
    end
end

%% Class
taxa = class;
u_taxa = unique(taxa);
u_Cla = u_taxa;
LL_Cla = zeros(size(data,1),length(u_taxa));
COEF_Cla = zeros(size(data,1),length(u_taxa));
for P = 1:length(u_taxa)
    P
    % reconstruct binary variable
    bin = cell(length(taxa),1);
    for P1 = 1:length(taxa)
        if strcmp(taxa{P1},u_taxa{P})==1
            bin{P1} = 'yes';
        else
            bin{P1} = 'no';
        end
    end
    for I = 1:size(data,1)
        tbl = table(data(I,:)',log_gSize,'VariableNames',{'m_PM','log_gSize'});
        tbl.P = nominal(bin);
        lm = fitlm(tbl,'m_PM~P+log_gSize+log_gSize^2');
        LL_Cla(I,P) = lm.LogLikelihood;
        COEF_Cla(I,P) = lm.Coefficients{3,1};
    end
end

%% Order
% taxa = order;
% u_taxa = unique(taxa);
% u_Ord = u_taxa;
% LL_Ord = zeros(size(data,1),length(u_taxa));
% COEF_Ord = zeros(size(data,1),length(u_taxa)); 
% for P = 1:length(u_taxa)
%     P
%     % reconstruct binary variable
%     bin = cell(length(taxa),1);
%     for P1 = 1:length(taxa)
%         if strcmp(taxa{P1},u_taxa{P})==1
%             bin{P1} = 'yes';
%         else
%             bin{P1} = 'no';
%         end
%     end
%     for I = 1:size(data,1)
%         tbl = table(data(I,:)',log_gSize,'VariableNames',{'m_PM','log_gSize'});
%         tbl.P = nominal(bin);
%         lm = fitlm(tbl,'m_PM~P+log_gSize+log_gSize^2');
%         LL_Ord(I,P) = lm.LogLikelihood;
%         COEF_Ord(I,P) = lm.Coefficients{3,1};
%     end
% end

%% Plot Data on Heatmap
LLMATRIX = ([LL_GSP,LL_GSN,LL_Phy,LL_Cla]-LL);
LLMATRIX(isnan(LLMATRIX)) = 0;

COEFMATRIX = ([COEF_GSP,COEF_GSN,COEF_Phy,COEF_Cla]);

% order phyla
phyOrd = [7,3,12,9,1,6,5,2,11,4,8,10];
phyNum = [2,2,1,1,1,5,1,5,1,1,1,1];
% order classes
claOrd = [7,1;2,3;9,1;6,3;3,1;2,2;12,1;6,1;7,2;2,5;2,4;6,5;3,2;5,1;2,1;1,1;6,4;6,2;11,1;4,1;8,1;10,1];

ORD = zeros(size(LLMATRIX,2),1);
ORD(1) = 1;
ORD(2) = 2;
for I = 1:length(phyOrd)
    sum = 0;
    for J = 1:length(phyOrd)
        if phyOrd(I) > phyOrd(J)
            sum = sum + phyNum(J) + 1;
        end
    end
    sum = sum + 1;
    ORD(I+2) = sum + 1;
end
for I = 1:length(claOrd)
    p = claOrd(I,1);
    c = claOrd(I,2);
    ind = find(phyOrd==p);
    ORD(I+2+12) = ORD(ind+2)+c;
end
[~,ORD1] = sort(ORD,'ascend');

COEFMATRIX_Ord = COEFMATRIX(:,ORD1);

colormap(flipud(gray))

imagesc(LLMATRIX(:,ORD1));
%colorbar
xticks(1:35);
yticks(1:88);
set(gca,'Xticklabel',[]);
set(gca,'Yticklabel',[]);
colorbar
hold on
% plot grid
for I = 1:88
    plot([1-0.5,36+0.5],[I+0.5,I+0.5],'color',[0.5 0.5 0.5])
    plot([1-0.5,36+0.5],[I-0.5,I-0.5],'color',[0.5 0.5 0.5])
end
for J = 1:36
    plot([J+0.5,J+0.5],[1-0.5,88+0.5],'color',[0.5 0.5 0.5])
    plot([J-0.5,J-0.5],[1-0.5,88+0.5],'color',[0.5 0.5 0.5])
end

labels = {'Gram Positive','Gram Negative',u_Phy{:},u_Cla{:}};
labels_ORD = labels(ORD1);

%% Likelihood ratio test
LL_2 = [LL_GSP,LL_GSN,LL_Phy,LL_Cla];
LL_2_Ord = LL_2(:,ORD1);

alpha = (10^-6)/length(LLMATRIX(:));

for I = 1:size(LL_2_Ord,2)
    [h(:,I),pValue(:,I)] = lratiotest(LL_2_Ord(:,I),LL,2,alpha);
end

for I = 1:size(h,1)
    for J = 1:size(h,2)
        if h(I,J) == true
            if COEFMATRIX_Ord(I,J) < 0
                plot(J,I,'color',[1 0 0],'marker','x','markersize',4)
            else
                plot(J,I,'color',[0 0 1],'marker','.','markersize',5)
            end
        end
    end
end


