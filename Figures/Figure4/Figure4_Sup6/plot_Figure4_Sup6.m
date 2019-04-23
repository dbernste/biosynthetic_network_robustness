clear all
close all
clc

load('Figure4_Sup6_Data.mat')

%% Metric correlations
[n1,n2] = size(PM_distance);
METRICS = zeros(n1,n2,7);
METRICS(:,:,1) = PM_distance;
METRICS(:,:,2) = PM_complementarity;
METRICS(:,:,3) = Seed_distance;
METRICS(:,:,4) = Seed_competition;
METRICS(:,:,5) = Seed_complementarity;
METRICS(:,:,6) = Rxn_distance;
METRICS(:,:,7) = Rxn_jaccard;
CORR1 = zeros(size(METRICS,3));
PVAL1 = zeros(size(METRICS,3));
for I = 1:size(METRICS,3)
    I
    for J = I+1:size(METRICS,3)
        J
        M1 = METRICS(:,:,I);
        M2 = METRICS(:,:,J);
        [CORR1(I,J),PVAL1(I,J)] = f_mantel(M1,M2,triu(ones(size(M1)),1)+tril(ones(size(M1)),-1),'spearman',1000);
    end
end

figure(1)
for I = 1:size(METRICS,3)
    I
    for J = I+1:size(METRICS,3)
        J
        M1 = METRICS(:,:,I);
        M2 = METRICS(:,:,J);
        subplot(size(METRICS,3),size(METRICS,3),(J-1)*size(METRICS,3)+I)
        plot(M1(:),M2(:),'k.','markersize',1)
        title(strcat('rho=',num2str(CORR1(I,J)),', pval=',num2str(PVAL1(I,J))),'fontsize',6)
    end
end

