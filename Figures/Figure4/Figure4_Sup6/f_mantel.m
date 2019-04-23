function [rho0,pval0] = f_mantel(matrix1,matrix2,matrixC,type,permutations)
% this function computes the correlation between two matricies using a
% permutation test

% output
% pval - the fraction of permutations with pvalue above the input matricies
% input
% matrix1 - first matrix to test
% matrix2 - second matrix to test
% matrixC - a matrix of ones and zeros with 1 for all elements that should
% be included in the correlation and 0 for all that should not.
% type - type of correlation to test
% permutations - number of permutations to use

[rho0,~] = corr(matrix1(matrixC==1),matrix2(matrixC==1),'type',type);
N = permutations; % total number of permutations
n = 0; % number of permutations with correlation greater than original
for I = 1:permutations
    p1 = randperm(size(matrix1,1));
    p2 = randperm(size(matrix1,2));
    matrix1p = matrix1(p1,p2); % permuted matrix 1
    [rho,~] = corr(matrix1p(matrixC==1),matrix2(matrixC==1),'type',type);
    if abs(rho) >= abs(rho0)
        n = n + 1;
    end
end
pval0 = (n+1)/(N+1);

end
