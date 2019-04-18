% David Bernstein

% This function is used to implement the combinatorial theory and calculate
% the output probability of a metabolite based on its minimal precursor set
% structure (S) and the input probability of each metabolite in the minimal
% precursor sets (Pin)

function [Pout] = f_calc_prob(S,Pin)
%INPUT
% S
% -S is a matrix of all of the sets of metabolite that allow for the
% production of the metabolite of interest. Each column corresponds to a
% particular metabolite and each row corresponds to a particular
% minimal precursor set, with 1 meaning that metabolite is included in that set.
% The rows can be though of as or relationships while the columns can be
% thought of as and relationships pertaining to the necesity of these
% metabolites to synthesize a particular metabolite

% Pin
% -Pin is a vector of the input probability of each metabolite it has 
% length of size(S,2)

%OUTPUT
% Pout
% -Pout is the output probability of synthesizing the metabolite of interest

SUM = 0;
for i = 1:size(S,1)
    mult = (-1).^(i+1); %define multiplier for inclusion exclusion principle
    jvect = zeros(1,i);
    SUM = SUM + mult * f_r_calc_prob(1,i,0,jvect,S,Pin); %recursive sum to enumerat binomial coeficient
end
Pout = SUM;

end




