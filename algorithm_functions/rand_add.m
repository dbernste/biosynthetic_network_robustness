%David Bernstein

function [added] = rand_add(bern)
% Function to give a random sample of input metabolites based on a 
% Bernoulli distribution for each metabolite

% INPUT
% bern: double, vector of length of metabolites with Bernoulli parameter 
% for each metabolite

% OUTPUT
% added: logical, vector of length of metabolites with 0 for metabolites 
% not added to the model and 1 for metabolites added to the model

added = zeros(length(bern),1);
for I = 1:length(bern)
    temp = rand;
    if temp < bern(I)
        added(I,1) = true;
    end
end

end
