%David Bernstein

function [prob] = prob(target,model1,inds,bern,ineq,samp)
% Function to determine the probability of a target metabolite set being 
% producible given the added metabolite probabilities and the given 
% metabolic network and specified constraints

% INPUT
% target: logical, vector of length of metabolites with 0 for non-target 
% metabolites and 1 for target metaboltes
% model1: structure, COBRA model after being run through prep_mod, exchange 
% reactions off, maintenance reaction off, additional exchange for all 
% metabolites added with indices recorded
% inds: int, vector of length of metabolites with the indices of metabolite 
% exchange reactions
% bern: double, vector of length of metabolites with Bernoulli parameter 
% for each metabolite
% ineq: logical, vector of length of metabolites with 0 for equality mass 
% balance (S*V>0) and 1 for inequality mass balance (S*V>=0)
% samp: int, number of samples taken to estimate probability of production

% OUTPUT
% prob: double, output of prob, probability of target metabolite being 
% producible given the metabolic network and input metabolite probabilities

%%
N = uint32(0); % variable for number of succesful trials
for J = 1:samp % for each sample
    added = rand_add(bern); % randomly add input metabolites based on p_in
    out = feas(target,model1,inds,added,ineq); % assess feasibility of producing target metabolite
    if out == 1
        N = N + 1; % count succesful trials
    end
end
prob = double(N)/double(samp); % estimate bernoilli parameter p_out

end
