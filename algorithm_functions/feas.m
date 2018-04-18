%David Bernstein

function [out] = feas(target,model1,inds,added,ineq)
% Function to determine if the production of a given target metabolite by a 
% given metabolic model with specified constraints is feasible or not

% INPUT
% target: logical, vector of length of metabolites with 0 for non-target 
% metabolites and 1 for target metaboltes
% model1: structure, COBRA model after being run through prep_mod, exchange 
% reactions off, maintenance reaction off, additional exchange for all 
% metabolites added with indices recorded
% inds: int, vector of length of metabolites with the indices of metabolite 
% exchange reactions
% added: logical, vector of length of metabolites with 0 for metabolites 
% not added to the model and 1 for metabolites added to the model
% ineq: logical, vector of length of metabolites with 0 for equality mass 
% balance (S*V>0) and 1 for inequality mass balance (S*V>=0)

% OUTPUT
% out: logical, output of feas, 0 if target metabolite production is not 
% feasible and 1 if target metabolite production is feasible

% HARD CODED VARIABLES
% thresh: double, threshold of target metabolite production above which the 
% model is said to be feasible

%%
% Update bounds of metabolite exchange reactions based on added and ineq
model1.ub(inds(added==1)) = 1000; %Allow import of metabolite
model1.lb(inds(ineq==1)) = -1000; %Allow export of metabolite

%Change objective function to target metabolites
% Reset old objective function
model1.c = zeros(length(model1.rxns),1);
% add new reaction consuming all target metabolites
Name = 'TARGET';
model1.rxns = [model1.rxns;Name];
S_rxn = zeros(size(model1.S,1),1);
S_rxn(target==1) = -1; %Reaction consums 1 unit of target metabolites
model1.S = [model1.S,S_rxn];
model1.lb(length(model1.rxns)) = 0;
model1.ub(length(model1.rxns)) = 1000;
model1.c(length(model1.rxns)) = 1;
model1.rev(length(model1.rxns)) = 0;
% Make new reaction objective function
model2 = changeObjective(model1,'TARGET'); % Set objective function to new reaction

% Optimize model to asses feasibility
out = false;
sol = optimizeCbModel(model2);
thresh = 0.001;
if sol.f > thresh
    out = true; % feasible
end

end