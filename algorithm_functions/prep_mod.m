% David Bernstein

function [model1,inds,e_m,i_m] = prep_mod(model)
% Function to prepare a model to be run by calc_PM_fit_nonlin: find and 
% turn off exchange reactions, find and turn off maintenance reactions, add
% an additional exchange reaction for all metabolites and record the index,
% record which metabolites are intracellular vs extracellular
% Note: prep_mod has not been optimized for all possible model types, it
% works best with models downloaded from KBase.us or the BiGG database

%INPUT
% model: structure, COBRA model

%OUTPUT
% model1: structure, COBRA model after being run through prep_mod, exchange
% reactions off, maintenance reaction off, additional exchange for all
% metabolites added with indices recorded
% inds: int, vector of length of metabolites with the indices of metabolite 
% exchange reactions
% e_m: logical, vector of length of metabolites with 0 for 
% non-extracellular and 1 for extracellular metabolites
% i_m: logical, vector of length of metabolites with 0 for 
% non-intracellular and 1 for intracellular metabolites

%%
% initialize model1
model1 = model;

% find and turn off exchange reactions
for I = 1:length(model1.rxns) % for each reaction
    index = find(model1.S(:,I)~=0); % find indicies of metabolites involved in reaction
    coeffnorm = zeros(length(index),1); % initialize variable to store normalized stoichiometric coefficients
    for J = 1:length(index) % for each metabolite involved in reaction
        coeffnorm(J) = model1.S(index(J),I)/abs(model1.S(index(J),I)); % record normalized stoichiometric coefficients (+ or - 1)
    end
    if length(unique(coeffnorm)) == 1 % if reaction only produces or only consumes metabolites assign it as exchange and turn it off, (note: this also includes reactions typically listed as source/sink reactions)
        % turn off reactions
        model1.lb(I) = 0;
        model1.ub(I) = 0;
    end
end

% find and turn off maintenance reactions
for I = 1:length(model1.rxns) % for each reaction
    if model1.lb(I) > 0 % if flux is forced through reaction in + direction
        model1.lb(I) = 0; % turn off reaction
    elseif model1.ub(I) < 0 % if flux is forced through reaction in - direction
        model1.ub(I) = 0; % turn off reaction
    end
end

% add an additional exchange reaction for all metabolites and keep track of the indicies
inds = uint32(zeros(length(model1.mets),1));%store the index of the exchange reaction associated with each metabolite
for I = 1:length(model1.mets)
    %add exchange reaction for that metabolite
    Name = I;
    model1.rxns = [model1.rxns;Name];
    S_rxn = zeros(size(model1.S,1),1);
    S_rxn(I) = 1; %Reaction adds 1 unit of metabolite to model (note: this is the oposite of typical exchange reactions)
    model1.S = [model1.S,S_rxn];
    model1.lb(length(model1.rxns)) = 0;
    model1.ub(length(model1.rxns)) = 0;
    model1.c(length(model1.rxns)) = 0;
    model1.rev(length(model1.rxns)) = 0;
    inds(I) = length(model1.rxns);
end

% find extracellular metabolites
% Identifiers (this code is set up to identify a common syntax at the end of the metabolite name)
ext_IDs = {'[e]',... M genetalium model iPS189
           '_e0',... KBase models
           '[e0]',...KBase new models
           '_e'};... BiGG models
ext_IDs_L = zeros(length(ext_IDs),1); % variable for ID lengths
for I = 1:length(ext_IDs); ext_IDs_L(I) = length(ext_IDs{I}); end % find lengths of IDs          
e_m = logical(false(length(model1.mets),1));% initialize logical variable for extracellular metabolites
for I = 1:length(model1.mets) % for each metabolite
    probe = model1.mets{I}; % probe metabolite name
    if length(probe)> min(ext_IDs_L) % avoid length errors
        for I1 = 1:length(ext_IDs) % for each candidate ID
            if length(probe) > ext_IDs_L(I1) % avoid length errors
                probe2 = probe(end-(ext_IDs_L(I1)-1):end); % probe the end of the metabolite name
                if strcmp(probe2,ext_IDs{I1}) % check for matches with probe and ID's
                    e_m(I) = true; % extracellular
                end
            end
        end
    end
end

% find intracellular metabolites
% Identifiers (this code is set up to identify a common syntax at the end of the metabolite name)
int_IDs = {'[c]',... M genetalium model iPS189
           '_c0',... KBase models
           '[c0]',...KBase new models
           '_c'};... BiGG models
int_IDs_L = zeros(length(int_IDs),1); % variable for ID lengths
for I = 1:length(int_IDs); int_IDs_L(I) = length(int_IDs{I}); end % find lengths of IDs
i_m = logical(false(length(model1.mets),1));% initialize logical variable for intracellular metabolites
for I = 1:length(model1.mets) % for each metabolite
    probe = model1.mets{I}; % probe metabolite name
    if length(probe)> min(int_IDs_L) % avoid length errors
        for I1 = 1:length(int_IDs) % for each candidate ID
            if length(probe) > int_IDs_L(I1) % avoid length errors
                probe2 = probe(end-(int_IDs_L(I1)-1):end); % probe the end of the metabolite name
                if strcmp(probe2,int_IDs{I1}) % check for matches with probe and ID's
                    i_m(I) = true; % intracellular
                end
            end
        end
    end
end

end
