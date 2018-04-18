% David Bernstein

function find_PM_mods_mets(results_location,models_location,mets,metnum,modnum,fixon,fixoff,a,s,samp,noise,n,thresh,runs)
% Function to run the producibility analysis for a set of models and
% metabolites, and allow for parallelization

% INPUT
% results_location: char, directory where final results will be saved
% models_location: char, directory containing model files
% mets: cell array - char, metabolite names that will be analyzed as target
% metabolites 
% metnum: int, index of target metabolite to analyze, used for
% parallelization 
% modnum: int, index of model to analyze, used for parallelization
% fixon: cell array - char, metabolite names that will always be added
% fixoff: cell array - char, metabolite names that will never be added
% a: int, indicator variable for metabolites randomly added during
% producibility analysis, a = 0: add intracellular metabolites (remove
% target), a = 1: add extracellular metabolites (remove target), a = 2: add
% intracellular metabolites (don't remove target)
% s: logical, Boolean variable for mass balance constraint on all
% metabolites, s = 0: equality mass balance constraint (S*V=0), s = 1:
% inequality mass balance constraints (S*V>=0)  
% samp: int, number of samples taken to estimate probability of production
% thresh: double, threshold variable used by nonlinear fitting
% n: int, number variable used by nonlinear fitting
% noise: double, noise variable used by nonlinear fitting
% runs: int, number of repeated runs used to estimate producibility metric

% OUTPUT (Saved to folder)
% results: double, vector of size (1, runs) containing the producibility
% metric results for a specific model and metabolite. One vector will be
% produced for each pair of model and metabolite analyzed.

%%
for I = 1:length(modnum) %for each model
    modelnum = modnum(I) %assign current model number
    
    % Load model
    files = dir(models_location);
    files = files(~[files.isdir]);
    fileName = files(modelnum).name;
    filePath = [models_location, filesep, fileName];
    S = load(filePath,'-mat');
    model = whos('-file',filePath);
    model = model.name;
    
    for J = 1:length(metnum) %for each metabolite
        metabnum = metnum(J) %assign current metabolite number
        
        % find target metabolites
        target = logical(false(length(S.(model).mets),1));
        mets_metnum = mets{metabnum};
        [~,target_id] = ismember(mets_metnum,S.(model).mets);
        if ismember(0,target_id) == 1
            % at least 1 of the target metabolites is not in the model
            PM = zeros(length(runs),1);
        else
            target(target_id) = true;
            
            % prepair model and extract relevant information
            [model1,inds,e_m,i_m] = prep_mod(S.(model));
            
            % set add variable
            add = zeros(length(model1.mets),1);
            if a == 0
                add = i_m; % add intracellular metabolites
                add(target) = 0; % don't add current target
            elseif a == 1
                add = e_m; % add extracellular metabolites
                add(target) = 0; % don't add current target
            elseif a == 2
                add = i_m; % add intracellular metabolites
            end
            
            % set fixed metabolites on or off
            % on
            [~,fixon_id] = ismember(fixon,model1.mets);
            add(fixon_id) = 0; % do not probabilistically add fixed metabolites
            model1.ub(inds(fixon_id)) = 1000;
            % off
            [~,fixoff_id] = ismember(fixoff,model1.mets);
            add(fixoff_id) = 0; % do not probabilistically add fixed metabolites
            
            % set ineq variable
            if s == 0
                ineq = logical(false(length(model1.mets),1)); % SV = 0 for all metabolites
            elseif s == 1
                ineq = logical(true(length(model1.mets),1)); % SV >= 0 for all metabolites
            end
            
            % Run aglorithm
            PM = zeros(length(runs),1); %initialize PM variable to store results
            for K = 1:runs
                PM(K) = calc_PM_fit_nonlin(target,model1,inds,add,ineq,samp,noise,n,thresh);
            end
        end
        % Save results
        cname = [results_location, filesep,'mod_',num2str(modnum(I)),'_met_',num2str(metnum(J)),'_PM.mat'];
        save(cname,'PM','fileName','mets_metnum')
    end
end