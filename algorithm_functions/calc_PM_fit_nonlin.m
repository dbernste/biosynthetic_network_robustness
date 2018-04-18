% David Bernstein

function [PM] = calc_PM_fit_nonlin(target,model1,inds,add,ineq,samp,noise,n,thresh)
% Function to caclulate the producibility metric of a given metabolic model 
% and target metabolite set using a nonlinear fitting algorithm

% INPUT
% target: logical, vector of length of metabolites with 0 for non-target 
% metabolites and 1 for target metaboltes
% model1: structure, COBRA model after being run through prep_mod, exchange 
% reactions off, maintenance reaction off, additional exchange for all 
% metabolites added with indices recorded
% inds: int, vector of length of metabolites with the indices of metabolite 
% exchange reactions
% add: logical, vector of length of metabolites with 0 for metabolites that 
% will not be randomly added and 1 for metabolites that will be randomly 
% added during the producibility analysis process
% ineq: logical, vector of length of metabolites with 0 for equality mass 
% balance (S*V>0) and 1 for inequality mass balance (S*V>=0)
% samp: int, number of samples taken to estimate probability of production
% thresh: double, threshold variable used by nonlinear fitting
% n: int, number variable used by nonlinear fitting
% noise: double, noise variable used by nonlinear fitting

% OUTPUT
% PM: double, output of calc_PM_fit_Nonlin, producibility metric

% HARD CODED VARIABLES
% limit: int, number of runs before algorithm cuts off and returns NaN for 
% producibility metric

%%
% Initialize by finding the p_out values for p_in = 0 (no mets added) and 1
% (all mets added). For a well defined metabolic network p_out(p_in=0)=0, 
% p_out(p_in=1)=1.

bern(add==1) = 0;
PROB0 = prob(target,model1,inds,bern,ineq,samp); %p_out(p_in=0)
bern(add==1) = 1;
PROB1 = prob(target,model1,inds,bern,ineq,samp); %p_out(p_in=1)

if PROB0 == 1 %if target metabolite is producible from no inputs
    PM = 1;
elseif PROB1 == 0 %if target metabolite is not producible from all inputs
    PM = 0;
else
    % Use nonlinear fitting algorithm to find PM
    
    % Initialize curve
    % Endpoints
    X = [0; 1];
    Y = [PROB0; PROB1];
    
    % Set parameters for nonlinear fitting
    % algorithm attempts to estimate the probability of adding metabolties 
    % (p_in) at which the probability of producing the target metabolite 
    % (p_out) = 0.5, (p_0.5)
    
    % Start in the middle
    Estimate_off = 0.5; % offset point at which to calculate the new probability value in order to improve the estimate
    Estimate_list = Estimate_off; % list of estimated points from nonlinear fitting
    count = 0; % count to keep track of algorithm runs
    Done = 0; % indicator to determine when algorithm converges
    cfit = [0.5 20]; % initial sigmoid parameters
    
    limit = 1000; % Set upper limit on algorithm in case convergence fails
    failed_limit = 0;
    
    while Done == 0
        count = count + 1;
        
        bern(add==1) = Estimate_off; % update p_in value
        PROB = prob(target,model1,inds,bern,ineq,samp); % calculate p_out value
        % update points to plot
        X = [X;Estimate_off];
        Y = [Y;PROB];
        
        %fit sigmoid function using MATLAB lsqnonlin function
        fun = @(c) 1./(1+exp((c(1)-X)*c(2)))-Y;
        c0 = cfit;
        opts1 =  optimset('display','off');
        cfit = lsqnonlin(fun,c0,[],[],opts1);
        
%        % display fitting figure for debugging
%         figure(1)
%         hold off
%        plot(X,Y,'mo')
%         xx = [0:0.00001:1];
%         yy = 1./(1+exp((cfit(1)-xx)*cfit(2)));
%         plot(xx,yy,'r-')
%        pause(eps)
        
        Est = cfit(1); % New estimate of p_0.5 value
        % check to make sure Est is between 0 and 1;
        if Est < 0
            Est = 0;
        elseif Est > 1
            Est = 1;
        end
        Estimate_list = [Estimate_list,Est];
        
        % check for convergence
        if length(Estimate_list)>=n
            Done = 1;
            % check the distance between all pairwise combinations of the
            % last n estimates to see if any is greater than thresh
            for I = 1:n
                for J = I+1:n
                    if Done == 1
                        if abs(Estimate_list(end-n+I) - Estimate_list(end-n+J)) > thresh
                            Done = 0; %if one of the last estimates is not within the threshold, then the algorithm has not converged
                        end
                    end
                end
            end
        end
        
        % find new value shifted slightly from the estimated value in order
        % to refit the nonlinear function with a new point near the
        % estimate the shifted value will alternate between being chosen
        % above and below the predicted value and will be offset based on
        % the noise parameter
        shift = rand(1)*noise;
        if mod(count,2) == 1 %odd count 
            %shift < 0
            shift_scaled = Est*shift;
            Estimate_off = Est - shift_scaled;
        else %even count
            %shift > 0
            shift_scaled = (1-Est)*shift;
            Estimate_off = Est + shift_scaled;
        end
        
        if count > limit %Number of runs is greater than limit, algorithm failed to converge
            Done = 1;
            failed_limit = 1;
        end
    end
    
    if failed_limit == 0
        PM = 1-Est; % PM = 1 - P_0.5
    else
        PM = NaN;
    end
    
end
end

