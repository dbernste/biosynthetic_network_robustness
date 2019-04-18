% David Bernstein

% This recursive function is used to be enumerate and sum binomial 
% coefficient in the theoretical calculation of the output probability

function [sumr] = f_r_calc_prob(s,i,sumr,jvect,S,Pin)
% Input
% s - (start) 1 + value of previous loop
% i - number of sets being considered as a component of the larger inclusion
% exclusion sum, this defines the number of "for loops" in the recursive sum
% sumr - value of ongoing recursive sum
% jvect - vector to keep track of the indicies of the current sets of interest
% S - A matrix of all of the sets of metabolite that allow for the
% production of the metabolite of interest. Each column corresponds to a
% particular metabolite and each row corresponds to a particular
% minimal precursor set, with 1 meaning that metabolite is included in that set.
% The rows can be though of as or relationships while the columns can be
% thought of as and relationships pertaining to the necesity of these
% metabolites to synthesize a particular metabolite
% Pin - a vector of the input probability of each metabolite it has length 
% of size(S,2)

if i >= 1 %condition for exiting recursive loop (count down from total number of sets being considered)
    for j = s:size(S,1)
        jvect(1,i) = j; %keep track of indicies of the for loops
        sumr = f_r_calc_prob(j+1,i-1,sumr,jvect,S,Pin); %Recursivelly call function
    end
else
    % Caclulate value to be summed
    prod = 1; %Define a product value
    % for all elements in the union of sets defined by jvect multiply prod
    % by their probability
    for k = 1:size(S,2)
        done = 0;
        for l = 1:length(jvect)
            if done == 0
                if S(jvect(l),k) == 1
                    prod = prod * Pin(k);
                    done = 1;
                end
            end
        end
    end
    sumr = sumr + prod;
end
