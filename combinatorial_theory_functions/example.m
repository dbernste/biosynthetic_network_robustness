% David Bernstein
% Example: Run combinatorial theory code

%% Define minimal precursor set structure.
% The set of minimal precursor sets can be described by a matrix: S
% The columns of S are different metabolites. 1 indicates that this 
% metabolite is in the minimal precursor set of the row, 0 indicates that 
% it is not. The Matrix contains only metabolites that are involved in at 
% least 1 minimal precursor set and could be extended with all 0's to 
% contain all metabolites, but this is not necessary.

% single metabolite minimal precursor set
% S = [1];

% single minimal precursor set of size n
% n = 10;
% S = ones(1,n);

% m minimal precursor sets of size 1
% m = 10;
% S = eye(m);

% arbitrary minimal precursor set structure
S = [1 0 0 1 0; 0 1 1 0 0; 1 0 1 0 1];

%% Calculate Pout
% define input probability vector (identical for all metabolties in the minimal precursor sets)
Pin = 0.5 *ones(size(S,2));

% calculate output probability
Pout = f_calc_prob(S,Pin);

%% Plot Producibility Curve
Pin_vect = [0:0.01:1];

Pout_vect = zeros(size(Pin_vect));
for I = 1:length(Pin_vect)
    Pin_I = Pin_vect(I)*ones(size(S,2));
    Pout_vect(I) = f_calc_prob(S,Pin_I);
end

figure()
plot(Pin_vect,Pout_vect,'k.')
axis([0 1 0 1])

