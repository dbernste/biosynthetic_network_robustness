% David Bernstein

% example_calculate_PM

% This example code goes through the calculation of the producibility
% metric (PM) from Bernstein et al.

% The code contains a step by step procedure for calculating an example PM
% for 2 metabolites on 2 different metabolic network models.

% Commented code is also provided for calculating the PM for 88 different
% essentail biomass components across 456 metabolic network models from 
% the oral microbiome. However, we do not reccomend running this code
% without parallelizing.

% 1) 
% Initialize the COBRA toolbox and set solver
% The COBRA toolbox is a popular compendium of metabolic network analysis
% methods. It contains the functions for running flux balance analysis
% (FBA) which are used to calculate the PM. Please visit the COBRA toolbox
% webpage for additional information on the toolbox including installation
% (https://opencobra.github.io/cobratoolbox/stable/)
initCobraToolbox
changeCobraSolver('glpk')

% 2) 
% Add to the MATLAB path the folder containing the algorithm funcitons
addpath('./algorithm_functions')

% 3)
% Set Parameters
% Define the parameters for the PM calculation: all parameters are
% described in the comments for the find_PM_mods_mets command which is also
% coppied here:
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
% a: char, indicator variable for metabolites randomly added during
% producibility analysis, a = int_nt: add intracellular metabolites (remove
% target), a = ext_nt: add extracellular metabolites (remove target), a = int: add
% intracellular metabolites (don't remove target), a = all_nt: add all
% metabolites (remove target)
% s: logical, Boolean variable for mass balance constraint on all
% metabolites, s = 0: equality mass balance constraint (S*V=0), s = 1:
% inequality mass balance constraints (S*V>=0)  
% samp: int, number of samples taken to estimate probability of production
% thresh: double, threshold variable used by nonlinear fitting
% n: int, number variable used by nonlinear fitting
% noise: double, noise variable used by nonlinear fitting
% runs: int, number of repeated runs used to estimate producibility metric
results_location = './example_results';
models_location = './example_models';
mets = {'cpd00118_c0','cpd00156_c0'}; % Putrescine, L-Valine
metnum = 1:2;
modnum = 1:2;
fixon = {};
fixoff = {};
a = 'int_nt';
s = 1;
samp = 50;
noise = 0.3;
n = 7;
thresh = 0.01;
runs = 10;

% 4)
% Calculate the PM
find_PM_mods_mets(results_location,models_location,mets,metnum,modnum,fixon,fixoff,a,s,samp,noise,n,thresh,runs)

% 5)
% Commented code for calculating PM for 88 essentail biomass components
% across 456 oral microbiome metabolic networks. 
% Running this code will take a long time
% The analysis for Bernstein et al. was conducted on the BU shared
% computing cluster where individual models and metabolites were
% parallelized to improve runtime. To parallelize runs, the metnum and
% modnum variables can be set to variables that change across computers.

% models_location = './metabolic_networks_mat';
% mets = {'cpd00118_c0','cpd00156_c0','cpd00322_c0','cpd00016_c0','cpd00132_c0','cpd00129_c0','cpd00051_c0','cpd00119_c0','cpd02229_c0','cpd00056_c0','cpd00039_c0','cpd00066_c0','cpd00069_c0','cpd00015_c0','cpd00220_c0','cpd00054_c0','cpd00357_c0','cpd00038_c0','cpd00241_c0','cpd00161_c0','cpd00107_c0','cpd00084_c0','cpd00060_c0','cpd00017_c0','cpd11493_c0','cpd00356_c0','cpd00115_c0','cpd00035_c0','cpd00010_c0','cpd00041_c0','cpd00001_c0','cpd00033_c0','cpd00052_c0','cpd00002_c0','cpd00062_c0','cpd00023_c0','cpd00053_c0','cpd00003_c0','cpd00006_c0','cpd00201_c0','cpd00087_c0','cpd00149_c0','cpd00166_c0','cpd15775_c0','cpd15776_c0','cpd15777_c0','cpd15668_c0','cpd15667_c0','cpd15669_c0','cpd15757_c0','cpd15758_c0','cpd15759_c0','cpd00557_c0','cpd11459_c0','cpd15432_c0','cpd00254_c0','cpd00034_c0','cpd00030_c0','cpd00058_c0','cpd00205_c0','cpd00063_c0','cpd15533_c0','cpd15695_c0','cpd15696_c0','cpd00065_c0','cpd00345_c0','cpd00042_c0','cpd15793_c0','cpd15794_c0','cpd15795_c0','cpd15540_c0','cpd15722_c0','cpd15723_c0','cpd00264_c0','cpd00099_c0','cpd00028_c0','cpd15767_c0','cpd15766_c0','cpd15768_c0','cpd15749_c0','cpd15748_c0','cpd15750_c0','cpd10516_c0','cpd10515_c0','cpd15665_c0','cpd15500_c0','cpd15560_c0','cpd15352_c0'};
% metnum = 1:88;
% modnum = 1:456;
% find_PM_mods_mets(results_location,models_location,mets,metnum,modnum,fixon,fixoff,a,s,samp,noise,n,thresh,runs)


