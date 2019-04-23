#!/bin/bash -l

#$ -N PercJobs
#$ -t 1-456
#$ -l h_rt=24:00:00
#$ -j y

module load cobratoolbox/809fa7367e

echo Start time is `date`

modnum=$SGE_TASK_ID
metnum=1:9

# MATLAB function to run
matlab -nodisplay -r "addpath('/algorithm_functions'); initCobraToolbox; find_PM_mods_mets('/Results','/Models',{'cpd00029_c0', 'cpd00159_c0', 'cpd00221_c0', 'cpd00211_c0', 'cpd01711_c0', 'cpd00036_c0', 'cpd00047_c0', 'cpd00141_c0', 'cpd05178_c0'},$metnum,$modnum,{},{},'int_nt',1,50,0.3,7,0.01,10); exit"

echo End time is `date`
