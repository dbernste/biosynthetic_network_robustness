#!/bin/bash -l

#$ -N PercJobs
#$ -t 1-52
#$ -l h_rt=03:00:00
#$ -j y

module load cobratoolbox/809fa7367e

echo Start time is `date`

metnum=$SGE_TASK_ID

# MATLAB function to run
matlab -nodisplay -r "addpath('/algorithm_functions'); initCobraToolbox; find_PM_mods_mets('/Results','/Model',{'13dpg[c]','2pg[c]','3pg[c]','6pgc[c]','6pgl[c]','ac[c]','acald[c]','accoa[c]','acon-C[c]','actp[c]','adp[c]','akg[c]','amp[c]','atp[c]','cit[c]','co2[c]','coa[c]','dhap[c]','e4p[c]','etoh[c]','f6p[c]','fdp[c]','for[c]','fum[c]','g3p[c]','g6p[c]','gln-L[c]','glu-L[c]','glx[c]','h2o[c]','h[c]','icit[c]','lac-D[c]','mal-L[c]','nad[c]','nadh[c]','nadp[c]','nadph[c]','nh4[c]','o2[c]','oaa[c]','pep[c]','pi[c]','pyr[c]','q8[c]','q8h2[c]','r5p[c]','ru5p-D[c]','s7p[c]','succ[c]','succoa[c]','xu5p-D[c]'},$metnum,1,{},{},'int_nt',1,50,0.3,7,0.01,10); exit"

echo End time is `date`
