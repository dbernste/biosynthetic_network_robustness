#!/bin/bash -l

#$ -N PercJobs
#$ -t 1-68
#$ -l h_rt=26:00:00
#$ -j y

module load cobratoolbox/809fa7367e

echo Start time is `date`

modnum=1:10
metnum=$SGE_TASK_ID

# MATLAB function to run
matlab -nodisplay -r "addpath('/algorithm_functions'); initCobraToolbox; find_PM_mods_mets('/Results','/Models',{'10fthf_c', '2fe2s_c', '2ohph_c', '4fe4s_c', 'ala__L_c', 'amet_c', 'arg__L_c', 'asn__L_c', 'asp__L_c', 'atp_c', 'bmocogdp_c', 'btn_c', 'ca2_c', 'cl_c', 'coa_c', 'cobalt2_c', 'ctp_c', 'cu2_c', 'cys__L_c', 'datp_c', 'dctp_c', 'dgtp_c', 'dttp_c', 'fad_c', 'fe2_c', 'fe3_c', 'gln__L_c', 'glu__L_c', 'gly_c', 'gtp_c', 'h2o_c', 'his__L_c', 'ile__L_c', 'k_c', 'leu__L_c', 'lys__L_c', 'met__L_c', 'mg2_c', 'mlthf_c', 'mn2_c', 'mobd_c', 'nad_c', 'nadp_c', 'nh4_c', 'ni2_c', 'pe160_c', 'pe161_c', 'phe__L_c', 'pheme_c', 'pro__L_c', 'pydx5p_c', 'ribflv_c', 'ser__L_c', 'sheme_c', 'so4_c', 'thf_c', 'thmpp_c', 'thr__L_c', 'trp__L_c', 'tyr__L_c', 'udcpdp_c', 'utp_c', 'val__L_c', 'zn2_c', 'kdo2lipid4_e', 'murein5px4p_p', 'pe160_p', 'pe161_p'},$metnum,$modnum,{},{},'int_nt',1,50,0.3,7,0.01,10); exit"

echo End time is `date`
