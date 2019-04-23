#!/bin/bash -l

#$ -N PercJobs
#$ -t 1-97
#$ -l h_rt=12:00:00
#$ -j y

module load cobratoolbox/809fa7367e

echo Start time is `date`

metnum=$SGE_TASK_ID
modnum=[1:3]

# MATLAB function to run
matlab -nodisplay -r "addpath('/algorithm_functions'); initCobraToolbox; find_PM_mods_mets('/Results_carb','/Models',{'cpd00118_c0','cpd00156_c0','cpd00322_c0','cpd00016_c0','cpd00132_c0','cpd00129_c0','cpd00051_c0','cpd00119_c0','cpd02229_c0','cpd00056_c0','cpd00039_c0','cpd00066_c0','cpd00069_c0','cpd00015_c0','cpd00220_c0','cpd00054_c0','cpd00357_c0','cpd00038_c0','cpd00241_c0','cpd00161_c0','cpd00107_c0','cpd00084_c0','cpd00060_c0','cpd00017_c0','cpd11493_c0','cpd00356_c0','cpd00115_c0','cpd00035_c0','cpd00010_c0','cpd00041_c0','cpd00001_c0','cpd00033_c0','cpd00052_c0','cpd00002_c0','cpd00062_c0','cpd00023_c0','cpd00053_c0','cpd00003_c0','cpd00006_c0','cpd00201_c0','cpd00087_c0','cpd00149_c0','cpd00166_c0','cpd15775_c0','cpd15776_c0','cpd15777_c0','cpd15668_c0','cpd15667_c0','cpd15669_c0','cpd15757_c0','cpd15758_c0','cpd15759_c0','cpd00557_c0','cpd11459_c0','cpd15432_c0','cpd00254_c0','cpd00034_c0','cpd00030_c0','cpd00058_c0','cpd00205_c0','cpd00063_c0','cpd15533_c0','cpd15695_c0','cpd15696_c0','cpd00065_c0','cpd00345_c0','cpd00042_c0','cpd15793_c0','cpd15794_c0','cpd15795_c0','cpd15540_c0','cpd15722_c0','cpd15723_c0','cpd00264_c0','cpd00099_c0','cpd00028_c0','cpd15767_c0','cpd15766_c0','cpd15768_c0','cpd15749_c0','cpd15748_c0','cpd15750_c0','cpd10516_c0','cpd10515_c0','cpd15665_c0','cpd15500_c0','cpd15560_c0','cpd15352_c0','cpd00029_c0','cpd00159_c0','cpd00221_c0','cpd00211_c0','cpd01711_c0','cpd00036_c0','cpd00047_c0','cpd00141_c0','cpd05178_c0'},$metnum,$modnum,{'cpd00027_c0'},{},'int_nt',1,50,0.3,7,0.01,10); exit"

echo End time is `date`
