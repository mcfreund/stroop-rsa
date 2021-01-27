#!/usr/bin/env bash

## get vars

glm_names=(pro_bias_acc-only_fmriprep pro_bias_acc-only_fmriprep pro_allstim_acc-only_fmriprep)
suffix_runwise=("" _runwise _runwise)
#filename="/data/nil-external/ccp/freund/ub55/in/ub55_subjects.txt"
#mapfile -t subjects < $filename
subjects=102008
runs=(1 2)
encoding_dir=(AP PA)

## directories

stimts=/data/nil-external/ccp/freund/stroop-rsa/glms/
out=/data/nil-external/ccp/freund/stroop-rsa/glms/
img=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS/
scripts=/data/nil-external/ccp/freund/stroop-rsa/code/glms/


## fit


for subject in ${subjects[@]}; do

	echo ${subject}
	
	for glm in ${glm_names[@]}; do
	
		source 3dDeconvolve_${glm}${suffix_runwise}.sh

	done

done

