#!/usr/bin/env bash


## get vars

glm_names=(pro_bias_acc-only_fmriprep pro_bias_acc-only_fmriprep_run1 pro_bias_acc-only_fmriprep_run2 pro_allstim_acc-only_fmriprep_run1 pro_allstim_acc-only_fmriprep_run2)
#glm_names=(pro_bias_acc-only_fmriprep)


filename="/data/nil-external/ccp/freund/stroop-rsa/in/subjects.txt"
mapfile -t subjects < $filename
#subjects=132017


## directories

stimts=/data/nil-external/ccp/freund/stroop-rsa/glms/
out=/data/nil-external/ccp/freund/stroop-rsa/glms/
img=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS/
scripts=/data/nil-external/ccp/freund/stroop-rsa/code/glms/

## fit


for subject in ${subjects[@]}; do
	#subject=132017
	
	echo ${subject}
	
	source 3dREMLfit_parallel.sh
	
	wait

done
