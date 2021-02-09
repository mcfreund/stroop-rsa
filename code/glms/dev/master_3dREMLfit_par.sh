#!/usr/bin/env bash


## get vars

glm_names=(pro_bias_acc-only_fmriprep pro_bias_acc-only_fmriprep pro_allstim_acc-only_fmriprep)
suffices=("" _runwise _runwise)


subjects=132017
#filename="/data/nil-external/ccp/freund/stroop-rsa/in/subjects.txt"
#mapfile -t subjects < $filename

runs=(1 2)
encoding_dir=(AP PA)
#hemis=(L R)

## directories

stimts=/data/nil-external/ccp/freund/stroop-rsa/glms/
out=/data/nil-external/ccp/freund/stroop-rsa/glms/
img=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS/
scripts=/data/nil-external/ccp/freund/stroop-rsa/code/glms/


## fit


for subject in ${subjects[@]}; do

	echo ${subject}
	
	for glm_i in ${!glm_names[@]}; do
	
		glm=${glm_names[$glm_i]}
		suffix=${suffices[$glm_i]}
		
		echo ${glm}${suffix}

		source 3dREMLfit${suffix}_par.sh

	done

	wait

done

