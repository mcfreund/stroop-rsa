#!/usr/bin/env bash


## get vars

glm_names=(pro_bias_acc-only_fmriprep pro_bias_acc-only_fmriprep_run1 pro_bias_acc-only_fmriprep_run2 pro_allstim_acc-only_fmriprep_run1 pro_allstim_acc-only_fmriprep_run2)

subjects=132017
#filename="/data/nil-external/ccp/freund/stroop-rsa/in/subjects.txt"
#mapfile -t subjects < $filename

## directories

stimts=/data/nil-external/ccp/freund/stroop-rsa/glms/
out=/data/nil-external/ccp/freund/stroop-rsa/glms/
img=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS/
scripts=/data/nil-external/ccp/freund/stroop-rsa/code/glms/



## fit


for subject in ${subjects[@]}; do
	#subject=132017
	
	echo ${subject}
	
	for glm_i in ${!glm_names[@]}; do	
		#glm_i=1
		
		glm=${glm_names[$glm_i]}
		
		echo ${glm}
		
		## inputs
		
		dir_stimts=${stimts}${subject}/input/pro
		dir_out=${out}${subject}/results/${glm}
		
		if [[ "$glm" == *"_run"* ]]; then
					
			if [[ "$glm" == *"_run1" ]]; then
			
				suffix=_run1
				encoding_dir=1_AP
				xmat=X_run1.xmat.1D
				
			else
				
				suffix=_run2
				encoding_dir=2_PA
				xmat=X_run2.xmat.1D
				
			fi
			
			name_img=${img}${subject}/INPUT_DATA/Stroop/proactive/lpi_scale_blur4_tfMRI_StroopPro${encoding_dir}.nii.gz
			
		else
		
			suffix=""
			
			name_img_run1=${img}${subject}/INPUT_DATA/Stroop/proactive/lpi_scale_blur4_tfMRI_StroopPro1_AP.nii.gz
			name_img_run2=${img}${subject}/INPUT_DATA/Stroop/proactive/lpi_scale_blur4_tfMRI_StroopPro2_PA.nii.gz
			name_img=${name_img_run1}" "${name_img_run2}
			xmat=X.xmat.1D
		
		fi
		
		## run
		
		source 3dREMLfit_par.sh

	done

	wait

done
