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


## glm function

function remlfit {
	
	cd ${dir_out}	

	/usr/local/pkg/afni_18/3dREMLfit \
	-matrix ${dir_out}/${xmat} \
	-input "${name_img}" \
	-Rvar ${dir_out}/stats_var_${subject}${suffix}.nii.gz \
	-Rbuck ${dir_out}/stats_${subject}${suffix}.nii.gz \
	-rwherr ${dir_out}/wherr_${subject}${suffix}.nii.gz \
	-rerrts ${dir_out}/errts_${subject}${suffix}.nii.gz \
	-GOFORIT \
	-fout \
	-tout \
	-nobout \
	-noFDR \
	-verb
	
}



## fit


for subject in ${subjects[@]}; do
	#subject=132017
	
	echo ${subject}
	
	for glm_i in ${!glm_names[@]}; do	
		#glm_i=1
		
		cd $scripts  ## write any output to this dir...
		
		glm=${glm_names[$glm_i]}
		
		echo ${glm}
		
		## inputs
		
		dir_stimts=${stimts}${subject}/input/pro
		dir_out=${out}${subject}/results/${glm}  ## will change to this dir when running remlfit
		
		if [[ "$glm" == *"_run"* ]]; then
					
			if [[ "$glm" == *"_run1" ]]; then
			
				suffix=_run1
				encoding_dir=1_AP
				xmat=X_1.xmat.1D
				
			else
				
				suffix=_run2
				encoding_dir=2_PA
				xmat=X_2.xmat.1D
				
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
		
		logpath=${out}${subject}/results/${glm}
		remlfit < /dev/null > ${logpath}/runtime.log 2>&1 &


	done

	wait

done
