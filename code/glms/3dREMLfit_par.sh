#!/usr/bin/env bash

# https://stackoverflow.com/questions/17621798/linux-process-in-background-stopped-in-jobs
# for writing stdin/stdout to file, for background running


wd=$(pwd)

function remlfit {
	
	
	## define paths and names
	dir_stimts=${stimts}${subject}/input/pro
	dir_out=${out}${subject}/results/${glm}
	name_img_run1=${img}${subject}/INPUT_DATA/Stroop/proactive/lpi_scale_blur4_tfMRI_StroopPro1_AP.nii.gz
	name_img_run2=${img}${subject}/INPUT_DATA/Stroop/proactive/lpi_scale_blur4_tfMRI_StroopPro2_PA.nii.gz

	cd ${dir_out}	

	/usr/local/pkg/afni_18/3dREMLfit \
	-matrix ${dir_out}/X.xmat.1D \
	-input ${name_img_run1} ${name_img_run1} \
	-Rvar ${dir_out}/stats_var_${subject}.nii.gz \
	-Rbuck ${dir_out}/stats_${subject}.nii.gz \
	-rwherr ${dir_out}/wherr_${subject}.nii.gz \
	-rerrts ${dir_out}/errts_${subject}.nii.gz \
	-GOFORIT \
	-fout \
	-tout \
	-nobout \
	-noFDR \
	-verb

}


logpath=${out}${subject}/results/${glm}

remlfit < /dev/null > ${logpath}/runtime.log 2>&1 &


cd ${wd}

