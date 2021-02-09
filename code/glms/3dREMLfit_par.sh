#!/usr/bin/env bash

# https://stackoverflow.com/questions/17621798/linux-process-in-background-stopped-in-jobs
# for writing stdin/stdout to file, for background running


wd=$(pwd)

function remlfit {
	
	cd ${dir_out}	

	/usr/local/pkg/afni_18/3dREMLfit \
	-matrix ${dir_out}/X.xmat.1D \
	-input ${name_img} \
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


logpath=${out}${subject}/results/${glm}

#echo ${logpath}

remlfit < /dev/null > ${logpath}/runtime.log 2>&1 &


cd ${wd}

