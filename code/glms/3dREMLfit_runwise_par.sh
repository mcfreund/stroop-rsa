#!/usr/bin/env bash

# https://stackoverflow.com/questions/17621798/linux-process-in-background-stopped-in-jobs
# for writing stdin/stdout to file, for background running


wd=$(pwd)

function remlfit {
	
	local run_i=$1
	
	## define paths and names
	dir_stimts=${stimts}${subject}/input/pro
	dir_out=${out}${subject}/results/${glm}
	name_img=${img}${subject}/INPUT_DATA/${task}/${sessions[$session_i]}/lpi_scale_tfMRI_${task}${sess}${runs[$run_i]}_${encoding_dir[$run_i]}_${hemi}.nii.gz

	cd ${dir_out}	

	/usr/local/pkg/afni_18/3dREMLfit \
	-matrix ${dir_out}/X.xmat.1D \
	-input ${name_img} \
	-Rvar ${dir_out}/stats_var_${subject}_run${runs[$run_i]}.nii.gz \
	-Rbuck ${dir_out}/stats_${subject}_run${runs[$run_i]}.nii.gz \
	-rwherr ${dir_out}/wherr_${subject}_run${runs[$run_i]}.nii.gz \
	-rerrts ${dir_out}/errts_${subject}_run${runs[$run_i]}.nii.gz \
	-GOFORIT \
	-fout \
	-tout \
	-nobout \
	-noFDR \
	-verb

}


for run_i in ${!runs[@]}; do

	logpath=${out}${subject}/results/${glm}

	#echo ${logpath}

	remlfit ${run_i} < /dev/null > ${logpath}/runtime.log 2>&1 &

done


cd ${wd}

