#!/usr/bin/env bash

wd=$(pwd)

for run_i in ${!runs[@]}; do

	## define paths and names
	dir_stimts=${stimts}${subject}/input/pro
	dir_out=${out}${subject}/results/${glm}_run${runs[$run_i]}
	name_img=${img}${subject}/INPUT_DATA/Stroop/proactive/lpi_scale_blur4_tfMRI_StroopPro${runs[$run_i]}_${encoding_dir[$run_i]}.nii.gz

	## make result dir
	mkdir -p ${dir_out}
	cd ${dir_out}

	## build xmat
	/usr/local/pkg/afni_18/3dDeconvolve \
	-local_times \
	-x1D_stop \
	-allzero_OK \
	-input ${name_img} \
	-polort A \
	-float \
	-censor ${dir_stimts}/movregs_FD_mask_run${runs[$run_i]}.txt \
	-num_stimts 21 \
	-stim_times 1 ${dir_stimts}/${subject}_pro_blueBLUE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)'  -stim_label 1 blueBLUE \
	-stim_times 2 ${dir_stimts}/${subject}_pro_purpleBLUE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)'  -stim_label 2 purpleBLUE \
	-stim_times 3 ${dir_stimts}/${subject}_pro_redBLUE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)'  -stim_label 3 redBLUE \
	-stim_times 4 ${dir_stimts}/${subject}_pro_whiteBLUE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)'  -stim_label 4 whiteBLUE \
	-stim_times 5 ${dir_stimts}/${subject}_pro_bluePURPLE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)'  -stim_label 5 bluePURPLE \
	-stim_times 6 ${dir_stimts}/${subject}_pro_purplePURPLE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)'  -stim_label 6 purplePURPLE \
	-stim_times 7 ${dir_stimts}/${subject}_pro_redPURPLE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)'  -stim_label 7 redPURPLE \
	-stim_times 8 ${dir_stimts}/${subject}_pro_whitePURPLE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)'  -stim_label 8 whitePURPLE \
	-stim_times 9 ${dir_stimts}/${subject}_pro_blueRED_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)'  -stim_label 9 blueRED \
	-stim_times 10 ${dir_stimts}/${subject}_pro_purpleRED_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)'  -stim_label 10 purpleRED \
	-stim_times 11 ${dir_stimts}/${subject}_pro_redRED_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)'  -stim_label 11 redRED \
	-stim_times 12 ${dir_stimts}/${subject}_pro_whiteRED_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)'  -stim_label 12 whiteRED \
	-stim_times 13 ${dir_stimts}/${subject}_pro_blueWHITE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)'  -stim_label 13 blueWHITE \
	-stim_times 14 ${dir_stimts}/${subject}_pro_purpleWHITE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)'  -stim_label 14 purpleWHITE \
	-stim_times 15 ${dir_stimts}/${subject}_pro_redWHITE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)'  -stim_label 15 redWHITE \
	-stim_times 16 ${dir_stimts}/${subject}_pro_whiteWHITE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)'  -stim_label 16 whiteWHITE \
	-stim_times 17 ${dir_stimts}/${subject}_pro_pc50_incon_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)'  -stim_label 17 pc50_i \
	-stim_times 18 ${dir_stimts}/${subject}_pro_pc50_congr_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)'  -stim_label 18 pc50_c \
	-stim_times 19 ${dir_stimts}/${subject}_pro_nuisance_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)'  -stim_label 19 nuisance \
	-stim_times_AM1 20 ${dir_stimts}/${subject}_pro_sustained_run${runs[$run_i]}.txt 'dmBLOCK(1)'  -stim_label 20 sustained \
	-stim_times 21 ${dir_stimts}/${subject}_pro_transient_run${runs[$run_i]}.txt 'TENTzero(0,16.8,8)'  -stim_label 21 transient \
	-ortvec ${dir_stimts}/Movement_Regressors_StroopPro${runs[$run_i]}_${encoding_dir[$run_i]}.1D movregs \
	-x1D ${dir_out}/X_run${runs[$run_i]}.xmat.1D \
	-xjpeg ${dir_out}/X_run${runs[$run_i]}.jpg \
	-nobucket
	
	

done

cd ${wd}

