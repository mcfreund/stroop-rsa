#!/usr/bin/env bash

wd=$(pwd)

for run_i in ${!runs[@]}; do

	## define paths and names
	dir_stimts=${stimts}${subject}/input/pro
	dir_out=${out}${subject}/results/${glm_names[$glm_i]}_run${runs[$run_i]}
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
	-num_stimts 34 \
	-stim_times 1 ${dir_stimts}/${subject}_pro_blackBLACK_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 1 blackBLACK \
	-stim_times 2 ${dir_stimts}/${subject}_pro_blackGREEN_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 2 blackGREEN \
	-stim_times 3 ${dir_stimts}/${subject}_pro_blackPINK_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 3 blackPINK \
	-stim_times 4 ${dir_stimts}/${subject}_pro_blackYELLOW_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 4 blackYELLOW \
	-stim_times 5 ${dir_stimts}/${subject}_pro_blueBLUE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 5 blueBLUE \
	-stim_times 6 ${dir_stimts}/${subject}_pro_bluePURPLE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 6 bluePURPLE \
	-stim_times 7 ${dir_stimts}/${subject}_pro_blueRED_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 7 blueRED \
	-stim_times 8 ${dir_stimts}/${subject}_pro_blueWHITE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 8 blueWHITE \
	-stim_times 9 ${dir_stimts}/${subject}_pro_greenBLACK_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 9 greenBLACK \
	-stim_times 10 ${dir_stimts}/${subject}_pro_greenGREEN_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 10 greenGREEN \
	-stim_times 11 ${dir_stimts}/${subject}_pro_greenPINK_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 11 greenPINK \
	-stim_times 12 ${dir_stimts}/${subject}_pro_greenYELLOW_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 12 greenYELLOW \
	-stim_times 13 ${dir_stimts}/${subject}_pro_pinkBLACK_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 13 pinkBLACK \
	-stim_times 14 ${dir_stimts}/${subject}_pro_pinkGREEN_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 14 pinkGREEN \
	-stim_times 15 ${dir_stimts}/${subject}_pro_pinkPINK_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 15 pinkPINK \
	-stim_times 16 ${dir_stimts}/${subject}_pro_pinkYELLOW_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 16 pinkYELLOW \
	-stim_times 17 ${dir_stimts}/${subject}_pro_purpleBLUE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 17 purpleBLUE \
	-stim_times 18 ${dir_stimts}/${subject}_pro_purplePURPLE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 18 purplePURPLE \
	-stim_times 19 ${dir_stimts}/${subject}_pro_purpleRED_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 19 purpleRED \
	-stim_times 20 ${dir_stimts}/${subject}_pro_purpleWHITE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 20 purpleWHITE \
	-stim_times 21 ${dir_stimts}/${subject}_pro_redBLUE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 21 redBLUE \
	-stim_times 22 ${dir_stimts}/${subject}_pro_redPURPLE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 22 redPURPLE \
	-stim_times 23 ${dir_stimts}/${subject}_pro_redRED_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 23 redRED \
	-stim_times 24 ${dir_stimts}/${subject}_pro_redWHITE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 24 redWHITE \
	-stim_times 25 ${dir_stimts}/${subject}_pro_whiteBLUE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 25 whiteBLUE \
	-stim_times 26 ${dir_stimts}/${subject}_pro_whitePURPLE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 26 whitePURPLE \
	-stim_times 27 ${dir_stimts}/${subject}_pro_whiteRED_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 27 whiteRED \
	-stim_times 28 ${dir_stimts}/${subject}_pro_whiteWHITE_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 28 whiteWHITE \
	-stim_times 29 ${dir_stimts}/${subject}_pro_yellowBLACK_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 29 yellowBLACK \
	-stim_times 30 ${dir_stimts}/${subject}_pro_yellowGREEN_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 30 yellowGREEN \
	-stim_times 31 ${dir_stimts}/${subject}_pro_yellowPINK_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 31 yellowPINK \
	-stim_times 32 ${dir_stimts}/${subject}_pro_yellowYELLOW_run${runs[$run_i]}_acc-only.txt 'BLOCK(1,1)' -stim_label 32 yellowYELLOW \
	-stim_times_AM1 33 ${dir_stimts}/${subject}_pro_sustained_run${runs[$run_i]}.txt 'dmBLOCK(1)'  -stim_label 33 sustained \
	-stim_times 34 ${dir_stimts}/${subject}_pro_transient_run${runs[$run_i]}.txt 'TENTzero(0,16.8,8)'  -stim_label 34 transient \
	-ortvec ${dir_stimts}/Movement_Regressors_StroopPro${runs[$run_i]}_${encoding_dir[$run_i]}.1D movregs \
	-x1D ${dir_out}/X_${runs[$run_i]}.xmat.1D \
	-xjpeg ${dir_out}/X_${runs[$run_i]}.jpg \
	-nobucket

done

cd ${wd}

