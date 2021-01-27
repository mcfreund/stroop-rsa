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
	-num_stimts 34 \
	-stim_times_AM1 1 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_block_run${runs[$run_i]}.txt 'dmBLOCK(1)' -stim_label 1 block \
	-stim_times 2 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_blockONandOFF_shifted_run${runs[$run_i]}.txt 'TENTzero(0,16.8,15)' -stim_label 2 blockONandOFF \
	-stim_times 3 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_blackBLACK_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 3 blackBLACK \
	-stim_times 4 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_blackGREEN_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 4 blackGREEN \
	-stim_times 5 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_blackPINK_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 5 blackPINK \
	-stim_times 6 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_blackYELLOW_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 6 blackYELLOW \
	-stim_times 7 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_blueBLUE_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 7 blueBLUE \
	-stim_times 8 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_bluePURPLE_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 8 bluePURPLE \
	-stim_times 9 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_blueRED_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 9 blueRED \
	-stim_times 10 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_blueWHITE_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 10 blueWHITE \
	-stim_times 11 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_greenBLACK_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 11 greenBLACK \
	-stim_times 12 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_greenGREEN_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 12 greenGREEN \
	-stim_times 13 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_greenPINK_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 13 greenPINK \
	-stim_times 14 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_greenYELLOW_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 14 greenYELLOW \
	-stim_times 15 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_pinkBLACK_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 15 pinkBLACK \
	-stim_times 16 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_pinkGREEN_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 16 pinkGREEN \
	-stim_times 17 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_pinkPINK_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 17 pinkPINK \
	-stim_times 18 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_pinkYELLOW_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 18 pinkYELLOW \
	-stim_times 19 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_purpleBLUE_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 19 purpleBLUE \
	-stim_times 20 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_purplePURPLE_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 20 purplePURPLE \
	-stim_times 21 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_purpleRED_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 21 purpleRED \
	-stim_times 22 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_purpleWHITE_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 22 purpleWHITE \
	-stim_times 23 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_redBLUE_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 23 redBLUE \
	-stim_times 24 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_redPURPLE_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 24 redPURPLE \
	-stim_times 25 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_redRED_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 25 redRED \
	-stim_times 26 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_redWHITE_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 26 redWHITE \
	-stim_times 27 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_whiteBLUE_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 27 whiteBLUE \
	-stim_times 28 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_whitePURPLE_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 28 whitePURPLE \
	-stim_times 29 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_whiteRED_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 29 whiteRED \
	-stim_times 30 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_whiteWHITE_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 30 whiteWHITE \
	-stim_times 31 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_yellowBLACK_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 31 yellowBLACK \
	-stim_times 32 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_yellowGREEN_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 32 yellowGREEN \
	-stim_times 33 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_yellowPINK_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 33 yellowPINK \
	-stim_times 34 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_yellowYELLOW_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 34 yellowYELLOW \
	-ortvec ${dir_stimts}/motion_demean_proactive_run${runs[$run_i]}.1D movregs \
	-x1D ${dir_out}/X_${runs[$run_i]}.xmat.1D \
	-xjpeg ${dir_out}/X_${runs[$run_i]}.jpg \
	-nobucket

done

cd ${wd}

