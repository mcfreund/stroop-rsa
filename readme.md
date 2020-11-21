# RSA of color-word Stroop

Repository for manuscript

Submitted for publication, 2020-10-21

---

## PIPELINE SUMMMARY

This summary is broken down by section of the manuscript.
The tables here connect analysis scripts to parts of the manuscript (Figures, Tables, text sections).

### Estimating Coding Strength \beta

 * scripts in **./code/behav/**

script  | manuscript | purpose | 
--------|------------|---------|
1_define_models.R | .csv files in rsa/mods/ | Figure 1 | creates model similarity matrices
2_estimate_rsms.R | .rds files in rsa/obsv | Method::Estimating Coding Strength \beta | estimates empirical similarity matrices per ROI
3_regress_run_effects.R |.rds files in rsa/obsv with suffix tag \_residual | Method::Estimating Coding Strength \beta | prewhitens empirical similarity matrices with "run-bias" model
4_model_rsms.R | .csv files in rsa\stats | Method::Estimating Coding Strength \beta | fits RSA models (RSA regression)


### Behavioral analyses

 * scripts in **./code/behav/**
 * master file: **\_behav.rmd** (this sources all scripts in order and generates report)
 * view report: [**\_behav.html**](https://htmlpreview.github.io/?https://github.com/mcfreund/stroop-rsa/blob/ms_v1/code/behav/_behav.html)

script | manuscript | purpose | 
-------|------------|---------|
prelim_behavioral_models.R | Method::Selection of Behavioral Measures... | models behavioral data of 'primary analysis set'
microphone_comparison.R | Method::Selection of Behavioral Measures... | estimates impact of microphone change on RT measures
prelim_behavioral_models_validation_set.R | Method::Selection of Behavioral Measures... | models RT data of 'validation set'


### Group-level analyses

 * scripts in **./code/group/**
 * master file: **\_group.rmd** (this sources all scripts in order and generates report)
 * view report: [**\_group.html**](https://htmlpreview.github.io/?https://github.com/mcfreund/stroop-rsa/blob/ms_v1/code/group/_group.html)

script | manuscript | purpose | 
-------|------------|---------|
fpc_dissoc.R |  Results::Group::Dorsolateral...; Figure 2A; Tables A1--A3 | tests hypotheses regarding group-level coding dissociations
mds.R | Figure 2B | dimensionality reduction to visualize DMFC (L), V1, and SomMot--Mouth geometries
visual_sm_dissoc.R | Results::Group::Sensitivity and control analyses; Figure A1 | test SomMot--mouth and V1 coding for positive control analysis
fpc_dissoc_parcel.R | Results::Group::Sensitivity and control analyses; Figure A2 | sensitivity test for group-level analysis (parcel-level)
fpc_dissoc_altdef.R | Results::Group::Sensitivity and control analyses; Table A4 | sensitivity test for group-level analysis (alternate ROI definitions for DMFC and DLPFC)
noise_ceiling.R | Results::Group::Sensitivity and control analyses, stats inline | estimate noise ceilings per ROI
noise_ceiling_tost.R | Results::Group::Sensitivity and control analyses, stats inline | contrast noise ceiling estimates across ROIs, provide two one-sided test for equivalence


### Individual (Difference) analyses

 * scripts in **./code/indiv/**
 * master file: **\_indiv.rmd** (this sources all scripts in order and generates report)
 * view report: [**\_indiv.html**](https://htmlpreview.github.io/?https://github.com/mcfreund/stroop-rsa/blob/ms_v1/code/indiv/_indiv.html)

script | manuscript | purpose | 
-------|------------|---------|
bivar_superparcel.R | Figure 3A; Figure A3; | create bivariate scatterplots Stroop ~ coding strength
single_roi.R |  Results::Individual::Better-performing...; Table A5 | test stroop\*coding-strength relationship in each ROI\*coding-scheme separately
wn_roi_contrast.R | Results::Individual::Better-performing...; Table A6 | contrast stroop\*coding-strength relationship between coding schemes (incongr., target), within-ROI
bn_roi_contrast.R | Results::Individual::Better-performing...;  Table A7 | contrast stroop\*coding-strength relationship between ROIs, within coding schemes (incongr., target)
model_selection.R | Results::Individual::Model selection...; Figure 3B,C | elastic net model selection and validation-set prediction
bivar_allcorrs_table.R | Results::Individual::Model selection...; Table A9 | largest 20% of stroop~coding-scheme correlations observed across all superparcels


### Exploratory whole-cortex analyses

 * scripts in **./code/explor/**

script | manuscript | purpose | 
-------|------------|---------|
explor.rmd | Results::Exploratory; Figure 4, Table A10--12, Figure A6 | test group-level coding of each scheme (target, distr, incongr.) in each MMP parcel
movregs_rsa.rmd | Results::Exploratory; Table A13 | negative control: are coding schemes 'encoded' within movement regressors?


### other scripts

#### 'helpers'
scripts for sourcing, all located in ./code
 * packages.R, strings.R, funs.R, read_atlases.R (depends: write_atlases.R), read_masks.R (depends: write_masks.R)

#### 'prelim scripts'
scripts that set up initial things
 * define_behav_subjects.R (depends:), write_atlases.R (depends: ), write_masks.R (depends: )

#### misc

---

## MORE ABOUT REPO

### ./in
* contains behavioral data assembled by stroop-rsa/write_behav_stroop_rsa.R.
* no other script writes to this directory.
* underlay surfaces downloaded from BALSA: https://balsa.wustl.edu/sceneFile/show/7qP5m

### ./code
* primary analysis pipeline (bash, .R scripts)
* dynamic reports (.rmd files)
* reads from ./data, nil-bluearc (for 3D+t images), ./out
* writes to ./out
* __./behav__: scripts / reports for generating behavioral analysis (RT and accuracy)
* __./group__: scripts / reports for generating group-level analyses
* __./indiv__: scripts / reports for generating individual difference analyses (brain~behavior)

### ./out
* all output of scripts are directed to this directory
* the one exeption is output related to fMRI GLMs, which is saved within __./glms__
* __./rsa__
  * __./mods__: representational similarity models generated by define_models.R
  * __./obsv__: observed similarity matrices. saved as arrays within .rds files.
  * __./stats__: subject level fits and group statistics, saved in long-form .csvs
* __./summaries__: various summary tables (.csvs) for QC, misc analyses, and things that might be read into a manuscript file.
* __./masks__: functionally or anatomically defined brain masks, created by scripts within __./code/masks__
* __./behav__: output of scripts / reports for generating behavioral analysis (RT and accuracy)
* __./group__: output of group-level analyses
* __./indiv__: output of individual difference analyses (brain~behavior)

### ./glms
* contains AFNI GLM input (e.g., stimtime and movreg files), shell scripts, and output (e.g., .nii brick files)
* 3D+t images read from nil-bluearc
* GLMs fit on ccplinux1, and results merged with local directory via rsync
