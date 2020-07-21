## about ----
## 
## reads .nii files that are masks for 'user-specified' ROIs.
## these ROIs were specified on the basis of neuroanatomy (e.g., VWFA) and on function (results from RSA).
## 
## mike freund, 2019-01-10

## TODO
## create masks.key file with definitions (same structure as atlas.key)
## embed reading in of masks in loop

## initialize ----

# mask.names <- c(
#   "vwfa"
# )
# 
# masks <- setNames(vector("list", length(mask.names)), mask.names)
# rm(mask.names)
# 
# 
# ## vwfa ----
# 
# masks$vwfa <- oro.nifti::readNIfTI(here::here("out", "masks", "vwfa.nii.gz"), reorient = FALSE)
# 

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
fnames <- list.files(here::here("out", "masks"), pattern = "\\.nii")
=======
fnames <- list.files(here::here("out", "masks"), pattern = "\\.nii\\.gz$")
>>>>>>> 4de56ea180c2a28fca16b7a04dea418f9609edc0
=======
fnames <- list.files(here::here("out", "masks"), pattern = "\\.nii\\.gz$")
>>>>>>> 4de56ea180c2a28fca16b7a04dea418f9609edc0
=======
fnames <- list.files(here::here("out", "masks"), pattern = "\\.nii\\.gz$")
>>>>>>> 4de56ea180c2a28fca16b7a04dea418f9609edc0

masks <- lapply(
  fnames, 
  function(x)
  oro.nifti::readNIfTI(here::here("out", "masks", x), reorient = FALSE)
)

names(masks) <- gsub("\\.nii.*$", "", fnames)
