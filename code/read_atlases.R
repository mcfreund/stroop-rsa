## about ----
##
## reads atlases and atlas info.
## Multi-modal parcellation: Glasser (2018)
## Gordon: Gordon (2014)
## 
## depends on successful completion of  _write_atlases.R script.

## atlases

atlas <- setNames(vector("list", 2), c("mmp", "gordon"))

atlas$mmp <- oro.nifti::readNIfTI(here::here("out", "atlases", "mmp.nii.gz"), reorient = FALSE)

atlas.key$mmp[atlas.key$mmp$hemi == "L", "num.roi"] <- atlas.key$mmp$num.roi[atlas.key$mmp$hemi == "L"] + 180
atlas.key$mmp[atlas.key$mmp$hemi == "R", "num.roi"] <- atlas.key$mmp$num.roi[atlas.key$mmp$hemi == "R"] - 180

atlas$gordon <- oro.nifti::readNIfTI(here::here("out", "atlases", "gordon.nii.gz"), reorient = FALSE)


## keys

atlas.key <- setNames(vector("list", 2),  c("mmp", "gordon"))

atlas.key$mmp <- read.csv(here::here("out", "atlases", "mmp.csv"), stringsAsFactors = FALSE)
atlas.key$gordon <- read.csv(here::here("out", "atlases", "gordon.csv"), stringsAsFactors = FALSE)

atlas.key$mmp$X <- NULL
atlas.key$gordon$X <- NULL

