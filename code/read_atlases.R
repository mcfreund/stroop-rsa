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

atlas$gordon <- oro.nifti::readNIfTI(here::here("out", "atlases", "gordon.nii.gz"), reorient = FALSE)


## keys

atlas.key <- setNames(vector("list", 2),  c("mmp", "gordon"))

atlas.key$mmp <- read.csv(here::here("out", "atlases", "mmp.csv"), stringsAsFactors = FALSE)
atlas.key$gordon <- read.csv(here::here("out", "atlases", "gordon.csv"), stringsAsFactors = FALSE)

atlas.key$mmp$X <- NULL
atlas.key$gordon$X <- NULL


## rearrange for workbench underlay

inds.left <- atlas.key$mmp$roi %>% grep("_L$", .)
inds.right <- atlas.key$mmp$roi %>% grep("_R$", .)
atlas.key$mmp$roi <- atlas.key$mmp$roi[c(inds.right, inds.left)]  ## correctly order mmp atlas

## for workbench underlay (get filename, write if non-existent)

fname.pscalar.mmp <- here(
  "out", "wb", 
  "Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.pscalar.nii"
)

if (!file.exists(fname.pscalar.mmp)) cifti.parcellate(name.atlas = "glasser", dir.to.write = here("out", "wb"))

