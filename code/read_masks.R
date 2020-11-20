## about ----
## 
## reads .nii files that are masks for ROIs.
## 
## mike freund, 2019-01-10

fnames <- list.files(here::here("out", "masks"), pattern = "\\.nii")
fnames <- fnames[-grep("pscalar", fnames)]

masks <- lapply(
  fnames, 
  function(x)
  oro.nifti::readNIfTI(here::here("out", "masks", x), reorient = FALSE)
)

names(masks) <- gsub("\\.nii.*$", "", fnames)

# sum(masks$vwfa@.Data)  ## nvox in VWFA mask
