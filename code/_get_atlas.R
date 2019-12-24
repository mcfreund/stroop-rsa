# if (!interactive()) stop("session not interactive")  ## commented for knitrs 2019-02-26

## about ----


## dependencies prompts and input validation ----

library(dplyr)

## do.mb

if (!exists("do.mb")) {
  do.mb <- c(
    mb4 = readline(prompt = "do mb 4? [y, *any key*]  (*any key* for no)    \n") == "y",
    mb8 = readline(prompt = "do mb 8? [y, *any key*]  (*any key* for no)    \n") == "y"
  )
}
if (!is.logical(do.mb)) stop("funky values of do.mb var")
if (!length(do.mb) > 1) stop("do.mb too long")
if (do.mb["mb4"]) {
  voxel.label <- "2p4"   ## for getting the correct ROI image below
} else if (do.mb["mb8"]) { voxel.label <- "222" }

## do.read.atlas

if (!exists("dir.atlas")) {
  source(here::here("..", "gen/funs/get_dirs_remote.R"))
  # print("sourcing ./gen/funs/get_dirs_remote.R")  ## commented for knitrs 2019-02-26
}
if (!exists("do.read.atlas")) {
  do.read.atlas <- c(
    mmp    = readline(prompt = "get mmp? [y, *any key*]  (*any key* for no)    \n") == "y",
    gordon = readline(prompt = "get gordon? [y, *any key*]  (*any key* for no)    \n") == "y"
  )
}
if (!is.logical(do.read.atlas) & !length(do.read.atlas) > 2) stop("funky values of do.read.atlas var")
if (any(do.read.atlas)) {
  if (!any(do.mb)) stop("!any(do.mb) == TRUE")
  atlas <- vector("list", sum(do.read.atlas)) %>% setNames(names(do.read.atlas)[do.read.atlas])
  for (ii in seq_along(atlas)) atlas[[ii]] <- vector("list", 2) %>% setNames(c("l", "r"))
  rm(ii)
}

## read atlases ----

if (do.read.atlas["mmp"]) {
  ## get atlas
  atlas$mmp$l <-oro.nifti::readNIfTI(
    paste0(dir.atlas, "/HCP-MMP1_L_on_MNI152_ICBM2009a_nlin_", voxel.label, ".nii.gz"),
    reorient = FALSE
    )
  atlas$mmp$r <- oro.nifti::readNIfTI(
    paste0(dir.atlas, "/HCP-MMP1_R_on_MNI152_ICBM2009a_nlin_", voxel.label, ".nii.gz"),
    reorient = FALSE
    )
}

if (do.read.atlas["gordon"]) {
  ## all this extra wrangling to make gordon atlas indices in same format as MMP: i.e.,
  ## same numbers for R and L hemispheres (one number per bilateral ROI), but separate objects
  ## for each hemisphere atlas.
  atlas$gordon$l <- oro.nifti::readNIfTI(
    file.path(dir.atlas, paste0("gordon/communities_", voxel.label, "_resampled_LPI.nii.gz")),
    reorient = FALSE
    )
  atlas$gordon$r <- atlas$gordon$l
  ## right hemi == even numbers; left hemi == odd numbers
  atlas$gordon$l[atlas$gordon$l %% 2 == 0] <- 0  ## true if even (i.e., sets right hemi to 0)
  atlas$gordon$r[atlas$gordon$r %% 2 == 1] <- 0  ## true if odd (i.e., sets left hemi to 0)
  atlas$gordon$r <- atlas$gordon$r - 1  ## shift by one (so r and l have same indices)
  atlas$gordon$r[atlas$gordon$r == -1 & !is.na(atlas$gordon$r)] <- 0  ## reset -1 to 0
  ## for checking:
  # unique(c(image.atlas.l))
  # unique(c(image.atlas.r))
  # c(image.atlas.l)[!is.na(c(image.atlas.l))] %*% c(image.atlas.r)[!is.na(c(image.atlas.r))]
}

## read keys ----

atlas.key <- vector("list", length(do.read.atlas)) %>% setNames(names(do.read.atlas))
atlas.key$mmp <- read.csv(file.path(dir.atlas, "mmp.csv"), stringsAsFactors = FALSE)
atlas.key$gordon <- read.table(
  file.path(dir.atlas, "communityKey_wsubcort.txt"),
  stringsAsFactors = FALSE
  ) %>%
  rename(num.roi = comm.num, roi = comm.lbl) %>%
  mutate(
    roi = gsub("(.*).$", "\\1", roi)
  ) %>%
  filter(num.roi %% 2 == 1)
