## about
## 
## this script standardizes atlas formats.
## atlases used in this project are the MMP (Glasser 2016), and the Gordon (2014)
## 
## the standard atlas format consists of:
##  - a single nifti for each atlas
##  - each parcel indicated by a unique integer
##  - left hemisphere parcels assigned integers 1, 2, ... N/2
##  - right hemisphere parcels assigned integers N/2+1, N/2+2, ... N
##  - voxels not assigned to any roi are given a value of 0
##  - atlas keys have matching structures
## 
## this file writes these standardized atlases and atlas keys to ./out/atlases/
## this script should be run as a preliminary script, as all other scripts read atlases from ./out
## 

## atlas path dependent upon location

library(here)
library(magrittr)

source(here("code", "strings.R"))

## mmp ----

## NB. mmp is arranged in two separate hemispheres (b/c was warped from surface files).
## these hemispheres contain rois with the same valued integers (same range: 1 thru 180).
## thus, the right hemisphere values need to be increased by 180, and to be combined with the left.

## atlas

mmp.l <- oro.nifti::readNIfTI(
  paste0(dir.atlas, "/HCP-MMP1_L_on_MNI152_ICBM2009a_nlin_2p4.nii.gz"),
  reorient = FALSE
)
mmp.r <- oro.nifti::readNIfTI(
  paste0(dir.atlas, "/HCP-MMP1_R_on_MNI152_ICBM2009a_nlin_2p4.nii.gz"),
  reorient = FALSE
)

# sort(unique(c(mmp.l)))
# sum(is.na(c(mmp.l)))


mmp.r[mmp.r > 0] <- mmp.r[mmp.r > 0] + max(c(mmp.r))  ## add 180 to all values > 0
stopifnot(all.equal(sort(unique(c(mmp.r))), c(0, 181:360)))
stopifnot(all.equal(sort(unique(c(mmp.l))), c(0:180)))
stopifnot(c(mmp.r) %*% c(mmp.c) == 0)  ## overlap between hemis?

mmp <- mmp.r + mmp.l  ## combine hemispheres

## key

mmp.key.l <- read.csv(file.path(dir.atlas, "mmp.csv"), stringsAsFactors = FALSE)
mmp.key.r <- mmp.key.l

mmp.key.l$hemi <- "L"
mmp.key.r$hemi <- "R"

mmp.key.r$num.roi <- mmp.key.r$num.roi + max(mmp.key.r$num.roi)

mmp.key <- rbind(mmp.key.l, mmp.key.r)
mmp.key$roi <- paste0(mmp.key$roi, "_", mmp.key$hemi)
mmp.key <- mmp.key[c("roi", "num.roi", "hemi", "community", "community.short")]
# head(mmp.key)


## gordon ----

## NB. Gordon is arranged in one file (bilateral).
## however, right hemisphere parcels are assingned even integers, and left, odd.
## additionally, voxels not assigned to any ROI are NA, not zero (unlike MMP).
## these features are changed.

## just create two temporary files:

gordon.l <- oro.nifti::readNIfTI(
  file.path(dir.atlas, paste0("gordon/communities_2p4_resampled_LPI.nii.gz")),
  reorient = FALSE
)
gordon.r <- gordon.l

# sort(unique(c(gordon.l)))
# sum(is.na(c(gordon.l)))

gordon.l[is.na(gordon.l)] <- 0
gordon.r[is.na(gordon.r)] <- 0

gordon.l[gordon.l %% 2 == 0] <- 0  ## true if even (i.e., sets right hemi to 0)
gordon.r[gordon.r %% 2 == 1] <- 0  ## true if odd (i.e., sets left hemi to 0)

gordon.l[gordon.l > 0] <- (gordon.l[gordon.l > 0] + 1) / 2
gordon.r[gordon.r > 0] <- (gordon.r[gordon.r > 0]) / 2 + max(c(gordon.l))

stopifnot(all.equal(sort(unique(c(gordon.l))), 0:14))
stopifnot(all.equal(sort(unique(c(gordon.r))), c(0, 15:28)))
stopifnot(c(gordon.r) %*% c(gordon.l) == 0)  ## fail if overlap between hemis

gordon <- gordon.l + gordon.r

## key

gordon.key <- read.table(
  file.path(dir.atlas, "communityKey_wsubcort.txt"),
  stringsAsFactors = FALSE
  ) %>%
  dplyr::rename(num.roi = comm.num, roi = comm.lbl) %>%
  dplyr::mutate(
    hemi = gsub("(.*)(.$)", "\\2", roi),
    roi = paste0(tolower(gsub("(.*)(.$)", "\\1", roi)), "_", hemi),
  )

## left hemisphere:
gordon.key.l <- gordon.key[gordon.key$num.roi %% 2 == 1, ]
gordon.key.l$num.roi <- (gordon.key.l[gordon.key.l$num.roi %% 2 == 1, "num.roi"] + 1) / 2

## right:
gordon.key.r <- gordon.key[gordon.key$num.roi %% 2 == 0, ]
gordon.key.r$num.roi <- gordon.key.r[gordon.key.r$num.roi %% 2 == 0, "num.roi"] / 2 + max(c(gordon.l))

gordon.key <- rbind(gordon.key.l, gordon.key.r)


## write ----

oro.nifti::writeNIfTI(mmp, here("out", "atlases", "mmp"))
oro.nifti::writeNIfTI(gordon, here("out", "atlases", "gordon"))

write.csv(mmp.key, here("out", "atlases", "mmp.csv"))
write.csv(gordon.key, here("out", "atlases", "gordon.csv"))

