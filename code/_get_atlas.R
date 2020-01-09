## about ----


## read atlases ----

atlases.to.read <- c("mmp", "gordon")
atlas <- setNames(vector("list", length(atlases.to.read)), atlases.to.read)
for (ii in seq_along(atlas)) atlas[[ii]] <- setNames(vector("list", 2), c("l", "r"))
rm(ii)


## mmp --

atlas$mmp$l <- oro.nifti::readNIfTI(
  paste0(dir.atlas, "/HCP-MMP1_L_on_MNI152_ICBM2009a_nlin_2p4.nii.gz"),
  reorient = FALSE
  )
atlas$mmp$r <- oro.nifti::readNIfTI(
  paste0(dir.atlas, "/HCP-MMP1_R_on_MNI152_ICBM2009a_nlin_2p4.nii.gz"),
  reorient = FALSE
  )

## gordon --

atlas$gordon$l <- oro.nifti::readNIfTI(
  file.path(dir.atlas, paste0("gordon/communities_2p4_resampled_LPI.nii.gz")),
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

## vwfa --



## read keys ----

atlas.key <- setNames(vector("list", length(atlases.to.read)),  atlases.to.read)
atlas.key$mmp <- read.csv(file.path(dir.atlas, "mmp.csv"), stringsAsFactors = FALSE)
atlas.key$gordon <- read.table(
  file.path(dir.atlas, "communityKey_wsubcort.txt"),
  stringsAsFactors = FALSE
  ) %>%
  dplyr::rename(num.roi = comm.num, roi = comm.lbl) %>%
  dplyr::mutate(
    roi = gsub("(.*).$", "\\1", roi)
  ) %>%
  dplyr::filter(num.roi %% 2 == 1)
