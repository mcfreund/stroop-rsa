## about ----
## 

library(here)
library(mikeutils)
library(magrittr)
library(dplyr)
library(data.table)
library(oro.nifti)
source(here("code", "strings.R"))
source(here("code", "read_atlases.R"))
source(here("code", "funs.R"))

write.masks <- function(l, atlas.name) {
  
  .atlas <- atlas[[atlas.name]]
  .atlas.key <- atlas.key[[atlas.name]]
  
  for (ii in seq_along(l)) {
    num.rois <- .atlas.key$num.roi[.atlas.key$roi %in% l[[ii]]]
    mask.i <- array(as.numeric(.atlas@.Data %in% num.rois), dim = c(75, 90, 75))
    
    writeNIfTI(as.nifti(mask.i), here::here("out", "masks", names(l)[ii]))
    
  }
  
}


## define superparcels ----

## from MMP:

superparcels <- list(
  ## lateral posterior parietal cortex (intraparietal sulcus):
  ## parcels tiling intraparietal sulcus (from IP0 to S1 [area 2]) and IPL (medial aspect) / SPL (lateral aspect)
  lppc     = c("IP0", "IPS1", "IP1", "MIP", "7PL", "7AL", "7PC", "VIP", "LIPv", "LIPd", "AIP", "IP2"),
  ## dorsomedial fc:
  ##  SMA/pre-sma: SCEF, 8BM
  ##  dACC: a32pr, p32pr, d32 (d32 is default mode in Cole-Anticevic)
  dmfc      = c("SCEF", "8BM",  "p32pr", "a32pr"),
  dmfc.alt  = c("SCEF", "8BM",  "p32pr", "a32pr", "d32"),
  ## dlpfc:
  ## parcels tiling middle frontal gyrus
  ## 46 is default mode in cole-anticevic
  ## a9-46v, 9-46d are rostral, can be considered fronto-polar cortex
  dlpfc     = c("p9-46v", "i6-8", "8Av", "8C"), ## along middle frontal gyrus
  dlpfc.alt = c("p9-46v", "i6-8", "8Av", "8C", "46", "a9-46v", "9-46d"),
  dpm       = c("6a", "FEF", "6ma"), ## dorsal premotor
  vpm       = c("6v", "6v", "PEF", "6r", "55b", "IFJp"),  ## ventral premotor
  vvis      = c("FFC", "VVC", "V8", "VMV3"),
  ins       = c("FOP4", "FOP5", "FOP3", "AVI", "MI"),
  ifc       = c("IFJa", "44", "45", "IFSa", "p47r", "IFSp", "a47r", "47l"),
  ofc       = c("47s", "47m", "11l", "13l", "OFC", "10pp", "10v", "10r", "s32", "p32", "a24", "25", "pOFC"),
  fpc       = c("46", "9-46d", "a9-46v", "a10p", "p10p", "a47r", "9a", "9p", "10d")
)

superparcels.r <- lapply(superparcels, function(x) paste0(x, "_R"))
superparcels.l <- lapply(superparcels, function(x) paste0(x, "_L"))
superparcels <- c(
  setNames(superparcels.l, paste0(names(superparcels.r), "_L")),
  setNames(superparcels.r, paste0(names(superparcels.r), "_R")),
  list(
    aud   = combo_paste(c("LBelt", "A1", "MBelt", "PBelt"), c("L", "R")),
    evis  = paste0("V1", c("_L", "_R"))
  )
)

write.masks(superparcels, "mmp")


## from gordon:

num.rois <- atlas.key$gordon$num.roi[atlas.key$gordon$roi %in% c("smmouth_L", "smmouth_R")]
smmouth <- array(as.numeric(atlas$gordon@.Data %in% num.rois), dim = c(75, 90, 75))
writeNIfTI(as.nifti(smmouth), here::here("out", "masks", "smmouth"))


## write parcel IDs ----

saveRDS(superparcels, here::here("out", "masks", "ids_superparcels.RDS"))


## create workbench files ----

## all superparcels

# superparcels <- atlas.key$mmp
# superparcels$id <- "none"
# for (superparcel.i in seq_along(superparcels)) {
#   superparcels$id[superparcels$roi %in% superparcels[[superparcel.i]]] <- names(superparcels)[superparcel.i]
# }
# superparcels$id <- relevel(as.factor(superparcels$id), "none")
# overlay <- superparcels %>% select(roi, id)
# overlay$id <- as.numeric(overlay$id) - 1
# 
# inds.left <- atlas.key$mmp$roi %>% grep("_L$", .)
# inds.right <- atlas.key$mmp$roi %>% grep("_R$", .)
# atlas.key$mmp$roi <- atlas.key$mmp$roi[c(inds.right, inds.left)]
# 
# overlay$roi.num <- match(overlay$roi, atlas.key$mmp$roi)
# overlay %<>% arrange(roi.num)
# 
# cifti.convert(
#   fname.overlay = "superparcels_superparcels",
#   values.overlay = overlay$id,
#   dir.template = here("out", "figs"),
#   name.atlas = "glasser",
#   dir.to.write = here("out", "figs", "ms_v1_2020-03", "indif_explor")
# )


## only lateral and medial pfc, ips

# superparcels$id <- "none"
# anatfunc.dissoc <- anatfunc[grep("lppc|mfc|dlpfc", names(anatfunc))]
# for (superparcel.i in seq_along(anatfunc.dissoc)) {
#   superparcels$id[superparcels$roi %in% anatfunc.dissoc[[superparcel.i]]] <- names(anatfunc.dissoc)[superparcel.i]
# }
# superparcels$id <- relevel(as.factor(superparcels$id), "none")
# overlay <- superparcels %>% select(roi, id)
# overlay$id <- as.numeric(overlay$id) - 1
# 
# inds.left <- atlas.key$mmp$roi %>% grep("_L$", .)
# inds.right <- atlas.key$mmp$roi %>% grep("_R$", .)
# atlas.key$mmp$roi <- atlas.key$mmp$roi[c(inds.right, inds.left)]
# 
# overlay$roi.num <- match(overlay$roi, atlas.key$mmp$roi)
# overlay %<>% arrange(roi.num)
# 
# cifti.convert(
#   fname.overlay = "superparcels_anatfunc_dissoc",
#   values.overlay = overlay$id,
#   dir.template = here("out", "figs"),
#   name.atlas = "glasser",
#   dir.to.write = here("out", "figs", "ms_v1_2020-03", "indif_dissoc")
# )
