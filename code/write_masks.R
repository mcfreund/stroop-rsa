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

write.masks <- function(l, prefix, atlas.name) {
  
  .atlas <- atlas[[atlas.name]]
  .atlas.key <- atlas.key[[atlas.name]]
  
  for (ii in seq_along(l)) {
    num.rois <- .atlas.key$num.roi[.atlas.key$roi %in% l[[ii]]]
    mask.i <- array(as.numeric(.atlas@.Data %in% num.rois), dim = c(75, 90, 75))
    
    writeNIfTI(as.nifti(mask.i), here::here("out", "masks", paste0(prefix, names(l)[ii])))
    
  }
  
}


## anatomical ----

## from MMP
## auditory (bil)
## visual (bil)
## dmpfc (l, r)
## dlpfc (l, r)
## ips (l, r)

anatomical <- list(
  # lppc   = c("IP0", "IPS1", "IP1", "MIP", "7PL", "7Pm", "7Am", "7PL", "7AL", "7PC", "VIP", "LIPv", "LIPd", "AIP", "IP2"),
  lppc   = c("IP0", "IPS1", "IP1", "MIP", "7PL", "7AL", "7PC", "VIP", "LIPv", "LIPd", "AIP", "IP2"),
  ## lateral posterior parietal cortex (intraparietal sulcus):
  ## parcels tiling intraparietal sulcus (from IP0 to S1 [area 2]) and IPL (medial aspect) / SPL (lateral aspect)
  mfc    = c("SCEF", "8BM",  "p32pr", "a32pr", "d32"),
  ## dorsomedial fc:
  ##  SCEF, 8BM: SMA/pre-sma (Assem 2019)
  ##  p32pr: incongruency coding (right)
  ##  a32pr, d32: dACC (Assem 2019)
  dlpfc  = c("p9-46v", "i6-8", "8Av", "8C", "46"),
  ## along middle frontal sulcus
  dpm    = c("6a", "FEF", "6ma"), ## dorsal premotor
  vpm    = c("6v", "6v", "PEF", "6r", "55b"),  ## ventral premotor
  ifj    = c("IFJp", "IFJa"),
  vvis   = c("FFC", "VVC", "V8", "VMV3"),
  ins    = c("FOP4", "FOP5", "FOP3", "AVI", "MI"),
  ifc    = c("44", "45", "IFSa", "p47r", "IFSp", "a47r", "47l"),
  ofc    = c("47s", "47m", "11l", "13l", "OFC", "10pp", "10v", "10r", "s32", "p32", "a24", "25", "pOFC"),
  fpc    = c("9-46d", "a9-46v", "a10p", "p10p", "a47r", "9a", "9p", "10d")
)

anatomical.r <- lapply(anatomical, function(x) paste0(x, "_R"))
anatomical.l <- lapply(anatomical, function(x) paste0(x, "_L"))
anatomical <- c(
  setNames(anatomical.l, paste0(names(anatomical.r), "_L")),
  setNames(anatomical.r, paste0(names(anatomical.r), "_R")),
  list(
    aud   = combo_paste(c("LBelt", "A1", "MBelt", "PBelt"), c("L", "R")),
    evis  = combo_paste(c("V1", "V2", "V3"), c("L", "R"))
  )
)

write.masks(anatomical, "anat_", "mmp")


## from gordon
## smmouth (bil)

num.rois <- atlas.key$gordon$num.roi[atlas.key$gordon$roi %in% c("smmouth_L", "smmouth_R")]
smmouth <- array(as.numeric(atlas$gordon@.Data %in% num.rois), dim = c(75, 90, 75))
writeNIfTI(as.nifti(smmouth), here::here("out", "masks", "anat_smmouth"))


## assem (mmp) ----

md <- list(
  core       = c("p9-46v", "a9-46v", "i6-8", "AVi", "8C", "IFJp", "IP2", "IP1", "PFm", "8BM", "SCEF"),
  extended   = c(
    "p9-46v", "a9-46v", "i6-8", "AVi", "8C", "IFJp", "IP2", "IP1", "PFm", "8BM",
    "TE1m", "TE1p", "PGs", "PFm", "AIP", "MIP", "LIPd", "IP1", "IP2", "s6-8", 
    "i6-8", "a9-46v", "FOP5", "AVI", "11l", "a10p", "p10p", "a47r", "p47r"
  )
)

md.r <- lapply(md, function(x) paste0(x, "_R"))
md.l <- lapply(md, function(x) paste0(x, "_L"))
md <- c(
  setNames(md.l, paste0(names(md.l), "_L")),
  setNames(md.r, paste0(names(md.r), "_R"))
)

write.masks(md, "md_", "mmp")


## anat * func ----

## get names of task-responsive ROIs

stats.subjs <- fread(
  here("out", "rsa", "stats", paste0("group_pro_bias_acc-only_mmp_pearson_residual_glm-tdic.csv"))
)

stats.group.sig <- stats.subjs %>%
  filter(y == "rank", param %in% c("target", "incongruency", "distractor"), measure == "beta", p.fdr < 0.05)


rois <- unique(stats.group.sig$roi)


## get intersection with anatomical / apriori ROIs

anatfunc <- lapply(anatomical, intersect, rois)
write.masks(anatfunc, "anatfunc_", "mmp")



## write parcel IDs ----

saveRDS(anatfunc, here::here("out", "masks", "ids_anatfunc.RDS"))
saveRDS(anatomical, here::here("out", "masks", "ids_anat.RDS"))
saveRDS(md, here::here("out", "masks", "ids_md.RDS"))


## create workbench files ----

## all superparcels

superparcels <- atlas.key$mmp
superparcels$id <- "none"
for (superparcel.i in seq_along(anatomical)) {
  superparcels$id[superparcels$roi %in% anatomical[[superparcel.i]]] <- names(anatomical)[superparcel.i]
}
superparcels$id <- relevel(as.factor(superparcels$id), "none")
overlay <- superparcels %>% select(roi, id)
overlay$id <- as.numeric(overlay$id) - 1

inds.left <- atlas.key$mmp$roi %>% grep("_L$", .)
inds.right <- atlas.key$mmp$roi %>% grep("_R$", .)
atlas.key$mmp$roi <- atlas.key$mmp$roi[c(inds.right, inds.left)]

overlay$roi.num <- match(overlay$roi, atlas.key$mmp$roi)
overlay %<>% arrange(roi.num)

cifti.convert(
  fname.overlay = "superparcels_anatomical",
  values.overlay = overlay$id,
  dir.template = here("out", "figs"),
  name.atlas = "glasser",
  dir.to.write = here("out", "figs", "ms_v1_2020-03", "indif_explor")
)


## only lateral and medial pfc, ips

superparcels$id <- "none"
anatfunc.dissoc <- anatfunc[grep("lppc|mfc|dlpfc", names(anatfunc))]
for (superparcel.i in seq_along(anatfunc.dissoc)) {
  superparcels$id[superparcels$roi %in% anatfunc.dissoc[[superparcel.i]]] <- names(anatfunc.dissoc)[superparcel.i]
}
superparcels$id <- relevel(as.factor(superparcels$id), "none")
overlay <- superparcels %>% select(roi, id)
overlay$id <- as.numeric(overlay$id) - 1

inds.left <- atlas.key$mmp$roi %>% grep("_L$", .)
inds.right <- atlas.key$mmp$roi %>% grep("_R$", .)
atlas.key$mmp$roi <- atlas.key$mmp$roi[c(inds.right, inds.left)]

overlay$roi.num <- match(overlay$roi, atlas.key$mmp$roi)
overlay %<>% arrange(roi.num)

cifti.convert(
  fname.overlay = "superparcels_anatfunc_dissoc",
  values.overlay = overlay$id,
  dir.template = here("out", "figs"),
  name.atlas = "glasser",
  dir.to.write = here("out", "figs", "ms_v1_2020-03", "indif_dissoc")
)
