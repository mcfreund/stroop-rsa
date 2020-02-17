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
  lppc   = c("IP0", "IPS1", "IP1", "MIP", "7PL", "7Pm", "7Am", "7PL", "7AL", "7PC", "VIP", "LIPv", "LIPd", "AIP", "IP2"),
  mfc    = c("SCEF", "8BM",  "p32pr", "a32pr"),  ## medial fc
  dlpfc  = c("p9-46v", "i6-8", "8Av", "8C", "46"),
  dpm    = c("6a", "FEF", "6ma"), ## dorsal premotor
  vpm    = c("6v", "6v", "PEF", "6r"),  ## ventral premotor
  ifj    = c("IFJp", "IFJa"),
  vvis   = c("FFC", "VVC", "V8", "VMV3"),
  ins    = c("FOP4", "FOP5", "FOP3", "AVI", "MI"),
  ifc    = c("44", "45", "IFSa", "p47r"),
  ofc    = c("47s", "47m", "11l", "13l", "OFC", "10pp"),
  fpc    = c("9-46d", "a9-46v", "a10p", "p10p", "a47r")
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


## func ----


tolower(c("IP0_r", "IP1_r", "IP2_r", "IPS1_r", "MIP_r", "LIPd_r", "LIPv_r", "AIP_r", "7PC_r", "7PL_r")) %>%
  setdiff(., tolower(anatomical$lppc_R))