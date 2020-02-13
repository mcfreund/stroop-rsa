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

atlas.key$mmp$X <- NULL
atlas.key$gordon$X <- NULL

atlas.key$mmp[atlas.key$mmp$hemi == "L", "num.roi"] <- atlas.key$mmp$num.roi[atlas.key$mmp$hemi == "L"] + 180
atlas.key$mmp[atlas.key$mmp$hemi == "R", "num.roi"] <- atlas.key$mmp$num.roi[atlas.key$mmp$hemi == "R"] - 180


write.masks <- function(l, prefix, atlas.name) {
  
  .atlas <- atlas.key[[atlas.name]]
  .atlas.key <- atlas.key[[atlas.name]]
  
  for (ii in seq_along(l)) {
    num.rois <- .atlas$num.roi[.atlas.key$roi %in% l[[ii]]]
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
  lppc  = c("IP0", "IPS1", "IP1", "MIP", "7PL", "7Pm", "7Am", "7PL", "7AL", "7PC", "VIP", "LIPv", "LIPd", "AIP", "IP2"),
  mfc   = c("SCEF", "8BM",  "p32pr", "a32pr"),  ## medial fc
  dlpfc = c("p9-46v", "i6-8", "8Av", "8C", "46"),
  dpm   = c("6a", "FEF", "6ma"), ## dorsal premotor
  vpm   = c("6v", "6v", "PEF"),  ## ventral premotor
  ifj   = c("IFJp", "IFJa"),
  vvis  = c("FFC", "PH", "VVC", "PIT", "V8", "LO2", "TE2p", "VMV3", "VMV2", "VMV1", "PHA3", "PHA2", "PHA1")
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
  core       = c("p9-46v", "a9-46v", "i6-8", "AVi", "8C", "IFJp", "IP2", "IP1", "PFm", "8BM"),
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
stats.subjs.tdic <- fread(
  here("out", "rsa", "stats", paste0("subjs_pro_bias_acc-only_mmp_pearson_residual_glm-tdic.csv"))
)
stats.subjs.tdic <- stats.subjs.tdic[is.analysis.group == TRUE & y == "rank", ]  ## EXCLUDE HELD OUT SUBJECTS!
stats.subjs.tdic <- stats.subjs.tdic[, "coef" := NULL]

stats.subjs.tdic %<>% full_join(atlas.key$mmp, by = "roi")

stats.group.tdic <- stats.subjs.tdic %>%
  group_by(num.roi, model, param) %>%
  summarize(
    v    = wilcox.test(beta, alternative = "greater")$statistic,
    p    = wilcox.test(beta, alternative = "greater")$p.value,
    beta = tanh(mean(atanh(beta))),  ## must be last in this summarize()
  ) %>%
  mutate(p.fdr = p.adjust(p, method = "fdr"), p.holm = p.adjust(p, method = "holm")) %>%
  ungroup

stats.group.tdic %<>% full_join(atlas.key$mmp, by = "num.roi") %>% as.data.table

rois <- unique(stats.group.tdic[p.fdr < 0.05]$roi)


## get intersection with anatomical / apriori ROIs

anatfunc <- lapply(anatomical, intersect, rois)
write.masks(anatfunc, "anatfunc_", "mmp")


## orig ----

orig <- list(
  evis.L     = paste0(c("V1", "V2", "V3"), "_L"),
  evis.R     = paste0(c("V1", "V2", "V3"), "_R"),
  vvis.L     = paste0(c("FFC", "VVC", "V8", "VMV3"), "_L"),
  vvis.R     = paste0(c("FFC", "VVC", "V8", "VMV3", "VMV2"), "_R"),
  ips.L      = paste0(c("IP0", "IP1", "LIPd", "VIP", "LIPv", "AIP", "7PC"), "_L"),
  ips.R      = paste0(c("IP0", "IP1", "IP2", "IPS1", "MIP", "LIPd", "LIPv", "AIP", "7PC"), "_R"),
  premotor.L = paste0(c("FEF", "55b", "PEF", "6v", "6r", "43", "6d"), "_L"),
  premotor.R = paste0(c("FEF", "55b", "PEF", "6v", "6r", "43", "6a"), "_R"),
  dlpfc.L    = paste0(c("p9-46v", "8C", "8Av", "i6-8"), "_L"),
  dlpfc.R    = paste0(c("p9-46v", "8C", "8Av", "i6-8"), "_R"),
  mpfc.L     = paste0(c("a32pr", "p32pr", "8BM", "SCEF"), "Ll"),
  mpfc.R     = paste0(c("a32pr", "p32pr", "8BM", "SCEF"), "_R"),
  ins.L      = paste0(c("FOP4", "FOP5", "FOP3"), "_L"),
  ins.R      = paste0(c("FOP4", "FOP5", "FOP3"), "_R"),
  ifc.L      = paste0(c("44", "45", "IFSa", "IFSp", "p47r", "p47l"), "_L"),
  ifc.R      = paste0(c("44", "45", "IFSa", "IFSp", "p47r", "p47l"), "_R"),
  lang.L     = paste0(c("PSL", "55b", "44", "45", "SFL"), "_L"),
  lang.R     = paste0(c("PSL", "55b", "44", "45", "SFL"), "_R")
)

write.masks(orig, "orig_", "mmp")


## func ----




