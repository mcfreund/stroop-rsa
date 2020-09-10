

library(here)
library(knitr)
library(magrittr)
library(data.table)
library(ggplot2)
library(mikeutils)
library(doParallel)
library(foreach)
library(gifti)
library(viridis)
library(colorspace)
library(boot)
library(vegan)
library(grid)
library(gridExtra)
library(lme4)
library(lmerTest)
library(multcomp)
library(dplyr)
# library(magick)
library(viridis)
library(cowplot)
library(lemon)
source(here("code", "strings.R"))
source(here("code", "funs.R"))
source(here("code", "read_atlases.R"))


## global settings


theme_set(theme_bw(base_size = 8))
n_cores <- detectCores()
nresamps <- 1E4
figwidth <- 4.2  ## cm

axis.text.size <- rel(1)
axis.title.size <- rel(1)
axis.line.size <- rel(1)
label.size <- rel(3)
p.value.size <- rel(2)
p.line.size <- rel(0.5)
geom.line.size <- rel(1)
geom.point.size <- rel(2)


## functions

boot_mean_ci <- function(x, R = 1E4, type = "bca", ...) {
  
  out <- boot::boot(x, statistic = function(x, ii) mean(x[ii]), R = R)
  ci <- boot::boot.ci(out, type = type)[[type]][4:5]
  
  data.frame(y = out$t0, ymin = ci[1], ymax = ci[2])
  
}


## strings


colors.model <- c(incongruency = "#d95f02", target = "#1b9e77", distractor = "#7570b3")
colors.region <- setNames(viridis(3), c("DLPFC", "LPPC", "DMFC"))

params.interest <- names(colors.model)

md <- list(
  core       = c("p9-46v", "a9-46v", "i6-8", "AVI", "8C", "IFJp", "IP2", "IP1", "PFm", "8BM", "SCEF"),
  extended   = c(
    "p9-46v", "a9-46v", "i6-8", "AVi", "8C", "IFJp", "IP2", "IP1", "PFm", "8BM",
    "TE1m", "TE1p", "PGs", "PFm", "AIP", "MIP", "LIPd", "IP1", "IP2", "s6-8", 
    "i6-8", "a9-46v", "FOP5", "AVI", "11l", "a10p", "p10p", "a47r", "p47r"
  )
)


## data

## subject-level data:
stats.subjs.mmp <- 
  fread(
    here("out", "rsa", "stats", paste0("subjs_pro_bias_acc-only_mmp_residual.csv"))
  )
stats.subjs.mmp <- stats.subjs.mmp[is.analysis.group == TRUE, ]  ## EXCLUDE HELD OUT SUBJECTS!
stats.subjs.mmp %<>% full_join(atlas.key$mmp, by = "roi")

stats.subjs.super <- 
  fread(
    here("out", "rsa", "stats", paste0("subjs_pro_bias_acc-only_masks_residual.csv"))
  ) %>%
  filter(is.analysis.group, param %in% c("target", "distractor", "incongruency")) %>%
  mutate(roi.hemi = roi, roi = gsub("_L|_R", "", roi), hemi = ifelse(grepl("_L", roi.hemi), "L", "R"))

## group-level data:
stats.group.mmp <- fread(here("out", "rsa", "stats", "group_pro_bias_acc-only_mmp_residual.csv"))
stats.group.mmp <- stats.group.mmp[measure == "beta" & y == "rank" & param %in% params.interest, ]
stats.group.mmp <- stats.group.mmp[, c("measure", "y") := NULL]

## similarity matrices:
rsarray.mmp <- readRDS(
  here("out", "rsa", "obsv", "rsarray_pro_bias_acc-only_mmp_residual-rank.rds")
)[, , unique(stats.subjs.mmp$subj), ]
rsarray.super <- readRDS(
  here("out", "rsa", "obsv", "rsarray_pro_bias_acc-only_masks_residual-rank.rds")
)[, , unique(stats.subjs.mmp$subj), ]


## underlay for workbench

inds.left <- atlas.key$mmp$roi %>% grep("_L$", .)
inds.right <- atlas.key$mmp$roi %>% grep("_R$", .)
atlas.key$mmp$roi <- atlas.key$mmp$roi[c(inds.right, inds.left)]  ## correctly order mmp atlas

fname.pscalar.mmp <- here(
  "out", "wb", 
  "Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.pscalar.nii"
)

if (!file.exists(fname.pscalar.mmp)) cifti.parcellate(name.atlas = "glasser", dir.to.write = here("out", "wb"))
