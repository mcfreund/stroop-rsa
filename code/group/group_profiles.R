library(here)
library(knitr)
library(magrittr)
library(dplyr)
library(data.table)
library(ggplot2)
library(mikeutils)
library(doParallel)
library(foreach)
library(gifti)
library(fANCOVA)
library(viridis)
library(colorspace)
library(boot)
library(vegan)
library(grid)
library(gridExtra)
library(lme4)
source(here("code", "strings.R"))
source(here("code", "funs.R"))
source(here("code", "read_atlases.R"))

## underlay images

hcp <- list(
  L  = readGIfTI(
    file.path(dir.atlas, "surf", "HCP_S1200_GroupAvg_v1", "S1200.L.very_inflated_MSMAll.32k_fs_LR.surf.gii")
  ),
  R = readGIfTI(
    file.path(dir.atlas, "surf", "HCP_S1200_GroupAvg_v1", "S1200.R.very_inflated_MSMAll.32k_fs_LR.surf.gii")
  )
)

## MMP template for defining masks

mmp <- list(
  L = read_gifti2matrix(file.path(dir.atlas, "MMP_surface", "mmpL.func.gii")) %>% c,
  R = read_gifti2matrix(file.path(dir.atlas, "MMP_surface", "mmpR.func.gii")) %>% c
)

mmp$L[mmp$L > 0] <- mmp$L[mmp$L > 0] - 180
mmp$R[mmp$R > 0] <- mmp$R[mmp$R > 0] + 180


## data

stats.subjs.tdic <- fread(
  here("out", "rsa", "stats", paste0("subjs_pro_bias_acc-only_mmp_pearson_residual_glm-tdic.csv"))
)
stats.subjs.tdic <- stats.subjs.tdic[is.analysis.group == TRUE & y == "rank", ]  ## EXCLUDE HELD OUT SUBJECTS!
stats.subjs.tdic <- stats.subjs.tdic[, "coef" := NULL]
stats.subjs.tdic %<>% full_join(atlas.key$mmp, by = "roi")

stats.group.tdic <- fread(here("out", "rsa", "stats", "group_pro_bias_acc-only_mmp_pearson_residual_glm-tdic.csv"))
stats.group.tdic <- stats.group.tdic[measure == "beta" & y == "rank" & param %in% params.interest, ]
stats.group.tdic <- stats.group.tdic[, c("measure", "y") := NULL]

stats.pairs.tdic <- fread(
  here("out", "group", paste0("group_pro_bias_acc-only_mmp_pearson_residual_glm_tdic_pairwise.csv"))
)
stats.tosts.tdic <- fread(
  here("out", "group", paste0("group_pro_bias_acc-only_mmp_pearson_residual_glm_tdic_tost.csv"))
)


## 'coding' results ----

rois.targt <- stats.group.tdic[param == "target" & p.fdr < 0.05, roi]
rois.distr <- stats.group.tdic[param == "distractor" & p.fdr < 0.05, roi]
rois.incon <- stats.group.tdic[param == "incongruency" & p.fdr < 0.05, roi]
rois.congr <- stats.group.tdic[param == "congruency" & p.fdr < 0.05, roi]

keep.these.bois <- c("m", "v", "p.fdr", "num.roi", "roi", "hemi")
stats.coding.targt <- stats.group.tdic[param == "target", ..keep.these.bois]
stats.coding.distr <- stats.group.tdic[param == "distractor", ..keep.these.bois]
stats.coding.incon <- stats.group.tdic[param == "incongruency", ..keep.these.bois]
stats.coding.congr <- stats.group.tdic[param == "congruency", ..keep.these.bois]

## intersections

rois.targt.incon <- intersect(rois.targt, rois.incon)
rois.targt.distr <- intersect(rois.targt, rois.distr)
rois.distr.incon <- intersect(rois.distr, rois.incon)  ## none
rois.targt.distr.incon <- Reduce(intersect, list(rois.targt, rois.incon, rois.distr))  ## none


## pairwise comparisons ----

rois.pref.targt.vs.distr <- stats.pairs.tdic[roi %in% rois.targt & p.targt.distr < 0.05 & b.targt.distr > 0, roi]
rois.pref.targt.vs.incon <- stats.pairs.tdic[roi %in% rois.targt & p.targt.incon < 0.05 & b.targt.incon > 0, roi]
rois.pref.distr.vs.targt <- stats.pairs.tdic[roi %in% rois.distr & p.targt.distr < 0.05 & b.targt.distr < 0, roi]
rois.pref.distr.vs.incon <- stats.pairs.tdic[roi %in% rois.distr & b.incon.distr < 0.05 & b.incon.distr < 0, roi]
rois.pref.incon.vs.distr <- stats.pairs.tdic[roi %in% rois.incon & b.incon.distr < 0.05 & b.incon.distr > 0, roi]
rois.pref.incon.vs.targt <- stats.pairs.tdic[roi %in% rois.incon & p.targt.incon < 0.05 & b.targt.incon < 0, roi]

rois.pref.targt <- intersect(rois.pref.targt.vs.incon, rois.pref.targt.vs.distr)
rois.pref.distr <- intersect(rois.pref.distr.vs.targt, rois.pref.distr.vs.incon)
rois.pref.incon <- intersect(rois.pref.incon.vs.targt, rois.pref.incon.vs.distr)

## equivalency tests ----

rois.tost.targt.vs.distr <- stats.tosts.tdic[roi %in% rois.pref.targt.vs.distr & param == "distractor", ] %>%
  ungroup %>%
  mutate(p.fdr = p.adjust(p, method = "fdr"), p.holm = p.adjust(p, method = "holm")) %>%
  filter(p.fdr < 0.05) %>%
  .$roi

rois.tost.targt.vs.incon <- stats.tosts.tdic[roi %in% rois.pref.targt.vs.incon & param == "incongruency", ] %>%
  ungroup %>%
  mutate(p.fdr = p.adjust(p, method = "fdr"), p.holm = p.adjust(p, method = "holm")) %>%
  filter(p.fdr < 0.05) %>%
  .$roi

rois.target.selective <- intersect(rois.tost.targt.vs.incon, rois.tost.targt.vs.distr)

rois.tost.distr.vs.incon <- stats.tosts.tdic[roi %in% rois.pref.distr.vs.incon & param == "incongruency", ] %>%
  ungroup %>%
  mutate(p.fdr = p.adjust(p, method = "fdr"), p.holm = p.adjust(p, method = "holm")) %>%
  filter(p.fdr < 0.05) %>%
  .$roi

rois.tost.incon.vs.distr <- stats.tosts.tdic[roi %in% rois.pref.incon.vs.distr & param == "distractor", ] %>%
  ungroup %>%
  mutate(p.fdr = p.adjust(p, method = "fdr"), p.holm = p.adjust(p, method = "holm")) %>%
  filter(p.fdr < 0.05) %>%
  .$roi

rois.tost.incon.vs.targt <- stats.tosts.tdic[roi %in% rois.pref.incon.vs.targt & param == "target", ] %>%
  ungroup %>%
  mutate(p.fdr = p.adjust(p, method = "fdr"), p.holm = p.adjust(p, method = "holm")) %>%
  filter(p.fdr < 0.05) %>%
  .$roi

stats.group.tdic.coding <- stats.group.tdic[roi %in% unique(c(rois.distr, rois.incon, rois.pref.targt.vs.distr))]


## figure, panel a: brains -----

## counts:

length(rois.distr)
length(rois.targt)
length(rois.incon)
length(rois.pref.targt.vs.distr)
length(rois.pref.targt.vs.incon)

rois.profiles <- list(
  distr       = rois.distr,
  targt0distr = setdiff(rois.tost.targt.vs.distr, rois.tost.targt.vs.incon),
  targt0incon = setdiff(rois.tost.targt.vs.incon, rois.tost.targt.vs.distr),
  targt0distrincon = intersect(rois.tost.targt.vs.incon, rois.tost.targt.vs.distr)
)

## reshape to wide, add profiles

stats.group.tdic.profiles.wide <- stats.group.tdic %>%
  select(param, m, roi, community.short) %>%
  filter(param %in% params.interest) %>%
  tidyr::spread(param, m) %>%
  filter(roi %in% unlist(rois.profiles)) %>%
  mutate(
    prof = ifelse(
      roi %in% rois.profiles$distr, "distr",
      ifelse(
       roi %in% rois.profiles$targt0distr, "targt0distr",
       ifelse(
        roi %in% rois.profiles$targt0incon, "targt0incon",
        ifelse(
          roi %in% rois.profiles$targt0distrincon, "targt0distrincon",
          NA
        )
      )
    )
  )
)

overlay <- rbind(
  stats.group.tdic.profiles.wide %>% select(roi, prof),
  data.frame(
    roi = setdiff(atlas.key$mmp$roi, unlist(rois.profiles)), 
    prof = "none", stringsAsFactors = FALSE
  )
)
overlay$profnum <- as.numeric(factor(overlay$prof, levels = c("none", unique(stats.group.tdic.coding.wide$prof)))) - 1

inds.left <- atlas.key$mmp$roi %>% grep("_L$", .)
inds.right <- atlas.key$mmp$roi %>% grep("_R$", .)
atlas.key$mmp$roi <- atlas.key$mmp$roi[c(inds.right, inds.left)]

overlay$roi.num <- match(overlay$roi, atlas.key$mmp$roi)
overlay %<>% arrange(roi.num)

cifti.parcellate(name.atlas = "glasser", dir.to.write = here("out", "figs"))
cifti.convert(
  fname.overlay = "profiles",
  values.overlay = overlay$profnum,
  dir.template = here("out", "figs"),
  name.atlas = "glasser",
  dir.to.write = here("out", "figs", "ms_v1_2020-03", "group_profiles")
)


## figure, panel b: mds ----


rsarray <- readRDS(
  here("out", "rsa", "obsv", "rsarray_pro_bias_acc-only_mmp_pearson_residual-rank.rds")
)[, , unique(stats.subjs.tdic$subj), ]
# dimnames(rsarray)

## define examples

eg <- list(
  
  distr = stats.group.tdic.coding %>%
    filter(roi %in% rois.coding$distr, param == "distractor") %>%
    filter(v == max(v)) %>%
    pull("roi"),
  
  targt0distr = stats.group.tdic.coding %>%
    filter(roi %in% rois.profiles$targt0distr) %>%
    select(param, roi, v) %>%
    tidyr::spread(param, v) %>%
    mutate(targt.incon = incongruency + target) %>%
    filter(targt.incon == max(targt.incon)) %>%
    pull("roi"),
  
  targt0incon = stats.group.tdic.coding %>%
    filter(roi %in% rois.coding$targt0incon, param == "target") %>%
    filter(v == max(v)) %>%
    pull("roi")
)


## mds

# eg.mds <- lapply(eg, nmds, arr = rsarray)

eg.mds <- lapply(eg, nmds.rank, arr = rsarray)

## plots

p.mds <- lapply(eg.mds, plot.mds)
# p.mds

## add lines

mds.lines <- list(
  distr = list(
    ## distractor coding
    ## BLUE
    mds.line("blueBLUE", "purpleBLUE"),
    mds.line("purpleBLUE", "redBLUE"),
    mds.line("redBLUE", "whiteBLUE"),
    mds.line("whiteBLUE", "blueBLUE"),
    ## WHITE
    mds.line("blueWHITE", "purpleWHITE"),
    mds.line("purpleWHITE", "whiteWHITE"),
    mds.line("whiteWHITE", "redWHITE"),
    mds.line("redWHITE", "blueWHITE"), 
    ## RED
    mds.line("blueRED", "whiteRED"),
    mds.line("whiteRED", "redRED"),
    mds.line("redRED", "purpleRED"),
    mds.line("purpleRED", "blueRED"),
    ## PURPLE
    mds.line("redPURPLE", "bluePURPLE"),
    mds.line("bluePURPLE", "purplePURPLE"), 
    mds.line("purplePURPLE", "whitePURPLE"),
    mds.line("whitePURPLE", "redPURPLE")
  ),
  targt0distr = list(
    ## congruency & target
    ## blue
    mds.line("blueBLUE", "blueWHITE"),
    mds.line("blueBLUE", "blueRED"),
    mds.line("blueBLUE", "bluePURPLE"),
    ## white
    mds.line("whiteWHITE", "whiteBLUE"),
    mds.line("whiteWHITE", "whiteRED"),
    mds.line("whiteWHITE", "whitePURPLE"),
    ## red
    mds.line("redRED", "redBLUE"),
    mds.line("redRED", "redWHITE"),
    mds.line("redRED", "redPURPLE"),
    ## purple
    mds.line("purplePURPLE", "purpleBLUE"),
    mds.line("purplePURPLE", "purpleWHITE"),
    mds.line("purplePURPLE", "purpleRED")
  ),
  targt0incon = list(
    ## target
    ## blue
    mds.line("blueRED", "blueWHITE"),
    mds.line("blueWHITE", "blueBLUE"),
    mds.line("blueBLUE", "bluePURPLE"),
    mds.line("bluePURPLE", "blueRED"),
    ## white
    mds.line("whiteWHITE", "whitePURPLE"),
    mds.line("whitePURPLE", "whiteBLUE"),
    mds.line("whiteBLUE", "whiteRED"),
    mds.line("whiteRED", "whiteWHITE"),
    ## red
    mds.line("redBLUE", "redRED"),
    mds.line("redRED", "redWHITE"),
    mds.line("redWHITE", "redPURPLE"),
    mds.line("redPURPLE", "redBLUE"),
    ## purple
    mds.line("purplePURPLE", "purpleRED"),
    mds.line("purpleRED", "purpleBLUE"),
    mds.line("purpleBLUE", "purpleWHITE"),
    mds.line("purpleWHITE", "purplePURPLE")
  )
)

## get grobs

grobs <- list(
  distr       = plot.mds(eg.mds$distr, .add.lines = mds.lines$distr, .size = rel(7)),
  targt0incon = plot.mds(eg.mds$targt0incon, .add.lines = mds.lines$targt0incon, .size = rel(7)),
  targt0distr = plot.mds(eg.mds$targt0distr, .add.lines = mds.lines$targt0distr, .size = rel(7))
) %>%
  lapply(ggplotGrob)

labs.roi <- unlist(eg)

## draw

for (ii in seq_along(grobs)) {
  
  name.ii <- names(grobs)[ii]
  
  cairo_pdf(
    here("out", "figs", "ms_v1_2020-03", "group_profiles", paste0("mds_distr_", eg[[name.ii]], "_.pdf")), 
    height = 5, width = 5
  )
  
  grid.draw(grobs[[name.ii]])
  
  dev.off()
  
}

