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

## variables

params.interest <- c("target", "distractor", "incongruency")
colors.profs <- c(
  incon = "#d95f02ff", targt = "#1b9e77ff", distr = "#7570b3ff", targt.incon = "#26190dff", targt.distr = "#4682b4ff"
)


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

## get p values

stats.pairs.tdic <- stats.subjs.tdic %>%
  filter(param %in% params.interest) %>%
  group_by(roi) %>%
  summarize(
    out = list(
      pairwise.wilcox.test(beta, param, paired = TRUE, alternative = "two.sided", p.adjust.method = "fdr")
    )
  )
pvals <- vapply(stats.pairs.tdic$out, function(.) .$p.value[lower.tri(diag(2), diag = TRUE)], numeric(3))
pvals <- t(pvals)
colnames(pvals) <- c("p.incon.distr", "p.targt.distr", "p.targt.incon")
pvals <- apply(pvals, 2, p.adjust, method = "fdr")
stats.pairs.tdic <- data.frame(roi = stats.pairs.tdic$roi, pvals, stringsAsFactors = FALSE)

## get mean diffs

stats.pairs.tdic %<>% full_join(
  stats.subjs.tdic %>%
    filter(param %in% params.interest) %>%
    group_by(roi) %>%
    summarize(
      b.incon.distr = tanh(mean(beta[param == "incongruency"] - beta[param == "distractor"])),
      b.targt.distr = tanh(mean(beta[param == "target"] - beta[param == "distractor"])),
      b.targt.incon = tanh(mean(beta[param == "target"] - beta[param == "incongruency"])),
    ),
  by = "roi"
)
stats.pairs.tdic %<>% full_join(atlas.key$mmp, by = "roi") %>% as.data.table


## results

rois.pref.targt.vs.distr <- stats.pairs.tdic[roi %in% rois.targt & p.targt.distr < 0.05 & b.targt.distr > 0, roi]
rois.pref.targt.vs.incon <- stats.pairs.tdic[roi %in% rois.targt & p.targt.incon < 0.05 & b.targt.incon > 0, roi]
rois.pref.distr.vs.targt <- stats.pairs.tdic[roi %in% rois.distr & p.targt.distr < 0.05 & b.targt.distr < 0, roi]
rois.pref.distr.vs.incon <- stats.pairs.tdic[roi %in% rois.distr & b.incon.distr < 0.05 & b.incon.distr < 0, roi]
rois.pref.incon.vs.distr <- stats.pairs.tdic[roi %in% rois.incon & b.incon.distr < 0.05 & b.incon.distr > 0, roi]
rois.pref.incon.vs.targt <- stats.pairs.tdic[roi %in% rois.incon & p.targt.incon < 0.05 & b.targt.incon < 0, roi]

rois.pref.targt <- intersect(rois.pref.targt.vs.incon, rois.pref.targt.vs.distr)
rois.pref.distr <- intersect(rois.pref.distr.vs.targt, rois.pref.distr.vs.incon)
rois.pref.incon <- intersect(rois.pref.incon.vs.targt, rois.pref.incon.vs.distr)




## figure ----


## (a) venn diagram ----

length(rois.distr)
length(rois.targt)
length(rois.incon)
length(intersect(rois.incon, rois.targt))
length(intersect(rois.distr, rois.targt))

# rois.coding <- list(
#   distr = setdiff(rois.distr, c(rois.targt, rois.incon)),
#   targt = setdiff(rois.targt, c(rois.distr, rois.incon)),
#   incon = setdiff(rois.incon, c(rois.targt, rois.distr)),
#   targt.distr = intersect(rois.targt, rois.distr),
#   targt.incon = intersect(rois.targt, rois.incon)
# )


rois.coding <- list(
  distr = setdiff(rois.distr, c(rois.pref.targt.vs.distr, rois.incon)),
  targt = setdiff(rois.pref.targt, c(rois.distr, rois.incon)),
  taskr = setdiff(rois.pref.targt.vs.distr, c(rois.distr, rois.incon, rois.pref.targt)),
  incon = setdiff(rois.incon, c(rois.pref.targt.vs.distr, rois.distr)),
  targt.distr = intersect(rois.targt, rois.distr),
  targt.incon = intersect(rois.pref.targt.vs.distr, rois.incon)
)

## (b) brains ----

## reshape to wide, add profiles

stats.group.tdic.coding.wide <- stats.group.tdic %>%
  select(param, m, roi, community.short) %>%
  filter(param %in% params.interest) %>%
  tidyr::spread(param, m) %>%
  filter(roi %in% unlist(rois.coding)) %>%
  mutate(
    prof = ifelse(
      roi %in% rois.coding$taskr, "taskr",
      # ifelse(
      #   roi %in% rois.coding$targt, "targt",
        ifelse(
          roi %in% rois.coding$targt.distr, "targt.distr",
          ifelse(
            roi %in% rois.coding$distr, "distr",
            ifelse(
              roi %in% rois.coding$incon, "incon",
              ifelse(
                roi %in% rois.coding$targt.incon, "targt.incon",
                NA
              )
            )
          )
        )
      ),
    targt = roi %in% rois.targt,
    distr = roi %in% rois.distr,
    incon = roi %in% rois.incon
    )
  # )

stats.subjs.tdic.coding <- stats.subjs.tdic %>%
  filter(roi %in% unlist(rois.coding), param != "congruency") %>%
  select(subj, param, beta, roi, community.short)

overlay <- rbind(
  stats.group.tdic.coding.wide %>% select(roi, prof),
  data.frame(
    roi = setdiff(atlas.key$mmp$roi, unlist(rois.coding)), 
    prof = "none", stringsAsFactors = FALSE
  )
)
overlay$profnum <- as.numeric(factor(overlay$prof, levels = c("none", unique(stats.group.tdic.coding.wide$prof)))) - 1
# overlay %<>% mutate(is.left = grepl("_L$", roi)) %>% arrange(is.left, roi.num)  ## weirdness with MMP atlas! flip hemis.

inds.left <- atlas.key$mmp$roi %>% grep("_L$", .)
inds.right <- atlas.key$mmp$roi %>% grep("_R$", .)
atlas.key$mmp$roi <- atlas.key$mmp$roi[c(inds.right, inds.left)]

overlay$roi.num <- match(overlay$roi, atlas.key$mmp$roi)
overlay %<>% arrange(roi.num)

cifti.parcellate(name.atlas = "glasser", dir.to.write = here("out", "figs"))
cifti.convert(
  fname.overlay = "group_coding",
  values.overlay = overlay$profnum,
  dir.template = here("out", "figs"),
  name.atlas = "glasser",
  dir.to.write = here("out", "figs", "ohbm")
)


## indiv models

# overlay.targt <- rbind(
#   data.frame(
#     roi = rois.targt, 
#     is.sig = TRUE
#   ),
#   data.frame(
#     roi = setdiff(atlas.key$mmp$roi, rois.targt), 
#     is.sig = FALSE
#   )
# )

overlay.targt <- stats.group.tdic %>% filter(param == "target") %>% select(roi, num.roi, hemi, v, p.fdr, m)
overlay.targt <- overlay.targt %>% mutate(val = ifelse(p.fdr < 0.05, v, 0))
overlay.targt %<>%  arrange(rev(hemi), num.roi)  ## weirdness with MMP atlas! flip hemis.
cifti.convert(
  fname.overlay = "group_target",
  values.overlay = overlay.targt$val,
  dir.template = here("out", "figs"),
  name.atlas = "glasser",
  dir.to.write = here("out", "figs", "ohbm")
)


overlay.distr <- stats.group.tdic %>% filter(param == "distractor") %>% select(roi, num.roi, hemi, v, p.fdr, m)
overlay.distr <- overlay.distr %>% mutate(val = ifelse(p.fdr < 0.05, v, 0))
overlay.distr %<>%  arrange(rev(hemi), num.roi)  ## weirdness with MMP atlas! flip hemis.
cifti.convert(
  fname.overlay = "group_distractor",
  values.overlay = overlay.distr$val,
  dir.template = here("out", "figs"),
  name.atlas = "glasser",
  dir.to.write = here("out", "figs", "ohbm")
)


overlay.incon <- stats.group.tdic %>% filter(param == "incongruency") %>% select(roi, num.roi, hemi, v, p.fdr, m)
overlay.incon <- overlay.incon %>% mutate(val = ifelse(p.fdr < 0.05, v, 0))
overlay.incon %<>%  arrange(rev(hemi), num.roi)  ## weirdness with MMP atlas! flip hemis.
cifti.convert(
  fname.overlay = "group_incongruency",
  values.overlay = overlay.incon$val,
  dir.template = here("out", "figs"),
  name.atlas = "glasser",
  dir.to.write = here("out", "figs", "ohbm")
)



## dissoc ----


unique(atlas.key$mmp$community.short)
atlas.key$mmp[atlas.key$mmp$community.short %in% c("dlPFC", "PM", "aCC and mPFC", "pCC") & atlas.key$mmp$roi %in% rois.coding$taskr, ]

# atlas.key$mmp[atlas.key$mmp$roi %in% rois.coding$taskr, ]
# atlas.key$mmp[atlas.key$mmp$roi %in% rois.coding$targt, ]

subjs <- unique(stats.subjs.tdic$subj)

# rois.lateral <- atlas.key$mmp[
#   atlas.key$mmp$community.short %in% c("dlPFC", "PM") & atlas.key$mmp$roi %in% c(rois.pref.targt.vs.distr, rois.distr, rois.incon), 
#   "roi"
#   ]
# rois.medial <- atlas.key$mmp[
#   atlas.key$mmp$community.short %in% c("aCC and mPFC") & atlas.key$mmp$roi %in% c(rois.pref.targt.vs.distr, rois.distr, rois.incon), 
#   "roi"
#   ]

# rois.lateral <- atlas.key$mmp[
#   atlas.key$mmp$community.short %in% c("dlPFC", "PM") & atlas.key$mmp$roi %in% rois.coding$taskr, 
#   "roi"
#   ]
# rois.medial <- atlas.key$mmp[
#   atlas.key$mmp$community.short %in% c("aCC and mPFC") & atlas.key$mmp$roi %in% c(rois.coding$incon, rois.coding$targt.incon), 
#   "roi"
#   ]
rois.coding$dist

rois.medial <- combo_paste(c("8BM", "SCEF", "p32pr", "a32pr"), c("L", "R"))
rois.lateral <- combo_paste(c("i6-8", "p9-46v", "8Av", "8C"), c("L", "R"))


# rois.medial <- c("8BM_R", "p32pr_R")
# rois.lateral <- c("p9-46v_R", "FEF_R")
rois.visual <- intersect(c("V1_L", "V2_L", "V3_L", "V1_R", "V2_R", "V3_R"), rois.distr)
rois.sommot <- intersect(c("3b_L", "4_L", "3b_R", "4_R"), rois.targt)
rois <- c(rois.medial, rois.lateral, rois.visual, rois.sommot)

bounds <- 0.96
crit <- qnorm(1 - (1 - bounds) / 2)
n.resamples <- 1E4
colors.model <- c(incongruency = "#d95f02", target = "#1b9e77", distractor = "#7570b3")

d <- stats.subjs.tdic %>% filter(roi %in% rois)
e <- d %>% filter(param != "congruency")


unique(e$roi)
unique(e$region)
e$region <- ""
e$region[e$roi %in% rois.sommot] <- "somatomotor"
e$region[e$roi %in% rois.visual] <- "visual"
e$region[e$roi %in% rois.medial] <- "mFC"
e$region[e$roi %in% rois.lateral] <- "lFC"


e.pfc.i <- e %>%
  filter(param == "incongruency", region %in% c("lFC", "mFC")) %>%
  # group_by(roi, subj) %>%
  group_by(subj, param, region) %>% summarize(beta = mean(beta)) %>%
  mutate(beta_jbar = mean(beta)) %>%
  ungroup %>%
  mutate(beta_bar = mean(beta)) %>%
  group_by(region, param) %>%
  summarize(
    beta = mean(beta),
    ci = sqrt((var(beta - beta_jbar + beta_bar) * 1/2) / n()) * crit
  )

e.pfc.t <- e %>%
  filter(param == "target", region %in% c("lFC", "mFC")) %>%
  # group_by(roi, subj) %>%
  group_by(subj, param, region) %>% summarize(beta = mean(beta)) %>%
  mutate(beta_jbar = mean(beta)) %>%
  ungroup %>%
  mutate(beta_bar = mean(beta)) %>%
  group_by(region, param) %>%
  summarize(
    beta = mean(beta),
    ci = sqrt((var(beta - beta_jbar + beta_bar) * 1/2) / n()) * crit
  )

e.pfc <- bind_rows(e.pfc.i, e.pfc.t)

e.sm.d <- e %>%
  filter(param == "distractor", region %in% c("somatomotor", "visual")) %>%
  # group_by(roi, subj) %>%
  group_by(subj, param, region) %>% summarize(beta = mean(beta)) %>%
  mutate(beta_jbar = mean(beta)) %>%
  ungroup %>%
  mutate(beta_bar = mean(beta)) %>%
  group_by(region, param) %>%
  summarize(
    beta = mean(beta),
    ci = sqrt((var(beta - beta_jbar + beta_bar) * 1/2) / n()) * crit
  )

e.sm.t <- e %>%
  filter(param == "target", region %in% c("somatomotor", "visual")) %>%
  # group_by(roi, subj) %>%
  group_by(subj, param, region) %>% summarize(beta = mean(beta)) %>%
  mutate(beta_jbar = mean(beta)) %>%
  ungroup %>%
  mutate(beta_bar = mean(beta)) %>%
  group_by(region, param) %>%
  summarize(
    beta = mean(beta),
    ci = sqrt((var(beta - beta_jbar + beta_bar) * 1/2) / n()) * crit
  )

e.sm <- bind_rows(e.sm.d, e.sm.t)


ylims <- e %>%
  group_by(param, roi) %>%
  summarize(b = mean(beta))

e %<>% mutate(region = factor(region, levels = c("visual", "somatomotor", "mFC", "lFC")))

# e %>% filter(roi %in% c(rois.lateral, rois.medial))


p.frontal <- e %>%  
  
  filter(region %in% c("lFC", "mFC"), param %in% c("incongruency", "target")) %>%
  
  ggplot(aes(param, beta, color = region)) +
  
  geom_point(
    data = . %>% group_by(roi, param, region) %>% summarize(beta = mean(beta))
    # position = position_dodge(width = 0.25)
  ) +
  
  geom_line(
    data = . %>% group_by(param, roi, region) %>% summarize(beta = mean(beta)),
    aes(group = roi),
    # position = position_dodge(width = 0.25),
    size = 0.5,
    linetype = "dashed"
  ) +

  geom_line(
    data = . %>% group_by(param, region) %>% summarize(beta = mean(beta)),
    aes(group = region),
    # position = position_dodge(width = 0.5),
    size = 2,
  ) +
    
  geom_point(
    data = . %>% group_by(param, region) %>% summarize(beta = mean(beta)),
    size = 10
  ) +
  
  geom_errorbar(
    data = e.pfc,
    aes(ymin = beta + ci, ymax = beta - ci),
    width = 0,
    size = 2
  ) +
  
  scale_y_continuous(
    breaks = c(0, ylims %>% pull(b) %>% max %>% round(2))
  ) +
  
  scale_color_manual(values = c(lFC = "grey30", mFC = "firebrick")) +
  
  labs(y = "RSA model fit (β)") +
  
  annotate(
    geom = "text", x = 1.5, y = 0.08, 
    label = "medial FC", 
    # label = paste0("medial FC (", paste0(rois.medial, collapse = ", "), ")"), 
    hjust = 0, fontface = "bold.italic", size = 7, color = "firebrick"
    ) +
  annotate(
    geom = "text", x = 1.5, y = 0.0725, 
    label = paste0("lateral FC"), 
    # label = paste0("lateral FC (all 'task-relevant preference' \nparcels in dlPFC and premotor)"), 
    hjust = 0, fontface = "bold.italic", size = 7, color = "grey30"
    ) +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    # axis.text.x = element_text(face = "bold", size = rel(4), angle = 15, hjust = 0.5, vjust = 0.5),
    axis.text.x = element_text(size = rel(2.5), color = colors.model[c("incongruency", "target")], face = "bold"),
    # axis.text.y = element_text(color = rev(labels$colors), face = "bold", size = rel(1)),
    axis.title = element_blank(),
    panel.border = element_blank(),
    # axis.line.x = element_line(size = 2),
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  )
  

p.vismot <- e %>%  
  
  filter(region %in% c("somatomotor", "visual"), param %in% c("distractor", "target")) %>%
  
  ggplot(aes(param, beta, color = region)) +
  
  geom_point(
    data = . %>% group_by(roi, param, region) %>% summarize(beta = mean(beta))
    # position = position_dodge(width = 0.25)
  ) +
  
  geom_line(
    data = . %>% group_by(param, roi, region) %>% summarize(beta = mean(beta)),
    aes(group = roi),
    # position = position_dodge(width = 0.25),
    size = 0.5,
    linetype = "dashed"
  ) +
  
  geom_line(
    data = . %>% group_by(param, region) %>% summarize(beta = mean(beta)),
    aes(group = region),
    # position = position_dodge(width = 0.5),
    size = 2,
  ) +
  
  geom_point(
    data = . %>% group_by(param, region) %>% summarize(beta = mean(beta)),
    size = 10
  ) +
  
  geom_errorbar(
    data = e.sm,
    aes(ymin = beta + ci, ymax = beta - ci),
    width = 0,
    size = 2
  ) +
  
  geom_segment(
    aes(
      x = -Inf, xend = -Inf, y = 0, 
      yend = ylims %>% pull(b) %>% max %>% round(2)
    ),
    size = 4,
    color = "black"
  ) +
  
  scale_y_continuous(
    breaks = c(0, ylims %>% pull(b) %>% max %>% round(2))
  ) +
  
  scale_color_manual(values = c(visual = "#1f78b4", somatomotor = "#a6cee3")) +
  
  labs(y = "RSA model fit (β)") +
  
  annotate(geom = "text", x = 1, y = 0.08, label = "somato-motor", hjust = 0, fontface = "bold.italic", size = 7, color = "#a6cee3") +
  annotate(geom = "text", x = 1, y = 0.075, label = "early visual", hjust = 0, fontface = "bold.italic", size = 7, color = "#1f78b4") +
  
  # annotate(geom = "text", x = 1, y = 0.08, label = "somatomotor (3b_L, 4_L)", hjust = 0, fontface = "bold.italic", size = 7, color = "#a6cee3") +
  # annotate(geom = "text", x = 1, y = 0.075, label = "visual (V1_L, V2_L)", hjust = 0, fontface = "bold.italic", size = 7, color = "#1f78b4") +
  
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "bold", size = rel(2)),
    # axis.text.x = element_text(face = "bold", size = rel(4), angle = 15, hjust = 0.5, vjust = 0.5),
    axis.text.x = element_text(size = rel(2.5), color = colors.model[c("distractor", "target")], face = "bold"),
    # axis.text.y = element_text(color = rev(labels$colors), face = "bold", size = rel(1)),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = rel(2)),
    panel.border = element_blank(),
    # axis.line.x = element_line(size = 2),
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  )
  

grid.arrange(p.vismot, p.frontal, nrow = 1)



ggsave(
  here("docs", "posters", "ohbm2020", "figs", "group_dissoc_lineplots.pdf"),
  plot = arrangeGrob(p.vismot, p.frontal, nrow = 1),
  units = "cm",
  device = "pdf",
  height = 15,
  width = 30
)



fit.lm.t <- lmer(beta ~ region + (1 | subj), e %>% filter(roi %in% c(rois.lateral, rois.medial) & param == "target"))
summary(fit.lm.t)

fit.lm.c <- lmer(beta ~ region + (1 | subj), e, subset = roi %in% c(rois.lateral, rois.medial) & param == "incongruency")
summary(fit.lm.c)

fit.lm.c <- lmer(
  beta ~ region + (1 | subj), 
  e %>% filter(roi %in% c(rois.lateral, rois.medial) & param == "incongruency") %>%
    group_by(region, subj) %>%
    summarize(beta = mean(beta))
  )
summary(fit.lm.c)

fit.lm.t <- lmer(
  beta ~ region + (1 | subj), 
  e %>% filter(roi %in% c(rois.lateral, rois.medial) & param == "target") %>%
    group_by(region, subj) %>%
    summarize(beta = mean(beta))
)
summary(fit.lm.t)


## within region test:
fit.mfc <- lmer(beta ~ param + (1 | subj), e %>% filter(roi %in% rois.medial, param %in% c("target", "incongruency")))
summary(fit.mfc)





fit.sm.d <- lmer(beta ~ region + (1 | subj), e, subset = roi %in% c(rois.sommot, rois.visual) & param == "distractor")
summary(fit.sm.d)

fit.sm.t <- lmer(beta ~ region + (1 | subj), e, subset = roi %in% c(rois.sommot, rois.visual) & param == "target")
summary(fit.sm.t)


## ----


library(colorspace)


plot.model <- function(df, .title, .axis.text.size = rel(2)) {
  
  axis.title.colors.x <- rep(ifelse(bias.colors == "white", "grey50", bias.colors), each = 4)
  axis.title.colors.y <- rev(axis.title.colors.x)
  
  df %<>% mutate(Var1 = factor(Var1, levels = rev(levels(Var1))))
  
  p <- df %>%
    ggplot(aes(x = Var1, y = Var2)) +
    geom_tile(aes(fill = value), color = "white") +
    scale_fill_gradientn(colors = rev(diverge_hcl(200)), limits = c(-1, 1)) +
    theme(
      panel.background = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none",
      axis.text   = element_text(face = "bold", size = .axis.text.size),
      axis.text.x = element_text(
        color = axis.title.colors.x,
        angle = 90, hjust = 1, vjust = 0.5, margin = margin(-2.5, 0, 0, 0)
      ),
      axis.text.y = element_text(color = axis.title.colors.y, margin = margin(0, -2.5, 0, 0))
    ) +
    scale_x_discrete(labels = rep(bias.words, 4)) +
    scale_y_discrete(labels = rev(rep(bias.words, 4)))
  
  if (.title != "") p <- p + labs(title = .title) + theme(plot.title = element_text(face = "bold.italic", size = rel(4)))

  p
  
}


model.target <- as.matrix(fread(here("out", "rsa", "mods", "rsm_bias_target.csv"), data.table = FALSE)[, -1])
rownames(model.target) <- colnames(model.target)

model.distractor <- as.matrix(fread(here("out", "rsa", "mods", "rsm_bias_distractor.csv"), data.table = FALSE)[, -1])
rownames(model.distractor) <- colnames(model.distractor)

model.incongruency <- as.matrix(fread(here("out", "rsa", "mods", "rsm_bias_incongruency.csv"), data.table = FALSE)[, -1])
rownames(model.incongruency) <- colnames(model.incongruency)


p.mod.target <- model.target[sort(bias.items), sort(bias.items)] %>% reshape2::melt() %>% plot.model(.title = "target")
p.mod.distractor <- model.distractor[sort(bias.items), sort(bias.items)] %>% reshape2::melt() %>% plot.model(.title = "distractor")
p.mod.incongruency <- model.incongruency[sort(bias.items), sort(bias.items)] %>% reshape2::melt() %>% plot.model(.title = "incongruency")


ggsave(
  here("docs", "posters", "ohbm2020", "figs", "model_target.pdf"),
  plot = p.mod.target,
  units = "cm",
  device = "pdf",
  height = 15,
  width = 14
)


ggsave(
  here("docs", "posters", "ohbm2020", "figs", "model_incongruency.pdf"),
  plot = p.mod.incongruency,
  units = "cm",
  device = "pdf",
  height = 15,
  width = 14
)


ggsave(
  here("docs", "posters", "ohbm2020", "figs", "model_distractor.pdf"),
  plot = p.mod.distractor,
  units = "cm",
  device = "pdf",
  height = 15,
  width = 14
)
