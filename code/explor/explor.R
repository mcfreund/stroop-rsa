
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
library(magick)
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
    here("out", "rsa", "stats", paste0("subjs_pro_bias_acc-only_mmp_pearson_residual_glm-tdic_tmp.csv"))
  )
stats.subjs.mmp <- stats.subjs.mmp[is.analysis.group == TRUE & y == "rank", ]  ## EXCLUDE HELD OUT SUBJECTS!
stats.subjs.mmp <- stats.subjs.mmp[, "coef" := NULL]
stats.subjs.mmp %<>% full_join(atlas.key$mmp, by = "roi")

stats.group.mmp <- stats.subjs.mmp %>%
  
  filter(is.analysis.group, param %in% c("target", "distractor", "incongruency")) %>%
  group_by(param, roi) %>%
  
  summarize(
    p = wilcox.test(beta, alternative = "greater")$p.value,
    b = mean(beta)
  ) %>%
  
  mutate(p.fdr = p.adjust(p, method = "fdr"))
  
  








# ## group-level data:
# stats.group.mmp <- fread(here("out", "rsa", "stats", "group_pro_bias_acc-only_mmp_pearson_residual_glm-tdic_tmp.csv"))
# stats.group.mmp <- stats.group.mmp[measure == "beta" & y == "rank" & param %in% params.interest, ]
# stats.group.mmp <- stats.group.mmp[, c("measure", "y") := NULL]


## build workbench plots ----


(rois.targt <- stats.group.mmp[param == "target" & p.fdr < 0.05, roi])
(rois.distr <- stats.group.mmp[param == "distractor" & p.fdr < 0.01, roi])
(rois.incon <- stats.group.mmp[param == "incongruency" & p.fdr < 0.05, roi])
(rois.congr <- stats.group.mmp[param == "congruency" & p.fdr < 0.05, roi])

# keep.these.bois <- c("m", "v", "p.fdr", "num.roi", "roi", "hemi")
# stats.coding.targt <- stats.group.mmp[param == "target", ..keep.these.bois]
# stats.coding.distr <- stats.group.mmp[param == "distractor", ..keep.these.bois]
# stats.coding.incon <- stats.group.mmp[param == "incongruency", ..keep.these.bois]
# stats.coding.congr <- stats.group.mmp[param == "congruency", ..keep.these.bois]

## intersections

(rois.targt.incon <- intersect(rois.targt, rois.incon))
(rois.targt.distr <- intersect(rois.targt, rois.distr))
(rois.distr.incon <- intersect(rois.distr, rois.incon))
(rois.targt.distr.incon <- Reduce(intersect, list(rois.targt, rois.incon, rois.distr)))

## read or write files for workbench ----

fname.overlay.coding <- c(
  here("out", "explor", "target_veryinflated.tiff"), 
  here("out", "explor", "distractor_veryinflated.tiff"),
  here("out", "explor", "incongruency_veryinflated.tiff")
)

if (all(file.exists(fname.overlay.coding))) {
  
  p.coding <- lapply(fname.overlay.coding, image_read)
  grobs <- lapply(p.coding, rasterGrob)
  
  grid.arrange(
    arrangeGrob(grobs[[1]], top = textGrob("target", gp = gpar(fontsize = 24))), 
    arrangeGrob(grobs[[2]], top = textGrob("distractor", gp = gpar(fontsize = 24))),
    arrangeGrob(grobs[[3]], top = textGrob("incongruency", gp = gpar(fontsize = 24))),
    ncol = 3
  )
  
} else {
  
  overlay.targt <- stats.group.mmp %>% 
    
    filter(param == "target") %>% select(roi, num.roi, hemi, v, p.fdr, m) %>% 
    mutate(val = ifelse(p.fdr < 0.05, v, 0)) %>% 
    arrange(rev(hemi), num.roi)  ## weirdness with MMP atlas! flip hemis.
  
  cifti.convert(
    fname.overlay  = "group_target",
    dir.to.write   = here("out", "explor"),
    values.overlay = overlay.targt$val,
    dir.template   = here("out", "wb"),
    name.atlas     = "glasser"
  )
  
  overlay.distr <- stats.group.mmp %>% 
    
    filter(param == "distractor") %>% select(roi, num.roi, hemi, v, p.fdr, m) %>% 
    mutate(val = ifelse(p.fdr < 0.05, v, 0)) %>%
    arrange(rev(hemi), num.roi)  ## weirdness with MMP atlas! flip hemis.
  
  cifti.convert(
    fname.overlay  = "group_distractor",
    dir.to.write   = here("out", "explor"),
    values.overlay = overlay.distr$val,
    dir.template   = here("out", "wb"),
    name.atlas     = "glasser"
  )
  
  overlay.incon <- stats.group.mmp %>% 
    
    filter(param == "incongruency") %>% select(roi, num.roi, hemi, v, p.fdr, m) %>%
    mutate(val = ifelse(p.fdr < 0.05, v, 0)) %>% 
    arrange(rev(hemi), num.roi)  ## weirdness with MMP atlas! flip hemis.
  
  cifti.convert(
    fname.overlay  = "group_incongruency",
    dir.to.write   = here("out", "explor"),
    values.overlay = overlay.incon$val,
    dir.template   = here("out", "wb"),
    name.atlas     = "glasser"
  )
  
  ## now build and save plots with workbench, and re-run this .rmd file
  
}

rois.top.targt <- stats.group.mmp %>% filter(param == "target", roi %in% rois.targt) %>% slice_max(v, prop = 0.1) %>% pull(roi)

means.top.target <- stats.subjs.mmp %>%
  
  filter(roi %in% rois.top.targt, param %in% "target") %>%
  group_by(roi, community) %>%
  summarize(res = list(boot_mean_ci(beta))) %>% 
  
  tidyr::unnest(cols = c(res))


p.top.target <- 
  
  means.top.target %>%  
  
  ungroup %>%
  mutate(roi = factor(roi, levels = roi[order(y)])) %>%
  
  ggplot(aes(roi, y, color = community)) +
  geom_point(size = geom.point.size) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, size = geom.line.size) +
  # coord_capped_cart(left = "both") +
  labs(y = bquote("Target model fit ("*bar(beta)*")"), x = "Parcel") +
  scale_color_brewer(type = "qual", palette = 3) +
  scale_y_continuous(limits = c(0, 0.12), breaks = c(0, 0.05, 0.1), position = "right") +
  coord_flip() +
  coord_capped_flip(top = "both", left = "both") +
  theme(
    # legend.position = "none", 
    panel.grid      = element_blank(), 
    panel.border    = element_blank(),
    axis.line     = element_line(size = axis.line.size),
    axis.text       = element_text(size = axis.text.size*3/4),
    # axis.text.x     = element_text(color = colors.model),
    axis.ticks.x    = element_line(size = axis.line.size),
    # axis.ticks.x    = element_blank(),
    axis.title      = element_text(size = axis.title.size),
    # axis.text.y = element_text(angle = 90, vjust = 0, hjust = 1)
  )

ggsave(here("out", "explor", "top_target.pdf"), p.top.target, device = "pdf", width = 5, height = 3)




## whole-brain HLM ----
stats.subjs.mmp %<>% filter(param %in% params.interest)

fit.mmp <- lmer(
  beta ~ 0 + interaction(roi, param) + (param | subj), 
  stats.subjs.mmp
  # control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1E9))
)
fit.mmp.coef <- as.data.frame(coef(summary(fit.mmp)))

fit.mmp.coef$term <- rownames(fit.mmp.coef)
fit.mmp.coef$model <- gsub(".*(target|distractor|incongruency)", "\\1", fit.mmp.coef$term)
fit.mmp.coef$roi <- gsub("interaction\\(roi, param\\)(.*)\\.(target|distractor|incongruency)", "\\1", fit.mmp.coef$term)
# fit.mmp.coef %<>% rename(Pr)
fit.mmp.coef <- fit.mmp.coef %>%
  group_by(model) %>%
  mutate(p.fdr = p.adjust(`Pr(>|t|)`, method = "fdr"))

fit.mmp.coef %>% filter(p.fdr < 0.05, model == "incongruency") %>% View
fit.mmp.coef %>% filter(p.fdr < 0.05, model == "target") %>% View
fit.mmp.coef %>% filter(p.fdr < 0.05, model == "distractor")








fit.mmp.het.param <- lme(
  beta ~ 0 + interaction(roi, param),
  random = ~ 1 | subj, 
  weights = varIdent(form = ~ 1 | param),
  stats.subjs.mmp,
  method = "ML"
)
fit.mmp.het.all <- lme(
  beta ~ 0 + interaction(roi, param),
  random = ~ 1 | subj, 
  weights = varIdent(form = ~ 1 | interaction(roi, param)),
  stats.subjs.mmp,
  method = "ML"
)









is.distractor <- grepl("distractor", names(fixef(fit.super)))
is.incongruency <- grepl("incongruency", names(fixef(fit.super)))
is.target <- grepl("target", names(fixef(fit.super)))

is.dlpfc <- grepl("dlpfc", names(fixef(fit.super)))
is.lppc <- grepl("lppc", names(fixef(fit.super)))
is.dmfc <- grepl("dmfc", names(fixef(fit.super)))

is.left <- grepl("_L", names(fixef(fit.super)))
is.right <- grepl("_R", names(fixef(fit.super)))


contrasts.super <- rbind(
  
  ## region*model means (over hemi)
  
  DLPFC.target = is.dlpfc * is.target / 2,
  LPPC.target = is.lppc * is.target / 2,
  DMFC.target = is.dmfc * is.target / 2,
  
  DLPFC.distractor = is.dlpfc * is.distractor / 2,
  LPPC.distractor = is.lppc * is.distractor / 2,
  DMFC.distractor = is.dmfc * is.distractor / 2,
  
  DLPFC.incongruency = is.dlpfc * is.incongruency / 2,
  LPPC.incongruency = is.lppc * is.incongruency / 2,
  DMFC.incongruency = is.dmfc * is.incongruency / 2,
  
  ## within-region contrasts (over hemi)
  
  "target-distractor|DLPFC" = is.dlpfc * is.target - is.dlpfc * is.distractor / 2,
  "target-distractor|LPPC"  = is.lppc * is.target - is.lppc * is.distractor / 2,
  "target-distractor|DMFC"  = is.dmfc * is.target - is.dmfc * is.distractor / 2,
  
  "target-incongruency|DLPFC" = is.dlpfc * is.target - is.dlpfc * is.incongruency / 2,
  "target-incongruency|LPPC"  = is.lppc * is.target - is.lppc * is.incongruency / 2,
  "target-incongruency|DMFC"  = is.dmfc * is.target - is.dmfc * is.incongruency / 2,
  
  ## within-region target-incongruency contrsts (by hemi)
  
  "target-incongruency|DLPFC_L" = (is.dlpfc * is.target - is.dlpfc * is.incongruency) * is.left,
  "target-incongruency|LPPC_L" = (is.lppc * is.target - is.lppc * is.incongruency) * is.left,
  "target-incongruency|DMFC_L" = (is.dmfc * is.target - is.dmfc * is.incongruency) * is.left,
  
  "target-incongruency|DLPFC_R" = (is.dlpfc * is.target - is.dlpfc * is.incongruency) * is.right,
  "target-incongruency|LPPC_R" = (is.lppc * is.target - is.lppc * is.incongruency) * is.right,
  "target-incongruency|DMFC_R" = (is.dmfc * is.target - is.dmfc * is.incongruency) * is.right
  
)

glht.super <- summary(glht(fit.super, contrasts.super), test = adjusted("none"))

p.fdr.tvi.dlpfc.super <- glht.super$test$pvalues[grep("dlpfc_", rownames(contrasts.super))] %>% p.adjust(method = "fdr")
p.fdr.tvi.lppc.super <- glht.super$test$pvalues[grep("lppc_", rownames(contrasts.super))] %>% p.adjust(method = "fdr")
p.fdr.tvi.dmfc.super <- glht.super$test$pvalues[grep("dmfc_", rownames(contrasts.super))] %>% p.adjust(method = "fdr")

contrasts.super.all <- rbind(
  contrasts.super,
  ## cross-region contrasts:
  "DLPFC-DMFC_L|incongruency" = (is.dlpfc/2 - (is.dmfc * is.left)) * is.incongruency,
  "DLPFC-LPPC_L|incongruency" = (is.dlpfc/2 - (is.lppc * is.left)) * is.incongruency,
  "DLPFC-DMFC_L|target" = (is.dlpfc/2 - (is.dmfc * is.left)) * is.target,
  "DLPFC-LPPC_L|target" = (is.dlpfc/2 - (is.lppc * is.left)) * is.target,
  
  ## interactions:
  "(DLPFC-DMFC_L)(target-incongruency)" = 
    (is.dlpfc/2 - (is.dmfc * is.left)) * is.target - (is.dlpfc/2 - (is.dmfc * is.left)) * is.incongruency,
  "(DLPFC-LPPC_L)(target-incongruency)" = 
    (is.dlpfc/2 - (is.lppc * is.left)) * is.target - (is.dlpfc/2 - (is.lppc * is.left)) * is.incongruency
  
)

(glht.super.all <- summary(glht(fit.super, contrasts.super.all), test = adjusted("none")))


## table ----

table.group <- data.frame(
  contrast = names(glht.super.all$test$coefficients),
  b = glht.super.all$test$coefficients,
  se = glht.super.all$test$sigma,
  t = glht.super.all$test$tstat,
  p = glht.super.all$test$pvalues
)
table.group$p.raw <- table.group$p
table.group[
  table.group$contrast %in% c(names(p.fdr.tvi.dlpfc.super), names(p.fdr.tvi.lppc.super), names(p.fdr.tvi.dmfc.super)),
  "p"
  ] <- c(p.fdr.tvi.dlpfc.super, p.fdr.tvi.lppc.super, p.fdr.tvi.dmfc.super)
table.group$contrast %<>%
  gsub("\\.", " ", .) %>%
  gsub("\\|", " | ", .) %>%
  gsub("-", " - ", .)
rownames(table.group) <- NULL

kable(table.group)

fwrite(table.group, here("out", "group", "superparcels.txt"))


## figure ----

set.seed(0)
means.super <- stats.subjs.super %>%
  
  mutate(roi.hemi = paste0(roi, "_", hemi)) %>%
  filter(roi.hemi %in% c("dlpfc_L", "dlpfc_R", "dmfc_L", "lppc_L")) %>%
  
  group_by(roi, param, subj) %>%
  summarize(beta = mean(beta)) %>% 
  
  group_by(roi, param) %>%
  summarize(res = list(boot_mean_ci(beta))) %>% 
  
  tidyr::unnest(cols = c(res))


p.means.super <- means.super %>%  
  
  mutate(param = factor(as.factor(param), levels = c("incongruency", "target", "distractor"))) %>%
  
  ggplot(aes(param, y, group = roi, color = roi)) +
  geom_point(size = geom.point.size, position = position_dodge(width = 1/2)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position = position_dodge(width = 1/2), width = 0, size = geom.line.size) +
  geom_line(size = geom.line.size, position = position_dodge(width = 1/2)) +
  
  scale_color_manual(values = colors.region %>% setNames(tolower(names(.)))) +
  lemon::coord_capped_cart(left = "both") +
  # scale_y_continuous(breaks = c(0, 0.18), limits = c(-0.0205, 0.18)) +
  scale_x_discrete(labels = c("incon.", "targ.", "distr.")) +
  
  annotate(
    geom = "text", x = Inf, y = 0.1, label = "DLPFC (bil.)", color = colors.region["DLPFC"],
    hjust = 1, vjust = 1, size = label.size, fontface = "bold"
  ) +
  annotate(
    geom = "text", x = Inf, y = 0.12, label = "LPPC (L)", color = colors.region["LPPC"],
    hjust = 1, vjust = 1, size = label.size, fontface = "bold"
  ) +
  annotate(
    geom = "text", x = Inf, y = 0.14, label = "DMFC (L)", color = colors.region["DMFC"],
    hjust = 1, vjust = 1, size = label.size, fontface = "bold"
  ) +
  
  labs(y = bquote("Model fit ("*bar(beta)*")"), x = "Model") +
  
  theme(
    legend.position = "none", 
    panel.grid      = element_blank(), 
    panel.border    = element_blank(),
    axis.line.y     = element_line(size = axis.line.size),
    axis.text       = element_text(size = axis.text.size),
    # axis.text.x     = element_text(color = colors.model),
    axis.ticks.y    = element_line(size = axis.line.size),
    axis.ticks.x    = element_blank(),
    axis.title    = element_text(size = axis.title.size)
  )

ggsave(here("out", "group", "crossplot_superparcels.pdf"), p.means.super, device = "pdf", width = 4, height = 3)








