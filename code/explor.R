
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
    here("out", "rsa", "stats", paste0("subjs_pro_bias_acc-only_mmp_pearson_residual_glm-tdic.csv"))
  )
stats.subjs.mmp <- stats.subjs.mmp[is.analysis.group == TRUE & y == "rank", ]  ## EXCLUDE HELD OUT SUBJECTS!
stats.subjs.mmp <- stats.subjs.mmp[, "coef" := NULL]
stats.subjs.mmp %<>% full_join(atlas.key$mmp, by = "roi")


## group-level data:
stats.group.mmp <- fread(here("out", "rsa", "stats", "group_pro_bias_acc-only_mmp_pearson_residual_glm-tdic.csv"))
stats.group.mmp <- stats.group.mmp[measure == "beta" & y == "rank" & param %in% params.interest, ]
stats.group.mmp <- stats.group.mmp[, c("measure", "y") := NULL]


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








