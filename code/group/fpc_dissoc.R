#+ fpc-dissoc_setup, include = FALSE

if (interactive()) source(here::here("code", "group", "mds.R"))

#' ### plot means from all ROIs

#+ fpc-dissoc_plot-all

stats.subjs.super %>%
  
  filter(param %in% params.interest, roi %in% c("dlpfc", "dmfc", "lppc")) %>%
  
  ggplot(aes(roi.hemi, beta, color = param)) +
  stat_summary(fun.data = mean_cl_boot, size = 2, position = position_dodge(width = 0.5)) +
  
  scale_color_manual(values = colors.model) +
  annotate(
    geom = "text", x = -Inf, y = 0.15, label = "incongruency", color = colors.model["incongruency"],
    hjust = 0, vjust = 1, size = rel(6)
  ) +
  annotate(
    geom = "text", x = -Inf, y = 0.14, label = "target", color = colors.model["target"],
    hjust = 0, vjust = 1, size = rel(6)
  ) +
  annotate(
    geom = "text", x = -Inf, y = 0.13, label = "distractor", color = colors.model["distractor"],
    hjust = 0, vjust = 1, size = rel(6)
  ) +
  
  theme(legend.position = "none") +
  labs(title = "superparcels")



#' ### model
#+ fpc-dissoc_model_primary, fig.height = 7, fig.width = 9

fit.super <- lmer(
  beta ~ 0 + interaction(roi.hemi, param) + (roi * param | subj), 
  stats.subjs.super %>% filter(roi %in% c("dlpfc", "dmfc", "lppc"))
)
summary(fit.super)

## build contrasts

is.distractor <- grepl("distractor", names(fixef(fit.super)))
is.incongruency <- grepl("incongruency", names(fixef(fit.super)))
is.target <- grepl("target", names(fixef(fit.super)))
is.dlpfc <- grepl("dlpfc", names(fixef(fit.super)))
is.lppc <- grepl("lppc", names(fixef(fit.super)))
is.dmfc <- grepl("dmfc", names(fixef(fit.super)))
is.left <- grepl("_L", names(fixef(fit.super)))
is.right <- grepl("_R", names(fixef(fit.super)))

contrasts.super <- rbind(
  
  ## region*model*hemi means
  
  DLPFC_L.target = is.dlpfc * is.target * is.left,
  LPPC_L.target = is.lppc * is.target * is.left,
  DMFC_L.target = is.dmfc * is.target * is.left,
  DLPFC_R.target = is.dlpfc * is.target * is.right,
  LPPC_R.target = is.lppc * is.target * is.right,
  DMFC_R.target = is.dmfc * is.target * is.right,
  
  DLPFC_L.distractor = is.dlpfc * is.distractor * is.left,
  LPPC_L.distractor = is.lppc * is.distractor * is.left,
  DMFC_L.distractor = is.dmfc * is.distractor * is.left,
  DLPFC_R.distractor = is.dlpfc * is.distractor * is.right,
  LPPC_R.distractor = is.lppc * is.distractor * is.right,
  DMFC_R.distractor = is.dmfc * is.distractor * is.right,
  
  DLPFC_L.incongruency = is.dlpfc * is.incongruency * is.left,
  LPPC_L.incongruency = is.lppc * is.incongruency * is.left,
  DMFC_L.incongruency = is.dmfc * is.incongruency * is.left,
  DLPFC_R.incongruency = is.dlpfc * is.incongruency * is.right,
  LPPC_R.incongruency = is.lppc * is.incongruency * is.right,
  DMFC_R.incongruency = is.dmfc * is.incongruency * is.right,
  
  ## within-region contrasts (over/across hemi)
  
  "target-distractor|DLPFC" = is.dlpfc * is.target - is.dlpfc * is.distractor / 2,
  "target-distractor|LPPC"  = is.lppc * is.target - is.lppc * is.distractor / 2,
  "target-distractor|DMFC"  = is.dmfc * is.target - is.dmfc * is.distractor / 2,
  
  "target-incongruency|DLPFC" = is.dlpfc * is.target - is.dlpfc * is.incongruency / 2,
  "target-incongruency|LPPC"  = is.lppc * is.target - is.lppc * is.incongruency / 2,
  "target-incongruency|DMFC"  = is.dmfc * is.target - is.dmfc * is.incongruency / 2,
  
  ## within-region target-incongruency contrsts (within/by hemi); these need correcting
  
  "target-incongruency|DLPFC_L" = (is.dlpfc * is.target - is.dlpfc * is.incongruency) * is.left,
  "target-incongruency|LPPC_L" = (is.lppc * is.target - is.lppc * is.incongruency) * is.left,
  "target-incongruency|DMFC_L" = (is.dmfc * is.target - is.dmfc * is.incongruency) * is.left,
  
  "target-incongruency|DLPFC_R" = (is.dlpfc * is.target - is.dlpfc * is.incongruency) * is.right,
  "target-incongruency|LPPC_R" = (is.lppc * is.target - is.lppc * is.incongruency) * is.right,
  "target-incongruency|DMFC_R" = (is.dmfc * is.target - is.dmfc * is.incongruency) * is.right
  
)

glht.super <- summary(glht(fit.super, contrasts.super), test = adjusted("none"))
## correct p-values from within-region contrats:
p.fdr.tvi.dlpfc.super <- glht.super$test$pvalues[grep("\\|DLPFC_", rownames(contrasts.super))] %>% p.adjust(method = "fdr")
p.fdr.tvi.lppc.super <- glht.super$test$pvalues[grep("\\|LPPC_", rownames(contrasts.super))] %>% p.adjust(method = "fdr")
p.fdr.tvi.dmfc.super <- glht.super$test$pvalues[grep("\\|DMFC_", rownames(contrasts.super))] %>% p.adjust(method = "fdr")

glht.super

## now add btw region contrasts
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

(glht.super.all <- summary(glht(fit.super, contrasts.super.all), test = adjusted("none")))  ## p-values uncorrected


#' ### build table
#+ fpc-dissoc_table_primary

table.group.type <- c(
  "condition mean",
  "condition mean",
  "condition mean",
  "condition mean",
  "condition mean",
  "condition mean",
  "condition mean",
  "condition mean",
  "condition mean",
  "condition mean",
  "condition mean",
  "condition mean",
  "condition mean",
  "condition mean",
  "condition mean",
  "condition mean",
  "condition mean",
  "condition mean",
  "w/i-region contrast",
  "w/i-region contrast",
  "w/i-region contrast",
  "w/i-region contrast",
  "w/i-region contrast",
  "w/i-region contrast",
  "w/i-region contrast",
  "w/i-region contrast",
  "w/i-region contrast",
  "w/i-region contrast",
  "w/in-region contrast",
  "wtn-region contrast",
  "b/w-region contrast",
  "b/w-region contrast",
  "b/w-region contrast",
  "b/w-region contrast",
  "2-way interaction",
  "2-way interaction"
)

table.group.contrast <- c(
  "$\\text{target DLPFC (L)}$",
  "$\\text{target LPPC (L)}$",
  "$\\text{target DMFC (L)}$",
  "$\\text{target DLPFC (R)}$",
  "$\\text{target LPPC (R)}$",
  "$\\text{target DMFC (R)}$",
  "$\\text{distractor DLPFC (L)}$",
  "$\\text{distractor LPPC (L)}$",
  "$\\text{distractor DMFC (L)}$",
  "$\\text{distractor DLPFC (R)}$",
  "$\\text{distractor LPPC (R)}$",
  "$\\text{distractor DMFC (R)}$",
  "$\\text{incongruency DLPFC (L)}$",
  "$\\text{incongruency LPPC (L)}$",
  "$\\text{incongruency DMFC (L)}$",
  "$\\text{incongruency DLPFC (R)}$",
  "$\\text{incongruency LPPC (R)}$",
  "$\\text{incongruency DMFC (R)}$",
  "$\\text{DLPFC: } \\text{target}-\\text{distractor}$",
  "$\\text{LPPC: } \\text{target}-\\text{distractor}$",
  "$\\text{DMFC: } \\text{target}-\\text{distractor}$",
  "$\\text{DLPFC: } \\text{target}-\\text{incongruency}$",
  "$\\text{LPPC: } \\text{target}-\\text{incongruency}$",
  "$\\text{DMFC: } \\text{target}-\\text{incongruency}$",
  "$\\text{DLPFC (L): } \\text{target}-\\text{incongruency}$",
  "$\\text{LPPC (L): } \\text{target}-\\text{incongruency}$",
  "$\\text{DMFC (L): } \\text{target}-\\text{incongruency}$",
  "$\\text{DLPFC (R): } \\text{target}-\\text{incongruency}$",
  "$\\text{LPPC (R): } \\text{target}-\\text{incongruency}$",
  "$\\text{DMFC (R): } \\text{target}-\\text{incongruency}$",
  "$\\text{incongruency: } \\text{DLPFC}-\\text{DMFC (L)}$",
  "$\\text{incongruency: } \\text{DLPFC}-\\text{LPPC (L)}$",
  "$\\text{target: } \\text{DLPFC}-\\text{DMFC (L)}$",
  "$\\text{target: } \\text{DLPFC}-\\text{LPPC (L)}$",
  "$(\\text{DLPFC}-\\text{DMFC (L)})\\times(\\text{target}-\\text{incongruency})$",
  "$(\\text{DLPFC}-\\text{LPPC (L)})\\times(\\text{target}-\\text{incongruency})$"
)

## bind
table.group <- data.frame(
  contrast = table.group.contrast,
  type = table.group.type,
  b = glht.super.all$test$coefficients,
  se = glht.super.all$test$sigma,
  t = glht.super.all$test$tstat,
  p = glht.super.all$test$pvalues
)
## correct p-values again
table.group$p.raw <- table.group$p
table.group[
  table.group$contrast %in% c(names(p.fdr.tvi.dlpfc.super), names(p.fdr.tvi.lppc.super), names(p.fdr.tvi.dmfc.super)),
  "p"
  ] <- c(p.fdr.tvi.dlpfc.super, p.fdr.tvi.lppc.super, p.fdr.tvi.dmfc.super)
rownames(table.group) <- NULL

kable(table.group)

fwrite(table.group, here("out", "group", "superparcels.txt"))



#' ### plot
#+ fpc-dissoc_plot_primary, fig.height = 7, fig.width = 9

set.seed(0)
means.super <- stats.subjs.super %>%
  
  mutate(roi.hemi = paste0(roi, "_", hemi)) %>%
  filter(roi.hemi %in% c("dlpfc_L", "dlpfc_R", "dmfc_L", "lppc_L")) %>%
  
  group_by(roi, param, subj) %>%
  summarize(beta = mean(beta), .groups = "drop_last") %>% 
  
  summarize(res = list(boot_mean_ci(beta)), .groups = "drop") %>% 
  
  tidyr::unnest(cols = c(res))


p.means.super <- means.super %>%  
  
  mutate(param = factor(as.factor(param), levels = c("incongruency", "target", "distractor"))) %>%
  
  ggplot(aes(param, y, group = roi, color = roi)) +
  geom_point(size = geom.point.size, position = position_dodge(width = 1/2)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position = position_dodge(width = 1/2), width = 0, size = geom.line.size) +
  geom_line(size = geom.line.size, position = position_dodge(width = 1/2)) +
  
  scale_color_manual(values = colors.region %>% setNames(tolower(names(.)))) +
  lemon::coord_capped_cart(left = "both") +
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
    axis.ticks.y    = element_line(size = axis.line.size),
    axis.ticks.x    = element_blank(),
    axis.title    = element_text(size = axis.title.size)
  )

ggsave(here("out", "group", "crossplot_superparcels.pdf"), p.means.super, device = "pdf", width = 4, height = 3)



#' ### arrange figure
#+ fpc-dissoc_arrange-fig_primary

## assumes p.mds object is in env (need to source mds.R script prior, or run _group.rmd)
plot.group <- plot_grid(
  p.mds,
  p.means.super,
  ncol = 2,
  vjust = c(1.25, 1.25),
  labels = "AUTO",
  label_size = 12
)

plot.group

ggsave(
  here("out", "group", "fig_group.pdf"), 
  plot.group, 
  device = "pdf", height = 4.25, width = 8.5, unit = "cm"
)
