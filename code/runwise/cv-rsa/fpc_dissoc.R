#+ fpc-dissoc_setup, include = FALSE

#' ### plot means from all ROIs

#+ fpc-dissoc_plot-all

stats.subjs.super <- d.super[[measure]]

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

#+


#' ### model
#+ fpc-dissoc_model_primary, fig.height = 7, fig.width = 9

fit.super <- lmer(
  beta ~ 0 + interaction(roi.hemi, param) + (roi + param | subj), 
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
  
  "target-distractor|DLPFC" = (is.dlpfc * is.target - is.dlpfc * is.distractor) / 2,
  "target-distractor|LPPC"  = (is.lppc * is.target - is.lppc * is.distractor) / 2,
  "target-distractor|DMFC"  = (is.dmfc * is.target - is.dmfc * is.distractor) / 2,
  
  "target-incongruency|DLPFC" = (is.dlpfc * is.target - is.dlpfc * is.incongruency) / 2,
  "target-incongruency|LPPC"  = (is.lppc * is.target - is.lppc * is.incongruency) / 2,
  "target-incongruency|DMFC"  = (is.dmfc * is.target - is.dmfc * is.incongruency) / 2,
  
  "distractor-incongruency|DLPFC" = (is.dlpfc * is.distractor - is.dlpfc * is.incongruency) / 2,
  "distractor-incongruency|LPPC"  = (is.lppc * is.distractor - is.lppc * is.incongruency) / 2,
  "distractor-incongruency|DMFC"  = (is.dmfc * is.distractor - is.dmfc * is.incongruency) / 2,
  
  ## within-region target-incongruency contrasts (within/by hemi)
  
  "target-distractor|DLPFC_L" = (is.dlpfc * is.target - is.dlpfc * is.distractor) * is.left,
  "target-distractor|LPPC_L" = (is.lppc * is.target - is.lppc * is.distractor) * is.left,
  "target-distractor|DMFC_L" = (is.dmfc * is.target - is.dmfc * is.distractor) * is.left,
  "target-distractor|DLPFC_R" = (is.dlpfc * is.target - is.dlpfc * is.distractor) * is.right,
  "target-distractor|LPPC_R" = (is.lppc * is.target - is.lppc * is.distractor) * is.right,
  "target-distractor|DMFC_R" = (is.dmfc * is.target - is.dmfc * is.distractor) * is.right,
  
  "target-incongruency|DLPFC_L" = (is.dlpfc * is.target - is.dlpfc * is.incongruency) * is.left,
  "target-incongruency|LPPC_L" = (is.lppc * is.target - is.lppc * is.incongruency) * is.left,
  "target-incongruency|DMFC_L" = (is.dmfc * is.target - is.dmfc * is.incongruency) * is.left,
  "target-incongruency|DLPFC_R" = (is.dlpfc * is.target - is.dlpfc * is.incongruency) * is.right,
  "target-incongruency|LPPC_R" = (is.lppc * is.target - is.lppc * is.incongruency) * is.right,
  "target-incongruency|DMFC_R" = (is.dmfc * is.target - is.dmfc * is.incongruency) * is.right,

  "distractor-incongruency|DLPFC_L" = (is.dlpfc * is.distractor - is.dlpfc * is.incongruency) * is.left,
  "distractor-incongruency|LPPC_L" = (is.lppc * is.distractor - is.lppc * is.incongruency) * is.left,
  "distractor-incongruency|DMFC_L" = (is.dmfc * is.distractor - is.dmfc * is.incongruency) * is.left,
  "distractor-incongruency|DLPFC_R" = (is.dlpfc * is.distractor - is.dlpfc * is.incongruency) * is.right,
  "distractor-incongruency|LPPC_R" = (is.lppc * is.distractor - is.lppc * is.incongruency) * is.right,
  "distractor-incongruency|DMFC_R" = (is.dmfc * is.distractor - is.dmfc * is.incongruency) * is.right,
  
  ## cross-region contrasts:
  
  "DMFC_L-DLPFC|incongruency" = (is.dmfc * is.left - is.dlpfc / 2) * is.incongruency,
  "DMFC_L-LPPC|incongruency"  = (is.dmfc * is.left - is.lppc / 2) * is.incongruency,
  "DMFC_L-DLPFC|target"       = (is.dmfc * is.left - is.dlpfc / 2) * is.target,
  "DMFC_L-LPPC|target"        = (is.dmfc * is.left - is.lppc / 2) * is.target,
  
  ## interactions:
  
  "(DLPFC-DMFC_L)(target-incongruency)" = 
    (is.dlpfc/2 - (is.dmfc * is.left)) * is.target - (is.dlpfc/2 - (is.dmfc * is.left)) * is.incongruency,
  "(DLPFC-LPPC_L)(target-incongruency)" = 
    (is.dlpfc/2 - (is.lppc * is.left)) * is.target - (is.dlpfc/2 - (is.lppc * is.left)) * is.incongruency
  
  
)

(glht.super <- summary(glht(fit.super, contrasts.super), test = adjusted("none")))


## correct p-values from within-region contrats:

contrs2correct <- c(

  "target-distractor\\|DLPFC_.",
  "target-distractor\\|LPPC_.",
  "target-distractor\\|DMFC_.",
  "target-incongruency\\|DLPFC_.",
  "target-incongruency\\|LPPC_.",
  "target-incongruency\\|DMFC_.",
  "distractor-incongruency\\|DLPFC_.",
  "distractor-incongruency\\|LPPC_.",
  "distractor-incongruency\\|DMFC_."
  
)

glht.super$test$pvalues.fdr <- glht.super$test$pvalues
for (contr.i in seq_along(contrs2correct)) {
  
  p <- glht.super$test$pvalues[grep(contrs2correct[contr.i], names(glht.super$test$pvalues))]
  glht.super$test$pvalues.fdr[grep(contrs2correct[contr.i], names(glht.super$test$pvalues.fdr))] <- p.adjust(p, "fdr")
  
}

#+



#' ### build table
#+ fpc-dissoc_table_primary

table.group.contrast <- c(
  
  "$\\text{target DLPFC (L)}$",
  "$\\text{target LPPC (L)}$",
  "$\\text{target DMFC (L)}$",
  "$\\text{target DLPFC (R)}$",
  "$\\text{target LPPC (R)}$",
  "$\\text{target DMFC (R)}$",
  "$\\text{distr. DLPFC (L)}$",
  "$\\text{distr. LPPC (L)}$",
  "$\\text{distr. DMFC (L)}$",
  "$\\text{distr. DLPFC (R)}$",
  "$\\text{distr. LPPC (R)}$",
  "$\\text{distr. DMFC (R)}$",
  "$\\text{incon. DLPFC (L)}$",
  "$\\text{incon. LPPC (L)}$",
  "$\\text{incon. DMFC (L)}$",
  "$\\text{incon. DLPFC (R)}$",
  "$\\text{incon. LPPC (R)}$",
  "$\\text{incon. DMFC (R)}$",
  
  "$\\text{target}-\\text{distr.: } \\text{DLPFC}$",
  "$\\text{target}-\\text{distr.: } \\text{LPPC}$",
  "$\\text{target}-\\text{distr.: } \\text{DMFC}$",
  "$\\text{target}-\\text{incon.: } \\text{DLPFC}$",
  "$\\text{target}-\\text{incon.: } \\text{LPPC}$",
  "$\\text{target}-\\text{incon.: } \\text{DMFC}$",
  "$\\text{distr.}-\\text{incon.: } \\text{DLPFC}$",
  "$\\text{distr.}-\\text{incon.: } \\text{LPPC}$",
  "$\\text{distr.}-\\text{incon.: } \\text{DMFC}$",
  
  "$\\text{target}-\\text{distr.: } \\text{DLPFC (L)}$",
  "$\\text{target}-\\text{distr.: } \\text{LPPC (L)}$",
  "$\\text{target}-\\text{distr.: } \\text{DMFC (L)}$",
  "$\\text{target}-\\text{distr.: } \\text{DLPFC (R)}$",
  "$\\text{target}-\\text{distr.: } \\text{LPPC (R)}$",
  "$\\text{target}-\\text{distr.: } \\text{DMFC (R)}$",
  
  "$\\text{target}-\\text{incon.: } \\text{DLPFC (L)}$",
  "$\\text{target}-\\text{incon.: } \\text{LPPC (L)}$",
  "$\\text{target}-\\text{incon.: } \\text{DMFC (L)}$",
  "$\\text{target}-\\text{incon.: } \\text{DLPFC (R)}$",
  "$\\text{target}-\\text{incon.: } \\text{LPPC (R)}$",
  "$\\text{target}-\\text{incon.: } \\text{DMFC (R)}$",
  
  "$\\text{distr.}-\\text{incon.: } \\text{DLPFC (L)}$",
  "$\\text{distr.}-\\text{incon.: } \\text{LPPC (L)}$",
  "$\\text{distr.}-\\text{incon.: } \\text{DMFC (L)}$",
  "$\\text{distr.}-\\text{incon.: } \\text{DLPFC (R)}$",
  "$\\text{distr.}-\\text{incon.: } \\text{LPPC (R)}$",
  "$\\text{distr.}-\\text{incon.: } \\text{DMFC (R)}$",

  "$\\text{incon.: } \\text{DMFC (L)}-\\text{DLPFC}$",
  "$\\text{incon.: } \\text{DMFC (L)}-\\text{LPPC}$",
  "$\\text{target: } \\text{DMFC (L)}-\\text{DLPFC}$",
  "$\\text{target: } \\text{DMFC (L)}-\\text{LPPC}$",

  "$[\\text{target}-\\text{incon.}]\\cdot[\\text{DLPFC}-\\text{DMFC (L)}]$",
  "$[\\text{target}-\\text{incon.}]\\cdot[\\text{DLPFC}-\\text{LPPC (L)}]$"
  
)

table.group <- data.frame(
  contrast = table.group.contrast,
  b  = glht.super$test$coefficients,
  se = glht.super$test$sigma,
  t  = glht.super$test$tstat,
  p  = glht.super$test$pvalues.fdr
)
rownames(table.group) <- NULL

table.group.means <- table.group[1:18, ]
table.group.wnregion <- table.group[19:45, ]
table.group.bnregion <- table.group[46:51, ]

kable(table.group.means)
kable(table.group.wnregion)
kable(table.group.bnregion)

# fwrite(table.group.means, here("out", "group", "superparcels_means.txt"))
# fwrite(table.group.wnregion, here("out", "group", "superparcels_wnregion.txt"))
# fwrite(table.group.bnregion, here("out", "group", "superparcels_bnregion.txt"))


#+


#' ### plot
#+ fpc-dissoc_plot_primary, fig.height = 7, fig.width = 9

set.seed(0)
means.super <- stats.subjs.super %>%

  mutate(roi.hemi = paste0(roi, "_", hemi)) %>%
  filter(roi.hemi %in% c("dlpfc_L", "dlpfc_R", "dmfc_L")) %>%

  group_by(roi, param, subj) %>%
  summarize(beta = mean(beta), .groups = "drop_last") %>%

  summarize(res = list(boot_mean_ci(beta)), .groups = "drop") %>%

  tidyr::unnest(cols = c(res))

means.super %<>% mutate(param = factor(as.factor(param), levels = c("incongruency", "target", "distractor")))
ymax <- max(means.super$ymax) + 0.015

p.means.super <-
  means.super %>%
  
  mutate(roi = relevel(factor(roi), ref = "dmfc")) %>%
  
  ggplot(aes(roi, y, color = param)) +
    
    geom_point(size = geom.point.size, position = position_dodge(width = 1/2)) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), position = position_dodge(width = 1/2), width = 0, size = geom.line.size) +
    
    scale_color_manual(values = colors.model) +
    lemon::coord_capped_cart(left = "both") +
    scale_x_discrete(labels = c("DMFC (L)", "DLPFC", "LPPC")) +
    scale_y_continuous(breaks = c(0, 0.2, 0.2)) +
    
    # annotate(
    #   geom = "text", x = 1.5, y = 0.14, label = "target", color = colors.model["target"],
    #   hjust = 0, vjust = 1, size = label.size, fontface = "bold"
    # ) +
    # annotate(
    #   geom = "text", x = 1.5, y = 0.16, label = "incon.", color = colors.model["incongruency"],
    #   hjust = 0, vjust = 1, size = label.size, fontface = "bold"
    # ) +
    # annotate(
    #   geom = "text", x = 1.5, y = 0.12, label = "distr.", color = colors.model["distractor"],
    #   hjust = 0, vjust = 1, size = label.size, fontface = "bold"
    # ) +
    # 
    # annotate(
    #   geom = "segment", color = "grey40", size = geom.line.size/2,
    #   y = ymax, yend = ymax,
    #   x = 0.825 + 1, xend = 1.175 + 1
    # ) +
    # annotate(
    #   geom = "segment", color = "grey40", size = geom.line.size/2,
    #   y = ymax - 0.005, yend = ymax - 0.005,
    #   x = 1 + 1, xend = 1.175 + 1
    # ) +
    # 
    # annotate(
    #   geom = "segment", color = "grey40", size = geom.line.size/2, 
    #   y = ymax, yend = ymax,
    #   x = 0.825, xend = 1.175
    # ) +
    # annotate(
    #   geom = "segment", color = "grey40", size = geom.line.size/2, 
    #   y = ymax - 0.005, yend = ymax - 0.005,
    #   x = 1, xend = 1.175
    # ) +
    # annotate(
    #   geom = "segment", color = "grey40", size = geom.line.size/2, 
    #   y = ymax - 0.01, yend = ymax - 0.01,
    #   x = 0.825, xend = 1
    # ) +
    
    # annotate(
    #   geom = "segment", color = "grey40", size = geom.line.size/2, 
    #   y = ymax, yend = ymax,
    #   x = 0.825 + 2, xend = 1.175 + 2
    # ) +
    # annotate(
    #   geom = "segment", color = "grey40", size = geom.line.size/2, 
    #   y = ymax - 0.005, yend = ymax - 0.005,
    #   x = 3, xend = 1.175 + 2
    # ) +
    
    labs(y = bquote("Model fit ("*bar(beta)*")"), x = "Region") +
    
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

p.means.super
# ggsave(here("out", "group", "crossplot_superparcels.pdf"), p.means.super, device = "pdf", width = 4, height = 3)


#+


#+
rm(stats.subjs.super)
rm(measure)
#+