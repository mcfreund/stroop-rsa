#+ fpc-dissoc-altdef_setup, include = FALSE

if (interactive()) source(here::here("code", "group", "fpc_dissoc.R"))

#' ### model and table
#+ fpc-dissoc-altdef_model, fig.height = 7, fig.width = 9

fit.super.alt <- update(
  fit.super,
  . ~ . -  (roi * param | subj) + (roi + param | subj), 
  data = stats.subjs.super %>% filter(roi %in% c("dlpfc.alt", "dmfc.alt", "lppc")),
  control = lmerControl(optimizer = "bobyqa")
)
summary(fit.super.alt)

glht.super.alt <- summary(glht(fit.super.alt, contrasts.super), test = adjusted("none"))


glht.super.alt$test$pvalues.fdr <- glht.super.alt$test$pvalues
for (contr.i in seq_along(contrs2correct)) {
  
  p <- glht.super.alt$test$pvalues[grep(contrs2correct[contr.i], names(glht.super.alt$test$pvalues))]
  glht.super.alt$test$pvalues.fdr[grep(contrs2correct[contr.i], names(glht.super.alt$test$pvalues.fdr))] <- p.adjust(p, "fdr")
  
}

## table


table.group.alt <- data.frame(
  contrast = table.group.contrast,
  b  = glht.super.alt$test$coefficients,
  se = glht.super.alt$test$sigma,
  t  = glht.super.alt$test$tstat,
  p  = glht.super.alt$test$pvalues.fdr
)
table.group.alt <- table.group.alt[grep("DLPFC|DMFC", table.group.alt$contrast), ]  ## subset important rows
rownames(table.group.alt) <- NULL

kable(table.group.alt)

fwrite(table.group.alt, here("out", "group", "superparcels_alt.txt"))

# table.group.alt.means <- table.group.alt[1:12, ]
# table.group.alt.wnregion <- table.group.alt[13:30, ]
# table.group.alt.bnregion <- table.group.alt[31:51, ]
# 
# kable(table.group.alt.means)
# kable(table.group.alt.wnregion)
# kable(table.group.alt.bnregion)
# 
# fwrite(table.group.alt.means, here("out", "group", "superparcels_alt_means.txt"))
# fwrite(table.group.alt.wnregion, here("out", "group", "superparcels_alt_wnregion.txt"))
# fwrite(table.group.alt.bnregion, here("out", "group", "superparcels_alt_bnregion.txt"))




#+ 

#' ### plot
#+ fpc-dissoc-altdef_plot, fig.height = 7, fig.width = 9

set.seed(0)
means.super.alt <- stats.subjs.super %>%
  
  mutate(roi.hemi = paste0(roi, "_", hemi)) %>%
  filter(roi.hemi %in% c("dlpfc.alt_L", "dlpfc.alt_R", "dmfc.alt_L", "lppc_L")) %>%
  
  group_by(roi, param, subj) %>%
  summarize(beta = mean(beta)) %>% 
  
  group_by(roi, param) %>%
  summarize(res = list(boot_mean_ci(beta))) %>% 
  
  tidyr::unnest(cols = c(res))

p.means.super.alt <- means.super.alt %>%  
  
  ungroup %>%
  mutate(
    roi = gsub(".alt", "", roi),
    param = factor(as.factor(param), levels = c("incongruency", "target", "distractor"))
  ) %>%
  
  ggplot(aes(param, y, group = roi, color = roi)) +
  geom_point(size = geom.point.size, position = position_dodge(width = 1/2)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position = position_dodge(width = 1/2), width = 0, size = geom.line.size) +
  geom_line(size = geom.line.size*3/4, position = position_dodge(width = 1/2)) +
  
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

p.means.super.alt

ggsave(here("out", "group", "crossplot_superparcels_alt.pdf"), p.means.super.alt, device = "pdf", width = 4, height = 3)

