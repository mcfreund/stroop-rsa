#+ fpc-dissoc-altdef_setup, include = FALSE

if (interactive()) source(here::here("code", "group", "mds.R"))

#' ### model and table
#+ fpc-dissoc-altdef_model, fig.height = 7, fig.width = 9

fit.super.alt <- update(
  fit.super,
  . ~ . -  (roi * param | subj) + (roi + param | subj), 
  data = stats.subjs.super %>% filter(roi %in% c("dlpfc.alt", "dmfc.alt", "lppc"))
)
summary(fit.super.alt)

glht.super.all.alt <- summary(glht(fit.super.alt, contrasts.super.all), test = adjusted("none"))

p.fdr.tvi.dlpfc.super.alt <- glht.super.all.alt$test$pvalues[grep("dlpfc_", rownames(contrasts.super))] %>% 
  p.adjust(method = "fdr")
p.fdr.tvi.lppc.super.alt <- glht.super.all.alt$test$pvalues[grep("lppc_", rownames(contrasts.super))] %>% 
  p.adjust(method = "fdr")
p.fdr.tvi.dmfc.super.alt <- glht.super.all.alt$test$pvalues[grep("dmfc_", rownames(contrasts.super))] %>% 
  p.adjust(method = "fdr")

## table ----

table.group.alt <- data.frame(
  contrast = table.group.contrast,
  type = table.group.type,
  b        = glht.super.all.alt$test$coefficients,
  se       = glht.super.all.alt$test$sigma,
  t        = glht.super.all.alt$test$tstat,
  p        = glht.super.all.alt$test$pvalues
)
table.group.alt$p.raw <- table.group.alt$p
table.group.alt[
  table.group.alt$contrast %in% 
    c(names(p.fdr.tvi.dlpfc.super.alt), names(p.fdr.tvi.lppc.super.alt), names(p.fdr.tvi.dmfc.super.alt)),
  "p"
  ] <- c(p.fdr.tvi.dlpfc.super.alt, p.fdr.tvi.lppc.super.alt, p.fdr.tvi.dmfc.super.alt)
rownames(table.group.alt) <- NULL
kable(table.group.alt)

fwrite(table.group.alt, here("out", "group", "superparcels_alt.txt"))

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

