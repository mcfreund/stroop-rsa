#+ visual-sm-dissoc_setup, include = FALSE

if (interactive()) {
  
  source(here::here("code", "packages.R"))
  source(here("code", "strings.R"))
  source(here("code", "funs.R"))
  source(here("code", "plots.R"))
  source(here("code", "read_atlases.R"))
  
  stats.subjs.super <- read_subj_stats(roi.set = "masks")
  
}

#' ### model
#+ visual-sm-dissoc-model, include = TRUE

fit.prelim <- lmer(
  beta ~ 0 + interaction(roi, param) + (1 | subj),
  stats.subjs.super %>% filter(param %in% c("distractor", "target"), roi %in% c("V1", "smmouth"))
)
summary(fit.prelim)

contrasts.prelim <- rbind(
  diag(4),  #' means
  ## cross-parcel contrasts:
  "V1-smmouth|distractor" = c(1, -1, 0, 0),
  "smmouth-V1|target"     = c(0, 0, 1, -1),
  ## within-parcel contrasts:
  "distractor-target|V1"     = c(0, 1, 0, -1),
  "target-distractor|smmouth"= c(-1, 0, 1, 0),
  ## interaction:
  "(smmouth-V1)(target-distractor)" = c(-1, -1, 1, 1)
)
rownames(contrasts.prelim)[1:4] <- c("V1.distractor", "smmouth.distractor", "V1.target", "smmouth.target")

glht.prelim <- summary(glht(fit.prelim, contrasts.prelim, alternative = "greater"), test = adjusted("none"))
glht.prelim$test$pvalues["(smmouth-V1)(target-distractor)"] <- 
  glht.prelim$test$pvalues["(smmouth-V1)(target-distractor)"]*2  ## interaction should be two-sided

glht.prelim


#' ### plot
#+ visual-sm-dissoc-plot, include = TRUE, fig.width = 6, fig.height = 4


means.prelim <- stats.subjs.super %>% 
  
  filter(param %in% c("distractor", "target"), roi %in% c("evis", "smmouth")) %>%
  
  group_by(roi, param) %>%
  summarize(res = list(boot_mean_ci(beta)), .groups = "drop") %>% 
  
  unnest(cols = c(res))

p.means.prelim <- means.prelim %>%
  
  ggplot(aes(param, y, group = roi, color = roi)) +
  
  geom_point(size = geom.point.size, position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position = position_dodge(width = 0.25), width = 0, size = geom.line.size) +
  geom_line(size = geom.line.size, position = position_dodge(width = 0.25)) +
  
  scale_color_manual(values = c(evis = "#b2182b", smmouth = "#2166ac")) +
  coord_capped_cart(left = "both") +
  
  annotate(
    geom = "text", x = 0.5, y = 0.12, label = "ventral S1/M1", color = "#2166ac",
    hjust = 0, vjust = 1, size = label.size*1.5, fontface = "bold"
  ) +
  annotate(
    geom = "text", x = 0.5, y = 0.09, label = "V1", color = "#b2182b",
    hjust = 0, vjust = 1, size = label.size*1.5, fontface = "bold"
  ) +
  
  annotate(
    geom = "text", x = 0.5, y = 0.08,
    label = paste0("p(distr.-targt.)>0 = ", format.pval(glht.prelim$test$pvalues["distractor-target|V1"])),
    color = "#b2182b",
    hjust = 0, vjust = 1, size = label.size, fontface = "bold"
  ) +
  
  annotate(
    geom = "text", x = 0.5, y = 0.11,
    label = paste("p(targt.-distr.)>0 ", format.pval(glht.prelim$test$pvalues["target-distractor|smmouth"])),
    color = "#2166ac",
    hjust = 0, vjust = 1, size = label.size, fontface = "bold"
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

p.means.prelim

ggsave(here("out", "group", "prelim.pdf"), p.means.prelim, device = "pdf", width = figwidth, height = 3)

