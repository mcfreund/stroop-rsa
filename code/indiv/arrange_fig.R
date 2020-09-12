#+ arrange-fig
set.seed(0)
cor.vvis.l.incongruency <- cor_ci(w.super[c("vvis_L_incongruency", "stroop")], R = 1E4)
cor.vvis.l.incongruency.rank <- w.super[c("vvis_L_incongruency", "stroop")] %>% mutate_all(rank) %>% cor_ci(R = 1E4)
cor.dlpfc.l.distractor <- cor_ci(w.super[c("dlpfc_L_distractor", "stroop")], R = 1E4)
cor.dlpfc.l.distractor.rank <- w.super[c("dlpfc_R_target", "stroop")] %>% mutate_all(rank) %>% cor_ci(R = 1E4)

w.super %>%
  
  filter(is.analysis.group) %>%
  
  ggplot() +
  
  stat_boot_ci(aes(vvis_L_incongruency, stroop), n = 1E4, alpha = 0.3, fill = colors.model["incongruency"]) +
  stat_smooth(aes(vvis_L_incongruency, stroop), color = colors.model["incongruency"], method = "lm", se = FALSE) +
  geom_point(
    aes(vvis_L_incongruency, stroop), 
    fill = colors.model["incongruency"], color = "white", shape = 21, size = geom.point.size
  ) +
  
  coord_capped_cart(left = "both", bottom = "both") +
  
  annotate(
    geom = "text", x = 0.1, y = 140, 
    label = paste0(
      "r = ", 
      round(cor.vvis.l.incongruency$t0, 2), ", [", 
      round(cor.vvis.l.incongruency$lower, 2), ", ", 
      round(cor.vvis.l.incongruency$upper, 2), "]\n\U03C1 = ",
      
      round(cor.vvis.l.incongruency.rank$t0, 2), ", [", 
      round(cor.vvis.l.incongruency.rank$lower, 2), ", ", 
      round(cor.vvis.l.incongruency.rank$upper, 2), "]"
    ),
    size = label.size/2,
    hjust = 0,
    fontface = "italic",
    color = colors.model["incongruency"]
  ) +
  
  theme_bw(base_size = 8) +
  
  theme(
    panel.grid      = element_blank(), 
    panel.border    = element_blank(),
    # plot.margin     = unit(c(0, 0, 0, 0), "cm") ,
    axis.line       = element_line(size = axis.line.size),
    axis.text       = element_text(size = axis.text.size),
    axis.ticks      = element_line(size = axis.line.size),
    axis.title      = element_text(size = axis.title.size)
  ) +
  
  labs(y = "Stroop effect (RT)", x = bquote("Vent. Vis. (L) "*beta["incon."]*""))

p.vvis <- last_plot()

w.super %>%
  
  filter(is.analysis.group) %>%
  
  ggplot() +
  
  stat_boot_ci(aes(dlpfc_L_distractor, stroop), n = 1E4, alpha = 0.3, fill = colors.model["distractor"]) +
  stat_smooth(aes(dlpfc_L_distractor, stroop), color = colors.model["distractor"], method = "lm", se = FALSE) +
  
  geom_point(
    aes(dlpfc_L_distractor, stroop), 
    fill = colors.model["distractor"], color = "white", shape = 21, size = geom.point.size
  ) +
  
  theme_bw(base_size = 8) +
  
  
  
  annotate(
    geom = "text", x = -0.135, y = 35, 
    label = paste0(
      "r = ", 
      round(cor.dlpfc.l.distractor$t0, 2), ", [", 
      round(cor.dlpfc.l.distractor$lower, 2), ", ", 
      round(cor.dlpfc.l.distractor$upper, 2), "]\n\U03C1 = ",
      
      round(cor.dlpfc.l.distractor.rank$t0, 2), ", [", 
      round(cor.dlpfc.l.distractor.rank$lower, 2), ", ", 
      round(cor.dlpfc.l.distractor.rank$upper, 2), "]"
    ),
    size = label.size/2,
    hjust = 0,
    fontface = "italic",
    # nudge_x = -0.1,
    color = colors.model["distractor"]
  ) +
  
  coord_capped_cart(left = "both", bottom = "both") +
  
  theme(
    panel.grid      = element_blank(), 
    panel.border    = element_blank(),
    # plot.margin     = unit(c(0, 0, 0, 0), "cm") ,
    axis.line       = element_line(size = axis.line.size),
    axis.text       = element_text(size = axis.text.size),
    axis.ticks      = element_line(size = axis.line.size),
    axis.title      = element_text(size = axis.title.size),
    axis.title.y    = element_blank(),
    axis.text.y     = element_blank()
  ) +
  
  labs(y = "Stroop effect (RT)", x = bquote("DLPFC (L) "*beta["distr."]*""))

p.dlpfc <- last_plot()

# plot_grid(p.vvis, p.dlpfc, axis = "")
# p.explor <- plot_grid(p.vvis, p.dlpfc, align = "v")

data.frame(x = cor.perm.super[!is.na(cor.perm.super)]) %>%
  
  ggplot(aes(x = x)) +
  geom_density(fill = "grey40", color = "black", size = geom.line.size) +
  geom_vline(xintercept = cor.obs.super, size = geom.line.size, color = "firebrick1") +
  
  coord_flip() +
  
  theme_bw(base_size = 8) +
  
  theme(
    panel.grid      = element_blank(), 
    panel.border    = element_blank(),
    axis.line       = element_line(size = axis.line.size),
    axis.text       = element_text(size = axis.text.size),
    axis.ticks      = element_line(size = axis.line.size),
    axis.title      = element_text(size = axis.title.size),
    axis.title.x    = element_blank(),
    axis.text.x     = element_blank(),
    axis.line.x     = element_blank(),
    axis.ticks.x    = element_blank()
  ) +
  
  annotate(
    geom = "text", y = 0, x = 0, label = "null", color = "black", size = label.size, fontface = "italic",
    hjust = 0
  ) +
  
  annotate(
    geom = "text", y = 0, x = cor.obs.super, vjust = 0.5, hjust = 0,
    label = paste0("r = ", round(cor.obs.super, 2), "\np = ", round(p.super, 4)), 
    color = "firebrick1", size = label.size, fontface = "italic"
  ) +
  
  labs(x = "Predictive accuracy")

p.heldout <- last_plot()


## arrange ----

ratio <- (1 + sqrt(5))/2
labelmargin <- 0.1
fig.width <- 11.6
unit.height <- fig.width/ratio  ## cm
fig.height <- unit.height*3/2 + unit.height*labelmargin*2

bottom_row <- plot_grid(
  p.vvis, p.dlpfc, p.heldout, nrow = 1, rel_widths = c(1, 5/6, 2/3),
  labels = c("B", "", "C"),
  vjust = c(0, 0, 0)
)

plot.indiv <- plot_grid(
  NULL, p.dissoc, NULL, bottom_row,
  nrow = 4,
  rel_heights = c(labelmargin, 2, labelmargin, 1), 
  labels = c("", "A"),
  vjust = c(0, 0.3, 0)
)

plot.indiv

ggsave(
  here("out", "indiv", "fig_indiv.pdf"),
  plot.indiv,
  device = "pdf", height = fig.height, width = fig.width, unit = "cm"
)
#+