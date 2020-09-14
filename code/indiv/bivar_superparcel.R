#' * create bivariate scatterplots:
#'    * fig-indiv panel A
#'    * fig-indiv-hyp-bivar
#' * create table of stats
#'    * table-indiv-hyp-sep
#'    * table-indiv-hyp-all


#+ bivar-superparcel_correlate

set.seed(0)

cors <- d.super %>%
  
  filter(is.analysis.group, id %in% hyps) %>%
  
  group_by(id) %>%
  
  summarize(
    r.line = cor(beta, stroop),
    r.rank = cor(beta, stroop, method = "spearman"),
    r2.line = r.line^2,
    r2.rank = r.rank^2
  )

## get confidence intervals

l <- d.super %>%
  
  filter(is.analysis.group, id %in% hyps) %>%
  group_by(id) %>%
  split(f = .$id) %>%
  purrr::map(~.[c("stroop", "beta")])

ci.line <- lapply(l, cor_ci) %>% bind_rows
ci.rank <- lapply(l, cor_ci, method = "spearman") %>% bind_rows

names(ci.line) %<>% paste0("_line")
names(ci.rank) %<>% paste0("_rank")

cors <- bind_cols(cors, ci.line, ci.rank)
cors$p.geq0_line <- 1 - cors$p.leq0_line
cors$p.geq0_rank <- 1 - cors$p.leq0_rank

kable(cors %>% arrange(-r2.line), digits = 2)

fwrite(cors, here("out", "indiv", "hyp_bivariate.txt"))
#+

#' ### plot bivariate relationships, write correlation stats
#+ bivar-superparcel_plot

p.allcors <- d.super %>%
  
  mutate(
    hemi   = gsub("(.*)_(L|R)", "\\2", roi),
    region = gsub("(.*)_(L|R)", "\\1", roi) %>% toupper
  ) %>%
  
  filter(region %in% c("DLPFC", "DMFC", "LPPC"), param %in% c("target", "incongruency")) %>%
  
  ggplot(aes(beta, stroop, fill = param, color = param)) +
  
  stat_boot_ci(n = 1E4, alpha = 0.3, color = "transparent") +
  geom_smooth(method = "lm", se = FALSE, size = geom.line.size/2) +
  geom_point(size = geom.point.size/3) +
  
  scale_fill_manual(values = colors.model) +
  scale_color_manual(values = colors.model) +
  
  facet_grid(vars(hemi), vars(region)) +
  
  coord_capped_cart(left = "both", bottom = "both") +
  scale_y_continuous(limits = c(0, 150)) +
  scale_x_continuous(limits = c(-0.45, 0.65), breaks = c(-0.4, 0, 0.6)) +
  
  theme(
    panel.grid      = element_blank(),
    panel.border    = element_blank(),
    # panel.background = element_rect(fill = "grey90"),
    # plot.margin     = unit(c(0, 0, 0, 0), "cm") ,
    axis.line       = element_line(size = axis.line.size),
    axis.text       = element_text(size = axis.text.size),
    axis.ticks      = element_line(size = axis.line.size),
    axis.title      = element_text(size = axis.title.size*1.5),
    strip.background = element_blank(),
    strip.text = element_text(size = axis.title.size),
    legend.position = "none"
  )

p.allcors

ggsave(here("out", "indiv", "all_cors.pdf"), p.allcors, width = 17.6, height = 17.6/2, units = "cm")


w.super %>%
  
  filter(is.analysis.group) %>% 
  
  {
    
    grid.arrange(
      
      ggplot(., aes(dlpfc_R_target, lppc_R_target)) +
        stat_boot_ci(n = 1E4, alpha = 0.3, fill = colors.model["target"]) +
        geom_point(fill = colors.model["target"], color = "white", shape = 21, size = 4) +
        ggtitle("DLPFC~LPPC (right) target"),
      
      ggplot(., aes(dlpfc_R_incongruency, lppc_R_incongruency)) +
        stat_boot_ci(n = 1E4, alpha = 0.3, fill = colors.model["incongruency"]) +
        geom_point(fill = colors.model["incongruency"], color = "white", shape = 21, size = 4) +
        ggtitle("DLPFC~LPPC (right) incongruency"),
      
      ncol = 2
      
    )
    
  }


cors.lfp <- w.super %>%
  filter(is.analysis.group) %>% 
  summarize(
    r_target = cor(dlpfc_R_target, lppc_R_target), 
    r_incongruency = cor(dlpfc_R_incongruency, lppc_R_incongruency)
  ) %>%
  as.data.table


fwrite(cors.lfp, here("out", "indiv", "cors_lfp.txt"))
#+

#' ### make LFP variable

#+ bivar-superparcel_combine

w.super$lfp_R_target <- (w.super$dlpfc_R_target + w.super$lppc_R_target) / 2
w.super$lfp_R_incongruency <- (w.super$dlpfc_R_incongruency + w.super$lppc_R_incongruency) / 2


d.dissoc.hlm$lfp_R_target <- (d.dissoc.hlm$dlpfc_R_target + d.dissoc.hlm$lppc_R_target) / 2
d.dissoc.hlm$lfp_R_incongruency <- (d.dissoc.hlm$dlpfc_R_incongruency + d.dissoc.hlm$lppc_R_incongruency) / 2

d.super %<>%
  
  {
    bind_rows(
      .,
      filter(., roi %in% c("dlpfc_R", "lppc_R"), param %in% c("target", "incongruency")) %>%
        group_by(subj, scheme, is.analysis.group, param) %>%
        summarize(
          congr = unique(congr), stroop = unique(stroop), 
          beta = mean(beta), beta.s = mean(beta.s)
        ) %>%
        ungroup %>%
        mutate(roi = "lfp_R", id = paste0(roi, "_", param))
    )
  }
#+

#+ plot-panel-a

set.seed(0)
cor.dmfc.l.incongruency <- cor_ci(w.super[c("dmfc_L_incongruency", "stroop")], R = 1E4)
cor.lfp.r.target <- cor_ci(w.super[c("lfp_R_target", "stroop")], R = 1E4)
cor.dmfc.l.incongruency.rank <- w.super[c("dmfc_L_incongruency", "stroop")] %>% mutate_all(rank) %>% cor_ci(R = 1E4)
cor.lfp.r.target.rank <- w.super[c("lfp_R_target", "stroop")] %>% mutate_all(rank) %>% cor_ci(R = 1E4)

w.super %>%
  
  filter(is.analysis.group) %>%
  
  ggplot() +
  
  stat_boot_ci(aes(lfp_R_target, stroop), n = 1E4, alpha = 0.3, fill = colors.model["target"]) +
  stat_smooth(aes(lfp_R_target, stroop), method = "lm", color = colors.model["target"], se = FALSE) +
  
  stat_boot_ci(aes(dmfc_L_incongruency, stroop), n = 1E4, alpha = 0.3, fill = colors.model["incongruency"]) +
  stat_smooth(aes(dmfc_L_incongruency, stroop), method = "lm", color = colors.model["incongruency"], se = FALSE) +
  
  geom_point(
    aes(lfp_R_target, stroop), 
    fill = colors.model["target"], color = "white", shape = 21, size = geom.point.size*2
  ) +
  geom_point(
    aes(dmfc_L_incongruency, stroop),
    fill = colors.model["incongruency"], color = "white", shape = 21, size = geom.point.size*2
  ) +
  
  annotate(
    geom = "text", x = 0.3, y = 40, 
    label = paste0(
      "MFC incon.:\nr = ", 
      round(cor.dmfc.l.incongruency$t0, 2), ", [", 
      round(cor.dmfc.l.incongruency$lower, 2), ", ", 
      round(cor.dmfc.l.incongruency$upper, 2), "]\n\U03C1 = ",
      
      round(cor.dmfc.l.incongruency.rank$t0, 2), ", [", 
      round(cor.dmfc.l.incongruency.rank$lower, 2), ", ", 
      round(cor.dmfc.l.incongruency.rank$upper, 2), "]"
    ),
    size = label.size,
    hjust = 0,
    fontface = "italic",
    color = colors.model["incongruency"]
  ) +
  
  annotate(
    geom = "text", x = -0.35, y = 140, 
    label = paste0(
      "LFP target:\nr = ", 
      round(cor.lfp.r.target$t0, 2), ", [", 
      round(cor.lfp.r.target$lower, 2), ", ", 
      round(cor.lfp.r.target$upper, 2), "]\n\U03C1 = ",
      
      round(cor.lfp.r.target.rank$t0, 2), ", [", 
      round(cor.lfp.r.target.rank$lower, 2), ", ", 
      round(cor.lfp.r.target.rank$upper, 2), "]"
    ),
    size = label.size,
    hjust = 0,
    fontface = "italic",
    color = colors.model["target"]
  ) +
  
  theme_bw(base_size = 8) +
  
  coord_capped_cart(left = "both", bottom = "both") +
  
  theme(
    panel.grid      = element_blank(), 
    panel.border    = element_blank(),
    # plot.margin     = unit(c(0, 0, 0, 0), "cm") ,
    axis.line       = element_line(size = axis.line.size),
    axis.text       = element_text(size = axis.text.size),
    axis.ticks      = element_line(size = axis.line.size),
    axis.title      = element_text(size = axis.title.size*1.5)
  ) +
  
  labs(y = "Stroop effect (RT)", x = bquote("Model fit ("*beta*")"))

p.dissoc <- last_plot()

#+