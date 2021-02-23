#+ echo = FALSE
stats.subjs <- data.super[[measure]]

stats.subjs <- stats.subjs[stats.subjs$param %in% c("target", "distractor", "incongruency"), ]
stats.subjs <- full_join(blups, stats.subjs, by = "subj")
stats.subjs %<>% ungroup %>% mutate(id = paste0(roi.hemi, "_", param))
stats.subjs %<>% group_by(id) %>% mutate(beta.s = c(scale(beta)))  ## scale betas

## to long-form data-frame

d.super <- stats.subjs[stats.subjs$subj %in% subjs.analysis, ]

#+

#' N subj: `r length(unique(d.super$subj))`

#+ echo = FALSE

set.seed(0)

cors <- d.super %>%
  
  filter(id %in% hyps) %>%
  
  group_by(id) %>%
  
  summarize(
    r.line = cor(beta, stroop),
    r.rank = cor(beta, stroop, method = "spearman"),
    r2.line = r.line^2,
    r2.rank = r.rank^2
  )

## get confidence intervals

l <- d.super %>%
  
  filter(id %in% hyps) %>%
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

#+

#+ echo = FALSE, fig.width = 12, fig.height = 5

p.allcors <- d.super %>%
  
  filter(id %in% hyps) %>%
  
  mutate(region = toupper(roi)) %>%
  
  filter(region %in% c("DLPFC", "DMFC", "LPPC"), param %in% c("target", "incongruency")) %>%
  
  ggplot(aes(beta, stroop, fill = param, color = param)) +
  
  stat_boot_ci(n = 1E4, alpha = 0.3, color = "transparent") +
  geom_smooth(method = "lm", se = FALSE, size = geom.line.size/2) +
  geom_point(size = geom.point.size) +
  
  scale_fill_manual(values = colors.model) +
  scale_color_manual(values = colors.model) +
  
  facet_grid(vars(hemi), vars(region)) +
  
  coord_capped_cart(left = "both", bottom = "both") +
  scale_y_continuous(limits = c(0, 150)) +
  # scale_x_continuous(limits = c(-0.45, 0.65), breaks = c(-0.4, 0, 0.6)) +
  
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

#+
