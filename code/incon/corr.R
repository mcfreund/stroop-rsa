source(here::here("code", "packages.R"))
source(here("code", "strings.R"))
source(here("code", "funs.R"))
source(here("code", "plots.R"))
source(here("code", "read_atlases.R"))

## get correlations among fits ----

stats.subj.cc <- read_subj_stats(suffix = "_intercept-cc_residual", params = "incongruency")
stats.subj <- read_subj_stats(params = c("congruency", "incongruency"))
stats.subj <- stats.subj %>% 
  pivot_wider(names_from = c("param"), values_from = c("beta")) %>%
  mutate(
    beta_cc = incongruency - congruency,
    beta_ic = incongruency
    ) %>%
  dplyr::select(-congruency, -incongruency)

stats.subj <- full_join(
  stats.subj, stats.subj.cc %>% select(-param, beta_iccc = beta), 
  by = c("subj", "is.analysis.group", "roi", "roi.hemi", "hemi")
  )


group.cor <- stats.subj %>%
  filter(roi %in% c("dlpfc", "lppc", "dmfc")) %>%
  group_by(roi.hemi) %>%
  summarize(
    r_iccc = cor(beta_ic, beta_iccc),
    r_cc = cor(beta_ic, beta_cc)
    )


p.roi <-
  
  stats.subj %>%
  filter(roi %in% c("dlpfc", "lppc", "dmfc")) %>%
  pivot_longer(cols = c("beta_cc", "beta_ic", "beta_iccc"), names_to = "modtype") %>%
    
  ggplot(aes(modtype, value)) +
    
  geom_line(aes(group = subj), color = "grey80") +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, size = 1, fun.args = list(B = 1E4)) +
  facet_grid(rows = vars(hemi), cols = vars(roi), labeller = labeller(roi = toupper)) +
  geom_hline(yintercept = 0) +
  
  scale_x_discrete(labels = c("II-CC", "II-IC\n(original)", "II-(CC+IC)")) +
  labs(x = "Incongruency model parameterization", y = bquote("Model fit ("*bar(beta)*")")) +
  theme_minimal() +
  lemon::coord_capped_cart(left = "both") +
  theme(
    legend.position = "none", 
    panel.grid      = element_blank(), 
    panel.border    = element_blank(),
    axis.line.y     = element_line(size = axis.line.size),
    axis.text       = element_text(size = axis.text.size*0.75),
    axis.ticks.y    = element_line(size = axis.line.size),
    axis.ticks.x    = element_blank(),
    axis.title    = element_text(size = axis.title.size)
  )


## now for all MMP ---

stats.subj.cc <- read_subj_stats(suffix = "_intercept-cc_residual", roi.set = "mmp", param = c("incongruency"))
stats.subj <- read_subj_stats(roi.set = "mmp", param = c("congruency", "incongruency"))
stats.subj <- stats.subj %>% 
  pivot_wider(names_from = c("param"), values_from = c("beta")) %>%
  mutate(
    beta_cc = incongruency - congruency,
    beta_ic = incongruency
  ) %>%
  dplyr::select(-congruency, -incongruency)


stats.subj <- full_join(
  stats.subj, stats.subj.cc %>% select(-param, beta_iccc = beta), 
  by = c("subj", "is.analysis.group", "roi", "num.roi", "hemi", "community", "community.short", "network", "md"),
  suffix = c("_ic", "_iccc")
  )

group.stats.mmp <- stats.subj %>%
  group_by(roi) %>%
  summarize(
    p_ic = wilcox.test(beta_ic, alternative = "greater")$p.value,
    b_ic = mean(beta_ic),
    p_cc = wilcox.test(beta_cc, alternative = "greater")$p.value,
    b_cc = mean(beta_cc),
    p_iccc = wilcox.test(beta_iccc, alternative = "greater")$p.value,
    b_iccc = mean(beta_iccc)
  ) %>%
  mutate(p_ic = p.adjust(p_ic, "fdr"), p_cc = p.adjust(p_cc, "fdr"), p_iccc = p.adjust(p_iccc, "fdr"))

group.stats.mmp %>% filter(p_iccc < 0.05)
group.stats.mmp %>% filter(p_ic < 0.05)
group.stats.mmp %>% filter(p_cc < 0.05)


p.mmp.ic <-
  
  group.stats.mmp %>%
  ggplot(aes(b_iccc, b_ic)) +
  geom_abline() +
  geom_point(shape = 21, fill = "black", color = "white", size = 2) +
    
  theme_minimal() +
  lemon::coord_capped_cart(left = "both", bottom = "both") +
  theme(
    legend.position = "none", 
    panel.grid      = element_blank(), 
    panel.border    = element_blank(),
    axis.line     = element_line(size = axis.line.size),
    axis.text       = element_text(size = axis.text.size*0.75),
    axis.ticks    = element_line(size = axis.line.size),
    axis.title    = element_text(size = axis.title.size)
  ) +
  labs(y = "II-IC (original)", x = "II-(IC+CC)")


p.mmp.cc <- 
  
  group.stats.mmp %>%
  ggplot(aes(b_cc, b_ic)) +
  geom_abline() +
  geom_point(shape = 21, fill = "black", color = "white", size = 2) +

  theme_minimal() +
  lemon::coord_capped_cart(left = "both", bottom = "both") +
  theme(
    legend.position = "none", 
    panel.grid      = element_blank(), 
    panel.border    = element_blank(),
    axis.line     = element_line(size = axis.line.size),
    axis.text       = element_text(size = axis.text.size*0.75),
    axis.ticks    = element_line(size = axis.line.size),
    axis.title    = element_text(size = axis.title.size)
  ) +
  labs(y = "II-IC (original)", x = "II-CC")

  
group.cor <- stats.subj %>%
  group_by(roi) %>%
  summarize(
    r_ic = cor(beta_ic, beta_iccc),
    r_cc = cor(beta_ic, beta_cc)
  )

range(group.cor$r_ic)
range(group.cor$r_cc)


R <- diag(16)
colnames(R) <- bias.items
rownames(R) <- bias.items
R[bias.items.incon, bias.items.con] <- 1
R[bias.items.con, bias.items.incon] <- 1
sum(R[lower.tri(R)]) / ((16^2-16)/2)
sum(R[lower.tri(R)])
(16^2-16)/2


p.roi
p.mmp <- plot_grid(p.mmp.cc, p.mmp.ic)


## save ----

ggsave(here("out", "incon", "group_roi.pdf"), p.roi, device = "pdf", width = 16, height = 9, unit = "cm")
ggsave(here("out", "incon", "group_mmp.pdf"), p.mmp, device = "pdf", width = 16, height = 9, unit = "cm")

p <- plot_grid(p.roi, p.mmp, ncol = 1, rel_heights = c(1, 0.75), labels = "AUTO")

ggsave(here("out", "incon", "group_incon.pdf"), p, device = "pdf", width = 16, height = 18, unit = "cm")
