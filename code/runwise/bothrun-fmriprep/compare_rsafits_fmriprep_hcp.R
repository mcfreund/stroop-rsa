source(here::here("code", "packages.R"))
source(here("code", "strings.R"))

stats.subjs.fmriprep <- fread(
  here("out", "rsa", "stats",  paste0("subjs_pro_bias_acc-only_fmriprep_mmp_residual.csv"))
)
stats.subjs.hcp <- fread(
  here("out", "rsa", "stats",  paste0("subjs_pro_bias_acc-only_mmp_residual.csv"))
)

subjs <- intersect(unique(stats.subjs.hcp$subj), unique(stats.subjs.fmriprep$subj))

stats.subjs.hcp <- stats.subjs.hcp[subj %in% subjs]
stats.subjs.fmriprep <- stats.subjs.fmriprep[subj %in% subjs]

d_subj <- full_join(
  stats.subjs.fmriprep, 
  stats.subjs.hcp,
  by = c("subj", "is.analysis.group", "param", "roi"),
  suffix = c("_fmriprep", "_hcp")
)
d_subj %>%
  
  group_by(param, roi) %>%
  
  summarize(r = cor(beta_fmriprep, beta_hcp)) %>%
  
  ggplot(aes(param, r)) +
  geom_boxplot()




stats.group.hcp <- stats.subjs.hcp %>%
  group_by(param, roi) %>%
  summarize(
    b = mean(beta),
    p = wilcox.test(beta, alternative = "greater")$p.value
  ) %>%
  ungroup %>%
  mutate(p.fdr = p.adjust(p, "fdr"))


stats.group.fmriprep <- stats.subjs.fmriprep %>%
  group_by(param, roi, .drop = "last") %>%
  summarize(
    b = mean(beta),
    p = wilcox.test(beta, alternative = "greater")$p.value
  ) %>%
  ungroup %>%
  mutate(p.fdr = p.adjust(p, "fdr"))

stats.group.fmriprep %>% filter(param == "target", p.fdr < 0.05) %>% nrow
stats.group.hcp %>% filter(param == "target", p.fdr < 0.05) %>% nrow

stats.group.fmriprep %>% filter(param == "distractor", p.fdr < 0.05) %>% nrow
stats.group.hcp %>% filter(param == "distractor", p.fdr < 0.05) %>% nrow

stats.group.fmriprep %>% filter(param == "incongruency", p.fdr < 0.05) %>% nrow
stats.group.hcp %>% filter(param == "incongruency", p.fdr < 0.05) %>% nrow

stats.group.fmriprep %>% 
  
  d <- full_join(
    stats.group.fmriprep, 
    stats.group.hcp,
    by = c("param", "roi"),
    suffix = c("_fmriprep", "_hcp")
  )


d %>%
  group_by(param) %>%
  summarize(r = cor(b_fmriprep, b_hcp))

d %>%
  ggplot(aes(b_hcp, b_fmriprep)) +
  geom_abline() +
  geom_point(alpha = 0.5) +
  facet_grid(cols = vars(param)) +
  theme_minimal()


intersect(
  d %>% filter(param == "target", p.fdr_fmriprep < 0.01) %>% pull(roi),
  d %>% filter(param == "target", p.fdr_hcp < 0.01) %>% pull(roi)
)

setdiff(
  d %>% filter(param == "target", p.fdr_fmriprep < 0.01) %>% pull(roi),
  d %>% filter(param == "target", p.fdr_hcp < 0.01) %>% pull(roi)
)

setdiff(
  d %>% filter(param == "target", p.fdr_hcp < 0.01) %>% pull(roi),
  d %>% filter(param == "target", p.fdr_fmriprep < 0.01) %>% pull(roi)
)




