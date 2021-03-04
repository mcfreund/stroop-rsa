source(here::here("code", "packages.R"))
source(here("code", "strings.R"))
source(here("code", "funs.R"))
source(here("code", "plots.R"))
source(here("code", "read_atlases.R"))

## get correlations among fits ----

stats.subj.cc <- read_subj_stats(suffix = "_intercept-cc_residual", params = "incongruency") %>%
  filter(!subj %in% "562345")
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
    r_ic_cc = cor(beta_ic, beta_cc),
    r_ic_iccc = cor(beta_ic, beta_iccc)
    )

## get correlations with behavior ----


hyps <- combo_paste(c("dlpfc", "lppc", "dmfc"), c("R", "L"), c("target", "incongruency"))
hyps.alt <- combo_paste(c("dlpfc.alt", "lppc", "dmfc.alt"), c("R", "L"), c("target", "incongruency"))

# behav.mod.objs <- readRDS(here("out", "behav", "mod_objs.RDS"))  ## behavioral model objects

blups <- 
  bind_rows(
    read.csv(here("out", "behav", "stroop_blups_rt_group201902.csv"), stringsAsFactors = FALSE)
  ) %>%
  filter(!subj %in% "562345")

## wrange data ----

stats.subj.beh <- full_join(blups, stats.subj, by = "subj")

group.beh <- stats.subj.beh %>%
  
  filter(roi %in% c("dlpfc", "lppc", "dmfc")) %>%
  group_by(roi.hemi) %>%
  summarize(
    r_ic = cor(beta_ic, stroop),
    r_cc = cor(beta_cc, stroop),
    r_iccc = cor(beta_iccc, stroop)
  )

## now for all schemes ---
stats.subj.cc <- read_subj_stats(suffix = "_intercept-cc_residual", roi.set = "mmp") %>%
  filter(!subj %in% "562345")
stats.subj <- read_subj_stats(roi.set = "mmp")

stats.subj <- full_join(
  stats.subj, stats.subj.cc,
  by = c("subj", "is.analysis.group", "roi", "param", "num.roi", "hemi", "community", "community.short", "network", "md"),
  suffix = c("_ic", "_iccc")
  )

group.cor <- stats.subj %>%
  # filter(roi %in% c("dlpfc", "lppc", "dmfc")) %>%
  # group_by(roi.hemi) %>%
  group_by(roi) %>%
  summarize(r = cor(beta_ic, beta_iccc))

range(group.cor$r)

R <- diag(16)
colnames(R) <- bias.items
rownames(R) <- bias.items
R[bias.items.con, bias.items.con] <- 1
sum(R[lower.tri(R)])