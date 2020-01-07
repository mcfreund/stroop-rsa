## about ----

## set up env. ----

library(here)
library(magrittr)
library(dplyr)
library(data.table)
source(here("..", "gen", "funs", "_get_dirs_local.R"))
source(here("..", "gen", "funs", "_funs.R"))
source(here("r", "group-201902", "_get_misc_vars.R"))
# source(here("r", "group-201902", "_read_files_local.R"))

## list of subjects for analysis

group201902.analysis <- fread(file.path(dir.box.stroop.sheets, "subj-sample_group-201902_for-analysis.csv"))$x
rsarray <- readRDS(
  file.path(dir.box.stroop.data , "rsarray_pearson_mmp_group-201902_pro_bias_acc-only_rank-residual.rds")
)
rsarray <- rsarray[, , group201902.analysis, , ]  ## get analysis sample.

fname.suffix <- "pcor_rsm-pearson_mmp_group-201902_pro_bias_acc-only_residual_congruency-schemes.csv"
coefs <- fread(file.path(dir.box.stroop.stats, paste0("subject-coefs_", fname.suffix)))
results <- fread(file.path(dir.box.stroop.stats, paste0("parcel-coefs_", fname.suffix)))

## check two:
results %>% filter(parcel == "9-46d", hemi == "l")
results %>% filter(parcel == "IP0", hemi == "r")

## check all
results %>% filter(rho.srtest.p < 0.001)

## paired comparison:
coefs %<>%
  select(subj, parcel, hemi, variable, rho) %>%
  tidyr::spread(variable, rho) %>%
  mutate(delta = tanh(atanh(incongruency) - atanh(congruency)))  ## POSITIVE MEANS INCON BETTER FIT

paired <- coefs %>%
  group_by(parcel, hemi) %>%
  summarize(
    mean.congruency     = mean.rzr(congruency),  ## tanh(mean(atanh(x)))
    congruency.srtest.v = wilcox.test(congruency, alternative = "greater")$statistic,
    congruency.srtest.p = wilcox.test(congruency, alternative = "greater")$p.value,
    
    mean.incongruen     = mean.rzr(incongruency),  ## tanh(mean(atanh(x)))
    incongruen.srtest.v = wilcox.test(incongruency, alternative = "greater")$statistic,
    incongruen.srtest.p = wilcox.test(incongruency, alternative = "greater")$p.value,
    
    mean.delta          = mean.rzr(delta),  ## tanh(mean(atanh(x)))
    delta.srtest.v      = wilcox.test(delta, alternative = "two.sided")$statistic,
    delta.srtest.p      = wilcox.test(delta, alternative = "two.sided")$p.value,
    
  ) %>%
  rename(congruency = mean.congruency, incongruen = mean.incongruen, delta = mean.delta) ## for consistency

## reshape
paired %<>%
  select(parcel, hemi, congruency, incongruen, delta) %>%
  melt(value.name = "rho") %>%
  full_join(
    paired %>%
      select(parcel, hemi, congruency.srtest.p, incongruen.srtest.p, delta.srtest.p) %>%
      melt(value.name = "p") %>%
      mutate(variable = gsub(".srtest.p", "", variable)),
    by = c("parcel", "hemi", "variable")
    ) %>%
  full_join(
    paired %>%
      select(parcel, hemi, congruency.srtest.v, incongruen.srtest.v, delta.srtest.v) %>%
      melt(value.name = "v") %>%
      mutate(variable = gsub(".srtest.v", "", variable))
  )

## correct p value
paired %<>%
  dplyr::full_join(atlas.key$mmp %>% rename(parcel = "roi"), by = "parcel")  %>%  ## need ROI numbers
  dplyr::mutate(parcel.hemi = paste0(parcel, "_", hemi)) %>%
  ## create p.adjusted col:
  dplyr::group_by(variable) %>%
  dplyr::mutate(
    p.adj.wb = p.adjust(p, method = "fdr"),
    is.roi   = parcel.hemi %in% rois.mmp$roi.hemi[rois.mmp$dmcc2.taskgen > 1]  ## dmcc2.taskgen == 2 -> most conserv.
  ) %>%
  dplyr::group_by(variable, is.roi) %>%
  dplyr::mutate(
    p.adj.roi = p.adjust(p, method = "fdr"),
    p.adj.roi = ifelse(is.roi, p.adj.roi, p.adj.wb)
  )


## examine results ----

## where incongruency is positive and sig:
parcels.i <- paired %>% filter(variable == "incongruen", p.adj.roi < 0.05) %>% .$parcel.hemi
## and where delta is sig & positive:
paired %>% filter(parcel.hemi %in% parcels.i, variable == "delta", rho > 0, p < 0.05)

## where congruency is positive and sig:
parcels.c <- paired %>% filter(variable == "congruency", p.adj.roi < 0.05) %>% .$parcel.hemi
## and where delta is sig & negative:
paired %>% filter(parcel.hemi %in% parcels.c, variable == "delta", rho < 0, p < 0.05)
## or where delta is neg:
paired %>% filter(parcel.hemi %in% parcels.c, variable == "delta", rho < 0)

library(ggplot2)
coefs %>%
  filter(parcel == "9-46d", hemi == "l") %>%
  ggplot(aes(congruency, incongruency)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

## where incongruency is positive and sig, and congruency is not sig
paired %>% filter(parcel.hemi %in% parcels.i, variable == "congruency", p > 0.05)

## where congruency is positive and sig, and incongruency is not sig
paired %>% filter(parcel.hemi %in% parcels.c, variable == "incongruency", p > 0.05)

## of the parcels sig coding for congruency, those that are better fit by incongruent model
incons <- paired %>% 
  filter(parcel.hemi %in% parcels.c, variable == "delta") %>%
  ungroup %>%
  mutate(p.adj.cong = p.adjust(p, "fdr")) %>%
  filter(p.adj.cong < 0.05) %>%
  .$parcel.hemi

paired %>% 
  filter(parcel.hemi %in% parcels.c, variable == "delta") %>%
  ungroup %>%
  mutate(p.adj.cong = p.adjust(p, "fdr")) %>%
  write.csv(
    file.path(
      dir.box.stroop.stats, 
      "incong-vs-cong_incon_pcor_rsm-pearson_mmp_group-201902_pro_bias_acc-only_residual.csv"
    )
  )


## pairwise comparisons with other coding schemes (distractor, target)

coefs.all <- data.table::fread(
  file.path(dir.box.stroop.stats, "subject-coefs_pcor_rsm-pearson_mmp_group-201902_pro_bias_acc-only_residual.csv")
) %>%
  dplyr::full_join(atlas.key$mmp %>% dplyr::rename(parcel = "roi"), by = "parcel")  %>%  ## need ROI numbers
  dplyr::mutate(parcel.hemi = paste0(parcel, "_", hemi))

coefs.all
coefs.incon <- coefs %>%
  select(congruency = incongruency, subj, parcel, hemi) %>%
  melt(value.name = "rho")
  
paircomps <- coefs.all %>%
  select(subj, parcel, hemi, variable, rho) %>%
  filter(variable %in% c("distractor", "target")) %>%
  bind_rows(coefs.incon) %>%
  tidyr::spread(variable, rho) %>%
  melt %>%
  group_by(parcel, hemi) %>%
  summarize(
    pairwise = list(
      pairwise.wilcox.test(value, variable, paired = TRUE, alternative = "two.sided", p.adjust.method = "fdr")
    )
  )  ## warnings indicate non-exact p-values
pvals <- sapply(paircomps$pairwise, function(.) .$p.value[lower.tri(diag(2), diag = TRUE)])
pvals <- t(pvals)
# paircomps$pairwise[[1]]$p.value
# paircomps$pairwise[[1]]$p.value[lower.tri(diag(2), diag = TRUE)]
colnames(pvals) <- c("cong.dist", "targ.cong", "targ.dist")
paircomps <- cbind(as.data.frame(paircomps[c("parcel", "hemi")]), as.data.frame(pvals))
paircomps$parcel.hemi <- paste0(paircomps$parcel, "_", paircomps$hemi)

write.csv(
  paircomps,
  file.path(
    dir.box.stroop.stats, 
    "parcel-coefs_pairwise-incon_pcor_rsm-pearson_mmp_group-201902_pro_bias_acc-only_residual.csv"
    )
  )

## interaction ----
paired %>% group_by(variable) %>% filter(rho == max(rho))

inter <- coefs %>%
  filter(parcel %in% c("IP1", "9-46d"), hemi == "l") %>%
  select(subj, parcel, congruency, incongruency) %>%
  melt %>%
  mutate(variable = as.factor(variable), parcel = as.factor(parcel))
fit <- lmer(atanh(value) ~ variable * parcel + (1 | subj), inter)
summary(fit)
plot(fit)

## possibly influenced by larger number of incongruenct versus congruent items?

## correlations with behavior and new congruency scheme
## first, assess correlation between coding readouts
coefs %>%
  mutate(parcel.hemi = paste0(parcel, "_", hemi)) %>%
  # filter(parcel.hemi %in% incons) %>%
  group_by(parcel.hemi) %>%
  summarize(
    r = cor(congruency, incongruency)
  )  ## not that different.
coefs %>%
  mutate(parcel.hemi = paste0(parcel, "_", hemi)) %>%
  filter(parcel.hemi %in% incons) %>%
  ggplot(aes(congruency, incongruency)) +
  geom_point(aes(color = parcel.hemi))

stroop.pro.rt <- stroop.pro %>% filter(acc == 1, !is.na(rt))
stroop.pro.rt <- stroop.pro.rt %>%
  group_by(subj) %>%
  mutate(
    is.sd3.5 = abs(rt - mean(rt)) > 3.5 * sd(rt)
  )
library(lme4)
mrt.base <- lmer(
  rt ~ trial.type + (trial.type | subj), 
  stroop.pro.rt, subset = !is.sd3.5, control = lmer.cl
)
library(nlme)
cl1 <- lmeControl(
  maxIter = 100000, msMaxIter = 100000, niterEM = 100000,
  msMaxEval = 100000, tolerance = 0.000001, msTol = 0.0000001, returnObject = TRUE,
  minAbsParApVar = 0.05, opt = c("nlminb"), optimMethod = "BFGS"
)
mrt.base.varident <- lme(
  rt ~ trial.type, random = ~ 1 + trial.type | subj, 
  correlation = varIdent(form = ~ 1 | subj),
  data = stroop.pro.rt %>% filter(!is.sd3.5), control = cl1
)
behav <- data.frame(
  rt.hv = coef(mrt.base.varident)$trial.typei,
  rt = coef(mrt.base)$subj$trial.typei,
  subj = rownames(coef(mrt.base.varident))
)

incons <- coefs %>%
  mutate(parcel.hemi = paste0(parcel, "_", hemi)) %>%
  filter(parcel.hemi %in% incons) %>%
  full_join(behav, by = "subj") %>%
  group_by(parcel.hemi)
incons %>% summarize(rc = cor(rt, congruency), ri = cor(rt, incongruency))
incons %>% summarize(rc = cor(rt, congruency), ri = cor(rt, incongruency))

coefs %>%
  mutate(parcel.hemi = paste0(parcel, "_", hemi)) %>%
  full_join(behav, by = "subj") %>%
  group_by(parcel.hemi) %>%
  summarize(
    rc = cor(rt, congruency), 
    ri = cor(rt, incongruency),
    d  = tanh(atanh(rc) - atanh(ri))
  ) %>% View