#' ## read data

#+ prelim-behav-vset_setup, include = FALSE

if (interactive()) {
  
  source(here::here("code", "packages.R"))
  source(here("code", "strings.R"))
  source(here("code", "funs.R"))
  source(here("code", "plots.R"))
  source(here("code", "read_atlases.R"))
  
  behav <- fread(here("in", "behavior-and-events_group201902.csv"))
  
}

#+


#' ## model and write

#+ model

## initial fit

behav.vset.rt <- behav %>% filter(acc == 1, !is.na(rt), rt < 3000, rt > 250, !is.analysis.set)

fit.vset <- lme(
  rt ~ trial.type, 
  random  = ~ trial.type | subj,
  data    = behav.vset.rt,
  weights = varIdent(form = ~ 1 | subj),
  control = lmeControl(maxIter = 1e5, msMaxIter = 1e5, niterEM = 1e5, msMaxEval = 1e5),
  method  = "REML"
)

## trim and re-fit

behav.vset.rt$resid.p <- resid(fit.vset, type = "p")
behav.vset.rt$is.far.out <- farout(behav.vset.rt$resid.p)
fit.vset.trim <- update(fit.vset, subset = !is.far.out)

## extract predictions

blups.vset <- as.data.frame(coef(fit.vset.trim))
blups.vset %<>% rename(congr = "(Intercept)", stroop = "trial.typeincon") %>% tibble::rownames_to_column("subj")

## write and save

write.csv(blups.vset, here("out", "behav", "stroop_blups_rt_group201902_validation.csv"),  row.names = FALSE)
saveRDS(fit.vset.trim, here("out", "behav", "fit1-het-trim_group201902_validation_set.RDS"))

#+