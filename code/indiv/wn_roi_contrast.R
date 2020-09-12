#+ wn-roi-contrasts_model

d.dissoc.hlm %<>% group_by(subj) %>% mutate(rt.s = c(scale(rt)))

mods.lme <- list(
  dlpfc_L_targ = update(fit1.het.trim, rt.s ~ . + trial.type * dlpfc_L_target, data = d.dissoc.hlm),
  dlpfc_R_targ = update(fit1.het.trim, rt.s ~ . + trial.type * dlpfc_R_target, data = d.dissoc.hlm),
  lppc_L_targ  = update(fit1.het.trim, rt.s ~ . + trial.type * lppc_L_target, data = d.dissoc.hlm),
  lppc_R_targ  = update(fit1.het.trim, rt.s ~ . + trial.type * lppc_R_target, data = d.dissoc.hlm),
  mfc_L_incon  = update(fit1.het.trim, rt.s ~ . + trial.type * mfc_L_incongruency, data = d.dissoc.hlm),
  mfc_R_incon  = update(fit1.het.trim, rt.s ~ . + trial.type * mfc_R_incongruency, data = d.dissoc.hlm)
) 
sums.lme <- lapply(mods.lme, summary)
pvals.lme <- lapply(sums.lme, function(.) coef(.)[4, "p-value"])
pvals.lme <- reshape2::melt(bind_rows(pvals.lme), value.name = "p", variable = "id")
pvals.lme$group  <- rep(1:3, each = 2)
pvals.lme %<>% group_by(group) %>% mutate(p.fdr = p.adjust(p, method = "fdr"))

lapply(sums.lme, coef)  ## print results
pvals.lme  ## p values (corrected)


## now for within-region dissociations:

mods$dlpfc_R <- update(mods$dlpfc_R_targ, . ~ . + trial.type * dlpfc_R_incongruency)
mods$lppc_R <- update(mods$lppc_R_targ, . ~ . + trial.type * lppc_R_incongruency)
mods$dmfc_L <- update(mods$dmfc_L_incon, . ~ . + trial.type * dmfc_L_target)

W.single <- rbind("[Btarg(I-C)-Bincon(I-C)]" = c(0, 0, 0, 0, 1, -1))

summary(mods$dlpfc_R)
(contrast.dlpfc_R <- summary(glht(mods$dlpfc_R, W.single), test = adjusted("none")))  ## predicted: negative

summary(mods$lppc_R)
(contrast.lppc_R <- summary(glht(mods$lppc_R, W.single), test = adjusted("none")))  ## predicted: negative

summary(mods$dmfc_L)
(contrast.dmfc_L <- summary(glht(mods$dmfc_L, W.single), test = adjusted("none")))  ## predicted: positive
#+

#+ fullmod
## finally, the 'full' model, to test for 4-way interaction:

mods$dmfc_L.lfp_R <- lme(
  
  rt ~ 
    trial.type * lfp_R_target +
    trial.type * lfp_R_incongruency +
    trial.type * dmfc_L_target +
    trial.type * dmfc_L_incongruency, 
  
  random  = ~ trial.type | subj,
  data    = d.dissoc.hlm,
  weights = varIdent(form = ~ 1 | subj),
  control = lmeControl(maxIter = 1e5, msMaxIter = 1e5, niterEM = 1e5, msMaxEval = 1e5),
  method  = "REML"
  
)

summary(mods$dmfc_L.lfp_R)

W.full <- rbind(
  "Btarg(I-C)-Bincon(I-C)|lfp"         = c(0, 0, 0, 0, 0, 0, 1, -1, 0, 0),  ## 3-way interaction within lfp
  "Btarg(I-C)-Bincon(I-C)|dmfc"        = c(0, 0, 0, 0, 0, 0, 0, 0, 1, -1),  ## 3-way interaction within dmfc
  "[Btarg(I-C)-Bincon(I-C)](lfp-dmfc)" = c(0, 0, 0, 0, 0, 0, 1, -1, -1, 1)  ## 4-way interaction
)

(contrast.dmfc_L.lfp_R <- summary(glht(mods$dmfc_L.lfp_R, W.full), test = adjusted("none")))

## suppression?


## plot
# 
# ## refit model to get level-2 residuals sans mfc_L_target:
# 
# mods.lme$sans.mfc.targ <- update(mods.lme$mfc_L.lfp_R, . ~ . - trial.type * mfc_L_target)
# 
# resids <- data.frame(
#   
#   stroop = w.super %>% filter(is.analysis.group) %>% .$stroop,
#   
#   stroop.resid = ranef(mods.lme$sans.mfc.targ)[, 2],
#   
#   mfc_L_targ.resid = lm(
#     mfc_L_target ~ 
#       lfp_R_target + lfp_R_incongruency + mfc_L_incongruency, 
#     w.super %>% filter(is.analysis.group)
#     )$residuals,
#   
#   mfc_L_targ = w.super %>% filter(is.analysis.group) %>% .$mfc_L_target
#   
# )
# 
# cor(resids)
# 
# resids %>%
#   
#   ggplot(., aes(mfc_L_targ.resid, stroop.resid)) + 
#   stat_boot_ci(n = 1E4, alpha = 0.3, fill = colors.model["target"]) + 
#   geom_point(color = "white", shape = 21, size = 4, fill = colors.model["target"]) +
#   
#   labs(title = "partial correlation: stroop~MFC_L_target")

## get OLS estimates:

fits.bivar <- list(
  
  dmfc_L_target       = lm(stroop ~ dmfc_L_target, w.super %>% filter(is.analysis.group)),
  dmfc_L_incongruency = lm(stroop ~ dmfc_L_incongruency, w.super %>% filter(is.analysis.group)),
  lfp_R_target       = lm(stroop ~ lfp_R_target, w.super %>% filter(is.analysis.group)),
  lfp_R_incongruency = lm(stroop ~ lfp_R_incongruency, w.super %>% filter(is.analysis.group))
  
)
r2.bivar <- purrr::map(fits.bivar, summary) %>% purrr::map_dbl("r.squared")


fit <- lm(
  stroop ~ dmfc_L_target + lfp_R_target + lfp_R_incongruency + dmfc_L_incongruency,
  w.super %>% filter(is.analysis.group) %>% mutate_if(is.numeric, scale)
)
summary(fit)
car::vif(fit)  ## not very collinear
kappa(fit)  ## relatively low...

r2.plus <-
  
  fits.bivar[-grep("dmfc_L_target", names(fits.bivar))] %>%
  
  purrr::map(function(x) update(x, . ~ . + dmfc_L_target)) %>%
  purrr::map(summary) %>%
  purrr::map_dbl("r.squared")

r2.delta <- r2.plus - r2.bivar[names(r2.plus)]

r2.delta.pev <- r2.delta / r2.bivar["dmfc_L_target"]  ## ~20-fold increase in variance explained....


## save results ----

indiv.mod.objs <- list(
  mods = mods, 
  sep.p = pvals, 
  contrasts.3way = list(
    dmfc_L = contrast.dmfc_L, 
    lppc_R = contrast.lppc_R,
    dlpfc_R = contrast.dlpfc_R
  ),
  contrasts.4way = contrast.dmfc_L.lfp_R
)
saveRDS(indiv.mod.objs, here("out", "indiv", "indiv.RDS"))

# indiv.mod.objs.sep <- readRDS(here("out", "indiv", "indiv_separate.RDS"))
# mods <- indiv.mod.objs.sep$mods





#+ 
