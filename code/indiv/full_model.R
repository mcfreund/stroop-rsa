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

## refit model to get level-2 residuals sans mfc_L_target:
# 
# mods$sans.mfc.targ <- update(mods$dmfc_L.lfp_R, . ~ . - trial.type * dmfc_L_target)
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
  lfp_R_target        = lm(stroop ~ lfp_R_target, w.super %>% filter(is.analysis.group)),
  lfp_R_incongruency  = lm(stroop ~ lfp_R_incongruency, w.super %>% filter(is.analysis.group))
  
)
r2.bivar <- map(fits.bivar, summary) %>% purrr::map_dbl("r.squared")


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
