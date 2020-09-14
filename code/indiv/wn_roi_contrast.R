#+ wn-roi-contrasts_model

d.dissoc.hlm %<>% group_by(subj) %>% mutate(rt.s = c(scale(rt)))

mods <- list(
  dlpfc_L_targ = update(fit1.het.trim, rt.s ~ . + trial.type * dlpfc_L_target, data = d.dissoc.hlm),
  dlpfc_R_targ = update(fit1.het.trim, rt.s ~ . + trial.type * dlpfc_R_target, data = d.dissoc.hlm),
  lppc_L_targ  = update(fit1.het.trim, rt.s ~ . + trial.type * lppc_L_target, data = d.dissoc.hlm),
  lppc_R_targ  = update(fit1.het.trim, rt.s ~ . + trial.type * lppc_R_target, data = d.dissoc.hlm),
  dmfc_L_incon  = update(fit1.het.trim, rt.s ~ . + trial.type * dmfc_L_incongruency, data = d.dissoc.hlm),
  dmfc_R_incon  = update(fit1.het.trim, rt.s ~ . + trial.type * dmfc_R_incongruency, data = d.dissoc.hlm)
) 
sums <- lapply(mods, summary)
pvals <- lapply(sums, function(.) coef(.)[4, "p-value"])
pvals <- reshape2::melt(bind_rows(pvals), value.name = "p", variable = "id")
pvals$group  <- rep(1:3, each = 2)
pvals %<>% group_by(group) %>% mutate(p.fdr = p.adjust(p, method = "fdr"))

lapply(sums, coef)  ## print results
pvals  ## p values (corrected)


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
