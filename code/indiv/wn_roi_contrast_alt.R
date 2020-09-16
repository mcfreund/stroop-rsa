#+ wn-roi-contrasts-alt_model

d.dissoc.hlm.alt <- full_join(
  fit1.het.trim$data %>% filter(!is.far.out), 
  w.super %>% filter(is.analysis.group) %>% ungroup %>% dplyr::select(subj, one_of(hyps.alt)),
  by = "subj"
)

names(d.dissoc.hlm.alt) <- gsub("\\.alt", "", names(d.dissoc.hlm.alt))
d.dissoc.hlm.alt %<>% group_by(subj) %>% mutate(rt.s = c(scale(rt)))

mods.alt <- list(
  dlpfc_L_targ  = update(fit1.het.trim, rt.s ~ . + trial.type * dlpfc_L_target, data = d.dissoc.hlm.alt),
  dlpfc_R_targ  = update(fit1.het.trim, rt.s ~ . + trial.type * dlpfc_R_target, data = d.dissoc.hlm.alt),
  lppc_L_targ   = update(fit1.het.trim, rt.s ~ . + trial.type * lppc_L_target, data = d.dissoc.hlm.alt),
  lppc_R_targ   = update(fit1.het.trim, rt.s ~ . + trial.type * lppc_R_target, data = d.dissoc.hlm.alt),
  dmfc_L_incon  = update(fit1.het.trim, rt.s ~ . + trial.type * dmfc_L_incongruency, data = d.dissoc.hlm.alt),
  dmfc_R_incon  = update(fit1.het.trim, rt.s ~ . + trial.type * dmfc_R_incongruency, data = d.dissoc.hlm.alt)
) 
sums.alt <- lapply(mods.alt, summary)
pvals.alt <- lapply(sums.alt, function(.) coef(.)[4, "p-value"])
pvals.alt <- reshape2::melt(bind_rows(pvals.alt), value.name = "p", variable = "id")
pvals.alt$group  <- rep(1:3, each = 2)
pvals.alt %<>% group_by(group) %>% mutate(p.fdr = p.adjust(p, method = "fdr"))

lapply(sums.alt, coef)  ## print results
pvals.alt  ## p values (corrected)


## now for within-region dissociations:

mods.alt$dlpfc_R <- update(mods.alt$dlpfc_R_targ, . ~ . + trial.type * dlpfc_R_incongruency)
mods.alt$lppc_R <- update(mods.alt$lppc_R_targ, . ~ . + trial.type * lppc_R_incongruency)
mods.alt$dmfc_L <- update(mods.alt$dmfc_L_incon, . ~ . + trial.type * dmfc_L_target)

summary(mods.alt$dlpfc_R)
(contrast.dlpfc_R.alt <- summary(glht(mods.alt$dlpfc_R, W.single), test = adjusted("none")))  ## predicted: negative

summary(mods.alt$lppc_R)
(contrast.lppc_R.alt <- summary(glht(mods.alt$lppc_R, W.single), test = adjusted("none")))  ## predicted: negative

summary(mods$dmfc_L.alt)
(contrast.dmfc_L.alt <- summary(glht(mods.alt$dmfc_L, W.single), test = adjusted("none")))  ## predicted: positive


## save

saveRDS(mods.alt, here("out", "indiv", "mods_wnroi_alt.RDS"))
saveRDS(pvals.alt, here("out", "indiv", "pvals_wnroi_alt.RDS"))
saveRDS(contrast.dlpfc_R.alt, here("out", "indiv", "contrast.dlpfc_R_alt.RDS"))
saveRDS(contrast.lppc_R.alt, here("out", "indiv", "contrast.lppc_R_alt.RDS"))
saveRDS(contrast.dmfc_L.alt, here("out", "indiv", "contrast.dmfc_L_alt.RDS"))

#+ 
