#+ single-roi-contrasts_setup, include = FALSE
# if (interactive()) source(here::here("code", "indiv", "bivar_superparcel.R"))
#+

#+ single-roi-contrasts_model

## fit (only if RDS object doesn't already exist... skip time-consuming step)

# fname.mods.single <- here("out", "indiv", "mods_single.RDS")

if (file.exists(fname.mods.single)) {
  
  mods.single <- readRDS(here("out", "indiv", "mods_single.RDS"))
  
} else {
  
  mods.single <- list(
    
    dlpfc_R  = update(fit1.het.trim, rt ~ . + trial.type * dlpfc_R_target, data = d.dissoc.hlm),
    dlpfc_L  = update(fit1.het.trim, rt ~ . + trial.type * dlpfc_L_target, data = d.dissoc.hlm),
    lppc_R   = update(fit1.het.trim, rt ~ . + trial.type * lppc_R_target, data = d.dissoc.hlm),
    lppc_L   = update(fit1.het.trim, rt ~ . + trial.type * lppc_L_target, data = d.dissoc.hlm),
    dmfc_R   = update(fit1.het.trim, rt ~ . + trial.type * dmfc_R_target, data = d.dissoc.hlm),
    dmfc_L   = update(fit1.het.trim, rt ~ . + trial.type * dmfc_L_target, data = d.dissoc.hlm),
    
    dlpfc_R  = update(fit1.het.trim, rt ~ . + trial.type * dlpfc_R_incongruency, data = d.dissoc.hlm),
    dlpfc_L  = update(fit1.het.trim, rt ~ . + trial.type * dlpfc_L_incongruency, data = d.dissoc.hlm),
    lppc_R   = update(fit1.het.trim, rt ~ . + trial.type * lppc_R_incongruency, data = d.dissoc.hlm),
    lppc_L   = update(fit1.het.trim, rt ~ . + trial.type * lppc_L_incongruency, data = d.dissoc.hlm),
    dmfc_R   = update(fit1.het.trim, rt ~ . + trial.type * dmfc_R_incongruency, data = d.dissoc.hlm),
    dmfc_L   = update(fit1.het.trim, rt ~ . + trial.type * dmfc_L_incongruency, data = d.dissoc.hlm)
  
  )
  
  # saveRDS(mods.single, here("out", "indiv", "mods_single.RDS"))
    
}

lapply(mods.single, summary) %>% lapply(coef)

## wrangle and correct p-values

tab.single <- lapply(mods.single, summary) %>% 
  lapply(function(.) data.frame(coef(.), term = rownames(coef(.)))) %>% 
  bind_rows(.id = "roi") %>% 
  filter(grepl(":", term))

tab.single$group <- rep(1:6, each = 2)  ## group across hemisphere for p-value correction

tab.single %<>% group_by(group) %>% mutate(p.fdr = p.adjust(p.value, method = "fdr"))

tab.single %<>% select(term, roi, b = Value, se = Std.Error, t = t.value, p = p.value, p.fdr)


## save

# fwrite(tab.single, here("out", "indiv", "table_single.csv"))

#+ 
