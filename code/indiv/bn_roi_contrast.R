#+ bn-roi-contrasts_setup, include = FALSE
if (interactive()) source(here::here("code", "indiv", "bivar_superparcel.R"))
#+

#+ bn-roi-contrasts_model

## fit (only if RDS object doesn't already exist... skip time-consuming step)

fname.mods.bnroi <- here("out", "indiv", "mods_bnroi.RDS")

if (file.exists(fname.mods.bnroi)) {
  
  mods.bnroi <- readRDS(here("out", "indiv", "mods_bnroi.RDS"))
  
} else {
  
  mods.bnroi <- list(
    
    incon_dmfc_dlpfc  = update(fit1.het.trim, rt ~ . + trial.type * dmfc_L_incongruency + trial.type * dlpfc_R_incongruency, data = d.dissoc.hlm),
    target_dmfc_dlpfc = update(fit1.het.trim, rt ~ . + trial.type * dmfc_L_target + trial.type * dlpfc_R_target, data = d.dissoc.hlm),
    
    incon_dmfc_lppc   = update(fit1.het.trim, rt ~ . + trial.type * dmfc_L_incongruency + trial.type * lppc_R_incongruency, data = d.dissoc.hlm),
    target_dmfc_lppc  = update(fit1.het.trim, rt ~ . + trial.type * dmfc_L_target + trial.type * lppc_R_target, data = d.dissoc.hlm)
    
  )
  
  saveRDS(mods.bnroi, here("out", "indiv", "mods_bnroi.RDS"))
  
}


lapply(mods.bnroi, summary) %>% lapply(coef)

## get contrasts

contrasts.bnroi <- lapply(mods.bnroi, glht, linfct = W.single) %>% lapply(summary, test = adjusted("none"))
tab.bnroi <- contrasts.bnroi %>% map("test") %>% map_df(~ .[c("coefficients", "sigma", "tstat", "pvalues")], .id = "roi")

tab.bnroi %<>% 
  rename(b = coefficients, se = sigma, t = tstat, p = pvalues) %>%
  mutate(
    model = c("incon.", "target", "incon.", "target"),
    contrast = c(
      "$\\text{DMFC (L)} - \\text{DLPFC (R) } | \\text{ incon.}$",
      "$\\text{DMFC (L)} - \\text{DLPFC (R) } | \\text{ target}$",
      
      "$\\text{DMFC (L)} - \\text{LPPC (R) } | \\text{ incon.}$",
      "$\\text{DMFC (L)} - \\text{LPPC (R) } | \\text{ target}$"
    )
  )


kable(tab.bnroi, escape = FALSE)

saveRDS(contrasts.bnroi, here("out", "indiv", "contrasts_bnroi.RDS"))
fwrite(tab.bnroi, here("out", "indiv", "tab_bnroi.csv"))


#+ 
