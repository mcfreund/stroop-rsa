#+ wn-roi-contrasts_setup, include = FALSE
if (interactive()) source(here::here("code", "indiv", "bivar_superparcel.R"))
#+

#+ wn-roi-contrasts_model

## fit (only if RDS object doesn't already exist... skip time-consuming step)

fname.mods.wnroi <- here("out", "indiv", "mods_wnroi.RDS")

if (file.exists(fname.mods.wnroi)) {
  
  mods.wnroi <- readRDS(here("out", "indiv", "mods_wnroi.RDS"))
  
} else {
  
  mods.wnroi <- list(
    
    dmfc_L   = update(fit1.het.trim, rt ~ . + trial.type * dmfc_L_target + trial.type * dmfc_L_incongruency, data = d.dissoc.hlm),
    dlpfc_R  = update(fit1.het.trim, rt ~ . + trial.type * dlpfc_R_target + trial.type * dlpfc_R_incongruency, data = d.dissoc.hlm),
    lppc_R   = update(fit1.het.trim, rt ~ . + trial.type * lppc_R_target + trial.type * lppc_R_incongruency, data = d.dissoc.hlm)
    
  )
  
  saveRDS(mods.wnroi, here("out", "indiv", "mods_wnroi.RDS"))
  
}


lapply(mods.wnroi, summary) %>% lapply(coef)

## get contrasts

W.single <- rbind(c(0, 0, 0, 0, 1, -1))
contrasts.wnroi <- lapply(mods.wnroi, glht, linfct = W.single) %>% lapply(summary, test = adjusted("none"))
tab.wnroi <- contrasts.wnroi %>% map("test") %>% map_df(~ .[c("coefficients", "sigma", "tstat", "pvalues")], .id = "roi")

tab.wnroi %<>% 
  rename(b = coefficients, se = sigma, t = tstat, p = pvalues) %>%
  mutate(
    model = c("DMFC (L)", "DLPFC (R)", "LPPC (R)"),
    contrast = c(
      "$\\text{target}-\\text{incon. } | \\text{ DMFC (L)}$",
      "$\\text{target}-\\text{incon. } | \\text{ DLPFC (R)}$",
      "$\\text{target}-\\text{incon. } | \\text{ LPPC (R)}$"
    )
  )


kable(tab.wnroi, escape = FALSE)


saveRDS(contrasts.wnroi, here("out", "indiv", "contrasts_wnroi.RDS"))
fwrite(tab.wnroi, here("out", "indiv", "tab_wnroi.csv"))


#+ 
