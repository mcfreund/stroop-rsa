## about ----
## 
## description...
## 
## mike freund, 2019-03-17
## adapted for new project directory 2019-12-27


## setup ----

library(here)
library(mikeutils)
library(dplyr)
library(data.table)
library(purrr)
source(here("code", "_strings.R"))
source(here("code", "_get_atlas.R"))


## atlas loop ----

for (atlas.i in seq_along(atlas.key)) {
  # atlas.i = 1
  
  name.atlas.i <- names(atlas.key)[atlas.i]
  
  ## reprisimil matrix
  
  rsarray <- readRDS(
    here(
      "out", "rsa", "obsv", paste0("rsarray_pro_bias_acc-only_", name.atlas.i, "_pearson.rds")
    )
  )
  
  run.rsm <- as.matrix(read.csv(here("out", "rsa", "mods", "rsm_bias_run.csv"), row.names = 1))
  
  
  ## regress ----
  
  ## unwrap run model to lower-tri vector:
  
  run.rsv <- mat2vec(run.rsm)
  names(run.rsv) <- c("r", "c", "run.model")
  
  ## initialize indices and lists
  
  subjs <- dimnames(rsarray)$subj
  parcels <- dimnames(rsarray)$roi
  hemis <- c("l", "r")
  rsarray.resid.rank <- array(NA, dim = dim(rsarray), dimnames = dimnames(rsarray))
  rsarray.resid.linear <- rsarray.resid.rank
  
  ## loop
  
  for (n.subj.i in seq_along(subjs)) {
    for (n.parcel.j in seq_along(parcels)) {
      for (n.hemi.k in seq_along(hemis)) {
        # n.hemi.k = 1; n.parcel.j = 1; n.subj.i = 1
        
        ## get matrix
        
        rsm <- rsarray[, , n.subj.i, n.parcel.j, n.hemi.k]
        
        ## matrix to vector:
        
        rsv <- mat2vec(rsm)
        rsv <- full_join(rsv, run.rsv, by = c("r", "c"))
        
        ## model:
        
        fit.linear <- lm(atanh(value) ~ run.model, rsv)
        fit.rank   <- update(fit.linear, rank(value) ~ .)
        resids <- cbind(
          rsv[c("r", "c")], 
          ## re-center (and transform back to r):
          residuals.linear = tanh(residuals(fit.linear) + coef(fit.linear)["(Intercept)"]),
          residuals.rank   = residuals(fit.rank) + coef(fit.rank)["(Intercept)"]
        )
        
        ## vector to matrix:
        
        rsm.linear.i <- vec2mat(resids[, "residuals.linear"], dnames = bias.items)
        rsm.rank.i   <- vec2mat(resids[, "residuals.rank"], dnames = bias.items)
        
        ## store:
        
        name.ijk <- paste(subjs[n.subj.i], parcels[n.parcel.j], hemis[n.hemi.k], sep = "_")
        
        rsarray.resid.linear[, , n.subj.i, n.parcel.j, n.hemi.k] <- rsm.linear.i
        rsarray.resid.rank[, , n.subj.i, n.parcel.j, n.hemi.k]   <- rsm.rank.i
        
      }
    }
  }
  
  
  ## save ----
  
  saveRDS(
    rsarray.resid.rank, 
    here(
      "out", "rsa", "obsv", 
      paste0("rsarray_pro_bias_acc-only_", name.atlas.i, "_pearson_residual-rank.rds")
    )
  )
  
  saveRDS(
    rsarray.resid.linear, 
    here(
      "out", "rsa", "obsv", 
      paste0("rsarray_pro_bias_acc-only_", name.atlas.i, "_pearson_residual-linear.rds")
    )
  )
  
  
}


