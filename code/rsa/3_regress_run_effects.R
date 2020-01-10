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
source(here("code", "strings.R"))
source(here("code", "read_atlases.R"))


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
  
  run.rsv <- mat2vec(run.rsm, value.name = "run.model")
  X <- cbind(1, run.rsv$run.model)  ## model: intercept, run regressor
  
  ## initialize indices and lists
  
  is.lower.tri <- lower.tri(diag(length(bias.items)))
  subjs <- dimnames(rsarray)$subj
  parcels <- dimnames(rsarray)$roi
  hemis <- c("l", "r")
  rsarray.resid.rank <- array(NA, dim = dim(rsarray), dimnames = dimnames(rsarray))
  rsarray.resid.line <- rsarray.resid.rank
  
  ## loop
  
  for (n.subj.i in seq_along(subjs)) {
    for (n.parcel.j in seq_along(parcels)) {
      for (n.hemi.k in seq_along(hemis)) {
        # n.hemi.k = 1; n.parcel.j = 1; n.subj.i = 1
        
        ## get matrix and unwrap to lower triangle
        
        rsm <- rsarray[, , n.subj.i, n.parcel.j, n.hemi.k]
        rsv <- as.matrix(rsm[is.lower.tri])
        
        ## transform and create response matrix
        
        Y <- cbind(atanh(rsv), rank(rsv))
        
        ## regress run component from rsv
        
        fits <- .lm.fit(X, Y)
        E <- fits$residuals  ## residuals
        B0 <- fits$coef[1, ]  ## intercepts
        regressed <- sweep(E, 2, B0, "+")  ## re-center residuals (add intercepts)
        
        ## extract
        
        regressed.r <- tanh(regressed[, 1])
        regressed.rank <- regressed[, 2]
        
        ## vector to matrix:
        
        rsm.line.i <- vec2mat(regressed.r, dnames = bias.items)  ## linear regressed
        rsm.rank.i <- vec2mat(regressed.rank, dnames = bias.items)  ## rank regressed
        
        ## store:
        
        rsarray.resid.line[, , n.subj.i, n.parcel.j, n.hemi.k] <- rsm.line.i
        rsarray.resid.rank[, , n.subj.i, n.parcel.j, n.hemi.k] <- rsm.rank.i
        
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
    rsarray.resid.line, 
    here(
      "out", "rsa", "obsv", 
      paste0("rsarray_pro_bias_acc-only_", name.atlas.i, "_pearson_residual-linear.rds")
    )
  )
  
  
}
