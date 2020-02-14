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
source(here("code", "read_masks.R"))


## loop over sets of ROIs ----

sets.of.rois <- c("mmp", "gordon", "masks")

for (set.i in sets.of.rois) {
  # set.i = "mmp"
  
  ## reprisimil matrix
  
  rsarray <- readRDS(
    here(
      "out", "rsa", "obsv", paste0("rsarray_pro_bias_acc-only_", set.i, "_pearson.rds")
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
  rois  <- dimnames(rsarray)$roi
  rsarray.resid.rank <- array(NA, dim = dim(rsarray), dimnames = dimnames(rsarray))
  rsarray.resid.line <- rsarray.resid.rank
  
  ## loop
  
  for (subj.i in seq_along(subjs)) {
    for (roi.j in seq_along(rois)) {
      # roi.j = 1; subj.i = 1
      
      ## get matrix and unwrap to lower triangle
      
      rsm <- rsarray[, , subj.i, roi.j]
      rsv <- as.matrix(rsm[is.lower.tri])
      
      ## transform and create response matrix
      
      Y <- cbind(rsv, rank(rsv))  ## don't atanh betas, as betas can range |beta > 1| !
      
      ## regress run component from rsv
      
      fits <- .lm.fit(X, Y)
      E <- fits$residuals  ## residuals
      B0 <- fits$coef[1, ]  ## intercepts
      regressed <- sweep(E, 2, B0, "+")  ## re-center residuals (add intercepts)
      
      ## extract, vector to matrix
      
      rsm.line.i <- vec2mat(regressed[, 1], dnames = bias.items)  ## linear regressed
      rsm.rank.i <- vec2mat(regressed[, 2], dnames = bias.items)  ## rank regressed
      
      ## store:
      
      rsarray.resid.line[, , subj.i, roi.j] <- rsm.line.i
      rsarray.resid.rank[, , subj.i, roi.j] <- rsm.rank.i
      
    }
  }
  
  
  ## save ----
  
  saveRDS(
    rsarray.resid.rank, 
    here(
      "out", "rsa", "obsv", 
      paste0("rsarray_pro_bias_acc-only_", set.i, "_pearson_residual-rank.rds")
    )
  )
  
  saveRDS(
    rsarray.resid.line, 
    here(
      "out", "rsa", "obsv", 
      paste0("rsarray_pro_bias_acc-only_", set.i, "_pearson_residual-linear.rds")
    )
  )
  
  
}
