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

sets.of.rois <- c("mmp", "gordon", "masks")
is.lower.tri <- lower.tri(diag(length(bias.items)))
run.rsm <- as.matrix(fread(here("out", "rsa", "mods", "rsm_pro_bias_run.csv")), rownames = 1)
run.rsv <- mat2vec(run.rsm, value.name = "run.model")  ## unwrap run model to lower-tri vector
U <- cbind(1, run.rsv$run.model)  ## model: intercept, run regressor


## loop over sets of ROIs ----


for (set.i in sets.of.rois) {
  # set.i = "mmp"
  
  ## reprisimil matrix
  
  rsarray <- readRDS(here("out", "rsa", "obsv", paste0("rsarray_pro_bias_acc-only_", set.i, "_pearson.rds")))
  
  
  ## regress ----
  
  ## initialize indices and lists
  
  subjs <- dimnames(rsarray)$subj
  rois  <- dimnames(rsarray)$roi
  rsarray.resid.rank <- array(NA, dim = dim(rsarray), dimnames = dimnames(rsarray))
  rsarray.resid.line <- rsarray.resid.rank
  
  for (subj.i in seq_along(subjs)) {
    for (roi.j in seq_along(rois)) {
      # roi.j = 1; subj.i = 1
      
      ## get matrix and unwrap to lower triangle
      
      rsm <- rsarray[, , subj.i, roi.j]
      rsv <- as.matrix(rsm[is.lower.tri])
      R <- cbind(atanh(rsv), rank(rsv))  ## one for linear, one for rank
      
      fits <- .lm.fit(U, R)  ## regress run component from rsv
      B1 <- coef(fits)[2, ]  ## slopes
      Y <- R - U[, 2] %*% t(B1)    ## unscaled RSA response vectors
      
      rsm.line.i <- vec2mat(tanh(Y[, 1]), dnames = bias.items)  ## extract, vector to matrix
      rsm.rank.i <- vec2mat(Y[, 2], dnames = bias.items, diag.val = sum(is.lower.tri) + 1)  ## diag to n. unique elms+1
      
      rsarray.resid.rank[, , subj.i, roi.j] <- rsm.rank.i  ## store
      rsarray.resid.line[, , subj.i, roi.j] <- rsm.line.i  ## store
      
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
      paste0("rsarray_pro_bias_acc-only_", set.i, "_pearson_residual-line.rds")
    )
  )
  

}
