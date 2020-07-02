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

session <- "pro"


## loop over sets of ROIs ----

if (session == "bas") {
  glmname <- "bas_bias_acc-only_downsamp"
  sets.of.rois <- "mmp"
} else if (session == "pro") {
  glmname <- "pro_bias_acc-only"
  # sets.of.rois <- c("mmp", "gordon", "masks")
  sets.of.rois <- "masks"
} else if (session == "pro_donwsamp") {
  glmname <- "pro_bias_acc-only_downsamp"
  sets.of.rois <- "mmp"
}




for (set.i in sets.of.rois) {
  # set.i = "mmp"
  
  ## reprisimil matrix
  
  rsarray <- readRDS(here("out", "rsa", "obsv", paste0("rsarray_", glmname, "_", set.i, "_pearson.rds")))
  
  if (session == "pro_downsamp") {
    run.rsm <- as.matrix(read.csv(here("out", "rsa", "mods", paste0("rsm_bas_bias_run.csv")), row.names = 1))
  } else {
    run.rsm <- as.matrix(read.csv(here("out", "rsa", "mods", paste0("rsm_", session, "_bias_run.csv")), row.names = 1))
  }
  
  ## regress ----
  
  ## unwrap run model to lower-tri vector:
  
  run.rsv <- mat2vec(run.rsm, value.name = "run.model")
  X <- cbind(1, run.rsv$run.model)  ## model: intercept, run regressor
  
  ## if baseline, remove any subjs with any NA
  
  has.all.data <- !apply(rsarray, "subj", function(x) any(is.na(x)))
  rsarray <- rsarray[, , has.all.data, ]
  
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
      paste0("rsarray_", glmname, "_", set.i, "_pearson_residual-rank.rds")
    )
  )
  
  saveRDS(
    rsarray.resid.line, 
    here(
      "out", "rsa", "obsv", 
      paste0("rsarray_", glmname, "_", set.i, "_pearson_residual-linear.rds")
    )
  )
  
  
}
