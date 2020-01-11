## about ----
## 
## fits representational models with kendall's tau-a.
## 
## mike freund, 06 jan 2020a
## takes about an hour to run in current form (2020-01-06)

## TODO
## about
## loop for statistic (euclidean)


## setup ----

library(here)
library(mikeutils)
library(dplyr)
library(purrr)
library(magrittr)
library(data.table)
library(doParallel)  ## parallelize tau-a calculation
library(foreach)
source(here("code", "strings.R"))
source(here("code", "read_atlases.R"))
source(here("code", "read_masks.R"))

## read models

rsv.models.ltri <- fread(here("out", "rsa", "mods", "rsv_bias_lower-triangles.csv"), data.table = FALSE)

## format

variables <- as.matrix(rsv.models.ltri[sapply(rsv.models.ltri, is.numeric)])
variables <- cbind(variables, conclust = variables[, "incongruency"] | variables[, "congruency"])

## orthogonalize the categorical models

project <- function(y, X) c(X %*% MASS::ginv(crossprod(X)) %*% t(X) %*% y) ## y onto column-space of X
X <- scale(variables[, c("target", "distractor", "incongruency")])
P <- cbind(
  project(X[, 1], X[, c(2, 3)]),
  project(X[, 2], X[, c(1, 3)]),
  project(X[, 3], X[, c(1, 2)])
)
Z <- X - P
colnames(Z) <- paste0(colnames(Z), ".orth")

variables <- cbind(variables, Z)

# crossprod(Z, X)
# cor(Z, X)
# qcor(vec2mat(Z[, 1] / max(abs(Z[, 1])), bias.items, diag.val = 0))
# qcor(vec2mat(Z[, 2] / max(abs(Z[, 2])), bias.items, diag.val = 0))
# qcor(vec2mat(Z[, 3] / max(abs(Z[, 3])), bias.items, diag.val = 0))
# 
## qr decomposition
#
# Q <- cbind(
#   incongruency.q = qr.Q(qr(X))[, 3],
#   target.q       = qr.Q(qr(X[, 3:1]))[, 3],
#   distractor.q   = qr.Q(qr(X[, c(1, 3, 2)]))[, 3]
# )
# qcor(vec2mat(Q[, "incongruency.q"] / max(abs(Q[, "incongruency.q"])), bias.items, diag.val = 0))
# qcor(vec2mat(Q[, "distractor.q"] / max(abs(Q[, "distractor.q"])), bias.items, diag.val = 0))
# qcor(vec2mat(Q[, "target.q"] / max(abs(Q[, "target.q"])), bias.items, diag.val = 0))


## loop over sets of ROIs ----

sets.of.rois <- c("mmp", "gordon", "masks")

for (set.i in sets.of.rois) {
  # set.i = "mmp"
  
  ## read observed similarity matrices (arrays)
  
  rsarray.rank <- readRDS(
    here(
      "out", "rsa", "obsv",
      paste0("rsarray_pro_bias_acc-only_", set.i, "_pearson_residual-rank.rds")
    )
  )
  
  
  ## prepare similarity matrices for regression ----
  
  ## check if rows and col names are equal (should be, but just to be sure...)
  
  are.rowcol.equal <- isTRUE(
    all.equal(dimnames(rsarray.rank)[[1]], dimnames(rsarray.rank)[[2]])
  )
  if(!are.rowcol.equal) stop("rsarray isn't rsarray!")
  
  ## get indices and values
  
  n.dim <- dim(rsarray.rank)[1]
  is.lower.tri <- lower.tri(diag(n.dim))
  subjs <- dimnames(rsarray.rank)$subj
  rois <- dimnames(rsarray.rank)$roi
  n.mods <- length(subjs) * length(rois)
  
  ## unwrap into lower-triangle vector
  
  rsvectors <- vector("list", n.mods)
  names(rsvectors) <- combo_paste(subjs, rois)
  
  for (subj.i in seq_along(subjs)) {
    for (roi.j in seq_along(rois)) {
      # subj.i = 1; roi.j = 1
      
      rsm.rank <- rsarray.rank[, , subj.i, roi.j]
      
      name.ij <- paste0(subjs[subj.i], "_", rois[roi.j])  ## to match name
      
      rsvectors[[name.ij]] <- rsm.rank[is.lower.tri]
      
    }
  }
  
  ## check numbers
  
  n.rows.rsvectors <- map_dbl(rsvectors, length)
  if(sum(n.rows.rsvectors != 120) > 0) stop("missing row somewhere in rsvectors!")
  
  
  ## fit models ----
  
  n.cores <- detectCores()
  cl <- makeCluster(n.cores - 1)
  registerDoParallel(cl)
  
  taua <- foreach(ii = seq_along(rsvectors), .combine = rbind) %dopar% {
    
    apply(variables, 2, function(x) DescTools::KendallTauA(x, rsvectors[[ii]]))
      
  }
  
  stopCluster(cl)
  
  ## format ----
  
  taua <- as.data.frame(taua)
  taua$id <- names(rsvectors)
  taua <- melt(as.data.table(taua), id.vars = "id", variable.name = "param", value.name = "taua")
  
  ## create subj, roi, and hemi cols from id col
  
  taua <- bind_cols(
    taua,
    reshape2::colsplit(taua$id, pattern = "_", names = c("subj", "roi"))
  )
  
  ## add is.analysis.group col
  
  stroop <- fread(here("data", "behavior-and-events_group201902.csv"))
  sample.analysis <- unique(filter(stroop, is.analysis.group)$subj)
  taua$is.analysis.group <- taua$subj %in% sample.analysis
  
  ## create roi.hemi col, rearrange cols (and drop id col)
  
  taua %<>% select(subj, is.analysis.group, roi, param, taua)
  
  
  ## write ----
  
  fwrite(
    taua,
    here("out", "rsa", "stats", paste0("subjs_pro_bias_acc-only_", set.i, "_pearson_residual_taua.csv"))
  )
  
}

