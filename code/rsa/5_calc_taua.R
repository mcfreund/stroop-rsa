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
source(here("code", "_strings.R"))
source(here("code", "_read_atlases.R"))

## read models

rsv.models.ltri <- fread(here("out", "rsa", "mods", "rsv_bias_lower-triangles.csv"), data.table = FALSE)

## format

variables <- as.matrix(rsv.models.ltri[sapply(rsv.models.ltri, is.numeric)])
variables <- cbind(variables, conclust = variables[, "incongruency"] | variables[, "congruency"])

## loop over atlases ----

for (atlas.i in seq_along(atlas.key)) {
  # atlas.i = 1
  
  name.atlas.i <- names(atlas.key)[atlas.i]
  
  ## read observed similarity matrices (arrays)
  
  rsarray.rank <- readRDS(
    here(
      "out", "rsa", "obsv",
      paste0("rsarray_pro_bias_acc-only_", name.atlas.i, "_pearson_residual-rank.rds")
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
  hemis <- dimnames(rsarray.rank)$hemi
  n.mods <- length(subjs) * length(rois) * length(hemis)
  
  ## unwrap into lower-triangle vector
  
  rsvectors <- vector("list", n.mods)
  names(rsvectors) <- combo_paste3(subjs, rois, hemis)
  
  for (subj.i in seq_along(subjs)) {
    for (roi.j in seq_along(rois)) {
      for (hemi.k in seq_along(hemis)) {
        # subj.i = 1; roi.j = 1; hemi.k = 1
        
        rsm.rank <- rsarray.rank[, , subj.i, roi.j, hemi.k]
        
        name.ijk <- paste0(subjs[subj.i], "_", rois[roi.j], "_", hemis[hemi.k])  ## to match name
        
        rsvectors[[name.ijk]] <- rsm.rank[is.lower.tri]
        
      }
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
    reshape2::colsplit(taua$id, pattern = "_", names = c("subj", "roi", "hemi"))
  )
  
  ## add is.analysis.group col
  
  stroop <- fread(here("data", "behavior-and-events_group201902.csv"))
  sample.analysis <- unique(filter(stroop, is.analysis.group)$subj)
  taua$is.analysis.group <- taua$subj %in% sample.analysis
  
  ## create roi.hemi col, rearrange cols (and drop id col)
  
  taua %<>%
    mutate(roi.hemi = paste0(roi, "_", hemi)) %>%
    select(subj, is.analysis.group, roi.hemi, roi, hemi, param, taua)
  
  
  ## write ----
  
  fwrite(
    taua,
    here("out", "rsa", "stats", paste0("subjs_pro_bias_acc-only_", name.atlas.i, "_pearson_residual_taua.csv"))
  )
  
}

