## about ----
## fits general linear models on each subject's RSMs from each parcel then writes the results to files.
## 
## mike freund, created 2019-02-24, updated 2019-03-05
## adapted for new project directory 2019-12-29

## TODO
## about
## add continuous models
## add tau a
## loop for statistic (euclidean)
## loop for atlas


## setup ----

library(here)
library(mikeutils)
library(dplyr)
library(purrr)
library(magrittr)
library(data.table)
source(here("code", "_strings.R"))
source(here("code", "_funs.R"))
source(here("code", "_get_atlas.R"))

## read models

rsv.models.ltri <- fread(here("out", "rsa", "mods", "rsv_bias_lower-triangles.csv"))

## create design matrices

## models will be variants upon these main effects:
tdc <- as.matrix(rsv.models.ltri[, c("target", "distractor", "congruency")])
tdi <- as.matrix(rsv.models.ltri[, c("target", "distractor", "incongruency")])
tdc.std <- scale(tdc)
tdi.std <- scale(tdi)

m <- list(
  tdc = cbind(b0 = 1, tdc),
  tdi = cbind(b0 = 1, tdi),
  txi = cbind(b0 = 1, tdi, ti = rsv.models.ltri$target * rsv.models.ltri$incongruency),
  dxi = cbind(b0 = 1, tdi, di = rsv.models.ltri$distractor * rsv.models.ltri$incongruency)
  # all = cbind(b0 = 1, tdi, )  ## continuous and categorical
  )
m.std <- list(
  tdc = tdc.std,
  tdi = tdi.std,
  txi = cbind(tdi.std, ti = tdi.std[, "incongruency"] * tdi.std[, "target"]),
  dxi = cbind(tdi.std, di = tdi.std[, "incongruency"] * tdi.std[, "distractor"])
  # all = cbind(b0 = 1, tdi, )  ## continuous and categorical
)  ## get interaction btw standardized MEs, not standardized interaction


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
  rsarray.line <- readRDS(
    here(
      "out", "rsa", "obsv",
      paste0("rsarray_pro_bias_acc-only_", name.atlas.i, "_pearson_residual-linear.rds")
    )
  )
  
  
  ## prepare similarity matrices for regression ----
  
  ## check if rows and col names are equal (should be, but just to be sure...)
  
  are.rowcol.equal <- isTRUE(
    all.equal(dimnames(rsarray.rank)[[1]], dimnames(rsarray.rank)[[2]]) & 
    all.equal(dimnames(rsarray.line)[[1]], dimnames(rsarray.line)[[2]]) &
    all.equal(dimnames(rsarray.line)[[1]], dimnames(rsarray.rank)[[2]])
  )
  if(!are.rowcol.equal) stop("rsarray isn't rsarray!")
  
  ## get indices and values
  
  n.dim <- dim(rsarray.line)[1]
  is.lower.tri <- lower.tri(diag(n.dim))
  subjs <- dimnames(rsarray.line)$subj
  rois <- dimnames(rsarray.line)$roi
  hemis <- dimnames(rsarray.line)$hemi
  n.mods <- length(subjs) * length(rois) * length(hemis)
  
  ## unwrap into lower-triangle vector
  
  rsvectors <- vector("list", n.mods)
  names(rsvectors) <- combo_paste3(subjs, rois, hemis)
  
  for (subj.i in seq_along(subjs)) {
    for (roi.j in seq_along(rois)) {
      for (hemi.k in seq_along(hemis)) {
        # subj.i = 1; roi.j = 1; hemi.k = 1
        
        rsm.line <- rsarray.line[, , subj.i, roi.j, hemi.k]  ## get slice
        rsm.rank <- rsarray.rank[, , subj.i, roi.j, hemi.k]

        z <- atanh(rsm.line[is.lower.tri])  ## get lower.triangle vector (and transform to fisher's z)
        rank <- rsm.rank[is.lower.tri]
        
        name.ijk <- paste0(subjs[subj.i], "_", rois[roi.j], "_", hemis[hemi.k])  ## to match name
        rsvectors[[name.ijk]] <- cbind(z, rank)
        
      }
    }
  }
  
  ## check numbers
  
  n.rows.rsvectors <- map_dbl(rsvectors, nrow)
  if(sum(n.rows.rsvectors != 120) > 0) stop("missing row somewhere in rsvectors!")
  
  ## scale
  
  rsvectors.std <- rsvectors %>% map(scale)
  
  
  ## fit models ----
  
  stats.subjs <- setNames(vector("list", length(m)), names(m))
  
  for (m.i in seq_along(m)) {
    # m.i = 1
    
    X <- m[[m.i]]
    X.std <- m.std[[m.i]]
    
    ## get unstandardized coefficients
    
    fits <- rsvectors %>% map(~ .lm.fit( x = X, y = .))  ## fit models (two-column y)
    coefs <- as.data.frame(do.call(rbind, lapply(fits, coef))) ## to single data.frame
    
    names(coefs) <- c("z", "rank")  ## same order as in models
    coefs$param <- rep(colnames(X), n.mods)  ## same orders as in models
    coefs$id <- rep(names(rsvectors), each = ncol(X))
    coefs <- coefs[coefs$param != "b0", ]  ## ditch intercept
    ## put rank and linear (z) coefs in single long-form column:
    coefs <- melt(as.data.table(coefs), id.vars = c("id", "param"), variable = "y", value.name = "coef")
    
    ## get betas
    
    fits.std <- rsvectors.std %>% map(~ .lm.fit( x = X.std, y = .))
    betas <- as.data.frame(do.call(rbind, lapply(fits.std, coef)))
    
    names(betas) <- c("z", "rank")
    betas$param <- rep(colnames(X.std), n.mods)
    betas$id <- rep(names(rsvectors.std), each = ncol(X.std))
    betas <- melt(as.data.table(betas), id.vars = c("id", "param"), variable = "y", value.name = "beta")

    ## bind and save
    
    stats.subjs[[m.i]] <- full_join(coefs, betas, by = c("y", "id", "param"))
    
  }
  
  
  ## format ----
  
  ## collate into data.frame
  
  stats.subjs <- bind_rows(stats.subjs, .id = "model")
  
  ## remove colons from values of param col (and shorten)
  
  stats.subjs$param[grepl("target:incongruency", stats.subjs$param)] <- "ti"
  stats.subjs$param[grepl("distractor:incongruency", stats.subjs$param)] <- "di"
  
  ## create subj, roi, and hemi cols from id col
  
  stats.subjs <- bind_cols(
    stats.subjs,
    reshape2::colsplit(stats.subjs$id, pattern = "_", names = c("subj", "roi", "hemi"))
    )
  
  ## add is.analysis.group col
  
  stroop <- fread(here("data", "behavior-and-events_group201902.csv"))
  sample.analysis <- unique(filter(stroop, is.analysis.group)$subj)
  stats.subjs$is.analysis.group <- stats.subjs$subj %in% sample.analysis
  
  ## create roi.hemi col, rearrange cols (and drop id col)
  
  stats.subjs %<>%
    mutate(roi.hemi = paste0(roi, "_", hemi)) %>%
    select(subj, is.analysis.group, roi.hemi, roi, hemi, model, y, param, coef, beta)
  
  
  ## write ----
  
  fwrite(
    stats.subjs,
    here("out", "rsa", "stats", paste0("subjs_pro_bias_acc-only_", name.atlas.i, "_pearson_residual.csv"))
  )

}

