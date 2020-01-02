## about ----
## fits general linear models on each subject's RSMs from each parcel then writes the results to files.
## 
## mike freund, created 2019-02-24, updated 2019-03-05
## adapted for new project directory 2019-12-29

## TODO
## add continuous models
## add tau a
## parallelize loops
## check warning on lm.betas (and speed up?)
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

## models

models.wide <- fread(here("out", "rsa", "mods", "bias_full_matrix.csv"))

## formulas
f <- list(
  tdc = ~ target + distractor + congruency,
  tdi = ~ target + distractor + incongruency,
  txi = ~ target + distractor + incongruency + target:incongruency,
  dxi = ~ target + distractor + incongruency + distractor:incongruency
  # all = ~ target + distractor + incongruency + cielab + tphon + silh + orth + dphon
)

## strings and funs


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
  
  
  ## format similarity matrices ----
  
  ## check if rows and col names are equal (should be, but just to be sure...)
  
  are.rowcol.equal <- isTRUE(
    all.equal(dimnames(rsarray.rank)[[1]], dimnames(rsarray.rank)[[2]]) & 
    all.equal(dimnames(rsarray.line)[[1]], dimnames(rsarray.line)[[2]]) &
    all.equal(dimnames(rsarray.line)[[1]], dimnames(rsarray.rank)[[2]])
  )
  if(!are.rowcol.equal) stop("rsarray isn't rsarray!")
  
  ## get values 
  
  n.dim <- dim(rsarray.line)[1]
  is.upper.tri <- upper.tri(diag(n.dim), diag = TRUE)
  n.tri <- sum(!is.upper.tri) ## num unique elements ((ndim^2 - ndim) / 2)
  subjs <- dimnames(rsarray.line)$subj
  rois <- dimnames(rsarray.line)$roi
  hemis <- dimnames(rsarray.line)$hemi
  
  ## unwrap into lower-triangle vector and create model regressors
  
  rsvectors <- vector("list", length(subjs) * length(rois) * length(hemis))
  names(rsvectors) <- combo_paste3(subjs, rois, hemis)
  
  ## NB: prefer loops here to vectorized code bc of lower memory demand.
  ## will take a moment, though.
  for (subj.i in seq_along(subjs)) {
    for (roi.j in seq_along(rois)) {
      for (hemi.k in seq_along(hemis)) {
        # subj.i = 1; roi.j = 1; hemi.k = 1
        
        rsm.line <- rsarray.line[, , subj.i, roi.j, hemi.k]
        rsv.line <- mat2vec(rsm.line)
        rsv.line <- rename(rsv.line, x = r, y = c, r = value)
        
        rsm.rank <- rsarray.rank[, , subj.i, roi.j, hemi.k]
        rsv.rank <- mat2vec(rsm.rank)
        rsv.rank <- rename(rsv.rank, x = r, y = c, rank = value)
        
        rsv <- full_join(rsv.line, rsv.rank, by = c("x", "y"))
        
        ## bind model regressors to observed vector
        
        rsv <- left_join(rsv, models.wide, by = c("x", "y"))  ## left_join (not full) bc models.wide is full matrix
        
        ## store in list (rsvectors)
        
        name.ijk <- paste0(subjs[subj.i], "_", rois[roi.j], "_", hemis[hemi.k])  ## to match name
        rsvectors[[name.ijk]] <- rsv
        
      }
    }
  }
  
  ## check numbers --
  n.rows.rsvectors <- map_dbl(rsvectors, nrow)
  if(sum(n.rows.rsvectors != 120) > 0) stop("missing row somewhere in rsvectors!")
  
  
  ## fit models ----

  stats.subjs <- setNames(vector("list", length(f)), names(f))
  
  for (f.i in seq_along(f)) {
    # f.i = 1

    fits.rank <- rsvectors %>% map(lm, formula = update(f[[f.i]], rank ~ .))  ## fit models
    fits.line <- rsvectors %>% map(lm, formula = update(f[[f.i]], atanh(r) ~ .))
    
    ## get coefficients
    
    coefs.rank <- fits.rank %>% map(~ coef(.)[-1]) %>% do.call(rbind, .)  ## get coefficients
    coefs.line <- fits.line %>% map(~ coef(.)[-1]) %>% do.call(rbind, .)
    
    coefs.rank <- tibble::rownames_to_column(as.data.frame(coefs.rank, stringsAsFactors = FALSE), "id")  ## add rownames
    coefs.line <- tibble::rownames_to_column(as.data.frame(coefs.line, stringsAsFactors = FALSE), "id")
    
    coefs <- bind_rows(
      rank   = melt(as.data.table(coefs.rank), id.vars = "id", variable = "param", value.name = "coef"), 
      linear = melt(as.data.table(coefs.line), id.vars = "id", variable = "param", value.name = "coef"),
      .id = "y"
    )  ## melt to long-form, stack results from rank and linear, and add indicator column for response variable type

    
    ## get betas
    
    betas.rank <- fits.rank %>% map(lm.beta) %>% do.call(rbind, .)
    betas.line <- fits.line %>% map(lm.beta) %>% do.call(rbind, .)
    
    betas.rank <- tibble::rownames_to_column(as.data.frame(betas.rank, stringsAsFactors = FALSE), "id")
    betas.line <- tibble::rownames_to_column(as.data.frame(betas.line, stringsAsFactors = FALSE), "id")
    
    betas <- bind_rows(
      rank   = melt(as.data.table(betas.rank), id.vars = "id", variable = "param", value.name = "beta"), 
      linear = melt(as.data.table(betas.line), id.vars = "id", variable = "param", value.name = "beta"),
      .id = "y"
    )
    
    ## bind and save
    
    stats.subjs[[f.i]] <- full_join(betas, coefs, by = c("y", "id", "param"))
    
  }
  
  
  ## format ----
  
  ## collate into data.frame:
  
  stats.subjs <- bind_rows(stats.subjs, .id = "model")  ## will coerce to character vector and give warning
  
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

