## about ----
## fits general linear models on each subject's RSMs from each parcel, computes simple summary statistics, 
## then writes the results to files.
## 
## updated, 2019-03-05: added / commented input / ouput lines to read an array contaning 'run-regressed' RSMs, 
## and write a file of 'run-regressed' stats.
## see "run_effects" knitr in this directory for clarification.
## edited, june 04, 2019: add z / t column in results df, for effect size
## 
## mike freund, 2019-02-24
## adapted for new project directory 2019-12-29


## set up env. ----

## dependencies

library(here)
library(dplyr)
library(purrr)
library(reshape2)
library(data.table)
source(here("code", "_strings.R"))
source(here("code", "_get_atlas.R"))

# source(here("..", "gen", "funs", "_get_dirs_local.R"))
# source(here("..", "gen", "funs", "_funs.R"))
# source(here("r", "group-201902", "_get_misc_vars.R"))

## strings and funs

file.suffix <- ifelse(do.clusters, "_multiparcel", "")
if (do.held.out) {
  file.subj.sample.suffix <- "_held-out"
  file.rds.suffix <- "_held-out"
} else {
  file.subj.sample.suffix <- "_for-analysis" 
  file.rds.suffix <- ""
}

## list of subjects for analysis
group201902.analysis <- read.csv(
  file.path(dir.box.stroop.sheets, paste0("subj-sample_group-201902", file.subj.sample.suffix, ".csv")),
  stringsAsFactors = FALSE
  )[, 1]
## commented 2019-03-05
# atlases <- c("mmp", "gordon")
atlases <- "mmp"
# atlases <- "gordon"
model.names <- c("target", "distractor", "congruency")
# model.names <- c("target", "distractor", "congruency", "incongruency")
model.method <- "pcor"
# stats <- "pearson"
# group <- "group-201902"
# method <- "pro_bias_acc-only"


## loop over atlases ----

for (n.atlas.i in seq_along(atlases)) {
  # n.atlas.i = 1
  
  ## commented 2019-03-05
  # fname.i <- paste0("rsarray_pearson_", atlases[n.atlas.i], "_group-201902_pro_bias_acc-only.rds")
  # rsarray <- readRDS(file.path(dir.box.stroop.data , fname.i))
  fname.i <- paste0(
    "rsarray_pearson_", atlases[n.atlas.i], "_group-201902_pro_bias_acc-only_rank-residual", file.suffix, ".rds"
    )
  rsarray <- readRDS(file.path(dir.box.stroop.data , fname.i))
  rsarray <- rsarray[, , group201902.analysis, , ]  ## get analysis sample.
  
  ## format similarity matrices ----
  
  ## check if rows and col names are equal (should be, but just to be sure...) --
  are.rowcol.equal <- isTRUE(all.equal(dimnames(rsarray)[[1]], dimnames(rsarray)[[2]]))
  if(!are.rowcol.equal) stop("rsarray isn't rsarray!")
  n.dim <- dim(rsarray)[1]
  is.upper.tri <- upper.tri(diag(n.dim), diag = TRUE)
  n.tri <- sum(!is.upper.tri) ## num unique elements ((ndim^2 - ndim) / 2)
  
  ## unwrap into lower-triangle vector and create model regressors --
  subjs <- dimnames(rsarray)$subj
  rois <- dimnames(rsarray)$roi
  hemis <- dimnames(rsarray)$hemi
  rsvectors <- vector("list", length(subjs) * length(rois) * length(hemis))
  names(rsvectors) <- combo_paste3(subjs, rois, hemis)
  ## NB: prefer loops here to vectorized code bc of lower memory demand.
  ## will take a moment, though.
  for (n.subj.i in seq_along(subjs)) {
    for (n.roi.j in seq_along(rois)) {
      for (n.hemi.k in seq_along(hemis)) {
        # n.subj.i = 1; n.roi.j = 1; n.hemi.k = 1
        slice.ijk <- rsarray[, , n.subj.i, n.roi.j, n.hemi.k]
        slice.ijk[is.upper.tri] <- NA
        ## unwrap to vector, make dimnames cols, and remove NAs (i.e., upper.tri and diag):
        vector.ijk <- data.table::melt(slice.ijk, na.rm = TRUE) 
        name.ijk <- paste0(subjs[n.subj.i], "_", rois[n.roi.j], "_", hemis[n.hemi.k])  ## to match name
        ## create model regressors (target, distractor, and congruency):
        row.info <- split.str.item(vector.ijk$r)  ## see _get_misc_vars.R
        col.info <- split.str.item(vector.ijk$c)
        rsvectors[[name.ijk]] <- cbind(
          vector.ijk,
          target      = as.numeric(row.info$color == col.info$color),
          distractor  = as.numeric(row.info$word == col.info$word),
          congruency  = as.numeric(row.info$congruency == col.info$congruency)
          # incongruency = as.numeric(row.info$congruency == col.info$congruency & col.info$congruency == "I")
          )
      }
    }
  }
  
  ## check numbers --
  n.rows.rsvectors <- map_dbl(rsvectors, nrow)
  if(sum(n.rows.rsvectors != 120) > 0) stop("missing row somewhere in rsvectors!")
  
  
  ## fit models and summarize ----
  
  if (model.method == "glm") {
    
    ## general linear and rank models:
    form.rankmod <- scale(rank(value)) ~ scale(target) + scale(distractor) + scale(congruency)
    form.linmod  <- update(form.rankmod, scale(atanh(value)) ~ .)  ## "value" is z-transformed, then normalized

    coef.rankmod <- rsvectors %>% map(lm, formula = form.rankmod) %>% map(coef)
    coef.linmod  <- rsvectors %>% map(lm, formula = form.linmod) %>% map(coef)

    ## get output:
    coef.rankmod <- as.data.frame(do.call(rbind, coef.rankmod))  ## numeric (w/in list) to df
    coef.linmod  <- as.data.frame(do.call(rbind, coef.linmod))

    coef.rankmod <- coef.rankmod %>%  ## remove intercept col and rename slope cols
      select(target = "scale(target)", distractor = "scale(distractor)", congruency = "scale(congruency)")
    coef.linmod  <- coef.linmod %>%
      select(target = "scale(target)", distractor = "scale(distractor)", congruency = "scale(congruency)")

    coef.rankmod <- coef.rankmod %>%  ## split rownames into multiple ID cols, then melt to long-form
      bind_cols(colsplit(rownames(.), pattern = "_", names = c("subj", "parcel", "hemi"))) %>%
      melt(id.vars = c("subj", "parcel", "hemi"), value.name = "rank.beta")
    coef.linmod <- coef.linmod %>%
      bind_cols(colsplit(rownames(.), pattern = "_", names = c("subj", "parcel", "hemi"))) %>%
      melt(id.vars = c("subj", "parcel", "hemi"), value.name = "linear.beta")

    ## combine in single DF:
    coefs <- full_join(coef.rankmod, coef.linmod, by = c("subj", "parcel", "hemi", "variable"))

    ## summarize:
    results <- coefs %>%
      group_by(parcel, hemi, variable) %>%
      summarize(
        mean.rank.beta       = mean(rank.beta),
        mean.linear.beta     = mean(linear.beta),
        rank.beta.srtest.p   = wilcox.test(rank.beta, alternative = "greater")$p.value,
        linear.beta.srtest.p = wilcox.test(linear.beta, alternative = "greater")$p.value
      ) %>%
      rename(rank.beta = mean.rank.beta, linear.beta = mean.linear.beta)  ## for consistency
    
  } else if (model.method == "pcor") {
    
    ## partial correlations:
    rsvectors <- rsvectors %>% map(~.[c("value", model.names)])  ## get relevant cols (regressors and correlation)
    ## commented 2019-03-05
    # rsvectors <- rsvectors %>% map(~ mutate(., value = atanh(value)))  ## z-transform correlation (for linear model)
    r   <- rsvectors %>% map(psych::partial.r, method = "pearson") %>% map(~.["value", model.names])  ## lin mod
    rho <- rsvectors %>% map(psych::partial.r, method = "spearman") %>% map(~.["value", model.names])  ## rank mod
    
    ## get output ---
    r    <- as.data.frame(do.call(rbind, r))  ## numeric (w/in list) to df
    rho  <- as.data.frame(do.call(rbind, rho))
    
    r <- r %>%  ## split rownames into multiple ID cols, then melt to long-form
      bind_cols(colsplit(rownames(.), pattern = "_", names = c("subj", "parcel", "hemi"))) %>%
      melt(id.vars = c("subj", "parcel", "hemi"), value.name = "r")
    rho <- rho %>%
      bind_cols(colsplit(rownames(.), pattern = "_", names = c("subj", "parcel", "hemi"))) %>%
      melt(id.vars = c("subj", "parcel", "hemi"), value.name = "rho")
    
    ## combine in single DF:
    coefs <- full_join(r, rho, by = c("subj", "parcel", "hemi", "variable"))
    
    ## summarize:
    results <- coefs %>%
      group_by(parcel, hemi, variable) %>%
      summarize(
        mean.rho     = mean.rzr(rho),  ## tanh(mean(atanh(x)))
        mean.r       = mean.rzr(r),
        rho.srtest.v = wilcox.test(rho, alternative = "greater")$statistic,
        r.srtest.v   = wilcox.test(r, alternative = "greater")$statistic,
        rho.srtest.p = wilcox.test(rho, alternative = "greater")$p.value,
        r.srtest.p   = wilcox.test(r, alternative = "greater")$p.value
      ) %>%
      rename(rho = mean.rho, r = mean.r)  ## for consistency
    
  }
  
    
  ## write ----
  
  ## commented 2019-03-05
  # fname.suffix <- paste0(model.method, "_rsm-pearson_", atlases[n.atlas.i], "_group-201902_pro_bias_acc-only.csv")
  fname.suffix <- paste0(
    model.method, "_rsm-pearson_", atlases[n.atlas.i], "_group-201902_pro_bias_acc-only_residual",
    file.rds.suffix, file.suffix, ".csv"
    )
  fwrite(coefs, file.path(dir.box.stroop.stats, paste0("subject-coefs_", fname.suffix)))
  fwrite(results, file.path(dir.box.stroop.stats, paste0("parcel-coefs_", fname.suffix)))
  
}  ## end atlas loop
