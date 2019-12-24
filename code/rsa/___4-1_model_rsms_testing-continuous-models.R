## about ----
## mike freund, 2019-03-19

## set up env. ----

do.clusters <- TRUE
do.held.out <- TRUE

## dependencies

## NB: BOX DRIVE

library(here)
library(dplyr)
library(purrr)
library(reshape2)
library(data.table)
source(here("..", "gen", "funs", "_get_dirs_local.R"))
source(here("..", "gen", "funs", "_funs.R"))
source(here("r", "group-201902", "_get_misc_vars.R"))
coding.schemes <- read.csv(file.path(dir.box.stroop, "_to-clean", "_coding-schemes", "coding-schemes.csv"))
coding.schemes <- coding.schemes %>% 
  filter(set == "bias") %>% 
  select(r = a, c = b, cielab = cielab.vec, image = image.vec, silhouette = silh.vec) %>%
  mutate(r = as.character(r), c = as.character(c))
# coding.schemes %>% select(contains("vec")) %>% cor

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
model.names <- c("target", "distractor", "congruency", "image", "cielab", "silhouette")
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
          # congruency = as.numeric(row.info$congruency == col.info$congruency & col.info$congruency == "I")
        )
        rsvectors[[name.ijk]] <- rsvectors[[name.ijk]] %>% 
          mutate(r = as.character(r), c = as.character(c)) %>% 
          left_join(coding.schemes, by = c("r", "c"))
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
    

    r   <- rsvectors %>% map(cor, method = "pearson") %>% map(~.["value", model.names])  ## lin mod
    rho <- rsvectors %>% map(cor, method = "spearman") %>% map(~.["value", model.names])  ## rank mod

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
        rho.srtest.p = wilcox.test(rho, alternative = "greater")$p.value,
        r.srtest.p   = wilcox.test(r, alternative = "greater")$p.value
      ) %>%
      rename(rho = mean.rho, r = mean.r)  ## for consistency
    
  }
  
  
  ## write ----
  
  ## commented 2019-03-05
  # fname.suffix <- paste0(model.method, "_rsm-pearson_", atlases[n.atlas.i], "_group-201902_pro_bias_acc-only.csv")
  fname.suffix <- paste0(
    model.method, "_rsm-pearson_", atlases[n.atlas.i], "_group-201902_pro_bias_acc-only_residual_continuous_models",
    file.rds.suffix, file.suffix, ".csv"
  )
  fwrite(coefs, file.path(dir.box.stroop.stats, paste0("subject-coefs_", fname.suffix)))
  fwrite(results, file.path(dir.box.stroop.stats, paste0("parcel-coefs_", fname.suffix)))
  
}  ## end atlas loop

if (do.clusters) {
  
  ## test dlpfc ----
  # c("p9-46v_l", "8C_l", "8Av_l", "i6-8_l")
  ## group tendencies
  results.dlpfc <- results %>% filter(parcel %in% "dlpfc")
  results.dlpfc %>%
    group_by(hemi) %>%
    mutate(p.adj = p.adjust(rho.srtest.p, method = "fdr"))
  
  
  coefs.dlpfc <- coefs %>% filter(parcel %in% "dlpfc")
  coefs.dlpfc %>%
    group_by(hemi) %>%
    summarize(ps = list(pairwise.wilcox.test(rho, variable, paired = TRUE, p.adjust.method = "fdr"))) %>% .$ps
  
  ## behavioral correlation:
  
  stroop.pro <- fread(file.path(dir.box.crosstask.sheets, "dmcc2_behavior-and-events-stroop.csv")) %>%
    filter(subj %in% group201902.analysis, session == "pro") %>%
    arrange(subj, session, run, trial.num) %>%
    group_by(subj, session, run) %>%
    mutate(
      er            = 1 - acc,
      post.er       = lag(er),
      n1.trial.type = lag(trial.type),
      n2.trial.type = lag(n1.trial.type)
    )
  stroop.pro.rt <- stroop.pro %>% 
    filter(acc == 1, !is.na(rt), rt < 3000, rt > 250)  %>%
    group_by(subj) %>%
    mutate(is.sd3.5 = abs(rt - mean(rt)) > 3.5 * sd(rt))
  mrt.base <- lme4::lmer(
    rt ~ trial.type + (trial.type | subj) + (trial.type | color), 
    stroop.pro.rt, subset = !is.sd3.5, control = lmer.cl
  )
  behav <- coef(mrt.base)$subj %>%
    tibble::rownames_to_column("subj") %>%
    select(subj, rt = trial.typei)
  
  coefs.dlpfc <- coefs.dlpfc %>% full_join(behav, by = "subj")
  
  coefs.dlpfc.cor <- coefs.dlpfc %>%
    select(-r) %>%
    tidyr::spread(variable, rho) %>%
    select(-subj)
  
  
  coefs.dlpfc.cor %>%
    filter(hemi == "l") %>%
    select(rt:silhouette) %>%
    cor
  
  coefs.dlpfc.cor %>%
    filter(hemi == "r") %>%
    select(rt:silhouette) %>%
    cor(method = "pearson")
  
  
  plot(coefs.dlpfc.cor$silhouette[coefs.dlpfc.cor$hemi == "l"], coefs.dlpfc.cor$rt[coefs.dlpfc.cor$hemi == "l"])
  plot(coefs.dlpfc.cor$distractor[coefs.dlpfc.cor$hemi == "l"], coefs.dlpfc.cor$rt[coefs.dlpfc.cor$hemi == "l"])
  
  plot(coefs.dlpfc.cor$target[coefs.dlpfc.cor$hemi == "r"], coefs.dlpfc.cor$rt[coefs.dlpfc.cor$hemi == "r"])
  plot(coefs.dlpfc.cor$image[coefs.dlpfc.cor$hemi == "r"], coefs.dlpfc.cor$rt[coefs.dlpfc.cor$hemi == "r"])

  
  
  
  
  
  rsvectors[grepl("V8_l", names(rsvectors))] %>%
    map(~ select(., target, value)) %>% 
    map(~ DescTools::KendallTauA(x = .$target, y = .$value)) -> target
  rsvectors[grepl("V8_l", names(rsvectors))] %>%
    map(~ select(., cielab, value)) %>% 
    map(~ DescTools::KendallTauA(x = .$cielab, y = .$value)) -> cielab
  wilcox.test(unlist(cielab), unlist(target), paired = TRUE)
  plot(unlist(cielab), unlist(target))
  abline(0, 1)
  # cor(unlist(cielab), unlist(target))
  names(cielab) <- gsub("(^.*)_.*_.$", "\\1", names(cielab))
  cielab <- data.frame(cielab = unlist(cielab), subj = names(unlist(cielab)))
  cielab <- behav %>% full_join(cielab, by = "subj")
  plot(cielab$rt, cielab$cielab)
  cor(cielab$rt, cielab$cielab)
  
  
  # rsvectors[grepl("V8_l", names(rsvectors))] %>%
  #   map(~ select(., silhouette, value)) %>% 
  #   map(~ DescTools::KendallTauA(x = .$silhouette, y = .$value)) -> silhouette
  names(silhouette) <- gsub("(^.*)_.*_.$", "\\1", names(silhouette))
  silhouette <- data.frame(unlist(silhouette), subj = names(unlist(silhouette)))
  silhouette <- behav %>% full_join(silhouette, by = "subj")
  plot(silhouette$rt, silhouette$unlist.silhouette.)
  
  
} else {
  
  ## test V8 ----
  
  results.V8 <- results %>% mutate(parcel.hemi = paste0(parcel, "_", hemi)) %>% filter(parcel %in% c("V8"))
  results.V8 %>%
    group_by(parcel.hemi) %>%
    mutate(p.adj = p.adjust(rho.srtest.p, method = "fdr"))
  
  coefs.V8 <- coefs  %>% mutate(parcel.hemi = paste0(parcel, "_", hemi)) %>%
    filter(parcel %in% c("V8"))
  coefs.V8 %>%
    group_by(parcel.hemi) %>%
    summarize(ps = list(pairwise.wilcox.test(rho, variable, paired = TRUE, p.adjust.method = "fdr"))) %>% .$ps
  
  ## bhavioral correlations
  stroop.pro <- fread(file.path(dir.box.crosstask.sheets, "dmcc2_behavior-and-events-stroop.csv")) %>%
    filter(subj %in% group201902.analysis, session == "pro") %>%
    arrange(subj, session, run, trial.num) %>%
    group_by(subj, session, run) %>%
    mutate(
      er            = 1 - acc,
      post.er       = lag(er),
      n1.trial.type = lag(trial.type),
      n2.trial.type = lag(n1.trial.type)
    )
  stroop.pro.rt <- stroop.pro %>% 
    filter(acc == 1, !is.na(rt), rt < 3000, rt > 250)  %>%
    group_by(subj) %>%
    mutate(is.sd3.5 = abs(rt - mean(rt)) > 3.5 * sd(rt))
  mrt.base <- lme4::lmer(
    rt ~ trial.type + (trial.type | subj) + (trial.type | color), 
    stroop.pro.rt, subset = !is.sd3.5, control = lmer.cl
  )
  behav <- coef(mrt.base)$subj %>%
    tibble::rownames_to_column("subj") %>%
    select(subj, rt = trial.typei)
  
  coefs.V8 <- coefs.V8 %>% full_join(behav, by = "subj")
  
  coefs.V8 %>%
    filter(hemi == "l", variable == "congruency") %>%
    select(rho, rt) %>%
    cor
  
  
  
  ## test dlpfc ----
  # c("p9-46v_l", "8C_l", "8Av_l", "i6-8_l")
  ## group tendencies
  results.dlpfc <- results %>% mutate(parcel.hemi = paste0(parcel, "_", hemi)) %>%
    filter(parcel.hemi %in% c("p9-46v_l", "8C_l", "8Av_l", "i6-8_l", "p9-46v_r", "8C_r", "8Av_r", "i6-8_r"))
  results.dlpfc %>%
    group_by(parcel.hemi) %>%
    mutate(p.adj = p.adjust(rho.srtest.p, method = "fdr")) %>% View
  
  coefs.dlpfc <- coefs  %>% mutate(parcel.hemi = paste0(parcel, "_", hemi)) %>%
    filter(parcel.hemi %in% c("p9-46v_l", "8C_l", "8Av_l", "i6-8_l", "p9-46v_r", "8C_r", "8Av_r", "i6-8_r"))
  coefs.dlpfc %>%
    group_by(parcel.hemi) %>%
    summarize(ps = list(pairwise.wilcox.test(rho, variable, paired = TRUE, p.adjust.method = "fdr"))) %>% .$ps
  
  ## bhavioral correlations
  stroop.pro <- fread(file.path(dir.box.crosstask.sheets, "dmcc2_behavior-and-events-stroop.csv")) %>%
    filter(subj %in% group201902.analysis, session == "pro") %>%
    arrange(subj, session, run, trial.num) %>%
    group_by(subj, session, run) %>%
    mutate(
      er            = 1 - acc,
      post.er       = lag(er),
      n1.trial.type = lag(trial.type),
      n2.trial.type = lag(n1.trial.type)
    )
  stroop.pro.rt <- stroop.pro %>% 
    filter(acc == 1, !is.na(rt), rt < 3000, rt > 250)  %>%
    group_by(subj) %>%
    mutate(is.sd3.5 = abs(rt - mean(rt)) > 3.5 * sd(rt))
  mrt.base <- lme4::lmer(
    rt ~ trial.type + (trial.type | subj) + (trial.type | color), 
    stroop.pro.rt, subset = !is.sd3.5, control = lmer.cl
  )
  behav <- coef(mrt.base)$subj %>%
    tibble::rownames_to_column("subj") %>%
    select(subj, rt = trial.typei)
  
  coefs.dlpfc <- coefs.dlpfc %>% full_join(behav, by = "subj")
  
  coefs.dlpfc.cor <- coefs.dlpfc %>%
    select(-r) %>%
    tidyr::spread(variable, rho) %>%
    select(-subj)
  
  coefs.dlpfc.cor %>%
    filter(parcel.hemi == "p9-46v_l") %>%
    select(rt:silhouette) %>%
    cor
  
  coefs.dlpfc.cor %>%
    filter(parcel.hemi == "p9-46v_r") %>%
    select(rt:silhouette) %>%
    cor
  
  coefs.dlpfc.cor %>%
    filter(parcel.hemi == "i6-8_l") %>%
    select(rt:silhouette) %>%
    cor
  
  coefs.dlpfc.cor %>%
    filter(parcel.hemi == "i6-8_r") %>%
    select(rt:silhouette) %>%
    cor
  
  coefs.dlpfc.cor %>%
    filter(parcel.hemi == "8Av_l") %>%
    select(rt:silhouette) %>%
    cor
  
  coefs.dlpfc.cor %>%
    filter(hemi == "r") %>%
    select(rt:silhouette) %>%
    cor(method = "pearson")
  
  
  plot(coefs.dlpfc.cor$silhouette[coefs.dlpfc.cor$hemi == "l"], coefs.dlpfc.cor$rt[coefs.dlpfc.cor$hemi == "l"])
  plot(coefs.dlpfc.cor$distractor[coefs.dlpfc.cor$hemi == "l"], coefs.dlpfc.cor$rt[coefs.dlpfc.cor$hemi == "l"])
  
  plot(coefs.dlpfc.cor$target[coefs.dlpfc.cor$hemi == "r"], coefs.dlpfc.cor$rt[coefs.dlpfc.cor$hemi == "r"])
  plot(coefs.dlpfc.cor$image[coefs.dlpfc.cor$hemi == "r"], coefs.dlpfc.cor$rt[coefs.dlpfc.cor$hemi == "r"])
  
  
  ## taua
  rho <- rsvectors %>% map(cor, method = "spearman") %>% map(~.["value", model.names])  ## rank mod
  
  rsvectors[grepl("p9-46v_l", names(rsvectors))] %>%
    map(~ select(., target, value)) %>% 
    map(~ DescTools::KendallTauA(x = .$target, y = .$value)) -> target
  rsvectors[grepl("p9-46v_l", names(rsvectors))] %>%
    map(~ select(., cielab, value)) %>% 
    map(~ DescTools::KendallTauA(x = .$cielab, y = .$value)) -> cielab
  plot(unlist(cielab), unlist(target))
  abline(0, 1)
  rsvectors[grepl("p9-46v_l", names(rsvectors))] %>%
    map(~ select(., silhouette, value)) %>% 
    map(~ DescTools::KendallTauA(x = .$silhouette, y = .$value)) -> silhouette
  names(silhouette) <- gsub("(^.*)_.*_.$", "\\1", names(silhouette))
  silhouette <- data.frame(unlist(silhouette), subj = names(unlist(silhouette)))
  silhouette <- behav %>% full_join(silhouette, by = "subj")
  plot(silhouette$rt, silhouette$unlist.silhouette.)
  
  

  }

