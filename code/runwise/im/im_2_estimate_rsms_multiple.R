## about ----
## 
## mike freund, 2021-02-17


## setup ----

source(here::here("code", "packages.R"))
source(here::here("code", "strings.R"))
source(here::here("code", "read_atlases.R"))
source(here::here("code", "read_masks.R"))


glm.name <- "pro_bias_acc-only_fmriprep_im_run"

## functions

read_betas <- function(.dir.glms, .regs = bias.items) {
  # .dir.glms = dir.glms
  
  image.betas <- list(NULL, NULL)
  
  for (run.i in 1:2) {
    # run.i = 1
    
    fname.nii <- file.path(paste0(.dir.glms, run.i), paste0("stats_", fit.subjs[subj.i], "_run", run.i, ".nii.gz"))
    if (file.exists(fname.nii)) {
      ## dims of image.run.full [i, j, k, ???, regressor]
      # image.full <- oro.nifti::readNIfTI(fname.nii, reorient = FALSE)
      
      image.full <- RNifti::readNifti(fname.nii)
      
    } else stop("file nonexistant! ", paste0(fname.nii))
    
    
    ## get labels
    
    labs <- mikeutils::afni("3dinfo", paste("-label",  fname.nii))
    labs <- unlist(strsplit(labs, "\\|"))
    is.reg <- grepl(paste0(.regs, "#[0-9]_Coef", collapse = "|"), labs)
    regs <- gsub("([a:Z]*)(#.*)", "\\1", labs[is.reg])
    
    ## put subset of images (only the relevant regressors) into one array
    
    ## image.betas with dims [i, j, k, regs]:
    
    image.betas[[run.i]] <- image.full[, , , 1, is.reg]
    dimnames(image.betas[[run.i]]) <- list(i = NULL, j = NULL, k = NULL, condition = regs)
    
  }
  
  image.betas
  
}



## paths, vars

dir.analysis <- here::here("glms")
files.dir.analysis <- list.files(dir.analysis, pattern = "stats_var", recursive = TRUE)  ## get fit.subjs
files.dir.analysis <- files.dir.analysis[grep(glm.name, files.dir.analysis)]
fit.subjs <- unique(gsub("/results/.*", "", files.dir.analysis))

## regs will be used to pull out (via string match) the statistic from the afni brick;
## thus must have the reg.suffix
n.subj <- length(fit.subjs)
n.bias.items <- length(bias.items)

# cmat <- mikeutils::contrast_matrix(length(bias.items), bias.items)  ## contrast matrix

## loop over sets of ROIs ----
## each iteration collates RSMs into a single array, and saves it as a single .rds file.

for (set.i in sets.of.rois) {
  # set.i = "masks"
  
  ## get numbers and create storage objects
  
  if (set.i == "mmp") {
    n.roi <- nrow(atlas.key[[set.i]])
    roi.names <- atlas.key[[set.i]]$roi
  } else {
    n.roi <- length(masks)
    roi.names <- names(masks)
  }
  
  rsarray <- array(  ## for representational similarity matrices
    NA,
    dim = c(n.bias.items, n.bias.items, n.subj, n.roi),
    dimnames = list(
      .row = bias.items,
      .col = bias.items,
      subj = fit.subjs,
      roi  = roi.names
    )
  )
  
  rsarray_down <- rsarray
  
  crun.glmagg <- setNames(vector("list", length(fit.subjs)), fit.subjs)
  crun.glmagg.wn <- setNames(vector("list", length(fit.subjs)), fit.subjs)
  
  ## loop over subjs and rois ----
  
  for (subj.i in seq_along(fit.subjs)) {
    # subj.i <- 1
    
    ## get afni images
    
    dir.glms <- file.path(dir.analysis, fit.subjs[subj.i], "results", glm.name)
    image.betas <- read_betas(dir.glms)  ## returns list, as arrays may have differing numbers of dimensions by run
    
    crun.glmagg.roi <- setNames(vector("list", n.roi), roi.names)
    crun.glmagg.roi.wn <- setNames(vector("list", n.roi), roi.names)
      
    ## generate rsm for each roi
    
    for (roi.i in seq_len(n.roi)) {  ## careful: roi.i is an index for as.numeric(rois)
      # roi.i <- 23
      
      ## get and apply mask for roi.i
      
      if (set.i == "mmp") {
        mask.i <- atlas[[set.i]] == roi.i 
      } else {
        mask.i <- masks[[roi.i]] == 1
      }
      
      roi.betas <- lapply(image.betas, function(x1) apply(x1, "condition", function(x2) x2[mask.i]))  ## dirty

      ## estimate similarity matrices:
      
      
      ## cross-run ----
      
      R <- cor(roi.betas[[1]], roi.betas[[2]])  ## cross-run: rows = run1, cols = run2
      
      ## aggregate by condition over trials
      
      A1 <- model.matrix(~ 0 + dimnames(R)[[1]])  ## build dummy/indicator matrices
      A2 <- model.matrix(~ 0 + dimnames(R)[[2]])
      
      A1 <- sweep(A1, 2, colSums(A1), "/")  ## scale cols
      A2 <- sweep(A2, 2, colSums(A2), "/")
      
      # R_ave <- tanh(t(A1) %*% atanh(R) %*% A2)  ## average
      
      # rownames(R_ave) <- gsub("dimnames(R)[[1]]", "", rownames(R_ave), fixed = TRUE)  ## remove weird prefix
      # colnames(R_ave) <- gsub("dimnames(R)[[2]]", "", colnames(R_ave), fixed = TRUE)
      
      
      # rsarray[, , subj.i, roi.i] <- R_ave[bias.items, bias.items]
      
      
      
      ## cross-run, agg in glm ----
      
      conds1 <- dimnames(image.betas[[1]])$condition
      conds2 <- dimnames(image.betas[[2]])$condition
      
      distr1 <- gsub("[a-z]", "", conds1)
      distr2 <- gsub("[a-z]", "", conds2)
      targe1 <- gsub("[A-Z]", "", conds1)
      targe2 <- gsub("[A-Z]", "", conds2)
      congr1 <- as.numeric(targe1 == tolower(distr1))
      congr2 <- as.numeric(targe2 == tolower(distr2))
      incon1 <- as.numeric(!congr1)
      incon2 <- as.numeric(!congr2)

      X_distr <- model.matrix(~ 0 + distr1) %*% t(model.matrix(~ 0 + distr2))
      X_targe <- model.matrix(~ 0 + targe1) %*% t(model.matrix(~ 0 + targe2))
      X_congr <- congr1 %*% t(congr2)
      X_incon <- incon1 %*% t(incon2)
      X_diago <- X_distr * X_targe
      
      y <- scale(rank(c(R)))
      X <- scale(cbind(c(X_distr), c(X_targe), c(X_congr), c(X_incon), c(X_diago)))
      fit <- .lm.fit(X, y)
      betas <- coef(fit)
      betas <- as.data.table(betas)
      betas$param <- c("distractor", "target", "congruency", "incongruency", "diagonal")
      
      
      crun.glmagg.roi[[roi.i]] <- betas
      
      
      
      ## cross-run, downsample ----
      
      inds1 <- setNames(vector("list", n.bias.items), bias.items)
      for (cond.i in seq_along(bias.items)) {
        # cond.i = 1
        if (sum(bias.items[cond.i] == conds1) < 4) {
          take.these <- which(bias.items[cond.i] == conds1)
        } else {
          take.these <- sample(which(bias.items[cond.i] == conds1), 3)
        }
        inds1[[cond.i]] <- take.these
      }
      inds1 <- unlist(inds1)
      names(inds1) <- gsub("[0-9]", "", names(inds1))
      
      inds2 <- setNames(vector("list", n.bias.items), bias.items)
      for (cond.i in seq_along(bias.items)) {
        # cond.i = 1
        if (sum(bias.items[cond.i] == conds2) < 4) {
          take.these <- which(bias.items[cond.i] == conds2)
        } else {
          take.these <- sample(which(bias.items[cond.i] == conds2), 3)
        }
        inds2[[cond.i]] <- take.these
      }
      inds2 <- unlist(inds2)
      names(inds2) <- gsub("[0-9]", "", names(inds2))
      
      R_down <- R[inds1, inds2]
      
      ## aggregate by condition over trials
      
      A1_down <- model.matrix(~ 0 + dimnames(R_down)[[1]])  ## build dummy/indicator matrices
      A2_down <- model.matrix(~ 0 + dimnames(R_down)[[2]])
      
      A1_down <- sweep(A1_down, 2, colSums(A1_down), "/")  ## scale cols
      A2_down <- sweep(A2_down, 2, colSums(A2_down), "/")
      
      R_ave_down <- tanh(t(A1_down) %*% atanh(R_down) %*% A2_down)  ## average
      
      rownames(R_ave_down) <- gsub("dimnames(R_down)[[1]]", "", rownames(R_ave_down), fixed = TRUE)  ## remove weird prefix
      colnames(R_ave_down) <- gsub("dimnames(R_down)[[2]]", "", colnames(R_ave_down), fixed = TRUE)
      
      rsarray_down[, , subj.i, roi.i] <- R_ave_down[bias.items, bias.items]
       
      
      
      ## within-run, agg in glm ----
      
      R1 <- cor(roi.betas[[1]], roi.betas[[1]])
      R2 <- cor(roi.betas[[2]], roi.betas[[2]])
      
      X_distr1 <- model.matrix(~ 0 + distr1) %*% t(model.matrix(~ 0 + distr1))
      X_targe1 <- model.matrix(~ 0 + targe1) %*% t(model.matrix(~ 0 + targe1))
      X_congr1 <- congr1 %*% t(congr1)
      X_incon1 <- incon1 %*% t(incon1)
      X_diago1 <- X_distr1 * X_targe1
      is.lt1 <- lower.tri(R1)
      y1 <- scale(rank(R1[is.lt1]))
      X1 <- scale(cbind(X_distr1[is.lt1], X_targe1[is.lt1], X_congr1[is.lt1], X_incon1[is.lt1], X_diago1[is.lt1]))
      fit1 <- .lm.fit(X1, y1)
      betas1 <- coef(fit1)
      betas1 <- as.data.table(betas1)
      betas1$param <- c("distractor", "target", "congruency", "incongruency", "diagonal")

      X_distr2 <- model.matrix(~ 0 + distr2) %*% t(model.matrix(~ 0 + distr2))
      X_targe2 <- model.matrix(~ 0 + targe2) %*% t(model.matrix(~ 0 + targe2))
      X_congr2 <- congr2 %*% t(congr2)
      X_incon2 <- incon2 %*% t(incon2)
      X_diago2 <- X_distr2 * X_targe2
      is.lt2 <- lower.tri(R2)
      y2 <- scale(rank(R2[is.lt2]))
      X2 <- scale(cbind(X_distr2[is.lt2], X_targe2[is.lt2], X_congr2[is.lt2], X_incon2[is.lt2], X_diago2[is.lt2]))
      fit2 <- .lm.fit(X2, y2)
      betas2 <- coef(fit2)
      betas2 <- as.data.table(betas2)
      betas2$param <- c("distractor", "target", "congruency", "incongruency", "diagonal")
      
      
      betas.wn <- full_join(betas1, betas2, by = "param")
      crun.glmagg.roi.wn[[roi.i]] <- betas.wn

      

    }
    
    
    crun.glmagg[[subj.i]] <- bind_rows(crun.glmagg.roi, .id = "roi")
    crun.glmagg.wn[[subj.i]] <- bind_rows(crun.glmagg.roi.wn, .id = "roi")
    
    print(paste0(subj.i, ": subj ", fit.subjs[subj.i], " done!"), quote = FALSE)
    
    
  }  ## end subject loop
  
  crun.glmagg <- bind_rows(crun.glmagg, .id = "subj")
  crun.glmagg.wn <- bind_rows(crun.glmagg.wn, .id = "subj")
  
  
  ## store ----
    
  ## RSA results
  
  # saveRDS(
  #   rsarray, 
  #   here::here(
  #     "out", "rsa", "obsv", 
  #     paste0("rsarray_", glm.name, "_", set.i, ".rds")
  #     )
  #   )
  
  
  saveRDS(
    rsarray_down, 
    here::here(
      "out", "rsa", "obsv", 
      paste0("rsarray_", glm.name, "_downsampled-glmagg_", set.i, ".rds")
    )
  )

  fwrite(
    crun.glmagg,
    here("out", "rsa", "stats",  paste0("subjs_pro_bias_acc-only_fmriprep_im_run_glmagg_", set.i, ".csv"))
  )
  
  fwrite(
    crun.glmagg.wn,
    here("out", "rsa", "stats",  paste0("subjs_pro_bias_acc-only_fmriprep_im_run_glmagg-wn_", set.i, ".csv"))
  )
  
  

}  ## end atlas loop

