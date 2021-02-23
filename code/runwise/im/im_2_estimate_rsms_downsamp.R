## about ----
## 
## mike freund, 2021-02-17


## setup ----

source(here::here("code", "packages.R"))
source(here::here("code", "strings.R"))
source(here::here("code", "read_atlases.R"))
source(here::here("code", "read_masks.R"))

n.cores <- detectCores()
glm.name <- "pro_bias_acc-only_fmriprep_im_run"
type <- c("crun", "wnr1", "wnr2", "both")

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
    dim = c(n.bias.items, n.bias.items, n.subj, n.roi, length(type)),
    dimnames = list(
      .row = bias.items,
      .col = bias.items,
      subj = fit.subjs,
      roi  = roi.names,
      type = type
    )
  )
  

  ## loop over subjs and rois ----
  
  image.betas <- setNames(vector("list", length(fit.subjs)), fit.subjs)
  
  for (subj.i in seq_along(fit.subjs)) {
    # subj.i <- 1
    
    ## get afni images
    
    dir.glms <- file.path(dir.analysis, fit.subjs[subj.i], "results", glm.name)
    image.betas[[subj.i]] <- read_betas(dir.glms)  ## returns list, as arrays may have differing numbers of dimensions by run
    
    ## get downsampled trial inds
    
    conds1 <- dimnames(image.betas[[subj.i]][[1]])$condition
    conds2 <- dimnames(image.betas[[subj.i]][[2]])$condition
    
    nsamp <- 1000
    set.seed(0)
    
    inds1 <- setNames(vector("list", n.bias.items), bias.items)
    for (cond.i in seq_along(bias.items)) {
      # cond.i = 2
      if (sum(bias.items[cond.i] == conds1) > 3) {  ## downsample to 3 trials if >3 trials, otherwise keep # trials as is
        n.trials1 <- 3
      } else {
        n.trials1 <- sum(bias.items[cond.i] == conds1)
      }
      if (n.trials1 == 1) {  ## for flexible behavior of sample
        inds1[[cond.i]] <- matrix(rep(which(bias.items[cond.i] == conds1), nsamp), nrow = n.trials1)
      } else {
        inds1[[cond.i]] <- replicate(nsamp, sample(which(bias.items[cond.i] == conds1), n.trials1), nsamp)
      }
      rownames(inds1[[cond.i]]) <- rep(bias.items[[cond.i]], n.trials1)
    }
    inds1 <- do.call(rbind, inds1)
    
    inds2 <- setNames(vector("list", n.bias.items), bias.items)
    for (cond.i in seq_along(bias.items)) {
      # cond.i = 2
      if (sum(bias.items[cond.i] == conds2) > 3) {  ## downsample to 3 trials if >3 trials, otherwise keep # trials as is
        n.trials2 <- 3
      } else {
        n.trials2 <- sum(bias.items[cond.i] == conds2)
      }
      if (n.trials2 == 1) {  ## for flexible behavior of sample
        inds2[[cond.i]] <- matrix(rep(which(bias.items[cond.i] == conds2), nsamp), nrow = n.trials2)
      } else {
        inds2[[cond.i]] <- replicate(nsamp, sample(which(bias.items[cond.i] == conds2), n.trials2), nsamp)
      }
      rownames(inds2[[cond.i]]) <- rep(bias.items[[cond.i]], n.trials2)
    }
    inds2 <- do.call(rbind, inds2)
    
    for (roi.i in seq_len(n.roi)) {  ## careful: roi.i is an index for as.numeric(rois)
      # roi.i <- 23
      
      ## get and apply mask for roi.i
      
      if (set.i == "mmp") {
        mask.i <- atlas[[set.i]] == roi.i 
      } else {
        mask.i <- masks[[roi.i]] == 1
      }
      
      roi.betas <- lapply(image.betas[[subj.i]], function(x1) apply(x1, "condition", function(x2) x2[mask.i]))  ## dirty
      
      ## averaging matrices:
      A1 <- model.matrix(~ 0 + dimnames(inds1)[[1]])  ## build dummy/indicator matrices
      A1 <- sweep(A1, 2, colSums(A1), "/")  ## scale cols
      colnames(A1) <- gsub("dimnames(inds1)[[1]]", "", colnames(A1), fixed = TRUE)  ## remove weird prefix
      
      ## averaging matrices:
      A2 <- model.matrix(~ 0 + dimnames(inds2)[[1]])  ## build dummy/indicator matrices
      A2 <- sweep(A2, 2, colSums(A2), "/")  ## scale cols
      colnames(A2) <- gsub("dimnames(inds2)[[1]]", "", colnames(A2), fixed = TRUE)  ## remove weird prefix
      
      rsarray.samp <- array(  ## for representational similarity matrices
        NA,
        dim = c(n.bias.items, n.bias.items, nsamp, length(type)),
        dimnames = list(
          .row = bias.items,
          .col = bias.items,
          samp = NULL,
          type = type
        )
      )
      
      # cl <- makeCluster(n.cores/4)
      # registerDoParallel(cl)
      time.start <- Sys.time()
      for (samp.i in seq_len(nsamp)) {
      # res <- foreach(
      #   samp.i = seq_len(nsamp),
      #   .inorder = FALSE
      # ) %dopar% {
        # samp.i = 1
      
        i1 <- inds1[, samp.i]
        roi.betas.i1 <- roi.betas[[1]][, i1]  ## downsample
        roi.betas.i1.bar <- roi.betas.i1 %*% A1  ## average by condition
        
        i2 <- inds2[, samp.i]
        roi.betas.i2 <- roi.betas[[2]][, i2]  ## downsample
        roi.betas.i2.bar <- roi.betas.i2 %*% A2  ## average by condition
        
        roi.betas.bar <- (roi.betas.i1.bar + roi.betas.i2.bar)/2
        
        rsarray.samp[, , samp.i, "crun"] <- cor(roi.betas.i1.bar, roi.betas.i2.bar)[bias.items, bias.items]
        rsarray.samp[, , samp.i, "wnr1"] <- cor(roi.betas.i1.bar)[bias.items, bias.items]
        rsarray.samp[, , samp.i, "wnr2"] <- cor(roi.betas.i2.bar)[bias.items, bias.items]
        rsarray.samp[, , samp.i, "both"] <- cor(roi.betas.bar)[bias.items, bias.items]
        
        # list(R_crun, R_wnr1, R_wnr2)
      }
      # stopCluster(cl)
      (time.end <- Sys.time() - time.start)
      
      rsarray[, , subj.i, roi.i, ] <- tanh(apply(atanh(rsarray.samp), c(".row", ".col", "type"), mean))

    
  }
    
    print(paste0(subj.i, ": subj ", fit.subjs[subj.i], " done!"), quote = FALSE)
    
    
  }  ## end subject loop
  
  
  ## store ----
    
  ## RSA results
  
  saveRDS(
    rsarray,
    here::here(
      "out", "rsa", "obsv",
      paste0("rsarray_", glm.name, "_downsampled_", set.i, ".rds")
      )
  )
  
  
}  ## end atlas loop

