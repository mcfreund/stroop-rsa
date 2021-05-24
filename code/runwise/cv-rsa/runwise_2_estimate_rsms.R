## about ----
## 
## reads in afni images (beta estimates from GLM) into a list.
## this list contains one element per subject per parcel per hemisphere.
## given a  mask, similarity measures are then calculated, and saved as
##  .RDS files to stroop-rsa/out/rsa/.
## additionaly saved are mean values for conducting a univariate analysis.
## 
## mike freund, 2019-12-24


## setup ----

source(here::here("code", "packages.R"))
source(here::here("code", "strings.R"))
source(here::here("code", "read_atlases.R"))
source(here::here("code", "read_masks.R"))

## functions

read_betas <- function(.dir.glms, .regs) {
  
  image.betas <- array(
    NA, 
    dim = c(75, 90, 75, length(.regs), 2), 
    dimnames = list(i = NULL, j = NULL, k = NULL, reg = .regs, run = 1:2)
    )
  
  for (run.i in 1:2) {
    # run.i = 1
    
    fname.nii <- file.path(paste0(.dir.glms, run.i), paste0("stats_", fit.subjs[subj.i], "_run", run.i, ".nii.gz"))
    if (file.exists(fname.nii)) {
      ## dims of image.run.full [i, j, k, ???, regressor]
      image.full <- oro.nifti::readNIfTI(fname.nii, reorient = FALSE)
    } else stop("file nonexistant! ", paste0(fname.nii))
    
    brick.nums <- rep(NA, length(.regs))  ## + 1 for sustained
    for (reg.i in .regs) {
      # reg.i <- regs[1]
      
      image.label <- paste0(reg.i, "#0_Coef")
      
      brick.str <- mikeutils::afni("3dinfo", paste("-label2index", image.label, fname.nii))

      ## for error checking
      
      has.error <- grepl("error", brick.str, ignore.case = TRUE)
      if (any(has.error)) stop("error loading brick nums: ", paste0(fit.subjs[subj.i], " ", reg.i))
      
      ## to remove function call that is included in output (when is.local.session)
      
      brick.str <- brick.str[!grepl("3dinfo: AFNI version", brick.str)]
      brick.num <- as.numeric(brick.str)
      brick.nums[which(reg.i == .regs)] <- brick.num
      
    }
    
    ## get brick numbers for regressors
    
    if (any(is.na(brick.nums))) stop("brick nums equal zero! ", paste0(fit.subjs[subj.i]))
    
    ## put subset of images (only the relevant regressors) into one array
    
    ## image.betas with dims [i, j, k, regs]:
    image.betas[, , , , run.i] <- image.full[, , , 1, brick.nums + 1]
    
  }
  
  image.betas
  
}



distance_cv <- function(B1, B2, m, regressors) {
  # B1 = B[, , 1]; B2 = B[, , 2]
  
  D <- colSums(t(m %*% B1 * m %*% B2)) / ncol(B1)
  matrix(D, ncol = length(regressors), dimnames = list(.row = regressors, .col = regressors))
  
}


## paths, vars

dir.analysis <- here::here("glms")
files.dir.analysis <- list.files(dir.analysis, pattern = "stats_var", recursive = TRUE)  ## get fit.subjs
files.dir.analysis <- files.dir.analysis[grep("pro_bias_acc-only_fmriprep_run[1-2]/", files.dir.analysis)]
fit.subjs <- unique(gsub("/results/.*", "", files.dir.analysis))

## regs will be used to pull out (via string match) the statistic from the afni brick;
## thus must have the reg.suffix
regs <- c(bias.items, "pc50_c", "pc50_i", "sustained", "transient", "nuisance")
n.regs <- length(regs)
n.subj <- length(fit.subjs)
n.bias.items <- length(bias.items)

cmat <- mikeutils::contrast_matrix(length(bias.items), bias.items)  ## contrast matrix

## loop over sets of ROIs ----
## each iteration collates RSMs into a single array, and saves it as a single .rds file.

for (set.i in sets.of.rois) {
  # set.i = "mmp"
  
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
  

  ## loop over subjs and rois ----
  
  for (subj.i in seq_along(fit.subjs)) {
    # subj.i <- 1
    
    ## get afni images
    
    dir.glms <- file.path(dir.analysis, fit.subjs[subj.i], "results", "pro_bias_acc-only_fmriprep_run")
    image.betas <- read_betas(dir.glms, regs)
    
    
    ## generate rsm for each roi
    
    for (roi.i in seq_len(n.roi)) {  ## careful: roi.i is an index for as.numeric(rois)
      # roi.i <- 1
      
      ## get and apply mask for roi.i
      
      if (set.i == "mmp") {
        mask.i <- atlas[[set.i]] == roi.i 
      } else {
        mask.i <- masks[[roi.i]] == 1
      }
      
      roi.betas <- apply(image.betas, c("reg", "run"), function(.) .[mask.i])
      
      ## estimate similarity matrices:
      
      B <- aperm(roi.betas, c(2, 1, 3))  
      names(dimnames(B)) <- c("condition", "vertex", "run")
      
      B1 <- B[bias.items, , 1]
      B2 <- B[bias.items, , 2]
      
      D             <- distance_cv(B1, B2, cmat, bias.items)
      # D       <- distance_cv(t(scale(t(B1))), t(scale(t(B2))), cmat, bias.items)
      
      rsarray[, , subj.i, roi.i] <- D

    }
    
    print(paste0(subj.i, ": subj ", fit.subjs[subj.i], " done!"), quote = FALSE)
    
  }  ## end subject loop

  
  ## store ----
    
  ## RSA results
  
  saveRDS(
    rsarray, 
    here::here(
      "out", "rsa", "obsv", 
      paste0("rsarray_pro_bias_acc-only_fmriprep_runwise_cv-euclidean_", set.i, ".rds")
      )
    )
  

}  ## end atlas loop

