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
  

  ## loop over subjs and rois ----
  
  for (subj.i in seq_along(fit.subjs)) {
    # subj.i <- 1
    
    ## get afni images
    
    dir.glms <- file.path(dir.analysis, fit.subjs[subj.i], "results", glm.name)
    image.betas <- read_betas(dir.glms)  ## returns list, as arrays may have differing numbers of dimensions by run
    
    
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

      ## estimate similarity matrix:
      
      R <- cor(roi.betas[[1]], roi.betas[[2]])  ## cross-run: rows = run1, cols = run2
      
      ## aggregate by condition over trials
      
      A1 <- model.matrix(~ 0 + dimnames(R)[[1]])  ## build dummy/indicator matrices
      A2 <- model.matrix(~ 0 + dimnames(R)[[2]])
      
      A1 <- sweep(A1, 2, colSums(A1), "/")  ## scale cols
      A2 <- sweep(A2, 2, colSums(A2), "/")
      
      R_ave <- tanh(t(A1) %*% atanh(R) %*% A2)  ## average
      
      rownames(R_ave) <- gsub("dimnames(R)[[1]]", "", rownames(R_ave), fixed = TRUE)  ## remove weird prefix
      colnames(R_ave) <- gsub("dimnames(R)[[2]]", "", colnames(R_ave), fixed = TRUE)
    
      
      rsarray[, , subj.i, roi.i] <- R_ave[bias.items, bias.items]

    }
    
    print(paste0(subj.i, ": subj ", fit.subjs[subj.i], " done!"), quote = FALSE)
    
  }  ## end subject loop

  
  ## store ----
    
  ## RSA results
  
  saveRDS(
    rsarray, 
    here::here(
      "out", "rsa", "obsv", 
      paste0("rsarray_", glm.name, "_", set.i, ".rds")
      )
    )
  

}  ## end atlas loop

