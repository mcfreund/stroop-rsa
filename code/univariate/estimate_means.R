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

## paths, vars

dir.analysis <- here::here("glms")
files.dir.analysis <- list.files(dir.analysis, pattern = "stats_var", recursive = TRUE)  ## get fit.subjs
files.dir.analysis <- files.dir.analysis[grep("pro_bias_acc-only/", files.dir.analysis)]
fit.subjs <- unique(gsub("/results/.*", "", files.dir.analysis))

## regs will be used to pull out (via string match) the statistic from the afni brick;
## thus must have the reg.suffix
regs <- c(bias.items, "pc50_c", "pc50_i", "sustained", "transient", "nuisance")
n.regs <- length(regs)
n.subj <- length(fit.subjs)
n.bias.items <- length(bias.items)


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
  
  means <- array(  ## for representational similarity matrices
    NA,
    dim = c(n.bias.items, n.subj, n.roi),
    dimnames = list(
      cond = bias.items,
      subj = fit.subjs,
      roi  = roi.names
    )
  )
  
  
  ## loop over subjs and rois ----
  
  for (subj.i in seq_along(fit.subjs)) {
    # subj.i <- 1
    
    ## get afni images
    
    dir.glms <- file.path(dir.analysis, fit.subjs[subj.i], "results", "pro_bias_acc-only")
    fname.nii <- file.path(dir.glms, paste0("stats_", fit.subjs[subj.i], ".nii.gz"))
    
    if (file.exists(fname.nii)) {
      ## dims of image.run.full [i, j, k, ???, regressor]
      image.full <- oro.nifti::readNIfTI(fname.nii, reorient = FALSE)
    } else stop("file nonexistant! ", paste0(fname.nii))
    
    ## get brick numbers for regressors
    
    brick.nums <- rep(NA, length(regs))  ## + 1 for sustained
    for (reg.i in regs) {
      # reg.i <- regs[1]
      
      image.label <- paste0(reg.i, "#0_Coef")
      
      brick.str <- mikeutils::afni("3dinfo", paste("-label2index", image.label, fname.nii))
      
      ## for error checking
      
      has.error <- grepl("error", brick.str, ignore.case = TRUE)
      if (any(has.error)) stop("error loading brick nums: ", paste0(fit.subjs[subj.i], " ", reg.i))
      
      ## to remove function call that is included in output (when is.local.session)
      
      brick.str <- brick.str[!grepl("3dinfo: AFNI version", brick.str)]
      brick.num <- as.numeric(brick.str)
      brick.nums[which(reg.i == regs)] <- brick.num
      
    }
    
    if (any(is.na(brick.nums))) stop("brick nums equal zero! ", paste0(fit.subjs[subj.i]))
    
    ## put subset of images (only the relevant regressors) into one array
    
    ## image.betas with dims [i, j, k, regs]:
    image.betas <- image.full[, , , 1, brick.nums + 1]
    rm(image.full)
    gc()  ## take out the garbage
    dimnames(image.betas) <- list(i = NULL, j = NULL, k = NULL, reg = regs)
    
    ## generate rsm for each roi
    
    for (roi.i in seq_len(n.roi)) {  ## careful: roi.i is an index for as.numeric(rois)
      # roi.i <- 1
      
      ## get and apply mask for roi.i
      
      if (set.i == "mmp") {
        mask.i <- atlas[[set.i]] == roi.i 
      } else {
        mask.i <- masks[[roi.i]] == 1
      }
      
      roi.betas <- apply(image.betas, "reg", function(.) .[mask.i])
      
      ## get mean

      means[, subj.i, roi.i] <- colMeans(roi.betas[, bias.items])
      
    }
    
    print(paste0(subj.i, ": subj ", fit.subjs[subj.i], " done!"), quote = FALSE)
    
  }  ## end subject loop
  
  
  ## store ----
  
  saveRDS(
    means, 
    here::here(
      "out", "rsa", "obsv",  ## not an RSA, but save in ./out/rsa/ for consistency...
      paste0("means_pro_bias_acc-only_", set.i, ".rds")
    )
  )
  
  
}  ## end atlas loop

