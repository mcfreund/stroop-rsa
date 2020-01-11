#!/usr/bin/env Rscript

## about ----
## 
## reads in afni images (beta estimates from GLM) into a list.
## this list contains one element per subject per parcel per hemisphere.
## given a parcellation atlas or mask, similarity measures are then calculated from this list, and saved as
##  .RDS files to stroop-rsa/out/rsa/.
## additionaly saved are mean values for conducting a univariate analysis.
## 
## this script is configured to run as an executable on a *nix system.
## however, it can also be run locally (on mike's lenovo) in an interactive session.
## 
## mike freund, 2019-12-24

## TODO
## function for matching brick string

doc <- 
"Usage:
   2_estimate_rsms.R [-a <do_atlases> -m <do_masks> -u <univariate>]

Options:
   -a Conduct analysis using atlases (Glasser's Multi Modal Parcellation, and Gordon's RSFC communities)? [default: 0]
   -m Conduct analysis using user-specified masks? [default: 0]
   -u Estimate univariate statistics? [default: 0]

 ]"

opts <- docopt::docopt(doc)

do.atlas <- as.logical(as.integer(opts$a))
do.masks <- as.logical(as.integer(opts$m))
do.univa <- as.logical(as.integer(opts$u))

## defaults for interactive use (e.g., debugging, ...):
if (interactive()) {
  do.atlas <- TRUE
  do.masks <- TRUE
  do.univa <- TRUE
}

if (!any(do.atlas, do.masks, do.univa)) stop(paste0("you must do something!"))

## setup ----

source(here::here("code", "strings.R"))
if (do.atlas) source(here::here("code", "read_atlases.R"))
if (do.masks) source(here::here("code", "read_masks.R"))

## paths, vars

dir.analysis <- here::here("glms")
glm.name <- "pro_bias_acc-only"
files.dir.analysis <- list.files(dir.analysis, pattern = "stats_var", recursive = TRUE)  ## get fit.subjs
files.dir.analysis <- files.dir.analysis[grep(glm.name, files.dir.analysis)]
fit.subjs <- unique(gsub("/results/.*", "", files.dir.analysis))

## regs will be used to pull out (via string match) the statistic from the afni brick;
## thus must have the reg.suffix
regs <- c(bias.items, "pc50_c", "pc50_i", "sustained", "transient", "nuisance")
n.regs <- length(regs)
n.subj <- length(fit.subjs)
n.bias.items <- length(bias.items)

## this variable defines the outermost loop.
## each iteration collates RSMs into a single array, and saves it as a single .rds file.
sets.of.rois <- character(0)
if (do.atlas) sets.of.rois <- c(sets.of.rois, names(atlas))
if (do.masks) sets.of.rois <- c(sets.of.rois, "masks")


## loop over sets of ROIs ----

for (set.i in sets.of.rois) {
  # set.i = "mmp"
  
  ## get numbers and create storage objects
  
  if (set.i != "masks") {
    n.roi <- nrow(atlas.key[[set.i]])
    roi.names <- atlas.key[[set.i]]$roi
  } else {
    n.roi <- length(masks)
    roi.names <- names(masks)
  }
  
  rsarray.pearson <- array(  ## for representational similarity matrices
    NA,
    dim = c(n.bias.items, n.bias.items, n.subj, n.roi),
    dimnames = list(
      .row = bias.items,
      .col = bias.items,
      subj = fit.subjs,
      roi  = roi.names
    )
  )
  rsarray.euclidean <- rsarray.pearson
  
  ## for tallying voxels:
  voxels.silent <- matrix(NA, nrow = n.roi, ncol = n.subj)  ## for num unresponsive / roi
  dimnames(voxels.silent) <- list(roi = roi.names, subj = fit.subjs)
  voxels.number <- numeric(n.roi)  ## for total number of voxels / roi
  
  if (do.univa) {  ### for saving unvariate stats (means)
    roi.means <- matrix(NA, nrow = n.roi, ncol = n.regs)
    dimnames(roi.means) <- list(roi = roi.names, reg = regs)
  }
  
  ## loop over subjs and rois ----
  
  for (subj.i in seq_along(fit.subjs)) {
    # subj.i <- 1
    
    ## get afni images
    
    dir.glms <- file.path(dir.analysis, fit.subjs[subj.i], "results", glm.name)
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
      
      ## TODO: move this into wrapper function?
      if (nodename == "CCP-FREUND") {
        
        brick.str <- system2(
          "wsl",
          args = paste("/home/mcf/abin/3dinfo", "-label2index", image.label, fname.nii),
          stdout = TRUE
        )
        
      } else if (nodename == "ccplinux1"){
        
        brick.str <- system2(
          "/usr/local/pkg/linux_openmp_64/3dinfo",
          args = paste("-label2index", image.label, fname.nii),
          stdout = TRUE
        )  ## doesn't matter which run.
        
      }
      
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
      
      if (set.i != "masks") mask.i <- atlas[[set.i]] == roi.i else mask.i <- masks[[roi.i]] == 1
      
      roi.betas <- apply(image.betas, "reg", function(.) .[mask.i])
      
      ## tally number of unresponsive voxels
      
      is.all.zero <- vapply(rowSums(roi.betas), function(.) isTRUE(all.equal(., 0)), logical(1))
      voxels.silent[n.roi, subj.i] <- sum(is.all.zero)
      
      ## tally number of voxels (but only do once; same for all subjects)
      
      n.voxels <- nrow(roi.betas)
      if (subj.i == 1) voxels.number[n.roi] <- n.voxels
      
      ## get rsm (pearson and euclidean)
      
      rsarray.pearson[, , subj.i, roi.i] <- cor(roi.betas[, bias.items])
      rsarray.euclidean[, , subj.i, roi.i] <- mikeutils::dist2mat(roi.betas[, bias.items]) / n.voxels
      
      ## get univariate stats (across-voxel means)
      
      if (do.univa) roi.means[roi.i, ] <- apply(roi.betas, "reg", mean)
    
    }
    
    print(paste0(subj.i, ": subj ", fit.subjs[subj.i], " done!"), quote = FALSE)
    
  }  ## end subject loop

  
  ## store ----
    
  ## RSA results
  
  saveRDS(
    rsarray.pearson, 
    here::here(
      "out", "rsa", "obsv", 
      paste0("rsarray_", glm.name, "_", set.i, "_pearson.rds")
      )
    )
  
  saveRDS(
    rsarray.euclidean, 
    here::here(
      "out", "rsa", "obsv", 
      paste0("rsarray_", glm.name, "_", set.i, "_euclidean.rds")
    )
  )
  
  ## univariate results
  
  if (do.univa) {

    saveRDS(
      roi.means, 
      here::here(
        "out", "rsa", "obsv",  ## not an RSA, but save in ./out/rsa/ for consistency...
        paste0("roi-means_", glm.name, "_", set.i, ".rds")
      )
    )
    
  }
  
  ## voxel information
  
  voxels.silent <- data.table::as.data.table(voxels.silent)
  voxels.silent$roi <- roi.names
  voxels.silent <- data.table::melt(voxels.silent, id.vars = "roi", variable.name = "subj", value.name = "n.silent")
  data.table::fwrite(
    voxels.silent,
    here::here(
      "out", "summaries", paste0("voxel-counts_unresponsive_", set.i, ".csv")
    )
  )
  
  data.table::fwrite(
    data.table::data.table(roi = roi.names, n.total = voxels.number),
    here::here(
      "out", "summaries", paste0("voxel-counts_total_", set.i, ".csv")
    )
  )
  

}  ## end atlas loop
