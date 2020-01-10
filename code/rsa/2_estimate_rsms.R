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

#!/usr/bin/env Rscript

doc <- 
"Usage:
   2_estimate_rsms.R [-a <do_atlases> -m <do_masks> -u <univariate>]

Options:
   -a Conduct analysis using atlases (Glasser's Multi Modal Parcellation, and Gordon's RSFC communities)? [default: 1]
   -m Conduct analysis using user-specified masks? [default: 1]
   -u Estimate univariate statistics? [default: 1]

 ]"

opts <- docopt::docopt(doc)

do.atlas <- as.logical(as.integer(opts$a))
do.masks <- as.logical(as.integer(opts$m))
do.univa <- as.logical(as.integer(opts$u))

## setup ----

# library(here)
# library(mikeutils)
# library(magrittr)
# library(dplyr)
# library(abind)
# library(data.table)
# library(oro.nifti)

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
n.subj <- length(fit.subjs)
n.bias.items <- length(bias.items)

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
    # n.roi <- nrow(atlas.key[[atlas.i]])
    # roi.names <- atlas.key[[atlas.i]]$roi
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
  
  
  ## loop over subjs and rois ----
  
  for (subj.i in seq_along(fit.subjs)) {
    # subj.i <- 1
    
    ## get afni images
    
    dir.glms <- file.path(dir.analysis, fit.subjs[subj.i], "results", glm.name)
    fname.nii <- file.path(dir.glms, paste0("stats_", fit.subjs[subj.i], ".nii.gz"))
    
    if (file.exists(fname.nii)) {
      ## dims of image.run.full [i, j, k, ???, regressor]
      image.full <- readNIfTI(fname.nii, reorient = FALSE)
    } else { stop("file nonexistant! ", paste0(fname.nii)) }
    
    ## get brick numbers for regressors:
    
    brick.nums <- rep(NA, length(regs))  ## + 1 for sustained
    for (reg.i in regs) {
      # reg.i <- regs[1]
      
      image.label <- paste0(reg.i, "#0_Coef")
      
      if (nodename == "CCP-FREUND") {
        
        brick.str <- system2(
          "wsl",
          args = paste("/home/mcf/abin/3dinfo", "-label2index", image.label, win2lin(fname.nii)),
          stdout = TRUE
        )
        
      } else if (nodename == "ccplinux1"){
        
        brick.str <- system2(
          "/usr/local/pkg/linux_openmp_64/3dinfo",
          args = paste("-label2index", image.label, fname.nii),
          stdout = TRUE
        )  ## doesn't matter which run.
        
      }
      
      ## for error checking:
      
      has.error <- grepl("error", brick.str, ignore.case = TRUE)
      if (any(has.error)) stop("error loading brick nums: ", paste0(fit.subjs[subj.i], " ", reg.i))
      
      ## to remove function call that is included in output (when is.local.session):
      
      brick.str <- brick.str[!grepl("3dinfo: AFNI version", brick.str)]
      brick.num <- as.numeric(brick.str)
      brick.nums[which(reg.i == regs)] <- brick.num
      
    }
    
    if (any(is.na(brick.nums))) stop("brick nums equal zero! ", paste0(fit.subjs[subj.i]))
    
    ## put subset of images (only the relevant regressors) into one array:
    
    ## image.betas with dims [i, j, k, regs]
    image.betas <- image.full[, , , 1, brick.nums + 1]
    rm(image.full)  ## take out the garbage
    gc()
    dimnames(image.betas) <- list(i = NULL, j = NULL, k = NULL, reg = regs)
    
    ## generate rsm for each roi:
    
    rois <- atlas.key[[atlas.i]]$num.roi
    for (roi.i in seq_along(rois)) {  ## careful: roi.i is an index for as.numeric(rois)
      # roi.i <- 1
      
      for (hemi.i in hemis) {
        # hemi.i <- "r"
        
        ## for reading ease:
        this.atlas.hemi <- atlas[[atlas.i]][[hemi.i]]
        this.subj.roi.hemi <- paste(fit.subjs[subj.i], atlas.key[[atlas.i]][roi.i, "roi"], hemi.i, sep = "_")
        
        ## get roi:
        mask <- array(this.atlas.hemi %in% rois[roi.i], dim = dim(this.atlas.hemi))
        roi.betas <- apply(image.betas, "reg", function(slice.i) slice.i[mask])
        
        ## get rsm (pearson and euclidean):
        rsarray.pearson[, , subj.i, roi.i, hemi.i] <- cor(roi.betas[, bias.items])
        rsarray.euclidean[, , subj.i, roi.i, hemi.i] <- dist2mat(roi.betas[, bias.items]) / nrow(roi.betas)
      
      }
      
    }
    
    print(paste0(subj.i, ": subj ", fit.subjs[subj.i], " done!"), quote = FALSE)
    
  }

    ## store ----
    
    saveRDS(
      rsarray.pearson, 
      here(
        "out", "rsa", "obsv", 
        paste0("rsarray_", glm.name, "_", atlas.i, "_pearson.rds")
        )
      )
    
    saveRDS(
      rsarray.euclidean, 
      here(
        "out", "rsa", "obsv", 
        paste0("rsarray_", glm.name, "_", atlas.i, "_euclidean.rds")
      )
    )

}  ## end atlas loop
