## about ----
## mike freund, 2019-02-23
## 
## reads in afni images (beta estimates from GLM) into a list.
## this list contains one element per subject per parcel per hemisphere.
## similarity measures are then calculated from this list, then saved as an .RData file via boxr::box_save() to
##  ./stroop-rsa/data.
## the beta list is also optionally saved to this location.
## this is done separately for mmp and gordon parcels, yielding two .RData files.
## the subsequent script, 4_write_rsms.R,  will load these files, calculate RSMs,
## and write the results to  ./stroop-rsa/data as csvs.
## the purpose of this script is to get the data to a location (box) and format (.RData) that can be 
## accessed via boxdrive (locally)---enabling many files to be (more) efficiently written compared to  via boxr.

## updated 2019-03-14 :: conducting RSA on multi-parcel "clusters"
## updated 2019-03-17 :: run only for held out subjects


## set up env. ----

do.clusters <- TRUE
# do.held.out <- TRUE

## dependencies

library(boxr)
box_auth()
library(here)
library(dplyr)
library(abind)
library(data.table)
library(oro.nifti)
library(gifti)
library(plot3D)

source(here("..", "gen", "funs", "_funs.R"))
source(here("..", "gen", "funs", "_get_dirs_remote.R"))
do.mb <- c("mb4" = TRUE, mb8 = FALSE)  ## to skip prompts in 
do.read.atlas <- c("mmp" = TRUE, gordon = TRUE)
source(here("..", "gen", "funs", "_get_atlas.R"))
source(here("r", "group-201902", "_get_misc_vars.R"))
source(here("r", "group-201902", "_read_sheets_remote.R"))


## vars

dir.analysis <- file.path(dir.freund.external, "stroop-rsa", "afni-analysis_group-201902")
glm.name <- "pro_bias_acc-only"
hemis <- c("l", "r")
stats <- c("pearson", "euclidean")  ## which similarity statistics should be calculated?
save.betas <- FALSE  ## should the betas be written?
n.bias.items <- length(bias.items)
n.hemi <- length(hemis)  ## duh, but for clarity
n.stats <- length(stats)

## roi clusters:
## (multi-parcel regions that show task dimension coding)
clusters <- list(
  evis.l     = paste0(c("V1", "V2", "V3"), "_l"),
  evis.r     = paste0(c("V1", "V2", "V3"), "_r"),
  vvis.l     = paste0(c("FFC", "VVC", "V8", "VMV3"), "_l"),
  vvis.r     = paste0(c("FFC", "VVC", "V8", "VMV3", "VMV2"), "_r"),
  ips.l      = paste0(c("IP0", "IP1", "LIPd", "VIP", "LIPv", "AIP", "7PC"), "_l"),
  ips.r      = paste0(c("IP0", "IP1", "IP2", "IPS1", "MIP", "LIPd", "LIPv", "AIP", "7PC"), "_r"),
  premotor.l = paste0(c("FEF", "55b", "PEF", "6v", "6r", "43", "6d"), "_l"),
  premotor.r = paste0(c("FEF", "55b", "PEF", "6v", "6r", "43", "6a"), "_r"),
  dlpfc.l    = paste0(c("p9-46v", "8C", "8Av", "i6-8"), "_l"),
  dlpfc.r    = paste0(c("p9-46v", "8C", "8Av", "i6-8"), "_r"),
  mpfc.l     = paste0(c("a32pr", "p32pr", "8BM", "SCEF"), "_l"),
  mpfc.r     = paste0(c("a32pr", "p32pr", "8BM", "SCEF"), "_r"),
  ins.l      = paste0(c("FOP4", "FOP5", "FOP3"), "_l"),
  ins.r      = paste0(c("FOP4", "FOP5", "FOP3"), "_r"),
  ifc.l      = paste0(c("44", "45", "IFSa", "IFSp", "p47r", "p47l"), "_l"),
  ifc.r      = paste0(c("44", "45", "IFSa", "IFSp", "p47r", "p47l"), "_r"),
  lang.l     = paste0(c("PSL", "55b", "44", "45", "SFL"), "_l"),
  lang.r     = paste0(c("PSL", "55b", "44", "45", "SFL"), "_r")
)
cluster.regions <- unique(gsub("(.*)\\..*", "\\1", names(clusters)))

## strings
## regs will be used to pull out (via string match) the statistic from the afni brick;
## thus must have the reg.suffix
regs <- c(bias.items, "pc50_c", "pc50_i", "sustained", "transient", "nuisance")
files.dir.analysis <- list.files(dir.analysis, pattern = "stats_var", recursive = TRUE)  ## get fit.subjs
files.dir.analysis <- files.dir.analysis[grep(glm.name, files.dir.analysis)]
fit.subjs <- unique(gsub("/results/.*", "", files.dir.analysis))
# if (do.held.out) fit.subjs <- group201902.held.out


## loop atlas ----


for (atlas.i in names(atlas.key)) {
  # atlas.i = "mmp"
  
  if (do.clusters) if(atlas.i == "gordon") next  ## multiparcel clusters only for mmp
  
  
  ## numbers
  n.roi <- ifelse(do.clusters, length(cluster.regions), nrow(atlas.key[[atlas.i]]))
  if (do.clusters) {
    roi.names <- cluster.regions
  } else {
    roi.names <- atlas.key[[atlas.i]]$roi
  }
  n.subj <- length(fit.subjs)
  
  ## initialize storage objects:
  if (save.betas) {
    betalist <- vector("list", n.roi * n.subj * n.hemi)   ## for beta values
    names(betalist) <- combo_paste3(fit.subjs, atlas.key[[atlas.i]]$roi, hemis)
  }
  rsarray <- array(  ## for representational similarity matrices
    NA,
    dim = c(n.bias.items, n.bias.items, n.subj, n.roi, n.hemi),
    dimnames = list(
      r    = bias.items,
      c    = bias.items,
      subj = fit.subjs,
      roi  = roi.names,
      hemi = hemis
    )
  )
  
  ## create rsarray for each similarity measure
  for (ii in seq_along(stats)) assign(paste0("rsarray", ".", stats[ii]), rsarray)
  
  ## loop over subjs and rois ----
  
  for (subj.i in seq_along(fit.subjs)) {
    # subj.i <- 1
    
    ## get afni images:
    dir.glms <- file.path(dir.analysis, fit.subjs[subj.i], "results", glm.name)
    fname.nii <- file.path(dir.glms, paste0("stats_", fit.subjs[subj.i], ".nii.gz"))
    if (file.exists(fname.nii)) {
      ## dims of image.run.full [i, j, k, ???, regressor]
      image.full <- readNIfTI(fname.nii, reorient = FALSE)
    } else {
      stop("file nonexistant! ", paste0(fname.nii))
    }
    
    ## get brick numbers for regressors:
    brick.nums <- rep(NA, length(regs))  ## + 1 for sustained
    for (reg.i in regs) {
      # reg.i <- regs[1]
      image.label <- paste0(reg.i, "#0_Coef")
      if (is.local.session) {
        brick.str <- system2(
          "wsl",
          args = paste(
            "/home/mcf/abin/3dinfo", "-label2index", image.label, win2lin(fname.nii)
          ),
          stdout = TRUE
        )
      } else {
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
    dimnames(image.betas) <- list(i = NULL, j = NULL, k = NULL, reg = regs)
    
    
    if (do.clusters) {
      for(cluster.i in seq_along(clusters)) {
        # cluster.i = 1
        hemi.i <- gsub(".*\\.(.*)", "\\1", names(clusters)[cluster.i])
        region.i <- gsub("(.*)\\..*", "\\1", names(clusters)[cluster.i])
        rois.i <- gsub("(.*)_.*", "\\1", clusters[[cluster.i]])
        num.rois.i <- atlas.key$mmp$num.roi[atlas.key$mmp$roi %in% rois.i]
        this.atlas.hemi <- atlas[[atlas.i]][[hemi.i]]
        # this.subj.roi.hemi <- paste(fit.subjs[subj.i], atlas.key[[atlas.i]][roi.i, "roi"], hemi.i, sep = "_")
        ## get roi:
        mask <- array(this.atlas.hemi %in% num.rois.i, dim = dim(this.atlas.hemi))
        roi.betas <- apply(image.betas, "reg", function(slice.i) slice.i[mask])
        if (save.betas) betalist[[this.subj.roi.hemi]] <- roi.betas
        ## get rsm (pearson and euclidean):
        if ("pearson" %in% stats)  rsarray.pearson[, , subj.i, region.i, hemi.i] <- cor(roi.betas[, bias.items])
        if ("euclidean" %in% stats) rsarray.euclidean[, , subj.i, region.i, hemi.i] <- tdist(roi.betas[, bias.items])
      }
      print(paste0(subj.i, ": subj ", fit.subjs[subj.i], " done!"), quote = FALSE)
    } else {
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
        if (save.betas) betalist[[this.subj.roi.hemi]] <- roi.betas
        ## get rsm (pearson and euclidean):
        if ("pearson" %in% stats)  rsarray.pearson[, , subj.i, roi.i, hemi.i] <- cor(roi.betas[, bias.items])
        if ("euclidean" %in% stats) rsarray.euclidean[, , subj.i, roi.i, hemi.i] <- tdist(roi.betas[, bias.items])
      }
    }
    print(paste0(subj.i, ": subj ", fit.subjs[subj.i], " done!"), quote = FALSE)
  }
    }
  
  ## store ----
  
  # NB:
  ## the goal is to store these data on box in .csv format in ./stroop-rsa/data, in subject and parcel-specific folders.
  ## if this subject and parcel directory structrue does not exist, it should be created.
  ## but box does not talk to linux, except via web.
  ## boxr is an api tool for R, but isn't great for creating or searching through large directory structures.
  ## so, i write a single .RData file to ./stroop-rsa/data using boxr.
  ## i'll then use another, local, R script (and box drive) to read this .RData file, and create and write 
  ## the necessary directory structure.
  
  file.suffix <- ifelse(do.clusters, "_multiparcel.RData", ".RData")
  # file.suffix <- ifelse(do.held.out, paste0("_held-out", file.suffix), file.suffix)
  
  if ("pearson" %in% stats) {
    box_save(
      rsarray.pearson,
      dir_id = boxid.strooprsa.data,
      file_name = paste0("rsarray_pearson_", atlas.i, "_group-201902_", glm.name, file.suffix)
    )
  }
  if ("euclidean" %in% stats) {
    box_save(
      rsarray.euclidean,
      dir_id = boxid.strooprsa.data,
      file_name = paste0("rsarray_euclidean_", atlas.i, "_group-201902_", glm.name, file.suffix)
    )
  }
  ## commented 2019-02-24: too slow to store betas; will just move forward with RSMs
  if (save.betas) {
    box_save(
      betalist,
      dir_id = boxid.strooprsa.data,
      file_name = paste0("betalist_", atlas.i, "_group-201902_", glm.name, file.suffix)
      )
  }
  
}  ## end atlas loop


