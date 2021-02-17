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
files.dir.analysis <- files.dir.analysis[grep("fmriprep/", files.dir.analysis)]
fit.subjs <- unique(gsub("/results/.*", "", files.dir.analysis))

files.dir.analysis.orig <- list.files(dir.analysis, pattern = "stats_var", recursive = TRUE)  ## get fit.subjs
files.dir.analysis.orig <- files.dir.analysis.orig[grep("pro_bias_acc-only/", files.dir.analysis.orig)]
fit.subjs.orig <- unique(gsub("/results/.*", "", files.dir.analysis.orig))
fit.subjs.intersect <- intersect(fit.subjs, fit.subjs.orig)

## regs will be used to pull out (via string match) the statistic from the afni brick;
## thus must have the reg.suffix
regs <- c(bias.items, "pc50_c", "pc50_i", "sustained", "transient", "nuisance")
n.regs <- length(regs)
n.subj <- length(fit.subjs.intersect)
n.bias.items <- length(bias.items)



read_betas <- function(fname.nii, .regs = regs) {
  
  if (file.exists(fname.nii)) {
    ## dims of image.run.full [i, j, k, ???, regressor]
    image.full <- oro.nifti::readNIfTI(fname.nii, reorient = FALSE)
  } else stop("file nonexistant! ", paste0(fname.nii))
  
  ## get brick numbers for regressors
  
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
  
  if (any(is.na(brick.nums))) stop("brick nums equal zero! ", paste0(fit.subjs[subj.i]))
  
  ## put subset of images (only the relevant regressors) into one array
  
  ## image.betas with dims [i, j, k, regs]:
  image.betas <- image.full[, , , 1, brick.nums + 1]
  dimnames(image.betas) <- list(i = NULL, j = NULL, k = NULL, reg = .regs)
  
  image.betas
  
}



r <- numeric(length(fit.subjs.intersect))

is.brain <- c(atlas$mmp > 0)

for (subj.i in seq_along(fit.subjs.intersect)) {
  # subj.i <- 1
  
  ## get afni images: fmriprep
  

  dir.glms      <- file.path(dir.analysis, fit.subjs.intersect[subj.i], "results", "pro_bias_acc-only_fmriprep")
  dir.glms.orig <- file.path(dir.analysis, fit.subjs.intersect[subj.i], "results", "pro_bias_acc-only")
  
  fname      <- file.path(dir.glms, paste0("stats_", fit.subjs.intersect[subj.i], ".nii.gz"))
  fname.orig <- file.path(dir.glms.orig, paste0("stats_", fit.subjs.intersect[subj.i], ".nii.gz"))
  
  image.betas      <- read_betas(fname, regs)
  image.betas.orig <- read_betas(fname.orig, regs)
  
  image.betas.v <- array(image.betas, dim = c(75*90*75, 21))[is.brain, ]
  image.betas.orig.v <- array(image.betas.orig, dim = c(75*90*75, 21))[is.brain, ]
  
  # qcor(cor(image.betas.v, image.betas.orig.v))
  r[subj.i] <- tanh(mean(atanh(diag(cor(image.betas.v, image.betas.orig.v)))))

  print(paste0(subj.i, ": subj ", fit.subjs.intersect[subj.i], " done!"), quote = FALSE)
  
}  ## end subject loop



