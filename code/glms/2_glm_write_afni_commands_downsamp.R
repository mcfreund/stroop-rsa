stop("don't source me")
## about ----
## writes 3dDeconvolve and 3dREMLfit commands for GLMs.
## 
## to be executed from ccplinux1.
## 
## mike freund, 2019-02-20
## adapted for new project directory 2019-12-24

## setup ----

library(here)
library(mikeutils)
library(magrittr)
library(dplyr)
library(data.table)

source(here("code", "strings.R"))
stroop <- fread(here("data", "behavior-and-events_group201902.csv"))
subj_sum <- fread(here("data", "summary_group201902.csv"))


## funs ----


write.afni.cmds <- function(
  subj, session, .dmat, code.dir, out.folder, .dir.image = dir.nil.dmcc2.afni, .dir.analysis
) {
  # subj <- "132017"; session <- "proactive";
  # .dir.image = dir.nil.dmcc2.afni
  # .dir.analysis = file.path(dir.freund.external,"stroop-rsa", "afni-analysis_group-201902")
  # out.folder = "pro_bias_acc-only"
  
  ## TODO: validate inputs
  
  ## set up env:
  session.short <- paste0(toupper(substr(session, 1, 1)), substr(session, 2, 3))
  session.short.lc <- tolower(session.short)
  run.long <-c("1_AP", "2_PA")
  image.filename <- paste0(
    .dir.image, "/", subj, "/INPUT_DATA/Stroop/", session,
    "/lpi_scale_blur4_tfMRI_Stroop", session.short, run.long, ".nii.gz",
    collapse = " "
  )
  dir.results <- paste0(.dir.analysis, "/", subj, "/results/", out.folder)
  dir.input <- paste0(.dir.analysis, "/", subj, "/input/", session.short.lc)
  path.initial <- paste0(dir.input, "/", subj, "_", session.short.lc, "_")
  dir.glm <- paste0(code.dir, out.folder)
  
  ## create dir:
  if (!dir.exists(dir.glm)) dir.create(dir.glm, recursive = TRUE)
  
  ## write file
  ## binary connection, to write for bash:
  glm.filename <- file(paste0(dir.glm, "/glm_", subj), "wb")
  writeLines(
    c(
      ## begin file
      "#!/bin/bash",
      "\n",
      "#----- verify that the results directory does not yet exist",
      paste0("if [ -d ", dir.results, " ]; then"),
      paste0("  echo \"WARNING : output dir ", dir.results, " already exists\" >&2"),
      "  exit",
      "fi",
      paste0("mkdir -p ", dir.results),
      paste0("cd ", dir.results),
      "\n",
      
      ## begin deconvolve
      "3dDeconvolve \\",
      "  -local_times -x1D_stop -allzero_OK \\",
      paste0("  -input ", image.filename, " \\"),
      "  -polort A -float \\",
      paste0("  -censor ", dir.input, "/movregs_FD_mask.txt \\"),
      paste0("  -num_stimts ", nrow(dmat), " \\"),  ## plus two for "block" regressors
      
      ## begin design mat
      paste0(
        "  -", dmat[, ".time.opt"]," ", 1:nrow(dmat), " ", path.initial, dmat[, ".filename"], ".txt ", dmat[, ".fun"], 
        "  -stim_label ", 1:nrow(dmat), " ", dmat[, ".label"], " \\"
      ),
      paste0("  -ortvec ", dir.input, "/motion_demean_", session, ".1D movregs \\"),
      "  -x1D X.xmat.1D -xjpeg X.jpg -nobucket",
      "\n",
      
      ## begin REMLfit
      "3dREMLfit \\",
      "  -matrix X.xmat.1D \\",
      "  -GOFORIT 5 \\",
      paste0("  -input \"", image.filename, "\" \\", collapse = " "),
      paste0("  -Rvar stats_var_", subj, ".nii.gz \\"),
      paste0("  -Rbuck stats_", subj, ".nii.gz \\"),
      paste0("  -fout -tout -nobout -verb $*"),
      "\n",
      
      ## end file
      "wait",
      "echo -e \"EXIT STATUS : $?\"",
      "\n"
    ),
    glm.filename
  )
  close(glm.filename)
  unlink(glm.filename)
}  # end function


## write ----


fname.bias.con <- paste0(bias.items.con, "_acc-only_downsamp")
fname.bias.incon <- paste0(bias.items.incon, "_acc-only")

dmat.events <- cbind(
  .time.opt = "stim_times",
  .label    = c(bias.items.incon, bias.items.con, "pc50_i", "pc50_c", "nuisance", "congr_modelout"),
  .filename = c(fname.bias.incon, fname.bias.con, paste0(c("pc50_i", "pc50_c", "nuisance", "congr_modelout"), "_acc-only")),
  .fun      = "'BLOCK(1,1)'"
)
dmat.blocks <-   rbind(
  c("stim_times_AM1", "sustained", "sustained", "'dmBLOCK(1)'"),
  c("stim_times", "transient", "transient", "'TENTzero(0,16.8,8)'")
)
dmat <- rbind(dmat.events, dmat.blocks)


## baseline


stroop.bas.subj.nums <- stroop %>% filter(session == "bas") %>% split(list(.$subj)) %>% purrr::map_dbl(nrow)
stroop.bas.subjs <- names(stroop.bas.subj.nums)[stroop.bas.subj.nums == 216]

lapply(
  stroop.bas.subjs,
  function(x)
    write.afni.cmds(
      subj          = x,
      session       = "baseline", 
      .dmat         = dmat,
      code.dir      = "code/glms/",
      out.folder    = "bas_bias_acc-only_downsamp",
      .dir.image    = dir.nil.dmcc2.afni,
      .dir.analysis = here("glms")
    )
)
