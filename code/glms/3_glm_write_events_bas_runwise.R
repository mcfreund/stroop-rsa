stop("don't source me")
## about ----
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
library(purrr)

source(here("code", "strings.R"))

## split all stimtimes files in bas dir

dir.to.write.in <- here("glms")

subjs <- list.dirs(dir.to.write.in, recursive = FALSE, full.names = FALSE)

dirs.input <- file.path(dir.to.write.in, subjs, "input", "bas")
dirs.input <- dirs.input[dir.exists(dirs.input)]


onsets2file <- function(.onsets, .fname) {
  .fname <- paste0(.fname, ".txt")
  .fout <- file(.fname, "wt")
  cat(.onsets, file = .fout)
  close(.fout)
  unlink(.fout)
}

## split stimtimes files, and copy movregs

was.copied <- setNames(logical(length(subjs)), subjs)

for (dir.i in seq_along(dirs.input)) {
  # dir.i = 1
  
  name.dir.i <- dirs.input[dir.i]
  
  fnames.i <- list.files(name.dir.i)
  fnames.stimtimes.i <- fnames.i[grep(paste0("^", subjs[dir.i]), fnames.i)]  ## discard movregs
  fnames.stimtimes.i <- fnames.i[grep("acc-only\\.txt$", fnames.i)]  ## discard run-wise times
  if (length(fnames.stimtimes.i)!= 19) stop("bad length")
  
  for (name.stimtime.i in fnames.stimtimes.i) {
    # name.stimtime.i <- fnames.stimtimes.i[1]
    
    fname.i <- file.path(name.dir.i, name.stimtime.i)
    
    run1 <- readLines(fname.i, 1)
    run2 <- readLines(fname.i, 2)
    
    onsets2file(run1, paste0(gsub(".txt", "", fname.i), "_run1"))
    onsets2file(run2, paste0(gsub(".txt", "", fname.i), "_run2"))
    
  }
  
  
  fname.movregs <- file.path(
    "/data/nil-external/ccp/freund/sub-subj-glms/freund/sub-subj-glms/out/AFNI_ANALYSIS",
    subjs[dir.i],
    "INPUT_DATA/Stroop/baseline"
  )

  fname.movregs <- c(
    file.path(fname.movregs, "movregs_FD_mask_run1.txt"),
    file.path(fname.movregs, "movregs_FD_mask_run2.txt"),
    file.path(fname.movregs, "motion_demean_baseline_run2.1D"),
    file.path(fname.movregs, "motion_demean_baseline_run2.1D")
  )

  was.copied[dir.i] <- all(file.copy(fname.movregs, name.dir.i, overwrite = TRUE))
  
}

sum(!was.copied)
sum(was.copied)  ## full baseline set
subjs.bas <- names(was.copied)[was.copied]

write.csv(subjs.bas)