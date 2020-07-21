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

dirs.input <- file.path(dir.to.write.in, subjs, "input", "pro")
dirs.input <- dirs.input[dir.exists(dirs.input)]


onsets2file <- function(.onsets, .fname) {
  .fname <- paste0(.fname, ".txt")
  .fout <- file(.fname, "wt")
  cat(.onsets, file = .fout)
  close(.fout)
  unlink(.fout)
}

onsets4afni <- function(x) {
  ## takes a vector of onsets and puts them in format for afni
  x <- x[!is.na(x)]
  if (class(x) != "numeric") stop("needs numeric!")
  if (length(x) < 1) {
    x <- "*\n"
  } else {
    x <- paste0(paste(x, collapse = " "), " \n")
  }
  return(x)
}

## split stimtimes files, and copy movregs

was.copied <- setNames(logical(length(subjs)), subjs)

for (dir.i in seq_along(dirs.input)) {
  # dir.i = 1
  
  name.dir.i <- dirs.input[dir.i]
  
  fnames.i <- list.files(name.dir.i)
  fnames.stimtimes.i <- fnames.i[grep(paste0("^", subjs[dir.i]), fnames.i)]  ## discard movregs
  fnames.stimtimes.i <- fnames.i[grep("sustained\\.txt|transient\\.txt|acc-only\\.txt$", fnames.i)]  ## discard run-wise times
  if (length(fnames.stimtimes.i) != 21) stop("bad length")
  
  for (name.stimtime.i in fnames.stimtimes.i) {
    # name.stimtime.i <- fnames.stimtimes.i[1]
    
    fname.i <- file.path(name.dir.i, name.stimtime.i)
    
    # run1 <- readLines(fname.i, 1)
    # run2 <- readLines(fname.i, 2)
    
    stimtimes <- readChar(fname.i, file.info(fname.i)$size)
    stimtimes <- strsplit(stimtimes, split = "\n")[[1]]
    
    if (length(stimtimes) != 2) stop("bad length stimtimes")
    
    onsets2file(stimtimes[1], paste0(gsub(".txt", "", fname.i), "_run1"))
    onsets2file(stimtimes[2], paste0(gsub(".txt", "", fname.i), "_run2"))
    
  }
  
  
  fname.movregs <- file.path(
    "/data/nil-external/ccp/freund/sub-subj-glms/freund/sub-subj-glms/out/AFNI_ANALYSIS",
    subjs[dir.i],
    "INPUT_DATA/Stroop/proactive"
  )

  fname.movregs <- c(
    file.path(fname.movregs, "movregs_FD_mask_run1.txt"),
    file.path(fname.movregs, "movregs_FD_mask_run2.txt"),
    file.path(fname.movregs, "motion_demean_proactive_run1.1D"),
    file.path(fname.movregs, "motion_demean_proactive_run2.1D")
  )

  was.copied[dir.i] <- all(file.copy(fname.movregs, name.dir.i, overwrite = TRUE))
  
}

sum(!was.copied)
sum(was.copied)  ## full proactive set
subjs.pro <- names(was.copied)[was.copied]

write.csv(subjs.pro, "glms/subjs_pro_runwise.csv")




## downsampling items ----

for (dir.i in seq_along(dirs.input)) {
  # dir.i = 1
  
  name.dir.i <- dirs.input[dir.i]
  
  fnames.i <- list.files(name.dir.i)
  
  items.congruent <- c("redRED", "blueBLUE", "whiteWHITE", "purplePURPLE")
  items.congruent1 <- c("blueBLUE", "purplePURPLE")
  items.incongruent1 <- c("blueRED", "blueWHITE", "purpleBLUE", "redPURPLE", "whitePURPLE", "whiteRED")
  # items.incongruent2 <- c("blueRED", "blueWHITE", "purpleBLUE", "redPURPLE", "whitePURPLE", "whiteRED")
  
  fnames.i <- fnames.i[grep(paste(bias.items, collapse = "|"), fnames.i)]
  fnames.i <- fnames.i[grep("acc-only\\.txt$", fnames.i)]  ## discard run-wise times
  
  if (length(fnames.i) != 16) stop("bad length")
  
  modeled.out.congr <- vector("list", length = 2)
  modeled.out.incon <- vector("list", length = 2)
  
  for (stimtime.i in seq_along(fnames.i)) {
    # stimtime.i <- 1
    
    name.stimtime.i <- fnames.i[stimtime.i]
    
    fname.i <- file.path(name.dir.i, name.stimtime.i)
    
    # run1 <- readLines(fname.i, 1)
    # run2 <- readLines(fname.i, 2)
    
    stimtimes <- readChar(fname.i, file.info(fname.i)$size)
    stimtimes <- strsplit(stimtimes, split = "\n")[[1]]
    
    stimtimes1 <- as.numeric(strsplit(stimtimes[1], " ")[[1]])
    stimtimes2 <- as.numeric(strsplit(stimtimes[2], " ")[[1]])
    
    is.run1.item <- grepl(paste0(c(items.congruent1, items.incongruent1), collapse = "|"), name.stimtime.i)
    
    if (is.run1.item) {
      
      keep <- sample.int(length(stimtimes1), min(c(3, length(stimtimes1))))
      stimtimes.keep <- stimtimes1[keep]
      stimtimes.burn1 <- stimtimes1[-keep]
      stimtimes.burn2 <- stimtimes2
      onsets.keep <- paste0(onsets4afni(sort(stimtimes.keep)), "*\n")  ## star in second row (run)
      
    } else {
      
      keep <- sample.int(length(stimtimes2), min(c(3, length(stimtimes2))))
      stimtimes.keep <- stimtimes2[keep]
      stimtimes.burn2 <- stimtimes2[-keep]
      stimtimes.burn1 <- stimtimes1
      onsets.keep <- paste0("*\n", onsets4afni(sort(stimtimes.keep)))  ## star in first row (run)
      
    }
    
    onsets2file(onsets.keep, paste0(gsub(".txt", "", fname.i), "_downsamp"))
      
    
    is.congr <- grepl(paste0(items.congruent, collapse = "|"), name.stimtime.i)
    
    if (is.congr) {
      
      modeled.out.congr[[1]] <- c(stimtimes.burn1, modeled.out.congr[[1]])  ## run 1
      modeled.out.congr[[2]] <- c(stimtimes.burn2, modeled.out.congr[[2]])  ## run 2
      
    } else {
      
      modeled.out.incon[[1]] <- c(stimtimes.burn1, modeled.out.incon[[1]])  ## run 1
      modeled.out.incon[[2]] <- c(stimtimes.burn2, modeled.out.incon[[2]])  ## run 2
      
    }
    
    rm(stimtimes.burn1, stimtimes.burn2, onsets.keep)  ## reset, just in case
    
  }
  
  burn.congr <- c(
    onsets4afni(sort(modeled.out.congr[[1]])),
    onsets4afni(sort(modeled.out.congr[[2]]))
  )
  
  burn.incon <- c(
    onsets4afni(sort(modeled.out.incon[[1]])),
    onsets4afni(sort(modeled.out.incon[[2]]))
  )
  
  onsets2file(burn.congr, gsub("(.*_pro_).*$", "\\1congr_modelout_acc-only", fname.i))
  onsets2file(burn.incon, gsub("(.*_pro_).*$", "\\1incon_modelout_acc-only", fname.i))
  
}

