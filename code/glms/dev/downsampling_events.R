stop("don't source me")

## setup ----

library(here)
# library(mikeutils)
library(magrittr)
library(dplyr)
library(data.table)
library(purrr)

source(here("code", "strings.R"))


stroop <- fread(here("in", "behavior-and-events_group201902.csv"))
sample.analysis <- unique(filter(stroop, is.analysis.group)$subj)
dir.to.write.in <- here("glms")

subjs <- list.dirs(dir.to.write.in, recursive = FALSE, full.names = FALSE)
subjs <- subjs[subjs != "results"]

dirs.input <- file.path(dir.to.write.in, subjs, "input", "pro")
dirs.input <- dirs.input[dir.exists(dirs.input)]

bias.items.run1 <- c("blueBLUE", "bluePURPLE", "blueRED", "purplePURPLE", "purpleRED", "purpleWHITE", "redWHITE", "whiteBLUE")
bias.items.run2 <- c("blueWHITE", "purpleBLUE", "redBLUE", "redPURPLE", "redRED", "whitePURPLE", "whiteRED", "whiteWHITE")

for (dir.i in seq_along(dirs.input)) {
  # dir.i = 1
  
  name.dir.i <- dirs.input[dir.i]
  
  fnames.i <- list.files(name.dir.i)
  
  fnames.i
  
  
  
  
  
  items.congruent <- c("redRED", "blueBLUE", "whiteWHITE", "purplePURPLE")
  items.congruent1 <- c("blueBLUE", "purplePURPLE")
  # items.congruent2 <- c("blueBLUE", "purplePURPLE")
  
  fnames.i <- fnames.i[grep(paste(items.congruent, collapse = "|"), fnames.i)]
  fnames.i <- fnames.i[grep("acc-only\\.txt$", fnames.i)]  ## discard run-wise times
  
  if (length(fnames.i) != 4) stop("bad length")
  
  modeled.out <- vector("list", length = 2)
  
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
    
    is.run1.item <- grepl(paste0(items.congruent1, collapse = "|"), name.stimtime.i)
    
    if (is.run1.item) {
      
      keep <- sample.int(length(stimtimes1), 3)
      stimtimes.keep <- stimtimes1[keep]
      stimtimes.burn1 <- stimtimes1[-keep]
      stimtimes.burn2 <- stimtimes2
      onsets.keep <- paste0(onsets4afni(sort(stimtimes.keep)), "*\n")  ## star in second row (run)
      
    } else {
      
      keep <- sample.int(length(stimtimes2), 3)
      stimtimes.keep <- stimtimes2[keep]
      stimtimes.burn2 <- stimtimes2[-keep]
      stimtimes.burn1 <- stimtimes1
      onsets.keep <- paste0("*\n", onsets4afni(sort(stimtimes.keep)))  ## star in first row (run)
      
    }
    
    
    modeled.out[[1]] <- c(stimtimes.burn1, modeled.out[[1]])  ## run 1
    modeled.out[[2]] <- c(stimtimes.burn2, modeled.out[[2]])  ## run 2
    
    onsets2file(onsets.keep, paste0(gsub(".txt", "", fname.i), "_downsamp"))
    
  }
  
  burn <- c(
    onsets4afni(sort(modeled.out[[1]])),
    onsets4afni(sort(modeled.out[[2]]))
  )
  
  
  
  onsets2file(burn, gsub("(.*_bas_).*$", "\\1congr_modelout_acc-only", fname.i))
  
}

