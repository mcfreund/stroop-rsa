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






split.movregs <- function(
  to.split,
  dir.to   = dir.analysis,
  dir.from  = "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS/"
) {
  
  ## reads movregs files from nil-bluearc, splits by run
  ## NB: THIS FUNCTION ASSUMES THE MOVREGS FILE HAS nrow == TR
  ## NB: ONLY WORKS FOR MB4 PEOPLE!
  ## 
  ## arg to.split should contain a data.frame of all subject*task*session files to split and copy.
  ## e.g., 
  # to.split <- expand.grid(
  #   subj = "132017",
  #   task = c("Axcpt", "Cuedts", "Stern", "Stroop"),
  #   session = c("baseline", "proactive", "reactive")
  #   )
  # dir.from = dl$nil.dmcc2.hcp.afni
  # dir.to   = file.path(dl$nil.external.freund, "AFNI_ANALYSIS_SUBSUBJECT")
  
  
  dir.from   <- file.path(dir.from, to.split$subj, "INPUT_DATA", to.split$task, to.split$session)
  # dir.to     <- file.path(dir.to, to.split$subj, "INPUT_DATA", to.split$task, to.split$session)
  # filenames  <- paste0("Movement_Regressors_", to.split$session, ".1D")
  # files.from <- file.path(dir.from, filenames)
  # files.to   <- file.path(dir.to, filenames)
  file1 <- file.path(
    dir.from, 
    paste0(
      "Movement_Regressors_", 
      to.split$task, 
      toupper(substr(to.split$session, 1, 1)), substr(to.split$session, 2, 3),
      "1_AP.txt"
    )
  )
  file2 <- file.path(
    dir.from, 
    paste0(
      "Movement_Regressors_", 
      to.split$task, 
      toupper(substr(to.split$session, 1, 1)), substr(to.split$session, 2, 3),
      "2_PA.txt"
    )
  )
  file1.new <- file.path(
    dir.to, 
    paste0(
      "Movement_Regressors_", 
      to.split$task, 
      toupper(substr(to.split$session, 1, 1)), substr(to.split$session, 2, 3),
      "1_AP.1D"
    )
  )
  file2.new <- file.path(
    dir.to, 
    paste0(
      "Movement_Regressors_", 
      to.split$task, 
      toupper(substr(to.split$session, 1, 1)), substr(to.split$session, 2, 3),
      "2_PA.1D"
    )
  )
  
  ## create dirs if they do not exist
  for (file.i in seq_along(dir.to)) if (!dir.exists(dir.to[file.i])) dir.create(dir.to[file.i], recursive = TRUE)
  is.missing.dir <- !(dir.exists(dir.from) & dir.exists(dir.to))
  
  ## get number of TRs (for file splitting)
  trs <- expand.grid(
    task    = c("Axcpt", "Cuedts", "Stern", "Stroop"),
    session = c("proactive", "reactive", "baseline")
  )
  trs <- trs[with(trs, order(task, session)), ]
  trs$num <- c(1220, 1220, 1220, 1300, 1300, 1300, 1200, 1200, 1200, 1080, 1180, 1080)  ## order corresponds!
  
  ## read files
  has.unexpected.nrow <- vector("logical", nrow(to.split))
  for (file.i in seq_len(nrow(to.split))) {
    # file.i = 62
    if (is.missing.dir[file.i]) next
    
    session.i <- to.split$session[file.i] # gsub(".*(baseline).*|.*(proactive).*|.*(reactive).*", "\\1\\2\\3", to.split$session[file.i])
    task.i    <- to.split$task[file.i] # gsub(".*(Axcpt).*|.*(Cuedts).*|.*(Stern).*|.*(Stroop).*", "\\1\\2\\3\\4", files.from[file.i])
    num.trs.i <- trs$num[trs$task == task.i & trs$session == session.i]
    
    movregs.i.run1 <- data.table::fread(
      file1[file.i], 
      sep = "\t", header = FALSE, colClasses = rep("numeric", 6), data.table = FALSE
    )
    movregs.i.run2 <- data.table::fread(
      file2[file.i], 
      sep = "\t", header = FALSE, colClasses = rep("numeric", 6), data.table = FALSE
    )
    
    mask.i <- data.table::fread(
      file.path(dir.from, "movregs_FD.txt")[file.i],
      header = FALSE, data.table = FALSE
    )[[1]]  ## all FD_mask files have same name
    
    if (any(is.na(mask.i)) | any(!complete.cases(movregs.i.run1)) | any(!complete.cases(movregs.i.run2))) stop("NAs")
    
    is.unexpected.nrow <- (nrow(movregs.i.run1) + nrow(movregs.i.run2)) != num.trs.i | length(mask.i) != num.trs.i
    if (is.unexpected.nrow) {
      has.unexpected.nrow[file.i] <- TRUE
      next
    }
    
    ## write
    
    movregs.i.run1 %>% data.table::fwrite(file1.new[file.i], sep = " ", col.names = FALSE)
    movregs.i.run2 %>% data.table::fwrite(file2.new[file.i], sep = " ", col.names = FALSE)
    
    mask.i <- as.numeric(mask.i < 0.9)  ## FD threshold
    mask.i.run1 <- mask.i[seq(1, num.trs.i / 2)]
    mask.i.run2 <- mask.i[seq(num.trs.i / 2 + 1, num.trs.i)]
    
    mask.i.run1 %>% 
      data.frame %>% 
      data.table::fwrite(file.path(dir.to[file.i], "movregs_FD_mask_run1.txt"), sep = " ", col.names = FALSE)
    mask.i.run2 %>% 
      data.frame %>% 
      data.table::fwrite(file.path(dir.to[file.i], "movregs_FD_mask_run2.txt"), sep = " ", col.names = FALSE)
    
  }
  
  data.frame(to.split, is.missing.dir, has.unexpected.nrow)
  
}





## split all movregs, FDmask, and transient/sustained  files in dir

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
  # dir.i = 47
  
  name.dir.i <- dirs.input[dir.i]
  
  fnames.i <- list.files(name.dir.i)
  fnames.stimtimes.i <- fnames.i[grep(paste0("^", subjs[dir.i]), fnames.i)]  ## discard movregs
  fnames.stimtimes.i <- fnames.i[grep("sustained\\.txt|transient\\.txt$", fnames.i)]  ## discard run-wise times
  if (length(fnames.stimtimes.i) != 2) stop("bad length")
  
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
  
  if (subjs[dir.i] == "DMCC4854075") {
    
    res <- split.movregs(
      to.split = data.frame(subj = "562345", task = "Stroop", session = "proactive", stringsAsFactors = FALSE),
      dir.to = name.dir.i
    )
    
  } else {
    
    res <- split.movregs(
      to.split = data.frame(subj = subjs[dir.i], task = "Stroop", session = "proactive", stringsAsFactors = FALSE),
      dir.to = name.dir.i
    )
    
  }
  
  was.copied[dir.i] <- all(!res$is.missing.dir, !res$has.unexpected.nrow)
  
}

sum(!was.copied)
sum(was.copied)  ## full proactive set
subjs.pro <- names(was.copied)[was.copied]



write.csv(subjs.pro, here("glms", "subjs_pro_runwise.csv"))


