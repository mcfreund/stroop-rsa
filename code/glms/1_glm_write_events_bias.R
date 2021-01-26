stop("don't source me")
## about ----
## contains strings that are commonly used by scripts in this project.
## 
## these strings are mostly names of factor levels.
## some are pointers (paths) to external directories, such as the dir that contains the BOLD timecourse images.
## but, i try to limit the number of these external pointers.
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


stroop <- fread(here("in", "behavior-and-events_group201902.csv"))
subj_sum <- fread(here("in", "summary_group201902.csv"))
stroop.pro <- subset(stroop, session == "pro")
stroop.bas <- subset(stroop, session == "bas")


## funs ----


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

onsets4afni.dm <- function(x) {
  ## analogous to onsets4afni(), but adapted for use with dmBLOCK().
  if (class(x) != "data.frame") stop("needs data.frame!")
  if (any(!complete.cases(x$dm))) stop("incomplete cases in dm!")
  x <- x[complete.cases(x$time.onset), ]
  if (nrow(x) < 1) {
    ## no valid times needs -1 option, e.g.:
    ## https://afni.nimh.nih.gov/afni/community/board/read.php?1,84354,85947#msg-85947
    x <- "-1:0\n"  
  } else {
    x <- apply(x, 1, paste, collapse = ":") %>% paste0(collapse = " ") %>% paste0(" \n")
  }
  return(x)
}

onsets2file <- function(.onsets, .fname) {
  .fname <- paste0(.fname, ".txt")
  .fout <- file(.fname, "wt")
  cat(.onsets, file = .fout)
  close(.fout)
  unlink(.fout)
}


write.events <- function(
  df, 
  dir.analysis, 
  grouping.var.name,
  grouping.var.values,  ## NULL if each subject should have ONLY events that they have stim times for
  onset.var.name,
  reg.suffix
) {
  ## this function should operate on a df that contains data from single subject and
  ## single session. best used when called from lapply on a list of many such
  ## dataframes. requires composite 'grouping.var' as a column.
  ## for DURATION FIXED events.
  df <- as.data.frame(df)
  if (!grouping.var.name %in% names(df)) stop("add grouping.var!")
  df <- df %>% rename(grouping.var = paste0(grouping.var.name))  ## for convenience
  subj <- unique(df$subj)
  session <- unique(df$session)
  if (length(c(subj, session)) != 2) stop("splitting error!")
  dir.input <- file.path(dir.analysis, subj, "input", session)
  if (!dir.exists(dir.input)) dir.create(dir.input, recursive = TRUE)
  ## loop across values
  if (is.null(grouping.var.values)) grouping.var.values <- unique(df$grouping.var) ## programmatic option
  num.onsets <- vector("numeric", length = length(grouping.var.values))
  for (ii in seq_along(grouping.var.values)) {
    # ii <- 1
    value.ii <- grouping.var.values[ii]
    df.ii    <- df[df$grouping.var == value.ii, ]
    onsets1  <- df.ii[df.ii$run == "run1", onset.var.name] %>% sort %>% onsets4afni
    onsets2  <- df.ii[df.ii$run == "run2", onset.var.name] %>% sort %>% onsets4afni
    onsets   <- paste0(onsets1, onsets2, collapse = "")
    onsets2file(
      onsets, .fname = file.path(dir.input, paste0(subj, "_", session, "_", value.ii, "_", reg.suffix))
    )
    num.onsets[ii] <- nrow(df.ii)  ## book-keeping
  }
  noquote(num.onsets)
}


write.blocks <- function(
  df,
  dir.analysis
) {
  block.events <- df %>%
    select(subj, session, run, contains("time.block")) %>%
    filter(!duplicated(.))
  num.extra <- block.events %>% 
    select(-contains("time.block")) %>%
    table %>% as.data.frame %>% filter(Freq > 1) %>% nrow
  if (num.extra > 0) stop("multiple block times for single session * subj")
  block.labels <- unique(block.events[, c("subj", "session")]) %>% as.data.frame
  for (ii in seq_len(nrow(block.labels))) {
    # ii <- 1
    subj.ii <- block.labels[ii, "subj"]
    session.ii <- block.labels[ii, "session"]
    block.events.ii <- block.events %>% filter(subj == subj.ii, session == session.ii)
    if (nrow(block.events.ii) > 2) stop("should've been caught earlier (weird...)")
    dir.input <- file.path(dir.analysis, subj.ii, "input", session.ii)
    block.events.ii %>%
      reshape2::melt(id = c("subj", "session", "run")) %>%
      mutate(
        block = gsub(".*([1-3]).*", "\\1", variable),  ## pull out block
        variable = gsub("time.block[1-3].", "", variable)  ## pull out on / off
      ) %>%
      tidyr::spread(variable, value) %>%
      mutate(
        dm = off - on,
        ## for blocks with no off value (e.g. truncated run), replace duration with mean:
        dm = ifelse(is.na(off) & !is.na(on), mean(dm, na.rm = TRUE), dm)
      ) %>%
      rename(time.onset = on, time.offset = off) -> block.events.ii
    sustained1 <- block.events.ii %>% 
      filter(run == "run1", complete.cases(time.onset)) %>% arrange(time.onset) %>%
      select(time.onset, dm)
    sustained2 <- block.events.ii %>% 
      filter(run == "run2", complete.cases(time.onset)) %>% arrange(time.onset) %>% 
      select(time.onset, dm)
    sustained <- paste0(onsets4afni.dm(sustained1), onsets4afni.dm(sustained2), collapse = "")
    transient1 <- block.events.ii %>% filter(run == "run1") %>% select(time.onset, time.offset) %>% unlist
    transient2 <- block.events.ii %>% filter(run == "run2") %>% select(time.onset, time.offset) %>% unlist
    transient <- paste0(
      onsets4afni(transient1 %>% sort), onsets4afni(transient2 %>% sort),
      collapse = ""
    )
    sustained %>% onsets2file(file.path(dir.input, paste0(subj.ii, "_", session.ii, "_sustained")))
    transient %>% onsets2file(file.path(dir.input, paste0(subj.ii, "_", session.ii, "_transient")))
  }
}


copy.movregs <- function(
  .dir.nil = dir.nil.dmcc2.afni,
  .dir.analysis, 
  .overwrite = FALSE
) {
  patterns <- c("movregs_FD_mask", "motion_demean")
  subjs.in.dir <- list.files(.dir.analysis) %>% .[grep("^[0-9]|*^DMCC", .)]
  # subjs.in.nil <- list.files(.dir.nil) %>% .[grep("^[0-9]|*^DMCC", .)]
  pattern4grep <- paste0("*", patterns, "*") %>% paste0(collapse = "|")
  pattern4gsub <- paste0(patterns, ".*$") %>% paste0(collapse = "|")
  files.source <- list.files(.dir.nil, full.names = TRUE) %>% 
    ## get subj folders in /scratch1/* that have folders in /data/nil-bluearc/*:
    .[grep(paste0(subjs.in.dir, collapse = "|"), .)] %>%
    paste0("/INPUT_DATA/Stroop") %>%
    ## NB, 2018-10-20: line following this one commented; no mb8 subjs in sample,
    ## giving error:
    list.files(recursive = TRUE, full.names = TRUE, pattern = pattern4grep) #%>%
  # .[-grep("*movregs_24*", .)]  ## rm these files (from MB 8 subjs, it seems)
  folders.dest <- gsub(".*AFNI_ANALYSIS/", "", files.source) %>%
    gsub(pattern4gsub, "", .) %>%
    gsub("INPUT_DATA/Stroop", "input", .) %>%
    paste0(.dir.analysis, "/", .)
  folders.dest <- gsub("baseline", "bas", folders.dest)
  folders.dest <- gsub("proactive", "pro", folders.dest)
  folders.dest <- gsub("reactive", "rea", folders.dest)
  ## build df with file info to loop over
  df <- data.frame(
    subj    = gsub(".*ANALYSIS/([[:alnum:]]{6,11})/INPUT.*", "\\1", files.source),
    session = gsub(".*(bas)eline.*|.*(pro)active.*|.*(rea)ctive.*",  "\\1\\2\\3", files.source),
    source  = files.source,
    dest    = folders.dest,
    stringsAsFactors = FALSE
  )
  df$was.copied <- FALSE
  for (ii in seq_along(unique(df$dest))) {
    # ii <- 44
    # print(ii)
    this.dest <- unique(df$dest)[ii]
    if (!dir.exists(this.dest)) next
    df[df$dest == this.dest, "was.copied"] <- file.copy(
      df[df$dest == this.dest, "source"], this.dest, overwrite = .overwrite
    )
  }
  return(df)
}








## use (proactive) ----

## make grouping.var:
is.nuisance <- stroop.pro$acc.final %in% c("0", "no.response", "unintelligible")
stroop.pro <- stroop.pro %>%
  ungroup %>%
  mutate(
    reg = item,
    reg = ifelse(pc == "pc50", paste0("pc50_", trial.type), reg),
    reg = ifelse(is.nuisance, "nuisance", reg)
  )
unique(stroop.pro$reg)
sum(is.na(stroop.pro$reg))

## NB: write events checks for "run1" and "run2" text in column "run"!!!
if (any(unique(stroop.pro$run) %in% 1:2)) stroop.pro$run <- ifelse(stroop.pro$run == 1, "run1", "run2")

dir.to.write.in <- here("glms")

## first, check:
stroop.pro %>% split(list(.$subj, .$session)) %>% map_dbl(nrow) %>% unique  ## should be length 1, of value 216
head(stroop.pro[, grep("time.block", names(stroop.pro))])
stroop.pro[, grep("time.block", names(stroop.pro))] %>% range

## write events
num.events.written <- stroop.pro %>% 
  split(list(.$subj, .$session)) %>%
  map(
    write.events, 
    grouping.var.name   = "reg",
    grouping.var.values = unique(stroop.pro$reg),
    dir.analysis        = dir.to.write.in,
    onset.var.name      = "time.target.onset",
    reg.suffix          = "acc-only"
  )

## write blocks
stroop.pro %>% write.blocks(dir.analysis = dir.to.write.in)  ## runs silently

## copy movregs
summary.movregs <- copy.movregs(.dir.analysis = dir.to.write.in)

## write summaries ----

## event files
df.num.events.written <- matrix(unlist(num.events.written), nrow = length(unique(stroop.pro$subj)), byrow = TRUE)
df.num.events.written <- data.frame(unique(stroop.pro$subj), df.num.events.written)
names(df.num.events.written) <- c("subj", unique(stroop.pro$reg))
write.csv(df.num.events.written, here("out", "summaries", "event-files_group201902.csv"))

## movregs
summary.movregs <- summary.movregs %>% filter(session == "pro")
write.csv(summary.movregs, here("out", "summaries", "moveregs_group201902.csv"))










## use (baseline) ----

## make grouping.var:
is.nuisance <- stroop.bas$acc.final %in% c("0", "no.response", "unintelligible")
stroop.bas <- stroop.bas %>%
  ungroup %>%
  mutate(
    reg = item,
    reg = ifelse(pc == "pc50", paste0("pc50_", trial.type), reg),
    reg = ifelse(is.nuisance, "nuisance", reg)
  )
unique(stroop.bas$reg)
sum(is.na(stroop.bas$reg))

## NB: write events checks for "run1" and "run2" text in column "run"!!!
if (any(unique(stroop.bas$run) %in% 1:2)) stroop.bas$run <- ifelse(stroop.bas$run == 1, "run1", "run2")

dir.to.write.in <- here("glms")

## first, check:
## should be length 1, of value 216:

stroop.bas.subj.nums <- stroop.bas %>% split(list(.$subj, .$session)) %>% map_dbl(nrow)
stroop.bas.subjs <- names(stroop.bas.subj.nums)[stroop.bas.subj.nums == 216] %>% gsub(".bas", "", .)
stroop.bas %<>% filter(subj %in% stroop.bas.subjs)

stroop.bas %>% split(list(.$subj, .$session)) %>% map_dbl(nrow) %>% unique  ## should be length 1, of value 216
head(stroop.bas[, grep("time.block", names(stroop.bas))])
stroop.bas[, grep("time.block", names(stroop.bas))] %>% range

## write events

num.events.written <- stroop.bas %>% 
  split(list(.$subj, .$session)) %>%
  map(
    write.events, 
    grouping.var.name   = "reg",
    grouping.var.values = unique(stroop.bas$reg),
    dir.analysis        = dir.to.write.in,
    onset.var.name      = "time.target.onset",
    reg.suffix          = "acc-only"
  )

## write blocks
stroop.bas %>% write.blocks(dir.analysis = dir.to.write.in)  ## runs silently

## copy movregs
summary.movregs <- copy.movregs(.dir.analysis = dir.to.write.in)


## write summaries ----

## event files
df.num.events.written <- matrix(unlist(num.events.written), nrow = length(stroop.bas.subjs), byrow = TRUE)
df.num.events.written <- data.frame(stroop.bas.subjs, df.num.events.written)
names(df.num.events.written) <- c("subj", unique(stroop.bas$reg))
write.csv(df.num.events.written, here("out", "summaries", "event-files_bas_group201902.csv"))

## movregs
summary.movregs <- summary.movregs %>% filter(session == "bas")
write.csv(summary.movregs, here("out", "summaries", "moveregs_bas_group201902.csv"))



## downsampling congruent items ----

subjs <- list.dirs(dir.to.write.in, recursive = FALSE, full.names = FALSE)

dirs.input <- file.path(dir.to.write.in, subjs, "input", "bas")
dirs.input <- dirs.input[dir.exists(dirs.input)]


for (dir.i in seq_along(dirs.input)) {
  # dir.i = 1
  
  name.dir.i <- dirs.input[dir.i]
  
  fnames.i <- list.files(name.dir.i)
  
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
