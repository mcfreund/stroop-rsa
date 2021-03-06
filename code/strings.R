## about ----
## contains strings that are commonly used by scripts in this project.
## 
## these strings are mostly names of factor levels.
## but, some are pointers (paths) to external directories, such as the directory that contains the BOLD timecourse images.
## i try to limit the number of these external pointers.
##
## mike freund, 2019-02-20
## adapted for new project directory 2019-12-24


## paths: pointers for atlases and for BOLD timeseries ----

nodename <- Sys.info()["nodename"]

dir.nil.dmcc2.afni <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/AFNI_ANALYSIS"

if (nodename == "ccplinux1") {
  
  dir.atlas <- "/data/nil-external/ccp/freund/atlases"
  
} else if (nodename == "CCP-FREUND") {
  ## mike freund's (i.e., ccp's) thinkpad
  ## reliant on box drive
  ## assumes box drive location at ./Users/mcf/Box
  
  dir.atlas <- "C:/local/atlases"

}


## factor levels ----

bias.colors <- c("blue", "purple", "red", "white")
bias.words  <- toupper(bias.colors)
bias.items  <- mikeutils::combo.paste(bias.colors, bias.words, sep = "")
bias.items.incon <- bias.items[!bias.items %in% c("blueBLUE", "redRED", "whiteWHITE", "purplePURPLE")]
bias.items.con <- bias.items[!bias.items %in% bias.items.incon]
pc50.colors <- c("black", "green", "pink", "yellow")
pc50.words  <- toupper(pc50.colors)
pc50.items  <- mikeutils::combo.paste(pc50.colors, pc50.words, sep = "")
pc50.items.incon <- pc50.items[!pc50.items %in% c("blackBLACK", "greenGREEN", "pinkPINK", "yellowYELLOW")]
pc50.items.con <- pc50.items[!pc50.items %in% pc50.items.incon]
subjs.analysis <- 
  unique(data.table::fread(here::here("in", "behavior-and-events_group201902.csv"))[is.analysis.group == TRUE]$subj)
subjs.validation <- 
  unique(data.table::fread(here::here("in", "behavior-and-events_group201902.csv"))[is.analysis.group == FALSE]$subj)
subjs.validation <- subjs.validation[!subjs.validation %in% c("DMCC4260551")]  ## coil error!
subj_sum <- data.table::fread(here::here("in", "summary_group201902.csv"))[session == "pro"]
twinpairs <- unique(subj_sum[subj %in% subjs.validation]$twin.pair)
subjs.analysis.red <- unique(subj_sum[is.analysis.group == TRUE & !twin.pair %in% twinpairs]$subj)


# test1 <- unique(subj_sum[is.analysis.group == TRUE & twin.pair %in% twinpairs]$subj)
# test2 <- unique(subj_sum[is.analysis.group == FALSE]$subj)

# identical(subjs.validation, test2)
# setdiff(subjs.analysis.red, subjs.analysis)


sets.of.rois <- c("masks", "mmp")
