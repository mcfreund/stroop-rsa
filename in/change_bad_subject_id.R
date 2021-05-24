source(here::here("code", "strings.R"))
source(here::here("code", "packages.R"))
source(here::here("code", "funs.R"))


## behavior-and-events_group

behav <- fread(here("in", "behavior-and-events_group201902.csv"))
behav[subj == "DMCC4854075"]$subj <- "562345"
fwrite(behav, here("in", "behavior-and-events_group201902.csv"))

## summary_group20

subj_sum <- fread(here("in", "summary_group201902.csv"))
subj_sum[subj == "DMCC4854075"]$subj <- "562345"
fwrite(subj_sum, here("in", "summary_group201902.csv"))

## glm subdirectories

file.rename(here("glms", "DMCC4854075"), here("glms", "562345"))
"DMCC4854075" %in% list.files(here("glms"))
"562345" %in% list.files(here("glms"))


change_dimname <- function(fname) {
  
  R <- readRDS(here::here("out", "rsa", "obsv", fname))
  
  if ("DMCC4854075" %in% dimnames(R)$subj) {
    
    dimnames(R)$subj[dimnames(R)$subj == "DMCC4854075"] <- "562345"
    saveRDS(R, here::here("out", "rsa", "obsv", fname))
    
    dimnames(R)
    
  } else {
    
    return("not in dimnames")
    
  }
  
}

fs <- c(
  "rsarray_pro_bias_acc-only_masks.rds", "rsarray_pro_bias_acc-only_masks_residual-line.rds", "rsarray_pro_bias_acc-only_masks_residual-rank.rds",
  "rsarray_pro_bias_acc-only_mmp.rds", "rsarray_pro_bias_acc-only_mmp_residual-line.rds", "rsarray_pro_bias_acc-only_mmp_residual-rank.rds"
  )

res <- lapply(fs, change_dimname)

# res