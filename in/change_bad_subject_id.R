source(here::here("code", "strings.R"))
source(here::here("code", "packages.R"))



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

