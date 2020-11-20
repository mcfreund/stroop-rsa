## define subject samples ----


## full sample ----


## subject list for pro_bias glm ----
## mike freund, 2019-02-19

library(dplyr)
library(data.table)
library(purrr)

## read sheets ----

subj_sum <- read.csv(
  "C:/Users/mcf/Box/global/proj/cross-task-rsa/sheets/dmcc2_all-runs-summary.csv",
  # "C:/Users/mcf/Box/global/wustl/proj/cross-task-rsa/sheets/dmcc2_all-runs-summary.csv",
  stringsAsFactors = FALSE
)
stroop <- read.csv(
  "C:/Users/mcf/Box/global/proj/cross-task-rsa/sheets/dmcc2_behavior-and-events-stroop.csv",
  # "C:/Users/mcf/Box/global/wustl/proj/cross-task-rsa/sheets/dmcc2_behavior-and-events-stroop.csv",
  stringsAsFactors = FALSE
)
stroop.rt <- stroop %>% filter(acc == 1, !is.na(rt))
stroop.er <- stroop %>% filter(!is.na(acc)) %>% mutate(error = 1 - acc)


sample.pro.str <- subj_sum %>% filter(session == "pro", task == "str", mb == "four")
(missing.nii <- unique(filter(sample.pro.str, num.nii < 1)$subj))  ## 6 subjs missing images
sample.pro.str <- sample.pro.str %>% filter(!subj %in% missing.nii)
missing.acc <- unique(filter(sample.pro.str, missing.rows)$subj)  ## 5 subj missing accuracy coding
# missing.acc <- unique(filter(sample.pro.str, per.missing > 0)$subj)  ## 5 subj missing accuracy coding
sample.pro.str <- sample.pro.str %>% filter(!subj %in% missing.acc)

## look at their behavioral data:
data.sample.pro.str <- stroop %>% filter(subj %in% unique(sample.pro.str$subj), session == "pro")
sum(is.na(data.sample.pro.str$acc.final)) > 0   ## to check if any NAs---would be case if subj missing data
updated.sample <- unique(sample.pro.str$subj)

## these subjs should be good to go to fit fMRI GLMs.
## but, first compare to previous analysis group:
prior.sample <- read.csv(
  "C:/Users/mcf/Box/global/proj/stroop-rsa/old/wustl_proj_stroop-rsa/sheets/subj-sample_group-201812.csv",
  # "C:/Users/mcf/Box/global/wustl/proj/stroop-rsa/sheets/subj-sample_group_201812.csv",
  stringsAsFactors = FALSE, header = TRUE
)[[1]]

setdiff(prior.sample, updated.sample)  ## subjs NOT in new sample that were in old (should be character(0)---and is!)
(new.subjs <- setdiff(updated.sample, prior.sample))  ## subjs ADDED
length(new.subjs)  ## TOTAL num subjs added: 23 (including twins, and sleepy subjs)
sample.pro.str %>% 
  filter(subj %in% new.subjs, run == 1, !is.na(twin.pair)) %>% 
  .$twin.pair %>% duplicated %>% sum  ## num twins in ADDED subjs---suprisingly only 3?!
total.num.paired <- sample.pro.str %>% 
  filter(run == 1, !is.na(twin.pair)) %>% 
  .$twin.pair %>% duplicated %>% sum  ## TOTAL num twins in updated sample: 17
total.num.cotwins <- sample.pro.str %>% 
  filter(run == 1, !is.na(twin.pair)) %>% 
  nrow  ## TOTAL num twins in updated sample: 17
total.num.unrel <- sample.pro.str %>% filter(run == 1, is.na(twin.pair)) %>% nrow
total.num.cotwins - total.num.paired + total.num.unrel  ## final N: 49 subjects.

## write 'final' subject list to file ----
write.csv(
  updated.sample,
  "C:/Users/mcf/Box/global/wustl/proj/stroop-rsa/sheets/subj-sample_group_201902.csv",
  row.names = FALSE
)


## analysis and validation sets ----


## about ---- 
## generating subject list for group201902 ANALYSES 


## set up env ----

library(here)
library(dplyr)
source(here("..", "gen", "funs", "_get_dirs_local.R"))
source(here("..", "gen", "funs", "_funs.R"))
source(here("r", "group-201902", "_get_misc_vars.R"))

## read subject summary / overview sheet:
subj.summary <- read.csv(
  "C:/Users/mcf/Box/global/proj/cross-task-rsa/sheets/dmcc2_all-runs-summary.csv",
  # file.path(dir.box.crosstask.sheets, "dmcc2_all-runs-summary.csv"),
  stringsAsFactors = FALSE
)
subj.summary.pro <- filter(subj.summary, session == "pro", task == "str")
stroop <- read.csv(
  "C:/Users/mcf/Box/global/proj/cross-task-rsa/sheets/dmcc2_behavior-and-events-stroop.csv",
  # file.path(dir.box.crosstask.sheets, "dmcc2_behavior-and-events-stroop.csv"),
  stringsAsFactors = FALSE
)
stroop.pro <- filter(stroop, session == "pro")

## read previous subject list:
## ( a totally unrelated group of subjs)
group201812 <- read.csv(
  "C:/Users/mcf/Box/global/proj/stroop-rsa/old/wustl_proj_stroop-rsa/sheets/subj-sample_group-201812.csv",
  # file.path(dir.box.stroop.sheets, "subj-sample_group-201812.csv"),
  stringsAsFactors = FALSE
)[, 1]
# subj.summary.pro %>% filter(subj %in% group201812, run == 1, !is.na(twin.pair)) %>% .$twin.pair %>% duplicated %>% sum
## and current subject list---i.e., all pepople with GLMs fit:
## (contains twins!)
group201902 <- read.csv(
  "C:/Users/mcf/Box/global/proj/stroop-rsa/old/wustl_proj_stroop-rsa/sheets/subj-sample_group-201902.csv",
  # file.path(dir.box.stroop.sheets, "subj-sample_group-201902.csv"),
  stringsAsFactors = FALSE
)[, 1]
# subj.summary.pro %>% filter(subj %in% group201902, run == 1, !is.na(twin.pair)) %>% .$twin.pair %>% duplicated %>% sum

## compare subjects in each:
setdiff(group201812, group201902)  ## should be 0: no subjs in group201812 that aren't in 201902
setdiff(group201812, subj.summary.pro$subj)  ## should be 0
setdiff(group201902, subj.summary.pro$subj)  ## should be 0
setdiff(subj.summary.pro$subj, group201902)  ## probably no image subjects; mb8 subjects; twins...
# subj.summary.pro %>% filter(subj %in% setdiff(subj.summary.pro$subj, group201902), run == 1)


## remove cotwins in group201902 that weren't in group201812---keep all members of group201812.
pairs.to.remove.from.new <- unique(filter(subj.summary.pro, subj %in% group201812, !is.na(twin.pair))$twin.pair)
group201902.no.cross.dupes <- filter(
  subj.summary.pro, 
  !twin.pair %in% pairs.to.remove.from.new & subj %in% group201902
)$subj %>% unique

## "new" subjects (new in quotes bc not necessarily new---just not included in group201812)
new.subjs <- setdiff(group201902.no.cross.dupes, group201812)
## see C:\Users\mcf\Box\global\wustl\proj\stroop-rsa\knitrs\group-201812\stroop_rsa_behavior.pdf for original
## exclusion criteria.
## must have: mb4, proactive both runs, more than 90% 'useable.'
## known somnolent subj: 233326 (25% error; intradb notes / acc coding notes show)
new.subjs.summary <- filter(subj.summary.pro, subj %in% new.subjs)
new.subjs.summary <- new.subjs.summary %>% filter(subj != "233326")
## eliminate duplicate rows and extraneous cols
new.subjs.list <- new.subjs.summary %>% filter(run == 1)  %>% select(subj, twin.type, twin.pair)

set.seed(32)
## random sample of cotwins:
new.subjs.list.unrl <- bind_rows(
  new.subjs.list %>%
    filter(!is.na(twin.pair)) %>%  group_by(twin.pair) %>% sample_n(1),  ## randomly sample one cotwin from each "pair"
  new.subjs.list %>% filter(is.na(twin.pair))  ## get only non-twins
)

group201902.analysis <- c(group201812, new.subjs.list.unrl$subj)

## sanity checks:
subj.summary.pro %>% filter(subj %in% group201902.analysis, run == 1) %>% View
## all subjs in original analysis in this group?
setdiff(group201812, group201902.analysis)  ## yes
## any twins?
any(duplicated(filter(subj.summary.pro, subj %in% group201902.analysis, run == 1, !is.na(twin.pair))$twin.pair))
length(group201902.analysis)  ## my original number. alright.

## held out set:
held.out <- subj.summary.pro %>% filter(
  !subj %in% c(unique(filter(subj.summary.pro, missing.rows)$subj), group201902.analysis)
)
held.out %>% filter(run == 1, !is.na(twin.pair)) %>% .$twin.pair%>% duplicated %>% sum
held.out %>% filter(run == 1, !is.na(twin.pair)) %>% .$twin.pair %>% sort 
stroop.pro %>% 
  filter(subj %in% c("198855", "233326")) %>%
  group_by(subj, session) %>%
  summarize(
    commission = mean(acc.final == "0", na.rm = TRUE) * 100, 
    omission = mean(acc.final == "no.response", na.rm = TRUE) * 100
  )  ## these people were probably falling asleep
ok.subjs <- unique(filter(held.out, num.nii > 0, per.error < 20 )$subj)
## one twinpair left:
held.out %>% filter(subj %in% ok.subjs, run == 1, !is.na(twin.pair)) %>% .$twin.pair %>% duplicated %>% sum
held.out %>% filter(subj %in% ok.subjs, run == 1, !is.na(twin.pair)) %>% .$twin.pair %>% sort
remove.subj <- sample(unique(filter(held.out, twin.pair == "9")$subj), 1)  ## which twin to remove?)
ok.subjs <- setdiff(ok.subjs, remove.subj)
held.out <- held.out %>% filter(subj %in% ok.subjs)
## should run 5 additional subjects, most of whom are mb8
filter(held.out, subj %in% setdiff(held.out$subj, group201902))


## write:
write.csv(
  group201902.analysis, 
  file.path(dir.box.stroop.sheets, "subj-sample_group-201902_for-analysis.csv"),
  row.names = FALSE
)

write.csv(
  unique(held.out$subj), 
  file.path(dir.box.stroop.sheets, "subj-sample_group-201902_held-out.csv"),
  row.names = FALSE
)

