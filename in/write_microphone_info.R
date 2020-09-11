## read data ----
library(here)
library(dplyr)

stroop <- read.csv(here("in", "behavior-and-events_group201902.csv"),  stringsAsFactors = FALSE)
dates1 <- read.csv(here("in", "dates.csv"), stringsAsFactors = FALSE)
dates2 <- read.csv(here("in", "dates_missing_from-gcal.csv"), stringsAsFactors = FALSE)
dates1[c("BAS.day", "PRO.day", "REA.day")] <- lapply(dates1[c("BAS.day", "PRO.day", "REA.day")], as.Date)
dates2[c("session1", "session2", "session3")] <- lapply(dates2[c("session1", "session2", "session3")], as.Date)

## identify microphone used by dates ----

## fomri arrive (range):
fomri.arrive.lb <- "2018-05-01"
fomri.arrive.ub <- "2018-05-31"

## definite fomri
fomri.switch <- "2018-06-28"
## exceptions:
##  DMCC515268_baseline DMCC3 (June 30th, 2018) and DMCC6904377_reactive DMCC2 (July 1st, 2018)

dates1.is.microoptics <- dates1[c("BAS.day", "PRO.day", "REA.day")] < fomri.arrive.lb
dates1.is.fomri <- dates1[c("BAS.day", "PRO.day", "REA.day")] > fomri.switch

all(dates1.is.microoptics, na.rm = TRUE)  ## all microoptics

dates2.is.microoptics <- dates2[c("session1", "session2", "session3")] < fomri.arrive.lb
dates2.is.fomri <- dates2[c("session1", "session2", "session3")] > fomri.switch

dates2$all.microoptics <- apply(dates2.is.microoptics, 1, function(x) all(x))
dates2$all.fomri <- apply(dates2.is.fomri, 1, function(x) all(x))

subj.microoptics <- dates2 %>% filter(all.microoptics & subj != "DMCC6904377") %>% .$subj
subj.microoptics <- union(subj.microoptics, dates1$subj)
subj.fomri <- dates2 %>% filter(all.fomri & subj != "DMCC6904377") %>% .$subj
subj.unknown <- setdiff(unique(stroop$subj), c(subj.microoptics, subj.fomri))


mic <- ifelse(
  unique(stroop$subj) %in% subj.fomri, 
  "fomri",
  ifelse(
    unique(stroop$subj) %in% subj.microoptics,
    "micro", 
    ifelse(
      unique(stroop$subj) %in% subj.unknown,
      "unknown",
      NA
     )
   )
  )

micinfo <- data.frame(subj = unique(stroop$subj), mic = mic)

data.table::fwrite(micinfo, here("in", "microphone-info_group201902.csv"))
