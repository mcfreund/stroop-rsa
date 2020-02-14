## about ----
## 
## estimates subjects' stroop effects as BLUPs and writes to file.

## setup ----

library(here)
library(mikeutils)
library(dplyr)
library(data.table)
library(nlme)

source(here("code", "strings.R"))


stroop <- read.csv(here("data", "behavior-and-events_group201902.csv"),  stringsAsFactors = FALSE)

stroop %<>%
  arrange(subj, session, run, trial.num) %>%
  group_by(subj, session, run) %>%
  mutate(
    trial.type    = ifelse(trial.type == "i", "incon", "congr"),
    error         = 1 - acc,
    post.error    = lag(error),
    n1.trial.type = lag(trial.type),
    n2.trial.type = lag(n1.trial.type)
  )

stroop.pro.rt <- stroop %>% 
  filter(
    session == "pro", is.analysis.group,
    acc == 1, !is.na(rt), rt < 3000, rt > 250
    )

stroop.pro.rt <- stroop.pro.rt[!(stroop.pro.rt$subj %in% c("849971", "161832") & stroop.pro.rt$rt < 500), ]

## model ----

fit1.het <- lme(
  rt ~ trial.type, random = ~ trial.type | subj,
  data = stroop.pro.rt,
  weights = varIdent(form = ~ 1 | subj),
  control = lmeControl(maxIter = 1e5, msMaxIter = 1e5, niterEM = 1e5, msMaxEval = 1e5),
  method = "REML"
)

## extract

blups <- as.data.frame(coef(fit1.het))
blups %<>% rename(congr = "(Intercept)", stroop = "trial.typeincon") %>% tibble::rownames_to_column("subj")

## write ----

write.csv(blups, here("data", "stroop_blups_rt_group201902.csv"),  row.names = FALSE)
