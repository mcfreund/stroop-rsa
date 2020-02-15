## about ----
## 
## estimates subjects' stroop effects as BLUPs and writes to csv.
## additionally saves the lme() object to RDS file in same location as csv.
##


## setup ----

library(here)
library(mikeutils)
library(dplyr)
library(data.table)
library(nlme)

source(here("code", "strings.R"))


farout <- function(x) {
  
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr3 <- IQR(x) * 3
  
  x < (q1 - iqr3) | x > (q3 + iqr3)
  
}

## read and subset

stroop.pro <- read.csv(here("data", "behavior-and-events_group201902.csv"),  stringsAsFactors = FALSE) %>%
  filter(session == "pro", is.analysis.group) %>%
  mutate(trial.type = ifelse(trial.type == "i", "incon", "congr"))


## model ----

## initial fit  

stroop.pro.rt <- stroop.pro %>% filter(acc == 1, !is.na(rt), rt < 3000, rt > 250)
is.weird.rt <- stroop.pro.rt$subj %in% c("849971", "161832") & stroop.pro.rt$rt < 500
stroop.pro.rt <- stroop.pro.rt[!is.weird.rt, ]

fit1.het <- lme(
  rt ~ trial.type, 
  random  = ~ trial.type | subj,
  data    = stroop.pro.rt,
  weights = varIdent(form = ~ 1 | subj),
  control = lmeControl(maxIter = 1e5, msMaxIter = 1e5, niterEM = 1e5, msMaxEval = 1e5),
  method  = "REML"
)

## trim and re-fit

stroop.pro.rt$resid.p <- resid(fit1.het, type = "p")
stroop.pro.rt$is.far.out <- farout(stroop.pro.rt$resid.p)

fit1.het.trim <- update(fit1.het, subset = !is.far.out)

## extract predictions

blups <- as.data.frame(coef(fit1.het.trim))
blups %<>% rename(congr = "(Intercept)", stroop = "trial.typeincon") %>% tibble::rownames_to_column("subj")


## write ----

write.csv(blups, here("out", "behav", "stroop_blups_rt_group201902.csv"),  row.names = FALSE)  ## blups
saveRDS(fit1.het.trim, here("out", "behav", "fit1-het-trim_group201902.RDS"))
