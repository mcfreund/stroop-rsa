library(lme4)
library(dplyr)
library(mikeutils)
library(ggplot2)

## what does heirarchical estimation and shrinkage buy us in terms of reducing expected test error?

## read data ----

stroop <- read.csv(here::here("data", "behavior-and-events_group201902.csv"))
stroop.pro <- stroop %>% 
  filter(session == "pro", rt > 0, !is.na(rt), acc == 1) %>% 
  mutate(run = as.factor(paste0("run", run))) %>%
  select(subj, run, rt, trial.type)


## fit models ----


## ols

ols <- stroop.pro %>%
  group_by(subj, trial.type, run) %>%
  summarize(rt = mean(log(rt))) %>%
  tidyr::spread(trial.type, rt) %>%
  transmute(stroop = i - c, run = run) %>%
  tidyr::spread(run, stroop) %>%
  as.data.frame
cor(ols[c("run1", "run2")])


## full

fit <- lmer(
  log(rt) ~ trial.type * run + (trial.type * run | subj),
  stroop.pro,
  control = lmerControl(optim = "bobyqa", optCtrl = list(maxfun = 1E5))
)
summary(fit)
tau <- VarCorr(fit)$subj
m <- rbind(c(0, 1, 0, 0), c(0, 1, 0, 1))
r <- cov2cor(m %*% tau %*% t(m))[1, 2]
u <- as.matrix(coef(fit)$subj) %*% t(m)
cor.u <- cor(u)[1, 2]


## split by run

fit1 <- lmer(
  log(rt) ~ trial.type + (trial.type | subj),
  stroop.pro,
  control = lmerControl(optim = "bobyqa", optCtrl = list(maxfun = 1E5)),
  subset = run == "run1"
)
summary(fit1)
tau1 <- VarCorr(fit1)$subj
u1 <- as.matrix(coef(fit1)$subj)


fit2 <- lmer(
  log(rt) ~ trial.type + (trial.type | subj),
  stroop.pro,
  control = lmerControl(optim = "bobyqa", optCtrl = list(maxfun = 1E5)),
  subset  = run == "run2"
)
summary(fit2)
tau2 <- VarCorr(fit2)$subj
u2 <- as.matrix(coef(fit2)$subj)
cor(u1, u2)


## subsampled two-run model

## if still high correlation, then amt of data is important factor
## if low correlation, then some sort of information exchange between runs?

stroop.pro.subsamp <- stroop.pro %>%
  group_by(subj, trial.type, run) %>%
  sample_n(size = n() / 2)

fit.subsamp <- update(fit, data = stroop.pro.subsamp)
summary(fit.subsamp)
tau.subsamp <- VarCorr(fit.subsamp)$subj
r.subsamp <- cov2cor(m %*% tau.subsamp %*% t(m))[1, 2]
u.subsamp <- as.matrix(coef(fit.subsamp)$subj) %*% t(m)
cor.u.subsamp <- cor(u.subsamp)[1, 2]

