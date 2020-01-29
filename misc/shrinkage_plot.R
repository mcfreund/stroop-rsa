library(lme4)
library(dplyr)

## read data ----

stroop <- read.csv(here::here("data", "behavior-and-events_group201902.csv"))
stroop.pro <- stroop %>% 
  filter(session == "pro", rt > 0, !is.na(rt), acc == 1) %>% 
  mutate(run = as.factor(paste0("run", run))) %>%
  select(subj, run, rt, trial.type)

## fit model ----

fit <- lmer(
  log(rt) ~ trial.type * run + (trial.type * run | subj),
  stroop.pro,
  control = lmerControl(optim = "bobyqa", optCtrl = list(maxfun = 1E5))
)
summary(fit)
tau <- VarCorr(fit)$subj
m <- rbind(c(0, 1, 0, 0), c(0, 1, 0, 1))
r <- cov2cor(m %*% tau %*% t(m))[1, 2]
u <- as.matrix(ranef(fit)$subj) %*% t(m)
cor.u <- cor(u)[1, 2]


## bootstrap confint ----

beg.confint <- Sys.time()
ci <- confint(fit, method = "boot")
time.confint <- beg.confint - Sys.time()


## randomization test ----

n.resamp <- 1E5 + 1000
resamples <- replicate(n.resamp, sample(unique(stroop.pro$subj)))
resamples <- unique(resamples, MARGIN = 2)  ## remove duplicate permutations (prob. none)

d1 <- stroop.pro[stroop.pro$run == "run1", ]  ## split data by run
d2 <- stroop.pro[stroop.pro$run == "run2", ]
d2 <- split(d2, d2$subj)  ## split run 2 by subject

cors <- matrix(NA, ncol = 2, nrow = n.resamp) ## storage
colnames(cors) <- c("r", "u")
ranefs <- vector("list", n.resamp)
converged <- integer(n.resamp)

beg.resamp <- Sys.time()
for (resamp.i in seq_len(n.resamp)) {
  
  .d2 <- d2  ## make copy of data
  
  for (subj.i in seq_along(.d2)) .d2[[subj.i]]$subj <- resamples[subj.i, resamp.i]
  .d2 <- do.call(rbind, .d2)
  .d <- rbind(d1, .d2)
  
  .fit <- update(fit, data = .d)  ## refit
  
  ## check for convergence
  
  converged[resamp.i] <- !isTRUE(all.equal(.fit@optinfo$conv$opt, 0))

  ## get random effects
  
  .tau <- VarCorr(.fit)$subj
  .r <- cov2cor(m %*% .tau %*% t(m))[1, 2]
  .u <- as.matrix(ranef(.fit)$subj) %*% t(m)
  .cor.u <- cor(.u)[1, 2]
  
  ## save
  
  cors[resamp.i, c("r", "u")] <- c(.r, .cor.u)
  ranefs[[resamp.i]] <- .u
  
  print(resamp.i)
  
}
time.resamp <- beg.resamp - Sys.time()

save.image()


## scratch ----

# fit1a <- lmer(
#   rt ~ trial.type * run + (trial.type * run | subj),
#   stroop.pro %>% filter(),
#   control = lmerControl(optim = "bobyqa", optCtrl = list(maxfun = 5E6))
#   )
# summary(fit1a)
# tau <- VarCorr(fit1a)$subj
# m <- rbind(c(0, 1, 0, 0), c(0, 1, 0, 1))
# cov2cor(m %*% tau %*% t(m))
# 
# fit1 <- lmer(
#   rt ~ trial.type:run + 0 + (0 + trial.type:run | subj),
#   stroop.pro,
#   control = lmerControl(optim = "bobyqa", optCtrl = list(maxfun = 5E6))
#   )
# summary(fit1)
# coef(fit1)
# m1 <- rbind(c(-1, 1, 0, 0), c(0, 0, -1, 1))
# tau1 <- VarCorr(fit1)$subj
# cov2cor(m1 %*% tau1 %*% t(m1))
# ests1 <- as.matrix(coef(fit1)$subj) %*% t(m1)
# cor(ests1)
# 
# 
# 
# ?confint
# boot.confint <- confint(fit1, method = "boot")
# 
# VarCorr(fit1)$subj
# 
# 
# 
# 
# ests <- as.matrix(coef(fit1)$subj) %*% t(m)
# plot(ests)
# cor(ests)
# 
# 
# ests.raw <- stroop.pro %>%
#   group_by(subj, trial.type, run) %>%
#   summarize(rt = mean(rt, na.rm = TRUE)) %>%
#   tidyr::spread(trial.type, rt) %>%
#   transmute(stroop = i - c, run = run) %>% 
#   tidyr::spread(run, stroop) %>%
#   rename(run1 = "1", run2 = "2") %>%
#   as.data.frame %>%
#   select(-subj)
# ests.raw
# # cor(ests.raw$run1, ests[, 1])
# # cor(ests.raw$run2, ests[, 2])
# 
# plot(ests.raw, pch = 16)
# points(ests, pch = 16, col = "firebrick")
# 
# cor(ests)
# cor(ests.raw)
# 
# fit2 <- glmer(
#   rt ~ trial.type * run + (trial.type * run | subj), 
#   stroop.pro, 
#   control = glmerControl(optCtrl = list(maxfun = 5E6)),
#   family  = inverse.gaussian("identity")
#   )
# summary(fit2)
# cov2cor(m %*% VarCorr(fit2)$subj %*% t(m))