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
library(lme4)
library(ggplot2)
library(magrittr)

source(here("code", "strings.R"))


farout <- function(x) {
  
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr3 <- IQR(x) * 3
  
  x < (q1 - iqr3) | x > (q3 + iqr3)
  
}

logit2prob <- function(x) exp(x) / (1 + exp(x))

## read and subset

stroop.pro <- read.csv(here("in", "behavior-and-events_group201902.csv"),  stringsAsFactors = FALSE) %>%
  filter(session == "pro", is.analysis.group) %>%
  mutate(trial.type = ifelse(trial.type == "i", "incon", "congr"))


## model RT ----

## initial fit  

stroop.pro.rt <- stroop.pro %>% filter(acc == 1, !is.na(rt), rt < 3000, rt > 250)
is.weird.rt <- stroop.pro.rt$subj %in% c("849971", "161832") & stroop.pro.rt$rt < 500
stroop.pro.rt <- stroop.pro.rt[!is.weird.rt, ]

stroop.pro.er <- stroop.pro %>% filter(acc == 0)

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

fit1.het.trim <- update(fit1.het, data = stroop.pro.rt %>% filter(!is.far.out))
fit1.het.trim.ml <- update(fit1.het.trim, method = "ML")  ## for model comparisons


## test heterogeneity of variance

fit1.hom.trim.ml <- lme(
  rt ~ trial.type, 
  random  = ~ trial.type | subj,
  data    = stroop.pro.rt %>% filter(!is.far.out),
  control = lmeControl(maxIter = 1e5, msMaxIter = 1e5, niterEM = 1e5, msMaxEval = 1e5),
  method  = "ML"
)

rt.hom.v.het <- anova(fit1.hom.trim.ml, fit1.het.trim.ml)

## test stroop effect variance

fit0.het.trim.ml <- update(fit1.het.trim.ml, random  = ~ 1 | subj)
rt.stroopvar <- anova(fit0.het.trim.ml, fit1.het.trim.ml)

## estimate split-half reliability (across run)

fit1.het.run.trim <- lme(
  rt ~ trial.type * run, stroop.pro.rt %>% filter(!is.far.out), ~ trial.type * run | subj,
  weights = varIdent(form = ~ 1 | subj),
  control = lmeControl(maxIter = 1e5, msMaxIter = 1e5, niterEM = 1e5, msMaxEval = 1e5),
  method = "REML"
)

summary(fit1.het.run.trim)

## conditional covariances

m <- rbind(c(0, 1, 0, 0), c(0, 1, 0, 1))  ## a contrast matrix that gives us what we want
u.trim <- as.matrix(coef(fit1.het.run.trim)) %*% t(m)
cov(u.trim)
(r.conditional.trim <- cor(u.trim)[1, 2])  ## correlation of conditional modes

## unconditional / marginal covariances

tau.trim <- getVarCov(fit1.het.run.trim)
(tau.trim <- m %*% tau.trim %*% t(m))

(r.marginal.trim <- cov2cor(tau.trim)[1, 2])

u.trim <- as.data.frame(u.trim)
names(u.trim) <- c("stroop.run1", "stroop.run2")



## model errors ---

stroop.pro.er <- stroop.pro %>% mutate(error = 1 - acc)

## test for stroop effect variance

fit.error0 <- glmer(
  error ~ trial.type + (1 | subj),
  stroop.pro.er, 
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1E9))
  )
fit.error1 <- glmer(
  error ~ trial.type + (trial.type | subj), 
  stroop.pro.er, 
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1E9))
  )
er.stroopvar <- anova(fit.error0, fit.error1)



## wrangle estimates and plot figures ----

## get blups
blups <- as.data.frame(coef(fit1.het.trim))
blups %<>% rename(congr = "(Intercept)", stroop = "trial.typeincon") %>% tibble::rownames_to_column("subj")


blups %<>%
  full_join(
    data.frame(
      subj = rownames(coef(fit.error1)$subj),
      er.logit.stroop = coef(fit.error1)$subj$trial.typei,
      er.logit.congru = coef(fit.error1)$subj[["(Intercept)"]]
    ) %>%
      mutate(
        er.logit.incon = er.logit.stroop + er.logit.congru,
        er.odds = exp(er.logit.stroop),
        er = (logit2prob(er.logit.incon) - logit2prob(er.logit.congru)) * 100  ## percent error
      ) %>%
      select(subj, er.mlm = er),
    by = "subj"
  )


## draw figure ----

plot.behav <- blups %>%
  
  mutate(subj = factor(subj, levels = subj[order(stroop, decreasing = TRUE)])) %>%
  select(subj, stroop, er.mlm) %>%
  reshape2::melt() %>%
  
  ggplot(aes(subj, value)) +
  facet_grid(
    rows = vars(variable), scales = "free", switch = "y",
    labeller = as_labeller(c(rt.mlm.hetvar = "response time", er.mlm = "% error"))
  ) +
  
  geom_segment(aes(x = subj, y = 0, xend = subj, yend = value), color = "grey50") +
  geom_point(fill = "grey30", color = "black", shape = 21, size = 1) +
  
  xlab("participant") +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "mm"),
    strip.background = element_blank(),
    strip.text = element_text(face = "italic", size = 5),
    strip.placement = "outside",
    axis.title.x = element_text(face = "italic", size = 5),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 4),
    panel.background = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(2))


u.trim %>%
  
  ggplot(aes(stroop.run1, stroop.run2)) +
  geom_abline() +
  stat_boot_ci(n = 1E4, alpha = 0.2) +
  # geom_smooth(method = "lm", se = FALSE, color = "grey30") +
  geom_point(shape = 21, color = "white", fill = "black", size = 2) +
  annotate(
    geom = "text", label = paste0("r(unconditional) = ", round(r.marginal.trim, 2)), 
    y = 150, x = 20, vjust = 1, hjust = 0
  ) +
  theme(panel.grid = element_blank())


## write ----

write.csv(blups, here("out", "behav", "stroop_blups_rt_group201902.csv"),  row.names = FALSE)  ## blups
behav.mod.objs <- list(
  fit1.het.trim = fit1.het.trim,
  er.stroopvar = er.stroopvar,
  rt.stroopvar = rt.stroopvar,
  rt.hom.v.het = rt.hom.v.het,
  r.marginal.trim = r.marginal.trim,
  u.trim = u.trim
)
saveRDS(behav.mod.objs, here("out", "behav", "mod_objs.RDS"))
