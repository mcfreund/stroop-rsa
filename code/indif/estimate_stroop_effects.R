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
source(here("code", "strings.R"))

## read and subset

stroop.pro <- read.csv(here("data", "behavior-and-events_group201902.csv"),  stringsAsFactors = FALSE) %>%
  filter(session == "pro", is.analysis.group) %>%
  mutate(trial.type = ifelse(trial.type == "i", "incon", "congr"))


## model ----

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

fit1.het.trim <- update(fit1.het, subset = !is.far.out)

## extract predictions

blups <- as.data.frame(coef(fit1.het.trim))
blups %<>% rename(congr = "(Intercept)", stroop = "trial.typeincon") %>% tibble::rownames_to_column("subj")


## estimate split-half reliability (cross-run) ----

fit1.het.run.reml.trim <- lme(
  rt ~ trial.type * run, stroop.pro.rt, ~ trial.type * run | subj,
  weights = varIdent(form = ~ 1 | subj),
  control = lmeControl(maxIter = 1e5, msMaxIter = 1e5, niterEM = 1e5, msMaxEval = 1e5),
  method = "REML",
  subset = !is.far.out
)

summary(fit1.het.run.reml.trim)

## conditional covariances

m <- rbind(c(0, 1, 0, 0), c(0, 1, 0, 1))  ## a contrast matrix that gives us what we want
u.trim <- as.matrix(coef(fit1.het.run.reml.trim)) %*% t(m)
cov(u.trim)
(r.conditional.trim <- cor(u.trim)[1, 2])  ## correlation of conditional modes

## unconditional / marginal covariances

tau.trim <- getVarCov(fit1.het.run.reml.trim)
(tau.trim <- m %*% tau.trim %*% t(m))

(r.marginal.trim <- cov2cor(tau.trim)[1, 2])

u.trim <- as.data.frame(u.trim)

names(u.trim) <- c("stroop.run1", "stroop.run2")

## errors
stroop.pro.er <- stroop.pro %>% mutate(error = 1 - acc)
fit.error0 <- glmer(error ~ trial.type + (1 | subj), stroop.pro.er, family = binomial)
fit.error1 <- glmer(error ~ trial.type + (trial.type | subj), stroop.pro.er, family = binomial)
anova()


behav %<>%
  full_join(
    data.frame(
      subj = rownames(coef(m.er1)$subj),
      er.logit.stroop = coef(m.er1)$subj$trial.typei,
      er.logit.congru = coef(m.er1)$subj[["(Intercept)"]]
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

plot.behav <- behav %>%
  mutate(subj = factor(subj, levels = subj[order(rt.mlm.hetvar, decreasing = TRUE)])) %>%
  select(subj, rt.mlm.hetvar, er.mlm) %>%
  melt %>%
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
saveRDS(fit1.het.trim, here("out", "behav", "fit1-het-trim_group201902.RDS"))
