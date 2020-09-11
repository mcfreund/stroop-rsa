#' ## RT analysis

#+ prelim-behav_setup, include = FALSE

if (interactive()) {
  
  source(here::here("code", "packages.R"))
  source(here("code", "strings.R"))
  source(here("code", "funs.R"))
  source(here("code", "plots.R"))
  source(here("code", "read_atlases.R"))
  
  behav <- fread(here("in", "behavior-and-events_group201902.csv"))
  
}


#' ### plot
#+ prelim-behav_plot-rt

plot(behav$rt)
abline(h = 250)
abline(h = 3000)  ## marks beginning of subsequent trial

behav$is.implausible.rt <- behav$rt > 3000 | behav$rt < 250

grid.arrange(
  
  behav %>%
    
    filter(acc == 1, !is.implausible.rt) %>%
    
    ggplot(aes(rt)) +
    geom_density(fill = "slateblue", alpha = 0.3) +
    
    labs(title = "all correct trials") +
    theme(legend.position = "none"),
  
  behav %>%
    
    filter(acc == 1, !is.implausible.rt) %>%
    
    ggplot(aes(rt)) +
    geom_density(aes(fill = trial.type), alpha = 0.3) +
    
    labs(title = "by congruency") +
    scale_fill_brewer() +
    scale_color_brewer() +
    theme(legend.position = c(0.5, 0.5)),

  ncol = 2
  
)
#+

#+ prelim-behav_plot-rt-qq, fig.height = 10

behav %>%
  
  filter(acc == 1, !is.implausible.rt) %>%
  
  full_join(group_by(., subj) %>% summarize(r2norm = qqr2(rt)), by = "subj") %>%
  
  ggplot(aes(sample = rt)) +
  stat_qq(alpha = 0.8, size = 1) +
  stat_qq_line(size = 0.25) +
  
  facet_wrap(vars(subj)) +
  geom_text(
    aes(
      x = -1.25, y = 2000,
      label = paste0("r^2 = ", round(r2norm, 3))
    ), size = 3, color = "grey50"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.title = element_blank()) +
  labs(
    title    = "QQ rt versus normal",
    subtitle = paste0("overall r^2 with normal = ", round(qqr2(behav$rt), 3))
  )

#' 849971 and 161832 have odd patterns in RT distribution.
#' Strong clustering at 500 ms.
#' This is highly consistent with an artifact we found in other RT data from the same microphone/preprocessing method.
#+

#' ### handling likely-artifactual RTs
#+ prelim-behav_plot-rt-artifacts

behav %>%
  
  filter(acc == 1, !is.implausible.rt) %>%
  
  ggplot(aes(x = rt, group = subj)) +
  geom_density(data = . %>% filter(!subj %in% c("849971", "161832")), color = rgb(0, 0, 0, 0.15), size = 2) +
  geom_density(data = . %>% filter(subj %in% c("849971", "161832")), color = "firebrick1", size = 2) +
  labs(
    title    = "subjects with bimodal distributions (likely artifactual) in red"
  )

behav$is.artifactual.rt <- FALSE
behav$is.artifactual.rt[behav$subj %in% c("849971", "161832") & behav$rt < 500] <- TRUE

behav %>%

  filter(acc == 1, !is.implausible.rt, subj %in% c("849971", "161832")) %>%
  full_join(group_by(., subj) %>% summarize(r2norm = qqr2(rt)), by = "subj") %>%
  
  ggplot(aes(sample = rt)) +
  stat_qq(alpha = 0.8, size = 0.4) +
  stat_qq_line(size = 0.25) +
  geom_hline(yintercept = 500) +
  facet_wrap(vars(subj)) +
  theme_minimal(base_size = 10)
#+


#' ### fit RT model
#+ prelim-behav_rt-model

## make sure to exclude validation set!
behav.rt.aset <- behav %>% filter(acc == 1, !is.implausible.rt, !is.artifactual.rt, is.analysis.group)

fit1.het <- lme(
  rt ~ trial.type, 
  random  = ~ trial.type | subj,
  data    = behav.rt.aset,
  weights = varIdent(form = ~ 1 | subj),
  control = lmeControl(maxIter = 1E9, msMaxIter = 1E9, niterEM = 1E9, msMaxEval = 1E9),
  method  = "REML"
)

## define "far" outliers as points with resids more extreme than 3 IQR

behav.rt.aset$resid.p <- resid(fit1.het, type = "p")
behav.rt.aset$is.far.out <- farout(behav.rt.aset$resid.p)

## exclude far outliers and refit model

fit1.het.trim <- update(fit1.het, data = behav.rt.aset %>% filter(!is.far.out))
fit1.het.trim.ml <- update(fit1.het.trim, method = "ML")  ## for model comparisons
#+

#' ### test heterogeneity of level-I variance
#+ prelim-behav_test-hetvar

fit1.hom.trim.ml <- lme(
  rt ~ trial.type, 
  random  = ~ trial.type | subj,
  data    = behav.rt.aset %>% filter(!is.far.out),
  control = lmeControl(maxIter = 1E9, msMaxIter = 1E9, niterEM = 1E9, msMaxEval = 1E9),
  method  = "ML"
)

(rt.hom.v.het <- anova(fit1.hom.trim.ml, fit1.het.trim.ml))
#+

#' ### test variance in stroop effect
#+ prelim-behav_test-stroopvar-rt

fit0.het.trim.ml <- update(fit1.het.trim.ml, random  = ~ 1 | subj)
rt.stroopvar <- anova(fit0.het.trim.ml, fit1.het.trim.ml)

#' ### estimate cross-run reliability
#+ prelim-behav_est-reliability-rt
fit1.het.run.trim <- lme(
  rt ~ trial.type * run, 
  random  = ~ trial.type * run | subj,
  data    = behav.rt.aset %>% filter(!is.far.out),
  weights = varIdent(form = ~ 1 | subj),
  control = lmeControl(maxIter = 1E9, msMaxIter = 1E9, niterEM = 1E9, msMaxEval = 1E9),
  method  = "REML"
)
summary(fit1.het.run.trim)

## get unconditional / marginal covariances

m <- rbind(c(0, 1, 0, 0), c(0, 1, 0, 1))  ## contrast matrix
tau.trim <- getVarCov(fit1.het.run.trim)
(tau.trim <- m %*% tau.trim %*% t(m))  ## covariance matrix
(r.marginal.trim <- cov2cor(tau.trim)[1, 2])  ## correlation
svd(getVarCov(fit1.het.trim), nv = 0L)$d  ## covmat not degenerate (eigvals > 0)
#+

#' ## error analysis
#' ### test variance in stroop effect

#+ prelim-behav_test-stroopvar-er
behav <- behav %>% mutate(error = 1 - acc)

fit.error0 <- glmer(
  error ~ trial.type + (1 | subj),
  behav %>% filter(response.final != "unintelligible", is.analysis.group), 
  family  = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1E9))
)
fit.error1 <- glmer(
  error ~ trial.type + (trial.type | subj), 
  behav %>% filter(response.final != "unintelligible", is.analysis.group), 
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1E9))
)
(er.stroopvar <- anova(fit.error0, fit.error1))
#+

#' ### estimate cross-run reliability
#+ prelim-behav_est-reliability-er
fit.error1.run <- glmer(
  error ~ trial.type * run + (trial.type * run | subj),
  behav %>% filter(response.final != "unintelligible", is.analysis.group),
  family = binomial,
  control = glmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 1E9))
)
summary(fit.error1.run)
plot(rePCA(fit.error1.run)$subj)
#+

#' ## extract blups and draw figures

#+ prelim-behav_blups
## get RT blups
blups <- as.data.frame(coef(fit1.het.trim))
blups %<>% rename(congr = "(Intercept)", stroop = "trial.typeincon") %>% tibble::rownames_to_column("subj")

## bind with error blups

blups %<>%
  
  full_join(
    
    data.frame(
      subj = rownames(coef(fit.error1)$subj),
      er.logit.stroop = coef(fit.error1)$subj$trial.typei,  ## extract logits
      er.logit.congr  = coef(fit.error1)$subj[["(Intercept)"]]
    ) %>%
      mutate(
        er.logit.incon = er.logit.stroop + er.logit.congr,  ## logit of error on incon trials
        ##  blup stroop effect in units percent error:
        stroop.er = (logit2prob(er.logit.incon) - logit2prob(er.logit.congr)) * 100
      ) %>%
      dplyr::select(subj, stroop.er),
    
    by = "subj"
    
  )


## draw figure

plot.behav <- 
  
  blups %>%
  
  mutate(subj = factor(subj, levels = subj[order(stroop, decreasing = TRUE)])) %>%
  dplyr::select(subj, stroop, stroop.er) %>%
  reshape2::melt() %>%
  
  
  filter(!is.na(value)) %>%
  
  ggplot(aes(subj, value)) +
  facet_grid(
    rows = vars(variable), scales = "free", switch = "y",
    labeller = as_labeller(c(stroop = "Response time", stroop.er = "% error"))
  ) +
  
  geom_segment(aes(x = subj, y = 0, xend = subj, yend = value), color = "grey50", size = geom.line.size) +
  geom_point(fill = "grey30", color = "black", shape = 21, size = geom.point.size) +
  coord_capped_cart(left = "both") +
  
  xlab("Subject") +
  theme(
    panel.grid       = element_blank(), 
    panel.border     = element_blank(),
    panel.background = element_blank(),
    strip.placement  = "outside",
    strip.background = element_blank(),
    strip.text       = element_text(size = axis.title.size),
    axis.line.y     = element_line(size = axis.line.size),
    axis.text.y     = element_text(size = axis.text.size),
    axis.text.x     = element_blank(),
    axis.ticks.y    = element_line(size = axis.line.size),
    axis.ticks.x    = element_blank(),
    axis.title      = element_text(size = axis.title.size),
    axis.title.y = element_blank()
  )

ggsave()

## write

write.csv(blups, here("out", "behav", "stroop_blups_rt_group201902.csv"),  row.names = FALSE)
write.csv(behav, here("out", "behav", "behavior-and-events_group201902_with-subset-cols.csv"),  row.names = FALSE)

behav.mod.objs <- list(
  
  fit1.het.trim = fit1.het.trim,
  er.stroopvar = er.stroopvar,
  rt.stroopvar = rt.stroopvar,
  rt.hom.v.het = rt.hom.v.het,
  r.marginal.trim = r.marginal.trim
  
)
saveRDS(behav.mod.objs, here("out", "behav", "mod_objs.RDS"))
#+
