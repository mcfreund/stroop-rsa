#' ## read data

#+ prelim-behav-vset_setup, include = FALSE

if (interactive()) {
  
  source(here::here("code", "packages.R"))
  source(here("code", "strings.R"))
  source(here("code", "funs.R"))
  source(here("code", "plots.R"))
  source(here("code", "read_atlases.R"))
  
  behav <- fread(here("in", "behavior-and-events_group201902.csv"))
  
}

#+


#' ## model, plot, and write

#+ model_RT_valid

## initial fit

behav.vset.rt <- behav %>% filter(acc == 1, !is.na(rt), rt < 3000, rt > 250, !is.analysis.group)

fit.vset <- lme(
  rt ~ trial.type, 
  random  = ~ trial.type | subj,
  data    = behav.vset.rt,
  weights = varIdent(form = ~ 1 | subj),
  control = lmeControl(maxIter = 1e5, msMaxIter = 1e5, niterEM = 1e5, msMaxEval = 1e5),
  method  = "REML"
)

## trim and re-fit

behav.vset.rt$resid.p <- resid(fit.vset, type = "p")
behav.vset.rt$is.far.out <- farout(behav.vset.rt$resid.p)
fit.vset.trim <- update(fit.vset, subset = !is.far.out)

## model error

fit.error.vset <- glmer(
  error ~ trial.type + (trial.type | subj), 
  behav %>% filter(response.final != "unintelligible", !is.analysis.group) %>% mutate(error = 1 - acc), 
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1E9))
)
summary(fit.error.vset)


## extract predictions

blups.vset <- as.data.frame(coef(fit.vset.trim))
blups.vset %<>% rename(congr = "(Intercept)", stroop = "trial.typeincon") %>% tibble::rownames_to_column("subj")

## bind with error blups

blups.vset %<>%
  
  full_join(
    
    data.frame(
      subj = rownames(coef(fit.error.vset)$subj),
      er.logit.stroop = coef(fit.error.vset)$subj$trial.typei,  ## extract logits
      er.logit.congr  = coef(fit.error.vset)$subj[["(Intercept)"]]
    ) %>%
      mutate(
        er.logit.incon = er.logit.stroop + er.logit.congr,  ## logit of error on incon trials
        ##  blup stroop effect in units percent error:
        stroop.er = (logit2prob(er.logit.incon) - logit2prob(er.logit.congr)) * 100
      ) %>%
      dplyr::select(subj, stroop.er),
    
    by = "subj"
    
  )


## plot

plot.behav.vset <- 
  
  blups.vset %>%
  
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
    axis.line.y     = element_line(size = axis.line.size*0.6),
    axis.text.y     = element_text(size = axis.text.size),
    axis.text.x     = element_blank(),
    axis.ticks.y    = element_line(size = axis.line.size),
    axis.ticks.x    = element_blank(),
    axis.title      = element_text(size = axis.title.size),
    axis.title.y = element_blank()
  )

plot.behav.vset

## write and save

ggsave(here("out", "behav", "stroop-blups-validation.pdf"), height = 2.5, width = figwidth/2)
write.csv(blups.vset, here("out", "behav", "stroop_blups_rt_group201902_validation.csv"),  row.names = FALSE)
saveRDS(fit.vset.trim, here("out", "behav", "fit1-het-trim_group201902_validation_set.RDS"))

#+