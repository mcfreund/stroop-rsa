#' ## read data

#+ microphone-comparison_setup, include = FALSE

if (interactive()) {
  
  source(here::here("code", "packages.R"))
  source(here("code", "strings.R"))
  source(here("code", "funs.R"))
  source(here("code", "plots.R"))
  source(here("code", "read_atlases.R"))
  
  # behav <- fread(here("in", "behavior-and-events_group201902.csv")) %>% mutate(error = 1 - acc)
  behav <- fread(here("out", "behavior-and-events_group201902_with-subset-cols.csv"))
  cl1 <- lmeControl(maxIter = 1E5, msMaxIter = 1E5, niterEM = 1E5, msMaxEval = 1E5)
  
}

micinfo <- fread(here("in", "microphone-info_group201902.csv"), data.table = FALSE)
behav <- full_join(behav, micinfo, by = "subj")

#+



#' ## preliminary look

#+ r rawvals_all, fig.height = 15, fig.width = 20

behav %>%
  
  mutate(subj = factor(as.factor(subj), levels = micinfo[order(micinfo$mic), "subj"])) %>%
  
  ggplot(aes(subj, rt, color = mic)) +
  geom_point(position = position_jitter(width = 0.1), size = 0.5, alpha = 0.3) +
  geom_boxplot(position = position_nudge(1/3), notch = TRUE, width = 0.2) +
  
  scale_color_brewer(type = "qual", palette = 2) +
  
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0),
    legend.position = c(0.8, 0.8)
  )

#+ 

#+ rawvals

behav %>%
  
  group_by(mic, subj) %>%
  summarize(freq = sum(rt == 0, na.rm = TRUE)) %>%
  
  ggplot(aes(mic, freq)) +
  geom_point(size = 2) +
  
  labs(y = "frequency of rt == 0 per subj")

behav %>%
  
  group_by(mic, subj, error) %>%
  summarize(freq = sum(error, na.rm = TRUE)) %>%
  
  ggplot(aes(mic, freq)) +
  geom_point(size = 2, position = position_jitter(width = 0.1), alpha = 0.4) +
  
  labs(y = "frequency of errors per subj")

behav %>%
  group_by(mic, subj, acc.final) %>%
  summarize(freq = n()) %>%
  ggplot(aes(mic, freq)) +
  facet_grid(cols = vars(acc.final)) +
  geom_point(size = 2, position = position_jitter(width = 0.1), alpha = 0.4) +
  labs(y = "frequency of acc.final codings per subj")

#+ 


#+ qq, fig.height = 15, fig.width = 20

behav %>%
  
  ggplot(aes(sample = rt, color = mic)) +
  
  stat_qq(alpha = 0.3, size = 0.15) +
  stat_qq_line(size = 0.25) +
  
  facet_wrap(vars(subj)) +
  scale_color_brewer(type = "qual", palette = 2)

#+ 


#+ summary-stats

behav.mean.sd <- behav %>%
  
  group_by(mic, subj, session) %>%
  
  summarize(rt.mean = mean(rt, na.rm = TRUE), rt.sd = sd(rt, na.rm = TRUE))

grid.arrange(
  
  behav.mean.sd %>%
    ggplot(aes(mic, rt.mean)) +
    facet_grid(vars(session)) +
    geom_point(position = position_jitter(width = 0.1)) +
    geom_boxplot(position = position_nudge(1/3), notch = FALSE, width = 0.2),
  
  behav.mean.sd %>%
    ggplot(aes(mic, rt.sd)) +
    facet_grid(vars(session)) +
    geom_point(position = position_jitter(width = 0.1)) +
    geom_boxplot(position = position_nudge(1/3), notch = FALSE, width = 0.2)
  
)

#+



#' ## model

#' ### differences in mean stroop effect btw mics?
#' * leave in validation-set data for additional power

#+ means

behav.rt.mic <- behav %>% filter(!is.implausible.rt, !is.artifactual.rt, acc == 1, mic != "unknown")

fit <- lmer(
  rt ~ trial.type * mic + (trial.type | subj),
  behav.rt.mic
)
summary(fit)

# cl1 <- lmeControl(
#   maxIter = 100000, msMaxIter = 100000, niterEM = 100000,
#   msMaxEval = 100000, tolerance = 0.000001, msTol = 0.0000001, returnObject = TRUE,
#   minAbsParApVar = 0.05, opt = c("nlminb"), optimMethod = "BFGS"
# )

fit.het <- lme(
  rt ~ trial.type * mic, 
  random  = ~ trial.type | subj,
  weights = varIdent(form = ~ 1 | subj),
  data    = behav.rt.mic,
  control = cl1
)
summary(fit.het)

#' * fomri microphone recorded faster RTs
#' * no observed impact of microphone on size of mean stroop effect


## differences in variance btw microphones?

#' ### differences in variance between mics?
#+ vars

fit.het.ml  <- update(m.het, . ~ ., method = "ML")
fit.het.mic.ml <- update(fit.het.ml, . ~ ., weights = varIdent(form = ~ 1 | mic))
fit.hom.mic.ml <- update(fit.het.ml, . ~ ., weights = NULL)
anova(fit.het.ml, fit.het.mic.ml, fit.hom.mic.ml)

#+
#' * microoptics microphone recorded more variable RTs

