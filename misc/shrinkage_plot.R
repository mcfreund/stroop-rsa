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
u <- as.matrix(coef(fit)$subj) %*% t(m)
cor.u <- cor(u)[1, 2]


## bootstrap confint ----

beg.confint <- Sys.time()
ci <- confint(fit, method = "boot")
time.confint <- beg.confint - Sys.time()


## randomization test ----

n.resamp <- 1E5 + 500
resamples <- replicate(n.resamp, sample(unique(stroop.pro$subj)))
resamples <- unique(resamples, MARGIN = 2)  ## remove duplicate permutations (prob. none)

d1 <- stroop.pro[stroop.pro$run == "run1", ]  ## split data by run
d2 <- stroop.pro[stroop.pro$run == "run2", ]
d2 <- split(d2, d2$subj)  ## split run 2 by subject

cors <- matrix(NA, ncol = 2, nrow = n.resamp) ## storage
colnames(cors) <- c("r", "u")
ranefs <- vector("list", n.resamp)
converged <- logical(n.resamp)

beg.resamp <- Sys.time()
for (resamp.i in seq_len(n.resamp)) {
  
  .d2 <- d2  ## make copy of data
  
  for (subj.i in seq_along(.d2)) .d2[[subj.i]]$subj <- resamples[subj.i, resamp.i]
  .d2 <- do.call(rbind, .d2)
  .d <- rbind(d1, .d2)
  
  .fit <- update(fit, data = .d)  ## refit
  
  ## check for convergence
  
  converged[resamp.i] <- is.null(fit@optinfo$conv$lme4$code)

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

## load and look ----

cors <- cors[!is.na(cors[, 1]), ]
ranefs <- ranefs[!is.na(cors[, 1])]


ols <- stroop.pro %>%
  group_by(subj, trial.type, run) %>%
  summarize(rt = mean(log(rt))) %>%
  tidyr::spread(trial.type, rt) %>%
  transmute(stroop = i - c, run = run) %>%
  tidyr::spread(run, stroop) %>%
  as.data.frame

u <- as.data.frame(u)
names(u) <- c("run1", "run2")
u <- tibble::rownames_to_column(u, "subj")

ests <- bind_rows(u %>% mutate(est = "eb"), ols %>% mutate(est = "ols"))

library(ggplot2)

scatter <- ests %>%
  ggplot(aes(run1, run2)) +
  geom_abline() +
  geom_line(aes(group = subj), color = "grey60") +
  geom_point(aes(fill = est), color = "grey40", pch = 21, size = 2) +
  scale_fill_manual(values = c(ols = "black", eb = "firebrick2")) +
  theme(
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    axis.title = element_text(face = 'bold', color = "grey40")
  ) +
  geom_segment(y = -Inf, yend = -Inf, x = min(ests$run1), xend = max(ests$run1), color = "grey40") +
  geom_segment(x = -Inf, xend = -Inf, y = min(ests$run2), yend = max(ests$run2), color = "grey40") +
  scale_x_continuous(breaks = c(round(min(ests$run1), 2), 0, round(max(ests$run1), 2))) +
  scale_y_continuous(breaks = c(round(min(ests$run2), 2), 0, round(max(ests$run2), 2))) +
  labs(x = "run 1 stroop (log RT)", y = "run 2 stroop (log RT)") +
  annotate(
    "text", x = -0.05, y = Inf, 
    label = paste0("ols (r = ", round(cor(ols$run1, ols$run2), 2), ")"),
    hjust = 0, vjust = 1, size = 5
    ) +
  annotate(
    "text", x = -0.05, y = 0.255, 
    label = paste0("empirical bayes (r = ", round(cor(u$run1, u$run2), 2), ")"),
    hjust = 0, vjust = 1, color = "firebrick2", size = 5
  )
  
cors <- cors %>% as.data.frame
hyptest <- cors %>%
  ggplot(aes(y = r)) +
  geom_boxplot(notch = TRUE, fill = "steelblue", notchwidth = 0.1, outlier.colour = "steelblue") +
  theme(
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = 'bold', color = "grey40")
  ) +
  # geom_segment(x = -Inf, xend = -Inf, y = min(cors$r), yend = max(cors$r)) +
  scale_y_continuous(breaks = c(round(min(cors$r), 2), 0, round(max(cors$r), 2))) +
  annotate("point", x = 0, y = cor.u, color = "firebrick2", size = 5, shape = 8) +
  annotate("text", x = 0, y = 0.75, color = "steelblue", label = "null") +
  labs(y = "split-half correlation")


p <- gridExtra::grid.arrange(scatter, hyptest, ncol = 2, widths = c(1, 0.15))

ggsave("ols-vs-eb.pdf", p, device = "pdf", width = 9, height = 6)


## response transformations
# 
# fit.raw <- lmer(
#   rt ~ trial.type + (trial.type | subj),
#   stroop.pro,
#   control = lmerControl(optim = "bobyqa", optCtrl = list(maxfun = 1E5))
# )
# summary(fit.raw)
# plot(fit.raw)
# 
# stroop.pro$resid <- resid(fit.raw)
# 
# stroop.pro %<>%
#   group_by(subj) %>%
#   mutate(
#     sdev = sd(resid),
#     rtbar = mean(rt),
#     rt.scale = rt / sdev,
#     rt.center = rt - rtbar,
#     rt.z = rt.center / sdev
#     )
# 
# fit.rtbar <- lmer(rt ~ trial.type * rtbar + (trial.type | subj), stroop.pro)
# summary(fit.rtbar)
# plot(fit.rtbar)
# 
# fit.center <- lmer(rt.center ~ trial.type + (trial.type | subj), stroop.pro)
# summary(fit.center)
# plot(fit.center)
# 
# fit.sdev <- lmer(rt ~ trial.type * sdev + (trial.type | subj), stroop.pro)
# summary(fit.sdev)
# plot(fit.sdev)
# 
# fit.scale <- lmer(
#   rt.scale ~ trial.type + (trial.type | subj), stroop.pro
# )
# fit.scale.red <- lmer(
#   rt.scale ~ trial.type + (1 | subj), stroop.pro
# )
# anova(fit.scale, fit.scale.red)
# summary(fit.scale)
# 
# 
# 
# plot(fit.scale)
# plot(
#   coef(fit.scale)$subj$trial.typei,
#   coef(fit.rtbar)$subj$trial.typei
# )
# 
# fit.z <- lmer(
#   rt.z ~ trial.type + (trial.type | subj), stroop.pro,
#   control = lmerControl(optim = "bobyqa", optCtrl = list(maxfun = 1E5))
#   )
# summary(fit.z)
# plot(fit.z)
# 
# fit.scale2 <- lmer(
#   rt.scale ~ trial.type * run + (trial.type * run | subj), stroop.pro,
#   control = lmerControl(optim = "bobyqa", optCtrl = list(maxfun = 1E5))
# )
# summary(fit.scale2)
# tau <- VarCorr(fit.scale2)$subj
# r <- cov2cor(m %*% tau %*% t(m))[1, 2]
# u <- as.matrix(coef(fit.scale2)$subj) %*% t(m)
# cor.u <- cor(u)[1, 2]
# 
# plot(u)



