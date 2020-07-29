library(mikeutils)
library(magrittr)
library(here)
library(knitr)
library(dplyr)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(colorspace)
library(viridis)
library(nlme)
library(caret)
library(gtools)
library(vegan)
library(foreach)
library(doParallel)
library(lme4)
library(nlme)

source(here("code", "strings.R"))
source(here("code", "funs.R"))

theme_set(theme_bw(base_size = 12))

logit2prob <- function(x) exp(x) / (1 + exp(x))


## read data

blups <- read.csv(here("out", "behav", "stroop_blups_rt_group201902.csv"), stringsAsFactors = FALSE)
stats.subjs.tdic <- fread(
  here("out", "rsa", "stats", paste0("subjs_pro_bias_acc-only_masks_pearson_residual_glm-tdic.csv"))
)

## subset and bind

stats.subjs.tdic <- stats.subjs.tdic[is.analysis.group == TRUE, ]  ## EXCLUDE HELD OUT SUBJECTS!
stats.subjs.tdic <- stats.subjs.tdic[y == "rank" & param %in% c("target", "distractor", "incongruency"), ]
stats.subjs.tdic <- stats.subjs.tdic[, c("coef", "y", "model") := NULL]  ## remove useless cols
stats.subjs.tdic <- full_join(blups, stats.subjs.tdic, by = "subj")

## format cols

stats.subjs.tdic <- cbind(
  stats.subjs.tdic,
  reshape2::colsplit(stats.subjs.tdic$roi, "_", c("roi.set", "superparcel"))
)
stats.subjs.tdic$superparcel[stats.subjs.tdic$superparcel == ""] <- "vwfa"

## strings

colors.model <- c(incongruency = "#d95f02", target = "#1b9e77", distractor = "#7570b3")
colors.congruency <- c(I = "#d01c8b", C = "#4dac26")
colors.targets <- c(blue = "#08519c", red = "#a50f15", white = "grey50", purple = "#54278f")
  

d <- stats.subjs.tdic %>% filter(
  roi.set == "anatfunc",
  superparcel %in% c("dlpfc_L", "dlpfc_R", "mfc_L", "mfc_R", "lppc_L", "lppc_R", "vvis_L", "smmouth")
)
d <- full_join(
  d, 
  stats.subjs.tdic %>% filter(
    roi.set == "anat",
    superparcel %in% c("smmouth")
  )
)
  

d %<>% mutate(id = as.character(interaction(superparcel, param)))

w <- d %>%
  ungroup %>%
  dplyr::select(subj, stroop, id, beta) %>%
  tidyr::spread(id, beta)


## plot stroop effs ----

## read and subset

stroop.pro <- read.csv(here("data", "behavior-and-events_group201902.csv"),  stringsAsFactors = FALSE) %>%
  filter(session == "pro", is.analysis.group) %>%
  mutate(trial.type = ifelse(trial.type == "i", "incon", "congr"))


## model errors ----

stroop.pro.er <- stroop.pro %>% mutate(error = 1 - acc)

fit.error0 <- glmer(error ~ trial.type + (1 | subj), stroop.pro.er, family = binomial)
fit.error1 <- glmer(error ~ trial.type + (trial.type | subj), stroop.pro.er, family = binomial)
anova(fit.error0, fit.error1)


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

fit1.het.trim <- update(fit1.het, subset = !is.far.out)

fit0.het.trim <- lme(
  rt ~ trial.type, 
  random  = ~ 1 | subj,
  data    = stroop.pro.rt,
  weights = varIdent(form = ~ 1 | subj),
  control = lmeControl(maxIter = 1e5, msMaxIter = 1e5, niterEM = 1e5, msMaxEval = 1e5),
  method  = "REML",
  subset = !is.far.out
)

anova(fit1.het.trim, fit0.het.trim)


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
      select(subj, stroop.error = er),
    by = "subj"
  )

p.behav <- blups %>%
  mutate(subj = factor(subj, levels = subj[order(stroop, decreasing = TRUE)])) %>%
  select(subj, stroop, stroop.error) %>%
  reshape2::melt() %>%
  ggplot(aes(subj, value)) +
  facet_grid(
    rows = vars(variable), scales = "free", switch = "y",
    labeller = as_labeller(c(stroop = "RT", stroop.error = "% error"))
  ) +
  geom_segment(aes(x = subj, y = 0, xend = subj, yend = value), color = "grey50", size = 1) +
  geom_point(fill = "grey30", color = "black", shape = 21, size = 3) +
  xlab("participant") +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "mm"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = rel(3)),
    strip.placement = "outside",
    axis.title.x = element_text(face = "bold", size = rel(3)),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = rel(2)),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(size = rel(4))
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(2))


ggsave(
  here("out", "figs", "ohbm", "behav.pdf"),
  plot = p.behav,
  units = "cm",
  device = "pdf",
  height = 12,
  width = 20
)



## (a) brains ----

## see write_masks.R


## (b) scatterplot ----


## bootstrapped confidence interval for bivariate correlation

skipcor <- function(allx, ally, nboot = 10000, ...) {
  ## for bivariate cases only
  # allx = df$stroop.er.mean
  # ally = df$silhouette
  # if (ncol(allx) > 1 | ncol(ally) > 1) stop("need column vecs")
  allm <- cbind(allx, ally)
  # keep.obs <- outpro(allm)$keep
  keep.obs <- 1:nrow(allm)
  # keep.obs <- outpro(allm, cop = 2, plotit = FALSE)$keep
  nkeep <- length(keep.obs)
  # print(noquote(paste0("excluded ", nrow(allm) - nkeep, " obs as outliers")))
  x <- allm[keep.obs, 1]
  y <- allm[keep.obs, 2]
  rho <- cor(x, y, method = "spearman")
  r <- cor(x, y)
  ev <- c(spearman = rho, pearson = r)
  if (nboot < 1) return(ev)
  samps <- matrix(
    sample.int(nkeep, nkeep * nboot, replace = TRUE),
    nrow = nkeep, ncol = nboot
  )
  dist.r <- apply(samps, 2, function(s) cor(x[s], y[s]))
  dist.rho <- apply(samps, 2, function(s) cor(x[s], y[s], method = "spearman"))
  ## frequentist's p:
  lb <- (0.04 * nboot) / 2
  # lb <- (1 / 20 * nboot) / 2
  ub <- nboot - lb
  dist.rho <- sort(dist.rho)
  dist.r <- sort(dist.r)
  ci <- rbind(
    spearman = c(lb = dist.rho[lb], ub = dist.rho[ub]),
    pearson  = c(lb = dist.r[lb], ub = dist.r[ub])
  )
  cbind(ev, ci)
}



w$lfp_R.target <- (w$lppc_R.target + w$dlpfc_R.target) / 2


cor.lfp_R.target <- round(skipcor(w$lfp_R.target, w$stroop), 2)
cor.mfc_L.incongruency <- round(skipcor(w$mfc_L.incongruency, w$stroop), 2)
cor.dlpfc_R.target <- round(skipcor(w$dlpfc_R.target, w$stroop), 2)
cor.lppc_R.target <- round(skipcor(w$lppc_R.target, w$stroop), 2)


p.dissoc.scatter <-
  w %>%
  
  ggplot(aes(y = stroop)) +
  
  stat_boot_ci(aes(x = lfp_R.target), alpha = 0.3, n = 1E4, percent = 96, fill = colors.model["target"]) +
  stat_boot_ci(aes(x = mfc_L.incongruency), alpha = 0.3, n = 1E4, percent = 96, fill = colors.model["incongruency"]) +
  geom_smooth(aes(x = lfp_R.target), alpha = 0.3, color = colors.model["target"], method = "lm", se = FALSE) +
  geom_smooth(aes(x = mfc_L.incongruency), alpha = 0.3, color = colors.model["incongruency"], method = "lm", se = FALSE) +
  geom_point(aes(x = lfp_R.target), shape = 21, fill = colors.model["target"], color = "white", size = 6) +
  geom_point(aes(x = mfc_L.incongruency), shape = 21, fill = colors.model["incongruency"], color = "white", size = 6) +
  
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_blank(),
    panel.border =  element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = rel(2)),
    axis.text = element_text(size = rel(2)),
    axis.title = element_text(size = rel(2), face = "bold"),
    axis.ticks = element_blank()
  ) +
  
  scale_x_continuous(breaks = c(min(w$mfc_L.incongruency), 0, max(w$mfc_L.incongruency)) %>% round(2)) +
  scale_y_continuous(breaks = c(min(w$stroop), max(w$stroop)) %>% round) +
  geom_segment(
    aes(y = round(min(w$stroop), 2), yend = round(max(w$stroop), 2), x = -Inf, xend = -Inf),
    color = "grey40", size = 1
    ) +
  geom_segment(
    aes(x = round(min(w$mfc_L.incongruency), 2), xend = round(max(w$mfc_L.incongruency), 2), y = -Inf, yend = -Inf),
    color = "grey40", size = 1
  ) +
  
  labs(x = "RSA model fit (\u03B2)", y = "Stroop (ms)") +
  
  annotate(
    geom = "text", x = -0.35, y = 140, 
    label = paste0(
      "dMFC incongruency: r = ", 
      cor.mfc_L.incongruency["pearson", "ev"], 
      ", 96% CI [", cor.mfc_L.incongruency["pearson", "lb"], ", ", cor.mfc_L.incongruency["pearson", "ub"], "]"
      ),
    size = 8,
    hjust = 0,
    fontface = "italic",
    color = colors.model["incongruency"]
    ) +
    
    annotate(
      geom = "text", x = -0.35, y = 130, 
      label = paste0(
        "dLPFC+IPS target: r = ", 
        cor.lfp_R.target["pearson", "ev"], 
        " [", cor.lfp_R.target["pearson", "lb"], ", ", cor.lfp_R.target["pearson", "ub"], "]"
      ),
      size = 8,
      hjust = 0,
      fontface = "italic",
      color = colors.model["target"]
    )

ggsave(
  here("out", "figs", "ohbm", "indiv_predicted_dissoc_main.pdf"),
  plot = p.dissoc.scatter,
  units = "cm",
  device = "pdf",
  height = 15,
  width = 20
)


l <- w %>%
  
  select(subj, dlpfc_R.incongruency, dlpfc_R.target, lppc_R.incongruency, lppc_R.target, mfc_L.target, mfc_L.incongruency) %>%
  reshape2::melt() %>%
  full_join(blups)

l <- bind_cols(
  l, 
  reshape2::colsplit(as.character(l$variable), "\\.", c("parcel", "param"))
) 


l$parcel[l$parcel == "dlpfc_R"] <- "dlPFC (R)"
l$parcel[l$parcel == "lppc_R"] <- "IPS (R)"
l$parcel[l$parcel == "mfc_L"] <- "mFC (L)"

p.latmed <- l %>%
  
  filter(param %in% c("target", "incongruency")) %>%
  
  ggplot(aes(value, y = stroop, fill = param, color = param)) +
  stat_boot_ci(alpha = 0.3, n = 1E4, percent = 96, color = "transparent") +
  geom_smooth(alpha = 0.3, method = "lm", se = FALSE) +
  geom_point(shape = 21, color = "white", size = 3) +
  
  facet_wrap(vars(parcel), nrow = 1) +
  
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_blank(),
    panel.border =  element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = rel(2)),
    axis.text = element_text(size = rel(2)),
    axis.title = element_text(size = rel(2), face = "bold"),
    axis.ticks = element_blank()
  ) +
  
  scale_x_continuous(breaks = c(min(w$mfc_L.incongruency), 0, max(w$mfc_L.incongruency)) %>% round(2)) +
  scale_y_continuous(breaks = c(min(w$stroop), max(w$stroop)) %>% round) +
  
  geom_segment(
    aes(y = round(min(w$stroop)), yend = round(max(w$stroop)), x = -Inf, xend = -Inf),
    color = "grey40", size = 1
  ) +
  geom_segment(
    aes(x = round(min(w$mfc_L.incongruency), 2), xend = round(max(w$mfc_L.incongruency), 2), y = -Inf, yend = -Inf),
    color = "grey40", size = 1
  ) +
  
  scale_fill_manual(values = colors.model) +
  scale_color_manual(values = colors.model) +
  
  labs(x = "RSA model fit (\u03B2)", y = "stroop (ms)")


ggsave(
  here("out", "figs", "ohbm", "indiv_predicted_dissoc_all.pdf"),
  plot = p.latmed,
  units = "cm",
  device = "pdf",
  height = 10,
  width = 27
)


## exploratory analysis ----


l1 <- w %>%
  
  select(subj, dlpfc_L.distractor, vvis_L.incongruency, smmouth.distractor) %>%
  reshape2::melt() %>%
  full_join(blups)

l1 <- bind_cols(
  l1, 
  reshape2::colsplit(as.character(l1$variable), "\\.", c("parcel", "param"))
) 


l1$parcel[l1$parcel == "dlpfc_L"] <- "dlPFC (L)"
l1$parcel[l1$parcel == "vvis_L"] <- "ventral visual (L)"
l1$parcel[l1$parcel == "smmouth"] <- "SM-mouth (bil.)"

p.explor <-
  l1 %>%
  
  ggplot(aes(value, y = stroop, fill = param, color = param)) +
  stat_boot_ci(alpha = 0.3, n = 1E4, percent = 96, color = "transparent") +
  geom_smooth(alpha = 0.3, method = "lm", se = FALSE) +
  geom_point(shape = 21, color = "white", size = 3) +
  
  facet_wrap(vars(parcel), nrow = 1) +
  
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_blank(),
    panel.border =  element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = rel(2)),
    axis.text = element_text(size = rel(2)),
    axis.title = element_text(size = rel(2), face = "bold"),
    axis.ticks = element_blank()
  ) +
  
  scale_x_continuous(breaks = c(min(w$vvis_L.incongruency), 0, max(w$vvis_L.incongruency)) %>% round(2)) +
  scale_y_continuous(breaks = c(min(w$stroop), max(w$stroop)) %>% round) +
  
  geom_segment(
    aes(y = round(min(w$stroop)), yend = round(max(w$stroop)), x = -Inf, xend = -Inf),
    color = "grey40", size = 1
  ) +
  geom_segment(
    aes(x = round(min(w$mfc_L.incongruency), 2), xend = round(max(w$mfc_L.incongruency), 2), y = -Inf, yend = -Inf),
    color = "grey40", size = 1
  ) +
  
  scale_fill_manual(values = colors.model) +
  scale_color_manual(values = colors.model) +
  
  labs(x = "RSA model fit (\u03B2)", y = "stroop (ms)")


ggsave(
  here("out", "figs", "ohbm", "indiv_explor_dissoc.pdf"),
  plot = p.explor,
  units = "cm",
  device = "pdf",
  height = 10,
  width = 27
)




cor.vvis_L.incongruency <- round(skipcor(w$vvis_L.incongruency, w$stroop), 2)
cor.dlpfc_L.distractor <- round(skipcor(w$dlpfc_L.distractor, w$stroop), 2)
cor.smmouth.distractor <- round(skipcor(w$smmouth.distractor, w$stroop), 2)



fit.selected.scale <- lm(
  scale(stroop) ~ scale(dlpfc_L.distractor) + scale(vvis_L.incongruency) + scale(lppc_R.target) + scale(mfc_L.incongruency), w
)
coef(summary(fit.selected.scale))



## ----


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

fit0.het.trim <- lme(
  rt ~ trial.type, 
  random  = ~ 1 | subj,
  data    = stroop.pro.rt,
  weights = varIdent(form = ~ 1 | subj),
  control = lmeControl(maxIter = 1e5, msMaxIter = 1e5, niterEM = 1e5, msMaxEval = 1e5),
  method  = "REML",
  subset = !is.far.out
)

anova(fit1.het.trim, fit0.het)

## errors
stroop.pro.er <- stroop.pro %>% mutate(error = 1 - acc)
fit.error0 <- glmer(error ~ trial.type + (1 | subj), stroop.pro.er, family = binomial)
fit.error1 <- glmer(error ~ trial.type + (trial.type | subj), stroop.pro.er, family = binomial)
anova()
