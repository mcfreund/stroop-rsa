library(here)
library(knitr)
library(magrittr)
library(dplyr)
library(data.table)
library(ggplot2)
library(mikeutils)
library(doParallel)
library(foreach)
library(gifti)
library(fANCOVA)
library(viridis)
library(colorspace)
library(boot)
library(vegan)
library(grid)
library(gridExtra)
library(lme4)
library(HLMdiag)

source(here("code", "strings.R"))
source(here("code", "read_atlases.R"))

## data

stats.subjs.tdic <- fread(
  here("out", "rsa", "stats", paste0("subjs_pro_bias_acc-only_mmp_pearson_residual_glm-tdic.csv"))
)
stats.subjs.tdic <- stats.subjs.tdic[is.analysis.group == TRUE & y == "rank", ]  ## EXCLUDE HELD OUT SUBJECTS!
stats.subjs.tdic <- stats.subjs.tdic[, "coef" := NULL]
stats.subjs.tdic %<>% full_join(atlas.key$mmp, by = "roi")

## variables

subjs <- unique(stats.subjs.tdic$subj)
rois.medial <- c("8BM_R", "p32pr_R", "d32_R")
rois.lateral <- c("p9-46v_R", "IFSp_R", "FEF_R")
rois <- c(rois.medial, rois.lateral)

bounds <- 0.96
crit <- qnorm(1 - (1 - bounds) / 2)
n.resamples <- 1E4
colors.model <- c(incongruency = "#d95f02", target = "#1b9e77", distractor = "#7570b3")


## subset

d <- stats.subjs.tdic %>%
  filter(roi %in% c(rois.medial, rois.lateral), param %in% c("incongruency", "target")) %>%
  mutate(saggital = ifelse(roi %in% rois.medial, "medial", "lateral"))


## noise ceilings ----

nc <- data.frame(
  roi = c(rois.medial, rois.lateral),
  ub = NA,
  lb = NA
)

rsarray <- readRDS(
  here("out", "rsa", "obsv", "rsarray_pro_bias_acc-only_mmp_pearson_residual-linear.rds")
  )[, , unique(stats.subjs.tdic$subj), c(rois.medial, rois.lateral)]

## upper bound noise ceiling

z <- atanh(rsarray)
zbar <- apply(z, c(".row", ".col", "roi"), mean)

z.df <- apply(z, c("roi", "subj"), mat2vec)  ## returns matrix of lists
z.df <- bind_rows(apply(z.df, 1, function(.) bind_rows(., .id = "subj")), .id = "roi")  ## unwrap into df

zbar.df <- bind_rows(apply(zbar, "roi", mat2vec), .id = "roi")  ## all-subject average
zbar.df %<>% rename(all = value)

## split by roi and loop over:
z.dfw <- z.df %>% tidyr::spread(subj, value)
z.l <- split(z.dfw, z.dfw$roi)
zbar.l <- split(zbar.df, zbar.df$roi)

for (roi.i in seq_along(zbar.l)) {
  # roi.i = 1
  
  name.roi.i <- names(zbar.l[roi.i])
  xandy <- full_join(zbar.l[[roi.i]], z.l[[roi.i]], by = c(".row", ".col", "roi"))
  
  Y <- xandy[, subjs]
  x <- xandy[, "all"]
  
  # fit <- lm.fit(x = x, y = Y)
  # nc[nc$roi == name.roi.i, "ub"] <- mean(fit$coefficients)
  
  nc[nc$roi == name.roi.i, "ub"] <- tanh(mean(atanh(cor(x, Y))))
  
}

## lower-bound noise ceiling (leave one subject out)


for (roi.i in seq_along(rois)) {
  
  name.roi.i <- names(zbar.l[roi.i])
  nc.i <- numeric(length(subjs))
  
  for (subj.i in seq_along(subjs)) {
    # subj.i = 1
    
    x <- mat2vec(apply(z[, , -subj.i, roi.i], c(".row", ".col"), mean))
    y <- mat2vec(z[, , subj.i, roi.i])
    xandy <- full_join(x, y, by = c(".row", ".col"))
    
    nc.i[subj.i] <- cor(xandy$value.y, xandy$value.x)
    
  }
  
  nc[nc$roi == name.roi.i, "lb"] <- mean(nc.i)
  
}


## plot ----

## means and noise ceilings (barplot)

e <- d %>%
  group_by(roi, subj) %>%
  mutate(beta_jbar = mean(beta)) %>%
  ungroup %>%
  mutate(beta_bar = mean(beta)) %>%
  group_by(roi, param, saggital) %>%
  summarize(
    beta = mean(beta),
    ci = sqrt((var(beta - beta_jbar + beta_bar) * 1/2) / n()) * crit
  )

grid.arrange(
  d %>%
    ggplot(aes(roi, beta, color = param)) +
    facet_grid(cols = vars(saggital), scales = "free") +
    stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 1), fun.args = list(B = 1E4)),
  e %>%
    ggplot(aes(roi, beta, color = param)) +
    geom_point(position = position_dodge(width = 1)) +
    facet_grid(cols = vars(saggital), scales = "free") +
    # stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 1), fun.args = list(B = 1E4)) +
    # stat_summary(geom = "point", fun.y = "mean", position = position_dodge(width = 1)) +
    geom_errorbar(aes(ymin = beta - ci, ymax = beta + ci), position = position_dodge(width = 1))
)

nc4plot <- rbind(cbind(nc, param = "incon"), cbind(nc, param = "targt"))
nc4plot$saggital <- ifelse(nc4plot$roi %in% rois.lateral, "lateral", "medial")

as.numeric(nc4plot$roi)
nc4plot$xpos[nc4plot$roi == "FEF_R"]    <- 1
nc4plot$xpos[nc4plot$roi == "IFSp_R"]   <- 2
nc4plot$xpos[nc4plot$roi == "p9-46v_R"] <- 3
nc4plot$xpos[nc4plot$roi == "8BM_R"]    <- 1
nc4plot$xpos[nc4plot$roi == "d32_R"]    <- 2
nc4plot$xpos[nc4plot$roi == "p32pr_R"]  <- 3

p.dissoc <- d %>%
  mutate(roi = as.factor(roi)) %>%
  ggplot(aes(roi, beta, color = param)) +
  facet_grid(cols = vars(saggital), scales = "free") +
  stat_summary(
    fun.data = "mean_cl_boot",
    position = position_dodge(width = 0.25), fun.args = list(B = n.resamples),
    size = 0.5, color = "transparent"
  ) +
  geom_rect(
    data = nc4plot,
    mapping = aes(xmin = xpos - 1/4, xmax = xpos + 1/4, ymin = lb, ymax = ub),
    inherit.aes = FALSE,
    fill = "grey75",
    alpha = 0.3
  ) +
  stat_summary(
    fun.data = "mean_cl_boot",
    position = position_dodge(width = 1/3), fun.args = list(B = n.resamples),
    size = 0.5
    ) +
  scale_color_manual(values = c(incongruency = "#1b9e77", target = "#7570b3")) +
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(breaks = c(0, 0.1)) +
  geom_segment(
    data = data.frame(y = 0, yend = 0.1, x = -Inf, xend = -Inf, saggital = "lateral"),
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey50", size = 0.25, inherit.aes = FALSE
  )

ggsave(
  here("out", "figs", "ms_v1_2020-03", "group_dissoc", "group_dissoc_bar.pdf"),
  plot = p.dissoc,
  units = "cm",
  device = "pdf",
  height = 5,
  width = 11
)

## cross plot

dbar <- d %>%
  group_by(param, roi, saggital) %>%
  summarize(beta = mean(beta))

# p.dissoc.cross <- 
dbar %>%
  ggplot(aes(saggital, beta, color = param)) +
  geom_point(
    data = dbar %>% group_by(param, saggital) %>% summarize(beta = mean(beta)),
    size = 4
  ) +
  geom_line(
    data = dbar %>% group_by(param, saggital) %>% summarize(beta = mean(beta)),
    aes(group = param)
  ) +
  geom_point(aes(group = roi), alpha = 0.5, position = position_dodge(width = 0.5)) +
  # geom_errorbar(
  #   data = d %>%
  #     group_by(subj, roi) %>%
  #     mutate(
  #       beta_j = mean(beta)
  #     ) %>%
  #     group_by(roi) %>%
  #     mutate(
  #       beta_adj = beta - beta_j + mean(beta)
  #     ) %>%
  #     group_by(param, roi, saggital) %>%
  #     summarize(
  #       beta = mean(beta),
  #       se = sqrt(var(beta_adj) / length(unique(subj)) * 2)
  #     ),
  #   aes(ymin = beta - se * 1.96, ymax = beta + se * 1.96, width = 0.1, group = roi),
  #   position = position_dodge(width = 0.5)
  # ) +
  # geom_line(aes(group = param), alpha = 0.5, position = position_dodge(width = 0.5))
  scale_color_manual(values = colors.model) +
  # scale_color_manual(values = c(medial = "firebrick", lateral = "black"))
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(breaks = c(0, 0.1)) +
  geom_segment(
    data = data.frame(y = 0, yend = 0.1, x = -Inf, xend = -Inf, saggital = "lateral"),
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey50", size = 0.25, inherit.aes = FALSE
  )

ggsave(
  here("out", "figs", "ms_v1_2020-03", "group", "group_dissoc_cross.pdf"),
  plot = p.dissoc.cross,
  units = "cm",
  device = "pdf",
  height = 6,
  width = 8
)


## model ----

fit.all <- lmer(beta ~ param * saggital + (1 | subj), d)

## diagnostics

plot(fit.all, subj ~ resid(., type = "pearson"), main = "pearson's resids by subject")
plot(fit.all, param ~ resid(., type = "pearson"), main = "pearson's resids by subject")
plot(fit.all, interaction(roi, param) ~ resid(., type = "pearson"), main = "pearson's resids by subject")
ggplot_qqnorm(resid(fit.all), line = "rlm")

## summary
summary(fit.all)
fit.all.confint <- confint(fit.all, level = 0.96)



fit.lateral <- lmer(beta ~ param + (1 | subj), d, subset = saggital == "lateral")
summary(fit.lateral)

fit.medial <- lmer(beta ~ param + (1 | subj), d, subset = saggital == "medial")
summary(fit.medial)


## diagnostics

plot(fit.all, subj ~ resid(., type = "pearson"), main = "pearson's resids by subject")
plot(fit.all, param ~ resid(., type = "pearson"), main = "pearson's resids by subject")
plot(fit.all, interaction(roi, param) ~ resid(., type = "pearson"), main = "pearson's resids by subject")
ggplot_qqnorm(resid(fit.all), line = "rlm")

## summary
summary(fit.all)
fit.all.confint <- confint(fit.all, level = 0.96)
