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

stats.group.tdic <- stats.subjs.tdic %>%
  
  filter(is.analysis.group) %>%
  
  group_by(param, roi) %>%
  summarize(p = wilcox.test(beta, alternative = "greater")$p.value) %>%
  
  mutate(p.fdr = p.adjust(p, method = "fdr"))

stats.group.tdic %>%
  filter(param == "target", p.fdr < 0.05)

## variables

subjs <- unique(stats.subjs.tdic$subj)
rois.medial <- c("8BM_R", "p32pr_R")
rois.lateral <- c("p9-46v_R", "FEF_R")
rois.visual <- c("V1_L", "V2_L")
rois.sommot <- c("3b_L", "4_L")
rois <- c(rois.medial, rois.lateral, rois.visual, rois.sommot)

bounds <- 0.95
crit <- qnorm(1 - (1 - bounds) / 2)
n.resamples <- 1E4
colors.model <- c(incongruency = "#d95f02", target = "#1b9e77", distractor = "#7570b3")


## subset

d <- stats.subjs.tdic %>% filter(roi %in% rois)


## noise ceilings ----

nc <- data.frame(
  roi = c(rois.medial, rois.lateral),
  ub = NA,
  lb = NA
)

rsarray <- readRDS(
  here("out", "rsa", "obsv", "rsarray_pro_bias_acc-only_mmp_pearson_residual-linear.rds")
)[, , unique(stats.subjs.tdic$subj), rois]

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
  filter(param != "congruency") 
  # group_by(roi, subj) %>%
  # mutate(beta_jbar = mean(beta)) %>%
  # ungroup %>%
  # mutate(beta_bar = mean(beta)) %>%
  # group_by(roi, param, community) %>%
  # summarize(
  #   beta = mean(beta),
  #   ci = sqrt((var(beta - beta_jbar + beta_bar) * 1/2) / n()) * crit
  # )

# grid.arrange(
#   d %>%
#     filter(param != "congruency") %>%
#     ggplot(aes(roi, beta, color = param)) +
#     facet_grid(cols = vars(community), scales = "free") +
#     stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 1), fun.args = list(B = 1E4)),
#   e %>%
#     ggplot(aes(roi, beta, color = param)) +
#     geom_point(position = position_dodge(width = 1)) +
#     facet_grid(cols = vars(community), scales = "free") +
#     # stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 1), fun.args = list(B = 1E4)) +
#     # stat_summary(geom = "point", fun.y = "mean", position = position_dodge(width = 1)) +
#     geom_errorbar(aes(ymin = beta - ci, ymax = beta + ci), position = position_dodge(width = 1))
# )
# 
# 
# d %>%
#   filter(param != "congruency") %>%
#   ggplot(aes(roi, beta, color = param)) +
#   facet_grid(cols = vars(community), scales = "free") +
#   stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 1), fun.args = list(B = 1E4))


unique(e$roi)
unique(e$region)
e$region <- ""
e$region[e$roi %in% rois.sommot] <- "somatomotor"
e$region[e$roi %in% rois.visual] <- "visual"
e$region[e$roi %in% rois.medial] <- "mFC"
e$region[e$roi %in% rois.lateral] <- "lFC"


ylims <- e %>%  
  group_by(param, roi) %>%
  summarize(b = mean(beta))

e %<>% mutate(region = factor(region, levels = c("visual", "somatomotor", "mFC", "lFC")))



p.tar <- e %>%  
  
  filter(param == "target") %>%
  
  ggplot(aes(roi, beta, color = param)) +
  # geom_point(size = 6) +
  # geom_errorbar(aes(ymin = beta - ci, ymax = beta + ci), width = 0, size = 2) +
  stat_summary(fun.data = "mean_cl_boot", size = 2.5) +
  
  # facet_wrap(vars(param), ncol = 1, scales = "free_y") +
  
  geom_segment(
    aes(
      x = Inf, xend = Inf, y = 0, 
      yend = ylims %>% filter(param == "target") %>% pull(b) %>% max %>% round(2),
    ),
    size = 4,
    color = colors.model["target"]
  ) +
  geom_segment(
    aes(
      x = -Inf, xend = Inf, y = -Inf, yend = -Inf
    ),
    size = 4,
    color = "black"
  ) +
  
  scale_color_manual(values = colors.model) +
  scale_y_continuous(
    breaks = c(0, ylims %>% filter(param == "target") %>% pull(b) %>% max %>% round(2)), position = "right"
  ) +
  # scale_x_discrete(
  #   labels = c("M1", "dmPFC/SMA", "dPM/FEF", "V1")
  # ) +
  
  facet_wrap(vars(region), nrow = 4, scales = "free_y", strip.position = "left") +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "bold", size = rel(2)),
    # axis.text.x = element_text(face = "bold", size = rel(4), angle = 15, hjust = 0.5, vjust = 0.5),
    axis.text.x = element_text(size = rel(2.5), color = colors.model["target"]),
    # axis.text.y = element_text(color = rev(labels$colors), face = "bold", size = rel(1)),
    axis.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = rel(2)),
    strip.placement = "outside",
    panel.border = element_blank(),
    # axis.line.x = element_line(size = 2),
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black", size = rel(2))
  ) +
  
  coord_flip()




p.dis <- e %>%
  
  filter(param == "distractor") %>%
  
  ggplot(aes(roi, beta, color = param)) +
  # geom_point(size = 6) +
  # geom_errorbar(aes(ymin = beta - ci, ymax = beta + ci), width = 0, size = 2) +
  stat_summary(fun.data = "mean_cl_boot", size = 2.5) +
  
  # facet_wrap(vars(param), ncol = 1, scales = "free_y") +
  
  geom_segment(
    aes(
      x = Inf, xend = Inf, y = 0, 
      yend = ylims %>% filter(param == "distractor") %>% pull(b) %>% max %>% round(2)
    ),
    size = 4,
    color = colors.model["distractor"]
  ) +
  
  scale_color_manual(values = colors.model) +
  scale_y_continuous(
    breaks = c(0, ylims %>% filter(param == "distractor") %>% pull(b) %>% max %>% round(2)), position = "right"
  ) +
  
  facet_wrap(vars(region), nrow = 4, scales = "free_y", strip.position = "left") +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = rel(2.5), color = colors.model["distractor"]),
    # axis.text.y = element_text(color = rev(labels$colors), face = "bold", size = rel(1)),
    axis.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  ) +
  
  coord_flip()
  



p.con <- e %>%
  
  filter(param == "incongruency") %>%
    
    ggplot(aes(roi, beta, color = param)) +
    # geom_point(size = 6) +
    # geom_errorbar(aes(ymin = beta - ci, ymax = beta + ci), width = 0, size = 2) +
    stat_summary(fun.data = "mean_cl_boot", size = 2.5) +
    
  # facet_wrap(vars(param), ncol = 1, scales = "free_y") +
    
    geom_segment(
      aes(
        x = Inf, xend = Inf, y = 0, 
        yend = ylims %>% filter(param == "incongruency") %>% pull(b) %>% max %>% round(2)
      ),
      size = 4,
      color = colors.model["incongruency"]
    ) +
    
    scale_color_manual(values = colors.model) +
    scale_y_continuous(
      breaks = c(0, ylims %>% filter(param == "incongruency") %>% pull(b) %>% max %>% round(2)), position = "right"
    ) +
    
    facet_wrap(vars(region), nrow = 4, scales = "free_y", strip.position = "left") +
    
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = rel(2.5), color = colors.model["incongruency"]),
      # axis.text.y = element_text(color = rev(labels$colors), face = "bold", size = rel(1)),
      axis.title = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank()
    ) +
    
    coord_flip()
  
ggsave(
  here("docs", "posters", "ohbm2020", "figs", "group_dissoc.pdf"),
  plot = arrangeGrob(p.tar, p.dis, p.con, ncol = 3, widths = c(2, 1.5, 1.5)),
  units = "cm",
  device = "pdf",
  height = 25,
  width = 30
)



## model ----


fit.lm.t <- lmer(beta ~ region + (1 | subj), e, subset = roi %in% c(rois.lateral, rois.medial) & param == "target")
summary(fit.lm.t)

fit.lm.c <- lmer(beta ~ region + (1 | subj), e, subset = roi %in% c(rois.lateral, rois.medial) & param == "incongruency")
summary(fit.lm.c)

fit.vm.t <- lmer(beta ~ region + (1 | subj), e, subset = roi %in% c(rois.sommot, rois.visual) & param == "target")
summary(fit.vm.t)

fit.vm.c <- lmer(beta ~ region + (1 | subj), e, subset = roi %in% c(rois.sommot, rois.visual) & param == "incongruency")
summary(fit.vm.c)

fit.vm.d <- lmer(beta ~ region + (1 | subj), e, subset = roi %in% c(rois.sommot, rois.visual) & param == "distractor")
summary(fit.vm.d)
