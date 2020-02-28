
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

source(here("code", "strings.R"))
source(here("code", "funs.R"))

theme_set(theme_bw(base_size = 12))

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


## (a) pca ----

d <- stats.subjs.tdic %>% filter(
  roi.set == "anatfunc", 
  superparcel %in% c("dlpfc_L", "dlpfc_R", "mfc_L", "mfc_R", "lppc_L", "lppc_R")
)


d %<>% mutate(id = as.character(interaction(superparcel, param)))

w <- d %>%
  ungroup %>%
  dplyr::select(subj, stroop, id, beta) %>%
  tidyr::spread(id, beta)
m <- w[, -c(1, 2)]
rownames(m) <- w$subj

pca <- prcomp(scale(m))

plot(pca)


p.pca <- pca$rotation[, 1:2] %>%
  as.data.frame %>%
  tibble::rownames_to_column("id") %>%
  mutate(
    roi = gsub("(.*_.).*", "\\1", id),
    param = gsub(".*_..(.*)", "\\1", id)
  ) %>%
  ggplot(aes(PC1, PC2, color = param)) +
  geom_point() +
  geom_text(aes(label = roi, color = param), nudge_y = 0.00, nudge_x = 0.02, fontface = "bold") +
  geom_segment(
    aes(x = 0, y = 0, xend = 0, yend = 0.1), 
    arrow = arrow(length = unit(0.01, "npc"), type = "closed"), color = "grey40", size = 1
  ) +
  geom_segment(
    aes(x = 0, y = 0, xend = 0.1, yend = 0),
    arrow = arrow(length = unit(0.01, "npc"), type = "closed"), color = "grey40", size = 1
  ) +
  scale_color_manual(values = colors.model) +
  annotate(geom = "text", label = "PC1", x = 0.05, y = 0, vjust = 1, fontface = "bold.italic", color = "grey40") +
  annotate(geom = "text", label = "PC2", x = 0, y = 0.05, angle = 90, vjust = 0, fontface = "bold.italic", color = "grey40") +
  annotate(geom = "text", label = "0.1", x = -0.0025, y = 0.1, hjust = 1) +
  annotate(geom = "text", label = "0.1", y = -0.005, x = 0.1, vjust = 1) +
  annotate(geom = "text", label = "0", y = -0.005, x = -0.0025, vjust = 1) +
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.ticks = element_blank(), 
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.border = element_blank()
  ) +
  annotate(geom = "text", y = 0.24, x = 0.05, color = colors.model["target"], fontface = "bold.italic", label = "target coding") +
  annotate(geom = "text", y = 0.22, x = 0.05, color = colors.model["distractor"], fontface = "bold.italic", label = "distractor coding") +
  annotate(geom = "text", y = 0.20, x = 0.05, color = colors.model["incongruency"], fontface = "bold.italic", label = "conflict coding") +
  labs(title = "a")

p.pca


## (b) scatterplot


cor(w$stroop, (w$lppc_R.target + w$dlpfc_R.target) / 2)
cor(w$stroop, (w$lppc_R.target + w$dlpfc_R.target) / 2, method = "spearman")
w$lfp_R.target <- (w$lppc_R.target + w$dlpfc_R.target) / 2
cor(w$stroop, w$mfc_L.incongruency)
cor(w$mfc_R.incongruency, w$stroop)

p.dissoc.scatter <- w %>%
  ggplot(aes(y = stroop)) +
  stat_boot_ci(aes(x = lfp_R.target), alpha = 0.3, n = 1E3, percent = 96, fill = colors.model["target"]) +
  stat_boot_ci(aes(x = mfc_L.incongruency), alpha = 0.3, n = 1E4, percent = 96, fill = colors.model["incongruency"]) +
  geom_smooth(aes(x = lfp_R.target), alpha = 0.3, color = colors.model["target"], method = "lm", se = FALSE) +
  geom_smooth(aes(x = mfc_L.incongruency), alpha = 0.3, color = colors.model["incongruency"], method = "lm", se = FALSE) +
  geom_point(aes(x = lfp_R.target), shape = 21, fill = colors.model["target"], color = "white", size = 3) +
  geom_point(aes(x = mfc_L.incongruency), shape = 21, fill = colors.model["incongruency"], color = "white", size = 3) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_blank(),
    panel.border =  element_blank(),
    panel.background = element_blank(),
    # axis.text.y = element_blank(),
    # axis.title.y = element_blank(),
    # axis.ticks.y = element_blank()
  ) +
  scale_x_continuous(breaks = c(min(w$mfc_L.incongruency), 0, max(w$mfc_L.incongruency)) %>% round(2)) +
  scale_y_continuous(breaks = c(min(w$stroop), max(w$stroop)) %>% round(2)) +
  geom_segment(
    aes(y = round(min(w$stroop), 2), yend = round(max(w$stroop), 2), x = -Inf, xend = -Inf),
    color = "grey40", size = 1
    ) +
  geom_segment(
    aes(x = round(min(w$mfc_L.incongruency), 2), xend = round(max(w$mfc_L.incongruency), 2), y = -Inf, yend = -Inf),
    color = "grey40", size = 1
  ) +
  labs(title = "b", x = "RSA model fit")



## (b) trichotomized MDS

## put together and write ----

p.dissoc <- grid.arrange(
  p.pca, p.dissoc.scatter, ncol = 2, widths = c(1, 1)
)

ggsave(
  here("out", "figs", "ms_v1_2020-03", "indiff_dissoc", "indif_dissoc.pdf"), 
  plot = p.dissoc,
  units = "cm",
  device = "pdf",
  height = 10,
  width = 20# * 1.75
)
