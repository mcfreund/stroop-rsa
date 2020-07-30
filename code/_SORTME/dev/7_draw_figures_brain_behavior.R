## about ----
## 

## setup ----

set.seed(0)

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

theme_set(theme_bw(base_size = 12))

## read data

blups <- read.csv(here("out", "behav", "stroop_blups_rt_group201902.csv"), stringsAsFactors = FALSE)
blups$incon <- blups$congr + blups$stroop

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


d <- stats.subjs.tdic %>% 
  filter(roi.set == "anatfunc" | superparcel %in% c("vwfa", "smmouth"))


## plot ----

d <- d %>%
  select(subj, congr, incon, superparcel, param, beta) %>%
  as.data.table %>%
  melt(value.name = "rt", id.vars = c("subj", "superparcel", "param", "beta"), measure.vars = c("congr", "incon"))

d1 <- d %>% filter(param == "distractor", superparcel == "dlpfc_L")
d1 %>%
  ggplot(aes(beta, rt, color = variable, fill = variable)) +
  geom_smooth(
    data = d1 %>% filter(variable == "congr"), method = "lm", se = FALSE
  ) +
  geom_smooth(
    data = d1 %>% filter(variable == "incon"), method = "lm", se = FALSE
  ) +
  stat_boot_ci(
    data = d1 %>% filter(variable == "incon"), alpha = 0.2
  ) +
  stat_boot_ci(
    data = d1 %>% filter(variable == "congr"), alpha = 0.2
  ) +
  geom_point() +
  geom_line(aes(group = subj), color = "grey20")