## about ----
## creates model similarity matrices (for RSA) and writes them to .csv files.
## 
## mike freund, 2019-02-20
## adapted for new project directory 2019-12-24

## setup ----

library(here)
library(mikeutils)
library(magrittr)
library(dplyr)

source(here("code", "_strings.R"))


## create ----

rsm.empty <- matrix(0, ncol = 16, nrow = 16, dimnames = list(bias.items, bias.items))
rsm.target <- rsm.empty
rsm.distractor <- rsm.empty
rsm.congruency <- rsm.empty
rsm.incongruency <- rsm.empty
for (color.i in bias.colors) rsm.target[grepl(color.i, bias.items), grepl(color.i, bias.items)] <- 1
for (word.i in bias.words) rsm.distractor[grepl(word.i, bias.items), grepl(word.i, bias.items)] <- 1
bias.items.congruency <- c(
  TRUE, rep(FALSE, 4),
  TRUE, rep(FALSE, 4),
  TRUE, rep(FALSE, 4),
  TRUE
)  ## in same order as bias.items
rsm.congruency[bias.items.congruency, bias.items.congruency] <- 1
rsm.congruency[!bias.items.congruency, !bias.items.congruency] <- 1
rsm.incongruency[!bias.items.congruency, !bias.items.congruency] <- 1

## check ----

qcor(rsm.incongruency)
qcor(rsm.congruency)
qcor(rsm.target) ## NB: sorted by distractor
qcor(rsm.distractor)  ## NB: sorted by distractor

## melt to data.frame

rsv.models.full <- cbind(
  rsm.target %>% mat2vec(full.matrix = TRUE, value.name = "target"),
  distractor   = c(rsm.distractor),
  congruency   = c(rsm.congruency),
  incongruency = c(rsm.incongruency)
)
# rsmodels.full %<>%
#   mutate(
#     .row = factor(.row, levels = sort(bias.items)),
#     .col = factor(.col, levels = rev(sort(bias.items)))
#   )

rsv.models.ltri <- cbind(
  rsm.target %>% mat2vec(value.name = "target"),
  distractor   = rsm.distractor[lower.tri(rsm.distractor)],
  congruency   = rsm.congruency[lower.tri(rsm.congruency)],
  incongruency = rsm.incongruency[lower.tri(rsm.incongruency)]
)
# rsmodels.ltri %<>%
#   mutate(
#     x = factor(x, levels = sort(bias.items)),
#     y = factor(y, levels = rev(sort(bias.items)))
#   )


## write ----

write.csv(rsm.target, here("out", "rsa", "mods", "rsm_bias_target.csv"))
write.csv(rsm.distractor, here("out", "rsa", "mods", "rsm_bias_distractor.csv"))
write.csv(rsm.congruency, here("out", "rsa", "mods", "rsm_bias_congruency.csv"))
write.csv(rsm.incongruency, here("out", "rsa", "mods", "rsm_bias_incongruency.csv"))
write.csv(rsv.models.full, here("out", "rsa", "mods", "rsv_bias_full-matrices.csv"), row.names = FALSE)
write.csv(rsv.models.ltri, here("out", "rsa", "mods", "rsv_bias_lower-triangles.csv"), row.names = FALSE)

## run model ----

counts <- data.table::fread(here("data", "summary_event-counts.csv"))
# View(counts)
counts <- counts[counts$stimulus %in% bias.items, c("proactive1", "proactive2", "stimulus")]
counts$proactive1 <- counts$proactive1 > 3
counts$proactive2 <- counts$proactive2 > 3
run1.items <- counts$stimulus[counts$proactive1]
run2.items <- counts$stimulus[counts$proactive2]
rsm.run <- matrix(0, ncol = 16, nrow = 16, dimnames = list(bias.items, bias.items))
rsm.run[run1.items, run1.items] <- 1
rsm.run[run2.items, run2.items] <- 1
qcor(rsm.run, "run model", tl.cex = 0.5)  ## all looks good with schemes:

write.csv(rsm.run, here("out", "rsa", "mods", "rsm_bias_run.csv"))

