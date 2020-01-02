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

empty.rsm <- matrix(0, ncol = 16, nrow = 16, dimnames = list(bias.items, bias.items))
target.rsm <- empty.rsm
distractor.rsm <- empty.rsm
congruency.rsm <- empty.rsm
incongruency.rsm <- empty.rsm
for (color.i in bias.colors) target.rsm[grepl(color.i, bias.items), grepl(color.i, bias.items)] <- 1
for (word.i in bias.words) distractor.rsm[grepl(word.i, bias.items), grepl(word.i, bias.items)] <- 1
bias.items.congruency <- c(
  TRUE, rep(FALSE, 4),
  TRUE, rep(FALSE, 4),
  TRUE, rep(FALSE, 4),
  TRUE
)  ## in same order as bias.items
congruency.rsm[bias.items.congruency, bias.items.congruency] <- 1
congruency.rsm[!bias.items.congruency, !bias.items.congruency] <- 1
incongruency.rsm[!bias.items.congruency, !bias.items.congruency] <- 1

## check ----

qcor(incongruency.rsm)
qcor(congruency.rsm)
qcor(target.rsm) ## NB: sorted by distractor
qcor(distractor.rsm)  ## NB: sorted by distractor

## melt to data.frame

models.wide <- cbind(
  target.rsm %>% reshape2::melt(value.name = "target") %>% rename(x = Var1, y = Var2),
  distractor = c(distractor.rsm),
  congruency = c(congruency.rsm),
  incongruency = c(incongruency.rsm)
)
models.wide %<>%
  mutate(
    x = factor(x, levels = sort(bias.items)),
    y = factor(y, levels = rev(sort(bias.items)))
  )

## write ----

write.csv(target.rsm, here("out", "rsa", "mods", "bias_target.csv"))
write.csv(distractor.rsm, here("out", "rsa", "mods", "bias_distractor.csv"))
write.csv(congruency.rsm, here("out", "rsa", "mods", "bias_congruency.csv"))
write.csv(incongruency.rsm, here("out", "rsa", "mods", "bias_incongruency.csv"))
write.csv(models.wide, here("out", "rsa", "mods", "bias_full_matrix.csv"), row.names = FALSE)

## run model ----

counts <- data.table::fread(here("data", "summary_event-counts.csv"))
# View(counts)
counts <- counts[counts$stimulus %in% bias.items, c("proactive1", "proactive2", "stimulus")]
counts$proactive1 <- counts$proactive1 > 3
counts$proactive2 <- counts$proactive2 > 3
run1.items <- counts$stimulus[counts$proactive1]
run2.items <- counts$stimulus[counts$proactive2]
run.rsm <- matrix(0, ncol = 16, nrow = 16, dimnames = list(bias.items, bias.items))
run.rsm[run1.items, run1.items] <- 1
run.rsm[run2.items, run2.items] <- 1
qcor(run.rsm, "run model", tl.cex = 0.5)  ## all looks good with schemes:

write.csv(run.rsm, here("out", "rsa", "mods", "bias_run.csv"))
