library(here)
library(dplyr)
library(abind)

stats.conclust <- read.csv(here("out", "rsa", "stats", "subjs_pro_bias_acc-only_mmp_pearson_residual_glm-tdclust.csv"))

stats.orig <- read.csv(here("old", "wustl_proj_stroop-rsa", "stats", "subject-coefs_pcor_rsm-pearson_mmp_group-201902_pro_bias_acc-only_residual.csv"))

stats.conclust %>% head
stats.orig %>% head


stats.orig <- stats.orig %>% select(-r) %>% mutate(roi = paste0(parcel, "_", toupper(hemi))) %>% select(-parcel, -hemi)
stats.conclust %<>% filter(y == "rank") %>% 
  select(subj, roi, param, partr)
stats.orig <- stats.orig %>% rename(param = variable, partr = rho)

full_join(stats.orig, stats.conclust)
stats.conclust <- stats.conclust %>% filter(subj %in% unique(stats.orig$subj))

all.equal(
  stats.conclust %>% arrange(subj, param, roi) %>% select(partr),
  stats.orig %>% arrange(subj, param, roi) %>% select(partr)
)


results.orig <- stats.orig %>%
  group_by(roi, param) %>%
  summarize(
    rho.srtest.v = wilcox.test(partr, alternative = "greater")$statistic,
    rho.srtest.p = wilcox.test(partr, alternative = "greater")$p.value,
  ) %>%
  dplyr::group_by(param) %>%
  dplyr::mutate(
    p.adj.wb = p.adjust(rho.srtest.p, method = "fdr"),
  )

results.orig %>% filter(param == "target", p.adj.wb < 0.05) %>% nrow



results.conclust <- stats.conclust %>%
  group_by(roi, param) %>%
  summarize(
    rho.srtest.v = wilcox.test(partr, alternative = "greater")$statistic,
    rho.srtest.p = wilcox.test(partr, alternative = "greater")$p.value,
  ) %>%
  dplyr::group_by(param) %>%
  dplyr::mutate(
    p.adj.wb = p.adjust(rho.srtest.p, method = "fdr"),
  )


results.conclust %>% filter(param == "conclust", p.adj.wb < 0.05) %>% nrow


stats.group.tdic <- stats.conclust %>%
  group_by(roi, param) %>%
  summarize(
    v    = wilcox.test(partr, alternative = "greater")$statistic,
    p    = wilcox.test(partr, alternative = "greater")$p.value
  ) %>%
  # group_by(param) %>%
  mutate(p.fdr = p.adjust(p, method = "fdr"), p.holm = p.adjust(p, method = "holm"))

stats.group.tdic %>% filter(param == "target", p.fdr < 0.05) %>% nrow





partr.rank[rownames(partr.rank) == "132017_V1_L", ]

x <- rsvectors[["132017_V1_L"]]

rsv <- x[, "rank"]
d <- cbind(rsv, X[, -1])
psych::partial.r(d, method = "spearman")[1, -1]


partr.rank <- do.call(rbind, lapply(rsvectors, get.partr, "rank"))

















# ----

## are the RDMs identical?

rdm.rank <- readRDS(
  here("out", "rsa", "obsv", "rsarray_pro_bias_acc-only_mmp_pearson_residual-rank.rds")
)
rdm.old <- readRDS(
  here("old", "wustl_proj_stroop-rsa", "data", "rsarray_pearson_mmp_group-201902_pro_bias_acc-only_rank-residual.rds")
)
here("old", "wustl_proj_stroop-rsa", "data") %>% list.files
rdm.rank %>% dimnames
rdm.old %>% dimnames


rdm.old <- rdm.old[, , dimnames(rdm.rank)$subj, , ]
dimnames.rdm.old <- dimnames(rdm.old)

dimnames.rdm.old$roi <- mikeutils::combo_paste(dimnames.rdm.old$roi, c("L", "R"))
dimnames.rdm.old$hemi <- NULL

rdm.old <- array(c(rdm.old), c(16, 16, 66, 180 * 2))
dimnames(rdm.old) <- dimnames.rdm.old

identical(c(rdm.old), c(rdm.rank))

## yes


## ----

## new

library(data.table)

rsv.models.ltri <- fread(here("out", "rsa", "mods", "rsv_bias_lower-triangles.csv"), data.table = FALSE)
variables <- as.matrix(rsv.models.ltri[sapply(rsv.models.ltri, is.numeric)])
variables <- cbind(variables, conclust = variables[, "incongruency"] | variables[, "congruency"])

structures <- list(
  tdic    = c("target", "distractor", "incongruency", "congruency"),
  tdi     = c("target", "distractor", "incongruency"),
  tdclust = c("target", "distractor", "conclust"),
  all     = c("target", "cielab", "tphono", "distractor", "silhou", "orthog", "dphono", "incongruency", "congruency"),
  continu = c("target", "cielab", "tphono", "distractor", "silhou", "dphono", "incongruency", "congruency")
)

m <- list(
  tdic    = cbind(b0 = 1, variables[, structures$tdic]),
  tdi     = cbind(b0 = 1, variables[, structures$tdi]),
  tdclust = cbind(b0 = 1, variables[, structures$tdclust]),
  all     = cbind(b0 = 1, variables[, structures$all]),  ## continuous and categorical
  continu = cbind(b0 = 1, variables[, structures$continu])  ## perceptual and phonological for target and distractor
)

set.i <- "mmp"
subj.i <- "132017"
roi.j <- "V1_L"


rsarray.rank <- readRDS(
  here(
    "out", "rsa", "obsv",
    paste0("rsarray_pro_bias_acc-only_", set.i, "_pearson_residual-rank.rds")
  )
)

rsm.rank <- rsarray.rank[, , subj.i, roi.j]
rsv.rank <- rsm.rank[lower.tri(rsm.rank)]

X <- variables[, structures$tdclust]
d2 <- cbind(rsv.rank, X)
psych::partial.r(d2, method = "spearman")[1, -1]


## original

fname.i <- paste0(
  "rsarray_pearson_", set.i, "_group-201902_pro_bias_acc-only_rank-residual.rds"
)
rsarray <- readRDS(here("old", "wustl_proj_stroop-rsa", "data", fname.i))

slice.ijk <- rsarray[, , subj.i, "V1", "l"]
slice.ijk[upper.tri(diag(16), diag = TRUE)] <- NA
vector.ijk <- reshape2::melt(slice.ijk, na.rm = TRUE) 


row.info <- split.str.item(vector.ijk$r)  ## see _get_misc_vars.R
col.info <- split.str.item(vector.ijk$c)
x <- cbind(
  vector.ijk,
  target      = as.numeric(row.info$color == col.info$color),
  distractor  = as.numeric(row.info$word == col.info$word),
  congruency  = as.numeric(row.info$congruency == col.info$congruency)
  # incongruency = as.numeric(row.info$congruency == col.info$congruency & col.info$congruency == "I")
)
library(mikeutils)
is.local.session <- TRUE
source(here("old", "R01_freund_stroop-rsa", "stroop-rsa", "r", "group-201902", "_get_misc_vars.R"))


psych::partial.r(x[c("value", "target", "distractor", "congruency")], method = "spearman")["value", c("target", "distractor", "congruency")]

