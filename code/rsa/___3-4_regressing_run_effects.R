## about ----
## 

do.clusters <- TRUE
do.gordon <- FALSE


## set up env ----

file.suffix <- ifelse(do.clusters, "_multiparcel.rds", ".rds")
this.atlas <- ifelse(do.gordon, "gordon", "mmp")

## dependencies -- 

library(here)
library(dplyr)
library(data.table)
library(purrr)
source(here("..", "gen", "funs", "_get_dirs_local.R"))
source(here("..", "gen", "funs", "_funs.R"))
source(here("r", "group-201902", "_get_misc_vars.R"))

## objs --

## reprisimil matrix

rsarray <- readRDS(
  file.path(
    dir.box.stroop.data, 
    paste0("rsarray_pearson_", this.atlas, "_group-201902_pro_bias_acc-only", file.suffix)
    )
)

## event counts and run model
counts <- readxl::read_xlsx(file.path(dir.box.stroop.sheets, "summary_event-counts.xlsx"))
# View(counts)
counts <- counts[17:32, c("proactive1", "proactive2", "stimulus")]
counts$proactive1 <- counts$proactive1 > 3
counts$proactive2 <- counts$proactive2 > 3
run1.items <- counts$stimulus[counts$proactive1]
run2.items <- counts$stimulus[counts$proactive2]
run.rsm <- matrix(0, ncol = 16, nrow = 16, dimnames = list(bias.items, bias.items))
run.rsm[run1.items, run1.items] <- 1
run.rsm[run2.items, run2.items] <- 1
qcor(run.rsm, "run model", tl.cex = 0.5)  ## all looks good with schemes:

write.csv(run.rsm, file.path(dir.box.stroop.sheets, "rsmodel_run.csv"))


## regress ----

## prep -- 

## unwrap run model to lower-tri vector:
run.rsv <- mat2vec(run.rsm)
names(run.rsv) <- c("r", "c", "run.model")

## initialize indices and lists
subjs <- dimnames(rsarray)$subj
parcels <- dimnames(rsarray)$roi
# parcels <- unique(atlas.key$mmp$roi)  ## commented 2019-03-15 :: for multiparcel analysis (should work for both)
hemis <- c("l", "r")
residuals.linear <- vector("list", length(subjs) * length(parcels) * length(hemis))
names(residuals.linear) <- combo_paste3(subjs, parcels, hemis)
residuals.rank <- residuals.linear

## loop --

for (n.subj.i in seq_along(subjs)) {
  for (n.parcel.j in seq_along(parcels)) {
    for (n.hemi.k in seq_along(hemis)) {
      # n.hemi.k = 1; n.parcel.j = 1; n.subj.i = 1
      
      ## get rsm:
      rsm <- rsarray[, , n.subj.i, n.parcel.j, n.hemi.k]
      rsv <- mat2vec(rsm)
      rsv <- full_join(rsv, run.rsv, by = c("r", "c"))
      
      ## model:
      fit.linear <- lm(atanh(value) ~ run.model, rsv)
      fit.rank <- lm(rank(value) ~ run.model, rsv)
      resids <- cbind(
        rsv[c("r", "c")], 
        ## re-center (and transform back to r):
        residuals.linear = tanh(residuals(fit.linear) + coef(fit.linear)["(Intercept)"]),
        residuals.rank   = residuals(fit.rank) + coef(fit.rank)["(Intercept)"]
        )
      
      ## store:
      name.ijk <- paste(subjs[n.subj.i], parcels[n.parcel.j], hemis[n.hemi.k], sep = "_")
      residuals.linear[[name.ijk]] <- resids[c("r", "c", "residuals.linear")]
      residuals.rank[[name.ijk]] <- resids[c("r", "c", "residuals.rank")]
      
    }
  }
}

## wrap vectors back into array --

## sanity checking:
# test <- residuals.linear[[1]]
# a <- paste0(as.character(test$r), "_", as.character(test$c))
# b <- diag(16)
# colnames(b) <- bias.items
# rownames(b) <- bias.items
# b[lower.tri(b, diag=FALSE)] <- a
# b <- t(b)
# b[lower.tri(b, diag=FALSE)] <- a
# b
# test <- residuals.linear[[1]]
# a <- test$residuals.linear
# b <- diag(16)
# colnames(b) <- bias.items
# rownames(b) <- bias.items
# b[lower.tri(b, diag=FALSE)] <- a
# b <- t(b)
# b[lower.tri(b, diag=FALSE)] <- a
# qcor(b)
# isSymmetric.matrix(b)
# b
# test[20:40, ]

## prep --

empty.rsm <- diag(16)
colnames(empty.rsm) <- bias.items
rownames(empty.rsm) <- bias.items
rsarray.resid.rank <- array(NA, dim = dim(rsarray), dimnames = dimnames(rsarray))
rsarray.resid.linear <- rsarray.resid.rank

## loop --

for (n.subj.i in seq_along(subjs)) {
  for (n.parcel.j in seq_along(parcels)) {
    for (n.hemi.k in seq_along(hemis)) {
      # n.hemi.k = 1; n.parcel.j = 1; n.subj.i = 1
      
      name.ijk <- paste(subjs[n.subj.i], parcels[n.parcel.j], hemis[n.hemi.k], sep = "_")
      
      ## build matrices:
      
      rsm.rank.i <- empty.rsm
      
      rsm.rank.i[lower.tri(rsm.rank.i, diag = FALSE)] <- residuals.rank[[name.ijk]]$residuals.rank
      rsm.rank.i <- t(rsm.rank.i)
      rsm.rank.i[lower.tri(rsm.rank.i, diag = FALSE)] <- residuals.rank[[name.ijk]]$residuals.rank
      
      rsm.linear.i <- empty.rsm
      rsm.linear.i[lower.tri(rsm.linear.i, diag = FALSE)] <- residuals.linear[[name.ijk]]$residuals.linear
      rsm.linear.i <- t(rsm.linear.i)
      rsm.linear.i[lower.tri(rsm.linear.i, diag = FALSE)] <- residuals.linear[[name.ijk]]$residuals.linear
      
      ## store:
      rsarray.resid.rank[, , n.subj.i, n.parcel.j, n.hemi.k] <- rsm.rank.i
      rsarray.resid.linear[, , n.subj.i, n.parcel.j, n.hemi.k] <- rsm.linear.i

    }
  }
}


## save ----

saveRDS(
  rsarray.resid.rank, 
  file.path(
    dir.box.stroop.data, 
    paste0("rsarray_pearson_", this.atlas, "_group-201902_pro_bias_acc-only_rank-residual", file.suffix)
    )
  )
saveRDS(
  rsarray.resid.linear, 
  file.path(
    dir.box.stroop.data, 
    paste0("rsarray_pearson_", this.atlas, "_group-201902_pro_bias_acc-only_linear-residual", file.suffix)
    )
  )
