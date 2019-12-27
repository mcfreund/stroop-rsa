## about ----
## 


file.suffix <- ifelse(do.clusters, "_multiparcel.rds", ".rds")


## setup

library(here)
library(mikeutils)
library(dplyr)
library(data.table)
library(purrr)
source(here("code", "_strings.R"))
source(here("code", "_get_atlas.R"))

## reprisimil matrix

rsarray <- readRDS(
  here(
    "out", "rsa", "obsv", 
    paste0("rsarray_", glm.name, "_", atlas.i, "_pearson", file.suffix, ".rds")
    )
)

run.rsm <- as.matrix(read.csv(here("data", "bias_run.csv"), row.names = 1))

## regress ----

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

## loop

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
  here(
    "out", "rsa", "obsv", 
    paste0("rsarray_", glm.name, "_", atlas.i, "_pearson_residual-rank", file.suffix, ".rds")
    )
  )
saveRDS(
  rsarray.resid.linear, 
  here(
    "out", "rsa", "obsv", 
    paste0("rsarray_", glm.name, "_", atlas.i, "_pearson_residual-linear", file.suffix, ".rds")
    )
  )
