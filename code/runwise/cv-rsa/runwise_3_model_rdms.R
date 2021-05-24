## about ----
## 
## fits general linear models on each subject's RSMs from each parcel then writes the results to files.
## 
## mike freund, created 2019-02-24, updated 2019-03-05
## adapted for new project directory 2019-12-29


## setup ----

source(here::here("code", "packages.R"))
source(here("code", "strings.R"))
source(here("code", "read_atlases.R"))
source(here("code", "read_masks.R"))


## analysis group:

stroop <- fread(here("in", "behavior-and-events_group201902.csv"))
sample.analysis <- unique(filter(stroop, is.analysis.group)$subj)

n.dim <- length(bias.items)
is.lower.tri <- lower.tri(diag(n.dim))
# subjs <- unique(stroop$subj)


## read models, convert to distances

rsv.models.ltri <- fread(here("out", "rsa", "mods", "rsv_bias_lower-triangles.csv"), data.table = FALSE)
rsv.models.ltri[c("target", "distractor", "congruency", "incongruency")] <- 
  1 - rsv.models.ltri[c("target", "distractor", "congruency", "incongruency")]

## create design matrices

X <- scale(as.matrix(rsv.models.ltri[sapply(rsv.models.ltri, is.numeric)]))


## loop over sets of ROIs ----

for (set.i in sets.of.rois) {
  # set.i = "masks"
  
  ## read observed similarity matrices (arrays)
  
  rsarray <- readRDS(
    here(
      "out", "rsa", "obsv",
      paste0("rsarray_pro_bias_acc-only_fmriprep_runwise_cv-euclidean_", set.i, ".rds")
    )
  )
  # sds <- apply(rsarray, c("subj", "roi"), function(x) sd(x[lower.tri(x)]))
  # rsarray <- sweep(rsarray, MARGIN = 3:4, STATS = sds, FUN = "/")  ## scale distances by off-diagonal standard dev.
  

  subjs <- dimnames(rsarray)$subj
  
  
  ## prepare similarity matrices for regression ----
  
  ## check if rows and col names are equal (should be, but just to be sure...)
  
  are.rowcol.equal <- isTRUE(all.equal(dimnames(rsarray)[[1]], dimnames(rsarray)[[2]]))
  if(!are.rowcol.equal) stop("rsarray isn't rsarray!")
  
  ## get indices and values
  
  rois <- dimnames(rsarray)$roi
  n.mods <- length(subjs) * length(rois)
  
  ## unwrap into lower-triangle vector
  
  rsvectors <- vector("list", n.mods)
  names(rsvectors) <- combo_paste(subjs, rois)
  
  for (subj.i in seq_along(subjs)) {
    for (roi.j in seq_along(rois)) {
      # subj.i = 1; roi.j = 1
      
      name.ij <- paste0(subjs[subj.i], "_", rois[roi.j])  ## to match name
      rsvectors[[name.ij]] <- scale(rsarray[, , subj.i, roi.j][is.lower.tri])  ## z-score
      
    }
  }
  
  ## check numbers
  
  lengths.rsvectors <- map_dbl(rsvectors, length)
  if(sum(lengths.rsvectors != 120) > 0) stop("missing row somewhere in rsvectors!")
  
  
  ## fit glms ----
  
  fits <- rsvectors %>% map(~ .lm.fit( x = X, y = .))
  betas <- as.data.frame(do.call(rbind, lapply(fits, coef)))
  
  names(betas) <- colnames(X)
  betas$id <- rownames(betas)
  betas <- melt(as.data.table(betas), id.vars = "id", value.name = "beta", variable.name = "param")
  
  
  ## format ----
  
  ## create subj, roi, and hemi cols from id col
  
  stats.subjs <- bind_cols(
    betas,
    reshape2::colsplit(betas$id, pattern = "_", names = c("subj", "roi"))
  )
  
  stats.subjs$is.analysis.group <- stats.subjs$subj %in% sample.analysis  ## add is.analysis.group col
  
  ## rearrange cols (and drop id col)
  
  stats.subjs %<>% select(subj, is.analysis.group, roi, param, beta)
  
  ## write ----
  
  fwrite(
    stats.subjs,
    here("out", "rsa", "stats",  paste0("subjs_pro_bias_acc-only_fmriprep_runwise_cv-euclidean_", set.i, ".csv"))
  )
  
  
}

