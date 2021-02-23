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


glm.name <- "pro_bias_acc-only_fmriprep_im_run_downsampled"
type <- c("crun", "wnrun", "both")

## analysis group:

stroop <- fread(here("in", "behavior-and-events_group201902.csv"))
sample.analysis <- unique(filter(stroop, is.analysis.group)$subj)

n.dim <- length(bias.items)

# subjs <- unique(stroop$subj)


## read models

rsv.models <- fread(here("out", "rsa", "mods", "rsv_bias_full-matrices.csv"), data.table = FALSE)

## create design matrices

rsv.models <- as.matrix(rsv.models[sapply(rsv.models, is.numeric)])
rsv.models <- cbind(rsv.models, diagonal = c(diag(length(bias.items))))
X <- scale(rsv.models)


## loop over sets of ROIs ----

for (set.i in sets.of.rois) {
  # set.i = "masks"
  
  ## read observed similarity matrices (arrays)
  
  rsarray <- readRDS(
    here(
      "out", "rsa", "obsv",
      paste0("rsarray_", glm.name, "_", set.i, ".rds")
    )
  )
  subjs <- dimnames(rsarray)$subj
  rois <- dimnames(rsarray)$roi
  n.mods <- length(subjs) * length(rois)
  rsvectors <- vector("list", n.mods)
  names(rsvectors) <- combo_paste(subjs, rois)
  are.rowcol.equal <- isTRUE(all.equal(dimnames(rsarray)[[1]], dimnames(rsarray)[[2]]))
  if(!are.rowcol.equal) stop("rsarray isn't rsarray!")
  
  
  for (type.i in type) {
    # type.i = "wnrun"
    
    if (type.i == "wnrun") {
      rsarray.i <- tanh((atanh(rsarray[, , , , "wnr1"]) + atanh(rsarray[, , , , "wnr2"]))/2)
    } else {
      rsarray.i <- rsarray[, , , , type.i]
    }
  
    for (subj.i in seq_along(subjs)) {
      for (roi.j in seq_along(rois)) {
        # subj.i = 1; roi.j = 1
        
        name.ij <- paste0(subjs[subj.i], "_", rois[roi.j])  ## to match name
        rsvectors[[name.ij]] <- scale(rank(c(rsarray.i[, , subj.i, roi.j])))  ## z-score
        
      }
    }
    
    ## check numbers
    
    lengths.rsvectors <- map_dbl(rsvectors, length)
    if(sum(lengths.rsvectors != 256) > 0) stop("missing row somewhere in rsvectors!")
    
    
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
      here("out", "rsa", "stats",  paste0("subjs_", glm.name, "_", type.i, "_", set.i, ".csv"))
    )
    
      
  }
  
  
  
  
}

