## about ----
## mike freund, 2019-02-23
##
## first reads in RData files from box (using box drive) which contains an array of rsms, with dims: row, column, 
## subject, parcel, and hemi. 
## (optionally reads lists of beta values, portioned by parcel, hemisphere, and subject.)
## (these RData files were written to the box folder via boxr::box_save() from ccplinux1 by the script 3_save_betas.R.)
## next, saves these RData files as RDS (which support assignment upon loading into R), and deletes RData files.
## finally, creates directory structure for data, and writes rsms (and opitionally, betas) to csv within directory 
## structure.
## this script should be run locally (because it requires box drive)
##
## why save as RDS opposed to csv or RData? 
## though csvs are preferred for reproducibility (supporting reading into software other than R), using RDS or RData
## files are considerably faster than loading nsubj*nparcel files and compiling to an array.
## RDS files are for single objects; RData files are for any number of ojbects.
## reading RDS files additionally allow the newly read object to be assigned to a name; RData files do not natively
## allow this.
## thus, the analyses in this project will rely upon RDS files.
## (Note that RDS files could not be uploaded via boxr / ccplinux1, hence this file.)


## TODO:
## 1. load .rdata
## 2. write RDS
## 3. remove .rdata
## 4. (in time) write csvs. ---- don't write betas!!!!!



## set up env. ----


# do.held.out <- TRUE
do.clusters <- TRUE
# do.combine.samps <- TRUE  ## NB: use only when certain is needed!!

## dependencies

## NB: BOX DRIVE

library(here)
library(reshape2)
library(data.table)
source(here("..", "gen", "funs", "_get_dirs_local.R"))
source(here("..", "gen", "funs", "_funs.R"))
source(here("r", "group-201902", "_get_misc_vars.R"))
should.write.csvs <- FALSE

## strings and funs

# data.type <- c("betalist", "rsarray")
data.type <- "rsarray"
atlases <- c("mmp", "gordon")
stats <- c("pearson", "euclidean")
group <- "group-201902"
method <- "pro_bias_acc-only"  ## glm name
file.suffix <- ".RData"
if (do.clusters) file.suffix <- paste0("_multiparcel", file.suffix)
# if (do.held.out) file.suffix <- paste0("_held-out", file.suffix)
fnames <- combo_paste3(data.type, stats, paste0(atlases, "_", group, "_", method, file.suffix))

## for use with rdata files.
## preferably, the data files would be saved as RDS (to enable base R functioning to load and return the 
## file: readRDS()).
## but, boxr does not seem to be able to write anything but an RData file. 
load.rdata <- function(.fname){
  ## stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file
  ## loads an RData file, and returns it
  load(.fname)
  get(ls()[ls() != ".fname"])
}

for (n.fname.i in seq_along(fnames)) {
  # n.fname.i = 1
  fname.rdata.i <- file.path(dir.box.stroop.data, fnames[n.fname.i])
  if (!file.exists(fname.rdata.i)) next
  rsarray <- load.rdata(fname.rdata.i)
  fname.rds.i <- gsub("RData", "rds", fname.rdata.i)
  saveRDS(rsarray, gsub("RData", "rds", fname.rdata.i))
  if (file.exists(fname.rds.i)) file.remove(fname.rdata.i)
}


# if (do.combine.samps) {
#   fnames.rds <- list.files(dir.box.stroop.data, pattern = "\\.rds")
#   fnames.split <- reshape2::colsplit(fnames.rds, "_", c("type", "statistic", "atlas", "rest"))
#   fnames.split$rest <- gsub("_held-out|_multiparcel", "", fnames.split$rest)
#   fnames.split$is.held.out <- grepl("held-out", fnames.rds)
#   fnames.split$is.mutliparcel <- grepl("multiparcel", fnames.rds)
#   fnames.split$is.residual <- grepl("residual", fnames.rds)
#   for (stat.i in stats) {
#     for (atlas.i in atlases) {
#       # stat.i = "pearson"; atlas.i = "mmp"
#       for (do.multiparcel.i in c(TRUE, FALSE)) {
#         do.multiparcel.i = TRUE
#         if ((stat.i == "gordon" | atlas.i == "euclidean") & do.multiparcel) next  ## these don't have a multiparcel yet
#         is.row <- fnames.split$statistic == stat.i & fnames.split$atlas == atlas.i & 
#           as.logical(do.multiparcel.i * fnames.split$is.mutliparcel) & !fnames.split$is.residual
#         if (sum(is.row) > 2) stop("sum(is.row) > 2")
#         
#         fnames.rds[is.row & !fnames.split$is.held.out]
#         rds.held.out <- readRDS(file.path(dir.box.stroop.data, fnames.rds[is.row & fnames.split$is.held.out]))
#         rds.analysis <- readRDS(file.path(dir.box.stroop.data, fnames.rds[is.row & !fnames.split$is.held.out]))
#         all.good.held.out <- identical(dimnames(rds.held.out)$subj, dimnames(rds.held.out)$subj)
#         all.good.analysis <- identical(dimnames(rds.analysis)$subj, dimnames(rds.analysis)$subj)
#         if (any(!all.good.analysis, !all.good.held.out)) stop("not all good")
#         rds.combined <- abind::abind(rds.held.out, rds.analysis, along = 3)
#         saveRDS(rds.combined, file.path(dir.box.stroop.data, fnames.rds[is.row & !fnames.split$is.held.out]))
#         file.remove(file.path(dir.box.stroop.data, fnames.rds[is.row & fnames.split$is.held.out]))
#         print(file.path(dir.box.stroop.data, fnames.rds[is.row & !fnames.split$is.held.out]))
#       }
#     }
#   }
# }

## read data ---
# 
# 
# for (atlas.i in seq_along(atlases)) {  ## loop over first-lev (atlases)
#   # atlas.i = 2
#     
#   ## save rsarray.RData as .rds (MUCH faster I/O):
#   fname.rdata.i <- file.path(dir.box.stroop.data, fnames[atlas.i])
#   if (!file.exists(fname.rdata.i)) next
#   rsarray <- load.rdata(fname.rdata.i)
#   saveRDS(rsarray, gsub("RData", "rds", fname.rdata.i))
#   if (file.exists(fname.rds.i)) file.remove(fname.rdata.i)
#   
#   # if (should.write.csvs) {
#   #   
#   #   ## NB: heavy lifting, depending on num and size of parcels:
#   #   betalist.i <- load.rdata(file.path(dir.box.stroop.data, fnames[atlas.i]))
#   #   
#   #   ## get subject lists from data files:
#   #   betalist.index <- colsplit(names(betalist.i), "_", names = c("subj", "parcel", "hemi"))  ## long-form
#   #   subjs.i <- unique(betalist.index[, "subj"])
#   #   parcels.i <- unique(betalist.index[, "parcel"])
#   #   hemis.i <- unique(betalist.index[, "hemi"])
#   #   
#   #   ## for the RDS file:
#   #   rsarray <- array(
#   #     NA,
#   #     dim = c(length(bias.items), length(bias.items), length(subjs.i), length(parcels.i), length(hemis.i)),
#   #     dimnames = list(
#   #       r    = bias.items,
#   #       c    = bias.items,
#   #       subj = subjs.i,
#   #       parc = parcels.i,
#   #       hemi = hemis.i
#   #     )
#   #   )
#   #   
#   #   
#   #   for (subj.j in seq_along(subjs.i)) {
#   #     # subj.j = 1
#   #     
#   #     ## create relevant dirs (only if don't exist)
#   #     dir.subj.atlas.j <- file.path(dir.box.stroop.data, subjs.i[subj.j], atlases[atlas.i])
#   #     if (!dir.exists(dir.subj.atlas.j)) dir.create(dir.subj.atlas.j, recursive = TRUE)
#   #     
#   #     for (parcel.k in seq_along(parcels.i)) {
#   #       # parcel.k = 1
#   #       for (hemi.l in seq_along(hemis.i)) {
#   #         # hemi.l = 1
#   #         
#   #         ## get element (parcel betas):
#   #         elem.name <- paste0(subjs.i[subj.j], "_", parcels.i[parcel.k], "_", hemis.i[hemi.l])
#   #         elem <- betalist.i[[elem.name]]  ## a single subj's betas for a single parcel & hemisphere
#   #         
#   #         ## make parcel dir:
#   #         dir.parcel.k <- file.path(dir.subj.atlas.j, paste0(parcels.i[parcel.k], "_", hemis.i[hemi.l]))
#   #         if (!dir.exists(dir.parcel.k)) dir.create(dir.parcel.k)
#   #         
#   #         ## write files:
#   #         file.suffix <- paste0(
#   #           paste0(subjs.i[subj.j], "_", atlases[atlas.i], "_", parcels.i[parcel.k], "_", hemis.i[hemi.l]), "_", method, 
#   #           ".csv"
#   #         )
#   #         fwrite(
#   #           as.data.table(elem),  ## to avoid coercion warning
#   #           file = file.path(dir.parcel.k, paste0("betas_", file.suffix))
#   #         )
#   #         rsm <- cor(elem[, bias.items])  ## get correlation matrix for relevant regressors
#   #         fwrite(
#   #           as.data.table(rsm),  ## to avoid coercion warning
#   #           file = file.path(dir.parcel.k, paste0("rsm-pearson_", file.suffix))
#   #         )
#   #         rsarray[, , subj.j, parcel.k, hemi.l] <- rsm  ## add to rsarray
#   #         
#   #       }  ## end hemi loop
#   #     }  ## end parcel loop
#   #     
#   #     print(paste0(atlases[atlas.i], "  ", subj.j))
#   #     
#   #   }  ## end subj loop
#   #   
#   # }
#   
# }  ## end atlas loop
# 
# 
# 
# ## NOT RUN {
# ## 
# # beg.time <- Sys.time()  
# # atlas.i <- 1
# # ## testing efficiency of reading all csvs...
# # rsarray <- array(  ## for representational similarity matrices
# #   NA,
# #   dim = c(length(bias.items), length(bias.items) + 1, length(subjs.i), length(parcels.i), length(hemis.i)),
# #   dimnames = list(
# #     r    = bias.items,
# #     c    = c("rn", bias.items),
# #     subj = subjs.i,
# #     roi  = parcels.i,
# #     hemi = hemis.i
# #   )
# # )
# # for (subj.j in seq_along(subjs.i)) {
# #   # subj.j = 1
# #   for (parcel.k in seq_along(parcels.i)) {
# #     # parcel.k = 1
# #     for (hemi.l in seq_along(hemis.i)) {
# #       # hemi.l = 1
# #       fname.jkl <- file.path(
# #         dir.box.stroop.data, subjs.i[subj.j], atlases[atlas.i],
# #         paste0(parcels.i[parcel.k], "_", hemis.i[hemi.l]),
# #         paste0(
# #           "rsm-pearson_", subjs.i[subj.j], "_", atlases[atlas.i], "_", parcels.i[parcel.k], "-", hemis.i[hemi.l], "_",
# #         method, ".csv"
# #         )
# #       )
# #       rsarray[, , subj.j, parcel.k, hemi.l] <- as.matrix(fread(fname.jkl))
# #     }  ## end parcel loop
# #   }  ## end hemi loop
# # }  ## end subj loop
# # end.time <- Sys.time()  
# # beg.time - end.time
# ## } END NOT RUN RUN