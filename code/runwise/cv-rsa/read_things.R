## read info ----

# subjsum <- fread(here("in", "summary_group201902.csv"))
# subjsum %>% filter(session == "pro", task == "stp", run == 1) %>%
#   filter(!is.na(twin.pair)) %>%
#   filter(duplicated(twin.pair))


behav <- fread(here("out", "behav", "behavior-and-events_group201902_with-subset-cols.csv"))


N.total <- length(unique(behav$subj))
N.valid <- length(unique(behav$subj[!behav$is.analysis.group]))
N.analy <- length(unique(behav$subj[behav$is.analysis.group]))

## behav

micinfo <- fread(here("in", "microphone-info_group201902.csv"))
## warnings due to text encoding (nonissue):
behav.mod.objs <- suppressWarnings(readRDS(here("out", "behav", "mod_objs.RDS")))  ## behavioral model objects
hyp.cors <- fread(here("out", "indiv", "hyp_bivariate.txt"))



## group

table.group.means <- fread(here("out", "group", "superparcels_means.txt"))
table.group.wnregion <- fread(here("out", "group", "superparcels_wnregion.txt"))
table.group.bnregion <- fread(here("out", "group", "superparcels_bnregion.txt"))

table.group.alt <- fread(here("out", "group", "superparcels_alt.txt"))
ceilings.p <- fread(here("out", "group", "ceilings.txt"))
ceilings.m <- fread(here("out", "group", "ceilings_means.txt"))

## indiv

lfp.cors <- fread(here("out", "indiv", "cors_lfp.txt"))
tab.single <- fread(here("out", "indiv", "table_single.csv"))
tab.wnroi <- fread(here("out", "indiv", "tab_wnroi.csv"))
tab.bnroi <- fread(here("out", "indiv", "tab_bnroi.csv"))



dem <- readxl::read_xls("C:/Users/mcf/Box/DMCC_Phase2(HCP)/DataAndDemographics.xls")[, 1:2]
names(dem) <- c("subj", "gender")
dem <- filter(dem, subj %in% behav$subj)
dem$is.analysis.group <- dem$subj %in% unique(behav$subj[behav$is.analysis.group])


## explor

`table.group.mmp.targt` <- fread(here("out", "explor", "table_target.csv"), colClasses = c(md = "character"))


mov0 <- fread(here("out", "explor", "movregs_rsa_vs0.csv"))
movd <- fread(here("out", "explor", "movregs_rsa_delta.csv"))

mov06_p <- min(mov0[nmovregs == "m06", p])
mov06_b <- max(mov0[nmovregs == "m06", b])
mov12_p <- min(mov0[nmovregs == "m12", p])
mov12_b <- max(mov0[nmovregs == "m12", b])
movd_p <- min(movd[, p])
movd_b <- abs(max(movd[, b]))