## about ----
## copies behavioral (RT, ACC) and event info (stim times) from Box/DMCC_Phase2(HCP)/Preprocessed_Data
## NOT meant to be part of the pipleine.
## I.e., this script is for ONE-TIME USE.
## by writing the data from BOX to the project stroop-rsa/data directory, external pointers are avoided.
## THIS IS THE ONLY SCRIPT THAT WRITES TO THE data DIRECTORY.
##
## mike freund, 2019-02-20
## adapted for new project directory 2019-12-24


subj_sum <- read.csv(
  "C:/Users/mcf/Box/DMCC_Phase2(HCP)/Preprocessed_Data/_wrangled/summary/dmcc2_all-runs-summary.csv",
  stringsAsFactors = FALSE
)
stroop <- read.csv(
  "C:/Users/mcf/Box/DMCC_Phase2(HCP)/Preprocessed_Data/_wrangled/dmcc2_behavior-and-events_stroop.csv",
  stringsAsFactors = FALSE
)
subjs <- read.csv(
  "C:/Users/mcf/Box/global/wustl/proj/stroop-rsa/sheets/subj-sample_group-201902.csv",
  stringsAsFactors = FALSE
)$x
subjs.foranalysis <- read.csv(
  "C:/Users/mcf/Box/global/wustl/proj/stroop-rsa/sheets/subj-sample_group-201902_for-analysis.csv",
  stringsAsFactors = FALSE
)$x

subjlist <- data.frame(
  subj = subjs,
  is.analysis.group = subjs %in% subjs.foranalysis,
  strings.as.factors = FALSE
)

stroop <- subset(stroop, subj %in% subjlist$subj)
subj_sum <- subset(subj_sum, subj %in% subjlist$subj)

stroop$is.analysis.group <- stroop$subj %in% subjlist$subj[subjlist$is.analysis.group]
subj_sum$is.analysis.group <- subj_sum$subj %in% subjlist$subj[subjlist$is.analysis.group]

write.csv(stroop, "data/behavior-and-events_group201902.csv", row.names = FALSE)
write.csv(subj_sum, "data/summary_group201902.csv", row.names = FALSE)

