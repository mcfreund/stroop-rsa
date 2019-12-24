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
write.csv(stroop, "data/summary_group201902.csv", row.names = FALSE)

