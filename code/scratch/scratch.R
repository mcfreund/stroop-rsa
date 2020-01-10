# ## subset by sample
# rsarray.analysis <- rsarray[, , sample.analysis, , ]
# rsarray.held.out <- rsarray[, , sample.held.out, , ]


## list of subjects for analysis
# 
# stroop <- fread(here("data", "behavior-and-events_group201902.csv"))
# sample.analysis <- unique(filter(stroop, is.analysis.group)$subj)
# sample.held.out <- unique(filter(stroop, !is.analysis.group)$subj)
# 
# model.names <- c("target", "distractor", "congruency")
# # model.names <- c("target", "distractor", "congruency", "incongruency")
# model.method <- "pcor"
# # stats <- "pearson"
# # group <- "group-201902"
# # method <- "pro_bias_acc-only"


# ## summarize:
# results <- coefs %>%
#   group_by(parcel, hemi, variable) %>%
#   summarize(
#     mean.rho     = mean.rzr(rho),  ## tanh(mean(atanh(x)))
#     mean.r       = mean.rzr(r),
#     rho.srtest.v = wilcox.test(rho, alternative = "greater")$statistic,
#     r.srtest.v   = wilcox.test(r, alternative = "greater")$statistic,
#     rho.srtest.p = wilcox.test(rho, alternative = "greater")$p.value,
#     r.srtest.p   = wilcox.test(r, alternative = "greater")$p.value
#   ) %>%
#   rename(rho = mean.rho, r = mean.r)  ## for consistency
