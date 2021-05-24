source(here::here("code", "packages.R"))
source(here("code", "strings.R"))

subjs <- c(subjs.analysis, subjs.validation)
subjs <- setdiff(subjs, "562345")
n.subj <- length(subjs)


X_glm <- setNames(vector("list", length(subjs)), subjs)

for (subj.i in seq_along(subjs)) {
  # subj.i = 1
  
  X1 <- read_xmat(here("glms", subjs[subj.i], "results", "pro_bias_acc-only_fmriprep_im_run1", "X_1.xmat.1D"))
  X2 <- read_xmat(here("glms", subjs[subj.i], "results", "pro_bias_acc-only_fmriprep_im_run2", "X_2.xmat.1D"))
  colnames(X1) <- gsub("\\[|\\]|#0|#", "", colnames(X1))
  colnames(X2) <- gsub("\\[|\\]|#0|#", "", colnames(X2))
  
  X_glm[[subj.i]] <- list(run1 = X1, run2 = X2)
  
}

saveRDS(X_glm, here("glms", "xmats_pro_bias_acc-only_fmriprep_im.rds"))
