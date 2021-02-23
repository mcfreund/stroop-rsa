source(here::here("code", "packages.R"))
source(here("code", "strings.R"))

subjs <- c(subjs.analysis, subjs.validation)
conds <- dimnames(read_xmat(here("glms", subjs[1], "results", "pro_bias_acc-only", "X.xmat.1D")))$regressor
conds <- gsub("\\[|\\]|#0|#", "", conds)

n.tr <- 1080
n.cond <- length(conds)
n.subj <- length(subjs)


X_glm <- array(
  NA,
  dim = c(n.tr, n.cond, n.subj),
  dimnames = list(tr = NULL, condition = conds, subj = subjs)
)


for (subj.i in seq_along(subjs)) {
  # subj.i = 1
  
  X_glm.i <- read_xmat(here("glms", subjs[subj.i], "results", "pro_bias_acc-only", "X.xmat.1D"))
  colnames(X_glm.i) <- gsub("\\[|\\]|#0|#", "", colnames(X_glm.i))
  X_glm[, , subj.i] <- X_glm.i[, conds]
  
}

saveRDS(X_glm, here("glms", "xmats_pro_bias_acc-only.rds"))
