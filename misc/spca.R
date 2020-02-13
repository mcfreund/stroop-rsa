library(here)
library(knitr)
library(magrittr)
library(dplyr)
library(data.table)
library(ggplot2)
library(mikeutils)
library(doParallel)
library(foreach)
library(gifti)
library(fANCOVA)
library(viridis)
library(colorspace)
library(elasticnet)
source(here("code", "strings.R"))
source(here("code", "read_atlases.R"))


stats.subjs.tdic <- fread(
  here("out", "rsa", "stats", paste0("subjs_pro_bias_acc-only_mmp_pearson_residual_glm-tdic.csv"))
)

stats.subjs.tdic <- stats.subjs.tdic[is.analysis.group == TRUE & y == "rank", ]  ## EXCLUDE HELD OUT SUBJECTS!
stats.subjs.tdic <- stats.subjs.tdic[, "coef" := NULL]

stats.subjs.tdic %<>% full_join(atlas.key$mmp, by = "roi")

params.interest <- c("target", "distractor", "incongruency")

## get stats ----

## coding

stats.group.tdic <- stats.subjs.tdic %>%
  group_by(num.roi, model, param) %>%
  summarize(
    v    = wilcox.test(beta, alternative = "greater")$statistic,
    p    = wilcox.test(beta, alternative = "greater")$p.value,
    beta = tanh(mean(atanh(beta))),  ## must be last in this summarize()
  ) %>%
  mutate(p.fdr = p.adjust(p, method = "fdr"), p.holm = p.adjust(p, method = "holm")) %>%
  ungroup
stats.group.tdic %<>% full_join(atlas.key$mmp, by = "num.roi") %>% as.data.table


rois.targt <- stats.group.tdic[param == "target" & p.fdr < 0.05]$roi
rois.distr <- stats.group.tdic[param == "distractor" & p.fdr < 0.05]$roi
rois.incon <- stats.group.tdic[param == "incongruency" & p.fdr < 0.05]$roi

## preference
stats.pairs.tdic <- stats.subjs.tdic %>%
  filter(param %in% params.interest) %>%
  group_by(roi) %>%
  summarize(
    out = list(
      pairwise.wilcox.test(beta, param, paired = TRUE, alternative = "two.sided", p.adjust.method = "fdr")
    )
  )
pvals <- vapply(stats.pairs.tdic$out, function(.) .$p.value[lower.tri(diag(2), diag = TRUE)], numeric(3))
pvals <- t(pvals)
colnames(pvals) <- c("p.incon.distr", "p.targt.distr", "p.targt.incon")
stats.pairs.tdic <- data.frame(roi = stats.pairs.tdic$roi, pvals, stringsAsFactors = FALSE)

rois.pref.targt.vs.distr <- stats.pairs.tdic[roi %in% rois.targt & p.targt.distr < 0.05 & b.targt.distr > 0, roi]
rois.pref.targt.vs.incon <- stats.pairs.tdic[roi %in% rois.targt & p.targt.incon < 0.05 & b.targt.incon > 0, roi]

rois.pref.distr.vs.targt <- stats.pairs.tdic[roi %in% rois.distr & p.targt.distr < 0.05 & b.targt.distr < 0, roi]
rois.pref.distr.vs.incon <- stats.pairs.tdic[roi %in% rois.distr & b.incon.distr < 0.05 & b.incon.distr < 0, roi]

rois.pref.incon.vs.distr <- stats.pairs.tdic[roi %in% rois.incon & b.incon.distr < 0.05 & b.incon.distr > 0, roi]
rois.pref.incon.vs.targt <- stats.pairs.tdic[roi %in% rois.incon & p.targt.incon < 0.05 & b.targt.incon < 0, roi]

rois.pref.targt <- intersect(rois.pref.targt.vs.incon, rois.pref.targt.vs.distr)
rois.pref.distr <- intersect(rois.pref.distr.vs.targt, rois.pref.distr.vs.incon)
rois.pref.incon <- intersect(rois.pref.incon.vs.targt, rois.pref.incon.vs.distr)

stats.pairs.tdic %<>% full_join(
  stats.subjs.tdic %>%
    filter(param %in% params.interest) %>%
    group_by(roi) %>%
    summarize(
      b.incon.distr = tanh(mean(atanh(beta[param == "incongruency"]) - atanh(beta[param == "distractor"]))),
      b.targt.distr = tanh(mean(atanh(beta[param == "target"]) - atanh(beta[param == "distractor"]))),
      b.targt.incon = tanh(mean(atanh(beta[param == "target"]) - atanh(beta[param == "incongruency"]))),
    ),
  by = "roi"
)

stats.pairs.tdic %<>% full_join(atlas.key$mmp, by = "roi") %>% as.data.table


## wrangle to matrix ---

# qcor(cor(targt))

stats.w <- stats.subjs.tdic[, c("subj", "beta", "param", "roi")] %>% tidyr::spread(param, beta)

targt <- stats.w %>% select(subj, roi, target) %>% tidyr::spread(roi, target)
rownames(targt) <- targt$subj
targt$subj <- NULL
targt <- as.matrix(targt)

distr <- stats.w %>% select(subj, roi, distractor) %>% tidyr::spread(roi, distractor)
rownames(distr) <- distr$subj
distr$subj <- NULL
distr <- as.matrix(distr)

incon <- stats.w %>% select(subj, roi, incongruency) %>% tidyr::spread(roi, incongruency)
rownames(incon) <- incon$subj
incon$subj <- NULL
incon <- as.matrix(incon)


## pca ----

plot.loadings <- function(d) {
  plot(d)
  text(d[, 1], d[, 2], labels = row.names(d), pos = 1)
}


pca.targt.coding <- prcomp(targt[, rois.targt])
plot(pca.targt.coding)
plot.loadings(pca.targt.coding$rotation[, c(1, 2)])
plot.loadings(pca.targt.coding$rotation[, c(2, 3)])
plot.loadings(pca.targt.coding$rotation[, c(1, 3)])


r <- t(pca.targt.coding$rotation[, 1:2] %*% t(pca.targt.coding$rotation[, 1:2])) %*% as.matrix(targt[1, rois.targt])
cor(r, targt[1, rois.targt])

pca.targt.pref.distr <- prcomp(targt[, rois.pref.targt.vs.distr])
plot(pca.targt.pref.distr)
plot.loadings(pca.targt.pref.distr$rotation[, c(1, 2)])
plot.loadings(pca.targt.pref.distr$rotation[, c(2, 3)])
plot.loadings(pca.targt.pref.distr$rotation[, c(1, 3)])


pca.incon.coding <- prcomp(scale(incon[, rois.incon]))
pca.incon.coding$sdev
plot(pca.incon.coding)
plot.loadings(pca.incon.coding$rotation[, c(1, 2)])
plot.loadings(pca.incon.coding$rotation[, c(2, 3)])
plot.loadings(pca.incon.coding$rotation[, c(1, 3)])


## cross-validaiton to select number of components
# https://stats.stackexchange.com/questions/93845/how-to-perform-cross-validation-for-pca-to-determine-the-number-of-principal-com
# loopca <- function(d, K = "max") {
#   
#   if (K == "max") K <- min(nrow(d) - 1, ncol(d))
# 
#   mse <- numeric(K)
#   
#   for (k in seq_len(K)) {
#     
#     mse.k <- numeric(nrow(d))
#     
#     for (ii in seq_len(nrow(d))) {
#       
#       train <- d[-ii, ]
#       test  <- d[ii, ]
#       rot <- prcomp(train)$rotation[, seq_len(k)]
#       # D <- rot %*% t(rot)
#       # test.hat <- D %*% test  ## project test onto train components
#       test.hat <- fitted(lm(test ~ rot))
#       mse.k[ii] <- mean((test - test.hat)^2)
#       
#     }
#     
#     mse[k] <- mean(mse.k)  ## average MSE across folds (? is this kosher ?)
#     
#   }
#   
#   mse
#   
# }
# NB: THIS IS WRONG!
# mse.pca.targt.rois <- loopca(targt[, rois.targt])

plot(mse.pca.targt.rois)


pca.targt.coding <- prcomp(targt[, rois.targt])
plot(pca.targt.coding)
plot.loadings(pca.targt.coding$rotation[, c(1, 2)])
plot.loadings(pca.targt.coding$rotation[, c(2, 3)])
plot.loadings(pca.targt.coding$rotation[, c(1, 3)])


r <- t(pca.targt.coding$rotation[, 1:2] %*% t(pca.targt.coding$rotation[, 1:2])) %*% as.matrix(targt[1, rois.targt])
cor(r, targt[1, rois.targt])




## sparse pca ----



n.feats <- 5
spca.targt <- spca(
  targt[, rois.targt],
  K = n.feats, type = "predictor", sparse = "penalty", para = rev(c(10, 0.5, 0.1, 0.01, 0.001)),
  use.corr = TRUE,
  trace = TRUE
)

spca.targt$var.all


spca.targt$loadings[abs(spca.targt$loadings[, 1]) > 0, 1] %>% sort
spca.targt$loadings[abs(spca.targt$loadings[, 2]) > 0, 2] %>% sort
spca.targt$loadings[abs(spca.targt$loadings[, 3]) > 0, 3] %>% sort
spca.targt$loadings[abs(spca.targt$loadings[, 4]) > 0, 4] %>% sort
spca.targt$loadings[abs(spca.targt$loadings[, 5]) > 0, 5] %>% sort
