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


pca.incon.coding <- prcomp(scale(incon[, rois.incon]))
pca.incon.coding$sdev
plot(pca.incon.coding)
plot.loadings(pca.incon.coding$rotation[, c(1, 2)])
plot.loadings(pca.incon.coding$rotation[, c(2, 3)])
plot.loadings(pca.incon.coding$rotation[, c(1, 3)])


## cross-validaiton to select number of components
# https://stats.stackexchange.com/questions/93845/how-to-perform-cross-validation-for-pca-to-determine-the-number-of-principal-com


loopca <- function(d, K = "max") {
  ## see here:
  ## https://stats.stackexchange.com/questions/93845/how-to-perform-cross-validation-for-pca-to-determine-the-number-of-principal-com
  
  if (K == "max") K <- min(nrow(d) - 1, ncol(d) - 2)
  if (K < 2) stop("K > 2!")
  
  mse <- numeric(K - 1)
  
  for (k in seq_len(K)) {  ## for a particular number of PCs...
    
    features <- seq_len(k)
    
    if (k == 1) next
    
    mse.k <- numeric(nrow(d))  ## holds mse per observation ii
    
    for (ii in seq_len(nrow(d))) {  ## drop obs ii
  
      train <- d[-ii, ]
      test  <- d[ii, ]
      
      U <- prcomp(train)$rotation[, features]
      
      testhat <- numeric(k)
      for (j in features) {  ## drop features
        testhat[j] <- (U %*% MASS::ginv(U[-j, ]) %*% test[-j])[j]  ## hmmm...
      }
      
      mse.k[ii] <- mean((test[features] - testhat)^2)  ## mse across features for obs ii

    }
  
    mse[k - 1] <- mean(mse.k)  ## average MSE across folds (? is this kosher ?)
    
  }
  
  mse
  
}
thresh <- function(x, tile) names(x)[x > quantile(x, tile)]

## target

mse.pca.targt.rois <- loopca(targt[, rois.targt])
plot(mse.pca.targt.rois)
which.min(mse.pca.targt.rois) + 1

mse.pca.targt.rois <- loopca(t(scale(t(scale(targt[, rois.targt])))))
plot(mse.pca.targt.rois)
which.min(mse.pca.targt.rois)
which.min(mse.pca.targt.rois) + 1

mse.pca.targt.rois <- loopca(t(scale(t(targt[, rois.targt]))))
plot(mse.pca.targt.rois)
which.min(mse.pca.targt.rois)
which.min(mse.pca.targt.rois) + 1


pca.targt.coding <- prcomp(t(scale(t(scale(targt[, rois.targt])))))$rotation[, 1:4]

targt.top <- c(
  thresh(abs(pca.targt.coding[, 1]), 0.99),
  thresh(abs(pca.targt.coding[, 2]), 0.99),
  thresh(abs(pca.targt.coding[, 3]), 0.99),
  thresh(abs(pca.targt.coding[, 4]), 0.99)
)

atlas.key$mmp %>% filter(roi %in% targt.top)

## incongruency

mse.pca.incon.rois <- loopca(incon[, rois.incon])
plot(2:(length(mse.pca.incon.rois) + 1), mse.pca.incon.rois)
which.min(mse.pca.targt.rois) + 1

mse.pca.incon.rois <- loopca(scale(incon[, rois.incon]))
plot(2:(length(mse.pca.incon.rois) + 1), mse.pca.incon.rois)
which.min(mse.pca.targt.rois) + 1

mse.pca.incon.rois <- loopca(t(scale(t(scale(incon[, rois.incon])))))
plot(2:(length(mse.pca.incon.rois) + 1), mse.pca.incon.rois)
which.min(mse.pca.targt.rois) + 1

mse.pca.incon.rois <- loopca(t(scale(t(incon[, rois.incon]))))
plot(2:(length(mse.pca.incon.rois) + 1), mse.pca.incon.rois)
which.min(mse.pca.targt.rois) + 1

pca.incon.coding <- prcomp(scale(incon[, rois.targt]))$rotation[, 1:4]
apply(pca.incon.coding, 2, range)

incon.top <- c(
  thresh(abs(pca.incon.coding[, 1]), 0.99),
  thresh(abs(pca.incon.coding[, 2]), 0.99),
  thresh(abs(pca.incon.coding[, 3]), 0.99),
  thresh(abs(pca.incon.coding[, 4]), 0.99)
)

atlas.key$mmp %>% filter(roi %in% incon.top)

