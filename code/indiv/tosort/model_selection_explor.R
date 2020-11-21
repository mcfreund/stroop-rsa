#+ model-selection-whole-cortex_fit

rm.these <- grep("stroop|^congr$", colnames(m.mmp))
X.mmp <- m.mmp[ids$is.analysis.group, -rm.these]  ## analysis set

lambdas.mmp <- tune_lambda(X.mmp, strooprt, alpha = 0.5, selection_crit = function(fit) fit$lambda.min)
hist(lambdas.mmp, breaks = 50, col = "grey50")
lambda.best.mmp <- min(Mode(lambdas.mmp))

net.mmp <- glmnet(x = X.mmp, y = strooprt, lambda = lambda.best.mmp, alpha = 0.5)
coef(net.mmp)
#+
#' None selected, even with liberal minimum lambda criterion.


#' ## elastic net model: MD regions (MMP)

#+ model-selection-md_fit

keep.these <- grep(
  paste0(atlas.key$mmp$roi[atlas.key$mmp$md != "none"], collapse = "|"),
  colnames(m.mmp)
)
X.md <- m.mmp[ids$is.analysis.group, keep.these]  ## analysis set

lambdas.md <- tune_lambda(X.md, strooprt, alpha = 0.5, selection_crit = function(fit) fit$lambda.min)
hist(lambdas.md, breaks = 50, col = "grey50")
lambda.best.md <- min(Mode(lambdas.md))
net.md <- glmnet(x = X.md, y = strooprt, lambda = lambda.best.md, alpha = 0.5)
coef(net.md)


## estimate test error ----

X.md_vset <- m.mmp[!ids$is.analysis.group, keep.these]  ## 'held out' / validation set

## observed error

ys.md <- cbind(
  y = c(strooprt_vset),
  yhat = c(predict(net.md, newx = X.md_vset))
)

(cor.obs.md <- cor(ys.md)["y", "yhat"])

ys.md %>%
  as.data.frame %>%
  ggplot(aes(yhat, y)) +
  stat_boot_ci(n = 1e4, alpha = 0.3, color = "grey50") +
  geom_point()

## permuted validation error

n_cores <- detectCores() - 1
nresamps <- 1E4
cl <- makeCluster(n_cores - 1)
registerDoParallel(cl)
cor.perm.md <- foreach(ii = seq_len(nresamps), .inorder = FALSE, .combine = "c", .packages = "glmnet") %dorng% {
  
  yperm <- strooprt[sample.int(length(strooprt))]
  
  fit.ii <- glmnet(x = X.md, y = yperm, lambda = lambda.best.md, alpha = 0.5)
  
  cor(strooprt_vset, predict(fit.ii, newx = X.md_vset))
  
}
stopCluster(cl)

(p.md <- sum(cor.perm.md > cor.obs.md, na.rm = TRUE) / (nresamps - sum(is.na(cor.perm.md))))

plot(density(cor.perm.md, na.rm = TRUE))
abline(v = cor.obs.md, col = "firebrick", lwd = 3)

sum(is.na(cor.perm.md))  ## num models with no predictors
#+


#+ model-selection-explor_save
## save ----

saveRDS(net.md, here("out", "indiv", "selected_model_md.RDS"))
saveRDS(
  list(ys.md = ys.md, p.md = p.md), 
  here("out", "indiv", "selected_model_md_validation.RDS")
)

#+