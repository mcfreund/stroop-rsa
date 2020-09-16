#+ model-selection-alt_fit
## estimate model ----

rm.these.alt <- grep("stroop|^congr$|V1|dlpfc_|dmfc_", colnames(m.super))
X.super.alt <- m.super[ids$is.analysis.group, -rm.these.alt]  ## analysis set

lambdas.super.alt <- tune_lambda(X.super.alt, strooprt, alpha = 0.5, selection_crit = function(fit) fit$lambda.min)
hist(lambdas.super.alt, breaks = 50, col = "grey50")
lambda.best.super.alt <- min(Mode(lambdas.super.alt))

net.super.alt <- glmnet(x = X.super.alt, y = strooprt, lambda = lambda.best.super.alt, alpha = 0.5)
coef(net.super.alt)


## estimate test error ----

X.super.alt_vset <- m.super[!ids$is.analysis.group, -rm.these.alt]  ## 'held out' / validation set

## observed error

ys.super.alt <- cbind(
  y = c(strooprt_vset),
  yhat = c(predict(net.super.alt, newx = X.super.alt_vset))
)

(cor.obs.super.alt <- cor(ys.super.alt)["y", "yhat"])

ys.super.alt %>%
  as.data.frame %>%
  ggplot(aes(yhat, y)) +
  stat_boot_ci(n = 1e4, alpha = 0.3, color = "grey50") +
  geom_point()

## permuted validation error

n_cores <- detectCores() - 1
nresamps <- 1E4
cl <- makeCluster(n_cores - 1)
registerDoParallel(cl)
cor.perm.super.alt <- foreach(ii = seq_len(nresamps), .inorder = FALSE, .combine = "c", .packages = "glmnet") %dorng% {
  
  yperm <- strooprt[sample.int(length(strooprt))]
  
  fit.ii <- glmnet(x = X.super.alt, y = yperm, lambda = lambda.best.super.alt, alpha = 0.5)
  
  cor(strooprt_vset, predict(fit.ii, newx = X.super.alt_vset))
  
}
stopCluster(cl)

(p.super.alt <- sum(cor.perm.super.alt > cor.obs.super.alt, na.rm = TRUE) / (nresamps - sum(is.na(cor.perm.super.alt))))

plot(density(cor.perm.super.alt, na.rm = TRUE))
abline(v = cor.obs.super.alt, col = "firebrick", lwd = 3)

sum(is.na(cor.perm.super.alt))  ## num models with no predictors

## save ----

saveRDS(net.super, here("out", "indiv", "selected_model_alt.RDS"))
saveRDS(list(ys.super = ys.super, p.super = p.super), here("out", "indiv", "selected_model_validation_alt.RDS"))

#+


