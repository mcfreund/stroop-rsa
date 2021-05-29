#+ model-selection_fit
## estimate model ----

rm.these <- grep("stroop|^congr$|V1|\\.alt", colnames(m.super))
X.super <- m.super[ids$is.analysis.group, -rm.these]  ## analysis set

lambdas.super <- tune_lambda(X.super, strooprt, alpha = 0.5, selection_crit = function(fit) fit$lambda.1se)
hist(lambdas.super, breaks = 50, col = "grey50")
lambda.best.super <- min(Mode(lambdas.super))

## seq lambdas b/c warm start:
net.super <- glmnet(x = X.super, y = strooprt, lambda = seq(0, lambda.best.super, 1), alpha = 0.5)
coef(net.super, s = lambda.best.super)





## estimate test error ----

X.super_vset <- m.super[!ids$is.analysis.group, -rm.these]  ## 'held out' / validation set

## observed error

ys.super <- cbind(
  y = c(strooprt_vset),
  yhat = c(predict(net.super, newx = X.super_vset, s = lambda.best.super))
)

(cor.obs.super <- cor(ys.super)["y", "yhat"])




ys.super %>%
  as.data.frame %>%
  ggplot(aes(yhat, y)) +
  stat_boot_ci(n = 1e4, alpha = 0.3, color = "grey50") +
  geom_point()

## permuted validation error

n_cores <- detectCores() - 1
nresamps <- 1E4
cl <- makeCluster(n_cores - 1)
registerDoParallel(cl)
cor.perm.super <- foreach(ii = seq_len(nresamps), .inorder = FALSE, .combine = "c", .packages = "glmnet") %dorng% {
  
  yperm <- strooprt[sample.int(length(strooprt))]
  
  fit.ii <- glmnet(x = X.super, y = yperm, lambda = seq(0, lambda.best.super, 1), alpha = 0.5)
  
  cor(strooprt_vset, predict(fit.ii, newx = X.super_vset, s = lambda.best.super))
  
}
stopCluster(cl)

(p.super <- sum(cor.perm.super > cor.obs.super, na.rm = TRUE) / (nresamps - sum(is.na(cor.perm.super))))

plot(density(cor.perm.super, na.rm = TRUE))
abline(v = cor.obs.super, col = "firebrick", lwd = 3)

sum(is.na(cor.perm.super))  ## num models with no predictors
#+


#+ model-selection_save
## save ----

saveRDS(net.super, here("out", "indiv", "selected_model.RDS"))
saveRDS(
  list(ys.super = ys.super, p.super = p.super, cor.perm.super = cor.perm.super), 
  here("out", "indiv", "selected_model_validation.RDS")
  )

#+