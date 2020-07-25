

## fit "null" model
# 
# lmer0 <- lmer(
#   rt.s ~ trial.type + (trial.type | subj) + (trial.type | color),
#   d.dissoc.hlm,
#   control = lmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 1E9))
#   )
# 
# mods.lmer <- list(
#   dlpfc_L_targ = update(lmer0, . ~ . + trial.type * dlpfc_L_target),
#   dlpfc_R_targ = update(lmer0, . ~ . + trial.type * dlpfc_R_target),
#   lppc_L_targ  = update(lmer0, . ~ . + trial.type * lppc_L_target),
#   lppc_R_targ  = update(lmer0, . ~ . + trial.type * lppc_R_target),
#   mfc_L_incon  = update(lmer0, . ~ . + trial.type * mfc_L_incongruency),
#   mfc_R_incon  = update(lmer0, . ~ . + trial.type * mfc_R_incongruency)
# )
# sums.lmer <- lapply(mods.lmer, summary)
# pvals.lmer <- lapply(sums.lmer, function(.) coef(.)[4, "Pr(>|t|)"])
# pvals.lmer <- reshape2::melt(bind_rows(pvals.lmer), value.name = "p", variable = "id")
# pvals.lmer$group  <- rep(1:3, each = 2)
# pvals.lmer %<>% group_by(group) %>% mutate(p.fdr = p.adjust(p, method = "fdr"))

# set.seed(0)
# cors <- w.super.s %>%
#   
#   filter(is.analysis.group) %>%
#   select(subj, matches("dlpfc|lppc|mfc")) %>%
#   
#   mutate(subj = as.factor(subj)) %>%
#   
#   reshape2::melt(value.name = "beta", variable.name = "id") %>%
#   
#   left_join(w.super.s[c("stroop", "subj")], by = "subj") %>%
# 
#   group_by(id) %>%
#   summarize(
#     r.line = cor(beta, stroop),
#     r.rank = cor(beta, stroop, method = "spearman"),
#     ci.line = cor_ci(cbind(beta, stroop), R = 1E4, conf = 0.95) %>% list,
#     ci.rank = cor_ci(cbind(rank(beta), rank(stroop)), R = 1E4, conf = 0.95) %>% list
#   ) %>%
#   
#   unnest(c(ci.line, ci.rank), names_sep = ".")
# 
# cors





W <- rbind(
  lfp_target_vs_incongruency = c(0, 0, 0, 0, 1, -1, 1, -1, 0, 0),
  mfc_target_vs_incongruency = c(0, 0, 0, 0, 0, 0, 0, 0, 1, -1),
  lfpmfc_target_incongruency = c(0, 0, 0, 0, 0, 0, -1, 1, -1, 1)
)

(contrast.lmer.dissoc.lfp.mfc <- summary(glht(lmer.dissoc.lfp.mfc, W), test = adjusted("none")))


## separate models

lmer.dissoc.mfc <- update(
  lmer0,
  . ~ . +
    trial.type * mfc_L_target +
    trial.type * mfc_L_incongruency
)

summary.lmer.dissoc.mfc <- summary(lmer.dissoc.mfc)

lme.dissoc.mfc <- update(
  fit1.het.trim,
  . ~ . +
    trial.type * mfc_L_target +
    trial.type * mfc_L_incongruency,
  data = d.dissoc.hlm
)
summary(lme.dissoc.mfc)

(contrast.lme.dissoc.mfc <- summary(glht(lme.dissoc.mfc, rbind(c(0, 0, 0, 0, 1, -1))), test = adjusted("none")))








<!-- ## bootstrap model comparison -->
  
  <!-- ```{r} -->
  
  <!-- ## validation set: ridge ---- -->
  
  <!-- w.super <- cbind(w.super, lfc_R_target = (w.super[, "dlpfc_R_target"] + w.super[, "lppc_R_target"]) / 2) %>% scale -->
  <!-- w.super <- cbind(w.super, lfc_R_incongruency = (w.super[, "dlpfc_R_incongruency"] + w.super[, "lppc_R_incongruency"]) / 2) %>% scale -->
  
  <!-- # X.super      <- w.super[ids$is.analysis.group, cols] -->
  <!-- # y.super      <- w.super[ids$is.analysis.group, "stroop"] -->
  <!-- # Xprime.super <- w.super[!ids$is.analysis.group, cols] -->
  <!-- # yprime.super <- w.super[!ids$is.analysis.group, "stroop"] -->
  
  <!-- ### bootstrap comparison -->
  
  <!-- d <- w.super[ids$is.analysis.group, c("mfc_L_incongruency", "mfc_L_target", "lfc_R_target", "stroop")] -->
  <!-- valid <- w.super[!ids$is.analysis.group, c("mfc_L_incongruency", "mfc_L_target", "lfc_R_target", "stroop")] -->
  
  <!-- Xp1 <- cbind(1, valid[, c("lfc_R_target", "mfc_L_target")]) -->
  <!-- Xp0 <- Xp1[, -3] -->
  <!-- yp = valid[, "stroop"] -->
  <!-- heldout = list(Xp1 = Xp1, Xp0 = Xp0, yp = yp) -->
  <!-- modcomp <- function(d, heldout = heldout, ii) { -->
      
      <!--   dstar <- d[ii, ] -->
        
        <!--   X1 <- cbind(1, dstar[ii, c("lfc_R_target", "mfc_L_target")]) -->
          <!--   X0 <- X1[, -3] -->
            <!--   y = dstar[ii, "stroop"] -->
              
              <!--   fit0 <- .lm.fit(x = X0, y = dstar[, "stroop"]) -->
                <!--   fit1 <- .lm.fit(x = X1, y = dstar[, "stroop"]) -->
                  
                  <!--   yhat0 <- c(heldout$Xp0 %*% as.matrix(fit0$coefficients)) -->
                    <!--   yhat1 <- c(heldout$Xp1 %*% as.matrix(fit1$coefficients)) -->
                      
                      <!--   mse0 <- mean((heldout$yp - yhat0)^2) -->
                        <!--   mse1 <- mean((heldout$yp - yhat1)^2) -->
                          
                          <!--   mse0 - mse1  ## positive if full model is better fit -->
                        
                        <!-- } -->
  
  <!-- results <- boot(d, modcomp, R = 5E4, heldout = heldout) -->
  <!-- boot.ci(results, type = "bca")[["bca"]]  ## negative if mfc_target model better -->

<!-- sum(0 < results$t) / length(results$t)  ## p-value -->



<!-- cols <- c("mfc_L_incongruency", "dlpfc_R_target", "lppc_R_target") -->
  <!-- # cols <- c("mfc_L_incongruency", "lppc_R_target") -->
  <!-- # cols <- c("mfc_L_target", "lfc_R_target") -->
  <!-- # cols <- c("mfc_L_incongruency", "lfc_R_target") -->
  
  <!-- X.super      <- w.super[ids$is.analysis.group, cols] -->
  <!-- y.super      <- w.super[ids$is.analysis.group, "stroop"] -->
  <!-- Xprime.super <- w.super[!ids$is.analysis.group, cols] -->
  <!-- yprime.super <- w.super[!ids$is.analysis.group, "stroop"] -->
  
  
  <!-- ``` -->
  
  <!-- ### 1.b is there a better model? testing all hemisphere\*\scheme\*superparcel rois -->
  
  
  <!-- ```{r superparcels-hyp-modcomp-misc} -->
  
  <!-- lmer0 <- lmer( -->
                        <!--   rt.s ~  -->
                        <!--     trial.type + -->
                        <!--     (trial.type | subj) + (trial.type | color), -->
                        <!--   d.dissoc.hlm, -->
                        <!--   control = lmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 1E9)) -->
                        <!--   ) -->
  <!-- lmer.mfc_L <- lmer( -->
                             <!--   rt.s ~  -->
                             <!--     trial.type * mfc_L_target +  trial.type * mfc_L_incongruency + -->
                             <!--     (trial.type | subj) + (trial.type | color), -->
                             <!--   d.dissoc.hlm, -->
                             <!--   control = lmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 1E9)) -->
                             <!--   ) -->
  <!-- anova(lmer0, lmer.mfc_L) -->
  
  
  <!-- d.super.group <- d.super %>% -->
  
  <!--   filter(is.analysis.group) %>% -->
  <!--   group_by(param, roi) %>% -->
  
  <!--   summarize(p = wilcox.test(beta, alternative = "greater")$p.value) %>% -->
  <!--   ungroup %>% -->
  <!--   mutate(p.fdr = p.adjust(p, method = "fdr")) -->
  
  <!-- d.super.group %>% filter(p < 0.05) %>% View -->
  
  <!-- ``` -->
  
  
  
  
  


## held-out error? ----

w.super.s$lfp_R_target <- (w.super.s$dlpfc_R_target + w.super.s$lppc_R_target) / 2
w.super.s$lfp_R_incongruency <- (w.super.s$dlpfc_R_incongruency + w.super.s$lppc_R_incongruency) / 2
cols <- c("dlpfc_R_target", "lppc_R_target", "mfc_L_target")

X.super      <- w.super.s[ids$is.analysis.group, cols] %>% as.matrix
y.super      <- w.super.s[ids$is.analysis.group, "stroop"] %>% as.matrix
Xprime.super <- w.super.s[!ids$is.analysis.group, cols] %>% as.matrix
yprime.super <- w.super.s[!ids$is.analysis.group, "stroop"] %>% as.matrix

lambdas <- tune_lambda(X.super, y.super, alpha = 0, nreps = 5E3)
hist(lambdas)
lambda.best <- Mode(lambdas)
if (length(lambda.best) > 1) lambda.best <- max(lambda.best)
fit.super <- glmnet(x = X.super, y = y.super, lambda = lambda.best, alpha = 0)
coef(fit.super)

ys.super <- cbind(
  y = c(yprime.super),
  yhat = c(predict(fit.super, newx = Xprime.super))
)

## observed error
lm(y.super ~ X.super)
lm(yprime.super ~ Xprime.super)

# fit.super <- .lm.fit(x = cbind(1, X.super), y = y.super)
# ys.super <- cbind(
#   y = c(yprime.super),
#   yhat = c(cbind(1, Xprime.super) %*% as.matrix(fit.super$coefficients))
#   )

(cor.obs.super <- cor(ys.super)[1, 2])
(mse.obs.super <- mean((ys.super[, "y"] - ys.super[, "yhat"])^2))

ys.super %>%
  as.data.frame %>%
  ggplot(aes(y, yhat)) +
  stat_boot_ci(n = 1e4, alpha = 0.5, color = "grey50") +
  geom_point()


## ols model selection ----

ols.hyp.modcomp <- lm.allcombs(
  as.data.frame(w.super.s[ids$is.analysis.group, c(hyps, "stroop")]),
  "stroop"
)
ols.hyp.modcomp.selected <- ols.hyp.modcomp[[purrr::map_dbl(ols.hyp.modcomp, BIC) %>% which.min]]  ## selected model
summary(ols.hyp.modcomp.selected)


## hlm ----

d.dissoc.hlm <- full_join(
  d.dissoc.hlm,
  cbind(subj = ids$subj[ids$is.analysis.group], as.data.frame(w.super.s[ids$is.analysis.group, hyps])),
  by = "subj"
)

lmer.hyp.selected <- lmer(
  rt.s ~
    trial.type * mfc_L_target +  trial.type * lppc_R_target + trial.type * dlpfc_R_target +
    (trial.type | subj) + (trial.type | color),
  d.dissoc.hlm,
  control = lmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 1E9))
)
summary(lmer.hyp.selected)

lme.hyp.selected <- update(
  fit1.het.trim,
  rt.s ~ . +
    trial.type * mfc_L_target +  trial.type * lppc_R_target + trial.type * dlpfc_R_target,
  data = d.dissoc.hlm
)
summary(lme.hyp.selected)




## estimate model ----

X.super.hyp  <- w.super[ids$is.analysis.group, hyps]
y  <- w.super[ids$is.analysis.group, "stroop"]

## find lambda and fit

set.seed(0)

cl <- makeCluster(n_cores - 1)
registerDoParallel(cl)
lambdas.super.hyp <- foreach(ii = seq_len(nreps), .inorder = FALSE, .combine = "c", .packages = "glmnet") %dorng% {
  fit.ii <- cv.glmnet(x = X.super.hyp, y)
  (fit.ii$lambda.min + fit.ii$lambda.1se) / 2
}
stopCluster(cl)

hist(lambdas.super.hyp)
lambda.best.super.hyp <- Mode(lambdas.super.hyp)

fit.super.hyp <- glmnet(x = X.super.hyp, y = y, lambda = lambda.best.super.hyp)
coef(fit.super.hyp)

## estimate test error ----

Xp.super.hyp <- w.super[!ids$is.analysis.group, hyps]  ## 'held out' / validation set
yp <- w.super[!ids$is.analysis.group, "stroop"]

## observed error

ys.super.hyp <- cbind(
  y = y,
  yhat = c(predict(fit.super.hyp, newx = Xp.super.hyp))
)

(cor.obs.super.hyp <- cor(ys.super.hyp)["y", "yhat"])
(mse.obs.super.hyp <- mean((ys.super.hyp[, "y"] - ys.super.hyp[, "yhat"])^2))

ys.super.hyp %>%
  as.data.frame %>%
  ggplot(aes(y, yhat)) +
  stat_boot_ci(n = 1e4, alpha = 0.5, color = "grey50") +
  geom_point()

## permuted validation error









## lasso ----

#
#
# ## estimate model ----
#
# X.super.hyp  <- w.super[ids$is.analysis.group, hyps]
# y  <- w.super[ids$is.analysis.group, "stroop"]
#
# ## find lambda and fit
#
# set.seed(0)
#
# cl <- makeCluster(n_cores - 1)
# registerDoParallel(cl)
# lambdas.super.hyp <- foreach(ii = seq_len(nreps), .inorder = FALSE, .combine = "c", .packages = "glmnet") %dorng% {
#   fit.ii <- cv.glmnet(x = X.super.hyp, y)
#   (fit.ii$lambda.min + fit.ii$lambda.1se) / 2
# }
# stopCluster(cl)
#
# hist(lambdas.super.hyp)
# lambda.best.super.hyp <- Mode(lambdas.super.hyp)
#
# fit.super.hyp <- glmnet(x = X.super.hyp, y = y, lambda = lambda.best.super.hyp)
# coef(fit.super.hyp)
#
# ## estimate test error ----
#
# Xp.super.hyp <- w.super[!ids$is.analysis.group, hyps]  ## 'held out' / validation set
# yp <- w.super[!ids$is.analysis.group, "stroop"]
#
# ## observed error
#
# ys.super.hyp <- cbind(
#   y = y,
#   yhat = c(predict(fit.super.hyp, newx = Xp.super.hyp))
#   )
#
# (cor.obs.super.hyp <- cor(ys.super.hyp)["y", "yhat"])
# (mse.obs.super.hyp <- mean((ys.super.hyp[, "y"] - ys.super.hyp[, "yhat"])^2))
#
# ys.super.hyp %>%
#   as.data.frame %>%
#   ggplot(aes(y, yhat)) +
#   stat_boot_ci(n = 1e4, alpha = 0.5, color = "grey50") +
#   geom_point()
#
# ## permuted validation error






<!-- ```{r md-core} -->
  
  <!-- names.md.core <- colnames(w.mmp)[grep(paste0(md$core, collapse = "|"), colnames(w.mmp))] -->
  <!-- if (length(names.md.core) != 66) stop("something wrong") -->
  
  <!-- ## estimate model ---- -->
  
  <!-- X.md.core <- w.mmp[ids$is.analysis.group, names.md.core] -->
  
  <!-- ## find lambda and fit -->
  
  <!-- set.seed(0) -->
  
  <!-- cl <- makeCluster(n_cores - 1) -->
  <!-- registerDoParallel(cl) -->
  <!-- lambdas.md.core <- foreach(ii = seq_len(nreps), .inorder = FALSE, .combine = "c", .packages = "glmnet") %dorng% { -->
      <!--   fit.ii <- cv.glmnet(x = X.md.core, y) -->
        <!--   # (fit.ii$lambda.min + fit.ii$lambda.1se) / 2 -->
        <!--   fit.ii$lambda.min -->
        <!-- } -->
  <!-- stopCluster(cl) -->
  
  <!-- hist(lambdas.md.core) -->
  <!-- lambda.best.md.core <- Mode(lambdas.md.core) -->
  
  <!-- fit.md.core <- glmnet(x = X.md.core, y = y, lambda = lambda.best.md.core) -->
  <!-- coef(fit.md.core) -->
  
  
  <!-- ## estimate test error ---- -->
  
  <!-- Xp.md.core <- w.mmp[!ids$is.analysis.group, names.md.core]  ## 'held out' / validation set -->

<!-- ## observed error -->
  
  <!-- ys.md.core <- cbind( -->
                              <!--   y = yp, -->
                              <!--   yhat = c(predict(fit.md.core, newx = Xp.md.core)) -->
                              <!--   ) -->
  
  <!-- (cor.obs.md.core <- cor(ys.md.core)["y", "yhat"]) -->
  <!-- (mse.obs.md.core <- mean((ys.md.core[, "y"] - ys.md.core[, "yhat"])^2)) -->
  
  <!-- ys.md.core %>% -->
  <!--   as.data.frame %>% -->
  <!--   ggplot(aes(y, yhat)) + -->
  <!--   stat_boot_ci(n = 1e4, alpha = 0.5, color = "grey50") + -->
  <!--   geom_point() -->
  
  <!-- ## permuted validation error -->
  
  <!-- cl <- makeCluster(n_cores - 1) -->
  <!-- registerDoParallel(cl) -->
  <!-- mse.perm.md.core <- foreach(ii = seq_len(nresamps), .inorder = FALSE, .combine = "c", .packages = "glmnet") %dorng% { -->
      
      <!--   yperm <- y[sample.int(length(y))] -->
        <!--   fit.ii <- glmnet(x = X.md.core, y = yperm, lambda = lambda.best.md.core) -->
          
          <!--   mean(((yp - predict(fit.ii, newx = Xp.md.core))^2)) -->
          
          <!-- } -->
  <!-- stopCluster(cl) -->
  
  <!-- (p.md.core <- sum(mse.perm.md.core < mse.obs.md.core) / nresamps) -->
  
  
  <!-- ``` -->
  
  
  
  <!-- ```{r} -->
  
  <!-- d <- data.frame( -->
                          
                          <!--   presma_R_incon = rowMeans(w.mmp[ids$is.analysis.group, c("SCEF_R_incongruency", "p32pr_R_incongruency")]), -->
                          <!--   presma_R_targ  = rowMeans(w.mmp[ids$is.analysis.group, c("SCEF_R_target", "p32pr_R_target")]), -->
                          
                          <!--   presma_L_incon = rowMeans(w.mmp[ids$is.analysis.group, c("SCEF_L_incongruency", "p32pr_L_incongruency")]), -->
                          <!--   presma_L_targ  = rowMeans(w.mmp[ids$is.analysis.group, c("SCEF_L_target", "p32pr_L_target")]), -->
                          
                          
                          <!--   mfc_R_incon = rowMeans(w.mmp[ids$is.analysis.group, c("8BM_R_incongruency", "a32pr_R_incongruency", "d32_R_incongruency")]), -->
                          <!--   mfc_R_targ = rowMeans(w.mmp[ids$is.analysis.group, c("8BM_R_target", "a32pr_R_target", "d32_R_target")]), -->
                          
                          <!--   mfc_L_incon = rowMeans(w.mmp[ids$is.analysis.group, c("8BM_L_incongruency", "a32pr_L_incongruency", "d32_L_incongruency")]), -->
                          <!--   mfc_L_targ = rowMeans(w.mmp[ids$is.analysis.group, c("8BM_L_target", "a32pr_L_target", "d32_L_target")]), -->
                          
                          <!--   w.super[ids$is.analysis.group, ] -->
                          
                          <!-- ) -->
  
  
  <!-- fit <- lmrob( -->
                       <!--   stroop ~   -->
                       <!--     rowMeans(w.mmp[ids$is.analysis.group, c("8BM_R_target", "p32pr_R_target")]) + -->
                       <!--     # rowMeans(w.mmp[ids$is.analysis.group, c("8BM_R_incongruency", "p32pr_R_incongruency")]) + -->
                       <!--     # I(presma_R_targ + presma_L_targ) + -->
                       <!--     # I(mfc_R_incon + mfc_L_incon) + -->
                       <!--     I(lppc_R_target + dlpfc_R_target + lppc_L_target + dlpfc_L_target),  -->
                       <!--   d -->
                       <!--   ) -->
  <!-- summary(fit) -->
  <!-- plot(fit) -->
  <!-- cor(d) -->
  
  
  
  <!-- ## models ---- -->
  
  <!-- y  <- w.mmp[ids$is.analysis.group, "stroop"] -->
  
  <!-- X.targpref  <- w.mmp[ids$is.analysis.group, c(paste0(rois.pref.targt.vs.distr, "_target"))] -->
  <!-- lambdas.targpref <- tune_lambda(X.targpref, y) -->
  <!-- hist(lambdas.targpref) -->
  <!-- lambda.best.targpref <- Mode(lambdas.targpref) -->
  <!-- fit.targpref <- glmnet(x = X.targpref, y = y, lambda = lambda.best.targpref) -->
  <!-- coef(fit.targpref) -->
  
  <!-- X.targ  <- w.mmp[ids$is.analysis.group, c(paste0(rois.targt, "_target"))] -->
  <!-- lambdas.targ <- tune_lambda(X.targ, y) -->
  <!-- hist(lambdas.targ) -->
  <!-- lambda.best.targ <- Mode(lambdas.targ) -->
  <!-- fit.targ <- glmnet(x = X.targ, y = y, lambda = lambda.best.targ) -->
  <!-- coef(fit.targ) -->
  
  <!-- X.incon  <- w.mmp[ids$is.analysis.group, c(paste0(rois.incon, "_incongruency"))] -->
  <!-- lambdas.incon <- tune_lambda(X.incon, y) -->
  <!-- hist(lambdas.incon) -->
  <!-- lambda.best.incon <- Mode(lambdas.incon) -->
  <!-- fit.incon <- glmnet(x = X.incon, y = y, lambda = lambda.best.incon) -->
  <!-- coef(fit.incon) -->
  
  <!-- X.distr  <- w.mmp[ids$is.analysis.group, c(paste0(rois.distr, "_distractor"))] -->
  <!-- lambdas.distr <- tune_lambda(X.distr, y) -->
  <!-- hist(lambdas.distr) -->
  <!-- lambda.best.distr <- Mode(lambdas.distr) -->
  <!-- fit.distr <- glmnet(x = X.distr, y = y, lambda = lambda.best.distr) -->
  <!-- coef(fit.distr) -->
  
  <!-- regs <- combo_paste(c(rois.pref.targt.vs.distr), c("target", "distractor", "incongruency")) -->
  <!-- X.all  <- w.mmp[ids$is.analysis.group, regs] -->
  <!-- lambdas.all <- tune_lambda(X.all, y) -->
  <!-- hist(lambdas.all) -->
  <!-- lambda.best.all <- Mode(lambdas.all) -->
  <!-- fit.all <- glmnet(x = X.all, y = y, lambda = lambda.best.all) -->
  <!-- coef(fit.all) -->
  
  
  <!-- ``` -->
  
  