fits.rank[[1]] -> fit


lm.beta(fit)


fit.raw <- lm(rank ~ target * incongruency, rsv)
summary(fit.raw)


rsv.scale <- lapply(rsv[c("rank", "target", "incongruency")], scale)
fit.scale <- lm(rank ~ target * incongruency, rsv.scale)

lm.beta(fit.raw)
coef(fit.scale)[-1]
coef(fit.raw)[-1]

cor(rsv$target * rsv$distractor, rsv$target)


lm.beta <- function(MOD) {
  ## assumes model has intercept and only one response
  b  <- coef(MOD)[-1]
  
  sx <- apply(model.matrix(MOD)[, -1], 2, sd)
  sy <- sd(MOD$model[[1]])
  
  b * sx / sy
  
}


