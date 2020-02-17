## about ----
## 


## setup ----

set.seed(0)

library(mikeutils)
library(magrittr)
library(here)
library(knitr)
library(dplyr)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(colorspace)
library(viridis)
library(nlme)
library(caret)
library(gtools)

theme_set(theme_bw(base_size = 12))

## read data

blups <- read.csv(here("out", "behav", "stroop_blups_rt_group201902.csv"), stringsAsFactors = FALSE)
stats.subjs.tdic <- fread(
  here("out", "rsa", "stats", paste0("subjs_pro_bias_acc-only_masks_pearson_residual_glm-tdic.csv"))
)

## subset and bind

stats.subjs.tdic <- stats.subjs.tdic[is.analysis.group == TRUE, ]  ## EXCLUDE HELD OUT SUBJECTS!
stats.subjs.tdic <- stats.subjs.tdic[y == "rank" & param %in% c("target", "distractor", "incongruency"), ]
stats.subjs.tdic <- stats.subjs.tdic[, c("coef", "y", "model") := NULL]  ## remove useless cols

stats.subjs.tdic <- full_join(blups, stats.subjs.tdic, by = "subj")

## format cols

stats.subjs.tdic <- cbind(
  stats.subjs.tdic,
  reshape2::colsplit(stats.subjs.tdic$roi, "_", c("roi.set", "superparcel"))
)
stats.subjs.tdic$superparcel[stats.subjs.tdic$superparcel == ""] <- "vwfa"



## (1) bivariate correlations ----

d <- stats.subjs.tdic %>% 
  filter(roi.set == "anatfunc" | superparcel %in% c("vwfa", "smmouth"))

# d %<>%
#   group_by(superparcel, param) %>%
#   mutate(
#     n = seq_len(n()),
#     is.out.line = n %in% outpro(partr, stroop)$out.id,
#     is.out.rank = n %in% outpro(rank(partr), rank(stroop))$out.id
#   )

topcors <- d %>%
  group_by(superparcel, param) %>%
  summarize(
    r.line = cor(beta, stroop),
    r.rank = cor(beta, stroop, method = "spearman"),
    r2.line = r.line^2,
    r2.rank = r.rank^2,
    # r.line.nout = sum(is.out.line),
    # r.rank.nout = sum(is.out.rank)
  ) %>%
  ungroup %>%
  filter(rank(-r2.line) < 21 & rank(-r2.rank) < 21) %>%
  arrange(-r2.line) %>%
  mutate(id = interaction(superparcel, param))

topcors

## plot

d %<>% 
  mutate(id = as.character(interaction(superparcel, param))) %>%
  filter(id %in% topcors$id)

d %>%
  group_by(id) %>%
  mutate(
    beta.z = scale(beta), 
    stroop.z = scale(stroop),
    r = cor(beta, stroop)
  ) %>%
  ggplot(aes(beta.z, stroop.z, fill = param)) +
  facet_grid(vars(param), vars(superparcel)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_boot_ci(fill = "grey50", alpha = 0.5) +
  geom_point(shape = 21, color = "white", size = 2)


## (2) model comparison ----

w <- d %>%
  ungroup %>%
  select(subj, congr, stroop, beta, id) %>%
  tidyr::spread(id, beta)

## build model combinations

varnames <- as.character(unique(d$id))
nvars <- length(varnames)
combs <- lapply(1:nvars, function(x) gtools::combinations(n = nvars, r = x))  ## row indicates columns of design matrix
nmods <- sum(sapply(combs, nrow))
mse <- vector("list", nmods)
X <- as.matrix(w[varnames])
strooprt <- w$stroop
mod <- mse

## build cross-validation folds

folds <- createFolds(strooprt, k = 10, list = TRUE, returnTrain = FALSE)

print(cbind(sapply(combs, nrow)))  ## number of models (inner loop) per combination (outer loop)

## estimate test error
starttime <- Sys.time()
a <- 1
for (comb.i in seq_along(combs)) {
  # comb.i = 1
  
  comb <- combs[[comb.i]]
  
  for (model.i in seq_len(nrow(comb))) {
    # model.i = 2
    
    ## build model
    
    regs.i <- comb[model.i, ]  ## indices for regressors
    Xi <- cbind(b0 = 1, X[, regs.i])
    colnames(Xi)[-1] <- varnames[regs.i]
    
    ## get response (folds)
    
    mse.i <- lapply(
      folds,
      function(., x, y) {
        xtrn <- x[-., ]
        xtst <- x[., ]
        ytrn <- y[-.]
        ytst <- y[.]
        fit  <- .lm.fit(xtrn, ytrn)
        yhat <- xtst %*% fit$coefficients
        
        mean((ytst - yhat)^2)
        
      },
      x = Xi,
      y = strooprt
    )
    
    mse[[a]] <- unlist(mse.i, use.names = FALSE)
    mod[[a]] <- varnames[regs.i]
    
    a <- a + 1
   
  }
  
  print(comb.i)
  
}
(endtime <- Sys.time() - starttime)

## examine: find 'best' model

nparams <- vapply(mod, length, numeric(1)) + 1

range(mse.bar)

mse.bar <- vapply(mse, mean, numeric(1))
m.best <- which.min(mse.bar)
se.best <- sd(mse[[m.best]]) / sqrt(length(mse[[m.best]]))

thresh <- mse.bar[m.best] + se.best
is.below <- mse.bar < thresh
is.parsim <- nparams <= min(nparams[is.below])
is.viable <- is.below & is.parsim

mod[m.best]
mod[is.viable]
mod[is.viable][which.min(mse.bar[is.viable])]

plot(nparams, mse.bar)