## about ----
## 


## setup ----

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
source(here("code", "strings.R"))

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


farout <- function(x) {
  
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr3 <- IQR(x) * 3
  
  x < (q1 - iqr3) | x > (q3 + iqr3)
  
}


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
set.seed(0)
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

## find 'best' model

mse.bar <- vapply(mse, mean, numeric(1))
m.best <- which.min(mse.bar)
se.best <- sd(mse[[m.best]]) / sqrt(length(mse[[m.best]]))
nparams <- vapply(mod, length, numeric(1)) + 1

range(mse.bar)
# plot(nparams, mse.bar)

thresh <- mse.bar[m.best] + se.best
is.below <- mse.bar < thresh
is.parsim <- nparams <= min(nparams[is.below])
is.viable <- is.below & is.parsim

mod[m.best]
mod[is.viable]
(mod.winning <- mod[is.viable][[which.min(mse.bar[is.viable])]])
mse[m.best] %>% unlist %>% mean / nrow(X)
lapply(mse[is.viable], function(.) mean(unlist(.)) / nrow(X))


## partial regression plot ----

# d <- d %>%
#   select(subj, congr, incon, superparcel, param, beta) %>%
#   as.data.table %>%
#   melt(value.name = "rt", id.vars = c("subj", "superparcel", "param", "beta"), measure.vars = c("congr", "incon"))
# 

fit <- readRDS(here("out", "behav", "fit1-het-trim_group201902.RDS"))
behav <- fit$data
w.scale <- data.frame(subj = w$subj, scale(w[mod.winning]))
behav %<>% full_join(w.scale, by = "subj")

w.resid <- data.frame(
  lppc_R.target       = lm(lppc_R.target ~ mfc_L.incongruency + vvis_L.incongruency, w.scale)$resid,
  mfc_L.incongruency  = lm(mfc_L.incongruency ~ vvis_L.incongruency + lppc_R.target, w.scale)$resid,
  vvis_L.incongruency = lm(vvis_L.incongruency ~ lppc_R.target + mfc_L.incongruency, w.scale)$resid
)

blups <- cbind(blups, w.resid)

cor(blups$stroop, blups$mfc_L.incongruency)
cor(blups$stroop, blups$lppc_R.target)
cor(blups$stroop, blups$vvis_L.incongruency)


## fit MLMs ----

## read and subset
# 
# stroop.pro <- read.csv(here("data", "behavior-and-events_group201902.csv"),  stringsAsFactors = FALSE) %>%
#   filter(session == "pro", is.analysis.group) %>%
#   mutate(trial.type = ifelse(trial.type == "i", "incon", "congr"))
# 
# ## initial fit  
# 
# stroop.pro.rt <- stroop.pro %>% filter(acc == 1, !is.na(rt), rt < 3000, rt > 250)
# is.weird.rt <- stroop.pro.rt$subj %in% c("849971", "161832") & stroop.pro.rt$rt < 500
# stroop.pro.rt <- stroop.pro.rt[!is.weird.rt, ]
# fit1.het <- lme(
#   rt ~ trial.type, 
#   random  = ~ trial.type | subj,
#   data    = stroop.pro.rt,
#   weights = varIdent(form = ~ 1 | subj),
#   control = lmeControl(maxIter = 1e5, msMaxIter = 1e5, niterEM = 1e5, msMaxEval = 1e5),
#   method  = "REML"
# )
# 
# ## trim and re-fit
# 
# stroop.pro.rt$resid.p <- resid(fit1.het, type = "p")
# stroop.pro.rt$is.far.out <- farout(stroop.pro.rt$resid.p)
# fit1.het.trim <- update(fit1.het, subset = !is.far.out)


fit.winning <- lme(
    rt ~ trial.type + lppc_R.target + mfc_L.incongruency + vvis_L.incongruency +
      trial.type:lppc_R.target +
      trial.type:mfc_L.incongruency +
      trial.type:vvis_L.incongruency,
    random  = ~ trial.type | subj,
    data    = behav,
    weights = varIdent(form = ~ 1 | subj),
    control = lmeControl(maxIter = 1e5, msMaxIter = 1e5, niterEM = 1e5, msMaxEval = 1e5),
    method  = "REML"
  )

summary(fit.winning)
intervals(fit.winning, which = "fixed")

## plot: MLM

blups$incon <- blups$congr + blups$stroop

blups.long <- full_join(
  blups %>%
    select(subj, congr, incon) %>%
    reshape2::melt(variable.name = "trial.type", value.name = "rt"),
  blups %>%
    select(subj, lppc_R.target:vvis_L.incongruency) %>%
    reshape2::melt(value.name = "beta"),
  by = "subj"
)

blups.long %>%
  ggplot(aes(beta, rt, color = trial.type)) +
  facet_grid(cols = vars(variable)) +
  geom_smooth(
    data = blups.long %>% filter(trial.type == "congr"), method = "lm", se = FALSE
  ) +
  geom_smooth(
    data = blups.long %>% filter(trial.type == "incon"), method = "lm", se = FALSE
  ) +
  # stat_boot_ci(
  #   data = blups.long %>% filter(trial.type == "incon"), alpha = 0.2
  # ) +
  # stat_boot_ci(
  #   data = blups.long %>% filter(trial.type == "congr"), alpha = 0.2
  # )
  geom_point(alpha = 0.2) +
  scale_color_brewer(type = "qual", palette = 2)
  # geom_line(aes(group = subj), color = "grey20")
blups %>%
  select(-congr, -incon) %>%
  reshape2::melt(id.var = c("subj", "stroop")) %>%
  ggplot(aes(value, stroop, color = variable)) +
  facet_grid(cols = vars(variable)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_boot_ci(alpha = 0.3) +
  # stat_boot_ci(
  #   data = blups.long %>% filter(trial.type == "congr"), alpha = 0.2
  # )
  geom_point() +
  scale_color_brewer(type = "qual", palette = 2)




## MDS ---- 

library(vegan)


split.str.item <- function(col.j, prefix = "") {
  ## takes a single "item" vector, e.g. "blueBLUE", and decomposes
  ## it into color ("blue") word ("BLUE"), congruency ("C"), and label 
  ## ("C.BLUE", for plotting). 
  # prefix <- paste0(col.j, ".")
  col.j      <- as.character(col.j)
  color      <- gsub("[A-Z]", "", col.j)
  word       <- gsub("[a-z]", "", col.j)
  congruency <- ifelse(color == tolower(word), "C", "I")
  label      <- paste0(congruency, ".", word)
  cols       <- as.data.frame(cbind(color, word, congruency, label))
  colnames(cols) <- paste0(prefix, c("color", "word", "congruency", "label"))
  return(cols)
}
mds.to.df <- function(mat) {
  mat %>%
    as.data.frame %>%
    tibble::rownames_to_column("stim") %>%
    bind_cols(., split.str.item(.$stim))
}
plot.mds <- function(df) {
  df %>%
    ggplot(aes(MDS1, MDS2)) +
    geom_label(aes(label = word, color = color), fill = "grey60", fontface = "bold", label.size = 0) +
    scale_color_manual(values = setNames(bias.colors, bias.colors)) +
    theme(
      panel.background = element_blank(), 
      axis.text = element_blank(), 
      legend.position = "none", 
      axis.ticks = element_blank()
    )
}

## get data

models <- c("lppc_R", "mfc_L", "vvis_L")
subjs <- blups$subj
n.subj <- nrow(blups)
n.model <- length(models)
n.stim <- length(bias.items)

R <- readRDS(here("out", "rsa", "obsv", "rsarray_pro_bias_acc-only_masks_pearson_residual-linear.rds"))
dimnames(R)
R <- R[, , blups$subj, c("anatfunc_lppc_R", "anatfunc_mfc_L", "anatfunc_vvis_L")]
R.lppc <- R[, , , "anatfunc_lppc_R"]
R.mfc <- R[, , , "anatfunc_mfc_L"]
R.vvc <- R[, , , "anatfunc_vvis_L"]

## lppc ----

M.lppc <- array(
  NA,
  dim = c(stim = n.stim, dim = 2, subj = n.subj),
  dimnames = list(stim = bias.items, dim = c("MDS1", "MDS2"), subj = blups$subj)
)
D.lppc <- 1 - R.lppc

for (subj.i in seq_len(n.subj)) M.lppc[, , subj.i] <- vegan::metaMDS(D.lppc[, , subj.i], k = 2, trace = FALSE)$points

R.bar.lppc <- apply(R.lppc, 1:2, function(x) tanh(mean(atanh(x))))
D.bar.lppc <- 1 - R.bar.lppc
M.bar.lppc <- vegan::metaMDS(D.bar.lppc, k = 2)$points
M.bar.lppc %>% mds.to.df %>% plot.mds


## align

M.rot.lppc <- apply(
  M.lppc, 3, function(.x) vegan::procrustes(M.bar.lppc, .x, scale = FALSE)$Yrot
  )
M.rot.lppc <- array(
  M.rot.lppc, c(n.stim, 2, n.subj), 
  dimnames = list(stim = bias.items, d = c("MDS1", "MDS2"), subj = subjs)
  )

## get aligned centroids


rot.lppc <- reshape2::melt(M.rot.lppc)
rot.lppc %<>% cbind(split.str.item(rot.lppc$stim))

centroids.rot.lppc <- rot.lppc %>%
  group_by(color, d, subj) %>%
  summarize(value = mean(value)) %>%
  tidyr::spread(d, value) %>%
  full_join(blups, by = "subj") %>%
  mutate(stroop.z = scale(stroop))
centroids.rot.lppc %>%
  ggplot(aes(MDS1, MDS2, fill = color)) +
  geom_point(shape = 21, size = 2) +
  scale_fill_identity()


## regression

library(robustbase)

l.lppc <- split(centroids.rot.lppc, centroids.rot.lppc$color)

yhat.lppc <- lapply(
  l.lppc,
  function(x) {
    # b.d1 <- coef(lmrob(MDS1 ~ stroop.z, x))
    # b.d2 <- coef(lmrob(MDS2 ~ stroop.z, x))
    b.d1 <- coef(lm(MDS1 ~ stroop.z, x))
    b.d2 <- coef(lm(MDS2 ~ stroop.z, x))
    data.frame(
      mu.MDS1   = b.d1[["(Intercept)"]],
      mu.MDS2   = b.d2[["(Intercept)"]],
      yhat.MDS1 = b.d1[["(Intercept)"]] - b.d1[["stroop.z"]],
      yhat.MDS2 = b.d2[["(Intercept)"]] - b.d2[["stroop.z"]]
    )
  }
  ) %>%
  do.call(rbind, .) %>%
  tibble::rownames_to_column("color")


## plot

centroids.rot.lppc %>%
  ggplot(aes(MDS1, MDS2, fill = color)) +
  # geom_point(shape = 21, size = 2, color = "white", alpha = 0.5) +
  scale_fill_manual(values = c(blue = "blue", white = "grey40", red = "firebrick", purple = "purple")) +
  scale_color_manual(values = c(blue = "blue", white = "grey40", red = "firebrick", purple = "purple")) +
  geom_segment(
    data = yhat.lppc, aes(x = mu.MDS1, xend = yhat.MDS1, y = mu.MDS2, yend = yhat.MDS2),
    size = 2,
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_point(data = yhat.lppc, aes(x = mu.MDS1, y = mu.MDS2, fill = color), size = 5, shape = 21, color = "black")




## mfc ----

M.mfc <- array(
  NA,
  dim = c(stim = n.stim, dim = 2, subj = n.subj),
  dimnames = list(stim = bias.items, dim = c("MDS1", "MDS2"), subj = blups$subj)
)
D.mfc <- 1 - R.mfc

for (subj.i in seq_len(n.subj)) M.mfc[, , subj.i] <- vegan::metaMDS(D.mfc[, , subj.i], k = 2, trace = FALSE)$points

R.bar.mfc <- apply(R.mfc, 1:2, function(x) tanh(mean(atanh(x))))
D.bar.mfc <- 1 - R.bar.mfc
M.bar.mfc <- vegan::metaMDS(D.bar.mfc, k = 2)$points
M.bar.mfc %>% mds.to.df %>% plot.mds


## align

M.rot.mfc <- apply(
  M.mfc, 3, function(.x) vegan::procrustes(M.bar.mfc, .x, scale = FALSE)$Yrot
)
M.rot.mfc <- array(
  M.rot.mfc, c(n.stim, 2, n.subj), 
  dimnames = list(stim = bias.items, d = c("MDS1", "MDS2"), subj = subjs)
)

## get aligned centroids


rot.mfc <- reshape2::melt(M.rot.mfc)
rot.mfc %<>% cbind(split.str.item(rot.mfc$stim))

centroids.rot.mfc <- rot.mfc %>%
  group_by(color, congruency, d, subj) %>%
  summarize(value = mean(value)) %>%
  tidyr::spread(d, value) %>%
  full_join(blups, by = "subj") %>%
  mutate(stroop.z = scale(stroop))
centroids.rot.mfc %>%
  ggplot(aes(MDS1, MDS2, color = color)) +
  geom_point(aes(shape = congruency)) +
  scale_color_brewer(type = "qual", palette = 2)


## regression

# centroids.rot.mfc
# lm()
stroop.resid <- scale(lm(stroop ~ lppc_R.target + vvis_L.incongruency, w)$resid)
stroop.resid <- data.frame(stroop.resid, subj = w$subj)
centroids.rot.mfc %<>% full_join(stroop.resid, by = "subj")


l.mfc <- split(centroids.rot.mfc, interaction(centroids.rot.mfc$congruency, centroids.rot.mfc$color))
# l.mfc[[1]] -> x
yhat.mfc <- lapply(
  l.mfc,
  function(x) {
    # b.d1 <- coef(lmrob(MDS1 ~ scale(stroop.resid), x))
    # b.d2 <- coef(lmrob(MDS2 ~ scale(stroop.resid), x))
    # b.d1 <- coef(lm(MDS1 ~ scale(stroop.resid), x))
    # b.d2 <- coef(lm(MDS2 ~ scale(stroop.resid), x))
    b.d1 <- coef(lm(MDS1 ~ stroop.z, x))
    b.d2 <- coef(lm(MDS2 ~ stroop.z, x))
    # b.d1 <- coef(lmrob(MDS1 ~ stroop.z, x))
    # b.d2 <- coef(lmrob(MDS2 ~ stroop.z, x))
        
    data.frame(
      mu.MDS1   = b.d1[["(Intercept)"]],
      mu.MDS2   = b.d2[["(Intercept)"]],
      # yhat.MDS1 = b.d1[["(Intercept)"]] - b.d1[["scale(stroop.resid)"]],
      # yhat.MDS2 = b.d2[["(Intercept)"]] - b.d2[["scale(stroop.resid)"]]
      yhat.MDS1 = b.d1[["(Intercept)"]] - b.d1[["stroop.z"]],
      yhat.MDS2 = b.d2[["(Intercept)"]] - b.d2[["stroop.z"]]
    )
  }
  ) %>%
  do.call(rbind, .) %>%
  tibble::rownames_to_column("congruency.color") %>%
  mutate(
    congruency = substr(congruency.color, 1, 1),
    color = substr(congruency.color, 3, nchar(congruency.color))
  )


## plot

centroids.rot.mfc %>%
  mutate(congruency.color = interaction(congruency, color)) %>%
  ggplot(aes(MDS1, MDS2, fill = color)) +
  # geom_point(aes(shape = congruency, color = color), size = 3, alpha = 0.2) +
  scale_fill_manual(values = c(blue = "blue", white = "grey40", red = "firebrick", purple = "purple")) +
  scale_color_manual(values = c(blue = "blue", white = "grey40", red = "firebrick", purple = "purple")) +
  geom_segment(
    data = yhat.mfc, aes(x = mu.MDS1, xend = yhat.MDS1, y = mu.MDS2, yend = yhat.MDS2),
    size = 2,
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_point(data = yhat.mfc, aes(x = mu.MDS1, y = mu.MDS2, color = color, shape = congruency), size = 5)



## vvc ----


M.vvc <- array(
  NA,
  dim = c(stim = n.stim, dim = 2, subj = n.subj),
  dimnames = list(stim = bias.items, dim = c("MDS1", "MDS2"), subj = blups$subj)
)
D.vvc <- 1 - R.vvc

for (subj.i in seq_len(n.subj)) M.vvc[, , subj.i] <- vegan::metaMDS(D.vvc[, , subj.i], k = 2, trace = FALSE)$points

R.bar.vvc <- apply(R.vvc, 1:2, function(x) tanh(mean(atanh(x))))
D.bar.vvc <- 1 - R.bar.vvc
M.bar.vvc <- vegan::metaMDS(D.bar.vvc, k = 2)$points
M.bar.vvc %>% mds.to.df %>% plot.mds


## align

M.rot.vvc <- apply(
  M.vvc, 3, function(.x) vegan::procrustes(M.bar.vvc, .x, scale = TRUE)$Yrot
)
M.rot.vvc <- array(
  M.rot.vvc, c(n.stim, 2, n.subj), 
  dimnames = list(stim = bias.items, d = c("MDS1", "MDS2"), subj = subjs)
)

## get aligned centroids

rot.vvc <- reshape2::melt(M.rot.vvc)
rot.vvc %<>% cbind(split.str.item(rot.vvc$stim))

centroids.rot.vvc <- rot.vvc %>%
  group_by(color, d, subj) %>%
  summarize(value = mean(value)) %>%
  tidyr::spread(d, value) %>%
  full_join(blups, by = "subj") %>%
  mutate(stroop.z = scale(stroop))
centroids.rot.vvc %>%
  ggplot(aes(MDS1, MDS2, fill = color)) +
  geom_point(shape = 21, size = 2) +
  scale_fill_identity()
 

R.bar.vvc <- apply(R.vvc, 1:2, function(x) tanh(mean(atanh(x))))
D.bar.vvc <- 1 - R.bar.vvc
M.bar.vvc <- vegan::metaMDS(D.bar.vvc, k = 2)$points
M.bar.vvc %>% mds.to.df %>% plot.mds


## align

M.rot.vvc <- apply(
  M.vvc, 3, function(.x) vegan::procrustes(M.bar.vvc, .x, scale = FALSE)$Yrot
)
M.rot.vvc <- array(
  M.rot.vvc, c(n.stim, 2, n.subj), 
  dimnames = list(stim = bias.items, d = c("MDS1", "MDS2"), subj = subjs)
)

## get aligned centroids


rot.vvc <- reshape2::melt(M.rot.vvc)
rot.vvc %<>% cbind(split.str.item(rot.vvc$stim))

centroids.rot.vvc <- rot.vvc %>%
  group_by(color, d, subj) %>%
  summarize(value = mean(value)) %>%
  tidyr::spread(d, value) %>%
  full_join(blups, by = "subj") %>%
  mutate(stroop.z = scale(stroop))
centroids.rot.vvc %>%
  ggplot(aes(MDS1, MDS2, fill = color)) +
  geom_point(shape = 21, size = 2) +
  scale_fill_identity()


## regression

l.vvc <- split(centroids.rot.vvc, centroids.rot.vvc$color)

yhat.vvc <- lapply(
  l.vvc,
  function(x) {
    # b.d1 <- coef(lmrob(MDS1 ~ stroop.z, x))
    # b.d2 <- coef(lmrob(MDS2 ~ stroop.z, x))
    b.d1 <- coef(lm(MDS1 ~ stroop.z, x))
    b.d2 <- coef(lm(MDS2 ~ stroop.z, x))
    data.frame(
      mu.MDS1   = b.d1[["(Intercept)"]],
      mu.MDS2   = b.d2[["(Intercept)"]],
      yhat.MDS1 = b.d1[["(Intercept)"]] - b.d1[["stroop.z"]],
      yhat.MDS2 = b.d2[["(Intercept)"]] - b.d2[["stroop.z"]]
    )
  }
) %>%
  do.call(rbind, .) %>%
  tibble::rownames_to_column("color")


## plot

centroids.rot.vvc %>%
  ggplot(aes(MDS1, MDS2, fill = color)) +
  geom_point(shape = 21, size = 2, color = "white", alpha = 0.5) +
  scale_fill_manual(values = c(blue = "blue", white = "grey40", red = "firebrick", purple = "purple")) +
  scale_color_manual(values = c(blue = "blue", white = "grey40", red = "firebrick", purple = "purple")) +
  geom_segment(
    data = yhat.vvc, aes(x = mu.MDS1, xend = yhat.MDS1, y = mu.MDS2, yend = yhat.MDS2),
    size = 2,
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_point(data = yhat.vvc, aes(x = mu.MDS1, y = mu.MDS2, fill = color), size = 5, shape = 21, color = "black")




### MDS/regression by items (not centroids) ----

## lppc

rot.lppc <- reshape2::melt(M.rot.lppc)
rot.lppc %<>% cbind(split.str.item(rot.lppc$stim)) %>% tidyr::spread(d, value) %>%
  full_join(blups, by = "subj") %>%
  mutate(stroop.z = scale(stroop))


## regression

l.lppc <- split(rot.lppc, rot.lppc$stim)

yhat.lppc <- lapply(
  l.lppc,
  function(x) {
    # b.d1 <- coef(lmrob(MDS1 ~ stroop.z, x))
    # b.d2 <- coef(lmrob(MDS2 ~ stroop.z, x))
    b.d1 <- coef(lm(MDS1 ~ stroop.z, x))
    b.d2 <- coef(lm(MDS2 ~ stroop.z, x))
    data.frame(
      mu.MDS1 = b.d1[["(Intercept)"]],
      mu.MDS2 = b.d2[["(Intercept)"]],
      minus1sd.MDS1 = b.d1[["(Intercept)"]] + b.d1[["stroop.z"]],
      minus1sd.MDS2 = b.d2[["(Intercept)"]] + b.d2[["stroop.z"]],
      plus1sd.MDS1  = b.d1[["(Intercept)"]] - b.d1[["stroop.z"]],
      plus1sd.MDS2  = b.d2[["(Intercept)"]] - b.d2[["stroop.z"]]
    )
  }
  ) %>%
  do.call(rbind, .) %>%
  tibble::rownames_to_column("stim")

yhat.lppc %<>% cbind(split.str.item(.$stim))

## plot

yhat.lppc %>%
  ggplot(aes(x = minus1sd.MDS1, y = minus1sd.MDS2, fill = color)) +
  scale_fill_manual(values = c(blue = "blue", white = "grey40", red = "firebrick", purple = "purple")) +
  scale_color_manual(values = c(blue = "blue", white = "grey40", red = "firebrick", purple = "purple")) +
  geom_segment(
    aes(xend = plus1sd.MDS1, yend = plus1sd.MDS2, color = color),
    size = 2,
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_label(
    data = yhat.lppc,
    aes(label = word, color = color), fill = "grey60", fontface = "bold", label.size = 0
    )
  # geom_point(
  #   data = rot.lppc %>% group_by(color) %>% summarize(MDS1 = mean(MDS1), MDS2 = mean(MDS2)),
  #   aes(MDS1, MDS2, fill = color), size = 5, shape = 21, color = "black"
  #   )



## mfc

rot.mfc <- reshape2::melt(M.rot.mfc)
rot.mfc %<>% cbind(split.str.item(rot.mfc$stim)) %>% tidyr::spread(d, value) %>%
  full_join(blups, by = "subj") %>%
  mutate(stroop.z = scale(stroop))


## regression

stroop.resid <- scale(lm(stroop ~ lppc_R.target + vvis_L.incongruency, w)$resid)
stroop.resid <- data.frame(stroop.resid, subj = w$subj)
rot.mfc %<>% full_join(stroop.resid, by = "subj")
l.mfc <- split(rot.mfc, rot.mfc$stim)

yhat.mfc <- lapply(
  l.mfc,
  function(x) {
    # b.d1 <- coef(lmrob(MDS1 ~ stroop.z, x))
    # b.d2 <- coef(lmrob(MDS2 ~ stroop.z, x))
    # b.d1 <- coef(lm(MDS1 ~ stroop.z, x))
    # b.d2 <- coef(lm(MDS2 ~ stroop.z, x))
    b.d1 <- coef(lmrob(MDS1 ~ stroop.resid, x))
    b.d2 <- coef(lmrob(MDS2 ~ stroop.resid, x))
    data.frame(
      mu.MDS1 = b.d1[["(Intercept)"]],
      mu.MDS2 = b.d2[["(Intercept)"]],
      minus1sd.MDS1 = b.d1[["(Intercept)"]] + b.d1[[2]],
      minus1sd.MDS2 = b.d2[["(Intercept)"]] + b.d2[[2]],
      plus1sd.MDS1  = b.d1[["(Intercept)"]] - b.d1[[2]],
      plus1sd.MDS2  = b.d2[["(Intercept)"]] - b.d2[[2]]
    )
  }
) %>%
  do.call(rbind, .) %>%
  tibble::rownames_to_column("stim")

yhat.mfc %<>% cbind(split.str.item(.$stim))

## plot

yhat.mfc %>%
  ggplot(aes(x = minus1sd.MDS1, y = minus1sd.MDS2, fill = color)) +
  scale_fill_manual(values = c(blue = "blue", white = "grey40", red = "firebrick", purple = "purple")) +
  scale_color_manual(values = c(blue = "blue", white = "grey40", red = "firebrick", purple = "purple")) +
  geom_segment(
    data = yhat.mfc,
    aes(xend = plus1sd.MDS1, yend = plus1sd.MDS2, color = color),
    size = 2,
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_label(
    data = yhat.mfc %>% filter(congruency == "I"),
    aes(label = word, color = color), fill = "grey60", fontface = "bold", label.size = 0
  ) +
  geom_label(
    data = yhat.mfc %>% filter(congruency == "C"),
    aes(x = mu.MDS1, y = mu.MDS2, label = word, color = color), fill = "grey60", fontface = "bold", label.size = 0
  )

# geom_point(
#   data = rot.mfc %>% group_by(color) %>% summarize(MDS1 = mean(MDS1), MDS2 = mean(MDS2)),
#   aes(MDS1, MDS2, fill = color), size = 5, shape = 21, color = "black"
#   )



## quantile split ----

topthird <- w$subj[w$stroop < quantile(w$stroop, 1/3)]
botthird <- w$subj[w$stroop > quantile(w$stroop, 2/3)]

## lppc

D.top.lppc <- 1 - apply(R.lppc[, , topthird], 1:2, function(x) tanh(mean(atanh(x))))
D.bot.lppc <- 1 - apply(R.lppc[, , botthird], 1:2, function(x) tanh(mean(atanh(x))))
M.top.lppc <- vegan::metaMDS(D.top.lppc, k = 2, trace = FALSE)$points
M.bot.lppc <- vegan::metaMDS(D.bot.lppc, k = 2, trace = FALSE)$points

M.split.rot.lppc <- rbind(
  vegan::procrustes(M.bar.lppc, M.top.lppc, scale = FALSE)$Yrot %>% reshape2::melt() %>% mutate(split = "topthird"),
  vegan::procrustes(M.bar.lppc, M.bot.lppc, scale = FALSE)$Yrot %>% reshape2::melt() %>% mutate(split = "botthird")
)
M.split.rot.lppc %<>% 
  tidyr::spread(Var2, value) %>%
  rename(stim = Var1, MDS1 = "1", MDS2 = "2") %>%
  cbind(split.str.item(.$stim))
M.split.rot.lppc %<>% 
  as.data.table %>% 
  dcast(stim + color + word + congruency + label ~ split, value.var = c("MDS1", "MDS2"))
M.split.rot.lppc %>%
  ggplot(aes(x = MDS1_botthird, y = MDS2_botthird, fill = color)) +
  scale_fill_manual(values = c(blue = "blue", white = "grey40", red = "firebrick", purple = "purple")) +
  scale_color_manual(values = c(blue = "blue", white = "grey40", red = "firebrick", purple = "purple")) +
  geom_segment(
    aes(xend = MDS1_topthird, yend = MDS2_topthird, color = color),
    size = 2,
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_label(
    aes(label = word, color = color), fill = "grey60", fontface = "bold", label.size = 0
  )


## mfc

D.top.mfc <- 1 - apply(R.mfc[, , topthird], 1:2, function(x) tanh(mean(atanh(x))))
D.bot.mfc <- 1 - apply(R.mfc[, , botthird], 1:2, function(x) tanh(mean(atanh(x))))
M.top.mfc <- vegan::metaMDS(D.top.mfc, k = 2, trace = FALSE)$points
M.bot.mfc <- vegan::metaMDS(D.bot.mfc, k = 2, trace = FALSE)$points

M.split.rot.mfc <- rbind(
  vegan::procrustes(M.bar.mfc, M.top.mfc, scale = FALSE)$Yrot %>% reshape2::melt() %>% mutate(split = "topthird"),
  vegan::procrustes(M.bar.mfc, M.bot.mfc, scale = FALSE)$Yrot %>% reshape2::melt() %>% mutate(split = "botthird")
)
M.split.rot.mfc %<>% 
  tidyr::spread(Var2, value) %>%
  rename(stim = Var1, MDS1 = "1", MDS2 = "2") %>%
  cbind(split.str.item(.$stim))
# M.split.rot.mfc %<>%
#   as.data.table %>%
#   dcast(stim + color + word + congruency + label ~ split, value.var = c("MDS1", "MDS2"))
# M.split.rot.mfc %>%
#   filter(congruency == "I") %>%
#   ggplot(aes(x = MDS1_botthird, y = MDS2_botthird, fill = color)) +
#   scale_fill_manual(values = c(blue = "blue", white = "grey40", red = "firebrick", purple = "purple")) +
#   scale_color_manual(values = c(blue = "blue", white = "grey40", red = "firebrick", purple = "purple")) +
#   # geom_segment(
#   #   aes(xend = MDS1_topthird, yend = MDS2_topthird, color = color),
#   #   size = 2,
#   #   arrow = arrow(length = unit(0.25, "cm"))
#   # ) +
#   geom_line(aes(group = color), size = 2, color = "grey80") +
#   geom_label(
#     aes(label = word, color = color), fill = "grey80", fontface = "bold", label.size = 0
#   ) +
#   geom_line(aes(x = MDS1_topthird, y = MDS2_topthird, group = color), size = 2, color = "black") +
#   geom_label(
#     aes(x = MDS1_topthird, y = MDS2_topthird, label = word, color = color), fill = "black", fontface = "bold", label.size = 0
#   )

M.split.rot.mfc %>%
  # filter(congruency == "I") %>%
  ggplot(aes(x = MDS1, y = MDS2, fill = color)) +
  scale_fill_manual(values = c(blue = "blue", white = "grey40", red = "firebrick", purple = "purple")) +
  scale_color_manual(values = c(blue = "blue", white = "grey40", red = "firebrick", purple = "purple")) +
  # geom_segment(
  #   aes(xend = MDS1_topthird, yend = MDS2_topthird, color = color),
  #   size = 2,
  #   arrow = arrow(length = unit(0.25, "cm"))
  # ) +
  # geom_line(aes(group = color), size = 2, color = "grey80") +
  geom_label(
    aes(label = word, color = color), fill = "grey80", fontface = "bold", label.size = 0
  ) +
  facet_grid(cols = vars(split))
  # geom_line(aes(x = MDS1_topthird, y = MDS2_topthird, group = color), size = 2, color = "black") +
  geom_label(
    aes(x = MDS1_topthird, y = MDS2_topthird, label = word, color = color), fill = "black", fontface = "bold", label.size = 0
  )



## trichotomized bootstrap ----

topthird <- w$subj[w$stroop < quantile(w$stroop, 1/3)]
midthird <- w$subj[w$stroop < quantile(w$stroop, 2/3) & w$stroop > quantile(w$stroop, 1/3)]
botthird <- w$subj[w$stroop > quantile(w$stroop, 2/3)]


n.resamples <- 1E3
boot.lppc <- vector("list", n.resamples)

for (sample.i in seq_along(boot.lppc)) {
  
  ## get resamples
  
  topthird.i <- sample(topthird, replace = TRUE)
  midthird.i <- sample(midthird, replace = TRUE)
  botthird.i <- sample(botthird, replace = TRUE)
  
  ## get stats 
  
  D.top.lppc <- 1 - apply(R.lppc[, , topthird.i], 1:2, function(x) tanh(mean(atanh(x))))
  D.mid.lppc <- 1 - apply(R.lppc[, , midthird.i], 1:2, function(x) tanh(mean(atanh(x))))
  D.bot.lppc <- 1 - apply(R.lppc[, , botthird.i], 1:2, function(x) tanh(mean(atanh(x))))
  
  M.top.lppc <- vegan::metaMDS(D.top.lppc, k = 2, trace = FALSE)$points
  M.mid.lppc <- vegan::metaMDS(D.mid.lppc, k = 2, trace = FALSE)$points
  M.bot.lppc <- vegan::metaMDS(D.bot.lppc, k = 2, trace = FALSE)$points
  
  M.split.rot.lppc <- rbind(
    M.mid.lppc %>% vegan::procrustes(M.bar.lppc, ., scale = TRUE) %>% .$Yrot %>% 
      data.frame(split = "mid") %>% 
      rename(MDS1 = "X1", MDS2 = "X2") %>%
      tibble::rownames_to_column("stim") %>%
      bind_cols(., split.str.item(.$stim)),
    M.top.lppc %>% vegan::procrustes(M.bar.lppc, ., scale = TRUE) %>% .$Yrot %>% 
      data.frame(split = "top") %>% 
      rename(MDS1 = "X1", MDS2 = "X2") %>%
      tibble::rownames_to_column("stim") %>%
      bind_cols(., split.str.item(.$stim)),
    M.bot.lppc %>% vegan::procrustes(M.bar.lppc, ., scale = TRUE) %>% .$Yrot %>% 
      data.frame(split = "bot") %>% 
      rename(MDS1 = "X1", MDS2 = "X2") %>%
      tibble::rownames_to_column("stim") %>%
      bind_cols(., split.str.item(.$stim))
  )
  
  M.split.rot.lppc$iter <- sample.i
  
  boot.lppc[[sample.i]] <- M.split.rot.lppc
  
}

boot.lppc.df <- do.call(rbind, boot.lppc)
boot.lppc.df %>%
  mutate(
    congruency = as.character(congruency),
    color = as.character(color),
    split = relevel(split, "top")
  ) %>%
  group_by(iter, color, split) %>%
  summarize(MDS1 = mean(MDS1), MDS2 = mean(MDS2)) %>%
  # filter(congruency == "I") %>%
  ggplot(aes(x = MDS1, y = MDS2, color = color)) +
  geom_point(alpha = 0.1) +
  scale_color_manual(values = c(blue = "blue", white = "grey40", red = "firebrick", purple = "purple")) +
  stat_ellipse(type = "norm", size = 1) +
  facet_grid(cols = vars(split))
  # geom_label(
  #   data = M.bar.lppc %>%
  #     as.data.frame %>%
  #     tibble::rownames_to_column("stim") %>%
  #     bind_cols(., split.str.item(.$stim)),
  #     # filter(congruency == "C"),
    # aes(label = word, color = color), fill = "grey80", fontface = "bold", label.size = 0
  # )
  # geom_line(aes(x = MDS1_topthird, y = MDS2_topthird, group = color), size = 2, color = "black") +
  # geom_label(
  #   aes(x = MDS1_topthird, y = MDS2_topthird, label = word, color = color), fill = "black", fontface = "bold", label.size = 0
  # )



boot.mfc <- vector("list", n.resamples)
# topthird <- stroop.resid$subj[stroop.resid$stroop.resid < quantile(stroop.resid$stroop.resid, 1/3)]
# midthird <- stroop.resid$subj[
#   stroop.resid$stroop < quantile(stroop.resid$stroop, 2/3) & stroop.resid$stroop > quantile(stroop.resid$stroop, 1/3)
#   ]
# botthird <- stroop.resid$subj[stroop.resid$stroop > quantile(stroop.resid$stroop, 2/3)]

for (sample.i in seq_along(boot.mfc)) {
  
  ## get resamples
  
  topthird.i <- sample(topthird, replace = TRUE)
  midthird.i <- sample(midthird, replace = TRUE)
  botthird.i <- sample(botthird, replace = TRUE)
  
  ## get stats 
  
  D.top.mfc <- 1 - apply(R.mfc[, , topthird.i], 1:2, function(x) tanh(mean(atanh(x))))
  D.mid.mfc <- 1 - apply(R.mfc[, , midthird.i], 1:2, function(x) tanh(mean(atanh(x))))
  D.bot.mfc <- 1 - apply(R.mfc[, , botthird.i], 1:2, function(x) tanh(mean(atanh(x))))
  
  M.top.mfc <- vegan::metaMDS(D.top.mfc, k = 2, trace = FALSE)$points
  M.mid.mfc <- vegan::metaMDS(D.mid.mfc, k = 2, trace = FALSE)$points
  M.bot.mfc <- vegan::metaMDS(D.bot.mfc, k = 2, trace = FALSE)$points
  
  M.split.rot.mfc <- rbind(
    M.mid.mfc %>% vegan::procrustes(M.bar.mfc, ., scale = TRUE) %>% .$Yrot %>% 
      data.frame(split = "mid") %>% 
      rename(MDS1 = "X1", MDS2 = "X2") %>%
      tibble::rownames_to_column("stim") %>%
      bind_cols(., split.str.item(.$stim)),
    M.top.mfc %>% vegan::procrustes(M.bar.mfc, ., scale = TRUE) %>% .$Yrot %>% 
      data.frame(split = "top") %>% 
      rename(MDS1 = "X1", MDS2 = "X2") %>%
      tibble::rownames_to_column("stim") %>%
      bind_cols(., split.str.item(.$stim)),
    M.bot.mfc %>% vegan::procrustes(M.bar.mfc, ., scale = TRUE) %>% .$Yrot %>% 
      data.frame(split = "bot") %>% 
      rename(MDS1 = "X1", MDS2 = "X2") %>%
      tibble::rownames_to_column("stim") %>%
      bind_cols(., split.str.item(.$stim))
  )
  
  M.split.rot.mfc$iter <- sample.i
  
  boot.mfc[[sample.i]] <- M.split.rot.mfc
  
}

boot.mfc.df <- do.call(rbind, boot.mfc)


boot.mfc.df %>%
  mutate(
    congruency = as.character(congruency),
    color.congruency = as.character(interaction(color, congruency)),
    color.congruency = ifelse(congruency == "C", "C", color.congruency),
    tri = relevel(split, "top")
  ) %>%
  group_by(iter, congruency, tri) %>%
  # mutate(
  #   MDS1.c = MDS1 - mean(MDS1),
  #   MDS2.c = MDS2 - mean(MDS2)
  # ) %>%
  summarize(MDS1 = mean(MDS1), MDS2 = mean(MDS2)) %>%
  # filter(congruency == "C") %>%
  ggplot(aes(x = MDS1, y = MDS2, color = congruency)) +
  geom_point(alpha = 0.1, aes(fill = congruency, shape = congruency)) +
  scale_shape_manual(values = c(C = 21, I = 4)) +
  scale_fill_brewer(type = "qual", palette = 2) +
  scale_color_brewer(type = "qual", palette = 2) +
  # scale_fill_manual(values = c(C = "grey50", I = "black")) +
  # scale_color_manual(
  #   values = c(blue.C = "blue", white.C = "grey40", red.C = "firebrick", purple.C = "purple", I = "black")
  #   ) +
  stat_ellipse(type = "norm", size = 1) +
  facet_grid(cols = vars(tri))
