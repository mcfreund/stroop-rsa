
source(here::here("code", "packages.R"))
source(here("code", "strings.R"))
source(here("code", "funs.R"))
source(here("code", "plots.R"))
source(here("code", "read_atlases.R"))


## get xmats, calc derivs, fit glms, fit rsa models ----

subjs <- c(subjs.analysis, subjs.validation)
movs <- paste0("movregs", 0:5)
derivs <- paste0(movs, "dt")

lt <- lower.tri(diag(16))
X_rsa <- fread(here("out", "rsa", "mods", "rsv_bias_lower-triangles.csv"))
X_rsa <- scale(as.matrix(X_rsa[, c("target", "distractor", "incongruency")]))

run.rsm <- as.matrix(fread(here("out", "rsa", "mods", "rsm_pro_bias_run.csv")), rownames = 1)
run.rsv <- mat2vec(run.rsm, value.name = "run.model")  ## unwrap run model to lower-tri vector
U <- cbind(1, run.rsv$run.model)  ## model: intercept, run regressor



X_glm <- setNames(vector("list", length(subjs)), subjs)
B_glm <- X_glm
R <- array(NA, dim = c(16, 16, length(subjs), 2))
stats.subjs.movregs <- X_glm
n.tr <- 1080

for (subj.i in seq_along(subjs)) {
  # subj.i = 1
  
  X_glm.i <- read_xmat(here("glms", subjs[subj.i], "results", "pro_bias_acc-only", "X.xmat.1D"))
  colnames(X_glm.i) <- gsub("\\[|\\]|#0|#", "", colnames(X_glm.i))
  X_glm.i <- X_glm.i[, c("Run1Pol", "Run2Pol", movs, bias.items)]
  
  ## calculate time derivs
  
  d <- rbind(
    rbind(0, diff(X_glm.i[1:(n.tr/2), movs])),
    rbind(0, diff(X_glm.i[(n.tr/2+1):n.tr, movs]))
  )
  colnames(d) <- derivs
  X_glm.i <- cbind(X_glm.i, d)
  
  X_glm[[subj.i]] <- X_glm.i  ## store
  
  ## get censored frames
  
  is.ok <- as.logical(fread(here("glms", subjs[subj.i], "input", "pro", "movregs_FD_mask.txt"))[[1]])
  
  ## glm
  
  B_glm <- coef(
    .lm.fit(
      X_glm.i[is.ok, c("Run1Pol", "Run2Pol", bias.items)], 
      y = X_glm.i[is.ok, c(movs, derivs)]
      )
    )[-(1:2), ]  ## remove intercepts
  
  ## rsa
  
  R06 <- cor(t(B_glm[, 1:6]))
  R12 <- cor(t(B_glm))
  
  R[, , subj.i, 1] <- R06
  R[, , subj.i, 2] <- R12
  
  Y <- scale(cbind(rank(R06[lt]), rank(R12[lt])))  ## rank transform and z-score normalize
  
  fits <- .lm.fit(U, Y)  ## regress run component from rsv
  B1 <- coef(fits)[2, ]  ## slopes
  Y_prw <- Y - U[, 2] %*% t(B1)    ## unscaled RSA response vectors
  
  B_rsa <- coef(.lm.fit(X_rsa, Y_prw))
  dimnames(B_rsa) <- list(model = c("target", "distractor", "incongruency"), nmovregs = c("m06", "m12"))
  stats.subjs.movregs[[subj.i]] <- reshape2::melt(B_rsa)
  
  
}

stats.subjs.movregs <- bind_rows(stats.subjs.movregs, .id = "subj")


## inferential stats ----

group.subjs.movregs <- stats.subjs.movregs %>%
  
  group_by(model, nmovregs) %>%
  
  summarize(
    b = mean(value),
    p = wilcox.test(value, alternative = "greater")$p.value
    )


group.subjs.movregs.delta <- stats.subjs.movregs %>%
  
  pivot_wider(names_from = nmovregs, values_from = value) %>%
  group_by(model) %>%
  summarize(
    b = mean(m12 - m06),
    p = wilcox.test(m06, m12, paired = TRUE)$p.value
  )

group.subjs.movregs
group.subjs.movregs.delta


stats.subjs.movregs %>%
  
  pivot_wider(names_from = nmovregs, values_from = value) %>%
  
  ggplot(aes(m06, m12)) +
  geom_point() +
  geom_abline() +
  facet_grid(cols = vars(model))

stats.subjs.movregs %>%
  
  ggplot(aes(model, value, color = nmovregs)) +
  geom_boxplot(width = 0.1) +
  scale_color_grey()

fwrite(group.subjs.movregs, here("out", "group", "movregs_rsa_vs0.csv"))
fwrite(group.subjs.movregs.delta, here("out", "group", "movregs_rsa_delta.csv"))


