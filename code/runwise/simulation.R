## setup ----

source(here::here("code", "packages.R"))
source(here("code", "strings.R"))
source(here("code", "funs.R"))
source(here("code", "plots.R"))

# https://nwfsc-timeseries.github.io/atsa-labs/sec-tslab-random-walks-rw.html


## functions

ar1_cor <- function(n, rho) {
  # https://www.r-bloggers.com/2020/02/generating-correlation-matrix-for-ar1-model/
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
  rho^exponent
}



sim_null <- function(
  n_sim,
  S_tr, 
  n_cores = n.cores / 2,
  run_model = TRUE
) {
  # n_sim = 3; S_tr = diag(n.tr); n_cores = n.cores / 2; run_model = TRUE
  
  rng <- RNGseq(n_sim * n.subj, 0)  ## set seed for each iteration:
  ## https://cran.r-project.org/web/packages/doRNG/vignettes/doRNG.pdf
  
  time.start <- Sys.time()
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  B_rsa <- 
    foreach(
      sim.i = seq_len(n_sim),
      .export = ls(globalenv()),
      .inorder = FALSE
    ) %:% 
    
    foreach(
      subj.i = seq_along(subjs),
      r = rng[(sim.i - 1) * n.subj + 1:n.subj],
      .inorder = FALSE
    ) %dopar% {
      # subj.i <- 1
      
      rngtools::setRNG(r)
      
      ## generate noise
      
      E1 <- mvnfast::rmvn(n.vox, mu = numeric(n.tr), sigma = S_tr)
      E2 <- mvnfast::rmvn(n.vox, mu = numeric(n.tr), sigma = S_tr)
      E <- t(cbind(E1, E2))
      
      ## time-series GLM, similarity matrix
      
      B <- coef(.lm.fit(X_glm[, , subj.i], E))
      
      ## RSA
      
      R <- cor(t(B[is.bias.item, ]))  ## estimate similarity matrix
      dimnames(R) <- list(conds[is.bias.item], conds[is.bias.item])
      
      y <- scale(rank(R[lt]))
      
      if (run_model) {
        
        fit <- .lm.fit(U, y)  ## regress run component from rsv
        b1 <- coef(fit)[2]  ## slope
        y <- y - U[, 2] * b1  ## prewhitened response vector
        
      }
      
      
      b_rsa <- c(coef(.lm.fit(X_rsa, y)), subj.i)
      names(b_rsa) <- c("target", "distractor", "incongruency", "congruency", "subj")
      
      ## save
      
      b_rsa
      
    }
  
  
  
  stopCluster(cl)
  # (time.stop <- Sys.time() - time.start)
  
  
  d <- dplyr::bind_rows(B_rsa)
  d <- reshape2::melt(d, value.name = "beta", variable.name = "param", id = "subj")
  
  d
  
  
  
}



sim_null_simil <- function(
  n_sim,
  n_cores = n.cores / 2,
  run_model = TRUE
) {
  # n_sim = 3; S_tr = diag(n.tr); n_cores = n.cores / 2; run_model = TRUE
  
  rng <- RNGseq(n_sim * n.subj, 0)  ## set seed for each iteration:
  ## https://cran.r-project.org/web/packages/doRNG/vignettes/doRNG.pdf
  
  # time.start <- Sys.time()
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  R_rsa <- 
    
    foreach(
      sim.i = seq_len(n_sim),
      .export = ls(globalenv()),
      .inorder = FALSE
    ) %:% 
    
    foreach(
      subj.i = seq_along(subjs),
      r = rng[(sim.i - 1) * n.subj + 1:n.subj],
      .inorder = FALSE
    ) %dopar% {
      # subj.i <- 1
      
      rngtools::setRNG(r)
      
      ## generate noise
      
      # E1 <- mvnfast::rmvn(n.vox, mu = numeric(n.tr), sigma = S_tr)
      # E2 <- mvnfast::rmvn(n.vox, mu = numeric(n.tr), sigma = S_tr)
      # E <- t(cbind(E1, E2))
      E <- matrix(rnorm(n.vox*n.tr*2), nrow = n.vox, ncol = n.tr*2)
      
      ## time-series GLM, similarity matrix
      
      B <- coef(.lm.fit(X_glm[, , subj.i], E))
      
      ## RSA
      
      R <- cor(t(B[is.bias.item, ]))  ## estimate similarity matrix
      dimnames(R) <- list(conds[is.bias.item], conds[is.bias.item])
      
      R
      
    }
  
  stopCluster(cl)
  # (time.stop <- Sys.time() - time.start)
  
  R_rsa
  
  
}



symmat4ggplot <- function(R, var.names = c("v1", "v2"), val.name = "value") {
  
  ## make factors for row and column labels
  dn <- dimnames(R)
  if (is.null(dn)) {
    dn <- setNames(list(paste0("cell_", 1:nrow(R)), paste0("cell_", 1:ncol(R))), var.names)
  } else {
    names(dn) <- var.names  
  }
  
  labels <- expand.grid(dn, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = TRUE)
  labels[[2]] <- factor(labels[[2]], levels = rev(levels(labels[[2]])))
  
  r <- c(R)
  
  cbind(labels, setNames(as.data.frame(c(R)), val.name))
  
}


matplot <- function(x) {
  
  ggplot(symmat4ggplot(x), aes(v1, v2, fill = value)) +
    geom_raster() +
    scale_fill_viridis_c(option = "inferno") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
      axis.title = element_blank(), 
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid = element_blank()
    )
  
}


## strings, numerics, ...


subjs <- subjs.analysis
n.subj <- length(subjs)
n.cores <- detectCores()
X_glm <- readRDS(here("glms", "xmats_pro_bias_acc-only.rds"))[, , subjs]
n.tr <- dim(X_glm)[1]/2
n.vox <- 100
conds <- dimnames(X_glm)$condition
n.cond <- conds
is.bias.item <- conds %in% bias.items

any(conds[is.bias.item] %in% bias.items)


## modeling objects

## build X_rsa, for RSA models

lt <- lower.tri(diag(16))
X_rsa <- fread(here("out", "rsa", "mods", "rsv_bias_lower-triangles.csv"))
X_rsa <- scale(as.matrix(X_rsa[, c("target", "distractor", "incongruency", "congruency")]))

## build U, run-bias model (for prewhitening of similarity matrices)

run.rsm <- as.matrix(fread(here("out", "rsa", "mods", "rsm_pro_bias_run.csv")), rownames = 1)
run.rsv <- mat2vec(run.rsm, value.name = "run.model")  ## unwrap run model to lower-tri vector
U <- cbind(1, run.rsv$run.model)  ## model: intercept, run regressor

## simulation params

ar.coefs <- c(seq(0, 0.9, 0.1), 0.99, 0.999)
sigmas <- lapply(ar.coefs, ar1_cor, n = n.tr)


n.sims <- 1000


## simulate across range of AR(1)s ----


if (!file.exists(here("out", "runwise", "simulation.rds"))) {
  
  l <- vector("list", length(ar.coefs))
  t0 <- Sys.time()
  for (ar.coef.i in seq_along(ar.coefs)) {
    # ar.coef.i = 1
    
    d <- sim_null(n.sims, sigmas[[ar.coef.i]])
  
    d$ar.coef <- ar.coefs[ar.coef.i]
    
    l[[ar.coef.i]] <- d
    
    print(paste0("done: ", ar.coef.i, "/", length(ar.coefs)))
    
  }
  (runtime <- Sys.time() - t0)
  
  res <- bind_rows(l)
  
  res %<>% 
    group_by(subj, ar.coef) %>%
    mutate(sim = 1:n())
  
  saveRDS(res, here("out", "runwise", "simulation.rds"))

} else {
  
  res <- readRDS(here("out", "runwise", "simulation.rds"))
  
}


## look

group.stats <- res %>%
  group_by(ar.coef, param, sim) %>%
  summarize(
    b = mean(beta),
    d = b / sd(beta),
    p = wilcox.test(beta, alternative = "greater")$p.value
  )

group.stats %>%
  ggplot(aes(as.factor(ar.coef), b, fill = param)) +
  geom_boxplot()

group.stats %>%
  ggplot(aes(as.factor(ar.coef), b, color = param)) +
  stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.5))


p_res <- group.stats %>%
  
  filter(param != "congruency") %>%
  group_by(ar.coef, param) %>%
  summarize(fpr = mean(p < 0.05)) %>%
  
  ggplot(aes(ar.coef, fpr, color = param)) +
  
  scale_color_manual(values = colors.model) +
  scale_y_continuous(limits = c(0, 0.5)) +
  scale_x_continuous(limits = c(0, 1), breaks = c(seq(0, 0.8, 0.2), 0.999)) +
  coord_capped_cart(left = "both", bottom = "both") +

  geom_hline(yintercept = 0.05, alpha = 0.25) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  
  annotate("text", x = 0, y = 0.25, label = "target", color = colors.model["target"], fontface = "bold", hjust = 0, size = 5) +
  annotate("text", x = 0, y = 0.2, label = "distractor", color = colors.model["distractor"], fontface = "bold", hjust = 0, size = 5) +
  annotate("text", x = 0, y = 0.15, label = "incongruency", color = colors.model["incongruency"], fontface = "bold", hjust = 0, size = 5) +
  
  
  theme(
    panel.grid      = element_blank(),
    panel.border    = element_blank(),
    axis.line       = element_line(size = axis.line.size),
    axis.text       = element_text(size = axis.text.size),
    axis.ticks      = element_line(size = axis.line.size),
    axis.title      = element_text(size = axis.title.size*1.5),
    strip.background = element_blank(),
    strip.text = element_text(size = axis.title.size),
    legend.position = "none"
  ) +
  
  labs(y = "false positive rate", x = bquote("AR(1) coefficient ("*phi*")"))


# ggsave(here("out", "runwise", "simulation_fpr.jpg"), dev = "jpg", width = 7, height = 5)




## simulate rsm ----

if (!file.exists(here("out", "runwise", "simulation_rsm.rds"))) {
  
  simil <- sim_null_simil(10000, diag(n.tr), n_cores = 12)
  saveRDS(simil, here("out", "runwise", "simulation_rsm.rds"))
  
} else {
  
  simil <- readRDS(here("out", "runwise", "simulation_rsm.rds"))
  
}

## look

R <- simil[[1]][[1]]
R[] <- 0
diag(R) <- 1
for (ii in seq_along(simil)) {
  simil.ii <- simil[[ii]]
  for (jj in seq_along(simil.ii)) {
    R <- R + atanh(simil.ii[[jj]]) / (length(simil.ii)*length(simil))
  }
}
R <- tanh(R)

# axis.text.x = element_text(color = labels$colors, face = "bold", size = rel(4), angle = 90, hjust = 1, vjust = 0.5),
# axis.text.y = element_text(color = rev(labels$colors), face = "bold", size = rel(4)),


bias.colors.plot <- ifelse(bias.colors == "white", "grey50", bias.colors)
R <- R[sort(bias.items), sort(bias.items)]
p_rsm <- matplot(R[sort(bias.items), sort(bias.items)]) + 
  theme(
    legend.position = "right", 
    axis.text.x = element_text(color = rep(bias.colors.plot, each = 4), face = "bold", size = axis.text.size*0.7),
    axis.text.y = element_text(color = rev(rep(bias.colors.plot, each = 4)), face = "bold", size = axis.text.size*0.7)
    ) + 
  labs(fill = NULL) +
  scale_x_discrete(labels = rep(bias.words, 4)) +
  scale_y_discrete(labels = rev(rep(bias.words, 4)))

Rna <- R
diag(Rna) <- NA
p_rsm_scaled <- matplot(Rna) + 
  theme(
    legend.position = "right", 
    axis.text.x = element_text(color = rep(bias.colors.plot, each = 4), face = "bold", size = axis.text.size*0.7),
    axis.text.y = element_text(color = rev(rep(bias.colors.plot, each = 4)), face = "bold", size = axis.text.size*0.7)
  ) + 
  labs(fill = NULL) +
  scale_x_discrete(labels = rep(bias.words, 4)) +
  scale_y_discrete(labels = rev(rep(bias.words, 4)))


change_legend_size <- function(myPlot, pointSize = 2, textSize = 5, spaceLegend = 0.75) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
p_rsm <- change_legend_size(p_rsm)
p_rsm_scaled <- change_legend_size(p_rsm_scaled)

p <- plot_grid(
  p_res, plot_grid(p_rsm, p_rsm_scaled, nrow = 2), nrow = 1, rel_widths = c(1, 2/3),
  labels = c("A", "B")
)

ggsave(here("out", "runwise", "simulation_rsm.pdf"), p, dev = "pdf", width = 17.6, height = 17.6/5*3, unit = "cm")



phi_min <- matplot(sigmas[[1]]) + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_text()) + 
  labs(x = "TR", y = "TR", title = bquote(Sigma[phi==0]))

phi_max <- matplot(sigmas[[length(sigmas)]]) + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_text()) + 
  labs(x = "TR", y = "TR", title = bquote(Sigma[phi==0.999]))

ggsave(here("out", "runwise", "phi_min.pdf"), phi_min, dev = "pdf", width = 7, height = 7.5, unit = "cm")
ggsave(here("out", "runwise", "phi_max.pdf"), phi_max, dev = "pdf", width = 7, height = 7.5, unit = "cm")


# S_tr <- sigmas[[1]]
# S_tr1 <- diag(n.tr*2)
# 
# mb <- microbenchmark::microbenchmark(
#   
#   "fast" = {
#     
#     E1 <- mvnfast::rmvn(n.vox, mu = numeric(n.tr), sigma = S_tr)
#     E2 <- mvnfast::rmvn(n.vox, mu = numeric(n.tr), sigma = S_tr)
#     E <- t(cbind(E1, E2))
#     
#   },
#   
#   
#   "rnorm" = {
#     
#     E <- matrix(rnorm(n.vox*n.tr*2), nrow = n.vox, ncol = n.tr*2)
#     
#   },
#   
#   
#   "fast2" = {
#     
#     E <- mvnfast::rmvn(n.vox, mu = numeric(n.tr*2), sigma = S_tr1)
#     
#   }
# 
#   
# )
# 
# autoplot(mb)





## ----
# res <- sim_null(n.sims, sigmas[[1]])
# res %<>%
#   group_by(subj) %>%
#   mutate(sim = 1:n())
# group.stats <- res %>%
#   group_by(param, sim) %>%
#   summarize(
#     b = mean(beta),
#     d = b / sd(beta),
#     p = wilcox.test(beta, alternative = "greater")$p.value
#   )
# group.stats %>%
#   group_by(param) %>%
#   summarize(fpr = mean(p < 0.05))
# res.nreg <- sim_null(n.sims, sigmas[[1]], run_model = FALSE)
# 
# res.reg$run_model <- TRUE
# res.nreg$run_model <- FALSE
# res <- full_join(res.reg, res.nreg)
# 
# res %<>%
#   group_by(subj, run_model) %>%
#   mutate(sim = 1:n())
# 
# 
# group.stats <- res %>%
#   group_by(run_model, param, sim) %>%
#   summarize(
#     b = mean(beta),
#     d = b / sd(beta),
#     p = wilcox.test(beta, alternative = "greater")$p.value
#   )
# 
# group.stats %>%
#   ggplot(aes(run_model, b, fill = param)) +
#   geom_boxplot()
# 
# 
# group.stats %>%
#   group_by(run_model, param) %>%
#   summarize(fpr = mean(p < 0.05))
#   ggplot(aes(ar.coef, fpr, color = param)) +
#   geom_line()
# 
# # A tibble: 4 x 2
# param          fpr
# <fct>        <dbl>
#   1 target       0.106
# 2 distractor   0.067
# 3 incongruency 0.044
# 4 congruency   0.094
# >
## scratch ----


## random walk instead of AR?

# group.stats <- res %>%
#   group_by(run_model, param, sim) %>%
#   summarize(
#     b = mean(beta),
#     d = b / sd(beta),
#     p = wilcox.test(beta, alternative = "greater")$p.value
#   )




# XX <- t(X_glm.i) %*% X_glm.i
# XXinv <- solve(XX)
# XXinv <- XXinv[bias.items, bias.items]


# B <- coef(.lm.fit(X_rsa, XX[bias.items, bias.items][lt]))
# names(B) <- c("target", "distractor", "incongruency")
# stats.subjs[[subj.i]] <- B


# r <- 0.5
# Y <- replicate(n.subj*n.vox*2, arima.sim(list(order = c(length(r), 0, 0), ar = r), n = n.tr))  ## 2 runs per subj*vox
# Y <- array(Y, dim = c(n.tr, n.vox, n.subj, 2))

# library(marima)
# marima.sim(n.vox, )


