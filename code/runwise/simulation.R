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



## simulate function


sim_null <- function(
  n_sim,
  S_tr, 
  n_cores = n.cores / 2,
  run_model = TRUE
  ) {
  # n_sim = 3; S_tr = diag(n.tr); run_model = TRUE

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




## simulate ----

ar.coefs <- c(0, 0.5, 0.9, 0.99, 1)
sigmas <- lapply(ar.coefs, ar1_cor, n = n.tr)
n.sims <- 1000

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



# res %>%
#   group_by(param) %>%
#   summarize(
#     b = mean(beta),
#     p = wilcox.test(beta, alternative = "greater")$p.value
#   )







## scratch ----


## random walk instead of AR?






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


