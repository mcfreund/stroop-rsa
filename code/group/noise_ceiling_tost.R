#+ noise-ceiling-tost_setup, include = FALSE

if (interactive()) source(here::here("code", "group", "noise_ceiling.R"))

#' ## select the "smallest effect-size of interest"
# * here, defined as the smallest effect size for which we have 80% power to detect.

#+ noise-ceiling-tost_sesoi, fig.height = 7, fig.width = 9

n.resamps <- 1E4
n.subj <- 49

fname.seoi <- here("out", "group", "noise_ceiling_sesoi.txt")
fname.curve <- here("out", "group", "noise_ceiling_effect-size-curve.txt")

if (!file.exists(fname.seoi) || !file.exists(fname.curve)) {
  
  es <- seq(0.3, 0.4, 0.001)  ## define grid of effect sizes; range obtained from rough first-pass
  n.es <- length(es)
  n.tests <- 1E4
  thresh <- 0.05  ## alpha
  pwr <- numeric(n.es)
  n.cores <- detectCores()
  
  set.seed(0)
  for (n.es.i in seq_along(es)) {
    # n.es.i = 1
    
    X <- matrix(rnorm(n.subj * n.tests, es[n.es.i]), nrow = n.tests, ncol = n.subj)
    
    cl <- makeCluster(n.cores - 1)
    registerDoParallel(cl)
    freqs <- foreach(ii = seq_len(n.resamps), .combine = "+", .inorder = FALSE) %dopar% {
      rowMeans(X[, sample.int(n.subj, replace = TRUE)]) < 0
    }
    stopCluster(cl)
    
    pwr[n.es.i] <- sum(freqs / n.resamps <= thresh) / n.tests
    
  }
  
  d <- data.frame(effect.size = es, hit.rate = pwr)
  fit <- fANCOVA::loess.as(d$effect.size, d$hit.rate)  ## smooth over resampling error with loess
  d$fitted <- fit$fitted
  (sesoi <- approx(d$fitted, d$effect.size, xout = 0.8)$y)  ## interpolate minimum effect size of interest
  
  sink(fname.seoi)
  cat(sesoi)
  sink()
  
  fwrite(d, here("out", "group", "noise_ceiling_effect-size-curve.txt"))
  

} else {
  
  sesoi <- as.numeric(readLines(fname.seoi, 1, warn = FALSE))
  d <- fread(fname.curve)
  
}

## plot

plot(d$effect.size, d$hit.rate, main = "effect size (cohen's d) by hit rate (power) of percentile bootstrap")
lines(d$effect.size, d$fitted, col = "firebrick")
abline(h = 0.8)
abline(v = sesoi)


## get standardized contrasts:

es.ceilings.contr.rois <- ceilings.rois %>%
  
  dplyr::select(subj, region, lb) %>%
  mutate(lb = atanh(lb)) %>%
  
  tidyr::pivot_wider(., names_from = "region", values_from = "lb") %>%
  ungroup %>%
  mutate_if(is.numeric, function(.) . / sd(.)) %>%
  
  transmute(
    DLPFC.DMFC  = DLPFC - DMFC,
    LPPC.DMFC   = LPPC - DMFC,
    DLPFC.LPPC = LPPC - DLPFC
  )

## perform two one-sided tests:

less <- matrix(ncol = 3, nrow = n.resamps, dimnames = list(NULL, names(es.ceilings.contr.rois)))
greater <- less
set.seed(0)
for (ii in seq_len(n.resamps)) {
  means.ii <- colMeans(es.ceilings.contr.rois[sample.int(n.subj, replace = TRUE), ])
  less[ii, ] <- sesoi - means.ii   ## less than sesoi?
  greater[ii, ] <- means.ii + sesoi  ## greater than negative sesoi?
}

pvalues <- rbind(colMeans(less < 0), colMeans(greater < 0))
(tost.pvalues <- apply(pvalues, 2, max))
fwrite(as.data.frame(tost.pvalues), here("out", "group", "noise_ceiling_tost_p_values.txt"))

