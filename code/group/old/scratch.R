

# rsarray.mmp.rois <- rsarray.mmp[, , , c(rois.mmp.dlpfc, rois.mmp.mfc, rois.mmp.lppc)]
# 
# boot.noise.ceiling <- function(x, rois, ii) {
#   # x = rsarray.mmp.rois; rois = rois.mmp.dlpfc[1]; ii = sample.int(49, replace = TRUE)
#   
#   xi <- x[, , ii, ]  ## resample subjs
#   dimnames(xi)
#   is.lt <- lower.tri(xi[, , 1, 1])  ## for extracting lower tri
#   nsubj <- length(ii)
#   
#   ub <- rep(NA, length(rois))  ## for storing mean upper-bound per ROI
#   lb <- ub  ## for storing mean lower-bound per ROI
#   lb.i <- rep(NA, length(nsubj))  ## for storing lower-bound per subj (for a given ROI)
#   for (roi.i in seq_along(rois)) {
#     # roi.i = 1
#     
#     v <- apply(xi[, , , roi.i], 3, "[", is.lt)  ## lower-triangle RSM
#     
#     centroid <- rowMeans(v)
#     ub[roi.i] <- mean(atanh(cor(centroid, v)))  ## mean upper-bound z for roi.i
#     
#     lb.ii <- lb.i  ## copy
#     for (subj.i in seq_len(nsubj)) lb.ii[subj.i] <- cor(rowMeans(v[, -subj.i]), v[, subj.i])
#     lb[roi.i] <- mean(atanh(lb.ii))  ## mean lower-bound z for roi.i
#     
#   }
#   
#   c(lb = mean(lb), ub = mean(ub))  ## average across ROIs
#   
# }
# 
# 
# boot.noise.ceiling.contrast <- function(x, rois1, rois2, ii) {
#   # x = rsarray.mmp.rois; rois1 = rois.mmp.dlpfc; rois2 = rois.mmp.lppc; ii = 1:49
#   
#   xi <- x[, , ii, ]  ## resample subjs
#   is.lt <- lower.tri(xi[, , 1, 1])  ## for extracting lower tri
#   nsubj <- length(ii)
#   
#   ## get mean for rois1
#   
#   xi1 <- xi[, , , rois1]
#   ub1 <- rep(NA, length(rois1))  ## for storing mean upper-bound per ROI
#   lb1 <- ub1  ## for storing mean lower-bound per ROI
#   lb.i <- rep(NA, length(nsubj))  ## for storing lower-bound per subj (for a given ROI)
#   
#   for (roi.i in seq_along(rois1)) {
#     # roi.i = 1
#     
#     v1 <- apply(xi1[, , , roi.i], 3, "[", is.lt)  ## lower-triangle RSM
#     
#     centroid1 <- rowMeans(v1)
#     ub1[roi.i] <- mean(atanh(cor(centroid1, v1)))  ## mean upper-bound z for roi.i
#     
#     lb.ii1 <- lb.i  ## copy
#     for (subj.i in seq_len(nsubj)) lb.ii1[subj.i] <- cor(rowMeans(v1[, -subj.i]), v1[, subj.i])
#     lb1[roi.i] <- mean(atanh(lb.ii1))  ## mean lower-bound z for roi.i
#     
#   }
#   
#   ## get mean for rois2
#   
#   xi2 <- xi[, , , rois2]
#   ub2 <- rep(NA, length(rois2))  ## for storing mean upper-bound per ROI
#   lb2 <- ub2  ## for storing mean lower-bound per ROI
#   lb.i <- rep(NA, length(nsubj))  ## for storing lower-bound per subj (for a given ROI)
#   
#   for (roi.i in seq_along(rois2)) {
#     # roi.i = 1
#     
#     v2 <- apply(xi2[, , , roi.i], 3, "[", is.lt)  ## lower-triangle RSM
#     
#     centroid2 <- rowMeans(v2)
#     ub2[roi.i] <- mean(atanh(cor(centroid2, v2)))  ## mean upper-bound z for roi.i
#     
#     lb.ii2 <- lb.i  ## copy
#     for (subj.i in seq_len(nsubj)) lb.ii2[subj.i] <- cor(rowMeans(v2[, -subj.i]), v2[, subj.i])
#     lb2[roi.i] <- mean(atanh(lb.ii2))  ## mean lower-bound z for roi.i
#     
#   }
#   
#   
#   ## average across ROIs (within set), and contrast (across sets):
#   
#   c(lb = mean(lb1) - mean(lb2), ub = mean(ub1) - mean(ub2))
#   
# }

















res.boot.ceiling.mmp.rois <- lapply(
  rois.mmp[1],
  function(., d, statistic, R, ...) boot::boot(d, statistic, roi = ., R = R, ...),
  d = rsarray.mmp.rois,
  statistic = boot.noise.ceiling,
  R = 2E4
)

boot.ci(res.boot.ceiling.mmp.rois[[1]], index = 2)

lapply(
  
  res.boot.ceiling.mmp.rois,
  
  function(.) {
    
    ci <- boot::boot.ci(., type = type)[[type]][4:5]
    data.frame(y = out$t0, ymin = ci[1], ymax = ci[2])
    
  }
  
)










boot.noise.ceiling <- function(x, rois, ii) {
  # x = rsarray.mmp; roi = c("V1_L", "V2_L")
  # roi.i = 1
  xi <- x[, , ii, roi]
  is.lt <- lower.tri(xi[, , , 1])
  
  for (roi.i in seq_along(rois)) {
    v <- apply(xi, 3, "[", is.lt)
    
  }
  
  nsubj <- ncol(v)
  subjs <- colnames(v)
  
  centroid <- rowMeans(v)
  ub <- mean(cor(centroid, v))
  
  lb <- rep(NA, nsubj)
  for (subj.i in seq_len(nsubj)) lb[subj.i] <- cor(rowMeans(v[, -subj.i]), v[, subj.i])
  lb <- mean(lb)
  
  cbind(ub = ub, lb = lb)
  
}

## estimate ----

ceilings.mmp <- lapply(dimnames(rsarray.mmp)$roi, noise.ceiling, x = rsarray.mmp)
names(ceilings.mmp) <- dimnames(rsarray.mmp)$roi

ceilings.mmp %<>% bind_rows(.id = "roi") %>% full_join(atlas.key$mmp, by = "roi")


ceilings.mmp.rois <- ceilings.mmp %>% 
  
  filter(roi %in% c(rois.mmp.mfc, rois.mmp.dlpfc, rois.mmp.lppc)) %>%
  
  mutate(
    region = ifelse(
      roi %in% rois.mmp.dlpfc, "DLPFC", 
      ifelse(
        roi %in% rois.mmp.mfc, "MFC",
        ifelse(roi %in% rois.mmp.lppc, "LPPC", NA)
      )
    )
  )


# fit.mmp.lb <- lmer(atanh(lb) ~ -1 + roi + (1 | subj), ceilings.mmp.rois)
# summary(fit.mmp.lb)
# regnames.mmp.ceilings <- gsub("roi", "", rownames(coef(summary(fit.mmp.lb))))
# W.mmp.ceilings <- rbind(
#   dlpfc = regnames.mmp.ceilings %in% rois.mmp.dlpfc,
#   mfc   = regnames.mmp.ceilings %in% rois.mmp.mfc,
#   lppc  = regnames.mmp.ceilings %in% rois.mmp.lppc
# )
# contr.ceilings <- rbind(
#   "dlpfc-mfc"  = c(1, -1, 0),
#   "lppc-mfc"   = c(0, -1, 1),
#   "dlpfc-lppc" = c(1, 0, -1)
#   )
# W.mmp.ceilings <- rbind(W.mmp.ceilings, contr.ceilings %*% W.mmp.ceilings)
# 
# 
# summary(glht(fit.mmp.lb, W.mmp.ceilings, test = adjusted("none")))
# 
# 
# fit.mmp.ub <- refit(fit.mmp.lb, ceilings.mmp.rois$ub)
# summary(fit.mmp.ub)








## estimate mean ceiliings ----

res.boot.ceiling.mmp.dlpfc <- boot(
  rsarray.mmp.rois[, , , rois.mmp.dlpfc], 
  statistic = boot.noise.ceiling, 
  rois = rois.mmp.dlpfc,
  R = 1E4
)

res.boot.ceiling.mmp.lppc <- boot(
  rsarray.mmp.rois[, , , rois.mmp.lppc], 
  statistic = boot.noise.ceiling, 
  rois = rois.mmp.lppc,
  R = 1E4
)

res.boot.ceiling.mmp.mfc <- boot(
  rsarray.mmp.rois[, , , rois.mmp.mfc], 
  statistic = boot.noise.ceiling, 
  rois = rois.mmp.mfc,
  R = 1E4
)


## estimate differences in ceiliings ----

res.boot.ceiling.mmp.dlpfc.mfc <- boot(
  rsarray.mmp.rois[, , , c(rois.mmp.dlpfc, rois.mmp.mfc)], 
  statistic = boot.noise.ceiling.contrast, 
  rois1 = rois.mmp.dlpfc,
  rois2 = rois.mmp.mfc,
  R = 1E4
)
res.boot.ceiling.mmp.lppc.mfc <- boot(
  rsarray.mmp.rois[, , , c(rois.mmp.lppc, rois.mmp.mfc)], 
  statistic = boot.noise.ceiling.contrast, 
  rois1 = rois.mmp.lppc,
  rois2 = rois.mmp.mfc,
  R = 1E4
)
res.boot.ceiling.mmp.dlpfc.lppc <- boot(
  rsarray.mmp.rois[, , , c(rois.mmp.dlpfc, rois.mmp.lppc)], 
  statistic = boot.noise.ceiling.contrast, 
  rois1 = rois.mmp.dlpfc,
  rois2 = rois.mmp.lppc,
  R = 1E4
)


## wrangle results ----

means.ceilings.mmp.rois <- lapply(
  list(
    dlpfc = res.boot.ceiling.mmp.dlpfc,
    lppc = res.boot.ceiling.mmp.lppc,
    mfc = res.boot.ceiling.mmp.mfc
  ),
  function(., ...) {
    m <- rbind(c(
      .$t0,
      boot.ci(., type = "perc", index = 1)[["percent"]][4:5],
      boot.ci(., type = "perc", index = 2)[["percent"]][4:5]
    ))
    setNames(as.data.frame(m), c("lb", "ub", "lb.ymin", "lb.ymax", "ub.ymin", "ub.ymax"))
  }
) %>%
  bind_rows(.id = "region")

p.ceilings.mmp.rois <- lapply(
  list(
    dlpfc.mfc  = res.boot.ceiling.mmp.dlpfc.mfc,
    lppc.mfc   = res.boot.ceiling.mmp.lppc.mfc,
    dlpfc.lppc = res.boot.ceiling.mmp.dlpfc.lppc
  ),
  ci2p
) %>%
  bind_rows(.id = "region")



means.ceilings.mmp.rois %>%  
  
  ggplot(aes(region, lb, color = region)) +
  geom_point(size = 6) +
  geom_errorbar(aes(ymin = lb.ymin, ymax = lb.ymax), width = 0, size = 2) +
  geom_line(size = 2)

scale_color_manual(values = colors.region) +
  scale_y_continuous(breaks = c(0.02, 0.1), limits = c(0.02, 0.1)) +
  
  annotate(
    geom = "text", x = -Inf, y = 0.095, label = "DLPFC", color = colors.region["DLPFC"],
    hjust = 0, vjust = 1, size = rel(14), fontface = "bold"
  ) +
  annotate(
    geom = "text", x = -Inf, y = 0.09, label = "LPPC", color = colors.region["LPPC"],
    hjust = 0, vjust = 1, size = rel(14), fontface = "bold"
  ) +
  annotate(
    geom = "text", x = -Inf, y = 0.085, label = "MFC", color = colors.region["MFC"],
    hjust = 0, vjust = 1, size = rel(14), fontface = "bold"
  ) +
  
  labs(y = bquote("mean model fit ("*beta*")"), x = "model") +
  
  theme(
    legend.position = "none", 
    panel.grid      = element_blank(), 
    panel.border    = element_blank(),
    axis.line.y     = element_line(size = rel(2)),
    axis.text       = element_text(size = rel(3)),
    axis.text.x     = element_text(color = colors.model),
    axis.ticks.y    = element_line(size = rel(2)),
    axis.ticks.x    = element_blank(),
    axis.title    = element_text(size = rel(4))
  )  ## will trip warning due to vectorized input to axis.text.x (for colors)



boot.ci(res.boot.ceiling.mmp.dlpfc, type = "perc", index = 1)[["percent"]]
boot.ci(res.boot.ceiling.mmp.dlpfc, type = "bca", index = 1)$t

lapply(
  rois.mmp.dlpfc,
  function(x) colMeans(atanh(noise.ceiling(rsarray.mmp.rois, x)[, 1:2]))
) %>%
  do.call(rbind, .) %>%
  colMeans

noise.ceiling









# fit.mmp <- lmer(
#   beta ~ region * param + (region * param | subj),
#   d.mmp %>% filter(region %in% c("DLPFC", "MFC", "LPPC"), param %in% c("target", "incongruency"))
#   )
# summary(fit.mmp)
# 
# ## ORDER:
# ## dlpfcincon, lppc-dlpfc|incon, mfc-dlpfc|incon, target-incon|dlpfc, (lppc-dlpfc)(target-incon), (mfc-dlpfc)(target-incon)
# dlpfc.incon <- c(1, 0, 0, 0, 0, 0)
# mfc.incon   <- c(1, 0, 1, 0, 0, 0)
# lppc.incon  <- c(1, 1, 0, 0, 0, 0)
# dlpfc.targt <- c(1, 0, 0, 1, 0, 0)
# mfc.targt   <- c(1, 0, 1, 1, 0, 1)
# lppc.targt  <- c(1, 1, 0, 1, 1, 0)
# 
# W <- rbind(
# 
#   ## means of each condition:
#   "dlpfc_incon" = dlpfc.incon,
#   "dlpfc_targt" = dlpfc.targt,
#   "mfc_incon"   = mfc.incon,
#   "mfc_targt"   = mfc.targt,
#   "lppc_incon"  = lppc.incon,
#   "lppc_targt"  = lppc.targt,
# 
#   ## within-parcel contrasts:
#   "incon-targt|mfc"   = mfc.incon - mfc.targt,
#   "incon-targt|dlpfc" = dlpfc.incon - dlpfc.targt,
#   "incon-targt|lppc"  = lppc.incon - lppc.targt,
# 
#   ## between-parcel contrasts:
#   "mfc-dlpfc|targt" = mfc.targt - dlpfc.targt,
#   "mfc-lppc|targt"  = mfc.targt - lppc.targt,
#   "mfc-dlpfc|incon" = mfc.incon - dlpfc.incon,
#   "mfc-lppc|incon"  = mfc.incon - lppc.incon,
# 
#   ## interactions:
#   "(incon-targt)(mfc-dlpfc)" = (mfc.targt - dlpfc.targt) - (mfc.incon - dlpfc.incon),
#   "(incon-targt)(mfc-lppc)"  = (mfc.targt - lppc.targt) - (mfc.incon - lppc.incon)
# 
# )
# 
# summary(glht(fit.mmp, W), test = adjusted("none"))
















# 
# ## within-region contrasts:
# "incon-targt|mfc"   = (is.incon - is.targt) * is.mfc,
# "incon-targt|dlpfc" = (is.incon - is.targt) * is.dlpfc,
# "incon-targt|lppc"  = (is.incon - is.targt) * is.lppc,
# 
# ## between-parcel contrasts:
# "mfc-dlpfc|targt" = (is.mfc - is.dlpfc) * is.targt,
# "mfc-lppc|targt"  = (is.mfc - is.lppc) * is.targt,
# "mfc-dlpfc|incon" = (is.mfc - is.dlpfc) * is.incon,
# "mfc-lppc|incon"  = (is.mfc - is.lppc) * is.incon,
# 
# ## interactions:
# "(incon-targt)(mfc-dlpfc)" = (is.incon - is.targt) * (is.mfc - is.dlpfc),
# "(incon-targt)(mfc-lppc)"  = (is.incon - is.targt) * (is.mfc - is.lppc)
# 
