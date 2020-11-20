

# subjmeans.super <- stats.subjs.super %>%
#   
#   mutate(roi.hemi = paste0(roi, "_", hemi)) %>%
#   filter(roi.hemi %in% c("dlpfc_L", "dlpfc_R", "dmfc_L", "lppc_L")) %>%
#   
#   tidyr::pivot_wider(names_from = "param", values_from = "beta") %>%
#   
#   group_by(roi, subj) %>%
#   summarize_if(is.numeric, mean, .groups = "drop_last")
# 
# means.super <- subjmeans.super %>% 
#   
#   split(.$roi) %>% 
#   lapply(morey08, c("target", "distractor", "incongruency"), "param") %>%
#   
#   bind_rows(.id = "roi")
# 
# means.super %>%
#   
#   mutate(param = factor(as.factor(param), levels = c("incongruency", "target", "distractor"))) %>%
#   
#   ggplot(aes(roi, mu, color = param)) +
#   geom_point(size = geom.point.size, position = position_dodge(width = 1/2)) +
#   geom_errorbar(aes(ymin = mu - ci, ymax = mu + ci), position = position_dodge(width = 1/2), width = 0, size = geom.line.size)
# geom_line(size = geom.line.size, position = position_dodge(width = 1/2)) +
# 
# 
# 
# x = subjmeans.super %>% filter(roi == "dmfc")
# cols = c("target", "distractor", "incongruency")
# morey08 <- function(x, cols, varname, p = 0.95) {
#   
#   x <- x[, cols]
#   
#   m <- ncol(x)
#   n <- nrow(x)
#   df <- n - m
#   p <- p + (1 - p) / 2  ## 1-tailed to 2-tailed
#   q <- qt(p, df)  ## quantile of t distribution
#   # q <- qt(p, Inf)  ## quantile of t distribution
#   
#   normalized <- x - rowMeans(x) + mean(colMeans(x))
#   # sweep(x, 1, rowMeans(x), "-") + mean(unlist(x))
#   correction <- m / (m - 1)
#   
#   s2 <- apply(normalized, 2, var) * correction
#   ci <- sqrt(s2 / n) * q
#   
#   d <- data.frame(
#     mu = colMeans(x),
#     ci = ci, 
#     variable = cols
#   )
#   
#   if (missing(varname)) return(d)
#   
#   names(d)[3] <- varname
#   d
#   
# }
# 
# y <- subjmeans.super %>% filter(roi == "dmfc", ) %>% melt
# a <- summary(papaja::wsci(y, "subj", "variable", "value"))
# d <- subjmeans.super %>% filter(roi == "dmfc") %>% morey08(c("target", "incongruency"), "param")
# a$lower_limit
# 
# 
# subjmeans.super %>% group_by(roi) %>%
#   summarize(
#     p = t.test(incongruency - target)$p.value
#   )
# 
# 
