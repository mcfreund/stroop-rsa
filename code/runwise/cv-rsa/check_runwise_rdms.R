
source(here::here("code", "packages.R"))
source(here("code", "strings.R"))


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
      axis.title = element_blank(), 
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
      legend.position = "right",
      panel.border = element_blank(), panel.grid = element_blank()
    )
  
}



rsarray <- readRDS(
  here::here("out", "rsa", "obsv", "rsarray_pro_bias_acc-only_fmriprep_runwise_cv-euclidean-stand_mmp.rds")
)


dimnames(rsarray)
rsarray <- rsarray[sort(bias.items), sort(bias.items), , ]
D <- apply(rsarray, c(".row", ".col", "roi"), mean)

sds <- apply(rsarray, c("subj", "roi"), function(x) sd(x[lower.tri(x)]))
rsarray.s <- sweep(rsarray, MARGIN = 3:4, STATS = sds, FUN = "/")
D.s <- apply(rsarray.s, c(".row", ".col", "roi"), mean)

matplot(D.s[, , "V1_L"])
matplot(D.s[, , "LIPd_L"])
matplot(D.s[, , "IP1_L"])
matplot(D.s[, , "IP1_L"])

# # matplot(D[, , "smmouth"])
# matplot(D.s[, , "smmouth"])
# 
# # matplot(D[, , "aud"])
# matplot(D.s[, , "aud"])
# 
# # matplot(D[, , "fpc_R"])
# matplot(D.s[, , "fpc_R"])
# 
# matplot(D.s[, , "dmfc_R"])
# matplot(D.s[, , "dmfc_L"])
# 
# matplot(D.s[, , "lppc_L"])
# matplot(D.s[, , "lppc_R"])
# 
# matplot(D.s[, , "dlpfc_R"])
# matplot(D.s[, , "dlpfc_L"])
# 
# matplot(D.s[, , "ins_L"])
# matplot(D.s[, , "ins_R"])
# 
# matplot(D.s[, , "V1"])
