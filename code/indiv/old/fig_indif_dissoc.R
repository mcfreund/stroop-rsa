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
library(vegan)
library(foreach)
library(doParallel)

source(here("code", "strings.R"))
source(here("code", "funs.R"))

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

## strings

colors.model <- c(incongruency = "#d95f02", target = "#1b9e77", distractor = "#7570b3")
colors.congruency <- c(I = "#d01c8b", C = "#4dac26")
colors.targets <- c(blue = "#08519c", red = "#a50f15", white = "grey50", purple = "#54278f")
  

d <- stats.subjs.tdic %>% filter(
  roi.set == "anat",
  superparcel %in% c("dlpfc_L", "dlpfc_R", "mfc_L", "mfc_R", "lppc_L", "lppc_R")
)


d %<>% mutate(id = as.character(interaction(superparcel, param)))

w <- d %>%
  ungroup %>%
  dplyr::select(subj, stroop, id, beta) %>%
  tidyr::spread(id, beta)

## (a) pca ----
# m <- w[, -c(1, 2)]
# rownames(m) <- w$subj
# 
# pca <- prcomp(scale(m))
# 
# plot(pca)
# 
# 
# p.pca <- pca$rotation[, 1:2] %>%
#   as.data.frame %>%
#   tibble::rownames_to_column("id") %>%
#   mutate(
#     roi = gsub("(.*_.).*", "\\1", id),
#     param = gsub(".*_..(.*)", "\\1", id)
#   ) %>%
#   ggplot(aes(PC1, PC2, color = param)) +
#   geom_point() +
#   geom_text(aes(label = roi, color = param), nudge_y = 0.025, nudge_x = 0, fontface = "bold", size = 2.5) +
#   geom_segment(
#     aes(x = 0, y = 0, xend = 0, yend = 0.1), 
#     arrow = arrow(length = unit(0.01, "npc"), type = "closed"), color = "grey40", size = 1
#   ) +
#   geom_segment(
#     aes(x = 0, y = 0, xend = 0.1, yend = 0),
#     arrow = arrow(length = unit(0.01, "npc"), type = "closed"), color = "grey40", size = 1
#   ) +
#   scale_color_manual(values = c(colors.model)) +
#   annotate(geom = "text", label = "PC1", x = 0.05, y = 0, vjust = 1, fontface = "bold.italic", color = "grey40") +
#   annotate(geom = "text", label = "PC2", x = 0, y = 0.05, angle = 90, vjust = 0, fontface = "bold.italic", color = "grey40") +
#   annotate(geom = "text", label = "0.1", x = -0.0025, y = 0.1, hjust = 1) +
#   annotate(geom = "text", label = "0.1", y = -0.005, x = 0.1, vjust = 1) +
#   annotate(geom = "text", label = "0", y = -0.005, x = -0.0025, vjust = 1) +
#   theme(
#     axis.ticks.x = element_blank(),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     legend.position = "none",
#     axis.ticks = element_blank(), 
#     axis.text = element_blank(),
#     axis.title = element_blank(),
#     panel.border = element_blank()
#   ) +
#   annotate(
#     geom = "text", y = 0.24, x = -Inf, color = colors.model["target"], fontface = "bold.italic", label = "target coding",
#     size = 3, hjust = 0) +
#   annotate(
#     geom = "text", y = 0.22, x = -Inf, color = colors.model["distractor"], fontface = "bold.italic", label = "distractor coding",
#     size = 3, hjust = 0) +
#   annotate(
#     geom = "text", y = 0.20, x = -Inf, color = colors.model["incongruency"], fontface = "bold.italic", label = "conflict coding",
#     size = 3, hjust = 0) +
#   labs(title = "a")
# 
# p.pca
# 

## (a) brains ----

## see write_masks.R

## (b) scatterplot ----


w$lfp_R.target <- (w$lppc_R.target + w$dlpfc_R.target) / 2

# write.csv(w, file.path("C:/Users/mcf/Box/global/docs", "papers", "tics_2020", "figs", "fig_box2", "indif_data.csv"))

p.dissoc.scatter <- w %>%
  ggplot(aes(y = stroop)) +
  stat_boot_ci(aes(x = lfp_R.target), alpha = 0.3, n = 1E3, percent = 96, fill = colors.model["target"]) +
  stat_boot_ci(aes(x = mfc_L.incongruency), alpha = 0.3, n = 1E4, percent = 96, fill = colors.model["incongruency"]) +
  geom_smooth(aes(x = lfp_R.target), alpha = 0.3, color = colors.model["target"], method = "lm", se = FALSE) +
  geom_smooth(aes(x = mfc_L.incongruency), alpha = 0.3, color = colors.model["incongruency"], method = "lm", se = FALSE) +
  geom_point(aes(x = lfp_R.target), shape = 21, fill = colors.model["target"], color = "white", size = 3) +
  geom_point(aes(x = mfc_L.incongruency), shape = 21, fill = colors.model["incongruency"], color = "white", size = 3) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_blank(),
    panel.border =  element_blank(),
    panel.background = element_blank(),
    # axis.text.y = element_blank(),
    # axis.title.y = element_blank(),
    # axis.ticks.y = element_blank()
  ) +
  scale_x_continuous(breaks = c(min(w$mfc_L.incongruency), 0, max(w$mfc_L.incongruency)) %>% round(2)) +
  scale_y_continuous(breaks = c(min(w$stroop), max(w$stroop)) %>% round(2)) +
  geom_segment(
    aes(y = round(min(w$stroop), 2), yend = round(max(w$stroop), 2), x = -Inf, xend = -Inf),
    color = "grey40", size = 1
    ) +
  geom_segment(
    aes(x = round(min(w$mfc_L.incongruency), 2), xend = round(max(w$mfc_L.incongruency), 2), y = -Inf, yend = -Inf),
    color = "grey40", size = 1
  ) +
  labs(title = "b", x = "RSA model fit")




## (c) trichotomized MDS

fname <- here("out", "indif", "mds_trichotomy_bootstrap.csv")

if (file.exists(fname)) {
  
  mrot.boot <- fread(fname)

} else {
  
  subjs <- blups$subj
  n.subj <- length(subjs)
  n.stim <- length(bias.items)
  
  ## get data
  
  R <- readRDS(here("out", "rsa", "obsv", "rsarray_pro_bias_acc-only_masks_pearson_residual-linear.rds"))
  dimnames(R)
  R <- R[, , subjs, c("anatfunc_lppc_R", "anatfunc_mfc_L", "anatfunc_dlpfc_R")]
  R <- abind::abind(
    R,
    apply(R[, , , c("anatfunc_lppc_R", "anatfunc_dlpfc_R")], 1:3, function(.) tanh(mean(atanh(.))))
  )
  names(dimnames(R)) <- c(".row", ".col", "subj", "roi")
  dimnames(R)[["roi"]] <- c("lppc_R", "mfc_L", "dlpfc_R", "lfp_R")
    
  ## mean matrices
  
  Rbar <- apply(R, c(1:2, 4), function(x) tanh(mean(atanh(x))))  ## cross-subject average RSM per region
  Dbar <- 1 - Rbar  ## same, but correlation distance
  
  rois <- dimnames(Dbar)$roi
  Mbar <- array(
    NA, 
    dim = c(stim = n.stim, dim = 2, roi = length(rois)), 
    dimnames = list(stim = bias.items, dims = c("MDS1", "MDS2"), roi = rois)
  )
  p.mds.ave <- vector("list", length(rois)) %>% setNames(rois)
  for (roi in rois) {
    
    Mbar[, , roi] <- Dbar[, , roi] %>% vegan::metaMDS(k = 2, trace = FALSE) %>% .$points
    
    p.mds.ave[[roi]] <- Dbar[, , roi] %>% 
      vegan::metaMDS(k = 2, trace = FALSE) %>% 
      .$points %>%
      as.data.frame %>%
      tibble::rownames_to_column("stim") %>%
      bind_cols(., split.str.item(.$stim)) %>% 
      ggplot(aes(MDS1, MDS2)) +
      geom_label(aes(label = word, color = color), fill = "grey60", fontface = "bold", label.size = 0) +
      scale_color_manual(values = setNames(bias.colors, bias.colors)) +
      theme(
        panel.background = element_blank(), 
        axis.text = element_blank(), 
        legend.position = "none", 
        axis.ticks = element_blank()
      ) +
      labs(title = roi)
    
  }
  grid.arrange(p.mds.ave[[1]], p.mds.ave[[2]], p.mds.ave[[3]], p.mds.ave[[4]], ncol = 4)
  
  
  topthird <- w$subj[w$stroop < quantile(w$stroop, 1/3)]
  midthird <- w$subj[w$stroop < quantile(w$stroop, 2/3) & w$stroop > quantile(w$stroop, 1/3)]
  botthird <- w$subj[w$stroop > quantile(w$stroop, 2/3)]
  
  par(mfrow = c(1, 1))
  plot(
    w$stroop,
    pch = 16,
    col = ifelse(w$subj %in% topthird, "black", ifelse(w$subj %in% midthird, "orange", "firebrick")),
    main = "stroop RTs, split into 3 quantiles"
  )
  
  
  ## resamples
  
  n.resamples <- 1E4
  l.mrot.boot <- vector("list", length(rois)) %>% setNames(rois)
  # set.seed(0)
  time.beg <- Sys.time()
  for (roi in rois) {
    
    n.cores <- detectCores()
    cl <- makeCluster(n.cores - 1)
    registerDoParallel(cl)
    
    mrot.boot <- foreach(sample.i = seq_len(n.resamples)) %dopar% {
      
      set.seed(sample.i)  ## set seed anew each iteration (on each worker)
      
      topthird.i <- sample(topthird, replace = TRUE)
      midthird.i <- sample(midthird, replace = TRUE)
      botthird.i <- sample(botthird, replace = TRUE)
      
      ## get stats
      
      D.top <- 1 - apply(R[, , topthird.i, roi], 1:2, function(x) tanh(mean(atanh(x))))
      D.mid <- 1 - apply(R[, , midthird.i, roi], 1:2, function(x) tanh(mean(atanh(x))))
      D.bot <- 1 - apply(R[, , botthird.i, roi], 1:2, function(x) tanh(mean(atanh(x))))
      
      M.top <- vegan::metaMDS(D.top, k = 2, trace = FALSE)$points
      M.mid <- vegan::metaMDS(D.mid, k = 2, trace = FALSE)$points
      M.bot <- vegan::metaMDS(D.bot, k = 2, trace = FALSE)$points
      
      mrot.boot <- rbind(
        data.frame(tri = "top", vegan::procrustes(Mbar[, , roi], M.top, scale = TRUE)$Yrot),
        data.frame(tri = "mid", vegan::procrustes(Mbar[, , roi], M.mid, scale = TRUE)$Yrot),
        data.frame(tri = "bot", vegan::procrustes(Mbar[, , roi], M.bot, scale = TRUE)$Yrot)
      )
      mrot.boot$iter <- sample.i
      
      mrot.boot
      
    }
    stopCluster(cl)
    
    l.mrot.boot[[roi]] <- do.call(rbind, mrot.boot)
    
    print(roi)
    
  }
  (time.run <- Sys.time() - time.beg)
  
  mrot.boot <- bind_rows(lapply(l.mrot.boot, tibble::rownames_to_column, var = "stim"), .id = "roi")
  mrot.boot$stim <- gsub("[0-9]", "", mrot.boot$stim)
  mrot.boot %<>% cbind(., split.str.item(.$stim))
  fwrite(mrot.boot, here("out", "indif", "mds_trichotomy_bootstrap.csv"))
  
}


## lppc
mrot.boot.4plot <- bind_rows(
  mrot.boot %>%
    filter(roi == "lfp_R") %>%
    # filter(roi %in% c("anatfunc_dlpfc_L", "anatfunc_lppc_R")) %>%
    group_by(iter, color, tri, roi) %>%
    summarize(MDS1 = mean(X1), MDS2 = mean(X2)) %>%
    rename(variable = color),
  mrot.boot %>%
    filter(roi == "mfc_L") %>%
    # filter(roi %in% c("anatfunc_dlpfc_L", "anatfunc_lppc_R")) %>%
    group_by(iter, congruency, tri, roi) %>%
    summarize(MDS1 = mean(X1), MDS2 = mean(X2)) %>%
    rename(variable = congruency)
)

p.mdstri.dot <- mrot.boot.4plot %>%
  ungroup %>%
  mutate(tri = relevel(as.factor(tri), "top")) %>%
  ggplot(aes(x = MDS1, y = MDS2, fill = variable, color = variable)) +
  # stat_ellipse(type = "norm", level = 0.96, size = 0.5) +
  geom_point(alpha = 0.01, shape = 21, color = "transparent") +
  # stat_density_2d(size = 1) +
  scale_fill_manual(values = c(colors.congruency, colors.targets)) +
  scale_color_manual(values = c(colors.congruency, colors.targets)) +
  facet_grid(
    rows = vars(tri), 
    cols = vars(roi),
    labeller = labeller(
      tri = c(top = "small stroop", mid = "mid stroop", bot = "large stroop")
    ) 
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text =  element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )

# p.mdstri.dot
# p.mdstri.dot + geom_density_2d(size = 1, alpha = 0.1)

p.mdstri <- mrot.boot.4plot %>%
  ungroup %>%
  mutate(tri = relevel(as.factor(tri), "top")) %>%
  # filter(roi == "lfp_R") %>%
  ggplot(aes(x = MDS1, y = MDS2)) +
  
  scale_alpha_continuous(range = c(0.1, 1)) +
  
  stat_density2d(
    data = mrot.boot.4plot %>% filter(variable == "blue"),
    aes(alpha = stat(nlevel)), geom = "polygon", fill = colors.targets["blue"]
  ) +
  stat_density2d(
    data = mrot.boot.4plot %>% filter(variable == "red"),
    aes(alpha = stat(nlevel)), geom = "polygon", fill = colors.targets["red"]
  ) +
  stat_density2d(
    data = mrot.boot.4plot %>% filter(variable == "purple"),
    aes(alpha = stat(nlevel)), geom = "polygon", fill = colors.targets["purple"]
  ) +
  stat_density2d(
    data = mrot.boot.4plot %>% filter(variable == "white"),
    aes(alpha = stat(nlevel)), geom = "polygon", fill = colors.targets["white"]
  ) +
  
  stat_density2d(
    data = mrot.boot.4plot %>% filter(variable == "C"),
    aes(alpha = stat(nlevel)), geom = "polygon", fill = colors.congruency["C"]
  ) +
  stat_density2d(
    data = mrot.boot.4plot %>% filter(variable == "I"),
    aes(alpha = stat(nlevel)), geom = "polygon", fill = colors.congruency["I"]
  ) +

  # scale_fill_continuous_sequential(palette = "Blues 3") +
  # stat_density2d(
  #   data = mrot.boot.4plot %>% filter(variable == "blue"),
  #   aes(fill = stat(nlevel), alpha = stat(nlevel)), geom = "polygon"
  #   ) +
  # 
  # new_scale_fill() +
  # scale_fill_continuous_sequential(palette = "Purples 3") +
  # stat_density2d(
  #   data = mrot.boot.4plot %>% filter(variable == "purple"),
  #   aes(fill = stat(nlevel), alpha = stat(nlevel)), geom = "polygon"
  # ) +
  # 
  # new_scale_fill() +
  # scale_fill_continuous_sequential(palette = "Reds 3") +
  # stat_density2d(
  #   data = mrot.boot.4plot %>% filter(variable == "red"),
  #   aes(fill = stat(nlevel), alpha = stat(nlevel)), geom = "polygon"
  # ) +
  # 
  # new_scale_fill() +
  # scale_fill_continuous_sequential(palette = "Grays") +
  # stat_density2d(
  #   data = mrot.boot.4plot %>% filter(variable == "white"),
  #   aes(fill = stat(nlevel), alpha = stat(nlevel)), geom = "polygon"
  # ) +

  # geom_point(alpha = 0.01, shape = 21, color = "transparent") +
  # scale_fill_manual(values = c(colors.congruency, colors.targets)) +
  
  # scale_color_manual(values = c(colors.congruency, colors.targets)) +
  
  facet_grid(
    rows = vars(tri), 
    cols = vars(roi),
    labeller = labeller(
      tri = c(top = "small stroop", mid = "mid stroop", bot = "large stroop")
    ),
    switch = "x"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text =  element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  ) +
  
  geom_text(
    data = data.frame(roi = "mfc_L", tri = "bot"),
    inherit.aes = FALSE, fontface = "bold.italic",
    hjust = 0, x = -Inf, y = Inf, label = "\ncongruents", color = colors.congruency["C"],
    size = 3
  ) +
  
  geom_text(
    data = data.frame(roi = "mfc_L", tri = "bot"),
    inherit.aes = FALSE, fontface = "bold.italic",
    hjust = 0, x = -Inf, y = Inf, label = "\n\nincongruents", color = colors.congruency["I"],
    size = 3
  ) +
  
  geom_text(
    data = data.frame(roi = "lfp_R", tri = "bot"),
    inherit.aes = FALSE, fontface = "bold.italic",
    hjust = 0, x = -Inf, y = Inf, label = "\nblues", color = colors.targets["blue"],
    size = 3
  ) +
  
  geom_text(
    data = data.frame(roi = "lfp_R", tri = "bot"),
    inherit.aes = FALSE, fontface = "bold.italic",
    hjust = 0, x = -Inf, y = Inf, label = "\n\npurples", color = colors.targets["purple"],
    size = 3
  ) +
  
  geom_text(
    data = data.frame(roi = "lfp_R", tri = "bot"),
    inherit.aes = FALSE, fontface = "bold.italic",
    hjust = 0, x = -Inf, y = Inf, label = "\n\n\nwhites", color = colors.targets["white"],
    size = 3
  ) +
  
  geom_text(
    data = data.frame(roi = "lfp_R", tri = "bot"),
    inherit.aes = FALSE, fontface = "bold.italic",
    hjust = 0, x = -Inf, y = Inf, label = "\n\n\n\nreds", color = colors.targets["red"],
    size = 3
  ) +
  labs(title = "c")

## put together and write ----

# p.dissoc.tryp <- grid.arrange(
#   p.pca, p.dissoc.scatter, p.mdstri, ncol = 3, widths = c(2/3, 1, 2/3)
# )
# 
# 
# ggsave(
#   here("out", "figs", "ms_v1_2020-03", "indiff_dissoc", "indif_dissoc_tryp.pdf"), 
#   plot = p.dissoc.tryp,
#   units = "cm",
#   device = "pdf",
#   height = 10,
#   width = 10 * 7/3
# )



p.dissoc <- grid.arrange(
  p.dissoc.scatter + labs(title = "a"), p.mdstri + labs(title = "b"), ncol = 2, widths = c(1, 2/3)
)

ggsave(
  here("out", "figs", "ms_v1_2020-03", "indif_dissoc", "indif_dissoc.pdf"), 
  plot = p.dissoc,
  units = "cm",
  device = "pdf",
  height = 9,
  width = 9 * 6/3
)


