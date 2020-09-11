#+ fpc-dissoc-parcel_setup, include = FALSE

if (interactive()) source(here::here("code", "group", "mds.R"))

#+ fpc-dissoc-parcel_plot

## format

rois.lppc <- combo_paste(
  c("IP0", "IPS1", "IP1", "MIP", "7PL", "7AL", "7PC", "VIP", "LIPv", "LIPd", "AIP", "IP2"), 
  c("L", "R")
)
rois.dmfc <- combo_paste(c("SCEF", "8BM",  "p32pr", "a32pr"), c("L", "R"))
rois.dlpfc <- combo_paste(c("p9-46v", "i6-8", "8Av", "8C"), c("L", "R"))

d.mmp <- stats.subjs.mmp %>% 
  
  filter(roi %in% c(rois.dmfc, rois.dlpfc, rois.lppc)) %>%
  
  mutate(
    region = ifelse(
      roi %in% rois.dlpfc, "DLPFC", 
      ifelse(
        roi %in% rois.dmfc, "DMFC",
        ifelse(roi %in% rois.lppc, "LPPC", NA)
      )
    ),
    region.hemi = paste0(region, "_", hemi)
  )


## plot

set.seed(0)
means.mmp <- d.mmp %>%
  
  filter(param %in% c("incongruency", "target", "distractor")) %>%
  group_by(roi, param, subj) %>%
  summarize(beta = mean(beta)) %>% 
  
  group_by(roi, param) %>%
  summarize(res = list(boot_mean_ci(beta))) %>% 
  
  tidyr::unnest(cols = c(res))


p.mmp.lppc <- means.mmp %>%
  
  filter(param %in% c("incongruency", "target", "distractor"), roi %in% rois.lppc) %>%
  
  ggplot(aes(roi, y, color = param)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_point(size = geom.point.size, position = position_dodge(width = 1/2)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position = position_dodge(width = 1/2), width = 0, size = geom.line.size) +
  lemon::coord_capped_cart(left = "both") +
  scale_y_continuous(limits = c(-0.05, 0.225)) +
  scale_color_manual(values = colors.model) +
  labs(y = bquote("Model fit ("*bar(beta)*")")) +
  theme(
    legend.position = "none", 
    panel.grid      = element_blank(), 
    panel.border    = element_blank(),
    axis.line.y     = element_line(size = axis.line.size),
    axis.text       = element_text(size = axis.text.size*3/4),
    # axis.text.x     = element_text(color = colors.model),
    axis.ticks.y    = element_line(size = axis.line.size),
    axis.ticks.x    = element_blank(),
    axis.title      = element_text(size = axis.title.size),
    axis.title.x    = element_blank()
  )

p.mmp.dlpfc <- means.mmp %>%
  
  filter(param %in% c("incongruency", "target", "distractor"), roi %in% rois.dlpfc) %>%
  
  ggplot(aes(roi, y, color = param)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_point(size = geom.point.size, position = position_dodge(width = 1/2)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position = position_dodge(width = 1/2), width = 0, size = geom.line.size) +
  lemon::coord_capped_cart(left = "both") +
  scale_y_continuous(limits = c(-0.05, 0.225)) +
  scale_color_manual(values = colors.model) +
  labs(y = bquote("Model fit ("*bar(beta)*")")) +
  theme(
    legend.position = "none", 
    panel.grid      = element_blank(), 
    panel.border    = element_blank(),
    axis.line.y     = element_line(size = axis.line.size),
    axis.text       = element_text(size = axis.text.size*3/4),
    # axis.text.x     = element_text(color = colors.model),
    axis.ticks.y    = element_line(size = axis.line.size),
    axis.ticks.x    = element_blank(),
    axis.title    = element_text(size = axis.title.size),
    axis.title.x    = element_blank()
  )

p.mmp.dmfc <- means.mmp %>%
  
  filter(param %in% c("incongruency", "target", "distractor"), roi %in% rois.dmfc) %>%
  
  ggplot(aes(roi, y, color = param)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_point(size = geom.point.size, position = position_dodge(width = 1/2)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position = position_dodge(width = 1/2), width = 0, size = geom.line.size) +
  lemon::coord_capped_cart(left = "both") +
  scale_y_continuous(limits = c(-0.05, 0.225)) +
  scale_color_manual(values = colors.model) +
  labs(y = bquote("Model fit ("*bar(beta)*")")) +
  theme(
    legend.position = "none", 
    panel.grid      = element_blank(), 
    panel.border    = element_blank(),
    axis.line.y     = element_line(size = axis.line.size),
    axis.text       = element_text(size = axis.text.size*3/4),
    # axis.text.x     = element_text(color = colors.model),
    axis.ticks.y    = element_line(size = axis.line.size),
    axis.ticks.x    = element_blank(),
    axis.title    = element_text(size = axis.title.size),
    axis.title.x    = element_blank()
  )

p.parcel <- plot_grid(
  plot_grid(p.mmp.dlpfc, p.mmp.dmfc, labels = c("A", "B"), label_size = 12), 
  p.mmp.lppc, 
  nrow = 2, labels = c("", "C"), label_size = 12
)

p.parcel

ggsave(here("out", "group", "fpc_parcels.pdf"), p.parcel, device = "pdf", width = 20, height = 10, unit = "cm")


#+ fpc-dissoc-parcel_model, include = FALSE, eval = FALSE

## model ----
# 
# d.mmp.rois <- d.mmp %>% 
#   filter(region %in% c("DLPFC", "MFC", "LPPC"), param %in% c("distractor", "incongruency", "target"))
# 
# fit.mmp <- lmer(
#   beta ~ 0 + interaction(roi, param) + (region * param | subj), 
#   d.mmp.rois,
#   control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1E9))
#   # REML = FALSE
#   )
# summary(fit.mmp)
# 
# is.distractor <- grepl("distractor", names(fixef(fit.mmp)))
# is.incongruency <- grepl("incongruency", names(fixef(fit.mmp)))
# is.target <- grepl("target", names(fixef(fit.mmp)))
# 
# is.dlpfc <- grepl(paste0(rois.dlpfc, collapse = "|"), names(fixef(fit.mmp)))
# is.lppc <- grepl(paste0(rois.lppc, collapse = "|"), names(fixef(fit.mmp)))
# is.dmfc <- grepl(paste0(rois.dmfc, collapse = "|"), names(fixef(fit.mmp)))
# 
# is.left <- grepl("_L", names(fixef(fit.mmp)))
# is.right <- grepl("_R", names(fixef(fit.mmp)))
# 
# 
# contrasts.mmp <- rbind(
#   
#   ## region*model means (over hemi)
#   
#   dlpfc.target = is.dlpfc * is.target / 2,
#   lppc.target  = is.lppc * is.target / 2,
#   dmfc.target  = is.dmfc * is.target / 2,
#   
#   dlpfc.distractor = is.dlpfc * is.distractor / 2,
#   lppc.distractor  = is.lppc * is.distractor / 2,
#   dmfc.distractor  = is.dmfc * is.distractor / 2,
#   
#   dlpfc.incongruency = is.dlpfc * is.incongruency / 2,
#   lppc.incongruency  = is.lppc * is.incongruency / 2,
#   dmfc.incongruency  = is.dmfc * is.incongruency / 2,
#   
#   ## within-region contrasts (over hemi)
#   
#   "target-distractor|dlpfc" = is.dlpfc * is.target - is.dlpfc * is.distractor / 2,
#   "target-distractor|lppc"  = is.lppc * is.target - is.lppc * is.distractor / 2,
#   "target-distractor|dmfc"  = is.dmfc * is.target - is.dmfc * is.distractor / 2,
#   
#   "target-incongruency|dlpfc" = is.dlpfc * is.target - is.dlpfc * is.incongruency / 2,
#   "target-incongruency|lppc"  = is.lppc * is.target - is.lppc * is.incongruency / 2,
#   "target-incongruency|dmfc"  = is.dmfc * is.target - is.dmfc * is.incongruency / 2,
#   
#   ## within-region target-incongruency contrsts (by hemi)
# 
#   "target-incongruency|dlpfc_L" = (is.dlpfc * is.target - is.dlpfc * is.incongruency) * is.left,
#   "target-incongruency|lppc_L" = (is.lppc * is.target - is.lppc * is.incongruency) * is.left,
#   "target-incongruency|dmfc_L" = (is.dmfc * is.target - is.dmfc * is.incongruency) * is.left,
#   
#   "target-incongruency|dlpfc_R" = (is.dlpfc * is.target - is.dlpfc * is.incongruency) * is.right,
#   "target-incongruency|lppc_R" = (is.lppc * is.target - is.lppc * is.incongruency) * is.right,
#   "target-incongruency|dmfc_R" = (is.dmfc * is.target - is.dmfc * is.incongruency) * is.right,
#   
#   ## cross-region contrasts:
#   "dlpfc-mfc_L|incongruency" = (is.dlpfc/2 - (is.dmfc * is.left)) * is.incongruency,
#   "dlpfc-lppc_L|incongruency" = (is.dlpfc/2 - (is.lppc * is.left)) * is.incongruency,
#   "dlpfc-mfc_L|target" = (is.dlpfc/2 - (is.dmfc * is.left)) * is.target,
#   "dlpfc-lppc_L|target" = (is.dlpfc/2 - (is.lppc * is.left)) * is.target,
#   
#   ## interactions:
#   "(dlpfc-mfc_L)(target-incon)" = 
#     (is.dlpfc/2 - (is.dmfc * is.left)) * is.target - (is.dlpfc/2 - (is.dmfc * is.left)) * is.incongruency,
#   "(dlpfc-lppc_L)(target-incon)" = 
#     (is.dlpfc/2 - (is.lppc * is.left)) * is.target - (is.dlpfc/2 - (is.lppc * is.left)) * is.incongruency
# 
# )
# 
# (glht.mmp <- summary(glht(fit.mmp, contrasts.mmp), test = adjusted("none")))
