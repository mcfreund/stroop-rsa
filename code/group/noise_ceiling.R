#+ noise-ceiling_setup, include = FALSE

if (interactive()) {
  
  source(here::here("code", "packages.R"))
  source(here("code", "strings.R"))
  source(here("code", "funs.R"))
  source(here("code", "plots.R"))
  source(here("code", "read_atlases.R"))
  
  rsarray.rank <- read_simil_mats(tfrm = "residual-rank")  ## for noise ceiling
  
}


#` ## estimate noise ceilings
#+ noise-ceiling_estimate, fig.height = 7, fig.width = 9

ceilings <- lapply(dimnames(rsarray.rank)$roi, noise.ceiling, x = rsarray.rank)
names(ceilings) <- dimnames(rsarray.rank)$roi
ceilings %<>% bind_rows(.id = "roi")

ceilings <- ceilings %>% 
  
  mutate(
    region = toupper(gsub("_L|_R", "", roi)), 
    hemi = ifelse(grepl("_L", roi), "L", "R"),
    region.hemi = paste0(region, "_", hemi)
  )


## quick plot of all parcels


ceilings %>%
  
  filter(region %in% c("DLPFC", "DMFC", "LPPC")) %>%
  
  {
    grid.arrange(
      
      
      ggplot(., aes(roi, lb, color = region)) +
        stat_summary(fun.data = mean_cl_boot, size = 2) +
        scale_color_viridis_d() +
        theme(legend.position = "none") +
        labs(title = "lower bound", x = "parcel", y = "across-subject mean"),
      
      
      ggplot(., aes(roi, ub, color = region)) +
        stat_summary(fun.data = mean_cl_boot, size = 2) +
        scale_color_viridis_d() +
        annotate(
          geom = "text", x = -Inf, y = 0.25, label = "DLPFC", color = colors.region["DLPFC"],
          hjust = 0, vjust = 1, size = rel(6)
        ) +
        annotate(
          geom = "text", x = -Inf, y = 0.24, label = "LPPC", color = colors.region["LPPC"],
          hjust = 0, vjust = 1, size = rel(6)
        ) +
        annotate(
          geom = "text", x = -Inf, y = 0.23, label = "MFC", color = colors.region["MFC"],
          hjust = 0, vjust = 1, size = rel(6)
        ) +
        theme(legend.position = "none") +
        labs(title = "upper bound", x = "parcel", y = "across-subject mean"),
      
      ncol = 2,
      
      top = textGrob("noise ceilings", gp = gpar(fontsize = 24))
      
    )
  }


## extract ROIs

ceilings.rois <- ceilings %>% 
  filter(region.hemi %in% c("DLPFC_L", "DLPFC_R", "LPPC_L", "DMFC_L")) %>%
  group_by(subj, region) %>%
  summarize_if(is.numeric, mean)  ## average across hemisphere


means.ceilings.rois <- ceilings.rois %>%
  
  group_by(region) %>%
  
  { 
    bind_rows(
      
      summarize(., res = list(boot_mean_ci(lb))) %>% 
        tidyr::unnest(cols = c(res)) %>%
        mutate(ceiling = "lower"),
      
      summarize(., res = list(boot_mean_ci(ub))) %>% 
        tidyr::unnest(cols = c(res)) %>%
        mutate(ceiling = "upper")
      
    )
  }



#` ## contrast noise ceilings
#+ noise-ceiling_contrast, fig.height = 7, fig.width = 9

## mean values

means.ceilings.contr.rois <- ceilings.rois %>%
  
  mutate(ub = atanh(ub), lb = atanh(lb)) %>%
  
  group_by(region) %>%
  
  { 
    bind_rows(
      
      
      tidyr::pivot_wider(., -ub, names_from = "region", values_from = "lb") %>%
        
        transmute(
          DLPFC.DMFC  = DLPFC - DMFC,
          LPPC.DMFC   = LPPC - DMFC,
          DLPFC.LPPC = LPPC - DLPFC
        ) %>%
        summarize(
          DLPFC.DMFC  = list(boot_mean_ci(DLPFC.DMFC)),
          LPPC.DMFC   = list(boot_mean_ci(LPPC.DMFC)),
          DLPFC.LPPC = list(boot_mean_ci(DLPFC.LPPC)),
        ) %>% 
        tidyr::unnest(cols = c(DLPFC.DMFC, LPPC.DMFC, DLPFC.LPPC), names_sep = ".") %>%
        mutate(ceiling = "lower"),
      
      
      tidyr::pivot_wider(., -lb, names_from = "region", values_from = "ub") %>%
        
        transmute(
          DLPFC.DMFC  = DLPFC - DMFC,
          LPPC.DMFC   = LPPC - DMFC,
          DLPFC.LPPC = LPPC - DLPFC
        ) %>%
        summarize(
          DLPFC.DMFC  = list(boot_mean_ci(DLPFC.DMFC)),
          LPPC.DMFC   = list(boot_mean_ci(LPPC.DMFC)),
          DLPFC.LPPC = list(boot_mean_ci(DLPFC.LPPC)),
        ) %>% 
        tidyr::unnest(cols = c(DLPFC.DMFC, LPPC.DMFC, DLPFC.LPPC), names_sep = ".") %>%
        mutate(ceiling = "upper")
      
    )
    
  }

## get p-values for contrasts....
res.ceilings.contr.rois <- 
  
  boot(
    
    data = 
      
      ceilings.rois %>%
      mutate(ub = atanh(ub), lb = atanh(lb)) %>%
      group_by(region) %>%
      tidyr::pivot_wider(names_from = "region", values_from = c("lb", "ub")) %>%
      transmute(
        
        DLPFC.DMFC_lb  = lb_DLPFC - lb_DMFC,
        LPPC.DMFC_lb   = lb_LPPC - lb_DMFC,
        DLPFC.LPPC_lb = lb_LPPC - lb_DLPFC,
        
        DLPFC.DMFC_ub  = ub_DLPFC - ub_DMFC,
        LPPC.DMFC_ub   = ub_LPPC - ub_DMFC,
        DLPFC.LPPC_ub = ub_LPPC - ub_DLPFC
      ),
    
    statistic = function(x, ii) colMeans(x[ii, ]), 
    
    R = 1E4
    
  )

p.ceilings.contr.rois <- 
  lapply(1:6, ci2p, bootobj = res.ceilings.contr.rois) %>% 
  do.call(rbind, .) %>% 
  as.data.frame %>% setNames("p")
p.ceilings.contr.rois$contrast <- c("dlpfc.dmfc", "lppc.dmfc", "dlpfc.lppc", "dlpfc.dmfc", "lppc.dmfc", "dlpfc.lppc")
p.ceilings.contr.rois$ceiling <- rep(c("lb", "ub"), each = 3)

#` ## test for equivalence
#' * two one-sided test for equivalence

#+ noise-ceiling_tost, fig.height = 7, fig.width = 9


#` ## plot & save
#+ noise-ceiling_plot, fig.height = 7, fig.width = 9

p.means.ceilings <- means.ceilings.rois %>%
  
  filter(ceiling == "upper") %>%
  
  ggplot(aes(region, y)) +
  
  geom_line(
    data = ceilings.rois,
    aes(x = region, y = ub, group = subj),
    alpha = 0.1, size = 1,
    inherit.aes = FALSE
  ) +
  
  geom_point(size = 8) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, size = 3) +
  
  # scale_color_manual(values = colors.region) +
  
  annotate(
    geom = "text", x = 1.5, y = 0.4, size = rel(8), vjust = 0, color = "grey40", fontface = "italic",
    label = paste0(
      "p = ", 
      p.ceilings.contr.rois %>% filter(contrast == "dlpfc.dmfc", ceiling == "ub") %>% pull("p") %>% round(2)
    )
  ) +
  annotate(geom = "segment", x = 1, xend = 3, y = 0.395, yend = 0.395, color = "grey40", size = 1) +
  
  annotate(
    geom = "text", x = 2.5, y = 0.35, size = rel(8), vjust = 0, color = "grey40", fontface = "italic",
    label = paste0(
      "p = ", 
      p.ceilings.contr.rois %>% filter(contrast == "lppc.dmfc", ceiling == "ub") %>% pull("p") %>% round(2)
    )
  ) +
  annotate(geom = "segment", x = 2, xend = 3, y = 0.345, yend = 0.345, color = "grey40", size = 1) +
  
  annotate(
    geom = "text", x = 1.5, y = 0.3, size = rel(8), vjust = 0, color = "grey40", fontface = "italic",
    label = paste0(
      "p = ", 
      p.ceilings.contr.rois %>% filter(contrast == "dlpfc.lppc", ceiling == "ub") %>% pull("p") %>% round(2)
    )
  ) +
  annotate(geom = "segment", x = 1, xend = 2, y = 0.295, yend = 0.295, color = "grey40", size = 1) +
  
  labs(y = bquote("mean noise ceiling"), x = "region") +
  
  theme(
    legend.position = "none",
    panel.grid      = element_blank(),
    panel.border    = element_blank(),
    axis.line.y     = element_line(size = rel(2)),
    axis.text       = element_text(size = rel(3)),
    axis.ticks.y    = element_line(size = rel(2)),
    axis.ticks.x    = element_blank(),
    axis.title    = element_text(size = rel(4))
  )

p.means.ceilings

fwrite(p.ceilings.contr.rois, here("out", "group", "ceilings.txt"))
fwrite(means.ceilings.contr.rois, here("out", "group", "ceilings_means.txt"))
