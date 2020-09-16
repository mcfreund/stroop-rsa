#+ setup, include = FALSE

knitr::opts_chunk$set(
  fig.align = 'center', fig.width = 11.5, fig.fullwidth = TRUE, cache = TRUE
)

source(here::here("code", "packages.R"))
source(here("code", "strings.R"))
source(here("code", "funs.R"))
source(here("code", "plots.R"))
source(here("code", "read_atlases.R"))


hyps <- combo_paste(c("dlpfc", "lppc", "dmfc"), c("R", "L"), c("target", "incongruency"))
hyps.alt <- combo_paste(c("dlpfc.alt", "lppc", "dmfc.alt"), c("R", "L"), c("target", "incongruency"))

# behav.mod.objs <- readRDS(here("out", "behav", "mod_objs.RDS"))  ## behavioral model objects

blups <- 
  bind_rows(
    read.csv(here("out", "behav", "stroop_blups_rt_group201902.csv"), stringsAsFactors = FALSE),
    read.csv(here("out", "behav", "stroop_blups_rt_group201902_validation.csv"), stringsAsFactors = FALSE) 
  )

stats.subjs <- bind_rows(
  mmp         = 
    read_subj_stats(
      subjs = c(subjs.analysis, subjs.validation), 
      roi.set = "mmp"
      ) %>% 
    mutate(roi.hemi = roi, roi = gsub("_L$|_R$", "", roi)),
  superparcel = read_subj_stats(subjs = c(subjs.analysis, subjs.validation)),
  gordon      =     
    read_subj_stats(
      subjs = c(subjs.analysis, subjs.validation), 
      roi.set = "gordon"
    ) %>% 
    mutate(roi.hemi = roi, roi = gsub("_L$|_R$", "", roi)),
  .id = "scheme"
  )

## wrange data ----

stats.subjs <- stats.subjs[stats.subjs$param %in% c("target", "distractor", "incongruency"), ]
stats.subjs <- full_join(blups, stats.subjs, by = "subj")
stats.subjs %<>% ungroup %>% mutate(id = paste0(roi.hemi, "_", param))
stats.subjs %<>% group_by(scheme, id) %>% mutate(beta.s = c(scale(beta)))  ## scale betas

## to long-form data-frame

d.mmp <- stats.subjs %>% filter(scheme == "mmp")
d.super <- stats.subjs %>% filter(scheme == "superparcel")
d.gordon <- stats.subjs %>% filter(scheme == "gordon")

## to wide-form matrix

w.mmp <- stats.subjs %>% 
  
  filter(scheme == "mmp") %>%
  
  dplyr::select(subj, is.analysis.group, congr, stroop, beta, id) %>%
  pivot_wider(names_from = "id", values_from = "beta")


w.super <- stats.subjs %>% 
  
  filter(scheme == "superparcel") %>%
  
  dplyr::select(subj, is.analysis.group, congr, stroop, beta, id) %>%
  pivot_wider(names_from = "id", values_from = "beta")

w.gordon <- stats.subjs %>% 
  
  filter(scheme == "gordon") %>%
  
  dplyr::select(subj, is.analysis.group, congr, stroop, beta, id) %>%
  pivot_wider(names_from = "id", values_from = "beta")


are.identical <- identical(w.mmp[c("subj", "is.analysis.group")], w.super[c("subj", "is.analysis.group")])
are.identical <- 
  are.identical && identical(w.mmp[c("subj", "is.analysis.group")], w.gordon[c("subj", "is.analysis.group")])
if (!are.identical) {
  stop("something wrong")
} else {
  
  ids <- w.mmp[c("subj", "is.analysis.group")]
  
  m.mmp    <- w.mmp %>% ungroup %>% dplyr::select(-subj, -is.analysis.group, -scheme) %>% as.matrix
  m.super  <- w.super %>% ungroup %>%  dplyr::select(-subj, -is.analysis.group, -scheme) %>% as.matrix
  m.gordon <- w.gordon %>% ungroup %>%  dplyr::select(-subj, -is.analysis.group, -scheme) %>% as.matrix
  
}

## response variables

strooprt <- as.matrix(w.super[ids$is.analysis.group, "stroop"])  ## analysis set
strooprt_vset <- m.super[!ids$is.analysis.group, "stroop"]  ## validation set

fit1.het.trim <- readRDS(here("out", "behav", "fit1-het-trim_group201902.RDS"))  ## behavioral model
d.dissoc.hlm <- full_join(
  fit1.het.trim$data %>% filter(!is.far.out), 
  w.super %>% filter(is.analysis.group) %>% ungroup %>% dplyr::select(subj, one_of(hyps)),
  by = "subj"
)

#+