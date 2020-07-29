library(here)
library(magrittr)
library(dplyr)
library(data.table)
library(ggplot2)
library(mikeutils)
library(grid)
library(gridExtra)
source(here("code", "strings.R"))
source(here("code", "funs.R"))
source(here("code", "read_atlases.R"))

## variables

params.interest <- c("target", "distractor", "incongruency")
colors.model <- c(incongruency = "#d95f02", target = "#1b9e77", distractor = "#7570b3")

## data

stats.subjs.bas <- fread(
  here("out", "rsa", "stats", paste0("subjs_bas_bias_acc-only_downsamp_mmp_pearson_residual_glm-tdic.csv"))
)
stats.subjs.pro <- fread(
  here("out", "rsa", "stats", paste0("subjs_pro_bias_acc-only_downsamp_mmp_pearson_residual_glm-tdic.csv"))
)

stats.subjs.bas <- stats.subjs.bas[is.analysis.group == TRUE & y == "rank", ]  ## EXCLUDE HELD OUT SUBJECTS!
stats.subjs.bas <- stats.subjs.bas[, "coef" := NULL]
stats.subjs.bas %<>% full_join(atlas.key$mmp, by = "roi")

stats.subjs.pro <- stats.subjs.pro[is.analysis.group == TRUE & y == "rank", ]  ## EXCLUDE HELD OUT SUBJECTS!
stats.subjs.pro <- stats.subjs.pro[, "coef" := NULL]
stats.subjs.pro %<>% full_join(atlas.key$mmp, by = "roi")


stats.group.bas <- fread(
  here("out", "rsa", "stats", paste0("group_bas_bias_acc-only_downsamp_mmp_pearson_residual_glm-tdic.csv"))
)
stats.group.pro <- fread(
  here("out", "rsa", "stats", paste0("group_pro_bias_acc-only_downsamp_mmp_pearson_residual_glm-tdic.csv"))
)

stats.group.bas <- stats.group.bas[y == "rank" & param %in% params.interest, ]
stats.group.pro <- stats.group.pro[y == "rank" & param %in% params.interest, ]



## proportions ----

rois.targt <- stats.group.pro[param == "target" & p < 0.05, roi]
rois.distr <- stats.group.pro[param == "distractor" & p < 0.05, roi]
rois.incon <- stats.group.pro[param == "incongruency" & p < 0.05, roi]

n.identified <- length(c(rois.targt, rois.distr, rois.incon))


## using original FDR correction within models

counts.bas <- stats.group.bas %>%
  top_n(n.identified, -p.fdr) %>%
  group_by(param) %>%
  summarize(counts = n(), props = counts / n.identified) %>%
  mutate(session = "bas")

counts.pro <- stats.group.pro %>%
  top_n(n.identified, -p.fdr) %>%
  group_by(param) %>%
  summarize(counts = n(), props = counts / n.identified) %>%
  mutate(session = "pro")

counts.orig <- bind_rows(counts.pro, counts.bas)

p.orig <- counts.orig %>%
  mutate(session = relevel(as.factor(session), "pro")) %>%
  ggplot(aes(session, props, fill = param)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors.model) +
  theme_bw(base_size = 12) +
  geom_text(
    data = counts.orig[counts.orig$param == "target" & counts.orig$session == "bas", ], 
    aes(label = paste0("target: ", round(props * 100), "%")), y = 0.05,
    color = "white", fontface = "bold"
  ) +
  geom_text(
    data = counts.orig[counts.orig$param == "incongruency" & counts.orig$session == "bas", ], 
    aes(label = paste0("conflict: ", round(props * 100), "%")), y = 0.5,
    color = "white", fontface = "bold"
  ) +
  geom_text(
    data = data.frame(param = "distractor", props = 0, session = "bas"), 
    aes(label = paste0("distractor: ", round(props * 100), "%")), y = 1.02,
    color = colors.model["distractor"], fontface = "bold"
  ) +
  geom_text(
    data = counts.orig[counts.orig$param == "target" & counts.orig$session == "pro", ], 
    aes(label = paste0("target: ", round(props * 100), "%")), y = 0.5,
    color = "white", fontface = "bold"
  ) +
  geom_text(
    data = counts.orig[counts.orig$param == "incongruency" & counts.orig$session == "pro", ], 
    aes(label = paste0("conflict: ", round(props * 100), "%")), y = 0.97,
    color = "white", fontface = "bold"
  ) +
  geom_text(
    data = counts.orig[counts.orig$param == "distractor" & counts.orig$session == "pro", ], 
    aes(label = paste0("distractor: ", round(props * 100), "%")), y = 1.02,
    color = colors.model["distractor"], fontface = "bold"
  ) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y =  element_blank(),
    legend.position = "none"
  ) +
  labs(y = "% of parcels identified") +
  scale_x_discrete(labels = c(bas = "mostly congruent\n(alt. session)", pro = "mostly incongruent\n(main session)"))


# ggsave(
#   here("out", "figs", "ms_v1_2020-03", "group", "group_session_comparison.pdf"), 
#   plot = p.orig,
#   units = "cm",
#   device = "pdf",
#   height = 10,
#   width = 7.5
# )


## across all p-values (all models)

counts.all.bas <- stats.group.bas %>%
  top_n(n.identified, -p) %>%
  group_by(param) %>%
  summarize(counts = n(), props = counts / n.identified) %>%
  mutate(session = "bas")

counts.all.pro <- stats.group.pro %>%
  top_n(n.identified, -p) %>%
  group_by(param) %>%
  summarize(counts = n(), props = counts / n.identified) %>%
  mutate(session = "pro")

counts.all <- bind_rows(counts.all.bas, counts.all.pro)

p.all <- counts.all %>%
  ggplot(aes(session, props, fill = param)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors.model) +
  theme_bw(base_size = 12)

ggsave(
  here("out", "figs", "ms_v1_2020-03", "group", "group_session_comparison_runwise.pdf"),
  plot = p.all,
  units = "cm",
  device = "pdf",
  height = 10,
  width = 7.5
)


## proportions across entire range of alpha

alphas <- seq(0.00001, 0.1, by = 0.00001)
counts.bas  <- data.frame(
  alpha = alphas,
  targt = 0,
  distr = 0,
  incon = 0,
  n.obsv = 0
) %>% as.data.table
counts.pro  <- counts.bas

for (a.i in seq_along(alphas)) {
  # a.i = 1000
  
  a <- alphas[a.i]
  
  ## baseline
  
  observed.bas <- stats.group.bas[p < a]$param
  table.bas <- table(observed.bas)
  
  counts.bas[alpha == a, "n.obsv"] <- sum(table.bas)
  
  if ("distractor" %in% names(table.bas))   counts.bas[alpha == a, "distr"] <- table.bas["distractor"]
  if ("target" %in% names(table.bas))       counts.bas[alpha == a, "targt"] <- table.bas["target"]
  if ("incongruency" %in% names(table.bas)) counts.bas[alpha == a, "incon"] <- table.bas["incongruency"]
  
  ## proactive
  
  observed.pro <- stats.group.pro[p < a]$param
  table.pro <- table(observed.pro)
  
  counts.pro[a.i, "n.obsv"] <- sum(table.pro)
  
  if ("distractor" %in% names(table.pro))   counts.pro[alpha == a, "distr"] <- table.pro["distractor"]
  if ("target" %in% names(table.pro))       counts.pro[alpha == a, "targt"] <- table.pro["target"]
  if ("incongruency" %in% names(table.pro)) counts.pro[alpha == a, "incon"] <- table.pro["incongruency"]
  
  
}


counts <- rbind(
  cbind(counts.bas, session = "bas"),
  cbind(counts.pro, session = "pro")
)

counts <- counts %>% 
  melt(
    id.vars = c("alpha", "n.obsv", "session"), 
    measure.vars = c("targt", "distr", "incon"),
    variable.name = "param", value.name = "freq"
  )

counts$prop <- counts$freq / counts$n.obsv 

counts %>%
  ggplot(aes(alpha, prop, fill = param)) +
  geom_area() +
  facet_grid(cols = vars(session)) +
  scale_fill_brewer(type = "qual", palette = 2) +
  theme_bw(base_size = 12)
