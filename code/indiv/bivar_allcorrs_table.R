#' * create table of stats
#'    * table-indiv-superparcel-allcors


#+ include = FALSE
if (interactive()) source(here::here("code", "indiv", "setup.R"))
#+

#+ bivar-superparcel_correlate_all

set.seed(0)

cors <- d.super %>%
  
  filter(is.analysis.group) %>%
  
  group_by(id) %>%
  
  summarize(
    r.line = cor(beta, stroop),
    r.rank = cor(beta, stroop, method = "spearman"),
    r2.line = r.line^2,
    r2.rank = r.rank^2
  )

## get confidence intervals

l <- d.super %>%
  
  filter(is.analysis.group) %>%
  group_by(id) %>%
  split(f = .$id) %>%
  purrr::map(~.[c("stroop", "beta")])

ci.line <- lapply(l, cor_ci) %>% bind_rows
ci.rank <- lapply(l, cor_ci, method = "spearman") %>% bind_rows

names(ci.line) %<>% paste0("_line")
names(ci.rank) %<>% paste0("_rank")

cors <- bind_cols(cors, ci.line, ci.rank)
cors$p.geq0_line <- 1 - cors$p.leq0_line
cors$p.geq0_rank <- 1 - cors$p.leq0_rank

kable(cors %>% arrange(-r2.line), digits = 2)

fwrite(cors, here("out", "indiv", "allcor_bivariate.txt"))
#+