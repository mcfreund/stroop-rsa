library(here)
library(knitr)
library(magrittr)
library(dplyr)
library(data.table)
library(ggplot2)
library(mikeutils)
library(doParallel)
library(foreach)
library(gifti)
library(fANCOVA)
library(viridis)
library(colorspace)
library(boot)
library(vegan)
library(grid)
library(gridExtra)
library(lme4)
source(here("code", "strings.R"))
source(here("code", "funs.R"))
source(here("code", "read_atlases.R"))

## underlay images

hcp <- list(
  L  = readGIfTI(
    file.path(dir.atlas, "surf", "HCP_S1200_GroupAvg_v1", "S1200.L.very_inflated_MSMAll.32k_fs_LR.surf.gii")
  ),
  R = readGIfTI(
    file.path(dir.atlas, "surf", "HCP_S1200_GroupAvg_v1", "S1200.R.very_inflated_MSMAll.32k_fs_LR.surf.gii")
  )
)

## MMP template for defining masks

mmp <- list(
  L = read_gifti2matrix(file.path(dir.atlas, "MMP_surface", "mmpL.func.gii")) %>% c,
  R = read_gifti2matrix(file.path(dir.atlas, "MMP_surface", "mmpR.func.gii")) %>% c
)

mmp$L[mmp$L > 0] <- mmp$L[mmp$L > 0] - 180
mmp$R[mmp$R > 0] <- mmp$R[mmp$R > 0] + 180

## variables

params.interest <- c("target", "distractor", "incongruency")
colors.profs <- c(
  incon = "#d95f02ff", targt = "#1b9e77ff", distr = "#7570b3ff", targt.incon = "#26190dff", targt.distr = "#4682b4ff"
  )


## data

stats.subjs.tdic <- fread(
  here("out", "rsa", "stats", paste0("subjs_pro_bias_acc-only_mmp_pearson_residual_glm-tdic.csv"))
)
stats.subjs.tdic <- stats.subjs.tdic[is.analysis.group == TRUE & y == "rank", ]  ## EXCLUDE HELD OUT SUBJECTS!
stats.subjs.tdic <- stats.subjs.tdic[, "coef" := NULL]
stats.subjs.tdic %<>% full_join(atlas.key$mmp, by = "roi")


stats.group.tdic <- fread(here("out", "rsa", "stats", "group_pro_bias_acc-only_mmp_pearson_residual_glm-tdic.csv"))
stats.group.tdic <- stats.group.tdic[measure == "beta" & y == "rank" & param %in% params.interest, ]
stats.group.tdic <- stats.group.tdic[, c("measure", "y") := NULL]


## 'coding' results ----

rois.targt <- stats.group.tdic[param == "target" & p.fdr < 0.05, roi]
rois.distr <- stats.group.tdic[param == "distractor" & p.fdr < 0.05, roi]
rois.incon <- stats.group.tdic[param == "incongruency" & p.fdr < 0.05, roi]
rois.congr <- stats.group.tdic[param == "congruency" & p.fdr < 0.05, roi]

keep.these.bois <- c("m", "v", "p.fdr", "num.roi", "roi", "hemi")
stats.coding.targt <- stats.group.tdic[param == "target", ..keep.these.bois]
stats.coding.distr <- stats.group.tdic[param == "distractor", ..keep.these.bois]
stats.coding.incon <- stats.group.tdic[param == "incongruency", ..keep.these.bois]
stats.coding.congr <- stats.group.tdic[param == "congruency", ..keep.these.bois]

## intersections

rois.targt.incon <- intersect(rois.targt, rois.incon)
rois.targt.distr <- intersect(rois.targt, rois.distr)
rois.distr.incon <- intersect(rois.distr, rois.incon)  ## none
rois.targt.distr.incon <- Reduce(intersect, list(rois.targt, rois.incon, rois.distr))  ## none


## pairwise comparisons ----

## get p values

stats.pairs.tdic <- stats.subjs.tdic %>%
  filter(param %in% params.interest) %>%
  group_by(roi) %>%
  summarize(
    out = list(
      pairwise.wilcox.test(beta, param, paired = TRUE, alternative = "two.sided", p.adjust.method = "fdr")
    )
  )
pvals <- vapply(stats.pairs.tdic$out, function(.) .$p.value[lower.tri(diag(2), diag = TRUE)], numeric(3))
pvals <- t(pvals)
colnames(pvals) <- c("p.incon.distr", "p.targt.distr", "p.targt.incon")
stats.pairs.tdic <- data.frame(roi = stats.pairs.tdic$roi, pvals, stringsAsFactors = FALSE)

## get mean diffs

stats.pairs.tdic %<>% full_join(
  stats.subjs.tdic %>%
    filter(param %in% params.interest) %>%
    group_by(roi) %>%
    summarize(
      b.incon.distr = tanh(mean(atanh(beta[param == "incongruency"]) - atanh(beta[param == "distractor"]))),
      b.targt.distr = tanh(mean(atanh(beta[param == "target"]) - atanh(beta[param == "distractor"]))),
      b.targt.incon = tanh(mean(atanh(beta[param == "target"]) - atanh(beta[param == "incongruency"]))),
    ),
  by = "roi"
)
stats.pairs.tdic %<>% full_join(atlas.key$mmp, by = "roi") %>% as.data.table


## results

rois.pref.targt.vs.distr <- stats.pairs.tdic[roi %in% rois.targt & p.targt.distr < 0.05 & b.targt.distr > 0, roi]
rois.pref.targt.vs.incon <- stats.pairs.tdic[roi %in% rois.targt & p.targt.incon < 0.05 & b.targt.incon > 0, roi]
rois.pref.distr.vs.targt <- stats.pairs.tdic[roi %in% rois.distr & p.targt.distr < 0.05 & b.targt.distr < 0, roi]
rois.pref.distr.vs.incon <- stats.pairs.tdic[roi %in% rois.distr & b.incon.distr < 0.05 & b.incon.distr < 0, roi]
rois.pref.incon.vs.distr <- stats.pairs.tdic[roi %in% rois.incon & b.incon.distr < 0.05 & b.incon.distr > 0, roi]
rois.pref.incon.vs.targt <- stats.pairs.tdic[roi %in% rois.incon & p.targt.incon < 0.05 & b.targt.incon < 0, roi]

rois.pref.targt <- intersect(rois.pref.targt.vs.incon, rois.pref.targt.vs.distr)
rois.pref.distr <- intersect(rois.pref.distr.vs.targt, rois.pref.distr.vs.incon)
rois.pref.incon <- intersect(rois.pref.incon.vs.targt, rois.pref.incon.vs.distr)

## equivalency tests ----

fname.es <- here("out", "rsa", "smallest-effect-size-of-interest.csv")

if (!file.exists(fname.es)) {

  powersim <- function(es, n = 49, a = 0.05, nsim = 5000) {

    out <- numeric(nsim)

    for (sim.i in seq_len(nsim)) {
      x <- rnorm(n, es)
      p <- wilcox.test(x, alternative = "greater")$p.value
      out[sim.i] <- p < a
    }

    sum(out) / nsim

  }

  es <- seq(0, 0.5, 0.001)  ## range of effect sizes to simulate

  n.cores <- detectCores()
  cl <- makeCluster(n.cores - 1)
  registerDoParallel(cl)
  res <- foreach(ii = es, .combine = c) %dopar% {
    set.seed(ii * 1000)  ## set seed anew each iteraion (on each worker)
    powersim(ii)  ## power
  }
  stopCluster(cl)

  sesoi <- data.frame(effect.size = es, hit.rate = res)

  write.csv(sesoi, fname.es, row.names = FALSE)

} else {

  sesoi <- read.csv(fname.es)

}

## fit loess to generate smooth predicted values:
fit <- loess.as(sesoi$hit.rate, sesoi$effect.size)  ## uses model selection procedure to select optimal smoothing param
plot(sesoi$hit.rate, sesoi$effect.size, main = "effect size (cohen's d) by power level of signed rank")
lines(sesoi$hit.rate, fit$fitted, col = "firebrick")

(sesoi8 <- predict(fit, 0.8))  ## minimum effect size of interest


stats.tosts <- stats.subjs.tdic %>%
  filter(param %in% params.interest) %>%
  group_by(param, roi) %>%
  summarize(
    du = wilcox.test(beta / sd(beta), alternative = "less", mu = sesoi8)$p.value,  ## less than upward bound?
    dl = wilcox.test(beta / sd(beta), alternative = "greater", mu = -sesoi8)$p.value,  ## greater than lower bound?
    p = max(du, dl)  ## take largest p-value
  ) %>%
  as.data.table

rois.tost.targt.vs.distr <- stats.tosts[roi %in% rois.pref.targt.vs.distr & param == "distractor", ] %>%
  ungroup %>%
  mutate(p.fdr = p.adjust(p, method = "fdr"), p.holm = p.adjust(p, method = "holm")) %>%
  filter(p.fdr < 0.05) %>%
  .$roi

rois.tost.targt.vs.incon <- stats.tosts[roi %in% rois.pref.targt.vs.incon & param == "incongruency", ] %>%
  ungroup %>%
  mutate(p.fdr = p.adjust(p, method = "fdr"), p.holm = p.adjust(p, method = "holm")) %>%
  filter(p.fdr < 0.05) %>%
  .$roi

rois.target.selective <- intersect(rois.tost.targt.vs.incon, rois.tost.targt.vs.distr)

rois.tost.distr.vs.incon <- stats.tosts[roi %in% rois.pref.distr.vs.incon & param == "incongruency", ] %>%
  ungroup %>%
  mutate(p.fdr = p.adjust(p, method = "fdr"), p.holm = p.adjust(p, method = "holm")) %>%
  filter(p.fdr < 0.05) %>%
  .$roi

rois.tost.incon.vs.distr <- stats.tosts[roi %in% rois.pref.incon.vs.distr & param == "distractor", ] %>%
  ungroup %>%
  mutate(p.fdr = p.adjust(p, method = "fdr"), p.holm = p.adjust(p, method = "holm")) %>%
  filter(p.fdr < 0.05) %>%
  .$roi

rois.tost.incon.vs.targt <- stats.tosts[roi %in% rois.pref.incon.vs.targt & param == "target", ] %>%
  ungroup %>%
  mutate(p.fdr = p.adjust(p, method = "fdr"), p.holm = p.adjust(p, method = "holm")) %>%
  filter(p.fdr < 0.05) %>%
  .$roi

stats.group.tdic.coding <- stats.group.tdic[roi %in% unique(c(rois.distr, rois.incon, rois.pref.targt.vs.distr))]






# ## define profiles
# 
# profiles <- list(
#   distr       = rois.distr,
#   incon       = rois.incon,
#   targtdistr  = rois.targt.distr,
#   targtincon  = rois.targt.incon,
#   targt0distr = rois.tost.targt.vs.distr,
#   targt0incon = rois.tost.targt.vs.incon
# )
# 
# profile.distr       <- profiles$distr %>% setdiff(unlist(profiles[-1]))
# profile.incon       <- profiles$incon %>% setdiff(unlist(profiles[-2]))
# profile.targt0distr <- profiles$targt0distr
# profile.targt0incon <- profiles$targt0incon
# profile.targt.incon <- rois.targt.incon
# profile.targt.distr <- rois.targt.distr
# 
# 

## figure ----


## (a) venn diagram ----

length(rois.distr)
length(rois.targt)
length(rois.incon)
length(intersect(rois.incon, rois.targt))
length(intersect(rois.distr, rois.targt))

rois.coding <- list(
  distr = setdiff(rois.distr, c(rois.targt, rois.incon)),
  targt = setdiff(rois.targt, c(rois.distr, rois.incon)),
  incon = setdiff(rois.incon, c(rois.targt, rois.distr)),
  targt.distr = intersect(rois.targt, rois.distr),
  targt.incon = intersect(rois.targt, rois.incon)
)


## (b) brains ----

## reshape to wide, add profiles

stats.group.tdic.coding.wide <- stats.group.tdic %>%
  select(param, m, roi, community.short) %>%
  filter(param %in% params.interest)
  tidyr::spread(param, m) %>%
  filter(roi %in% unlist(rois.coding)) %>%
  mutate(
    prof = ifelse(
      roi %in% rois.coding$distr, "distr",
      ifelse(
        roi %in% rois.coding$incon, "incon",
        ifelse(
          roi %in% rois.coding$targt, "targt",
          ifelse(
            roi %in% rois.coding$targt.distr, "targt.distr",
            ifelse(
              roi %in% rois.coding$targt.incon, "targt.incon",
              NA
            )
          )
        )
      )
    )
  )

stats.subjs.tdic.coding <- stats.subjs.tdic %>%
  filter(roi %in% unlist(rois.coding), param != "congruency") %>%
  select(subj, param, beta, roi, community.short)

overlay <- rbind(
  stats.group.tdic.coding.wide %>% select(roi, prof),
  data.frame(
    roi = setdiff(atlas.key$mmp$roi, unlist(rois.coding)), 
    prof = "none", stringsAsFactors = FALSE
  )
)
overlay$profnum <- as.numeric(factor(overlay$prof, levels = c("none", unique(stats.group.tdic.coding.wide$prof)))) - 1
# overlay %<>% mutate(is.left = grepl("_L$", roi)) %>% arrange(is.left, roi.num)  ## weirdness with MMP atlas! flip hemis.

inds.left <- atlas.key$mmp$roi %>% grep("_L$", .)
inds.right <- atlas.key$mmp$roi %>% grep("_R$", .)
atlas.key$mmp$roi <- atlas.key$mmp$roi[c(inds.right, inds.left)]

overlay$roi.num <- match(overlay$roi, atlas.key$mmp$roi)
overlay %<>% arrange(roi.num)

cifti.parcellate(name.atlas = "glasser", dir.to.write = here("out", "figs"))
cifti.convert(
  fname.overlay = "group_coding",
  values.overlay = overlay$profnum,
  dir.template = here("out", "figs"),
  name.atlas = "glasser",
  dir.to.write = here("out", "figs", "ms_v1_2020-03", "group")
)




## (c) MDS ----


rsarray <- readRDS(
  here("out", "rsa", "obsv", "rsarray_pro_bias_acc-only_mmp_pearson_residual-linear.rds")
  )[, , unique(stats.subjs.tdic$subj), ]
# dimnames(rsarray)

## define examples

eg <- list(
  
  distr = stats.group.tdic.coding %>%
    filter(roi %in% rois.coding$distr, param == "distractor") %>%
    filter(v == max(v)) %>%
    pull("roi"),
  
  targt.distr = stats.group.tdic.coding %>%
    filter(roi %in% rois.coding$targt.distr) %>%
    select(param, roi, v) %>%
    tidyr::spread(param, v) %>%
    mutate(targt.distr = distractor + target) %>%
    filter(targt.distr == max(targt.distr)) %>%
    pull("roi"),
  
  targt.incon = stats.group.tdic.coding %>%
    filter(roi %in% rois.coding$targt.incon) %>%
    select(param, roi, v) %>%
    tidyr::spread(param, v) %>%
    mutate(targt.incon = incongruency + target) %>%
    filter(targt.incon == max(targt.incon)) %>%
    pull("roi"),

  incon = stats.group.tdic.coding %>%
    filter(roi %in% rois.coding$incon, param == "incongruency") %>%
    filter(v == max(v)) %>%
    pull("roi"),
  
  targt = stats.group.tdic.coding %>%
    filter(roi %in% rois.coding$targt, param == "target") %>%
    filter(v == max(v)) %>%
    pull("roi")
)


## mds

eg.mds <- lapply(eg, nmds, arr = rsarray)

## plots
p.mds <- lapply(eg.mds, plot.mds)
# p.mds

## add lines

mds.lines <- list(
  incon = list(
    ## congruency
    mds.line("whiteBLUE", "whiteRED"),
    mds.line("whiteRED", "redBLUE"),
    mds.line("redBLUE", "whitePURPLE"),
    mds.line("whitePURPLE", "blueWHITE"),
    mds.line("blueWHITE", "purpleWHITE"),
    mds.line("purpleWHITE", "redPURPLE"),
    mds.line("redPURPLE", "bluePURPLE"),
    mds.line("bluePURPLE", "redWHITE"),
    mds.line("redWHITE", "purpleBLUE"),
    mds.line("purpleBLUE", "purpleRED"),
    mds.line("purpleRED", "blueRED"),
    mds.line("blueRED", "whiteBLUE")
  ),
  # congr = list(
  #   ## congruency
  #   mds.line("whiteWHITE", "redRED"),
  #   mds.line("redRED", "blueBLUE"),
  #   mds.line("blueBLUE", "purplePURPLE"),
  #   mds.line("purplePURPLE", "whiteWHITE")
  # ),
  targt.distr = list(
    ## distractor coding
    ## BLUE
    mds.line("blueBLUE", "purpleBLUE"),
    mds.line("purpleBLUE", "whiteBLUE"),
    mds.line("whiteBLUE", "redBLUE"),
    mds.line("redBLUE", "blueBLUE"),
    ## WHITE
    mds.line("blueWHITE", "purpleWHITE"),
    mds.line("purpleWHITE", "whiteWHITE"),
    mds.line("whiteWHITE", "redWHITE"),
    mds.line("redWHITE", "blueWHITE"), 
    ## RED
    mds.line("blueRED", "whiteRED"),
    mds.line("whiteRED", "purpleRED"),
    mds.line("purpleRED", "redRED"),
    mds.line("redRED", "blueRED"),
    ## PURPLE
    mds.line("bluePURPLE", "redPURPLE"),
    mds.line("redPURPLE", "whitePURPLE"), 
    mds.line("whitePURPLE", "purplePURPLE"),
    mds.line("purplePURPLE", "bluePURPLE")
  ),
  distr = list(
    ## distractor coding
    ## BLUE
    mds.line("blueBLUE", "purpleBLUE"),
    mds.line("purpleBLUE", "whiteBLUE"),
    mds.line("whiteBLUE", "redBLUE"),
    mds.line("redBLUE", "blueBLUE"),
    ## WHITE
    mds.line("blueWHITE", "purpleWHITE"),
    mds.line("purpleWHITE", "whiteWHITE"),
    mds.line("whiteWHITE", "redWHITE"),
    mds.line("redWHITE", "blueWHITE"), 
    ## RED
    mds.line("blueRED", "whiteRED"),
    mds.line("whiteRED", "purpleRED"),
    mds.line("purpleRED", "redRED"),
    mds.line("redRED", "blueRED"),
    ## PURPLE
    mds.line("bluePURPLE", "redPURPLE"),
    mds.line("redPURPLE", "purplePURPLE"), 
    mds.line("purplePURPLE", "whitePURPLE"),
    mds.line("whitePURPLE", "bluePURPLE")
  ),
  targt.incon = list(
    ## congruency & target
    ## blue
    mds.line("blueBLUE", "blueWHITE"),
    mds.line("blueBLUE", "blueRED"),
    mds.line("blueBLUE", "bluePURPLE"),
    ## white
    mds.line("whiteWHITE", "whiteBLUE"),
    mds.line("whiteWHITE", "whiteRED"),
    mds.line("whiteWHITE", "whitePURPLE"),
    ## red
    mds.line("redRED", "redBLUE"),
    mds.line("redRED", "redWHITE"),
    mds.line("redRED", "redPURPLE"),
    ## purple
    mds.line("purplePURPLE", "purpleBLUE"),
    mds.line("purplePURPLE", "purpleWHITE"),
    mds.line("purplePURPLE", "purpleRED")
  ),
  targt = list(
    ## target
    ## blue
    mds.line("blueRED", "blueBLUE"),
    mds.line("blueBLUE", "blueWHITE"),
    mds.line("blueWHITE", "bluePURPLE"),
    mds.line("bluePURPLE", "blueRED"),
    ## white
    mds.line("whiteWHITE", "whitePURPLE"),
    mds.line("whitePURPLE", "whiteBLUE"),
    mds.line("whiteBLUE", "whiteRED"),
    mds.line("whiteRED", "whiteWHITE"),
    ## red
    mds.line("redRED", "redBLUE"),
    mds.line("redBLUE", "redWHITE"),
    mds.line("redWHITE", "redPURPLE"),
    mds.line("redPURPLE", "redRED"),
    ## purple
    mds.line("purplePURPLE", "purpleRED"),
    mds.line("purpleRED", "purpleBLUE"),
    mds.line("purpleBLUE", "purpleWHITE"),
    mds.line("purpleWHITE", "purplePURPLE")
  )
)

## get grobs

grobs <- list(
    targt       = plot.mds(eg.mds$targt, .add.lines = mds.lines$targt),
    distr       = plot.mds(eg.mds$distr, .add.lines = mds.lines$distr),
    incon       = plot.mds(eg.mds$incon, .add.lines = mds.lines$incon),
    targt.incon = plot.mds(eg.mds$targt.incon, .add.lines = mds.lines$targt.incon),
    targt.distr = plot.mds(eg.mds$targt.distr, .add.lines = mds.lines$targt.distr)
  ) %>%
  lapply(ggplotGrob)

# labs.profile <- c(
#   "target, distractor-insensitive",
#   "target & distractor",
#   "distractor",
#   "conflict",
#   "target, conflict-insensitive",
#   "target & conflict"
# )
labs.roi <- unlist(eg)

## create layout

radius.rel <- 0.1  ## relative size of each MDS circle (used to build positions of circle centers)
radius.abs <- unit(2, "cm")  ## absolute size of circle
a <- 1.2  ## scale factor overall plot
b <- 1.4  ## scale factor of MDS plot relative to surrounding circle (can be expanded)
h <- 10
w <- 10

coordinates <- data.frame(
  x = seq(radius.rel, radius.rel * 5, by = radius.rel),
  y = c(rep(c(sqrt(3) * radius.rel * 2, sqrt(3) * radius.rel), 2), sqrt(3) * radius.rel * 2)
)
coordinates <- as.data.frame(scale(coordinates, scale = FALSE) + 0.5)  ## center
coordinates <- coordinates / a
coordinates <- as.data.frame(scale(coordinates, scale = FALSE) + 0.5)  ## center

## order everything
mds.order <- c("distr", "targt.distr", "targt", "targt.incon", "incon")
colors.profs <- colors.profs[mds.order]
grobs <- grobs[mds.order]

## draw

cairo_pdf(
  here("out", "figs", "ms_v1_2020-03", "group", "mds_from_R.pdf"), 
  height = h, width = h
)

# plot.new()

grid.circle(coordinates$x, coordinates$y, r = radius.abs, gp = gpar(lwd = 4, fill = "transparent", col = colors.profs))

vps <- vector("list", length(grobs))
for (ii in seq_along(grobs)) {
  vp <- viewport(x = coordinates$x[ii], y = coordinates$y[ii], width = radius.abs * b, height = radius.abs * b)
  pushViewport(vp)
  grid.draw(grobs[[ii]])
  upViewport()
}

dev.off()


## (d) session comparison ----



## see fig_group_session_comparison.R

## (e) venn 2 ----

length(rois.tost.targt.vs.incon)
length(rois.tost.targt.vs.distr)



## (f) brains ----


rois.targt0distr <- setdiff(rois.tost.targt.vs.distr, rois.tost.targt.vs.incon)
rois.targt0incon <- setdiff(rois.tost.targt.vs.incon, rois.tost.targt.vs.distr)
rois.targt.selec <- intersect(rois.tost.targt.vs.incon, rois.tost.targt.vs.distr)

stats.group.tdic.coding.wide <- stats.group.tdic.coding %>%
  select(param, m, roi, community.short) %>%
  tidyr::spread(param, m) %>%
  filter(roi %in% c(rois.tost.targt.vs.distr, rois.tost.targt.vs.incon)) %>%
  mutate(
    prof = ifelse(
      roi %in% rois.targt0distr, "targt0distr",
      ifelse(
        roi %in% rois.targt0incon, "targt0incon",
        ifelse(
          roi %in% rois.targt.selec, "targt.selec",
          NA
          )
        )
    )
  )

overlay <- rbind(
  stats.group.tdic.coding.wide %>% select(roi, prof),
  data.frame(
    roi = setdiff(atlas.key$mmp$roi, c(rois.tost.targt.vs.distr, rois.tost.targt.vs.incon)),
    prof = "none", stringsAsFactors = FALSE
    )
)
overlay$profnum <- as.numeric(factor(overlay$prof, levels = c("none", unique(stats.group.tdic.coding.wide$prof)))) - 1

inds.left <- atlas.key$mmp$roi %>% grep("_L$", .)
inds.right <- atlas.key$mmp$roi %>% grep("_R$", .)
atlas.key$mmp$roi <- atlas.key$mmp$roi[c(inds.right, inds.left)]

overlay$roi.num <- match(overlay$roi, atlas.key$mmp$roi)
overlay %<>% arrange(roi.num)

cifti.parcellate(name.atlas = "glasser", dir.to.write = here("out", "figs"))
cifti.convert(
  fname.overlay = "group_selective_v1",
  values.overlay = overlay$profnum,
  dir.template = here("out", "figs"),
  name.atlas = "glasser",
  dir.to.write = here("out", "figs", "ms_v1_2020-03", "group")
)




## SCRATCH -----



# length(rois.tost.targt.vs.distr)
# length(rois.tost.targt.vs.incon)
# length(intersect(rois.tost.targt.vs.distr, rois.tost.targt.vs.incon))
# 
# length(profile.targt.incon)



## see inkscape


## (c) profile scatterplot ----
# 
# 
# ## reshape to wide, add profiles
# 
# stats.group.tdic.coding.wide <- stats.group.tdic.coding %>%
#   select(param, m, roi, community.short) %>%
#   tidyr::spread(param, m) %>%
#   filter(roi %in% unlist(profiles)) %>%
#   mutate(
#     prof = ifelse(
#       roi %in% profile.distr, "distr",
#       ifelse(
#         roi %in% profile.incon, "incon",
#         ifelse(
#           roi %in% profile.targt0distr, "targt0distr",
#           ifelse(
#             roi %in% profile.targt0incon, "targt0incon",
#             NA
#           )
#         )
#       )
#     ),
#     prof = ifelse(roi %in% profile.targt.incon, "targt.incon", prof),
#     prof = ifelse(roi %in% profile.targt.distr, "targt.distr", prof)
#   )
# 
# stats.subjs.tdic.coding <- stats.subjs.tdic %>%
#   filter(roi %in% unlist(profiles), param != "congruency") %>%
#   select(subj, param, beta, roi, community.short)
# 
# ## bootstrap
# 
# bootmean <- function(data, indices) mean(data[indices])
# n.samp <- 1E4
# 
# g <- interaction(stats.subjs.tdic.coding$roi, stats.subjs.tdic.coding$param)
# stats.subjs.tdic.coding$id <- g
# l <- split(stats.subjs.tdic.coding, g)
# 
# stats.group.tdic.coding.wide$lb.incongruency <- NA
# stats.group.tdic.coding.wide$ub.incongruency <- NA
# stats.group.tdic.coding.wide$lb.target <- NA
# stats.group.tdic.coding.wide$ub.target <- NA
# stats.group.tdic.coding.wide$lb.distractor <- NA
# stats.group.tdic.coding.wide$ub.distractor <- NA
# 
# for (roi.i in seq_along(l)) {
#   
#   l.i <- l[[roi.i]]
#   ci <- boot.ci(boot(l.i$beta, bootmean, R = n.samp), conf = 0.96)$bca[4:5]
#   
#   
#   cols <- paste0(c("lb.", 'ub.'), unique(l.i$param))
#   stats.group.tdic.coding.wide[
#     stats.group.tdic.coding.wide$roi == gsub(".target|.distractor|.incongruency", "", names(l)[roi.i]),
#     cols
#     ] <- ci
#   
#   # print(roi.i / length(l))
#   
# }
# 
# ## tables for plot
# 
# note.profs <- data.frame(
#   prof = c("targt0incon", "targt0distr", "incon", "targt.incon", "distr", "targt.distr"),
#   lab = c(
#     "target coding, distractor-insensitive", 
#     "target coding, conflict-insentive",
#     "conflict",
#     "target & conflict",
#     "distractor",
#     "target & distractor"
#   )
# )
# 
# colors.profs <- c(
#   targt0incon = '#a6cee3',
#   targt0distr = '#1f78b4',
#   incon = '#b2df8a',
#   targt.incon = '#33a02c',
#   distr = '#fb9a99',
#   targt.distr = '#e31a1c'
# )
# 
# plot.range <- stats.group.tdic.coding.wide %>% 
#   summarize(
#     max.targt = max(target),
#     min.targt = min(target),
#     max.distr = max(distractor),
#     min.distr = min(distractor),
#     max.incon = max(incongruency),
#     min.incon = min(incongruency)
#   )
# 
# labs.incon.targt <- c(
#   "3b_L", "4_L", "A5_R", "IP1_L", "LIPd_L", "MIP_L", "LO2_R", "9-46d_L", "6ma_L", "RSC_L", "IP0_R", "8BM_R", "s6-8_R",
#   "FEF_R", "PBelt_R"
#   )
# 
# labs.targt.distr <- c(
#   "3b_L", "4_L", "V1_L", "V2_L", "9-46d_L", "LO2_R", "IP0_R", "LO2_R", "9-46d_L", "IP1_L", "MIP_L", "6ma_L", "LIPd_L"
# )
# 
# 
# ## plot
# 
# 
# p.incon.targt <- stats.group.tdic.coding.wide %>%
#   ggplot(aes(incongruency, target, color = prof)) +
#   scale_color_manual(values = colors.profs) +
#   theme(
#     legend.position = "none",
#     panel.grid = element_blank(),
#     panel.background = element_blank()
#   ) +
#   geom_errorbar(aes(ymin = lb.target, ymax = ub.target), width = 0, alpha = 0.3) +
#   geom_errorbarh(aes(xmin = lb.incongruency, xmax = ub.incongruency, alpha = 0.3)) +
#   geom_point(size = 4) +
#   # geom_text(
#   #   data = note.profs[note.profs$prof == "targt0incon", ], aes(label = lab),
#   #   color = colors.profs["targt0incon"], 
#   #   x = -Inf, y = 0.01, fontface = "bold", vjust = 0, hjust = 0
#   # ) +
#   # geom_text(
#   #   data = note.profs[note.profs$prof == "targt0distr", ], aes(label = lab),
#   #   color = colors.profs["targt0distr"], 
#   #   x = -Inf, y = 0.005, fontface = "bold", vjust = 0, hjust = 0
#   # ) +
#   # geom_text(
#   #   data = note.profs[note.profs$prof == "incon", ], aes(label = lab),
#   #   color = colors.profs["incon"], 
#   #   x = -Inf, y = 0, fontface = "bold", vjust = 0, hjust = 0
#   # ) +
#   # geom_text(
#   #   data = note.profs[note.profs$prof == "targt.incon", ], aes(label = lab),
#   #   color = colors.profs["targt.incon"], 
#   #   x = -Inf, y = -0.005, fontface = "bold", vjust = 0, hjust = 0
#   # ) +
#   # geom_text(
#   #   data = note.profs[note.profs$prof == "distr", ], aes(label = lab),
#   #   color = colors.profs["distr"], 
#   #   x = -Inf, y = -0.01 , fontface = "bold", vjust = 0, hjust = 0
#   # ) +
#   # geom_text(
#   #   data = note.profs[note.profs$prof == "targt.distr", ], aes(label = lab),
#   #   color = colors.profs["targt.distr"], 
#   #   x = -Inf, y = -0.015 , fontface = "bold", vjust = 0, hjust = 0
#   # ) +
#   # geom_text(aes(label = roi), fontface = "bold", nudge_x = -0.005, nudge_y = 0.001) +
#   geom_text(
#     aes(label = ifelse(roi %in% labs.incon.targt, roi, "")),
#     fontface = "bold", nudge_x = -0.01, nudge_y = 0.004
#     ) +
#   scale_x_continuous(breaks = c(0, plot.range$max.incon) %>% round(2)) +
#   scale_y_continuous(breaks = c(0, plot.range$max.targt) %>% round(2)) +
#   geom_segment(
#     aes(y = 0, yend = plot.range$max.targt, x = -Inf, xend = -Inf),
#     color = "grey50", size = 0.25
#   ) +
#   geom_segment(
#     aes(x = 0, xend = plot.range$max.incon, y = -Inf, yend = -Inf),
#     color = "grey50", size = 0.25
#   )
# 
# p.distr.targt <- stats.group.tdic.coding.wide %>%
#   ggplot(aes(distractor, target, color = prof)) +
#   scale_color_manual(values = colors.profs) +
#   theme(
#     legend.position = "none",
#     panel.grid = element_blank(),
#     panel.background = element_blank(),
#     axis.text.y = element_blank(),
#     axis.title.y = element_blank(),
#     axis.ticks.y = element_blank()
#   ) +
#   geom_errorbar(aes(ymin = lb.target, ymax = ub.target), width = 0, alpha = 0.3) +
#   geom_errorbarh(aes(xmin = lb.distractor, xmax = ub.distractor), alpha = 0.3) +
#   geom_point(size = 3) +
#   geom_text(
#     aes(label = ifelse(roi %in% labs.targt.distr, roi, "")),
#     fontface = "bold", nudge_x = 0.005, nudge_y = -0.0025
#     ) +
#   # geom_text(aes(label = roi), fontface = "bold") +
#   scale_x_continuous(breaks = c(0, plot.range$max.distr) %>% round(2)) +
#   # scale_y_continuous(breaks = c(0, plot.range$max.targt) %>% round(2)) +
#   # geom_segment(
#   #   aes(y = 0, yend = plot.range$max.targt, x = -Inf, xend = -Inf),
#   #   color = "grey50", size = 0.25
#   #   ) +
#   geom_segment(
#     aes(x = 0, xend = plot.range$max.distr, y = -Inf, yend = -Inf),
#     color = "grey50", size = 0.25
#   )
# 
# p.profiles <- gridExtra::grid.arrange(p.incon.targt, p.distr.targt, ncol = 2)
# 
# ggsave(
#   here("out", "figs", "ms_v1_2020-03", "group_profiles", "profiles_from_R.pdf"), 
#   plot = p.profiles,
#   units = "cm",
#   device = "pdf",
#   height = 15,
#   width = 30
# )



## (b) brains ----
# 
# overlay <- rbind(
#   stats.group.tdic.coding.wide %>% select(roi, prof),
#   data.frame(
#     roi = setdiff(atlas.key$mmp$roi, unlist(profiles)), 
#     prof = "none", stringsAsFactors = FALSE
#     )
# )
# overlay$profnum <- as.numeric(factor(overlay$prof, levels = c("none", unique(stats.group.tdic.coding.wide$prof)))) - 1
# overlay %<>% mutate(is.left = grepl("_L$", roi)) %>% arrange(is.left, roi.num)  ## weirdness with MMP atlas! flip hemis.
# 
# inds.left <- atlas.key$mmp$roi %>% grep("_L$", .)
# inds.right <- atlas.key$mmp$roi %>% grep("_R$", .)
# atlas.key$mmp$roi <- atlas.key$mmp$roi[c(inds.right, inds.left)]
# 
# overlay$roi.num <- match(overlay$roi, atlas.key$mmp$roi)
# overlay %<>% arrange(roi.num)
# 
# cifti.parcellate(name.atlas = "glasser", dir.to.write = here("out", "figs"))
# cifti.convert(
#   fname.overlay = "group_profiles_v1",
#   values.overlay = overlay$profnum,
#   dir.template = here("out", "figs"),
#   name.atlas = "glasser",
#   dir.to.write = here("out", "figs", "ms_v1_2020-03", "group_profiles")
# )
# 
# 
# ## (d) MDS ----
# 
# rsarray <- readRDS(
#   here("out", "rsa", "obsv", "rsarray_pro_bias_acc-only_mmp_pearson_residual-linear.rds")
#   )[, , unique(stats.subjs.tdic$subj), ]
# # dimnames(rsarray)
# 
# ## define examples
# 
# eg <- list(
#   distr = stats.group.tdic.coding %>% 
#     filter(roi %in% profile.distr, param == "distractor") %>%
#     filter(v == max(v)) %>%
#     pull("roi"),
#   targt.distr = stats.group.tdic.coding %>% 
#     filter(roi %in% profile.targt.distr) %>%
#     select(param, roi, v) %>%
#     tidyr::spread(param, v) %>%
#     mutate(targt.distr = distractor + target) %>%
#     filter(targt.distr == max(targt.distr)) %>%
#     pull("roi"),
#   incon = stats.group.tdic.coding %>% 
#     filter(roi %in% profile.incon, param == "incongruency") %>%
#     filter(v == max(v)) %>%
#     pull("roi"),
#   targt0distr = stats.group.tdic.coding %>% 
#     filter(roi %in% profile.targt0distr, param == "target") %>%
#     filter(v == max(v)) %>%
#     pull("roi"),
#   targt0incon = stats.group.tdic.coding %>% 
#     filter(roi %in% profile.targt0incon, param == "target") %>%
#     filter(v == max(v)) %>%
#     pull("roi"),
#   targt.incon = stats.group.tdic.coding %>% 
#     filter(roi %in% profile.targt.incon) %>%
#     select(param, roi, v) %>%
#     tidyr::spread(param, v) %>%
#     mutate(targt.incon = incongruency + target) %>%
#     filter(targt.incon == max(targt.incon)) %>%
#     pull("roi")
# )


## mds
# 
# eg.mds <- lapply(eg, nmds, arr = rsarray)
# 
# ## plots
# # p.mds <- lapply(eg.mds, plot.mds)
# # p.mds
# 
# ## add lines
# 
# mds.lines <- list(
#   incon = list(
#     ## congruency
#     mds.line("whiteBLUE", "whiteRED"),
#     mds.line("whiteRED", "redBLUE"),
#     mds.line("redBLUE", "whitePURPLE"),
#     mds.line("whitePURPLE", "blueWHITE"),
#     mds.line("blueWHITE", "purpleWHITE"),
#     mds.line("purpleWHITE", "redPURPLE"),
#     mds.line("redPURPLE", "bluePURPLE"),
#     mds.line("bluePURPLE", "redWHITE"),
#     mds.line("redWHITE", "purpleBLUE"),
#     mds.line("purpleBLUE", "purpleRED"),
#     mds.line("purpleRED", "blueRED"),
#     mds.line("blueRED", "whiteBLUE")
#   ),
#   # congr = list(
#   #   ## congruency
#   #   mds.line("whiteWHITE", "redRED"),
#   #   mds.line("redRED", "blueBLUE"),
#   #   mds.line("blueBLUE", "purplePURPLE"),
#   #   mds.line("purplePURPLE", "whiteWHITE")
#   # ),
#   targt.distr = list(
#     ## distractor coding
#     ## BLUE
#     mds.line("blueBLUE", "purpleBLUE"),
#     mds.line("purpleBLUE", "whiteBLUE"),
#     mds.line("whiteBLUE", "redBLUE"),
#     mds.line("redBLUE", "blueBLUE"),
#     ## WHITE
#     mds.line("blueWHITE", "purpleWHITE"),
#     mds.line("purpleWHITE", "whiteWHITE"),
#     mds.line("whiteWHITE", "redWHITE"),
#     mds.line("redWHITE", "blueWHITE"), 
#     ## RED
#     mds.line("blueRED", "whiteRED"),
#     mds.line("whiteRED", "purpleRED"),
#     mds.line("purpleRED", "redRED"),
#     mds.line("redRED", "blueRED"),
#     ## PURPLE
#     mds.line("bluePURPLE", "redPURPLE"),
#     mds.line("redPURPLE", "whitePURPLE"), 
#     mds.line("whitePURPLE", "purplePURPLE"),
#     mds.line("purplePURPLE", "bluePURPLE")
#   ),
#   distr = list(
#     ## distractor coding
#     ## BLUE
#     mds.line("blueBLUE", "purpleBLUE"),
#     mds.line("purpleBLUE", "whiteBLUE"),
#     mds.line("whiteBLUE", "redBLUE"),
#     mds.line("redBLUE", "blueBLUE"),
#     ## WHITE
#     mds.line("blueWHITE", "purpleWHITE"),
#     mds.line("purpleWHITE", "whiteWHITE"),
#     mds.line("whiteWHITE", "redWHITE"),
#     mds.line("redWHITE", "blueWHITE"), 
#     ## RED
#     mds.line("blueRED", "whiteRED"),
#     mds.line("whiteRED", "purpleRED"),
#     mds.line("purpleRED", "redRED"),
#     mds.line("redRED", "blueRED"),
#     ## PURPLE
#     mds.line("bluePURPLE", "redPURPLE"),
#     mds.line("redPURPLE", "purplePURPLE"), 
#     mds.line("purplePURPLE", "whitePURPLE"),
#     mds.line("whitePURPLE", "bluePURPLE")
#   ),
#   targt.incon = list(
#     ## congruency & target
#     ## blue
#     mds.line("blueBLUE", "blueWHITE"),
#     mds.line("blueBLUE", "blueRED"),
#     mds.line("blueBLUE", "bluePURPLE"),
#     ## white
#     mds.line("whiteWHITE", "whiteBLUE"),
#     mds.line("whiteWHITE", "whiteRED"),
#     mds.line("whiteWHITE", "whitePURPLE"),
#     ## red
#     mds.line("redRED", "redBLUE"),
#     mds.line("redRED", "redWHITE"),
#     mds.line("redRED", "redPURPLE"),
#     ## purple
#     mds.line("purplePURPLE", "purpleBLUE"),
#     mds.line("purplePURPLE", "purpleWHITE"),
#     mds.line("purplePURPLE", "purpleRED")
#   ),
#   targt0incon = list(
#     ## target
#     ## blue
#     mds.line("blueRED", "blueBLUE"),
#     mds.line("blueBLUE", "blueWHITE"),
#     mds.line("blueWHITE", "bluePURPLE"),
#     mds.line("bluePURPLE", "blueRED"),
#     ## white
#     mds.line("whiteWHITE", "whitePURPLE"),
#     mds.line("whitePURPLE", "whiteBLUE"),
#     mds.line("whiteBLUE", "whiteRED"),
#     mds.line("whiteRED", "whiteWHITE"),
#     ## red
#     mds.line("redRED", "redBLUE"),
#     mds.line("redBLUE", "redWHITE"),
#     mds.line("redWHITE", "redPURPLE"),
#     mds.line("redPURPLE", "redRED"),
#     ## purple
#     mds.line("purplePURPLE", "purpleRED"),
#     mds.line("purpleRED", "purpleBLUE"),
#     mds.line("purpleBLUE", "purpleWHITE"),
#     mds.line("purpleWHITE", "purplePURPLE")
#   ),
#   targt0distr = list(
#     ## target
#     ## blue
#     mds.line("blueRED", "blueBLUE"),
#     mds.line("blueBLUE", "blueWHITE"),
#     mds.line("blueWHITE", "bluePURPLE"),
#     mds.line("bluePURPLE", "blueRED"),
#     ## white
#     mds.line("whiteWHITE", "whitePURPLE"),
#     mds.line("whitePURPLE", "whiteBLUE"),
#     mds.line("whiteBLUE", "whiteRED"),
#     mds.line("whiteRED", "whiteWHITE"),
#     ## red
#     mds.line("redRED", "redBLUE"),
#     mds.line("redBLUE", "redWHITE"),
#     mds.line("redWHITE", "redPURPLE"),
#     mds.line("redPURPLE", "redRED"),
#     ## purple
#     mds.line("purplePURPLE", "purpleRED"),
#     mds.line("purpleRED", "purpleWHITE"),
#     mds.line("purpleWHITE", "purpleBLUE"),
#     mds.line("purpleBLUE", "purplePURPLE")
#   )
# )
# 
# ## save (finish in inkscape)
# # 
# # p.mds <- grid.arrange(
# #   plot.mds(eg.mds$targt0distr, .add.lines = mds.lines$targt0distr, .size = rel(1.75)),
# #   plot.mds(eg.mds$targt.distr, .add.lines = mds.lines$targt.distr, .size = rel(1.75)),
# #   plot.mds(eg.mds$distr,       .add.lines = mds.lines$distr,        .size = rel(1.75)),
# #   plot.mds(eg.mds$incon,       .add.lines = mds.lines$incon,        .size = rel(1.75)),
# #   plot.mds(eg.mds$targt0incon, .add.lines = mds.lines$targt0incon, .size = rel(1.75)),
# #   plot.mds(eg.mds$targt.incon, .add.lines = mds.lines$targt.incon, .size = rel(1.75)),
# #   ncol = 3
# # )
# # 
# # ggsave(
# #   here("out", "figs", "ms_v1_2020-03", "group_profiles", "mds_from_R.pdf"), 
# #   plot = p.mds,
# #   units = "cm",
# #   device = "pdf",
# #   height = 10,
# #   width = 15
# # )
# # 
# 
# ## get grobs
# 
# grobs <- list(
#     targt0incon = plot.mds(eg.mds$targt0incon, .add.lines = mds.lines$targt0incon),
#     distr       = plot.mds(eg.mds$distr,       .add.lines = mds.lines$distr),
#     incon       = plot.mds(eg.mds$incon,       .add.lines = mds.lines$incon),
#     targt.incon = plot.mds(eg.mds$targt.incon, .add.lines = mds.lines$targt.incon),
#     targt.distr = plot.mds(eg.mds$targt.distr, .add.lines = mds.lines$targt.distr),
#     targt0distr = plot.mds(eg.mds$targt0distr, .add.lines = mds.lines$targt0distr)
#   ) %>%
#   lapply(ggplotGrob)
# 
# # labs.profile <- c(
# #   "target, distractor-insensitive",
# #   "target & distractor",
# #   "distractor",
# #   "conflict",
# #   "target, conflict-insensitive",
# #   "target & conflict"
# # )
# labs.roi <- unlist(eg)
# 
# ## create layout
# 
# radius.rel <- 0.1  ## relative size of each MDS circle (used to build positions of circle centers)
# radius.abs <- unit(2, "cm")  ## absolute size of circle
# a <- 1.2  ## scale factor overall plot
# b <- 1.4  ## scale factor of MDS plot relative to surrounding circle (can be expanded)
# h <- 10
# w <- 10
# 
# coordinates <- data.frame(
#   x = c(seq(radius.rel, radius.rel * 5, by = radius.rel), radius.rel * 3),
#   y = c(rep(c(sqrt(3) * radius.rel * 2, sqrt(3) * radius.rel), 2), sqrt(3) * radius.rel * 2, 0)
# )
# coordinates <- as.data.frame(scale(coordinates, scale = FALSE) + 0.5)  ## center
# coordinates <- coordinates / a
# coordinates <- as.data.frame(scale(coordinates, scale = FALSE) + 0.5)  ## center
# 
# ## order everything
# mds.order <- c("targt0incon", "targt.distr", "distr", "targt.incon", "incon", "targt0distr")
# colors.profs <- colors.profs[mds.order]
# grobs <- grobs[mds.order]
# 
# ## draw
# 
# cairo_pdf(
#   here("out", "figs", "ms_v1_2020-03", "group_profiles", "mds_from_R.pdf"), 
#   height = h, width = h
#   )
# 
# # plot.new()
# 
# grid.circle(coordinates$x, coordinates$y, r = radius.abs, gp = gpar(lwd = 4, fill = "transparent", col = colors.profs))
# 
# vps <- vector("list", length(grobs))
# for (ii in seq_along(grobs)) {
#   vp <- viewport(x = coordinates$x[ii], y = coordinates$y[ii], width = radius.abs * b, height = radius.abs * b)
#   pushViewport(vp)
#   grid.draw(grobs[[ii]])
#   upViewport()
# }
# 
# dev.off()


## draw labels
# circ.lab <- circ / 1.25  ## for labels
# grid.text(
#   labs.name, x = circ.lab$x, y = circ.lab$y, rot = c(0, -72, 36, -36, 72),
#   gp = gpar(fontface = "bold.italic", fontsize = 8, col = hues)
# )
# 
# grid.text(
#   labs.symb,
#   x = convertWidth(unit(circ.lab$x, "npc"), "npc"),
#   y = convertHeight(unit(circ.lab$y, "npc"), "npc"),
#   rot = c(0, -72, 36, -36, 72),
#   gp = gpar(fontface = "bold", fontsize = 6, col = hues, fontfamily = "Lucida Sans Unicode")
# )
# 
# grid.text(
#   titles, x = title.pos.x, y = title.pos.y,
#   gp = gpar(fontface = "bold.italic", fontsize = 6, col = hues)
# )



## visual stats ----
# 
# stats.subjs.tdic.masks <- fread(
#   here("out", "rsa", "stats", paste0("subjs_pro_bias_acc-only_masks_pearson_residual_glm-tdic.csv"))
# )
# stats.subjs.tdic.masks <- stats.subjs.tdic.masks[is.analysis.group == TRUE & y == "rank", ]  ## EXCLUDE HELD OUT SUBJECTS!
# stats.subjs.tdic.masks <- stats.subjs.tdic.masks[, "coef" := NULL]
# 
# stats.subjs.vis <- stats.subjs.tdic %>% 
#   filter(
#     roi %in% c("V1_L", "V2_L", "V8_L", "VVC_L"), 
#     param %in% c("target", "distractor")
#     ) %>%
#   full_join(
#     stats.subjs.tdic.masks %>%
#       filter(roi == "vwfa", param %in% c("target", "distractor"))
#   )
# 
# 
# stats.group.vis <- stats.subjs.vis %>%
#   group_by(param, roi) %>%
#   summarize(
#     p = wilcox.test(beta, alternative = "greater")$p.value,
#     beta = mean(beta)
#   )
#   
# 
# p.vis <- stats.subjs.vis %>%
#   ggplot(aes(roi, beta, color = param)) +
#   stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.25)) +
#   geom_line(
#     data = stats.group.vis,
#     aes(group = param), position = position_dodge(width = 0.25), linetype = "dashed"
#   )
# 
# 
# fit.v1v2 <- lmer(beta ~ roi + (1 | subj), stats.subjs.vis, subset = roi %in% c("V1_L", "V2_L") & param == "target")
# summary(fit.v1v2)
# confint(fit.v1v2)
# 
# 
# 
# 
# ## coding parcels ----
# 
# 
# 
# 
# 
# 
# 
# 
# ## hclust ----
# 
# 
# stats.subjs.tdic.coding %>%
#   filter(roi %in% profile.targt0distr)
# 
# rsarray.targ0distr <- rsarray[, , , profile.targt0distr]
# 
# for (roi.i in seq_along(profile.targt0distr)) {
#   for (subj.i in seq_along(dimnames(rsarray)$subj)) {
# 
#     rsm <- rsarray.targ0distr[, , subj.i, roi.i]
#     rsm[upper.tri(rsm, diag = TRUE)] <- NA
#     rsarray.targ0distr[, , subj.i, roi.i] <- rsm
# 
#   }
# }
# rsd.targ0distr <- reshape2::melt(rsarray.targ0distr, na.rm = TRUE)
# rsd.targ0distr <- rsd.targ0distr %>%
#   mutate(id = paste0(.row, "_", .col)) %>%
#   select(-.row, -.col)
# 
# rsdbar.targ0dist <- rsd.targ0distr %>%
#   group_by(roi, id) %>%
#   summarize(value = mean(value)) %>%
#   tidyr::spread(roi, value)
# R <- cor(rsdbar.targ0dist[, -1])
# D <- as.dist(1 - R)
# qcor(R)
# ?hclust
# 
# hclust.targ0distr <- hclust(D, method = "ward.D")
# 
# plot(hclust.targ0distr)
# 
# 
# heatmap(R)
# 
# 
# 
# 
# pca.targ0distr <- prcomp(t(scale(t(rsdbar.targ0dist[, -1]))), scale = TRUE)
# plot(pca.targ0distr)
# 
# pca <- pca.targ0distr$rotation %>%
#   as.data.frame %>%
#   # cbind(., roi = rownames(R)) %>%
#   tibble::rownames_to_column("roi") %>%
#   left_join(atlas.key$mmp, by = "roi")
#   ggplot(aes(PC1, PC2)) +
#   geom_label(aes(label = roi, color = community.short))
# str(hclust.targ0distr)
# 
# 
# 
# 
# 
# 
# 
# 
# d <- rsdbar.targ0dist[, -1]
# 
# 
# 
# 
# 
# clustk2 <- kmeans(d, centers = 4)
# plot(clustk2)
# 
# 
# 
# 
# setdiff(profile.targt.incon, profile.targt0distr)
# 
# 
# install.packages("kernlab")
# library(kernlab)
# 
# sc <- specc(scale(d), centers = 5)
# plot(pca[, c("PC1", "PC2")], col = sc, pch = 16)            # estimated classes (x)
# text(pca[, c("PC1", "PC2")], labels = pca$roi, col = sc)
# 
# 
# km <- kmeans(scale(d), centers = 4)
# plot(pca[, c("PC1", "PC2")], col = km$cluster, pch = 16)            # estimated classes (x)
# text(pca[, c("PC1", "PC2")], labels = pca$roi, col = km$cluster)
# 
# 
# 
# intersect(profile.targt.incon, profile.targt0distr)
# 
# 
# 
# 
# 
# 
