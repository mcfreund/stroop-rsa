#+ mds_setup, include = FALSE

if (interactive()) {
  
  source(here::here("code", "packages.R"))
  source(here("code", "strings.R"))
  source(here("code", "funs.R"))
  source(here("code", "plots.R"))
  source(here("code", "read_atlases.R"))
  
  stats.subjs.super <- read_subj_stats(roi.set = "masks")
  rsarray.line <- read_simil_mats(tfrm = "residual-line")  ## for MDS
  
}

#' ## fit

#+ mds_fit

rois.mds <- c(distr = "evis", incon = "dmfc_L", targt = "smmouth")
set.seed(0)
mds <- lapply(rois.mds, nmds, arr = rsarray.line)

#' ## plot

#+ mds_plot, fig.width = 4, fig.height = 4

## add lines

mds.lines <- list(
  
  distr = list(
    ## distractor coding
    ## BLUE
    mds.line("blueBLUE", "purpleBLUE"),
    mds.line("purpleBLUE", "redBLUE"),
    mds.line("whiteBLUE", "blueBLUE"),
    mds.line("redBLUE", "whiteBLUE"),
    ## WHITE
    mds.line("purpleWHITE", "whiteWHITE"),
    mds.line("purpleWHITE", "blueWHITE"),
    mds.line("redWHITE", "blueWHITE"),
    mds.line("whiteWHITE", "redWHITE"), 
    ## RED
    mds.line("blueRED", "whiteRED"),
    mds.line("whiteRED", "redRED"),
    mds.line("purpleRED", "redRED"),
    mds.line("purpleRED", "blueRED"),
    ## PURPLE
    mds.line("bluePURPLE", "redPURPLE"),
    mds.line("redPURPLE", "purplePURPLE"), 
    mds.line("purplePURPLE", "whitePURPLE"),
    mds.line("whitePURPLE", "bluePURPLE")
  ),
  
  incon = list(
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
    mds.line("blueWHITE", "bluePURPLE"),
    mds.line("blueRED", "bluePURPLE"),
    mds.line("blueWHITE", "blueBLUE"),
    mds.line("blueBLUE", "blueRED"),
    ## white
    mds.line("whiteWHITE", "whitePURPLE"),
    mds.line("whitePURPLE", "whiteRED"),
    mds.line("whiteBLUE", "whiteRED"),
    mds.line("whiteBLUE", "whiteWHITE"),
    ## red
    mds.line("redRED", "redBLUE"),
    mds.line("redBLUE", "redPURPLE"),
    mds.line("redWHITE", "redPURPLE"),
    mds.line("redRED", "redWHITE"),
    ## purple
    mds.line("purplePURPLE", "purpleRED"),
    mds.line("purpleRED", "purpleBLUE"),
    mds.line("purpleBLUE", "purpleWHITE"),
    mds.line("purpleWHITE", "purplePURPLE")
  )
  
)

## get grobs

grobs <- list(
  target       = plot.mds(mds$targt, .add.lines = mds.lines$targt, .size = rel(2), .fill = "transparent"),
  distractor   = plot.mds(mds$distr, .add.lines = mds.lines$distr, .size = rel(2), .fill = "transparent"),
  incongruency = plot.mds(mds$incon, .add.lines = mds.lines$incon, .size = rel(2), .fill = "transparent")
) %>%
  lapply(ggplotGrob)

## create layout

radius.rel <- 1.2  ## relative size of each MDS circle (used to build positions of circle centers)
radius.abs <- unit(figwidth/4.5, "cm")  ## absolute size of circle
a <- 0.2  ## scale factor of coordinates (controls distances among MDS circles)
b <- 2  ## scale factor of MDS embedding
h <- 2.54*figwidth

coordinates <- data.frame(
  x = seq(radius.rel, radius.rel * 3, by = radius.rel),
  y = c(sqrt(3) * radius.rel * 2, sqrt(3) * radius.rel, sqrt(3) * radius.rel * 2)
)
x.center <- 0.5
y.center <- 0.575
coordinates <- coordinates * a
coordinates <- sweep(scale(coordinates, scale = FALSE), 2, c(x.center, y.center), "+") %>% as.data.frame

## order everything
mds.order <- c("target", "incongruency", "distractor")
colors.profs <- colors.model[mds.order]
grobs <- grobs[mds.order]

## draw to file:

cairo_pdf(
  here("out", "group", "mds.pdf"),
  height = h, width = h
)

grid.circle(
  coordinates$x, coordinates$y, 
  r = radius.abs, gp = gpar(lwd = 3, fill = "transparent", col = colors.profs)
)

vps <- vector("list", length(grobs))
for (ii in seq_along(grobs)) {
  vp <- viewport(x = coordinates$x[ii], y = coordinates$y[ii], width = radius.abs * b, height = radius.abs * b)
  pushViewport(vp)
  grid.draw(grobs[[ii]])
  upViewport()
}

p.mds <- grid.grab(wrap.grobs = TRUE)

dev.off()

## draw for rmd:
rois.mds
grid.arrange(p.mds)
