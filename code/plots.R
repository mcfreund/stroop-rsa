colors.model <- c(incongruency = "#d95f02", target = "#1b9e77", distractor = "#7570b3")
colors.region <- setNames(viridis(3), c("DLPFC", "LPPC", "DMFC"))
params.interest <- names(colors.model)

theme_set(theme_bw(base_size = 8))
figwidth <- 4.2  ## cm
axis.text.size <- rel(1)
axis.title.size <- rel(1)
axis.line.size <- rel(1)
label.size <- rel(3)
p.value.size <- rel(2)
p.line.size <- rel(0.5)
geom.line.size <- rel(1)
geom.point.size <- rel(2)
