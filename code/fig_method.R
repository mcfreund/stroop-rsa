library(ggplot2)
library(here)
library(magrittr)
library(ellipse)

theme_set(theme_bw(base_size = 12))
colors.model <- c(incon = "#d95f02", target = "#1b9e77", distractor = "#7570b3")
text.size <- 7

shade <- function(x, ...) {
  dnorm(x, ...)
}

group.mfc <- data.frame(x = c(-3, 4)) %>%
  
  ggplot(aes(x = x)) +
  
  stat_function(fun = dnorm, color = colors.model["target"], size = 1.5) +
  stat_function(fun = dnorm, color = colors.model["incon"], size = 1.5, args = list(1)) +
  
  stat_function(fun = shade, geom = "area", fill = colors.model["target"], alpha = 0.8) +
  stat_function(fun = shade, geom = "area", fill = colors.model["incon"], alpha = 0.8, args = list(mean = 1)) +
  
  scale_x_continuous(labels = 0, breaks = 0) +
  
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),
    axis.line.x = element_line(), axis.title.y = element_blank(), axis.text = element_text(size = rel(text.size), color = "black"), 
    axis.title = element_text(size = rel(text.size))
    ) +
  
  labs(x = expression(beta))

group.dlpfc <- data.frame(x = c(-3, 4)) %>%
  
  ggplot(aes(x = x)) +
  
  stat_function(fun = dnorm, color = colors.model["incon"], size = 1.5) +
  stat_function(fun = dnorm, color = colors.model["target"], size = 1.5, args = list(1)) +
  
  stat_function(fun = shade, geom = "area", fill = colors.model["incon"], alpha = 0.8) +
  stat_function(fun = shade, geom = "area", fill = colors.model["target"], alpha = 0.8, args = list(mean = 1)) +
  
  scale_x_continuous(labels = 0, breaks = 0) +
  
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),
    axis.line.x = element_line(), axis.title.y = element_blank(), axis.text = element_text(size = rel(text.size), color = "black"), 
    axis.title = element_text(size = rel(text.size))
  ) +
  
  labs(x = expression(beta))


# https://stackoverflow.com/questions/36221596/plot-multivariate-gaussian-contours-with-ggplot2
m <- c(0, 0)
sigma <- matrix(c(1,.5,.5,1), nrow=2)
# sigma_neg <- matrix(c(1,-.5,-.5,1), nrow=2)
alpha_levels <- seq(0.5, 0.95, by = 0.05) ## or whatever you want
names(alpha_levels) <- alpha_levels ## to get id column in result
contour_data <- plyr::ldply(alpha_levels,ellipse,x=sigma,
                      scale=c(1,1),  ## needed for positional matching
                      centre=m)

indif.mfc <- contour_data %>%
  ggplot(aes(x, y, group = .id)) + geom_polygon(aes(fill = as.numeric(.id)), alpha = 0.1, fill = colors.model["incon"]) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),
    axis.line = element_line(), axis.title = element_text(size = rel(text.size))
  ) +
  labs(x = expression(beta["incongruency"]), y = "Stroop\neffect")
  
  
indif.dlpfc <- contour_data %>%
  ggplot(aes(x, -y, group = .id)) + geom_polygon(aes(fill = as.numeric(.id)), alpha = 0.1, fill = colors.model["target"]) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),
    axis.line = element_line(), axis.title = element_text(size = rel(text.size))
  ) +
  labs(x = expression(beta["target"]), y = "Stroop\neffect")



ggsave(here("out/figs/ms_v1_2020-03/method/group_mfc.pdf"), group.mfc, device = "pdf", width = 10, height = 6, units = "in")
ggsave(here("out/figs/ms_v1_2020-03/method/group_dlpfc.pdf"), group.dlpfc, device = "pdf", width = 10, height = 6, units = "in")
ggsave(here("out/figs/ms_v1_2020-03/method/indif_mfc.pdf"), indif.mfc, device = "pdf", width = 10, height = 7, units = "in")
ggsave(here("out/figs/ms_v1_2020-03/method/indif_dlpfc.pdf"), indif.dlpfc, device = "pdf", width = 10, height = 7, units = "in")

