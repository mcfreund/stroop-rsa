library(ggplot2)
library(viridis)
library(colorspace)
library(grid)
library(gridExtra)
library(cowplot)
library(reshape2)
library(dplyr)
library(magrittr)
library(here)
library(MASS)
library(mikeutils)


source(here::here("code/packages.R"))

rsarray <- read_simil_mats(subjs = c(subjs.analysis, subjs.validation), tfrm = "residual-line")  ## for MDS


theme_set(theme_bw(base_size = 12))
colors.model <- c(incon = "#d95f02", target = "#1b9e77", distractor = "#7570b3")
text.size <- 2.75



## create models ----

colors <- c("blue", "red", "purple", "white")
# colors <- c("blue", "red")
words <- toupper(colors)
labels <- expand.grid(lapply(list(words, colors), sort))
names(labels) <- c("word", "color")
is.level <- model.matrix(~ . + 0, labels, contrasts.arg = lapply(labels, contrasts, contrasts = FALSE)) == 1

empty <- diag(nrow(labels))
colnames(empty) <- apply(labels, 1, paste, collapse = "_")
rownames(empty) <- apply(labels, 1, paste, collapse = "_")

hue <- empty  ## hue
wor <- empty  ## word

for (col.i in grep("color", colnames(is.level))) hue[is.level[, col.i], is.level[, col.i]] <- 1
for (col.i in grep("word", colnames(is.level))) wor[is.level[, col.i], is.level[, col.i]] <- 1

is.inc <- labels$color != tolower(labels$word)
inc <- empty
diag(inc) <- 0
inc[is.inc, is.inc] <- 1




## model ----



symmat4ggplot <- function(R, var.names = c("v1", "v2"), val.name = "value") {
  
  ## make factors for row and column labels
  dn <- dimnames(R)
  if (is.null(dn)) {
    dn <- setNames(list(paste0("cell_", 1:nrow(R)), paste0("cell_", 1:ncol(R))), var.names)
  } else {
    names(dn) <- var.names  
  }
  
  labels <- expand.grid(dn, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = TRUE)
  labels[[2]] <- factor(labels[[2]], levels = rev(levels(labels[[2]])))
  
  r <- c(R)
  
  cbind(labels, setNames(as.data.frame(c(R)), val.name))
  
}

labels$color <- as.character(labels$color)
labels$color[labels$color == "white"] <- "grey50"

plotmat <- function(x, title, axis.text.size = 1.5, strip.text.size = 3) {
  
  
  # x %>%
  #   
  # symmat4ggplot %>%
  #   
  #   ggplot(aes(v1, v2, fill = value)) +
  #   # geom_tile(aes(fill = value), color = "transparent") +
  #   geom_raster(aes(color = value)) +
  #   
  #   scale_fill_continuous_diverging(palette = "Blue-Red", rev = TRUE) +
  #   scale_x_discrete(labels = as.character(labels$word)) +
  #   scale_y_discrete(labels = rev(as.character(labels$word))) +
  #   
  #   theme(
  #     legend.position = "none",
  #     panel.grid = element_blank(),
  #     axis.text.x = element_text(color = as.character(labels$color), face = "bold", size = rel(axis.text.size), angle = 90, hjust = 1, vjust = 0.5),
  #     axis.text.y = element_text(color = rev(as.character(labels$color)), face = "bold", size = rel(axis.text.size)),
  #     axis.title = element_blank(),
  #     plot.title = element_text(size = rel(5))
  #   ) +
  #   
  #   labs(title = title)
  x %>%
    
    symmat4ggplot %>%
    separate(v1, c("word1", "color1", "rule1", "modality1")) %>%
    separate(v2, c("word2", "color2", "rule2", "modality2")) %>%
    
    ggplot(aes(interaction(word1, color1), rev(interaction(word2, color2)), fill = value)) +
    geom_raster() +
    
    scale_fill_continuous_diverging(palette = "Blue-Red", rev = TRUE) +
    scale_x_discrete(labels = as.character(labels$word)) +
    scale_y_discrete(labels = rev(as.character(labels$word))) +
    
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text.x = element_text(color = as.character(labels$color), face = "bold", size = rel(axis.text.size), angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(color = rev(as.character(labels$color)), face = "bold", size = rel(axis.text.size)),
      axis.title = element_blank(),
      plot.title = element_text(size = rel(5)),
      panel.spacing = unit(0, "lines"), 
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(size = rel(strip.text.size), face = "bold")
    ) +
    
    labs(title = title)
  
  
}


p.target <- hue %>% plotmat("target")
p.distractor <- wor %>% plotmat("distractor")
p.incongruency <- inc %>% plotmat("incongruency")


m.smmouth <- apply(atanh(rsarray[, , , "smmouth"]), c(".row", ".col"), mean) %>% tanh
# m.smmouth <- apply(rsarray[, , , "smmouth"], c(".row", ".col"), mean)

# dimnames(m.smmouth) <- dimnames(hue)


d <- m.smmouth %>% rank.mat %>% symmat4ggplot
# d <- m.smmouth %>% symmat4ggplot
d1 <- split.str.item(d$v1)
names(d1) <- paste0(names(d1), "1")
d2 <- split.str.item(d$v2)
names(d2) <- paste0(names(d2), "2")

axis.text.size = 1.5
strip.text.size = 3

smmouth <- d %>% 
  bind_cols(d1, d2) %>%
  ggplot(aes(interaction(word1, color1), rev(interaction(word2, color2)), fill = value)) +
  geom_raster() +
  
  scale_fill_continuous_diverging(palette = "Blue-Red", rev = TRUE) +
  scale_x_discrete(labels = as.character(labels$word)) +
  scale_y_discrete(labels = rev(as.character(labels$word))) +
  
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(color = as.character(labels$color), face = "bold", size = rel(1.9), angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(color = rev(as.character(labels$color)), face = "bold", size = rel(1.9)),
    axis.title = element_blank(),
    plot.title = element_text(size = rel(5)),
    panel.spacing = unit(0, "lines"), 
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = rel(strip.text.size), face = "bold")
  ) +
  
  labs(title = "smmouth")
  


ggsave(here("out/method_fig/target.pdf"), p.target, device = "pdf", width = 5, height = 5.5, units = "in")
ggsave(here("out/method_fig/distractor.pdf"), p.distractor, device = "pdf", width = 5, height = 5.5, units = "in")
ggsave(here("out/method_fig/incongruency.pdf"), p.incongruency, device = "pdf", width = 5, height = 5.5, units = "in")
ggsave(here("out/method_fig/observed.pdf"), smmouth, device = "pdf", width = 5, height = 5.5, units = "in")
ggsave(here("out/method_fig/target_withledge.pdf"), p.target + theme(legend.position = "top"), device = "pdf", width = 5, height = 5.5, units = "in")


