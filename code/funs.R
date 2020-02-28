
split.str.item <- function(col.j, prefix = "") {
  ## takes a single "item" vector, e.g. "blueBLUE", and decomposes
  ## it into color ("blue") word ("BLUE"), congruency ("C"), and label 
  ## ("C.BLUE", for plotting). 
  # prefix <- paste0(col.j, ".")
  col.j      <- as.character(col.j)
  color      <- gsub("[A-Z]", "", col.j)
  word       <- gsub("[a-z]", "", col.j)
  congruency <- ifelse(color == tolower(word), "C", "I")
  label      <- paste0(congruency, ".", word)
  cols       <- as.data.frame(cbind(color, word, congruency, label))
  colnames(cols) <- paste0(prefix, c("color", "word", "congruency", "label"))
  return(cols)
}


nmds <- function(arr, roi, ...) {
  
  D <- 1 - apply(rsarray[, , , roi], c(".row", ".col"), function(.) tanh((mean(atanh(.)))))
  D <- as.dist(D)
  metaMDS(D, trace = FALSE, ...)
  
}

split.str.item <- function(col.j, prefix = "") {
  ## takes a single "item" vector, e.g. "blueBLUE", and decomposes
  ## it into color ("blue") word ("BLUE"), congruency ("C"), and label 
  ## ("C.BLUE", for plotting). 
  # prefix <- paste0(col.j, ".")
  col.j      <- as.character(col.j)
  color      <- gsub("[A-Z]", "", col.j)
  word       <- gsub("[a-z]", "", col.j)
  congruency <- ifelse(color == tolower(word), "C", "I")
  label      <- paste0(congruency, ".", word)
  cols       <- as.data.frame(cbind(color, word, congruency, label))
  colnames(cols) <- paste0(prefix, c("color", "word", "congruency", "label"))
  return(cols)
}

plot.mds <- function(.df, .title = "", .add.lines = "nolines", .fill = "white", .size = rel(2.5), .margins = rep(0, 4)) {
  
  .df <- data.frame(scale(.df$points), item = row.names(.df$points), stringsAsFactors = FALSE)
  .df <- cbind(.df, split.str.item(.df$item))
  .df$word <- substr(.df$word, 1, 1)
  
  p <- .df %>% ggplot(aes(MDS1, MDS2))
  
  if (any(.add.lines == "nolines")) {
    
  } else { p <- p + .add.lines }
  
  p <- p +
    geom_label(
      aes(label = word, color = color), fill = .fill,
      fontface = "bold", label.padding =  unit(0, "lines"), label.size = 0, size = .size
    ) +
    scale_color_manual(values = c(blue = "blue", purple = "purple", red = "red", white = "grey50")) +
    theme_void(base_size = 6) +
    theme(
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      # title      = element_blank(),
      title      = element_text(face = "italic", size = rel(0.8)),
      panel.background = element_blank(),
      legend.position = "none",
      axis.line = element_blank(),
      plot.margin = unit(.margins, "cm")
    ) 
  # scale_x_continuous(expand = c(-0.5, 0.5)) +
  # scale_y_continuous(expand = c(-0.5, 0.5))
  
  if (.title != "") p <- p + labs(title = .title)
  
  p
  
}

mds.line <- function(b, e) {
  geom_segment(
    aes(x = MDS1[item == b], xend = MDS1[item == e], y = MDS2[item == b], yend = MDS2[item == e]), 
    size = 0.15, color = "grey70"
  )
}



cifti.parcellate <- function(
  path.wb = "C:/Program Files/workbench-windows64-v1.2.3/workbench/bin_windows64",
  dir.to.write = here::here("out", "runwise", "wb"),
  name.atlas = "schaefer400"
) {
  ## adapted from Jo Etzel, 2019-02-26
  ## see: http://mvpa.blogspot.com/2017/11/assigning-arbitrary-values-to-surface.html
  
  if (!dir.exists(dir.to.write)) stop(paste0("missing: ", dir.to.write))
  
  ## location of wb_command.exe, to call wb_command functions via system() within R:
  dir.origwd <- getwd()
  on.exit(setwd(dir.origwd))  ## return to original wd upon exit
  setwd(path.wb)
  
  ## file names and dirs
  name.atlas <- tolower(name.atlas)
  dir.s1200 <- file.path(dir.atlas, "surf/HCP_S1200_GroupAvg_v1") ## HCP S1200 Group Average Data Release
  
  if (name.atlas == "schaefer400") {
    
    fname.dlabel <- "Schaefer2018_400Parcels_7Networks_order"
    fname.dlabel.full <- file.path(
      dir.atlas, "Schaefer2018_Parcellations", "HCP", "fslr32k", "cifti", paste0(fname.dlabel, ".dlabel.nii")
    )
    
  } else if (name.atlas == "glasser") {
    
    fname.dlabel <- "Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR"
    fname.dlabel.full <- file.path(dir.s1200, paste0(fname.dlabel, ".dlabel.nii"))
    
  } else if (name.atlas == "gordon") {
    
    fname.dlabel <- "Parcels_LR"
    fname.dlabel.full <- file.path(
      dir.atlas, "gordon", "gordon_parcels", "Parcels", paste0(fname.dlabel, ".dlabel.nii")
    )
    
  } else {stop("did not find name.atlas")}
  
  fname.template <- file.path(dir.to.write, paste0(fname.dlabel, ".pscalar.nii"))
  
  stdout.template <- system(
    paste0(
      "wb_command -cifti-parcellate ", dir.s1200, "/S1200.thickness_MSMAll.32k_fs_LR.dscalar.nii ", 
      fname.dlabel.full, " COLUMN ",
      fname.template
    )
  )
  
  if (stdout.template != 0) stop("problem writing template")
  if (!file.exists(fname.template)) stop(paste0("missing: ", fname.template))
  
  return(paste0("wrote pscalar for ", name.atlas))
  
}

# cifti.parcellate()
# cifti.parcellate(name.atlas = "gordon")
# cifti.parcellate(name.atlas = "glasser", dir.to.write = here("out", "figs"))

cifti.convert <- function(
  path.wb = "C:/Program Files/workbench-windows64-v1.2.3/workbench/bin_windows64",
  dir.to.write = here::here("out", "runwise", "wb"),
  fname.overlay = "asdf",
  values.overlay = as.numeric(as.factor(parcellation$key)),
  name.atlas = "schaefer400",
  dir.template
) {
  ## adapted from Jo Etzel, 2019-02-26
  ## see: http://mvpa.blogspot.com/2017/11/assigning-arbitrary-values-to-surface.html
  
  if (!dir.exists(dir.to.write)) stop(paste0("missing: ", dir.to.write))
  
  ## location of wb_command.exe, to call wb_command functions via system() within R:
  dir.origwd <- getwd()
  on.exit(setwd(dir.origwd))  ## return to original wd upon exit
  setwd(path.wb)
  
  ## file names and dirs
  dir.s1200 <- file.path(dir.atlas, "surf/HCP_S1200_GroupAvg_v1") ## HCP S1200 Group Average Data Release
  fname.text <- file.path(dir.to.write, paste0(fname.overlay, ".txt"))
  fname.cifti <- file.path(dir.to.write, paste0(fname.overlay, ".pscalar.nii"))
  
  if (name.atlas == "schaefer400") {
    
    fname.template <- file.path(dir.template, "Schaefer2018_400Parcels_7Networks_order.pscalar.nii")
    
  } else if (name.atlas == "glasser") {
    
    fname.template <- file.path(
      dir.template, 
      "Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.pscalar.nii"
    )
    
  } else if (name.atlas == "gordon") {
    
    fname.template <- file.path(dir.template, "Parcels_LR.pscalar.nii")
    
  } else {stop("did not find name.atlas")}
  
  if (!file.exists(fname.template)) stop(paste0("missing: ", fname.template))
  
  ## the text file needs to be arranged with the right hemisphere parcels in the first N/2 rows (in order),   
  ## then the N/2 parcels for the left hemisphere.
  write.table(values.overlay, fname.text, col.names = FALSE, row.names = FALSE)
  
  ## create a CIFTI from the text file for viewing in Workbench  
  stdout.cifti <- system(
    paste0("wb_command -cifti-convert -from-text ", fname.text, " ", fname.template, " ", fname.cifti)
  )
  if (stdout.cifti != 0) stop("problem writing cifti")
  if (!file.exists(fname.cifti)) stop(paste("missing:", fname.cifti))
  
  ## remove text file
  was.removed <- file.remove(file.path(dir.to.write, paste0(fname.overlay, ".txt")))
  if (!was.removed) stop("trouble removing text file")
  
}
