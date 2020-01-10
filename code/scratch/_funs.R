`%>%` <- magrittr::`%>%`

rank.to.rt <- function(df, subj, n, from.top = TRUE) {
  rt <- df[df$subj == subj, "rt"]
  if (from.top) n <- length(rt) - n
  rt[rank(rt) == n]
}

split.str.cell <- function(df, colname = "cell", sep = "_", new.names = c("a", "b")) {
  ## takes a df and "cell" vector, e.g., "blueBLUE_redBLUE", splits it at sep into
  ## two cols, and returns cols bound to df
  cols <- as.character(df[, colname])
  cols <- dplyr::bind_cols(strsplit(cols, split = sep))
  cols <- as.data.frame(t(cols))
  names(cols) <- new.names
  dplyr::bind_cols(df[, -grep(colname, names(df))], cols)
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

lm.beta <- function(MOD) {
  ## see ?QuantPsyc::lm.beta()
  
  b  <- summary(MOD)$coef[-1, 1]
  sx <- sapply(MOD$model[-1], sd)
  sy <- sapply(MOD$model[1], sd)
  
  b * sx / sy
  
}
