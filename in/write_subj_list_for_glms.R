source(here::here("code", "strings.R"))

subjs <- c(subjs.analysis, subjs.validation)


write.table(subjs, here::here("in", "subjects.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)