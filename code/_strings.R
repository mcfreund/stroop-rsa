## paths: pointers to nil-bluearc
dir.nil.dmcc2.afni <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/AFNI_ANALYSIS"
# dir.freund.external <- "/data/nil-external/ccp/freund"
# dir.freund.scratch1 <- "/scratch1/ccp/freundm"
# dir.git <- file.path(dir.freund.scratch1, "R01")
# dir.atlas <- "/data/nil-external/ccp/freund/atlases"

## factor levels

bias.colors <- c("blue", "purple", "red", "white")
bias.words  <- toupper(bias.colors)
bias.items  <- mikeutils::combo.paste(bias.colors, bias.words, sep = "")
bias.items.incon <- bias.items[!bias.items %in% c("blueBLUE", "redRED", "whiteWHITE", "purplePURPLE")]
bias.items.con <- bias.items[!bias.items %in% bias.items.incon]
pc50.colors <- c("black", "green", "pink", "yellow")
pc50.words  <- toupper(pc50.colors)
pc50.items  <- mikeutils::combo.paste(pc50.colors, pc50.words, sep = "")
pc50.items.incon <- pc50.items[!pc50.items %in% c("blackBLACK", "greenGREEN", "pinkPINK", "yellowYELLOW")]
pc50.items.con <- pc50.items[!pc50.items %in% pc50.items.incon]