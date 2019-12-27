## get afni images:
## TODO:
## runwise functionality (?)

read.subbrick <- function(
  fname, 
  afni.path = "/usr/local/pkg/linux_openmp_64/", 
  .label = "asdf", 
  .suffix = "+A", 
  wsl = FALSE
  ) {
  
  if (wsl) {
    
    stdout.3dinfo <- system2(
      "wsl",
      args = paste(afni.path, "/3dinfo", "-label2index", fname, .label, .suffix),
      stdout = TRUE
    )
    
  } else {
    
    stdout.3dinfo <- system2(
      file.path(afni.path, "3dinfo"),
      args = paste("-label2index", fname),
      stdout = TRUE
    )
    
  }
  
  ## for error checking:
  
  has.error <- grepl("error", stdout.3dinfo, ignore.case = TRUE)
  if (any(has.error)) stop("error loading brick nums: ", paste0(fname))
  
  ## to remove function call that is included in output (when is.local.session):
  
  stdout.3dinfo <- stdout.3dinfo[!grepl("3dinfo: AFNI version", stdout.3dinfo)]
  subbrick.num <- as.numeric(stdout.3dinfo)
  if (is.na(subbrick.num)) stop("subbrick num is NA: ", paste0(fname, " ", .label, .suffix))
  
  
}

