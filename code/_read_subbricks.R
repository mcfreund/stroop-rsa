## get afni images:
## TODO:
## -> REMOVE LOOP call for single function
## gifti functionality
## runwise functionality

read.subbricks <- function(
  image.full, 
  afni.path, 
  .label = "asdf", 
  .suffix = "+A", 
  runwise = FALSE, 
  wsl = FALSE
  ) {
  
  ## input validation
  if (!length(.suffix) %in% c(1, length(.label))) stop(".suffix wrong length (not 1 or length(.label)")
  if (!file.exists(fname)) stop("file nonexistant! ", paste0(fname))
  path.3dinfo <- file.path(afni.path, "3dinfo")
  if (!dir.exists(path.3dinfo)) stop(paste0(path.3dinfo), " does not exist")

  ## read image and create objects for looping
  image.full <- readNIfTI(fname, reorient = FALSE)  ## dims of image.run.full [i, j, k, ???, regressor]
  strings <- cbind(label = .labels, suffix = .suffixes, name = paste0(.labels, .suffixes))
  n.subbricks <- nrow(strings)
  subbrick.indices <- rep(NA, n.subbricks)

  for (subbrick.i in seq_len(n.subbricks)) {
    
    image.name <- strings[subbrick.i, "name"]
    
    if (wsl) {
      stdout.3dinfo <- system2(
        "wsl",
        args = paste("/home/mcf/abin/3dinfo", "-label2index", image.name, fname),
        stdout = TRUE
      )
    } else {
      stdout.3dinfo <- system2(
        "/usr/local/pkg/linux_openmp_64/3dinfo",
        args = paste("-label2index", image.name, fname),
        stdout = TRUE
      )
    }

    ## for error checking:
    has.error <- grepl("error", stdout.3dinfo, ignore.case = TRUE)
    if (any(has.error)) stop("error loading subbrick nums: ", paste0(fname))
    
    ## to remove function call that (may) be included in output (when is.local.session):
    stdout.3dinfo <- stdout.3dinfo[!grepl("3dinfo: AFNI version", stdout.3dinfo)]
    
    subbrick.nums[subbrick.i] <- as.numeric(stdout.3dinfo)
    
  }
  
  if (any(is.na(subbrick.nums))) stop("brick nums equal zero! ", paste0(fname))

  ## put subbricks into one array:
  image.stats <- image.full[, , , 1, brick.nums + 1]  ## +1 for 0 -> 1-based index
  dimnames(image.stats) <- list(i = NULL, j = NULL, k = NULL, reg = strings[, "name"])
  
  image.stats

}

