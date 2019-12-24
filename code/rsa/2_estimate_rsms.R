## about ----
## creates model similarity matrices (for RSA) and writes them to .csv files.
## 
## mike freund, 2019-02-20
## adapted for new project directory 2019-12-24

## setup ----

library(here)
library(mikeutils)
library(magrittr)
library(dplyr)
library(abind)
library(data.table)
library(oro.nifti)

source(here("code", "_strings.R"))
source(here("code", "_funs.R"))
do.clusters <- TRUE
# do.held.out <- TRUE
do.mb <- c("mb4" = TRUE, mb8 = FALSE)  ## to skip prompts in 
do.read.atlas <- c("mmp" = TRUE, gordon = TRUE)
source(here("code", "_get_atlas.R"))
# source(here("r", "group-201902", "_get_misc_vars.R"))

# library(gifti)
# library(plot3D)

# source(here("..", "gen", "funs", "_funs.R"))
# source(here("..", "gen", "funs", "_get_dirs_remote.R"))
# source(here("r", "group-201902", "_read_sheets_remote.R"))
