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
do.read.atlas <- c("mmp" = TRUE, gordon = TRUE)
source(here("code", "_get_atlas.R"))

