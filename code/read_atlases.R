## about ----
##
## reads atlases and atlas info.
## Multi-modal parcellation: Glasser (2018)
## Gordon: Gordon (2014)
## 
## depends on successful completion of  _write_atlases.R script.

## atlases

atlas <- setNames(vector("list", 2), c("mmp", "gordon"))

atlas$mmp <- oro.nifti::readNIfTI(here::here("out", "atlases", "mmp.nii.gz"), reorient = FALSE)

atlas$gordon <- oro.nifti::readNIfTI(here::here("out", "atlases", "gordon.nii.gz"), reorient = FALSE)


## keys

atlas.key <- setNames(vector("list", 2),  c("mmp", "gordon"))

atlas.key$mmp <- read.csv(here::here("out", "atlases", "mmp.csv"), stringsAsFactors = FALSE)
atlas.key$gordon <- read.csv(here::here("out", "atlases", "gordon.csv"), stringsAsFactors = FALSE)

atlas.key$mmp$X <- NULL
atlas.key$gordon$X <- NULL


## rearrange for workbench underlay

inds.left <- atlas.key$mmp$roi %>% grep("_L$", .)
inds.right <- atlas.key$mmp$roi %>% grep("_R$", .)
atlas.key$mmp$roi <- atlas.key$mmp$roi[c(inds.right, inds.left)]  ## correctly order mmp atlas

## for workbench underlay (get filename, write if non-existent)

fname.pscalar.mmp <- here(
  "out", "wb", 
  "Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.pscalar.nii"
)

if (!file.exists(fname.pscalar.mmp)) cifti.parcellate(name.atlas = "glasser", dir.to.write = here("out", "wb"))

## add network and MD info ----

# coleanticevic <- RCurl::getURL(
#   "https://raw.githubusercontent.com/ColeLab/ColeAnticevicNetPartition/master/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR_LabelKey.txt"
# )
# coleanticevic <- fread(text = coleanticevic)


coleanticevic <- read.table(
  url("https://raw.githubusercontent.com/ColeLab/ColeAnticevicNetPartition/master/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR_LabelKey.txt"),
  header = TRUE
  )
coleanticevic <- as.data.table(coleanticevic)
coleanticevic <- coleanticevic[!is.na(GLASSERLABELNAME), c("NETWORK", "GLASSERLABELNAME")]
coleanticevic$GLASSERLABELNAME <- gsub("(^.)_(.*)_ROI", "\\2_\\1", coleanticevic$GLASSERLABELNAME)
coleanticevic %<>% rename(roi = GLASSERLABELNAME, network = NETWORK)
atlas.key$mmp <- full_join(atlas.key$mmp, coleanticevic, by = "roi")

## MD assignments (Assem et al)

md <- list(
  core       = c("p9-46v", "a9-46v", "i6-8", "AVI", "8C", "IFJp", "IP2", "IP1", "PFm", "8BM", "SCEF"),
  extended   = c(
    "a9-46v", "p10p", "a10p", "11l", "a47r", "p47r", "FOP5", "AVI", "p9-46v", "8C", "IFJp", "6r", "s6-8", "i6-8",
    "SCEF", "8BM", "a32pr", "d32",
    "TE1m", "TE1p",
    "AIP", "IP2", "LIPd", "MIP", "IP1", "PGs", "PFm", "POS2"
    # "p9-46v", "a9-46v", "i6-8", "AVi", "8C", "IFJp", "IP2", "IP1", "PFm", "8BM",
    # "TE1m", "TE1p", "PGs", "PFm", "AIP", "MIP", "LIPd", "IP1", "IP2", "s6-8", 
    # "i6-8", "a9-46v", "FOP5", "AVI", "11l", "a10p", "p10p", "a47r", "p47r"
  )
)
# length(setdiff(md$extended, md$core))
# length(md$core)
# md$extended[!md$extended %in% gsub("_R|_L", "", atlas.key$mmp$roi)]

atlas.key$mmp$md <- ifelse(
  atlas.key$mmp$roi %in% combo_paste(md$core, c("R", "L")), 
  "core",
  ifelse(
    atlas.key$mmp$roi %in% combo_paste(md$extended, c("R", "L")), 
    "extended", "none"
  )
)
