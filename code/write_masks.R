## about ----
## 
## writes masks, .nii files, that define ROIs
## 
## warning: this script is spagett

library(here)
library(mikeutils)
library(magrittr)
library(dplyr)
library(data.table)
library(oro.nifti)
source(here("code", "strings.R"))
source(here("code", "read_atlases.R"))
source(here("code", "funs.R"))

write.masks <- function(l, atlas.name) {
  
  .atlas <- atlas[[atlas.name]]
  .atlas.key <- atlas.key[[atlas.name]]
  
  for (ii in seq_along(l)) {
    num.rois <- .atlas.key$num.roi[.atlas.key$roi %in% l[[ii]]]
    mask.i <- array(as.numeric(.atlas@.Data %in% num.rois), dim = c(75, 90, 75))
    
    writeNIfTI(as.nifti(mask.i), here::here("out", "masks", names(l)[ii]))
    
  }
  
}


## define superparcels ----

## from MMP:

superparcels <- list(
  ## lateral posterior parietal cortex (intraparietal sulcus):
  ## parcels tiling intraparietal sulcus (from IP0 to S1 [area 2]) and IPL (medial aspect) / SPL (lateral aspect)
  lppc     = c("IP0", "IPS1", "IP1", "MIP", "7PL", "7AL", "7PC", "VIP", "LIPv", "LIPd", "AIP", "IP2"),
  ## dorsomedial fc:
  ##  SMA/pre-sma: SCEF, 8BM
  ##  dACC: a32pr, p32pr, d32 (d32 is default mode in Cole-Anticevic)
  dmfc      = c("SCEF", "8BM",  "p32pr", "a32pr"),
  dmfc.alt  = c("SCEF", "8BM",  "p32pr", "a32pr", "d32"),
  ## dlpfc:
  ## parcels tiling middle frontal gyrus
  ## 46 is default mode in cole-anticevic
  ## a9-46v, 9-46d are rostral, can be considered fronto-polar cortex
  dlpfc     = c("p9-46v", "i6-8", "8Av", "8C"), ## along middle frontal gyrus
  dlpfc.alt = c("p9-46v", "i6-8", "8Av", "8C", "46", "a9-46v", "9-46d"),
  dpm       = c("6a", "FEF", "6ma"), ## dorsal premotor
  vpm       = c("6v", "6v", "PEF", "6r", "55b", "IFJp"),  ## ventral premotor
  vvis      = c("FFC", "VVC", "V8", "VMV3"),
  ins       = c("FOP4", "FOP5", "FOP3", "AVI", "MI"),
  ifc       = c("IFJa", "44", "45", "IFSa", "p47r", "IFSp", "a47r", "47l"),
  ofc       = c("47s", "47m", "11l", "13l", "OFC", "10pp", "10v", "10r", "s32", "p32", "a24", "25", "pOFC"),
  fpc       = c("46", "9-46d", "a9-46v", "a10p", "p10p", "a47r", "9a", "9p", "10d")
)

superparcels.r <- lapply(superparcels, function(x) paste0(x, "_R"))
superparcels.l <- lapply(superparcels, function(x) paste0(x, "_L"))
superparcels <- c(
  setNames(superparcels.l, paste0(names(superparcels.r), "_L")),
  setNames(superparcels.r, paste0(names(superparcels.r), "_R")),
  list(
    aud   = combo_paste(c("LBelt", "A1", "MBelt", "PBelt"), c("L", "R")),
    evis  = combo_paste(c("V1", "V2", "V3"), c("L", "R")),
    V1    = paste0("V1", c("_L", "_R"))  ## for prelim. validation analysis (group-level)
  )
)

write.masks(superparcels, "mmp")


## from gordon:

num.rois <- atlas.key$gordon$num.roi[atlas.key$gordon$roi %in% c("smmouth_L", "smmouth_R")]
smmouth <- array(as.numeric(atlas$gordon@.Data %in% num.rois), dim = c(75, 90, 75))
writeNIfTI(as.nifti(smmouth), here::here("out", "masks", "smmouth"))


## write parcel IDs ----

saveRDS(superparcels, here::here("out", "masks", "ids_superparcels.RDS"))


## create workbench files ----



## all superparcels

atlas.key$mmp$id <- "none"
superparcels4wb <- superparcels[names(superparcels)[-grep("V1|\\.alt", names(superparcels))]]
for (superparcel.i in seq_along(superparcels4wb)) {
  atlas.key$mmp$id[atlas.key$mmp$roi %in% superparcels4wb[[superparcel.i]]] <- names(superparcels4wb)[superparcel.i]
}
atlas.key$mmp$id <- relevel(as.factor(atlas.key$mmp$id), "none")
overlay <- atlas.key$mmp %>% select(roi, id)
overlay$id <- as.numeric(overlay$id) - 1
overlay$roi.num <- match(overlay$roi, atlas.key$mmp$roi)
overlay %<>% arrange(roi.num)
cifti.convert(
  fname.overlay  = "superparcels_all",
  values.overlay = overlay$id,
  dir.template   = here("out", "wb"),
  name.atlas     = "glasser",
  dir.to.write   = here("out", "masks")
)


## only ROIs
overlay.rois <- overlay
inds.rois <- which(levels(atlas.key$mmp$id) %in% c("lppc_L", "lppc_R", "dlpfc_L", "dlpfc_R", "dmfc_L", "dmfc_R"))-1
overlay.rois$id[!overlay.rois$id %in% inds.rois] <- 0
cifti.convert(
  fname.overlay  = "superparcels_rois",
  values.overlay = overlay.rois$id,
  dir.template   = here("out", "wb"),
  name.atlas     = "glasser",
  dir.to.write   = here("out", "masks")
)


## for gordon smmouth

fname.overlay <- "gordon"
dir.to.write <- here("out", "masks")
dir.origwd <- getwd()
setwd("C:/Program Files/workbench-windows64-v1.2.3/workbench/bin_windows64")
## for reading and writing files
s1200.path <- "C:/local/atlases/surf/HCP_S1200_GroupAvg_v1/" ## HCP S1200 Group Average Data Release  
fname.template <- "gordon-template.pscalar.nii"  ## filename for the template we'll make

## make a template pscalar.nii from the GORDON atlas  
system(
  paste0(
    "wb_command -cifti-parcellate ", s1200.path, "S1200.thickness_MSMAll.32k_fs_LR.dscalar.nii ", 
    "C:/local/atlases/gordon/gordon_parcels/Parcels/Parcels_LR.dlabel.nii COLUMN ",  
    file.path(dir.to.write, fname.template)
  )
)
if (!file.exists(file.path(dir.to.write, fname.template))) stop(paste0("missing: ", dir.to.write, fname.template))

## make a text version of the template, which has 360 rows, but the values in each row aren't the parcel numbers.  
system(
  paste0(
    "wb_command -cifti-convert -to-text ", file.path(dir.to.write, fname.template), " ", file.path(dir.to.write, fname.overlay), "_text.txt"
  )
)


## the text file needs to be arranged with the 180 right hemisphere parcels in the first 180 rows (in order),   
## then the 180 parcels for the left hemisphere.  

gordon.parcels <- read.csv("C:/local/atlases/gordon/Parcels.csv", stringsAsFactors = FALSE)
# gordon.parcels <- arrange(gordon.parcels, desc(Hem), ParcelID)
# gordon.parcels$comm.num <- as.numeric(as.factor(gordon.parcels$Community))
# write.table(gordon.parcels$comm.num, file.path(dir.to.write, paste0(fname.overlay, ".txt")), col.names = FALSE, row.names = FALSE)
write.table(
  as.numeric(gordon.parcels$Community == "SMmouth"),
  file.path(dir.to.write, paste0(fname.overlay, ".txt")), col.names = FALSE, row.names = FALSE
)

## create a CIFTI from the text file for viewing in Workbench  
system(
  paste0(
    "wb_command -cifti-convert -from-text ", file.path(dir.to.write, paste0(fname.overlay, ".txt ")), 
    file.path(dir.to.write, fname.template), " ", file.path(dir.to.write, paste0(fname.overlay, ".pscalar.nii"))
  )
)
if (!file.exists(file.path(dir.to.write, paste0(fname.overlay, ".pscalar.nii")))) {
  stop(paste("missing:", file.path(dir.to.write, paste0(fname.overlay, ".pscalar.nii"))))
}
setwd(dir.origwd)




## make table for masks ----

## get network info for each MMP parcel:

coleanticevic <- RCurl::getURL(
  "https://raw.githubusercontent.com/ColeLab/ColeAnticevicNetPartition/master/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR_LabelKey.txt"
)
coleanticevic <- fread(text = coleanticevic)
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



## LPPC
lppc <- atlas.key$mmp %>% filter(roi %in% c(superparcels$lppc_L, superparcels$lppc_R))

## DLPFC
dlpfc <- atlas.key$mmp %>%  filter(roi %in% c(superparcels$dlpfc_L, superparcels$dlpfc_R))

## DMFC
dmfc <- atlas.key$mmp %>% filter(roi %in% c(superparcels$dmfc_L, superparcels$dmfc_R))

## DLPFC---alt
dlpfc.alt <- atlas.key$mmp %>%  filter(roi %in% c(superparcels$dlpfc.alt_L, superparcels$dlpfc.alt_R))

## DMFC---alt
dmfc.alt <- atlas.key$mmp %>% filter(roi %in% c(superparcels$dmfc.alt_L, superparcels$dmfc.alt_R))

lppc
dmfc
dlpfc.alt
dmfc.alt

## superparcel table


superparcel.info <- data.frame(region = names(superparcels4wb), nvox = NA)
for (n.superparcel.i in seq_along(superparcels4wb)) {
  # n.superparcel.i = 1
  
  superparcel.i <- names(superparcels4wb)[n.superparcel.i]
  is.in.superparcel <- atlas.key$mmp$roi %in% superparcels4wb[[superparcel.i]]
  inds.superparcel <- atlas.key$mmp$num.roi[is.in.superparcel]
  nvox <- sum(atlas$mmp %in% inds.superparcel)
  
  superparcel.info[n.superparcel.i, "nvox"] <- nvox
  
}

superparcel.info$nparc <- sapply(superparcels4wb, length)
superparcel.info$hemi <- substr(
  superparcel.info$region,
  sapply(as.character(superparcel.info$region), nchar),
  sapply(as.character(superparcel.info$region), nchar)
)

nvox.smmouth <- sum(atlas$gordon %in% atlas.key$gordon$num.roi[atlas.key$gordon$roi %in% c("smmouth_L", "smmouth_R")])
superparcel.info %<>% rbind(
  data.frame(region = "smmouth", nvox = nvox.smmouth, nparc = NA, hemi = "bil.")
)


parcels <- lapply(superparcels4wb[-grep("_R$", names(superparcels4wb))], function(x) gsub("_L|_R", "", x) %>% unique)
parcels <- sapply(parcels, function(x) paste0(x, collapse = ", "))

tab.sp <- superparcel.info %>% filter(hemi != "R", region != "V1")
# tab.sp$hemi[tab.sp$hemi != "L"] <- "bil."
# tab.sp$hemi[tab.sp$hemi == "L"] <- "L, R"
tab.sp$parcels <- c(parcels, "")


tab.sp$abbreviation <- NA
tab.sp$abbreviation[tab.sp$region == "lppc_L"] <- "LPPC"
tab.sp$abbreviation[tab.sp$region == "dmfc_L"] <- "DMFC"
tab.sp$abbreviation[tab.sp$region == "dlpfc_L"] <- "DLPFC"
tab.sp$abbreviation[tab.sp$region == "dpm_L"] <- "DPM"
tab.sp$abbreviation[tab.sp$region == "vpm_L"] <- "VPM"
tab.sp$abbreviation[tab.sp$region == "vvis_L"] <- "VVis"
tab.sp$abbreviation[tab.sp$region == "ins_L"] <- "AIns"
tab.sp$abbreviation[tab.sp$region == "ifc_L"] <- "IFC"
tab.sp$abbreviation[tab.sp$region == "ofc_L"] <- "OFC"
tab.sp$abbreviation[tab.sp$region == "fpc_L"] <- "FPC"
tab.sp$abbreviation[tab.sp$region == "evis"] <- "V1--V3"
tab.sp$abbreviation[tab.sp$region == "smmouth"] <- "VSM"
tab.sp$abbreviation[tab.sp$region == "aud"] <- "Aud"

tab.sp$nvox[tab.sp$hemi == "L"] <- paste0(
  tab.sp %>% filter(hemi == "L") %>% pull(nvox), ", ",
  superparcel.info %>% filter(hemi == "R", region != "V1") %>% pull(nvox)
  )

tab.sp %<>% select(superparcel = abbreviation, "constituent MMP parcels" = parcels, "# voxels (L, R)" = nvox)
fwrite(tab.sp, here("out", "masks", "superparcel_all.csv"))
