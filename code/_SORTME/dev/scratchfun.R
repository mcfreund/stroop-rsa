
# colorspace::sequential_hcl()
# range(mask.tdic.congr$L)
# 
# overlay <- mask.tdic.incon
# underlay <- hcp

# 
# plot.surface <- function(gii, ttl, pos.lims=c(NA,NA), neg.lims=c(NA,NA), print.scaling.label=TRUE, which.surface="HCP") {
#   # gii <- tmp; ttl <- paste(lbls[i], "L", sess.ids[ssid]); pos.lims <- c(2,6); neg.lims <- c(-6,-2);
#   full.min <- min(gii$values, na.rm=TRUE);  # store original minimum and maximum for later text label
#   full.max <- max(gii$values, na.rm=TRUE);
#   sc.text <- "";   # start scaling label with a blank string
# 
#   # first fix the values for proper plotting
#   if (!is.na(pos.lims[1])) {    # showing positive values, so fix too-big values (so don't disappear when plotting)
#     inds <- which(gii$values > pos.lims[2]);
#     if (length(inds) > 0) { gii$values[inds] <- pos.lims[2]; }
#   }
#   if (!is.na(neg.lims[1])) {    # showing negative values, so fix too-small values
#     inds <- which(gii$values < neg.lims[1]);
#     if (length(inds) > 0) { gii$values[inds] <- neg.lims[1]; }
#   }
#   if (!is.na(pos.lims[1]) & !is.na(neg.lims[1])) {    # both positive and negtive, so fix middle (zero-ish) values
#     inds <- which(gii$values > neg.lims[2] & gii$values < pos.lims[1]);   # NA out middle values so don't get plotted
#     if (length(inds) > 0) { gii$values[inds] <- NA; }
#     c.lims <- c(neg.lims[1], pos.lims[2]);   # set color scaling to smallest negative to biggest positive
#     cols <- c(rev(cols.cool), cols.warm);   # both hot and cool colors
#     sc.text <- paste("colors:", neg.lims[1], "to", neg.lims[2], "&", pos.lims[1], "to", pos.lims[2]);
#   }
#   if (!is.na(pos.lims[1]) & is.na(neg.lims[1])) {   # only positive values
#     c.lims <- pos.lims;
#     cols <- cols.warm;
#     sc.text <- paste("colors:", pos.lims[1], "to", pos.lims[2]);
#   }
#   if (is.na(pos.lims[1]) & !is.na(neg.lims[1])) {   # only negative values
#     c.lims <- neg.lims;
#     cols <- rev(cols.cool);
#     sc.text <- paste("colors:", neg.lims[1], "to", neg.lims[2]);
#   }
# 
#   # slightly different z-axis scaling and title positioning seems to work better for HCP and fsaverage5 (fMRIPrep) surface anatomies
#   if (which.surface == "HCP") { z.lim <- c(-60,90); ttl.line <- -2.5; }
#   if (which.surface == "fsaverage5") { z.lim <- c(-130,90); ttl.line <- -1; }
# 
#   # actually plot
#   triangle3D(tri=gii$coords, theta=90, phi=0, ltheta=90, lphi=0, bty='n', colkey=FALSE, zlim=z.lim, d=6, colvar=gii$values, col=cols, clim=c.lims, facets=FALSE, resfac=0.01);
#   mtext(text=ttl, side=3, line=ttl.line, cex=0.6, adj=0);
# 
# 
#   triangle3D(tri=gii$coords, theta=270, phi=0, ltheta=270, lphi=0, bty='n', colkey=FALSE, zlim=z.lim, d=6, colvar=gii$values, col=cols, clim=c.lims, facets=FALSE, resfac=0.01);
#   if (print.scaling.label == TRUE) {
#     mtext(text=paste0("range: ", round(full.min,3), " to ", round(full.max,3), "\n", sc.text), side=1, line=-2.5, cex=0.5, adj=1);
#   }
#   if (print.scaling.label == "flipped") {
#     mtext(text=paste0(sc.text, "\nrange: ", round(full.min,3), " to ", round(full.max,3)), side=1, line=-2.5, cex=0.5, adj=1);
#   }
# }

## SEE:
# source("/data/nil-external/ccp/JosetAEtzel/DMCC_files/giftiPlottingFunctions.R");  # surface plotting functions 
# 
# 
# plot.surface <- function(
#   overlay, underlay, 
#   facet.names = c("left_lateral", "left_medial", "right_lateral", "right_medial"),
#   dims = c(1, 4),
#   facet.params
#   ){
#   # underlay <- hcp; overlay <- mmp
#   # facet.names <- c("left_lateral", "left_medial", "right_lateral", "right_medial")
#   # dims = c(1, 4)
#   
#   ## TODO
#   ## limits
#   ## titles + text
#   ## superior & inferior views?
# 
#   ## input validation
#   
#   is.bad.type <- !all(is.list(overlay) | is.list(underlay) | is.character(facet.names) | is.numeric(dims))
#   if (is.bad.type) stop("input(s) of wrong type")
#   if (any(length(overlay) !=2, length(underlay) != 2)) stop("over/underlay not of length 2")
#   if (any(vapply(c(overlay, underlay), class, character(1)) != "gifti")) stop("bad image classes")
#   if (any(!names(c(overlay, underlay)) %in% c("L", "R"))) stop("over / underlay names not %in% c('L', 'R')")
#   if (any(!facet.names %in% c("left_lateral", "left_medial", "right_lateral", "right_medial"))) stop("bad facet.names")
#   if (length(dims) != 2) stop("length(dims) != 2")
#   if (prod(dims) != length(facet.names)) stop("prod(dims) != length(facet.names)")
#   if (!missing(facet.params)) {
#     are.bad.facet.params <- any(
#       !is.list(facet.params),  ## type ok?
#       !unique(vapply(facet.params, length, numeric(1))) == 2,  ## lengths ok?
#       !length(facet.params) == length(facet.names),  ## length ok?
#       !names(facet.params) %in% c("left_lateral", "left_medial", "right_lateral", "right_medial"),  ## names ok?
#       !unique(c(vapply(facet.params, names, character(2)))) %in% c("theta", "phi")  ## element names ok?
#     )
#     if (are.bad.facet.params) stop("bad facet.params")
#   }
# 
#   ## default orientations
#   
#   if (missing(facet.params)) {
#     
#     facet.params <- list(
#       right_medial  = c(theta = 270, phi = 0),
#       right_lateral = c(theta = 90, phi = 0),
#       left_medial   = c(theta = 90, phi = 0),
#       left_lateral  = c(theta = 270, phi = 0)
#     )
#     
#   }
#   
#   ## extract data from gifti objects and format
#   
#   triangles <- lapply(underlay, function(.) c(t(.[["data"]][["triangle"]])) + 1)
#   indices   <- lapply(triangles, function(.) .[seq(1, length(.), 3)])
#   pointsets <- lapply(underlay, function(.) .[["data"]][["pointset"]])
#   values    <- lapply(overlay, function(.) .[["data"]][[1]][, 1])
#   
#   ## assign colors to triangles
#   
#   coords <- list(L = pointsets$L[triangles$L, ], R = pointsets$R[triangles$R, ])
#   values <- list(L = values$L[indices$L], R = values$R[indices$R])
#   
#   ## draw
#   
#   par(mfrow = dims)
#   
#   for (facet.i in facet.names) {
#     # facet.i = "right_medial"
#     
#     hemi.i <- switch(grepl("right", facet.i) + 1, "L", "R")
#     
#     plot3D::triangle3D(
#       tri    = coords[[hemi.i]],  ## overlay
#       colvar = values[[hemi.i]],  ## underlay
#       theta  = facet.params[[facet.i]]["theta"],  ## "polar angle" / "colatitute"
#       phi    = facet.params[[facet.i]]["phi"],  ## "azithmuthal angle"
#       ## (see https://en.wikipedia.org/wiki/Spherical_coordinate_system)
#       ltheta = facet.params[[facet.i]]["theta"],  ## lighting source
#       lphi   = facet.params[[facet.i]]["phi"],
#       zlim   =  c(-60, 90), 
#       d      = 6, ## ?
#       resfac = 0.01,  ## resolution factor
#       bty    = "n",  ## background type
#       colkey = FALSE,  ## color key
#       # col    = cols,  ## colors
#       # clim   = c.lims,  ## color limits
#       facets = FALSE
#     )
#     
#   }
# }
# 
# plot.surface(mmp, hcp)

