overlay      = ov
underlay     = hcp
hues         = viridis(500)
pos.lims     = Inf
neg.lims     = -Inf
thresh0      = TRUE
facet.names  = c("left_lateral", "left_medial", "right_lateral", "right_medial")
dims         = c(1, 4)
anno         = FALSE
# facet.params  ## optional




plot_surface <- function (
  overlay, 
  underlay       = hcp, 
  hues           = viridis(500),
  pos.max        = NULL,
  pos.min        = NULL, 
  neg.max        = NULL,
  neg.min        = NULL,
  thresh0        = TRUE,
  facet.names    = c("left_lateral", "left_medial", "right_lateral", "right_medial"), 
  dims           = c(1, 4), 
  anno           = FALSE,
  facet.params,  ## optional
  anno.cex       = 0.5,
  not.rstudio.gd = TRUE,
  ...
) {
  
  
  ## input validation ----
  
  is.bad.type <- !all(is.list(overlay) | is.list(underlay) | is.character(facet.names) | is.numeric(dims))
  if (is.bad.type) stop("input(s) of wrong type")
  if (any(length(overlay) != 2, length(underlay) != 2)) stop("over/underlay not of length 2")
  if (any(vapply(underlay, class, character(1)) != "gifti")) stop("bad image classes")
  if (any(!names(c(overlay, underlay)) %in% c("L", "R"))) stop("over / underlay names not %in% c('L', 'R')")
  if (any(!facet.names %in% c("left_lateral", "left_medial", "right_lateral", "right_medial"))) {
    stop("bad facet.names")
  }
  if (length(dims) != 2) stop("length(dims) != 2")
  if (prod(dims) != length(facet.names)) stop("prod(dims) != length(facet.names)")
  if (!missing(facet.params)) {
    are.bad.facet.params <- any(
      !is.list(facet.params), 
      !unique(vapply(facet.params, length, numeric(1))) == 2,
      !length(facet.params) == length(facet.names), 
      !names(facet.params) %in% c("left_lateral", "left_medial", "right_lateral", "right_medial"), 
      !unique(c(vapply(facet.params, names, character(2)))) %in% c("theta", "phi")
    )
    if (are.bad.facet.params) stop("bad facet.params")
  }
  if (missing(facet.params)) {
    facet.params <- list(
      right_medial  = c(theta = 270, phi = 0), 
      right_lateral = c(theta = 90, phi = 0), 
      left_medial   = c(theta = 90, phi = 0),
      left_lateral  = c(theta = 270, phi = 0)
    )
  }
  if (any(!c(pos.max, pos.min) >= 0) | any(!c(neg.min, neg.max) >= 0)) stop("bad lims")
  
  
  ## threshold ----
  
  ## pos max
  
  if (length(pos.max) == 0) {
    ## if null, don't threshold
  } else if (pos.max == 0) {
    ## if 0, plot no positive
    overlay$L[overlay$L > 0] <- NA
    overlay$R[overlay$R > 0] <- NA
  } else {
    overlay$L[overlay$L > pos.max] <- pos.max  ## greater than pos.max, set to max
    overlay$R[overlay$R > pos.max] <- pos.max
  }
  
  ## pos min
  
  if (length(pos.min) == 0) {
    ## if null, don't threshold
  } else {
    overlay$L[overlay$L <= pos.min & overlay$L > 0] <- NA  ## less than pos.max, set to 0
    overlay$R[overlay$R <= pos.min & overlay$R > 0] <- NA
  }
  
  ## neg min
  
  if (length(neg.min) == 0) {
    ## if null, don't threshold
  } else if (neg.min == 0) {
    ## plot no negative
    overlay$L[overlay$L < 0] <- NA
    overlay$R[overlay$R < 0] <- NA
  } else {
    overlay$L[overlay$L < neg.min] <- neg.min  ## more extreme than neg.min, set to neg.min
    overlay$R[overlay$R < neg.min] <- neg.min
  }
  
  ## neg max
   
  if (length(neg.max) == 0) {
    ## if null, don't threshold
  } else {
    overlay$L[overlay$L > neg.max & overlay$L < 0] <- NA  ## less extreme than neg.max, set to 0
    overlay$R[overlay$R > neg.max & overlay$R < 0] <- NA
  }

  ## handling zero values (typ sub-cortical regions)
  
  if (thresh0) {
    overlay$L[overlay$L == 0] <- NA
    overlay$R[overlay$R == 0] <- NA
  }
  
  
  ## prepare objects for plot3D ----
  
  triangles <- lapply(underlay, function(.) c(t(.[["data"]][["triangle"]])) + 1)
  indices   <- lapply(triangles, function(.) .[seq(1, length(.), 3)])
  pointsets <- lapply(underlay, function(.) .[["data"]][["pointset"]])
  coords    <- list(L = pointsets$L[triangles$L, ], R = pointsets$R[triangles$R, ])
  values    <- list(L = overlay$L[indices$L], R = overlay$R[indices$R])
  
  if (not.rstudio.gd) if (interactive()) dev.new(noRStudioGD = TRUE)  ## don't draw to rstudio device (slow!)
  
  
  ## plot ----
  
  par(mfrow = dims)
  
  for (facet.i in facet.names) {
    
    hemi.i <- switch(grepl("right", facet.i) + 1, "L", "R")
    
    plot3D::triangle3D(
      tri    = coords[[hemi.i]], 
      colvar = values[[hemi.i]], 
      theta  = facet.params[[facet.i]]["theta"], 
      phi    = facet.params[[facet.i]]["phi"], 
      ltheta = facet.params[[facet.i]]["theta"], 
      lphi   = facet.params[[facet.i]]["phi"], 
      zlim   = c(-60, 90), 
      d      = 6, 
      resfac = 0.01, 
      bty    = "n", 
      colkey = switch((facet.i == "left_lateral") + 1, FALSE, list(side = 1, cex.axis = 1.5)),
      col    = hues,
      facets = FALSE,
      ...
    )
    
  }
  
  if (anno) {
    mtext(
      text = paste0("range displayed: ", paste0(round(range(overlay, na.rm = TRUE), 2), collapse = ", ")),
      side = 1, line = -2.5, cex = 0.5, adj = 1
    )
  }
  
}


values <- stats.coding.targt$beta


# dat <- stats.coding.targt[, c("roi", "beta", "num.roi")]
# dat$hemi <- substr(dat$roi, nchar(dat$roi), nchar(dat$roi))
# col.values = "beta"
# col.roi.num = "num.roi"
# col.hemi = "hemi"
# template = mmp


build_overlay <- function (
  d, col.values, col.roi.num = "num.roi", col.hemi = "hemi", template = mmp
  ) {
  
  if (data.table::is.data.table(d) | tibble::is.tibble(d)) d <- as.data.frame(d)
  if (!is.data.frame(d)) stop("d is not data.frame")
  if (any(!c(col.roi.num, col.hemi, col.values) %in% names(d))) stop("names not in names(d)")
  if (!identical(sort(names(template)), c("L", "R"))) stop("template names wrong (not L, R)")
  
  overlay <- setNames(vector("list", 2), c("L", "R"))
  
  for (hemi.i in seq_along(overlay)) {
    
    ## setup objs for inner loop
    
    name.hemi.i <- names(overlay)[hemi.i]
    di <- d[d[[col.hemi]] == name.hemi.i, ]
    template.i <- template[[name.hemi.i]]
    
    for (roi.j in seq_len(nrow(di))) {
      
      num.roi.j <- di[[col.roi.num]][roi.j]
      template.i[template.i == num.roi.j] <- di[roi.j, col.values]
      
    }
    
    overlay[[name.hemi.i]] <- template.i
    
  }
  
  overlay
  
}
  


rm(d, col.values, col.roi.num, col.hemi, overlay, hemi.i, di, template.i, num.roi.j)

dat <- stats.coding.targt[, c("roi", "beta", "num.roi")]
dat$hemi <- substr(dat$roi, nchar(dat$roi), nchar(dat$roi))

build_overlay(dat, col.values = "beta")

dat$beta2 <- NA
build_overlay(dat, col.values = "beta2")


dat %>% build_overlay(col.values = "beta") %>% plot_surface(pos.min = 0.05)


v <- dat %>% build_overlay(col.values = "beta") %>% unlist

sum(v < 0.05 & v > 0)

ov <- dat %>% build_overlay(col.values = "beta")

View(dat)

ov$L[is.na(ov$L)]


unique(unlist(ov, use.names = FALSE))
sum(is.na(unlist(ov, use.names = FALSE)))


## original ----

# plot_surface <- function (
#   overlay, 
#   underlay     = hcp, 
#   hues         = viridis(500),
#   pos.lims     = Inf,
#   neg.lims     = -Inf, 
#   thresh0      = TRUE,
#   facet.names  = c("left_lateral", "left_medial", "right_lateral", "right_medial"), 
#   dims         = c(1, 4), 
#   anno         = FALSE,
#   facet.params  ## optional
#   ) {
#   
#     is.bad.type <- !all(is.list(overlay) | is.list(underlay) | is.character(facet.names) | is.numeric(dims))
#     if (is.bad.type) stop("input(s) of wrong type")
#     if (any(length(overlay) != 2, length(underlay) != 2)) stop("over/underlay not of length 2")
#     if (any(vapply(underlay, class, character(1)) != "gifti")) stop("bad image classes")
#     if (any(!names(c(overlay, underlay)) %in% c("L", "R"))) stop("over / underlay names not %in% c('L', 'R')")
#     if (any(!facet.names %in% c("left_lateral", "left_medial", "right_lateral", "right_medial"))) {
#       stop("bad facet.names")
#     }
#     if (length(dims) != 2) stop("length(dims) != 2")
#     if (prod(dims) != length(facet.names)) stop("prod(dims) != length(facet.names)")
#     if (!missing(facet.params)) {
#         are.bad.facet.params <- any(
#           !is.list(facet.params), 
#           !unique(vapply(facet.params, length, numeric(1))) == 2,
#           !length(facet.params) == length(facet.names), 
#           !names(facet.params) %in% c("left_lateral", "left_medial", "right_lateral", "right_medial"), 
#           !unique(c(vapply(facet.params, names, character(2)))) %in% c("theta", "phi")
#         )
#         if (are.bad.facet.params) stop("bad facet.params")
#     }
#     if (missing(facet.params)) {
#         facet.params <- list(
#           right_medial  = c(theta = 270, phi = 0), 
#           right_lateral = c(theta = 90, phi = 0), 
#           left_medial   = c(theta = 90, phi = 0),
#           left_lateral  = c(theta = 270, phi = 0)
#           )
#     }
#     
#     
#     if (pos.lims == max(Inf)) {  ## plot all data
#     } else if (length(pos.lims) == 0) {  ## plot no positive
#       overlay$L[overlay$L > 0] <- NA
#       overlay$R[overlay$R > 0] <- NA
#     } else {  ## plot within limits
#       overlay$L[overlay$L > pos.lims[2] | (overlay$L <= pos.lims[1] & overlay$L > 0)] <- NA
#       overlay$R[overlay$R > pos.lims[2] | (overlay$R <= pos.lims[1] & overlay$R > 0)] <- NA
#     }
#     if (neg.lims == min(-Inf)) {
#     } else if (length(neg.lims) == 0) {
#       overlay$L[overlay$L < 0] <- NA
#       overlay$R[overlay$R < 0] <- NA
#     } else {
#       overlay$L[overlay$L < neg.lims[2] | (overlay$L > neg.lims[1] & overlay$L < 0)] <- NA
#       overlay$R[overlay$R < neg.lims[2] | (overlay$R > neg.lims[1] & overlay$R < 0)] <- NA
#     }
#     if (thresh0) {
#       overlay$L[overlay$L == 0] <- NA
#       overlay$R[overlay$R == 0] <- NA
#     }
#     
#     lim <- vapply(overlay, range, na.rm = TRUE, numeric(2))
#     lim <- c(min(lim[, 1]), max(lim[, 2]))
#     
#     triangles <- lapply(underlay, function(.) c(t(.[["data"]][["triangle"]])) + 1)
#     indices <- lapply(triangles, function(.) .[seq(1, length(.), 3)])
#     pointsets <- lapply(underlay, function(.) .[["data"]][["pointset"]])
#     coords <- list(L = pointsets$L[triangles$L, ], R = pointsets$R[triangles$R, ])
#     values <- list(L = overlay$L[indices$L], R = overlay$R[indices$R])
#     
#     if (interactive()) dev.new(noRStudioGD = TRUE)
#     
#     par(mfrow = dims)
#     
#     for (facet.i in facet.names) {
#         
#       hemi.i <- switch(grepl("right", facet.i) + 1, "L", "R")
#       
#       plot3D::triangle3D(
#         tri    = coords[[hemi.i]], 
#         colvar = values[[hemi.i]], 
#         theta  = facet.params[[facet.i]]["theta"], 
#         phi    = facet.params[[facet.i]]["phi"], 
#         ltheta = facet.params[[facet.i]]["theta"], 
#         lphi   = facet.params[[facet.i]]["phi"], 
#         zlim   = c(-60, 90), 
#         d      = 6, 
#         resfac = 0.01, 
#         bty    = "n", 
#         colkey = switch((facet.i == "left_lateral") + 1, FALSE, list(side = 1, cex.axis = 1.5)),
#         col    = hues,
#         clim   = lim,
#         facets = FALSE
#         )
#       
#     }
#     
#   if (anno) {
#     mtext(
#       text = paste0("range displayed: ", paste0(round(range(overlay, na.rm = TRUE), 2), collapse = ", ")),
#       side = 1, line = -2.5, cex = 0.5, adj = 1
#       )
#   }
#     
# }
# 
# build.overlay <- function (values, overlay = mmp) {
#   
#   out <- setNames(vector("list", 2), c("L", "R"))
#   
#   for (hemi.i in c("L", "R")) {
#     
#     rois <- atlas.key$mmp[atlas.key$mmp$hemi == hemi.i, "num.roi"]
#     overlay.i <- overlay[[hemi.i]]
#     
#     for (roi.i in seq_along(rois)) {
#       overlay.i[overlay.i == rois[roi.i]] <- values[rois[roi.i]]
#       out[[hemi.i]] <- overlay.i
#     }
#         
#   }
#   
#   out
#   
# }
