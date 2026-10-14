# ==============================================================================
# Helper Functions for Geostatistical Analysis & Spatial Variography
# Author: Beth (MSc Spatial Epidemiology)
# Project: Spatial Epidemiology in Public Health - Côte d'Ivoire Assessment
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(sf)
  if (!requireNamespace("geoR", quietly = TRUE)) {
    warning("Package 'geoR' is required for variogram envelopes.")
  }
})

#' Convert an sf object's CRS units from meters to kilometers
#'
#' @param x An sf or sfc object with a projected coordinate reference system (CRS).
#' @return A proj4string character vector with units changed from meters to kilometers.
#' @examples
#' # epsgKM(adm1_rdt_sp)
epsgKM <- function(x) {
  crs <- sf::st_crs(x)
  proj4KM <- gsub(pattern = "+.units=m", replacement = "+units=km", crs$proj4string)
  return(proj4KM)
}

#' Compute Monte Carlo Envelopes for an Empirical Variogram under CSR
#'
#' Permutes spatial data across sample coordinates under the null hypothesis of 
#' Complete Spatial Randomness (CSR) to compute pointwise 95% confidence intervals 
#' for an empirical semivariogram.
#'
#' @param geodata Optional geodata object containing coords and data.
#' @param coords Two-column matrix or dataframe of spatial coordinates.
#' @param data Numeric vector of observations at coords.
#' @param obj.variog A variogram object returned by `geoR::variog`.
#' @param nsim Integer; number of Monte Carlo simulations (default: 99).
#' @param save.sim Logical; whether to retain raw simulation results (default: FALSE).
#' @param messages Logical; whether to print progress messages to the console.
#' @return A list of class 'variogram.envelope' containing lag distances (`u`) 
#'   and lower/upper 95% envelope boundaries (`v.lower`, `v.upper`).
variog_envelope <- function(geodata, 
                            coords = geodata$coords, 
                            data = geodata$data, 
                            obj.variog, 
                            nsim = 99, 
                            save.sim = FALSE, 
                            messages) {
  call.fc <- match.call()
  if (missing(geodata)) {
    geodata <- list(coords = coords, data = data)
  }
  if (missing(messages)) {
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), 
                                         TRUE, getOption("geoR.messages")))
  } else {
    messages.screen <- messages
  }
  
  obj.variog$v <- NULL
  if ((is.matrix(data) || is.data.frame(data))) {
    if (ncol(data) > 1) {
      stop("Envelopes can only be computed for one dataset at a time.")
    }
  }
  
  estimator.type <- if (!is.null(obj.variog$estimator.type)) obj.variog$estimator.type else "classical"
  
  if (abs(obj.variog$lambda - 1) > 1e-04) {
    if (abs(obj.variog$lambda) < 1e-04) {
      data <- log(data)
    } else {
      data <- ((data^obj.variog$lambda) - 1) / obj.variog$lambda
    }
  }
  
  xmat <- unclass(geoR::trend.spatial(trend = obj.variog$trend, geodata = geodata))
  if (obj.variog$trend != "cte") {
    if (is.vector(data)) {
      data <- stats::lm(data ~ xmat + 0)$residuals
      names(data) <- NULL
    } else {
      only.res <- function(y, x) stats::lm(y ~ xmat + 0)$residuals
      data <- apply(data, 2, only.res, x = xmat)
    }
  }
  
  if (messages.screen) {
    cat(paste("variog.env: generating", nsim, "simulations by permuting data values\n"))
  }
  
  simula <- list(coords = coords)
  n.data <- length(data)
  perm.f <- function(i, data, n.data) data[sample(1:n.data)]
  simula$data <- apply(as.matrix(1:nsim), 1, perm.f, data = data, n.data = n.data)
  
  if (messages.screen) {
    cat(paste("variog.env: computing empirical variograms for", nsim, "simulations\n"))
  }
  
  nbins <- length(obj.variog$bins.lim) - 1
  if (obj.variog$direction == "omnidirectional") {
    bin.f <- function(sim) {
      cbin <- vbin <- sdbin <- rep(0, nbins)
      temp <- .C("binit", as.integer(obj.variog$n.data), 
                 as.double(as.vector(coords[, 1])), 
                 as.double(as.vector(coords[, 2])), 
                 as.double(as.vector(sim)), 
                 as.integer(nbins), 
                 as.double(as.vector(obj.variog$bins.lim)), 
                 as.integer(estimator.type == "modulus"), 
                 as.double(max(obj.variog$u)), 
                 as.double(cbin), 
                 vbin = as.double(vbin), 
                 as.integer(FALSE), 
                 as.double(sdbin), 
                 PACKAGE = "geoR")$vbin
      return(temp)
    }
    simula.bins <- apply(simula$data, 2, bin.f)
  } else {
    variog.vbin <- function(x, ...) {
      geoR::variog(geodata = geodata, 
                   data = x, 
                   uvec = obj.variog$uvec, 
                   estimator.type = obj.variog$estimator.type, 
                   nugget.tolerance = obj.variog$nugget.tolerance, 
                   max.dist = obj.variog$max.dist, 
                   pairs.min = obj.variog$pairs.min, 
                   direction = obj.variog$direction, 
                   tolerance = obj.variog$tolerance, 
                   messages.screen = FALSE, ...)$v
    }
    simula.bins <- apply(simula$data, 2, variog.vbin)
  }
  
  simula.bins <- simula.bins[obj.variog$ind.bin, ]
  if (!save.sim) {
    simula$data <- NULL
  }
  if (messages.screen) {
    cat("variog.env: computing confidence envelopes\n")
  }
  
  limits <- apply(simula.bins, 1, stats::quantile, prob = c(0.025, 0.975))
  res.env <- list(u = obj.variog$u, v.lower = limits[1, ], v.upper = limits[2, ])
  if (save.sim) {
    res.env$simulations <- simula$data
  }
  res.env$call <- call.fc
  oldClass(res.env) <- "variogram.envelope"
  return(res.env)
}

#' Plot an Empirical Semivariogram with Monte Carlo Randomness Envelope
#'
#' Computes empirical binned semivariance estimates and superimposes a 
#' 95% simulation envelope generated by permuting point observations.
#'
#' @param coords A 2-column matrix or dataframe of projected spatial coordinates (x, y).
#' @param data Numeric vector of attribute values at `coords` (e.g., empirical logit prevalence).
#' @param bins Integer; number of distance lag bins (default: 15).
#' @param maxdist Numeric; maximum lag distance to evaluate (default: 1/3 of max pairwise distance).
#' @param uvec Optional numeric vector specifying user-defined lag bin cutpoints.
#' @param nsim Integer; number of Monte Carlo simulations for envelope (default: 999).
#' @param color Character; fill color for the envelope ribbon and points (default: "royalblue1").
#' @param xlab Character; label for x-axis (default: "Distance (m)").
#' @param show_nbins Logical; whether to display point-pair counts beneath each point (default: TRUE).
#' @return A ggplot2 plot object.
ggvario <- function(coords, 
                    data, 
                    bins = 15, 
                    maxdist = max(dist(coords)) / 3, 
                    uvec = NULL, 
                    nsim = 999,
                    color = "royalblue1", 
                    xlab = "Distance (m)", 
                    show_nbins = TRUE) {
  if (!requireNamespace("geoR", quietly = TRUE)) {
    stop("Package 'geoR' is required to run ggvario.")
  }
  
  coords <- as.matrix(coords)
  min_dist <- min(stats::dist(coords))
  if (is.null(uvec)) {
    uvec <- seq(min_dist, maxdist, length.out = bins)
  }
  
  empvario <- geoR::variog(coords = coords, data = data, uvec = uvec, messages = FALSE)
  envmc <- variog_envelope(coords = coords, data = data, 
                           obj.variog = empvario, nsim = nsim, messages = FALSE)
  
  dfvario <- data.frame(
    distance = empvario$u, 
    empirical = empvario$v,
    lowemp = envmc$v.lower, 
    upemp = envmc$v.upper, 
    nbins = empvario$n
  )
  
  p1 <- ggplot(dfvario, aes(x = distance, y = empirical, label = nbins)) +
    geom_ribbon(aes(ymin = lowemp, ymax = upemp), fill = color, alpha = 0.25) +
    geom_point(aes(y = empirical), col = "black", fill = color, shape = 21, size = 3) +
    scale_x_continuous(name = xlab, limits = c(0, uvec[length(uvec)]),
                       breaks = round(seq(0, uvec[length(uvec)], length.out = 6))) +
    scale_y_continuous(name = "Semivariance", 
                       limits = c(0, max(c(dfvario$upemp, dfvario$empirical), na.rm = TRUE) * 1.05)) +
    labs(
      title = "Empirical Semivariogram with 95% Monte Carlo Envelope",
      subtitle = paste("Null envelope generated from", nsim, "permutations under CSR")
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      axis.title = element_text(face = "bold", size = 13)
    )
    
  p2 <- p1 + geom_text(vjust = 1.4, size = 3.5, color = "grey20", 
                       nudge_y = - diff(range(dfvario$empirical)) / 22)
  
  if (show_nbins) p2 else p1
}
