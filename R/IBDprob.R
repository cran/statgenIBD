#' Summary function for objects of class IBDprob
#'
#' Prints a short summary for objects of class \code{IBDprob}. The
#' summary consists of the population type, number of evaluation points,
#' number of individuals and names of the parents in the object.
#'
#' @param object An object of class \code{IBDprob}.
#' @param ... Not used.
#'
#' @return No return value, a summary is printed.
#'
#' @examples
#' ## Compute IBD probabilities for Steptoe Morex.
#' SxMIBD <- calcIBD(popType = "DH",
#'                   markerFile = system.file("extdata/SxM", "SxM_geno.txt",
#'                                           package = "statgenIBD"),
#'                   mapFile = system.file("extdata/SxM", "SxM_map.txt",
#'                                         package = "statgenIBD"))
#'
#' ## Print summary
#' summary(SxMIBD)
#'
#' @export
summary.IBDprob <- function(object,
                            ...) {
  cat("population type: ", object$popType, "\n")
  cat("Number of evaluation points: ", nrow(object$markers), "\n")
  cat("Number of individuals: ", ncol(object$markers),"\n")
  cat("Parents: ", object$parents, "\n")
}

#' Concatenate function for objects of class IBDprob
#'
#' Concatenates objects of class \code{IBDprob}. All objects that are
#' concatenated  should have the same population type and the same map. The
#' function is mainly meant for combining information for multiple crosses
#' with overlapping parents.
#'
#' @param ... Objects of class \code{IBDprob}.
#'
#' @return An object of class \code{IBDprob} containing data for all
#' concatenated objects.
#'
#' @examples
#' ## Compute IBD probabilties for AxB.
#' AB <- calcIBD(popType = "F4DH",
#'              markerFile = system.file("extdata/multipop", "AxB.txt",
#'                                   package = "statgenIBD"),
#'              mapFile = system.file("extdata/multipop", "mapfile.txt",
#'                                   package = "statgenIBD"))
#' ## Compute IBD probabilties for Axc.
#' AC <- calcIBD(popType = "F4DH",
#'               markerFile = system.file("extdata/multipop", "AxC.txt",
#'                                       package = "statgenIBD"),
#'               mapFile = system.file("extdata/multipop", "mapfile.txt",
#'                                    package = "statgenIBD"))
#'
#' ## Combine results.
#' ABC <- c(AB, AC)
#'
#' ## Check summary.
#' summary(ABC)
#'
#' @export
c.IBDprob <- function(...) {
  args <- list(...)
  if (length(args) == 1) {
    return(args[[1]])
  }
  for (arg in args) {
    if (!inherits(arg, "IBDprob")) {
      stop("All inputs should be of class IBDprob.\n")
    }
  }
  pops <- unique(sapply(X = args, FUN = `[[`, "popType"))
  if (!length(pops) == 1) {
    stop("All inputs should have the same population type.\n")
  }
  parents <- unique(c(sapply(X = args, FUN = `[[`, "parents")))
  maps <- unique(lapply(X = args, FUN = `[[`, "map"))
  if (!length(maps) == 1) {
    stop("All inputs should have the same map.\n")
  } else {
    map <- maps[[1]]
  }
  markerLst <- lapply(X = args, FUN = `[[`, "markers")
  parentsNw <- unique(unlist(sapply(X = markerLst, FUN = function(mrk) {
    dimnames(mrk)[3]
  })))
  pedLst <- lapply(X = args, FUN = `[[`, "pedigree")
  pedTot <- do.call(rbind, pedLst)
  pedParNw <- unique(pedTot[pedTot[["par1"]] == 0 & pedTot[["par2"]] == 0, ])
  pedOffNw <- pedTot[pedTot[["par1"]] != 0 | pedTot[["par2"]] != 0, ]
  pedNw <- rbind(pedParNw, pedOffNw)
  genoNw <- unlist(sapply(X = markerLst, FUN = colnames))
  nGeno <- sapply(X = markerLst, FUN = ncol)
  genoCross <- data.frame(cross = paste0("cross",
                                         rep(seq_along(nGeno), times = nGeno)),
                          geno = genoNw)
  markersNw <- array(dim = c(nrow(markerLst[[1]]), length(genoNw),
                             length(parentsNw)),
                     dimnames = list(rownames(markerLst[[1]]), genoNw, parentsNw))
  for (i in 1:nrow(markerLst[[1]])) {
    markersNw[i, , ] <- as.matrix(dfBind(lapply(X = markerLst, FUN = function(mrk) {
      as.data.frame(mrk[i, , ])
    })))
  }
  res <- structure(list(map = map,
                        markers = markersNw,
                        popType = pops,
                        parents = parents,
                        pedigree = pedNw,
                        multiCross = TRUE),
                   class = "IBDprob",
                   genoCross = genoCross)
  return(res)
}

#' Plot function for objects of class IBDprob
#'
#' Creates a plot for an object of class \code{IBDprob}. Three types of plot can
#' be made:
#' \itemize{
#' \item{\code{singleGeno}}{ A plot for a single genotype showing the IBD
#' probabilities for all parents across the genome.}
#' \item{\code{allGeno}}{ A plot showing for all genotypes the IBD
#' probabilities of the parent with the highest probability per marker.}
#' \item{\code{pedigree}}{ A plot showing the structure of the pedigree of
#' the population.}
#' }
#'
#' @param x An object of class \code{IBDprob}.
#' @param ... Further arguments. Unused.
#' @param plotType A character string indicating the type of plot that should
#' be made.
#' @param genotype A character string indicating the genotype for which the
#' plot should be made. Only for \code{plotType = "singleGeno"}.
#' @param title A character string, the title of the plot.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE}, only a ggplot object is invisibly returned.
#'
#' @return A ggplot object is invisibly returned.
#'
#' @examples
#' \dontrun{
#' ## Compute IBD probabilities for Steptoe Morex.
#' ## Add extra evaluation positions in between exiting marker positions
#' ## to assure evaluation positions are at most 2 cM apart.
#' SxMIBD_Ext <- calcIBD(popType = "DH",
#'                       markerFile = system.file("extdata/SxM", "SxM_geno.txt",
#'                                               package = "statgenIBD"),
#'                       mapFile = system.file("extdata/SxM", "SxM_map.txt",
#'                                            package = "statgenIBD"),
#'                       evalDist = 2)
#'
#' ## Plot results for genotype dh005.
#' plot(SxMIBD_Ext,
#'      plotType = "singleGeno",
#'      genotype = "dh005")
#'
#' ## Plot results for all genotypes.
#' plot(SxMIBD_Ext,
#'      plotType = "allGeno")
#'
#' ## Plot sturcture of the pedigree.
#' plot(SxMIBD_Ext,
#'      plotType = "pedigree")
#' }
#'
#' @export
plot.IBDprob <- function(x,
                         ...,
                         plotType = c("singleGeno", "allGeno", "pedigree"),
                         genotype,
                         title = NULL,
                         output = TRUE) {
  map <- x$map
  markers <- x$markers
  parents <- x$parents
  pedigree <- x$pedigree
  multiCross <- x$multiCross
  popType <- x$popType
  genoCross <- attr(x = x, which = "genoCross")
  ## Input checks.
  plotType <- match.arg(plotType)
  if (!is.null(title) && (!inherits(title, "character") || length(title) > 1)) {
    stop("title should be a character string.\n")
  }
  if (plotType == "singleGeno") {
    p <- singleGenoPlot(markers = markers, map = map, parents = parents,
                        genotype = genotype, title = title)
  } else if (plotType == "allGeno") {
    p <- allGenoPlot(markers = markers, map = map, parents = parents,
                     title = title)
  } else if (plotType == "pedigree") {
    p <- pedPlot(pedigree = pedigree, offSpring = colnames(markers),
                 popType = popType, multiCross = multiCross,
                 genoCross = genoCross, title = title)
  }
  if (output) {
    plot(p)
  }
  invisible(p)
}
