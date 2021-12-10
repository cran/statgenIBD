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
                        multicross = TRUE),
                   class = "IBDprob",
                   genoCross = genoCross)
  return(res)
}

#' Plot function for objects of class IBDprob
#'
#' Creates a plot for an object of class \code{IBDprob}. Two types of plot can
#' be made, a plot for a single genotype showing the IBD probabilities for all
#' parents across the genome (\code{plotType = "singleGeno"}), or a plot showing
#' for all genotypes the probabilities of the parent with the highest
#' probability per marker (\code{plotType = "allGeno"}).
#'
#' @param x An object of class \code{IBDprob}.
#' @param ... Further arguments. Unused.
#' @param plotType A character string indicating the type of plot that should
#' be made.
#' @param genotype A character string indicating the genotype for which the
#' plot should be made. Ignored if \code{plotType = "allGeno"}.
#' @param title A character string, the title of the plot.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE}, only a ggplot object is invisibly returned.
#'
#' @return A ggplot object is invisibly returned.
#'
#' @examples
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
#' @export
plot.IBDprob <- function(x,
                         ...,
                         plotType = c("singleGeno", "allGeno"),
                         genotype,
                         title = NULL,
                         output = TRUE) {
  map <- x$map
  markers <- x$markers
  parents <- x$parents
  ## Input checks.
  plotType <- match.arg(plotType)
  if (!is.null(title) && (!inherits(title, "character") || length(title) > 1)) {
    stop("title should be a character string.\n")
  }
  if (plotType == "singleGeno") {
    if (!inherits(genotype, "character") || length(genotype) > 1) {
      stop("genotype should be a character string.\n")
    }
    if (!genotype %in% dimnames(markers)[[2]]) {
      stop(paste("genotype", genotype, "not defined\n"))
    }
    ## Convert to long format for plotting.
    markersLong <- markers3DtoLong(x)
    ## Merge map info to probabilities.
    plotDat <- merge(markersLong, map, by.x = "snp", by.y = "row.names")
    ## Restrict to selected genotype.
    plotDat <- plotDat[plotDat[["genotype"]] == genotype, ]
    ## Construct title.
    if (is.null(title)) {
      title <- genotype
    }
    p <- ggplot2::ggplot(plotDat,
                         ggplot2::aes_string(x = "pos", y = "parent",
                                             fill = "prob")) +
      ggplot2::geom_tile(width = 3) +
      ggplot2::facet_grid(". ~ chr", scales = "free", space = "free",
                          switch = "both") +
      ggplot2::scale_fill_gradient(low = "white", high = "black") +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::labs(title = title) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        panel.background = ggplot2::element_blank(),
        panel.spacing.x = ggplot2::unit(4, "mm"))
  } else if (plotType == "allGeno") {
    ## Get max IBD value and parent per marker-position combination.
    maxVals <- apply(X = markers, MARGIN = 1:2, FUN = max)
    maxPars <- parents[apply(X = markers, MARGIN = 1:2, FUN = which.max)]
    nGeno <- nrow(maxVals)
    ## Create plot data.
    plotDat <- data.frame(marker = factor(dimnames(maxVals)[[1]],
                                          levels = dimnames(maxVals)[[1]]),
                          genotype = rep(dimnames(maxVals)[[2]], each = nGeno),
                          maxVal = as.vector(maxVals),
                          maxPar = as.vector(maxPars))
    ## Construct title.
    if (is.null(title)) {
      title <- "IBD probabilities across the genome for all genotypes"
    }
    ## Get positions of start of new chromosomes.
    newChrs <- which(!duplicated(map[["chr"]]))[-1] - 0.5
    p <- ggplot2::ggplot(data = plotDat,
                         ggplot2::aes_string(x = "marker", y = "genotype",
                                             alpha = "maxVal",
                                             fill = "maxPar"))+
      ggplot2::geom_raster() +
      ggplot2::labs(title = title, x = "Genome", y = "Genotypes",
                    fill = "Parent", alpha = "Probability") +
      ggplot2::geom_vline(xintercept = newChrs, linetype = "dashed",
                          color = "black") +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(fill = NA),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )
  }
  if (output) {
    plot(p)
  }
  invisible(p)
}
