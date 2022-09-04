#' Helper function for plotting probabilities for all genotypes
#'
#' @keywords internal
allGenoPlot <- function(markers,
                        map,
                        parents,
                        title) {
  ## Get max IBD value and parent per marker-position combination.
  maxVals <- apply(X = markers, MARGIN = 2:1, FUN = max)
  maxPars <- parents[apply(X = markers, MARGIN = 2:1, FUN = which.max)]
  nMarkers <- nrow(maxVals)
  ## Create plot data.
  plotDat <- data.frame(marker = factor(rownames(maxVals),
                                        levels = rownames(maxVals)),
                        genotype = rep(colnames(maxVals), each = nMarkers),
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
  return(p)
}
