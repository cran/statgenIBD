#' Helper function for plotting probabilities for all genotypes
#'
#' @keywords internal
allGenoPlot <- function(markers,
                        map,
                        parents,
                        title) {
  maxVals <- maxPars <- matrix(nrow = dim(markers)[1], ncol = dim(markers)[2],
                               dimnames = list(rownames(markers),
                                               colnames(markers)))
  for (mrk in colnames(markers)) {
    ## Use markers3DtoMat for summing homozygeous and heterozygeous probs.
    mrkMat <- markers3DtoMat(markers = markers, parents = parents,
                             markerSel = mrk)
    ## Get max IBD value and parent per marker-position combination.
    maxVals[, mrk] <- apply(X = mrkMat, MARGIN = 1, FUN = max)
    maxPars[, mrk] <- parents[apply(X = mrkMat, MARGIN = 1, FUN = which.max)]
  }
  nGeno <- nrow(maxVals)
  ## Create plot data.
  plotDat <- data.frame(genotype = factor(rownames(maxVals),
                                          levels = rownames(maxVals)),
                        marker = factor(rep(colnames(maxVals), each = nGeno),
                                        levels = colnames(maxVals)),
                        maxVal = as.vector(maxVals),
                        maxPar = factor(as.vector(maxPars),
                                        levels = parents))
  ## Construct title.
  if (is.null(title)) {
    title <- "IBD probabilities across the genome for all genotypes"
  }
  ## Get positions of start of new chromosomes.
  newChrs <- which(!duplicated(map[["chr"]]))[-1] - 0.5
  ## Get mid points of the chromosomes.
  xMarks <- c(0, newChrs) +
    diff(c(which(!duplicated(map[["chr"]])) - 0.5, nrow(map) - 0.5)) / 2
  chrs <- unique(map[["chr"]])
  p <- ggplot2::ggplot(data = plotDat,
                       ggplot2::aes(x = .data[["marker"]],
                                    y = .data[["genotype"]],
                                    alpha = .data[["maxVal"]],
                                    fill = .data[["maxPar"]])) +
    ggplot2::geom_raster() +
    ggplot2::labs(title = title, x = "Chromosome", y = "Genotypes",
                  fill = "Parent", alpha = "Probability") +
    ggplot2::geom_vline(xintercept = newChrs, linetype = "dashed",
                        color = "black") +
    ggplot2::scale_x_discrete(breaks = levels(plotDat$marker)[round(xMarks)],
                              labels = chrs) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  return(p)
}
