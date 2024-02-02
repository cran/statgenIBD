#' Helper function for plotting mean of probabilities.
#'
#' @noRd
#' @keywords internal
totalCoveragePlot <- function(markers,
                              map,
                              parents,
                              chr = NULL,
                              title = NULL) {
  if (!all(chr %in% map[["chr"]])) {
    stop("chr not found in map.\n")
  }
  ## Restrict map to selected chromosomes.
  if (!is.null(chr)) {
    map <- map[map[["chr"]] %in% chr, ]
    markers <- markers[, colnames(markers) %in% rownames(map), ]
  }
  ## Compute means.
  mrkDat <- markers3DtoLong(markers = markers, parents = parents)
  plotDat <- tapply(X = mrkDat$prob, INDEX = list(mrkDat$snp, mrkDat$parent),
                    FUN = mean)
  plotDat <- data.frame(parent = parents,
                        coverage = 100 * colMeans(plotDat))
  ## Construct title.
  if (is.null(title)) {
    title <- "Total coverage per parent"
  }
  p <- ggplot2::ggplot(plotDat,
                       ggplot2::aes(x = .data[["parent"]],
                                    y = .data[["coverage"]],
                                    fill = .data[["parent"]])) +
    ggplot2::geom_col() +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = "parent", y = "Total coverage (%)", title = title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_line(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank())
  return(p)
}
