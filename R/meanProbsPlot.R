#' Helper function for plotting mean of probabilities.
#'
#' @importFrom stats aggregate
#'
#' @noRd
#' @keywords internal
meanProbsPlot <- function(markers,
                          map,
                          parents,
                          chr = NULL,
                          title = NULL) {
  if (!all(chr %in% map[["chr"]])) {
    stop("chr not found in map.\n")
  }
  ## Compute means.
  mrkDat <- markers3DtoLong(markers = markers, parents = parents)
  plotDat <- tapply(X = mrkDat$prob, INDEX = list(mrkDat$snp, mrkDat$parent),
                    FUN = mean)
  colnames(plotDat) <- parents
  ## Restrict map to selected chromosomes.
  if (!is.null(chr)) {
    map <- map[map[["chr"]] %in% chr, ]
    markers <- markers[, colnames(markers) %in% rownames(map), ]
  }
  map$chr <- factor(map$chr, levels = unique(map$chr))
  ## Get the boundary for each of the chromosomes.
  ## Has to be converted to numeric to avoid integer overflow in the next step.
  chrBnd <- aggregate(x = map$pos, by = list(map$chr),
                      FUN = function(p) {as.numeric(max(p))})
  ## Compute cumulative positions.
  addPos <- data.frame(chr = chrBnd[, 1],
                       add = c(0, cumsum(chrBnd[, 2]))[1:nrow(chrBnd)],
                       stringsAsFactors = FALSE)
  ## Compute cumulative position.
  map[["snp"]] <- rownames(map)
  map <- merge(map, addPos, by = "chr", sort = FALSE)
  map$cumPos <- map$pos + map$add
  ## Extract central chromosome positions from map.
  chrs <- levels(map$chr)
  xMarks <- aggregate(x = map$cumPos, by = list(map$chr), FUN = function(x) {
    min(x) + (max(x) - min(x)) / 2
  })[, 2]
  ## convert to data.frame.
  plotDat <- as.data.frame.table(plotDat)
  ## Merge map to get positions.
  plotDat <- merge(map, plotDat, by.x = "snp", by.y = "Var1", sort = FALSE)
  ## Construct title.
  if (is.null(title)) {
    title <- "Coverage per parent"
  }
  p <- ggplot2::ggplot(plotDat,
                       ggplot2::aes(x = .data[["cumPos"]],
                                    y = .data[["Freq"]],
                                    color = .data[["Var2"]],
                                    group = .data[["Var2"]])) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Chromosome", y = "Representation", title = title) +
    ggplot2::scale_x_continuous(breaks = xMarks, labels = chrs,
                                expand = c(0, 0)) +
    ggplot2::scale_color_discrete(name = "parent") +
    ggplot2::geom_vline(data = addPos[-1, ],
                        ggplot2::aes(xintercept = .data[["add"]]),
                        linetype = "dashed") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      axis.line = ggplot2::element_line())
  return(p)
}
