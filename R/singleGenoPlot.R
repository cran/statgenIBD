#' Helper function for plotting probabilities for a single genotype
#'
#' @keywords internal
singleGenoPlot <- function(markers,
                           map,
                           parents,
                           genotype,
                           title) {
  if (!inherits(genotype, "character") || length(genotype) > 1) {
    stop("genotype should be a character string.\n")
  }
  if (!genotype %in% rownames(markers)) {
    stop(paste("genotype", genotype, "not defined\n"))
  }
  ## Convert to long format for plotting.
  markersLong <- markers3DtoLong(markers = markers, parents = parents)
  ## Merge map info to probabilities.
  plotDat <- merge(markersLong, map, by.x = "snp", by.y = "row.names")
  ## Restrict to selected genotype.
  plotDat <- plotDat[plotDat[["genotype"]] == genotype, ]
  ## Convert chr to factor to keep ordering intact.
  plotDat[["chr"]] <- factor(plotDat[["chr"]],
                             levels = unique(map[["chr"]]))
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
  return(p)
}
