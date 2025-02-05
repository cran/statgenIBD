#' Write to Flapjack format
#'
#' Export the results of an IBD calculation to Flapjack format so it can be
#' visualized there.
#'
#' @param IBDprob An object of class \code{IBDprob}.
#' @param outFileMap A character string, the full path to the output map file.
#' @param outFileGeno A character string, the full path to the output genotype
#' file.
#'
#' @returns No output. Output files are written to a folder.
#'
#' @examples
#' ## Compute IBD probabilities for Steptoe Morex.
#' SxMIBD <- calcIBD(popType = "DH",
#'                   markerFile = system.file("extdata/SxM", "SxM_geno.txt",
#'                                           package = "statgenIBD"),
#'                   mapFile = system.file("extdata/SxM", "SxM_map.txt",
#'                                        package = "statgenIBD"))
#'
#' ## Write output in Flapjack format to temporary files.
#' writeFlapjack(SxMIBD,
#'               outFileMap = tempfile(fileext = ".txt"),
#'               outFileGeno = tempfile(fileext = ".txt"))
#'
#' @importFrom stats reshape
#' @importFrom utils write.table
#' @export
writeFlapjack <- function(IBDprob,
                          outFileMap = "ibd_map.txt",
                          outFileGeno = "ibd_geno.txt") {
  if (!inherits(IBDprob, "IBDprob")) {
    stop("IBDprob should be an object of class IBDprob.\n")
  }
  chkFile(outFileMap, fileType = "txt")
  chkFile(outFileGeno, fileType = "txt")
  map <- IBDprob$map
  markers <- IBDprob$markers
  parents <- IBDprob$parents
  nPar <- length(parents)
  nGeno <- dim(markers)[[1]]
  nMarkers <- dim(markers)[[2]]
  ## Use first marker to get parent combinations present in output.
  mrk1 <- colnames(markers)[1]
  parVals <- colnames(getProbs(IBDprob, mrk1))[-1]
  parVals <- gsub(pattern = paste0(mrk1, "_"), replacement = "", x = parVals)
  for (i in seq_along(parents)) {
    parVals <- gsub(pattern = parents[i],
                    replacement = paste0(i, "/"), x = parVals)
  }
  parVals <- substring(parVals, first = 1, last = nchar(parVals) - 1)
  ## Create output.
  res <- sapply(X = colnames(markers), FUN = function(marker) {
    mrkProbs <- getProbs(IBDprob, marker)[-1]
    apply(mrkProbs, MARGIN = 1, FUN = function(x) {
      if (length(which(x > (0.5 + 0.15 / nPar))) == 1) {
        parVals[which.max(x)]
      } else {
        "-"
      }
    })
  })
  ## Convert to matrix.
  res <- matrix(data = res, nrow = nGeno, dimnames = dimnames(markers)[1:2])
  resPar <- matrix(data = seq_along(parents), nrow = nPar, ncol = nMarkers,
                   dimnames = list(parents, colnames(markers)))
  res <- rbind(resPar, res)
  ## Write map file.
  cat(file = outFileMap, "# fjFile = MAP\n")
  write.table(map, file = outFileMap,
              quote = FALSE, sep = "\t", na = "-", row.names = TRUE,
              col.names = FALSE, append = TRUE)
  ## Write marker file.
  cat(file = outFileGeno, "# fjFile = GENOTYPE\n\t")
  cat(colnames(res), file = outFileGeno, sep = "\t", append = TRUE)
  cat("\n", file = outFileGeno, append = TRUE)
  write.table(res, file = outFileGeno,
              quote = FALSE, sep = "\t", na = "-", row.names = TRUE,
              col.names = FALSE, append = TRUE)
}
