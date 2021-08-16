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
#' @return No output. Output files are written to a folder.
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
  nGeno <- dim(markers)[[2]]
  nMarkers <- dim(markers)[[1]]
  ## Convert to long format.
  markersLong <- markers3DtoLong(IBDprob)
  markersWide <- reshape(data = markersLong, idvar = c("genotype", "snp"),
                         timevar = "parent", direction = "wide")
  for (i in 1:nrow(markersWide)) {
     parentsI <- parents[markersWide[i , 3:(2 + nPar)] > (0.85 / nPar)]
     markersWide[i, "res"] <- paste0(parentsI, collapse = "/")
  }
  res <- matrix(data = markersWide[["res"]], nrow = nGeno, byrow = TRUE,
                dimnames = dimnames(markers)[2:1])
  ## Write map file.
  cat(file = outFileMap, "# fjFile = MAP\n")
  write.table(map, file = outFileMap,
              quote = FALSE, sep = "\t", na = "-", row.names = TRUE,
              col.names = FALSE, append = TRUE)
  ## Write marker file.
  cat(file = outFileGeno, "# fjFile = GENOTYPE\n\t")
  cat(colnames(res), file = outFileGeno, append = TRUE)
  write.table(res, file = outFileGeno,
              quote = FALSE, sep = "\t", na = "-", row.names = TRUE,
              col.names = FALSE, append = TRUE)
}
