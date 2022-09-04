#' Write pedigree to file.
#'
#' Writes pedigree information to a plain text, tab-delimited file.
#'
#' @param IBDprob An object of class \code{IBDprob} containing the IBD
#' probabilities.
#' @param outFile A character string specifying the path of the output file.
#'
#' @return No output. The output file is created as a result of calling this
#' function.
#'
#' @examples
#' ## Compute IBD probabilities for Steptoe Morex.
#' SxMIBD <- calcIBD(popType = "DH",
#'                  markerFile = system.file("extdata/SxM", "SxM_geno.txt",
#'                                          package = "statgenIBD"),
#'                  mapFile = system.file("extdata/SxM", "SxM_map.txt",
#'                                       package = "statgenIBD"))
#'
#' ## Write predigree information to temporary file.
#' writePedigree(IBDprob = SxMIBD, outFile = tempfile(fileext = ".txt"))
#'
#' @export
writePedigree <- function(IBDprob,
                          outFile) {
  if (!inherits(IBDprob, "IBDprob")) {
    stop("IBDprob should be an object of class IBDprob.\n")
  }
  chkFile(outFile, fileType = "txt")
  pedigree <- IBDprob$pedigree
  if (is.null(pedigree)) {
    stop("IBDprob contains no pedigree information.\n")
  }
  write.table(pedigree, file = outFile, quote = FALSE, sep = "\t",
              na = "-", row.names = FALSE, col.names = TRUE)
}
