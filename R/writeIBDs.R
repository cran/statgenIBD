#' Write IBD probabilities to file.
#'
#' Writes IBD probabilities to a plain text, tab-delimited .txt or .ibd file.
#' Information about the file format can be found in the vignette (
#' \code{vignette("IBDFileFormat", package = "statgenIBD")}).
#'
#' @param IBDprob An object of class \code{IBDprob} containing the IBD
#' probabilities.
#' @param outFile A character string specifying the path of the output file.
#' @param decimals An integer value specifying the number of decimals to include
#' in writing the output file.
#' @param minProb A numerical value between zero and 1 / number of parents,
#' specifying the minimum probability cutoff value. Probabilities below this
#' cutoff are set to zero and other probabilities are rescaled to make sure that
#' the probabilities sum up to one.
#' @param compress Should the output be compressed to .gz format?
#'
#' @returns No output. The output file is created as a result of calling this
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
#' ## Write IBDs to temporary files.
#' writeIBDs(IBDprob = SxMIBD, outFile = tempfile(fileext = ".txt"))
#'
#' ## Write IBDs to file, set values <0.05 to zero and only print 3 decimals.
#' writeIBDs(IBDprob = SxMIBD, outFile = tempfile(fileext = ".txt"),
#'          decimals = 3, minProb = 0.05)
#'
#' @export
writeIBDs <- function(IBDprob,
                      outFile,
                      decimals = 6,
                      minProb = 0,
                      compress = FALSE) {
  if (!inherits(IBDprob, "IBDprob")) {
    stop("IBDprob should be an object of class IBDprob.\n")
  }
  chkFile(outFile, fileType = c("txt", "ibd"))
  if (!is.numeric(decimals) || length(decimals) > 1) {
    stop("decimals should be a single numerical value.\n")
  }
  if (!is.numeric(minProb) || length(minProb) > 1 || minProb < 0 ||
      minProb >= 1 / length(IBDprob$parents)) {
    stop("minProb should be a numerical value between 0 and ",
         "1 / number of parents.\n")
  }
  markers <- IBDprob$markers
  markerNames <- colnames(markers)
  genoNames <- rownames(markers)
  parents <- IBDprob$parents
  fmt <- if (decimals > 0) paste0("%#.", decimals, "f") else "%f"
  if (minProb > 0) {
    ## Set values < minProb to zero and rescale.
    markers[markers < minProb] <- 0
    markers <- simplify2array(apply(X = markers, MARGIN = 1, FUN = function(x) {
      x / rowSums(x)
    }, simplify = FALSE))
    markers <- aperm(markers, c(1, 3, 2))
  }
  ## Use markers3DtoLong for summing homozygeous and heterozygeous probs.
  mrkDat <- markers3DtoLong(markers = markers, parents = parents)
  ## Create data for export.
  exportDat <- reshape(mrkDat, direction  = "wide", v.names = "prob",
                       timevar = "parent", idvar = c("snp", "genotype"))
  exportDat <- exportDat[c("snp", "genotype", paste0("prob.", parents))]
  colnames(exportDat) <- c("Marker", "Genotype", parents)
  for (parent in parents) {
    ## Format output.
    ## Trailing zeros are removed and number of decimals is adjusted.
    exportDat[[parent]] <-
      sub("\\.$", "",
          sub("0+$", "",
              sprintf(exportDat[[parent]], fmt = fmt)
          )
      )
  }
  if (compress) {
    outFile <- paste0(outFile, ".gz")
  }
  data.table::fwrite(exportDat, file = outFile,
                     quote = FALSE, sep = "\t", na = "-", row.names = FALSE,
                     col.names = TRUE)
}
