#' Extract Probabilities for markers
#'
#' Extract IBD probabilities for one or more markers from an object of class
#' \code{IBDprob}.
#'
#' @param IBDprob An object of class \code{IBDprob}.
#' @param markers A character vector of markers that should be extracted.
#' @param sumProbs Should the probabilities by summed per parent. If \code{TRUE}
#' the probability for e.g. parent A in a cross with parent B will be calculated
#' as pA + 0.5 * pAB. If \code{FALSE} both pA and pAB will be output without
#' further calculations.
#'
#' @return A data.frame with IBD probabilities for the extracted markers in the
#' column and genotypes in the rows.
#'
#' @examples
#' ## Compute IBD probabilities for Steptoe Morex.
#' SxMIBD <- calcIBD(popType = "DH",
#'                   markerFile = system.file("extdata/SxM", "SxM_geno.txt",
#'                                         package = "statgenIBD"),
#'                   mapFile = system.file("extdata/SxM", "SxM_map.txt",
#'                                         package = "statgenIBD"))
#'
#' ## Get probabilities for a single marker.
#' probOne <- getProbs(IBDprob = SxMIBD,
#'                     markers = "plc")
#' head(probOne)
#'
#' ## Get probabilities for a multiple markers.
#' probMult <- getProbs(IBDprob = SxMIBD,
#'                      markers = c("plc", "tuba1"))
#' head(probMult)
#'
#' @export
getProbs <- function(IBDprob,
                     markers,
                     sumProbs = FALSE) {
  ## Checks.
  if (!inherits(IBDprob, "IBDprob")) {
    stop(deparse(substitute(IBDprob)),
         " should be an object of class IBDprob\n")
  }
  if (is.null(markers) || !is.character(markers)) {
    stop("markers should be a character vector\n")
  }
  missMrk <- markers[!markers %in% rownames(IBDprob$markers)]
  if (length(missMrk) > 0) {
    stop("The following markers are not in ", deparse(substitute(IBDprob)), ": ",
         paste(missMrk, collapse = ", "), "\n")
  }
  if (sumProbs) {
    ## probabilities are summed, e.g. pA = pA + 0.5 * pAB
    probs <- as.data.frame(markers3DtoMat(markers = IBDprob$markers,
                                          parents = IBDprob$parents,
                                          markerSel = markers))
  } else {
    ## probabilites are taken directly from array, so leaving pA and PAB in.
    probs <- lapply(X = markers, FUN = function(marker) {
      prob <- IBDprob$markers[marker, , ]
      colnames(prob) <- paste0(marker, "_",
                               substring(colnames(prob), first = 2))
      return(prob)
    })
    probs <- as.data.frame(do.call(cbind, probs))
  }
  genoCross <- attr(x = IBDprob, which = "genoCross")
  if (!is.null(genoCross)) {
    ## Add a cross column only if relevant.
    probs <- merge(genoCross, probs, by.x = "geno", by.y = "row.names")
    probs <- probs[c("cross", "geno",
                     setdiff(colnames(probs), c("cross", "geno")))]
  } else {
    probs[["geno"]] <- rownames(probs)
    rownames(probs) <- NULL
    probs <- probs[c("geno", setdiff(colnames(probs), "geno"))]
  }
  return(probs)
}
