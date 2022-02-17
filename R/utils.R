#' Row bind data.frames
#'
#' Helper function for row binding data.frames with different columns.
#'
#' @param dfList A list of data.frames.
#'
#' @noRd
#' @keywords internal
dfBind <- function(dfList) {
  ## Filter empty data.frames from dfList
  dfList <- Filter(f = function(x) nrow(x) > 0, x = dfList)
  if (length(dfList) == 0) {
    return(data.frame())
  }
  ## Get variable names from all data.frames.
  allNms <- unique(unlist(lapply(dfList, names)))
  ## rbind all data.frames setting values for missing columns to 0.
  do.call(rbind,
          c(lapply(X = dfList, FUN = function(x) {
            nwDat <- sapply(X = setdiff(allNms, names(x)), FUN = function(y) {
              0
            })
            data.frame(c(x, nwDat), check.names = FALSE,
                       stringsAsFactors = FALSE)
          }), make.row.names = FALSE)
  )
}

#' @noRd
#' @keywords internal
chkFile <- function(outFile,
                    fileType = "csv") {
  if (!is.character(outFile) || length(outFile) > 1 ||
      tools::file_ext(outFile) != fileType) {
    stop("outFile should be a single character string ending in .",
         fileType, ".\n")
  }
  if (file.access(dirname(outFile), 2)) {
    stop("No permission to write to ", outFile, ".\n")
  }
}

#' Helper function for converting 3D probability matrix to df.
#'
#' Helper function for converting 3D probability matrix to df.
#'
#' @noRd
#' @keywords internal
markers3DtoLong <- function(markers,
                            parents,
                            markerSel = NULL) {
  ## Restrict markers to selected markers
  if (!is.null(markerSel)) {
    markers <- markers[markerSel, , , drop = FALSE]
  }
  markerCols <- dimnames(markers)[[3]]
  ## Create base data.frame for storing long format data.
  markersLongBase <- expand.grid(snp = dimnames(markers)[[1]],
                                 genotype = dimnames(markers)[[2]])
  markersLong <- NULL
  for (parent in parents) {
    ## Construct parent column.
    parentCol <- paste0("p", parent)
    ## Get other columns containing parent.
    parentSubCols <- markerCols[grep(pattern = parent, x = markerCols)]
    parentSubCols <- parentSubCols[-which(parentSubCols == parentCol)]
    ## Add values for parent to base.
    markersParent <- markersLongBase
    markersParent[["parent"]] <- parent
    ## Compute probability for parent.
    ## (2 * pPar + psubPar) / 2
    markersParent[["prob"]] <-
      c(markers[, , parentCol] +
          apply(X = markers[, , parentSubCols, drop = FALSE],
                MARGIN = 1:2, FUN = sum) / 2)
    ## Add to markersLong
    markersLong <- rbind(markersLong, markersParent)
  }
  return(markersLong)
}

#' Helper function for converting 3D probability matrix to 2D IBDMatrix.
#'
#' Helper function for converting 3D probability matrix to 2D IBDMatrix.
#'
#' @noRd
#' @keywords internal
markers3DtoMat <- function(markers,
                           parents,
                           markerSel = NULL) {
  ## Use markers3DtoLong for summing homozygeous and heterozygeous probs.
  markersLong <- markers3DtoLong(markers = markers, parents = parents,
                                 markerSel = markerSel)
  markersLong[["snpPar"]] <-
    paste0(markersLong[["snp"]], "_", markersLong[["parent"]])
  markersLong[["snpPar"]] <- factor(markersLong[["snpPar"]],
                                    levels = unique(markersLong[["snpPar"]]))
  ## Convert to matrix.
  markerMat <- tapply(X = markersLong[["prob"]],
                      INDEX = list(markersLong[["genotype"]],
                                   markersLong[["snpPar"]]), FUN = I)
  return(markerMat)
}
