#' Read IBD probabilities from file
#'
#' Reads IBD probabilities from a plain text, tab-delimited .txt or .ibd file.
#' Information about the file format can be found in the vignette (
#' \code{vignette("IBDFileFormat", package = "statgenIBD")}). A data.frame with
#' the map must be specified as well.
#'
#' @param infile A character string specifying the path of the input file.
#' Compressed files with extension ".gz" or ".bz2" are supported as well.
#' @param map A data.frame with columns \code{chr} for chromosome and
#' \code{pos} for position. Positions should be in centimorgan (cM). They
#' should not be cumulative over the chromosomes. Other columns are ignored.
#' Marker names should be in the row names. These should match the marker names
#' in the input file.
#'
#' @returns An object of class \code{IBDprob}.
#'
#' @examples
#' ## Read map for Steptoe Morex.
#' SxMmap <- read.delim(system.file("extdata/SxM", "SxM_map.txt",
#'                                package = "statgenIBD"), header = FALSE)
#' rownames(SxMmap) <- SxMmap$V1
#' SxMmap <- SxMmap[, -1]
#' colnames(SxMmap) <- c("chr", "pos")
#'
#'
#' ## Read IBD probabilities for Steptoe Morex.
#' SxMIBD <- readIBDs(system.file("extdata/SxM", "SxM_IBDs.txt",
#'                                package = "statgenIBD"),
#'                   map = SxMmap)
#'
#' ## Print summary.
#' summary(SxMIBD)
#'
#' @importFrom utils hasName read.table
#' @export
readIBDs <- function(infile,
                     map) {
  if (missing(infile) || !is.character(infile) || length(infile) > 1 ||
      file.access(infile, mode = 4) == -1 ||
      ## Compressed .csv files can be read by fread and should be
      ## allowed as inputs as well.
      !(tools::file_ext(infile) %in% c("txt", "ibd") ||
        (tools::file_ext(infile) %in% c("gz", "bz2") &&
         tools::file_ext(tools::file_path_sans_ext(infile)) %in% c("txt", "ibd")))) {
    stop("infile should be a character string indicating a readable ",
         ".txt or .ibd file.\n")
  }
  if (!is.data.frame(map)) {
    stop("map should be a data.frame.\n")
  }
  if (!all(hasName(x = map, name = c("chr", "pos")))) {
    ## chr and pos are obligatory cols.
    stop("chr and pos should be columns in map.\n")
  }
  if (!is.numeric(map[["pos"]])) {
    stop("pos should be a numeric column in map.\n")
  }
  ## Extract columns and order.
  chkNum <- tryCatch(as.numeric(map[["chr"]]), warning = function(w) w)
  if (is.numeric(map[["chr"]]) || inherits(chkNum, "warning")) {
    map <- map[order(map[["chr"]], map[["pos"]]), c("chr", "pos")]
  } else {
    map <- map[order(chkNum, map[["pos"]]), c("chr", "pos")]
  }
  if (!is.character(map[["chr"]])) {
    map[["chr"]] <- as.character(map[["chr"]])
  }
  ## Read file.
  inDat <- data.table::fread(infile, sep = "\t", header = TRUE)
  ## Check that data has required columns.
  if (!all(colnames(inDat)[1:2] == c("Marker", "Genotype"))) {
    stop("First two columns in infile should be named Marker and Genotype.\n")
  }
  if (ncol(inDat) < 4) {
    stop("At least two parent columns should be present in input.\n")
  }
  inDat[, -c(1,2)] <- apply(inDat[, -c(1,2)], 2, function(x) as.numeric(x))
  ## Check that inDat and map contain the same markers.
  missMrkMap <- rownames(map)[!rownames(map) %in% inDat[["Marker"]]]
  if (length(missMrkMap) > 0) {
    warning(length(missMrkMap), " markers in map are not in infile.\n",
            call. = FALSE)
    map <- map[!rownames(map) %in% missMrkMap, ]
  }
  missMrkDat <- unique(inDat[["Marker"]])[!unique(inDat[["Marker"]]) %in%
                                            rownames(map)]
  if (length(missMrkDat) > 0) {
    warning(length(missMrkDat), " markers in infile are not in map.\n",
            call. = FALSE)
    inDat <- inDat[!inDat[["Marker"]] %in% missMrkDat, ]
  }
  genoNamesIn <- unique(inDat[["Genotype"]])
  markerNamesIn <- unique(inDat[["Marker"]])
  ## Sort input data to get everything in expected order.
  inDat <- inDat[order(inDat[["Marker"]], inDat[["Genotype"]]), ]
  genoNames <- unique(inDat[["Genotype"]])
  markerNames <- unique(inDat[["Marker"]])
  parents <- colnames(inDat)[-c(1, 2)]
  nGeno <- length(genoNames)
  nMarkers <- length(markerNames)
  nPar <- length(parents)
  markers <- array(NA_real_,
                   dim = c(nGeno, nMarkers, nPar),
                   dimnames = list(genoNames, markerNames, parents)
  )
  for (geno in genoNames) {
    for (parent in parents) {
      markers[geno, , parent] <- inDat[inDat$Genotype == geno, parent]
    }
    markers[geno, , ] <- markers[geno, , ] / rowSums(markers[geno, , ])
  }
  markers <- markers[genoNamesIn, markerNamesIn, ]
  res <- structure(list(map = map,
                        markers = markers,
                        popType = NULL,
                        parents = parents),
                   class = c("IBDprob", "list"))
  return(res)
}
