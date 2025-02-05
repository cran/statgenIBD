#' Read IBD probabilities
#'
#' Read a file with IBD probabilities computed by the RABBIT software package,
#' see [RABBIT](https://github.com/Biometris/RABBIT) for details.
#' It is possible to additionally read the pedigree file that is also used by
#' RABBIT. Reading this file allows for plotting the pedigree.
#'
#' @param infile A character string, a link to a .csv file with IBD
#' probabilities. Compressed .csv files with extension ".gz" or ".bz2" are
#' supported as well.
#' @param pedFile A character string, a link to a .csv file with pedigree
#' information as used by RABBIT as input. Compressed .csv files with extension
#' ".gz" or ".bz2" are supported as well.
#'
#' @returns An \code{IBDprob} object with map and markers corresponding to the
#' imported information in the imported .csv file.
#'
#' @examples
#' \dontrun{
#' ## Read RABBIT data for barley.
#' genoFile <- system.file("extdata/barley", "barley_magicReconstruct.zip",
#'                        package = "statgenIBD")
#' barleyIBD <- readRABBIT(unzip(genoFile, exdir = tempdir()))
#' }
#'
#' @references Zheng, Chaozhi, Martin P Boer, and Fred A Van Eeuwijk.
#' “Recursive Algorithms for Modeling Genomic Ancestral Origins in a Fixed
#' Pedigree.” G3 Genes|Genomes|Genetics 8 (10): 3231–45.
#' https://doi.org/10.1534/G3.118.200340.
#'
#' @importFrom utils hasName unzip
#' @export
readRABBIT <- function(infile,
                       pedFile = NULL) {
  if (missing(infile) || !is.character(infile) || length(infile) > 1 ||
      file.access(infile, mode = 4) == -1 ||
      ## Compressed .csv files can be read by fread and should be
      ## allowed as inputs as well.
      !(tools::file_ext(infile) == "csv" ||
        (tools::file_ext(infile) %in% c("gz", "bz2") &&
         tools::file_ext(tools::file_path_sans_ext(infile)) == "csv"))) {
    stop("infile should be a character string indicating a readable .csv file")
  }
  if (!is.null(pedFile) &&
      (!is.character(pedFile) || length(pedFile) > 1 ||
       file.access(pedFile, mode = 4) == -1 ||
       ## Compressed .csv files can be read by fread and should be
       ## allowed as inputs as well.
       !(tools::file_ext(pedFile) == "csv" ||
         (tools::file_ext(pedFile) %in% c("gz", "bz2") &&
          tools::file_ext(tools::file_path_sans_ext(pedFile)) == "csv")))) {
    stop("pedFile should be a character string indicating a readable .csv file")
  }
  ## Read first header to determine whether infile is created
  ## with julia or with mathematica
  con <- file(infile, "r")
  firstHeader <- readLines(con = con, n = 1)
  close(con)
  if (firstHeader == "RABBIT,designinfo") {
    rabbitRes <- readRABBITJulia(infile)
  } else {
    rabbitRes <- readRABBITMathematica(infile)
  }
  map <- rabbitRes$map
  founderProbs <- rabbitRes$founderProbs
  genoNames <- rownames(founderProbs)
  parents <- dimnames(founderProbs)[[3]]
  popType <- "RABBIT"
  if (is.null(pedFile) && is.null(rabbitRes$pedDat)) {
    genoCross <- NULL
  } else {
    if (is.null(rabbitRes$pedDat)) {
      ## Read pedigree.
      pedDat <- data.table::fread(pedFile, skip = "Generation",
                                  data.table = FALSE, fill = TRUE)
    } else {
      ## Pedigree is included and read from julia input.
      pedDat <- rabbitRes$pedDat
    }
    ## Get offspring.
    offDat <- pedDat[pedDat[["Generation"]] %in% genoNames, ]
    offDat[["ID"]] <- offDat[["Generation"]]
    ## Construct genoCross.
    genoCross <- offDat[c("MemberID", "Generation")]
    colnames(genoCross) <- c("cross", "geno")
  }
  ## Create IBDprob object.
  res <- structure(list(map = map,
                        markers = founderProbs,
                        popType = popType,
                        parents = parents,
                        pedigree = NULL),
                   class = c("IBDprob", "list"),
                   genoCross = genoCross)
  return(res)
}

#' Read infile for Mathematica
#'
#' @noRd
#' @keywords internal
readRABBITMathematica <- function(infile) {
  ## Read map and marker probabilities.
  markMap <- data.table::fread(infile, skip = "haploprob", fill = TRUE,
                               data.table = FALSE)
  ## Extract map.
  map <- data.frame(chr = as.numeric(markMap[3, -1]),
                    pos = as.numeric(markMap[4, -1]),
                    row.names = as.character(markMap[2, -1]))
  markMap <- markMap[5:(which(markMap[[2]] == "ibdprob") - 1), ]
  ## Get names of genotypes and compute number of founder alleles per genotype.
  genoNames <- unique(sapply(X = strsplit(x = markMap[, 1],
                                          split = "_haplotype"), FUN = "[[", 1))
  nAlleles = nrow(markMap) / length(genoNames)
  ## Convert markers to 3D array.
  founderProbs <- array(dim = c(length(genoNames), nrow(map), nAlleles))
  for (i in 1:nrow(map)) {
    founderProbs[, i, ] <- matrix(as.numeric(markMap[, i + 1]),
                                  ncol = nAlleles, byrow = TRUE)
  }
  ## Read parent names from file.
  parents <- data.table::fread(infile, header = FALSE, nrows = nAlleles + 2,
                               skip = "haplotypes in order",
                               data.table = FALSE, select = 3)
  parents <- as.character(parents[3:nrow(parents), ])
  ## Add dimnames to markers: genotypes x markers x parents.
  dimnames(founderProbs) <- list(genoNames, rownames(map), parents)
  res <- list(map = map, founderProbs = founderProbs)
  return(res)
}


#' Read infile for Julia
#'
#' @noRd
#' @keywords internal
readRABBITJulia <- function(infile) {
  ## Get the line numbers of the RABBIT csv table blocks.
  rabbitHeaderLineIndexes <- getRabbitHeaderLines(infile)
  ## Get the pedigree.
  pedDat <- data.table::fread(infile,
                              skip = rabbitHeaderLineIndexes$designinfo$lineNr,
                              nrows = rabbitHeaderLineIndexes$designinfo$numLines,
                              data.table = FALSE,
                              header = TRUE)
  colnames(pedDat) <- c("MemberID", "MotherID", "FatherID")
  pedDat[["Generation"]] <- 0
  ## Get the founders.
  founderInfo <-
    data.table::fread(infile,
                      skip = rabbitHeaderLineIndexes$founderinfo$lineNr,
                      nrows = rabbitHeaderLineIndexes$founderinfo$numLines,
                      data.table = FALSE,
                      header = TRUE)
  ## Get the offspring.
  offspringInfo <-
    data.table::fread(infile,
                      skip = rabbitHeaderLineIndexes$offspringinfo$lineNr,
                      nrows = rabbitHeaderLineIndexes$offspringinfo$numLines,
                      data.table = FALSE,
                      header = TRUE)
  offPed <- offspringInfo[c("individual", "member")]
  colnames(offPed) <- c("Generation", "MemberID")
  offPed[["MotherID"]] <- offPed[["FatherID"]] <- NA
  pedDat <- rbind(pedDat, offPed)
  ## Get the markers.
  founderGeno <-
    data.table::fread(infile,
                      skip = rabbitHeaderLineIndexes$foundergeno$lineNr,
                      nrows = rabbitHeaderLineIndexes$foundergeno$numLines,
                      data.table = FALSE,
                      header = TRUE)
  ## Extract map.
  map <- data.frame(chr = founderGeno[, 2],
                    pos = as.numeric(founderGeno[, 3]),
                    row.names = as.character(founderGeno[, 1]))
  ## Prepare 3D array.
  nMarkers <- nrow(founderGeno)
  nOffspring <- nrow(offspringInfo)
  nFounders <- nrow(founderInfo)
  founderProbs <- array(data = NA_real_,
                        dim = c(nOffspring, nMarkers, nFounders),
                        dimnames = list(offspringInfo[["individual"]],
                                        founderGeno[["marker"]],
                                        founderInfo[["individual"]]))
  ## Fill the 3D array (offspring, markers, founders).
  chrStart <- 1
  chrLength <- 0
  for (i in 1:rabbitHeaderLineIndexes$haploprob$numLines) {
    ## Read the probabilities.
    ## reading everything in a table leads to memory problems for large files.
    haploProb <-
      data.table::fread(infile, sep = ",",
                        skip = rabbitHeaderLineIndexes$haploprob$lineNr + i,
                        nrows = 1,
                        header = FALSE)
    markerindex <- as.numeric(stringi::stri_split_regex(haploProb[[5]], "[|]",
                                                        simplify = TRUE))
    haploindex <- as.numeric(stringi::stri_split_regex(haploProb[[6]], "[|]",
                                                       simplify = TRUE))
    probabilities <- as.numeric(stringi::stri_split_regex(haploProb[[7]], "[|]",
                                                          simplify = TRUE))
    mat <- Matrix::sparseMatrix(i = markerindex, j = haploindex,
                                x = probabilities,
                                dims = c(haploProb[[3]], haploProb[[4]]),
                                check = FALSE)
    if (i %% nOffspring == 1) {
      chrStart <- chrStart + chrLength
      chrLength <- dim(mat)[1]
      chrEnd <- chrStart + chrLength - 1
    }
    founderProbs[if (i %% nOffspring == 0) nOffspring else i %% nOffspring,
                 chrStart:chrEnd, ] <- as.numeric(mat)
  }
  res <- list(map = map, founderProbs = founderProbs, pedDat = pedDat)
  return(res)
}

#' @noRd
#' @keywords internal
getRabbitHeaderLines <- function(filepath) {
  con <- file(filepath, "r")
  magicHeaders <- list()
  i <- 0
  curHeaderLineCount <- 0
  curHeader <- NULL
  while (TRUE) {
    line <- readLines(con, n = 1)
    curHeaderLineCount <- curHeaderLineCount + 1
    if (!length(line)) {
      magicHeaders[[curHeader]]$numLines = curHeaderLineCount - 2
      break
    } else if (!is.null(line) & startsWith(line, "RABBIT,")) {
      if (!is.null(curHeader)) {
        magicHeaders[[curHeader]]$numLines = curHeaderLineCount - 2
      }
      header <- gsub("RABBIT,", "", line)
      curHeader <- header
      curHeaderLineCount <- 0
      magicHeaders[[header]] <- list(lineNr = i + 1, numLines = 0)
    }
    i <- i + 1
  }
  close(con)
  return(magicHeaders);
}

