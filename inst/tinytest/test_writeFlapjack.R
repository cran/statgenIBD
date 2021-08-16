### Test writeFlapjack function.

## Define file locations.
SxMloc <- system.file("extdata/SxM", "SxM_geno.txt", package = "statgenIBD")
SxMmap <- system.file("extdata/SxM", "SxM_map.txt", package = "statgenIBD")

## IBD calculations .
SxMIBD <- calcIBD(popType = "DH", markerFile = SxMloc, mapFile = SxMmap)

mapOut <- tempfile(fileext = ".txt")
genoOut <- tempfile(fileext = ".txt")

## Check that input checks are working correctly.
expect_error(writeFlapjack(1),
             "IBDprob should be an object of class IBDprob")
expect_error(writeFlapjack(SxMIBD, outFileMap = 1),
             "outFile should be a single character string ending in .txt")
expect_error(writeFlapjack(SxMIBD, outFileMap = "a/b.txt"),
             "No permission to write to")
expect_error(writeFlapjack(SxMIBD, outFileMap = mapOut, outFileGeno = 1),
             "outFile should be a single character string ending in .txt")
expect_error(writeFlapjack(SxMIBD, outFileMap = mapOut, outFileGeno = "a/b.txt"),
             "No permission to write to")

## There is very little that can actually be tested.
## Just checking silent execution.
expect_silent(writeFlapjack(SxMIBD, outFileMap = mapOut, outFileGeno = genoOut))
