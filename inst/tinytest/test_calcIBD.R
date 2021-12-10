### Test calcIBD function.

## Define file locations.
SxMmarkers <- system.file("extdata/SxM", "SxM_geno.txt", package = "statgenIBD")
SxMmap <- system.file("extdata/SxM", "SxM_map.txt", package = "statgenIBD")

## Check the input check are working correctly.

expect_error(calcIBD(popType = "tst", markerFile = SxMmarkers, mapFile = SxMmap),
             "unknown type tst")
expect_error(calcIBD(popType = "DH", markerFile = "tst", mapFile = SxMmap),
             "Cannot open file tst")
expect_error(calcIBD(popType = "DH", markerFile = SxMmarkers, mapFile = "tst"),
             "Cannot read file tst")

## Check that the output structure is correct.
expect_silent(SxMIBD <- calcIBD(popType = "DH", markerFile = SxMmarkers,
                                mapFile = SxMmap))

expect_inherits(SxMIBD, "IBDprob")
expect_equal(names(SxMIBD),
             c("map", "markers", "popType", "parents", "multiCross"))
expect_inherits(SxMIBD$map, "data.frame")
expect_inherits(SxMIBD$markers, "array")
expect_equal(dim(SxMIBD$markers), c(116, 150, 2))
expect_inherits(SxMIBD$popType, "character")
expect_inherits(SxMIBD$parents, "character")
expect_inherits(SxMIBD$multiCross, "logical")

## Check that output content is correct.

expect_equal_to_reference(SxMIBD$map, "SxMIBD_map")
expect_equal_to_reference(SxMIBD$markers, "SxMIBD_markers", tolerance = 10e-6)
expect_equal(SxMIBD$popType, "DH")
expect_false(SxMIBD$multiCross)

## Check that option verbose works correctly.

outMsg <- capture.output(msg0 <- calcIBD(popType = "DH", markerFile = SxMmarkers,
                                         mapFile = SxMmap))
outMsg2 <- capture.output(msg1 <- calcIBD(popType = "DH", markerFile = SxMmarkers,
                                          mapFile = SxMmap, verbose = TRUE))
expect_equal(outMsg, character())
expect_true(any(grepl(pattern = "reading data", x = outMsg2)))
expect_true(any(grepl(pattern = "analysis of family", x = outMsg2)))

## Check that option evalPos works correctly.
evalPos <- read.table(system.file("extdata/SxM", "SxM_eval.txt",
                                  package = "statgenIBD"), header = TRUE)

expect_silent(SxMIBD_evalPos <- calcIBD(popType = "DH", markerFile = SxMmarkers,
                                        mapFile = SxMmap, evalPos = evalPos))

expect_equal_to_reference(SxMIBD_evalPos, "SxMIBD_evalPos")

## Check that calcIBD works correctly with character values chromosomes.
evalPos[["chr"]][4:9] <- rep(c("1b", "1a"), each = 3)

expect_silent(SxMIBD_evalPosChr <-
                calcIBD(popType = "DH", markerFile = SxMmarkers,
                        mapFile = SxMmap, evalPos = evalPos))

expect_equal_to_reference(SxMIBD_evalPosChr, "SxMIBD_evalPosChr")


