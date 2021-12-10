### Test IBDprob

## Test summary.

## Define file locations.
SxMmarkers <- system.file("extdata/SxM", "SxM_geno.txt", package = "statgenIBD")
SxMmap <- system.file("extdata/SxM", "SxM_map.txt", package = "statgenIBD")

SxMIBD <- calcIBD(popType = "DH", markerFile = SxMmarkers, mapFile = SxMmap)

SxMSumm <- capture.output(summary(SxMIBD))
expect_true(any(grepl("population type:  DH", SxMSumm)))
expect_true(any(grepl("Number of evaluation points:  116", SxMSumm)))
expect_true(any(grepl("Number of individuals:  150", SxMSumm)))
expect_true(any(grepl("Parents:  Morex Steptoe", SxMSumm)))

## Test plot for single genotype.
expect_error(plot(SxMIBD, genotype = 1),
             "should be a character string")
expect_error(plot(SxMIBD, genotype = "a"),
             "genotype a not defined")
expect_error(plot(SxMIBD, genotype = "dh001", title = 1),
             "title should be a character string")

expect_silent(p <- plot(SxMIBD, genotype = "dh001"))

expect_inherits(p, "ggplot")

## Check that option title functions correctly.

p1 <- plot(SxMIBD, genotype = "dh001", title = "tst")
expect_equal(p1$labels$title, "tst")

## Test plot for all genotypes.

expect_silent(p2 <- plot(SxMIBD, plotType = "allGeno"))

expect_inherits(p2, "ggplot")
