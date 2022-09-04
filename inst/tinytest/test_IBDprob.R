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

## Test pedigree plot.

expect_silent(p3 <- plot(SxMIBD, plotType = "pedigree"))

expect_inherits(p3, "ggplot")

## Test plot for mean probabilities

expect_silent(p4 <- plot(SxMIBD, plotType = "meanProbs"))

expect_inherits(p4, "ggplot")

## Check that option chr functions correctly.

expect_error(plot(SxMIBD, plotType = "meanProbs", chr = "a"),
             "chr not found in map")

expect_silent(p5 <- plot(SxMIBD, plotType = "meanProbs", chr = 3))

expect_equal(nrow(p5$data), 28)

## Test plot for total coverage

expect_silent(p6 <- plot(SxMIBD, plotType = "totalCoverage"))

expect_inherits(p6, "ggplot")

## Check that option chr functions correctly.

expect_error(plot(SxMIBD, plotType = "totalCoverage", chr = "a"),
             "chr not found in map")

expect_silent(p7 <- plot(SxMIBD, plotType = "totalCoverage", chr = 3))

expect_equal(p7$data[["coverage"]], c(46.8607505542283, 53.1392494457717))

