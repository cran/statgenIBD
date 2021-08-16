### Test getProbs function.

## Define file locations.
ABmarkers <- system.file("extdata/multipop", "AxB.txt", package = "statgenIBD")
ACmarkers <- system.file("extdata/multipop", "AxC.txt", package = "statgenIBD")
ABCmap <- system.file("extdata/multipop", "mapfile.txt", package = "statgenIBD")

## IBD calculations for two populations separately.
AB <- calcIBD(popType = "F4DH", markerFile = ABmarkers, mapFile = ABCmap)
AC <- calcIBD(popType = "F4DH", markerFile = ACmarkers, mapFile = ABCmap)
ABC <- c(AB, AC)

## Check that input checks are working correctly.
expect_error(getProbs(IBDprob = "tst", markers = "M1_1"),
             "should be an object of class IBDprob")
expect_error(getProbs(IBDprob = AB, markers = 1),
             "markers should be a character vector")
expect_error(getProbs(IBDprob = AB, markers = "M1"),
             "The following markers are not in AB: M1")

## Check that output is correct.
AB_M1_1 <- getProbs(IBDprob = AB, markers = "M1_1")

expect_inherits(AB_M1_1, "data.frame")
expect_equal_to_reference(AB_M1_1, "AB_M1_1")

AB_M1_1_M3_3 <- getProbs(IBDprob = AB, markers = c("M1_1", "M3_3"))

expect_equal(AB_M1_1, AB_M1_1_M3_3[, 1:3])

## For multicross there should be an extra column cross.

ABC_M1_1 <- getProbs(IBDprob = ABC, markers = "M1_1")
expect_equal(colnames(ABC_M1_1),
             c("cross", "geno", "M1_1_A", "M1_1_B", "M1_1_C"))

## Check that option sumProbs functions correctly.

## Define file locations.
F4markers <- system.file("extdata/popF4", "cross.txt", package = "statgenIBD")
F4map <- system.file("extdata/popF4", "mapfile.txt", package = "statgenIBD")

F4 <- calcIBD(popType = "F4", markerFile = F4markers, mapFile = F4map)

F4_M1_1a <- getProbs(IBDprob = F4, markers = "M1_1", sumProbs = FALSE)
F4_M1_1b <- getProbs(IBDprob = F4, markers = "M1_1", sumProbs = TRUE)

expect_equal(F4_M1_1b[["M1_1_A"]],
             F4_M1_1a[["M1_1_A"]] + 0.5 * F4_M1_1a[["M1_1_AB"]])
expect_equal(F4_M1_1b[["M1_1_B"]],
             F4_M1_1a[["M1_1_B"]] + 0.5 * F4_M1_1a[["M1_1_AB"]])

