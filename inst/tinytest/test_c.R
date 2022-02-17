### Test IBDprob concatenation.

## Define file locations.
ABmarkers <- system.file("extdata/multipop", "AxB.txt", package = "statgenIBD")
ACmarkers <- system.file("extdata/multipop", "AxC.txt", package = "statgenIBD")
ABCmap <- system.file("extdata/multipop", "mapfile.txt", package = "statgenIBD")

## IBD calculations for two populations separately.
AB <- calcIBD(popType = "F4DH", markerFile = ABmarkers, mapFile = ABCmap)
AC <- calcIBD(popType = "F4DH", markerFile = ACmarkers, mapFile = ABCmap)

## Alternative calculations for AC to test input checks.
ACalt1 <- calcIBD(popType = "F4", markerFile = ACmarkers, mapFile = ABCmap)
ACalt2 <- calcIBD(popType = "F4DH", markerFile = ACmarkers, mapFile = ABCmap,
                  evalDist = 5)

## Check that input checks are working correctly.
expect_error(c(AB, "tst"),
             "All inputs should be of class IBDprob")
expect_error(c(AB, ACalt1),
             "All inputs should have the same population type")
expect_error(c(AB, ACalt2),
             "All inputs should have the same map")

## c with single argument should return input.
expect_equal(c(AB), AB)

## Multiple outputs should be combined correctly.
ABC <- c(AB, AC)

expect_inherits(ABC, "IBDprob")
expect_equal(ABC$map, AB$map)
expect_equal_to_reference(ABC$markers, "ABC_markers")
expect_equal(ABC$parents, c("A", "B", "C"))
expect_equal(ABC$popType, AB$popType)
expect_true(ABC$multiCross)

## Check that genoCross attribute is added correctly.
genoCross <- attr(x = ABC, which = "genoCross")

expect_inherits(genoCross, "data.frame")
expect_equal(genoCross[["cross"]],
             rep(c("cross1", "cross2"), times = c(100, 80)))
expect_equal(genoCross[["geno"]],
             c(dimnames(AB$markers)[[2]], dimnames(AC$markers)[[2]]))






