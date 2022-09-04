### Test readRABBIT

## Define input files.

genoFile <- system.file("extdata/barley", "barley_magicReconstruct.zip",
                        package = "statgenIBD")
pedFile <- system.file("extdata/barley", "barley_pedInfo.csv",
                       package = "statgenIBD")

## Checks for correct input.
expect_error(readRABBIT(infile = 1),
             "infile should be a character string indicating a readable")
expect_error(readRABBIT(infile = "tst"),
             "infile should be a character string indicating a readable")
expect_error(readRABBIT(infile = "tst.csv"),
             "infile should be a character string indicating a readable")

expect_error(readRABBIT(infile = unzip(genoFile, exdir = tempdir()),
                        pedFile = 1),
             "pedFile should be a character string indicating a readable")
expect_error(readRABBIT(infile = unzip(genoFile, exdir = tempdir()),
                        pedFile = "tst"),
             "pedFile should be a character string indicating a readable")
expect_error(readRABBIT(infile = unzip(genoFile, exdir = tempdir()),
                        pedFile = "tst.csv"),
             "pedFile should be a character string indicating a readable")


## Different combinations of inputs should give similar output.
expect_silent(barleyMPP <-
                readRABBIT(infile = unzip(genoFile, exdir = tempdir())))
expect_silent(barleyMPP2 <-
                readRABBIT(infile = unzip(genoFile, exdir = tempdir()),
                           pedFile = pedFile))

## General structure.
expect_inherits(barleyMPP, "IBDprob")

expect_inherits(barleyMPP$map, "data.frame")
expect_equal(dim(barleyMPP$map), c(355, 2))

expect_inherits(barleyMPP$markers, "array")
expect_equal(dim(barleyMPP$markers), c(916, 355, 5))

expect_inherits(attr(barleyMPP2, "genoCross"), "data.frame")
expect_inherits(barleyMPP2$pedigree, "data.frame")

## map and markers should be the same for all.
expect_equal(barleyMPP$map, barleyMPP2$map)
expect_equal(barleyMPP$markers, barleyMPP2$markers)

