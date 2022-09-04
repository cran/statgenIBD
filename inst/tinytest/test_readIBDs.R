### Test readIBDs function.

## Define input file.
ibdFile <- system.file("extdata/SxM", "SxM_IBDs.txt", package = "statgenIBD")
mapFile <- system.file("extdata/SxM", "SxM_map.txt", package = "statgenIBD")
map <- read.delim(mapFile, header = FALSE)
rownames(map) <- map$V1
map <- map[,-1]
colnames(map) <- c("chr", "pos")

map2 <- map3 <- map
map2[["pos"]] <- as.character(map2[["pos"]])
rownames(map3)[1] <- "a"

## Checks for correct input.
expect_error(readIBDs(infile = 1),
             "should be a character string indicating a readable .txt or .ibd file")
expect_error(readIBDs(infile = "a/b.txt"),
             "should be a character string indicating a readable .txt or .ibd file")
expect_error(readIBDs(infile = ibdFile, map = 1),
             "map should be a data.frame")
expect_error(readIBDs(infile = ibdFile, map = map[, 1, drop = FALSE]),
             "chr and pos should be columns in map")
expect_error(readIBDs(infile = ibdFile, map = map2),
             "pos should be a numeric column in map")

expect_warning(readIBDs(infile = ibdFile, map = map[-1, ]),
               "1 markers in infile are not in map")
expect_warning(readIBDs(infile = ibdFile, map = map3),
               "1 markers in map are not in infile")

## Check for successful read and returned structure.
expect_silent(ibds <- readIBDs(infile = ibdFile, map = map))
expect_inherits(ibds, "IBDprob")
expect_inherits(ibds$markers, "array")
expect_equal(dim(ibds$markers), c(150, 116, 2))
expect_equal(length(ibds$parents), 2)
