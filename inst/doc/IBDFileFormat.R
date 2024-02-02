## ----setup, include=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.dim = c(7, 4)
)
library(statgenIBD)
op <- options(width = 90)

## ----SxMwriteIBD------------------------------------------------------------------------
## Compute IBD probabilities for Steptoe Morex.
SxMIBD <- calcIBD(popType = "DH",
                  markerFile = system.file("extdata/SxM", "SxM_geno.txt",
                                           package = "statgenIBD"),
                  mapFile = system.file("extdata/SxM", "SxM_map.txt",
                                        package = "statgenIBD"))

## Write IBDs to tab-delimited .txt file.
writeIBDs(SxMIBD, "SxM-IBD.txt")

## ----echo=FALSE-------------------------------------------------------------------------
ibdFile <- system.file("extdata/SxM", "SxM_IBDs.txt", package = "statgenIBD")
ibds <- read.delim(ibdFile)
knitr::kable(head(ibds))

## ----eval=FALSE-------------------------------------------------------------------------
#  ## Write IBDs to file, set values <0.05 to zero and only print 3 decimals.
#  writeIBDs(IBDprob = SxMIBD, outFile = tempfile(fileext = ".txt"),
#            decimals = 3, minProb = 0.05)

## ----SxMreadIBD-------------------------------------------------------------------------
## Get map.
SxMMap <- SxMIBD$map

## Read IBDs to tab-delimited .txt file.
SxMIBD <- readIBDs("SxM-IBD.txt", map = SxMMap)
summary(SxMIBD)

## ----echo=FALSE, results='hide'---------------------------------------------------------
unlink("SxM-IBD.txt")

