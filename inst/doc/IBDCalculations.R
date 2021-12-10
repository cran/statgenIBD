## ----setup, include = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.dim = c(7, 4)
)
library(statgenIBD)
op <- options(width = 90)

## ----inspectMap-------------------------------------------------------------------------
## Read the map and display the first tows.
map <- read.table(system.file("extdata/SxM", "SxM_map.txt", package = "statgenIBD"))
head(map)

## ----inspectLoc-------------------------------------------------------------------------
## Read the genotypic file and display the first rows and columns.
geno <- read.table(system.file("extdata/SxM", "SxM_geno.txt", package = "statgenIBD"),
                   header = TRUE)
head(geno[, 1:5])

## ----sxmIBD-----------------------------------------------------------------------------
## Compute IBD probabilities for Steptoe Morex.
SxMIBD <- calcIBD(popType = "DH",
                  markerFile = system.file("extdata/SxM", "SxM_geno.txt",
                                           package = "statgenIBD"),
                  mapFile = system.file("extdata/SxM", "SxM_map.txt",
                                        package = "statgenIBD"))

## Print summary.
summary(SxMIBD)

## ----check_dh001_plc--------------------------------------------------------------------
SxMIBD$markers["plc", "dh001", ]

## ----check_dh005_abg313b----------------------------------------------------------------
SxMIBD$markers["abg313b", "dh005", ]

## ----plotsxmIBD-------------------------------------------------------------------------
## Visualize IBD probabilities for dh005.
plot(SxMIBD, 
     plotType = "singleGeno",
     genotype = "dh005")

## ----plotsxmIBDAll----------------------------------------------------------------------
## Visualize IBD probabilities for all genotypes.
plot(SxMIBD, 
     plotType = "allGeno",
     genotype = "dh005")

## ----sxmIBD_extGrid---------------------------------------------------------------------
## Compute IBD probabilities for Steptoe Morex.
## Add extra evaluation positions on dense grid.
SxMIBD_Ext_grid <- calcIBD(popType = "DH",
                           markerFile = system.file("extdata/SxM", "SxM_geno.txt",
                                                    package = "statgenIBD"),
                           mapFile = system.file("extdata/SxM", "SxM_map.txt",
                                                 package = "statgenIBD"),
                           evalDist = 1,
                           grid = TRUE)

## Print summary.
summary(SxMIBD_Ext_grid)

## ----SxMIBD_Ext_grid--------------------------------------------------------------------
## Visualize IBD probabilities for dh005
plot(SxMIBD_Ext_grid, 
     genotype = "dh005")

## ----sxmIBD_ext-------------------------------------------------------------------------
## Compute IBD probabilities for Steptoe Morex.
## Add extra evaluation positions between existing markers.
SxMIBD_Ext <- calcIBD(popType = "DH",
                      markerFile = system.file("extdata/SxM", "SxM_geno.txt",
                                               package = "statgenIBD"),
                      mapFile = system.file("extdata/SxM", "SxM_map.txt",
                                            package = "statgenIBD"),
                      evalDist = 5,
                      grid = FALSE)

## Print summary.
summary(SxMIBD_Ext)

## ----sxmIBD_ext_map---------------------------------------------------------------------
## Show first rows of map in output.
head(SxMIBD_Ext$map)

## ----inspectPosFile---------------------------------------------------------------------
## Read the evalPos file and display the first rows.
evalPos <- read.table(system.file("extdata/SxM", "SxM_eval.txt", package = "statgenIBD"), 
                      header = TRUE)
head(evalPos)

## ----sxmIBD_posFile---------------------------------------------------------------------
SxMIBD_evalPos <- calcIBD(popType = "DH",
                          markerFile = system.file("extdata/SxM", "SxM_geno.txt",
                                                   package = "statgenIBD"),
                          mapFile = system.file("extdata/SxM", "SxM_map.txt",
                                                package = "statgenIBD"),
                          evalPos = evalPos)

## Print summary.
summary(SxMIBD_evalPos)

## ----SxMProbs_one-----------------------------------------------------------------------
## Extract marker probabilities for markers plc and ABG053.
SxM_probs <- getProbs(SxMIBD, markers = c("plc", "ABG053"))
head(SxM_probs)

## ----SxMwriteFlapjack, eval=FALSE-------------------------------------------------------
#  ## Write results to Flapjack format.
#  writeFlapjack(SxMIBD_Ext,
#                outFileMap = "map.txt",
#                outFileGeno = "geno.txt")

## ----F4IBD------------------------------------------------------------------------------
## Compute IBD probabilities for simulated F4 population.
F4IBD <- calcIBD(popType = "F4",
                 markerFile = system.file("extdata/popF4", "cross.txt",
                                          package = "statgenIBD"),
                 mapFile = system.file("extdata/popF4", "mapfile.txt",
                                       package = "statgenIBD"))

## Print summary.
summary(F4IBD)

## ----F4Probs----------------------------------------------------------------------------
## Extract marker probabilities for markers M1_2 and M1_4.
F4_probs <- getProbs(F4IBD, markers = c("M1_2", "M1_4"))
head(F4_probs)

## ----F4Probs_sum------------------------------------------------------------------------
## Extract marker probabilities for markers M1_2 and M1_4.
## Sum the probabilities to probabilities per parent.
F4_probs_sum <- getProbs(F4IBD, markers = c("M1_2", "M1_4"), sumProbs = TRUE)
head(F4_probs_sum)

## ----C3S4DHIBD--------------------------------------------------------------------------
## Compute IBD probabilities for simulated C4S3 population.
C4S3IBD <- calcIBD(popType = "C4S3",
                   markerFile = system.file("extdata/popC4S3", "cross.txt",
                                            package = "statgenIBD"),
                   mapFile = system.file("extdata/popC4S3", "mapfile.txt",
                                         package = "statgenIBD"))

## Print summary.
summary(C4S3IBD)

## ----multiIBD---------------------------------------------------------------------------
## Compute IBD probabilties for AxB.
AB <- calcIBD(popType = "F4DH",
              markerFile = system.file("extdata/multipop", "AxB.txt",
                                       package = "statgenIBD"),
              mapFile = system.file("extdata/multipop", "mapfile.txt",
                                    package = "statgenIBD"),
              evalDist = 1)

## Print summary.
summary(AB)

## Compute IBD probabilties for AxC.
AC <- calcIBD(popType = "F4DH",
              markerFile = system.file("extdata/multipop", "AxC.txt",
                                       package = "statgenIBD"),
              mapFile = system.file("extdata/multipop", "mapfile.txt",
                                    package = "statgenIBD"),
              evalDist = 1)

## Print summary.
summary(AC)

## ----combineIBD-------------------------------------------------------------------------
ABC <- c(AB, AC)
summary(ABC)

## ----checkMultiIBD----------------------------------------------------------------------
## Extract probabilities for markers EXT_1_1 and EXT_1_3.
ABCProbs <- getProbs(ABC, markers = c("EXT_1_1", "EXT_1_3"))

## Print probabilities for genotypes AxB0001 and AxC0001.
ABCProbs[ABCProbs$geno %in% c("AxB0001", "AxC0001"), ]

## ----plotMultiIBD, fig.show="hold", out.width="47%"-------------------------------------
plot(ABC, genotype = "AxB0001")
plot(ABC, genotype = "AxC0001")

## ----plotMultiIBDAll--------------------------------------------------------------------
plot(ABC, plotType = "allGeno")

## ----winddown, include = FALSE------------------------------------------------
options(op)

