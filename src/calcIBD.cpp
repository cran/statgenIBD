#include <string>
#include <RcppArmadillo.h>

#include "mainR.h"
#include "matvec.h"
#include "Loc.h"
#include "read_map.h"

#include "popt.h"
#include "util_genetics.h"

using namespace Rcpp;
using namespace std;
using namespace ibd;

//' Calculate IBD probabilities
//'
//' Calculate IBD probabilities for different types of populations.
//'
//' IBD probabilities can be calculated for many different types of populations.
//' In the following table all supported populations are listed. Note that the
//' value of x in the population types is variable, with its maximum value
//' depicted in the last column.
//'
//' | __Population type__ | __Cross__ | __Description__ | __max. x__ |
//' | ------ | ----------------- | -------------------------------------- | --- |
//' | DH | biparental | doubled haploid population | |
//' | Fx | biparental | Fx population (F1, followed by x-1 generations of selfing) | 8 |
//' | FxDH | biparental | Fx, followed by DH generation | 8 |
//' | BCx | biparental | backcross, second parent is recurrent parent | 9 |
//' | BCxDH | biparental | BCx, followed by DH generation | 9 |
//' | BC1Sx | biparental | BC1, followed by x generations of selfing | 7 |
//' | BC1SxDH | biparental | BC1, followed by x generations of selfing and DH | 6 |
//' | C3 | three-way | three way cross: (AxB) x C |  |
//' | C3DH | three-way | C3, followed by DH generation |  |
//' | C3Sx | three-way | C3, followed by x generations of selfing | 7 |
//' | C3SxDH | three-way | C3, followed by x generations of selfing and DH generation | 6 |
//' | C4 | four-way | four-way cross: (AxB) x (CxD)	| |
//' | C4DH | four-way | C4, followed by DH generation |  |
//' | C4Sx | four-way | C4, followed by x generations of selfing | 6 |
//' | C4SxDH | four-way | C4, followed by x generations of selfing and DH generation | 6 |
//'
//' @param popType A character string indicating the type of population. One of
//' DH, Fx, FxDH, BCx, BCxDH, BC1Sx, BC1SxDH, C3, C3DH, C3Sx, C3SxDH, C4, C4DH,
//' C4Sx, C4SxDH (see Details).
//' @param markerFile A character string indicating the location of the file with
//' genotypic information for the population. The file should be in
//' tab-delimited format with a header containing marker names.
//' @param mapFile A character string indicating the location of the map file
//' for the population. The file should be in tab-delimited format. It should
//' consist of exactly three columns, marker, chromosome and position. There
//' should be no header. The positions in the file should be in centimorgan.
//' @param evalPos A data.frame with evaluation positions to which the
//' calculations should be limited.
//' @param evalDist An optional numerical value indicating the maximum
//' distance for marker positions. Extra markers will be added based on the
//' value of \code{grid}.
//' @param grid Should the extra markers that are added to assure the a
//' maximum distince of \code{evalDist} be on a grid (\code{TRUE}) or in between
//' marker existing marker positions (\code{FALSE}).
//' @param verbose Should messages indicating the progress of the process be
//' printed?
//'
//' @return An object of class \code{IBDprob}, a \code{list} with five elements,
//' \describe{
//' \item{map}{a \code{data.frame} with chromosome and position of the markers.}
//' \item{markers}{a 3-dimensional \code{array} of IBD probabilities with
//' markers, genotypes and  parents as array dimensions.}
//' \item{parents}{the parents.}
//' \item{popType}{the population type.}
//' \item{multiCross}{a logical value indicating if multiple crosses have been
//' combined in the \code{IBDprob} object.}
//' }
//'
//' @examples
//' ## Compute IBD probabilities for Steptoe Morex.
//' SxMIBD <- calcIBD(popType = "DH",
//'                   markerFile = system.file("extdata/SxM", "SxM_geno.txt",
//'                                         package = "statgenIBD"),
//'                   mapFile = system.file("extdata/SxM", "SxM_map.txt",
//'                                         package = "statgenIBD"))
//'
//' ## Check summary.
//' summary(SxMIBD)
//'
//' ## Compute IBD probabilities for Steptoe Morex.
//' ## Add extra evaluation positions in between exiting marker positions
//' ## to assure evaluation positions are at most 5 cM apart.
//' SxMIBD_Ext <- calcIBD(popType = "DH",
//'                       markerFile = system.file("extdata/SxM", "SxM_geno.txt",
//'                                               package = "statgenIBD"),
//'                       mapFile = system.file("extdata/SxM", "SxM_map.txt",
//'                                             package = "statgenIBD"),
//'                       evalDist = 5)
//'
//' ## Check summary.
//' summary(SxMIBD_Ext)
//'
//' @export
// [[Rcpp::export]]
List calcIBD(CharacterVector& popType,
             CharacterVector& markerFile,
             CharacterVector& mapFile,
             Nullable<DataFrame&> evalPos = R_NilValue,
             Nullable<NumericVector&> evalDist = R_NilValue,
             const bool& grid = true,
             const bool& verbose = false)
{
  string _poptype = Rcpp::as<std::string>(popType);
  // only to check poptype has correct format:
  const pop_base *popt = init_pop(_poptype);

  int x;
  bool isDH = _poptype.find("DH") != std::string::npos;
  bool isBC = match(x, _poptype, "BCx");
  LinkageMap positions;
  arma::cube prob;
  vector<string> parents, offspring;
  double max_step_size = -1;
  if (evalDist.isNotNull())
    max_step_size = Rcpp::as<double>(evalDist);
  try
  {
    main_pedigreeR(prob, parents, offspring, positions,
                   _poptype,
                   Rcpp::as<std::string>(markerFile),
                   Rcpp::as<std::string>(mapFile),
                   evalPos,
                   max_step_size,
                   grid,
                   verbose);
  }
  catch (ibd_error& e)
  {
    forward_exception_to_r(e);
  }
  catch (std::exception& e)
  {
    forward_exception_to_r(e);
  }
  catch (...)
  {
    ::Rf_error("c++ exception (unknown reason)");
  }
  const int npar = parents.size();
  const int M = positions.size();
  unsigned int nSlices = prob.n_slices;
  // Remove slices with only zero
  for (arma::uword i = 0; i < nSlices; i++) {
    if (all(vectorise(prob.slice(nSlices - i - 1)) == 0)) {
      prob.shed_slice(nSlices - i - 1);
    }
  }
  // Construct vector of names for parents.
  CharacterVector parentNames (0);
  for (int i = 0; i < npar; i++)
  {
    if (!(isBC && i == 0))
    {
      parentNames.push_back("p" + parents[i]);
    }
  }
  if (!isDH)
  {
    if (npar == 2)
    {
      parentNames.push_back("p" + parents[0] + parents[1]);
    } else if (npar == 3)
    {
      parentNames.push_back("p" + parents[0] + parents[1]);
      parentNames.push_back("p" + parents[0] + parents[2]);
    } else if (npar == 4) {
      parentNames.push_back("p" + parents[0] + parents[2]);
      parentNames.push_back("p" + parents[0] + parents[3]);
      parentNames.push_back("p" + parents[1] + parents[2]);
      parentNames.push_back("p" + parents[1] + parents[3]);
    }
  }
  // Construct map file from positions.
  CharacterVector posNames = CharacterVector(M);
  IntegerVector chr = IntegerVector(M);
  NumericVector pos = NumericVector(M);
  for (int m = 0; m < M; m++) {
    posNames(m) = positions[m].GetName();
    chr(m) = positions[m].GetChr();
    pos(m) = positions[m].GetPosition();
  }
  DataFrame map = DataFrame::create(Named("chr") = chr, Named("pos") = pos);
  map.attr("row.names") = posNames;
  // Reshape prob to 3D array and add names to dimensions.
  NumericVector P = wrap(prob);
  P.attr("dimnames") = List::create(posNames, offspring, parentNames);
  // Create result list: map + markers.
  List res = List::create(Named("map") = map,
                          Named("markers") = P,
                          Named("popType") = popType,
                          Named("parents") = parents,
                          Named("multiCross") = false);
  res.attr("class") = "IBDprob";
  return res;
}
