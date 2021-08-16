#ifndef HEADER_MAIN_PEDIGREE
#define HEADER_MAIN_PEDIGREE

#include <string>
#include <RcppArmadillo.h>

#include "matvec.h"
#include "Loc.h"

const std::string version = "2.75";
const std::string date    = "april 13, 2020";

int main_pedigreeR(arma::cube& Z,
                   std::vector<std::string>& parents,
                   std::vector<std::string>& offspring,
                   LinkageMap& positions,
                   const std::string& poptype,
                   const std::string& locfile,
                   const std::string& mapfile,
                   const Rcpp::DataFrame& eval_pos_df,
                   const double& max_step_size,
                   const bool& grid,
                   const bool& verbose);

#endif

