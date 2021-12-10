#ifndef READ_MAP_PEDIGREE_HEADER
#define READ_MAP_PEDIGREE_HEADER

#include <string>
#include <vector>

// library file
#include "Loc.h"

// output: LinkageMap with only the markers which are in the loc file
LinkageMap reduce_markermap(const LinkageMap& markermap, const std::vector<std::string>& markers);

// if distance between markers is less than eps=1.0e-6 then distance is
// assumed to be delta = 0.001
LinkageMap adjust_markermap(const LinkageMap& markermap);

LinkageMap read_map_file(const std::string& mapfile);

void print_marker_warnings(const LinkageMap& markermap, const std::vector<std::string>& markers);

LinkageMap read_eval_pos_df(const Rcpp::DataFrame& evalposdf);

LinkageMap select_chr(const LinkageMap& markermap, std::string sel_chr);

#endif

