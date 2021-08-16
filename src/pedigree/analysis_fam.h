#ifndef ANALYSIS_FAMILIES_HEADER
#define ANALYSIS_FAMILIES_HEADER

#include <string>
#include <vector>
#include <map>

#include "matvec.h"
#include "popt.h"
#include "markerscore.h"
#include "OrdGeno.h"
#include "Loc.h"

//! calculation of IBD probabilities in families at QTL positions
/*!
 Function Object for calculation of IBD probabilities at putative QTL positions
 */
class IBD_fam
{
public:
  IBD_fam(const ibd::matrix<ibd::OrdGeno>& P,
          const std::vector<score>& Offspring,
          const LinkageMap& MarkerMap,
          const std::string& poptype);
  ~IBD_fam() { delete popt; }
  std::map<ibd::OrdGeno,double> operator()(const Locus& QTLpos) const;
  std::vector<double> check_scores(const std::vector<ibd::OrdGeno>& geno, const score& sc_off) const;
private:
  pop_base *popt;
  int len_inh; // length of inheritance vector
  std::vector<ibd::OrdGeno> par;
  LinkageMap markermap;
  ibd::matrix<double> l_cond, r_cond;
};

#endif


