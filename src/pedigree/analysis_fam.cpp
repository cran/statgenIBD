#include "misc.h"
#include "HMMalgo.h"
#include "TransMatSym2D.h"
#include "InhVector.h"
#include "analysis_fam.h"

using namespace ibd;
using namespace std;

IBD_fam::IBD_fam(const ibd::matrix<OrdGeno>& P,
                 const std::vector<score>& offspring,
                 const LinkageMap& MarkerMap,
                 const std::string& poptype)
  : markermap(MarkerMap)
{
  popt = init_pop(poptype);
  len_inh = popt->get_len();

  const unsigned int npar = P.NrRows();
  const unsigned int M = markermap.size();
  const unsigned int N = pow2(len_inh); // dimension of state space

  for (unsigned int i=0;i<npar;i++)
    par.push_back(OrdGeno(i,i));

  // init pi0
  vector<double> pi0(N,1.0/N);

  // init T
  vector<double> r = make_rec_map(markermap);

  int Nintervals = r.size();
  vector<TransMatSym2D> T(Nintervals); // Transition matrices between the loci.
  for (int k=0;k<Nintervals;k++)
    T[k] = TransMatSym2D(len_inh,r[k]);

  // init Q
  matrix<double> Q(M,N);
  for (unsigned int m=0;m<M;m++)
  {
    vector<OrdGeno> geno(npar);
    for (unsigned int par=0;par<npar;par++)
      geno[par] = P[par][m];
    Q[m] = check_scores(geno,offspring[m]);
  }

  // calculate left- and right- conditional probabilities.
  l_cond = calc_prob_left(pi0,Q,T);
  r_cond = calc_prob_right(Q,T);
}

std::vector<double> IBD_fam::check_scores(const std::vector<OrdGeno>& geno,
                                          const score& sc_off) const
{
  bool all_inconsistent = true;
  const unsigned int N = pow2(len_inh); // dimension of state space
  vector<double> q(N);
  for (InhVector u(len_inh);!u.end();u++)
  {
    OrdGeno g = popt->gen_off(geno,u);
    if (check_score(g,sc_off))
    {
      all_inconsistent = false;
      q[u] = 1.0;
    }
    else
      q[u] = 0.0;
  }
  if (all_inconsistent)
    return vector<double>(N,1.0);
  else
    return q;
}

std::map<OrdGeno,double> IBD_fam::operator()(const Locus& QTLpos) const
{
  const int left = pos_qtl(markermap,QTLpos);
  const int right = left+1;
  double r_left = recomb(markermap[left],QTLpos);
  double r_right = recomb(QTLpos,markermap[right]);
  const vector<double>& L = l_cond[left];
  const vector<double>& R = r_cond[right];
  TransMatSym2D T1(len_inh,r_left);
  TransMatSym2D T2(len_inh,r_right);
  vector<double> p = elem_prod(L*T1,T2*R);
  make_conditional(p);

  map<OrdGeno,double> IBDprob;
  for (InhVector u(len_inh);!u.end();u++)
  {
    OrdGeno g = popt->gen_off(par,u);
    IBDprob[g] += p[u];
  }
  return IBDprob;
}





