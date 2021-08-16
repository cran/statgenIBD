#include <RcppArmadillo.h>

#include "HMMalgo.h"

using namespace std;
using namespace ibd;

vector<double> ibd::forward_equation(const vector<double>& p_prev,
                                     const TransMatSym2D& T,
                                     const vector<double>& q)
{
  vector<double> p = p_prev*T;
  size_t N = p_prev.size();
  for (unsigned int i=0;i<N;i++)
    p[i] *= q[i];
  make_conditional(p);
  return p;
}

vector<double> ibd::backward_equation(const vector<double>& p_next,
                                      const TransMatSym2D& T,
                                      const vector<double>& q)
{
  vector<double> p = T*p_next;
  size_t N = p_next.size();
  for (unsigned int i=0;i<N;i++)
    p[i] *= q[i];
  make_conditional(p);
  return p;
}

matrix<double> ibd::calc_prob_left(const std::vector<double>& pi0,
                                   const matrix<double>& q,
                                   const vector<TransMatSym2D>& T)
{
  int M = q.NrRows();
  int N = q.NrCols();
  matrix<double> L(M,N,0.0);

  for (int i=0;i<N;i++)
    L[0][i] = pi0[i]*q[0][i];
  make_conditional(L[0]);

  for (int loc=1;loc<M;loc++)
    L[loc] = forward_equation(L[loc-1],T[loc-1],q[loc]);
  return L;
}

matrix<double> ibd::calc_prob_right(const matrix<double>& q,
                                    const vector<TransMatSym2D>& T)
{
  int M = q.NrRows();
  int N = q.NrCols();
  matrix<double> R(M,N,0.0);

  R[M-1] = q[M-1];
  make_conditional(R[M-1]);

  for (int loc=M-2;loc>=0;loc--)
    R[loc] = backward_equation(R[loc+1],T[loc],q[loc]);
  return R;
}
