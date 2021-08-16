#include "ibdexcept.h"
#include "TransMatSym2D.h"

using namespace std;
using namespace ibd;

// Efficient calculation of y = A*x,
// where A is the kronecker product of n symmetric 2x2 matrices B[k] of the form
//
// B[k] =  ( 1-r[k]   r[k]   )  =  ( s[k]  r[k] )
//         (  r[k]    1-r[k] )  =  ( r[k]  s[k] )  , with s[k] = 1-r[k].
//
// The calculations use a recursive formula, based on divide and conquer:
//
// Let C[k] be the kronecker product of B[k],B[k-1],...,B[2],B[1], which
// implicates that C[n] = A, and C[1] = B[1].
//
// Using the definition of the Kronecker product we can derive that:
//
//   y = A*x = C[n]*x --->
//
//  ( y[0] )   ( s[n]*C[n-1]  r[n]*C[n-1] )   ( x[0] )
//  (      ) = (                          ) * (      )   --->
//  ( y[1] )   ( r[n]*C[n-1]  s[n]*C[n-1] )	  ( x[1] )
//
//  y[0] =  s[n]*u[0] + r[n]*u[1]
//  y[1] =  r[n]*u[0] + s[n]*u[1]
//
//  with u[i] = C[n-1]*x[i]  (i = 0,1)
//
//  Note that the dimension of the vectors x and y is 2^n, the dimension
//  of the vectors x[i], y[i], and u[i] is 2^(n-1) .
//
vector<double> product(int k,
                       unsigned int p,
                       const TransMatSym2D& A,
                       const vector<double>& x)
{
  if (k<0)
    return vector<double>(1,x[p]);

  double r = A[k];
  double s = 1.0-r;
  unsigned int M = pow2(k);
  vector<double> u0 = product(k-1,p  ,A,x); // u0 = C[k-1]*x0 (see comments above)
  vector<double> u1 = product(k-1,p+M,A,x); // u1 = C[k-1]*x1 (see comments above)
  vector<double> y(2*M);
  for (unsigned int i=0;i<M;i++)
  {
    y[i]   = s*u0[i] + r*u1[i];	// y0 = s*u0 + r*u1 (see comments above)
    y[i+M] = r*u0[i] + s*u1[i]; // y1 = r*u0 + s*u1 (see comments above)
  }
  return y;
}

// calculates y = A*x
// where A is the kronecker product of N symmetric 2x2 matrices B[k] of the form
//
// B[k] =  ( 1-r[k]   r[k]   )
//         (  r[k]    1-r[k] )
//
vector<double> ibd::operator*(const TransMatSym2D& A,
                              const vector<double>& x)
{
  if (x.size() != A.Dimension())
    throw ibd_error("Error in multiplication A*y");
  return product(A.size()-1,0,A,x);
}

// calculates y^T = x^T*A
// where A is the kronecker product of N symmmetric 2x2 matrices of the form
//
// B[k] =  ( 1-r[k]   r[k]   )
//         (  r[k]    1-r[k] )
//
// Note that the matrix A is symmetric because:
//   kroneck(B[1],...,B[N])^T = kroneck(B[1]^T,...,B[N]^T)= kroneck(B[1],...,B[N])
// This implicates that:
// y^T = x^T*A -> y = (x^T*A)^T = A^T * x = A*x
vector<double> ibd::operator*(const vector<double>& x,
                              const TransMatSym2D& A)
{
  if (x.size() != A.Dimension())
    throw ibd_error("Error in multiplication y*A");
  return product(A.size()-1,0,A,x);
}


