#ifndef TRANSMATRIX_SYMMETRIC_2D_HEADER
#define TRANSMATRIX_SYMMETRIC_2D_HEADER

#include "misc.h"

namespace ibd
{

class TransMatSym2D : public std::vector<double>
{
public:
	TransMatSym2D() {}
	TransMatSym2D(const std::vector<double>& R) : std::vector<double>(R) {}
	TransMatSym2D(int N, double r) : std::vector<double>(N,r) { }

	unsigned int Dimension() const { return pow2(this->size()); }
};

std::vector<double> operator*(const TransMatSym2D& A, const std::vector<double>& x);
std::vector<double> operator*(const std::vector<double>& x,const TransMatSym2D& A);

}

#endif



