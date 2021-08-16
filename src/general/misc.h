#ifndef MISCELLANEOUS_HEADER
#define MISCELLANEOUS_HEADER

#include <cmath>
#include <string>
#include <vector>
#include <map>

#include "ibdexcept.h"

namespace ibd
{

class Interval
{
public:
	Interval() {}
	Interval(double l, double r) : left(l), right(r) {}
	double GetLeft() const { return left; }
	double GetRight() const { return right; }
	double Length() const { return right-left;}
private:
	double left,right;
};

unsigned int pow2(int n);

inline double sqr(double x) { return x*x; }
inline double cub(double x) { return x*x*x; }
inline double pow10(double x) { return pow(10.0,x); }

double round(double x, int precision = 0);

void make_conditional(std::vector<double>& p);
std::vector<double> elem_prod(const std::vector<double>& a,
							  const std::vector<double>& b);

std::istream& skip_rest_of_line(std::istream& inp);
std::istream& eatcomment(std::istream& inp);

void OpenFile(std::ofstream& outp, std::string filename);
void OpenFile(std::ifstream& inp, std::string filename);

template <class T>
std::map<T,int> make_index(const std::vector<T>& x)
{
	std::map<T,int> M;
	const int N = x.size();
	for (int i=0;i<N;i++)
		M[x[i]] = i;
	return M;
}

}

#endif

