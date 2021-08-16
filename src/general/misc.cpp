#include <algorithm>
#include <numeric>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cfloat>
#include <float.h> // for DBL_MIN

#include "misc.h"
#include "ibdexcept.h"

using namespace ibd;
using namespace std;

double ibd::round(double x,
                  int precision)
{
	double fac = pow10(precision);
	return floor(fac*x+0.5)/fac;
}

void ibd::make_conditional(vector<double>& p)
{
	double sum = accumulate(p.begin(),p.end(),0.0);
	if (sum <= DBL_MIN)
		throw ibd_error("function make_conditional()");
	for (vector<double>::iterator it = p.begin(); it != p.end(); ++it)
		*it/=sum;
}

vector<double> ibd::elem_prod(const vector<double>& a, const vector<double>& b)
{
	size_t dim = a.size();
	vector<double> result(dim);
	for (unsigned int i=0;i<dim;i++)
		result[i] = a[i]*b[i];
	return result;
}

istream& ibd::skip_rest_of_line(istream& inp)
{
	inp.ignore(numeric_limits<int>::max(),'\n');
	return inp;
}

istream& ibd::eatcomment(istream& inp)
{
    char c;
    while (inp.get(c))
	{
		if (c == ';' || c == '#')
			skip_rest_of_line(inp);
		else if (!isspace(c))
		{
			inp.putback(c);
			return inp;
		}
	}
	return inp;
}

void ibd::OpenFile(ofstream& outp, string filename)
{
	outp.open(filename.c_str());
	if (!outp)
		throw ibd_error("Cannot open file " + filename);
	outp.setf(ios::fixed, ios::floatfield);
}

void ibd::OpenFile(ifstream& inp, string filename)
{
	inp.open(filename.c_str());
	if (!inp)
		throw ibd_error("Cannot open file " + filename);
}

unsigned int ibd::pow2(int n)
{
	unsigned int x=1;
	return (x << n);
}
