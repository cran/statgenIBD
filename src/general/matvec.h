#ifndef MATVEC_HEADER
#define MATVEC_HEADER

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "ibdexcept.h"

namespace ibd
{

template<class T>
std::vector<T> init_vector(const T& x)
{
	std::vector<T> result(1,x);
	return result;
}

template<class T>
std::vector<T> init_vector(const T& x1, const T& x2)
{
	std::vector<T> result(2);
	result[0] = x1;
	result[1] = x2;
	return result;
}

template<class T>
std::vector<T> init_vector(const T& x1, const T& x2, const T& x3)
{
	std::vector<T> result(3);
	result[0] = x1;
	result[1] = x2;
	result[2] = x3;
	return result;
}

template<class T>
std::vector<T> init_vector(const T& x1, const T& x2, const T& x3, const T& x4)
{
	std::vector<T> result(4);
	result[0] = x1;
	result[1] = x2;
	result[2] = x3;
	result[3] = x4;
	return result;
}

template<class T1, class T2>
double inner_product(const std::vector<T1>& a, const std::vector<T2>& b)
{
	typedef std::vector<T1> vec_t1;
	typedef std::vector<T2> vec_t2;

	typename vec_t1::const_iterator it1;
	typename vec_t2::const_iterator it2;

	double sum = 0.0;
	for (it1 = a.begin(), it2 = b.begin(); it1 != a.end(); it1++, it2++)
		sum += (*it1)*(*it2);
	return sum;
}

template<class T>
class matrix  : public std::vector< std::vector<T> >
{
private:
	typedef std::vector< std::vector<T> > mat_t;
public:
	matrix() {}
	matrix(int r, int c) : std::vector< std::vector<T> >(r,std::vector<T>(c)){}
	matrix(int r, int c, const T& init) :
		std::vector< std::vector<T> >(r,std::vector<T>(c,init)) {}
	typename mat_t::size_type NrRows() const { return this->size(); }
	typename mat_t::size_type NrCols() const
		{ return (this->size() == 0) ? 0 : this->front().size();}
	void resize(int nr, int nc, const T& val = T());
};

template <class T>
void matrix<T>::resize(int nr, int nc, const T& val)
{
	if (this->size() < nr)
		insert(this->end(), nr - this->size(), std::vector<T>());
    else if (nr < this->size())
		erase(this->begin() + nr, this->end());

	typedef std::vector< std::vector<T> > vec_vec_t;
	for (typename vec_vec_t::iterator it = this->begin(); it != this->end(); ++it)
		it->resize(nc,val);
}


template<class T>
void chk_equal_size(const std::vector<T>& c1,
					const std::vector<T>& c2,
					const std::string& str)
{
	if (c1.size() != c2.size())
		throw ibd_error("unequal size of vectors in " + str);
}

template<class T>
std::vector<T>& operator-=(std::vector<T>& x1, const std::vector<T>& x2)
{
	typedef std::vector<T> vec_t;
	typename vec_t::iterator it1;
	typename vec_t::const_iterator it2;
	chk_equal_size(x1,x2,"operator-=");

	it1 = x1.begin();
	it2 = x2.begin();
	for (; it1 != x1.end(); it1++, it2++)
		*it1 -= *it2;
	return x1;
}

template<class T>
std::vector<T>& operator+=(std::vector<T>& x1, const std::vector<T>& x2)
{
	typedef std::vector<T> vec_t;
	typename vec_t::iterator it1;
	typename vec_t::const_iterator it2;

	chk_equal_size(x1,x2,"operator+=");
	it1 = x1.begin();
	it2 = x2.begin();
	for (; it1 != x1.end(); it1++, it2++)
		*it1 += *it2;
	return x1;
}

template<class T>
std::vector<T>& operator*=(std::vector<T>& x, const T& lambda)
{
	typedef std::vector<T> vec_t;
	typename vec_t::iterator it;
	for (it = x.begin(); it != x.end(); it++)
		*it *= lambda;
	return x;
}

template<class T>
std::vector<T>& operator/=(std::vector<T>& x, const T& lambda)
{
	typedef std::vector<T> vec_t;
	typename vec_t::iterator it;
	for (it = x.begin(); it != x.end(); it++)
		*it /= lambda;
	return x;
}

template<class T>
std::vector<T> operator-(const std::vector<T>& x1, const std::vector<T>& x2)
{
	chk_equal_size(x1,x2,"operator-");
	std::vector<T> y = x1;
	y-=x2;
	return y;
}

template<class T>
std::vector<T> operator+(const std::vector<T>& x1, const std::vector<T>& x2)
{
	chk_equal_size(x1,x2,"operator+");
	std::vector<T> y = x1;
	y+=x2;
	return y;
}

template<class T>
std::vector<T> operator*(const std::vector<T>& x, const T& lambda)
{
	std::vector<T> y = x;
	y*=lambda;
	return y;
}

template<class T>
std::vector<T> operator*(const T& lambda,const std::vector<T>& x)
{
	std::vector<T> y = x;
	y*=lambda;
	return y;
}

template<class T>
std::vector<T> operator/(const std::vector<T>& x, const T& lambda)
{
	std::vector<T> y = x;
	y/=lambda;
	return y;
}

template<class T>
std::vector<T> operator- (const std::vector<T>& a)
{
    typedef std::vector<T> vec_t;
    typename vec_t::iterator iter_c;
    std::vector<T> c = a;
    iter_c = c.begin();

    for(; iter_c != c.end(); ++iter_c)
		*iter_c = -(*iter_c);
    return c;
}

template<class T>
matrix<T> operator*(const matrix<T>& A, const matrix<T>& B)
{
	const int A_col = A.NrCols();
	const int A_row = A.NrRows();
	const int B_col = B.NrCols();
	matrix<T> C(A_row,B_col);
    for (int i=0;i<A_row;i++)
		for (int j=0;j<B_col;j++)
		{
			T sum = 0;
			for (int k=0;k<A_col;k++)
				sum += A[i][k]*B[k][j];
			C[i][j] = sum;
		}
	return C;
}

template<class T>
matrix<T> operator+(const matrix<T>& A, const matrix<T>& B)
{
	const int ncol = A.NrCols();
	const int nrow = A.NrRows();
	matrix<T> C(nrow,ncol);
    for (int i=0;i<nrow;i++)
		for (int j=0;j<ncol;j++)
			C[i][j] = A[i][j] + B[i][j];
	return C;
}

template<class T>
matrix<T> operator-(const matrix<T>& A, const matrix<T>& B)
{
	const int ncol = A.NrCols();
	const int nrow = A.NrRows();
	matrix<T> C(nrow,ncol);
    for (int i=0;i<nrow;i++)
		for (int j=0;j<ncol;j++)
			C[i][j] = A[i][j] - B[i][j];
	return C;
}

template<class T>
matrix<T> operator*(double lambda, const matrix<T>& A)
{
	const int nr = A.NrRows();
	const int nc = A.NrCols();
	matrix<double> B(nr,nc);
	for (int r=0;r<nr;r++)
		for (int c=0;c<nc;c++)
			B[r][c] = lambda*A[r][c];

	return B;
}

template<class T>
std::vector<T> operator*(const matrix<T>& A, const std::vector<T>& b)
{
	const int dim = A.NrRows();
	std::vector<T> c(dim);
    for (int i=0;i<dim;i++)
		c[i] = inner_product(A[i],b);
	return c;
}

template<class T>
std::ostream& operator<<(std::ostream& outp, const std::vector<T>& v)
{
	int width = outp.width();
	outp.width(0);
	typedef std::vector<T> vec_t;
	for (typename vec_t::const_iterator iter = v.begin(); iter != v.end(); iter++)
		outp << std::setw(width) << *iter;
	return outp;
}

template<class T>
std::ostream& operator<<(std::ostream& outp, const matrix<T>& m)
{
	int width = outp.width();
	outp.width(0);
	typedef matrix<T> mat_t;
	for (typename mat_t::const_iterator iter = m.begin(); iter != m.end(); iter++)
		outp << std::setw(width) << *iter << std::endl;
	return outp;
}

template<class T>
std::istream& operator>>(std::istream& inp, std::vector<T>& v)
{
	if (v.empty())
		throw ibd_error("operator>> vector: empty vector");
	typedef std::vector<T> vec_t;
	for (typename vec_t::iterator iter = v.begin(); iter != v.end(); iter++)
	{
		if (inp.eof())
			throw ibd_error("operator>> vector: end of file");
		inp >> *iter;
		if (inp.fail())
			throw ibd_error("operator>> vector: format error");
	}
	return inp;
}

template<class T>
std::istream& operator>>(std::istream& inp, matrix<T>& m)
{
	if (m.NrCols() == 0 || m.NrRows() == 0)
		throw ibd_error("operator >> matrix: empty matrix");
	typedef matrix<T> mat_t;
	for (typename mat_t::iterator iter = m.begin(); iter != m.end(); iter++)
		inp >> *iter;
	return inp;
}

template<class T>
matrix<T> transpose(const matrix<T>& X)
{
	const int nr = X.NrRows();
	const int nc = X.NrCols();
	matrix<T> Xt(nc,nr);
	for (int r=0;r<nr;r++)
	{
		const std::vector<T>& Xr = X[r];
		for (int c=0;c<nc;c++)
	  		Xt[c][r] = Xr[c];
	}
	return Xt;
}

template<class T>
std::vector<T> column(const matrix<T>& X, int c)
{
	const int nr = X.NrRows();
	std::vector<T> col(nr);
	for (int r=0;r<nr;r++)
		col[r] = X[r][c];
	return col;
}

template<class T>
void set_column(matrix<T>& X, const std::vector<T>& a, int c)
{
	const int nr = X.NrRows();
	for (int r=0;r<nr;r++)
		X[r][c] = a[r];
}


template<class T>
T * conversion_vector(const std::vector<T>& x)
{
	int n = x.size();
	T *ptr = new T [n];
	for (int i=0;i<n;i++)
		ptr[i] = x[i];
	return ptr;
}

template<class T>
void freevector(T *ptr)
{
	delete[] ptr;
}


template<class T>
T ** make_C_matrix(int nr, int nc)
{
	T **ptr = new T *[nr];
	T *block = new T[nr*nc];
	for (int i=0;i<nr;i++)
		*(ptr+i) = block+i*nc;
	return ptr;
}

template<class T>
T ** conversion_matrix(const matrix<T>& A)
{
	const int nr = A.NrRows();
	const int nc = A.NrCols();
	T **ptr = new T *[nr];
	T *block = new T[nr*nc];
	for (int i=0;i<nr;i++)
		*(ptr+i) = block+i*nc;
	for (int i=0;i<nr;i++)
		for (int j=0;j<nc;j++)
			ptr[i][j] = A[i][j];
	return ptr;
}

template<class T>
void freematrix(T **ptr)
{
	delete[] *ptr;
	delete[] ptr;
}

}

#endif


