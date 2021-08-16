#ifndef INHERITANCEVECTOR_HEADER
#define INHERITANCEVECTOR_HEADER

#include <iostream>
#include <vector>

namespace ibd
{

//! Vector with binary inheritance indicators
/*!
   Maximum length of the inheritance vector is 15 (or 31 for VC8), because inheritance vector
   is implemented as an unsigned int.
*/

class InhVector
{
public:
	InhVector(int len, unsigned int init = 0);
	InhVector operator++(int) {InhVector r=*this;this->x++;return r;}
	bool next_indic();
	operator unsigned int() const { return x; }
	bool end() const { return x >= max; }
	void print(std::ostream& outp) const;
	unsigned int length() const { return n; }
private:
	unsigned int x;
	unsigned int max;
	unsigned int n;
};

std::ostream& operator<<(std::ostream& outp, const InhVector& h);

}

#endif


