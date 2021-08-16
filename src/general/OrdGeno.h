#ifndef ORDERED_GENOTYPE_HEADER
#define ORDERED_GENOTYPE_HEADER

#include <string>
#include <iostream>

#include "InhVector.h"

namespace ibd
{

//! Ordered Genotype
/*!
   The ordered genotype is represented by a pair of integers, where
     this->first is the allele from the first parent and
	 this->second is the allele from the second parent.
*/

class OrdGeno : public std::pair<int,int>
{
public:
	OrdGeno() {}
	OrdGeno(int a, int b) : std::pair<int,int>(a,b) {}
	bool homozygous() const { return this->first == this->second; }
	std::string print_string() const;
};

std::ostream& operator<<(std::ostream& outp, const OrdGeno& g);

int get_haplo(const OrdGeno& g, bool h);
OrdGeno cross(const OrdGeno& g1, bool h1, const OrdGeno& g2, bool h2);

OrdGeno selfing(const OrdGeno& g, InhVector& u);
OrdGeno selfing(OrdGeno g, InhVector& u, int Ngen);
OrdGeno DH(const OrdGeno& g, InhVector& u);
OrdGeno BC(const OrdGeno& donor, const OrdGeno& background, InhVector& u, int ngen);

// RC: Recurrent cross with background (background assumed to be inbred)
OrdGeno RC(const OrdGeno& A, const OrdGeno& background, InhVector& u, int ngen);

}

#endif

