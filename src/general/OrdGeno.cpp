#include "OrdGeno.h"
#include "convert.h"

using namespace std;
using namespace ibd;

string OrdGeno::print_string() const
{
	string str = "[" + stringify(this->first) + "," + stringify(this->second) + "]";
	return str;
}

int ibd::get_haplo(const OrdGeno& g, bool h)
{
	return (h == 0) ? g.first : g.second;
}

OrdGeno ibd::cross(const OrdGeno& g1, bool h1, const OrdGeno& g2, bool h2)
{
	int first  = get_haplo(g1,h1);
	int second = get_haplo(g2,h2);
	return OrdGeno(first,second);
}

OrdGeno ibd::selfing(const OrdGeno& g, InhVector& u)
{
	bool h1 = u.next_indic();
	bool h2 = u.next_indic();
	return cross(g,h1,g,h2);
}

OrdGeno ibd::selfing(OrdGeno g, InhVector& u, int Ngen)
{
	for (int i=0;i<Ngen;i++)
		g = selfing(g,u);
	return g;
}

OrdGeno ibd::DH(const OrdGeno& g, InhVector& u)
{
	bool h = u.next_indic();
	int a = get_haplo(g,h);
	return OrdGeno(a,a);
}

OrdGeno ibd::BC(const OrdGeno& donor, const OrdGeno& background, InhVector& u, int ngen)
{
	OrdGeno g = cross(donor,0,background,0);
	for (int i=0;i<ngen;i++)
	{
		bool h = u.next_indic();
		g = cross(g,h,background,0);
	}
	return g;
}

// RC: Recurrent cross with background (background assumed to be inbred)
OrdGeno ibd::RC(const OrdGeno& A, const OrdGeno& background, InhVector& u, int ngen)
{
	OrdGeno g = A;
	for (int i=0;i<ngen;i++)
	{
		bool h = u.next_indic();
		g = cross(g,h,background,0);
	}
	return g;
}

std::ostream& ibd::operator<<(std::ostream& outp, const OrdGeno& g)
{
	string str = g.print_string();
	int nr_whitespaces = outp.width() - str.length();
	outp.width(0);
	if (nr_whitespaces > 0)
		outp << string(nr_whitespaces,' ');
	outp << str;
	return outp;
}
