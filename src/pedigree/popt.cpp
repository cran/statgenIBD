#include "ibdexcept.h"
#include "popt.h"
#include "util_genetics.h"

using namespace ibd;
using namespace std;

pop_base * init_pop(const std::string& poptype)
{
	int x,y;
	// biparental crosses:
	if (poptype == "DH")				return new popDH();
	if (match(x,poptype,"Fx"))			return new popFx(x);
	if (match(x,poptype,"FxDH"))		return new popFxDH(x);
	if (match(x,poptype,"BCx"))			return new popBCx(x);
	if (match(x,poptype,"BCxDH"))		return new popBCxDH(x);

	if (match(x,y,poptype,"BCxSy"))		return new popBCxSy(x,y);
	if (match(x,y,poptype,"BCxSyDH"))	return new popBCxSyDH(x,y);

	// three-way crosses:
	if (poptype == "C3")				return new popC3Sx(0);
	if (poptype == "C3DH")				return new popC3SxDH(0);
	if (match(x,poptype,"C3Sx"))		return new popC3Sx(x);
	if (match(x,poptype,"C3SxDH"))		return new popC3SxDH(x);
	if (match(x,y,poptype,"C3RCxSy"))	return new popC3RCxSy(x,y);
	if (match(x,y,poptype,"C3RCxSyDH"))	return new popC3RCxSyDH(x,y);

	// four-way crosses:
	if (poptype == "C4")				return new popC4Sx(0);
	if (poptype == "C4DH")				return new popC4SxDH(0);
	if (match(x,poptype,"C4Sx"))		return new popC4Sx(x);
	if (match(x,poptype,"C4SxDH"))		return new popC4SxDH(x);

	throw ibd_error("unknown type " + poptype);
	return 0;
}

OrdGeno popDH::gen_off(const std::vector<OrdGeno>& p,
                       InhVector x) const
{
	OrdGeno F1 = cross(p[0],0,p[1],0);
	OrdGeno g = DH(F1,x);
	return g;
}

OrdGeno popFx::gen_off(const std::vector<OrdGeno>& p,
                       InhVector x) const
{
	OrdGeno F1 = cross(p[0],0,p[1],0);
	OrdGeno g = selfing(F1,x,ngen_self);
	return g;
}

OrdGeno popFxDH::gen_off(const std::vector<OrdGeno>& p,
                         InhVector x) const
{
	OrdGeno F1 = cross(p[0],0,p[1],0);
	OrdGeno g = selfing(F1,x,ngen_self);
	g = DH(g,x);
	return g;
}

// backcross with second parent
OrdGeno popBCx::gen_off(const std::vector<OrdGeno>& p,
                        InhVector x) const
{
	OrdGeno g = BC(p[0],p[1],x,ngen);
	return g;
}

// backcross with second parent
OrdGeno popBCxDH::gen_off(const std::vector<OrdGeno>& p,
                          InhVector x) const
{
	OrdGeno g = BC(p[0],p[1],x,ngen);
	return DH(g,x);
}

// BCxSy population
OrdGeno popBCxSy::gen_off(const std::vector<OrdGeno>& p,
                          InhVector x) const
{
	OrdGeno g1 = BC(p[0],p[1],x,ngen_BC);
	OrdGeno g2 = selfing(g1,x,ngen_self);
	return g2;
}

// BCxSyDH population
OrdGeno popBCxSyDH::gen_off(const std::vector<OrdGeno>& p,
                            InhVector x) const
{
	OrdGeno g1 = BC(p[0],p[1],x,ngen_BC);
	OrdGeno g2 = selfing(g1,x,ngen_self);
	return DH(g2,x);
}

// three way cross: (AxB)xC, followed by x generations of selfing.
OrdGeno popC3Sx::gen_off(const std::vector<OrdGeno>& p,
                         InhVector x) const
{
	const OrdGeno& A = p[0];
	const OrdGeno& B = p[1];
	const OrdGeno& C = p[2];
	OrdGeno AB = cross(A,0,B,0);
	bool h = x.next_indic();
	OrdGeno ABxC = cross(AB,h,C,0);
	OrdGeno g = selfing(ABxC,x,ngen_self);
	return g;
}

// three way cross: (AxB)xC, followed by x generations of selfing and DH
OrdGeno popC3SxDH::gen_off(const std::vector<OrdGeno>& p,
                           InhVector x) const
{
	const OrdGeno& A = p[0];
	const OrdGeno& B = p[1];
	const OrdGeno& C = p[2];
	OrdGeno AB = cross(A,0,B,0);
	bool h = x.next_indic();
	OrdGeno ABxC = cross(AB,h,C,0);
	OrdGeno g = selfing(ABxC,x,ngen_self);
	return DH(g,x);
}

// 3 parents: (AxB), followed by x gen of RC with C, y gen selfing
OrdGeno popC3RCxSy::gen_off(const std::vector<OrdGeno>& p,
                            InhVector x) const
{
	const OrdGeno& A = p[0];
	const OrdGeno& B = p[1];
	const OrdGeno& C = p[2];
	OrdGeno AB = cross(A,0,B,0);
	OrdGeno g = RC(AB,C,x,ngen_RC);
	g = selfing(g,x,ngen_self);
	return g;
}

// 3 parents: (AxB), followed by x gen of RC with C, y gen selfing,DH
OrdGeno popC3RCxSyDH::gen_off(const std::vector<OrdGeno>& p,
                              InhVector x) const
{
	const OrdGeno& A = p[0];
	const OrdGeno& B = p[1];
	const OrdGeno& C = p[2];
	OrdGeno AB = cross(A,0,B,0);
	OrdGeno g = RC(AB,C,x,ngen_RC);
	g = selfing(g,x,ngen_self);
	return DH(g,x);
}

// four way cross: (AxB)x(CxD), followed by x generations of selfing
OrdGeno popC4Sx::gen_off(const std::vector<OrdGeno>& p,
                         InhVector x) const
{
	const OrdGeno& A = p[0];
	const OrdGeno& B = p[1];
	const OrdGeno& C = p[2];
	const OrdGeno& D = p[3];
	OrdGeno AB = cross(A,0,B,0);
	OrdGeno CD = cross(C,0,D,0);
	bool h1 = x.next_indic();
	bool h2 = x.next_indic();
	OrdGeno ABxCD = cross(AB,h1,CD,h2);
	OrdGeno g = selfing(ABxCD,x,ngen_self);
	return g;
}

// four way cross: (AxB)x(CxD), followed by x generations of selfing, and DH
OrdGeno popC4SxDH::gen_off(const std::vector<OrdGeno>& p,
                           InhVector x) const
{
	const OrdGeno& A = p[0];
	const OrdGeno& B = p[1];
	const OrdGeno& C = p[2];
	const OrdGeno& D = p[3];
	OrdGeno AB = cross(A,0,B,0);
	OrdGeno CD = cross(C,0,D,0);
	bool h1 = x.next_indic();
	bool h2 = x.next_indic();
	OrdGeno ABxCD = cross(AB,h1,CD,h2);
	OrdGeno g = selfing(ABxCD,x,ngen_self);
	return DH(g,x);
}
