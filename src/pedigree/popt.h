#ifndef POPULATION_TYPE_HEADER
#define POPULATION_TYPE_HEADER

#include <string>
#include <vector>

#include "OrdGeno.h"
#include "InhVector.h"

//! abstract virtual base class for population types
class pop_base
{
public:
	pop_base(int len_inh_vector) : n(len_inh_vector) {}
	virtual ~pop_base(){}

	unsigned int get_len() const { return n; }
	virtual ibd::OrdGeno gen_off(const std::vector<ibd::OrdGeno>& p,
                              ibd::InhVector x) const = 0;
private:
	unsigned int n;
};

//! population type DH (doubled haploids)
class popDH : public pop_base
{
public:
	popDH() : pop_base(1) {}
	virtual ~popDH(){}

  ibd::OrdGeno gen_off(const std::vector<ibd::OrdGeno>& p, ibd::InhVector x) const;
};

//! population type Fx (where x is number of generations)
class popFx : public pop_base
{
public:
	popFx(int x) : pop_base(2*x-2), ngen_self(x-1) {}
	virtual ~popFx(){}

  ibd::OrdGeno gen_off(const std::vector<ibd::OrdGeno>& p, ibd::InhVector x) const;
private:
	int ngen_self;
};

//! population type FxDH (x is a variable)
class popFxDH : public pop_base
{
public:
	popFxDH(int x) : pop_base(2*x-1), ngen_self(x-1) {}
	virtual ~popFxDH(){}

  ibd::OrdGeno gen_off(const std::vector<ibd::OrdGeno>& p, ibd::InhVector x) const;
private:
	int ngen_self;
};

//! population type BCx (x is a variable)
class popBCx : public pop_base
{
public:
	popBCx(int x) : pop_base(x), ngen(x) {}
	virtual ~popBCx(){}

  ibd::OrdGeno gen_off(const std::vector<ibd::OrdGeno>& p, ibd::InhVector x) const;
private:
	int ngen;
};

//! population type BCxDH (x is a variable)
class popBCxDH : public pop_base
{
public:
	popBCxDH(int x) : pop_base(x+1), ngen(x) {}
	virtual ~popBCxDH(){}

  ibd::OrdGeno gen_off(const std::vector<ibd::OrdGeno>& p, ibd::InhVector x) const;
private:
	int ngen;
};

//! population type BCxSy (x gen. backcross, followed by y gen. of selfing)
class popBCxSy : public pop_base
{
public:
	popBCxSy(int x, int y) : pop_base(x+2*y), ngen_BC(x), ngen_self(y) {}
	virtual ~popBCxSy(){}

  ibd::OrdGeno gen_off(const std::vector<ibd::OrdGeno>& p, ibd::InhVector x) const;
private:
	int ngen_BC;
	int ngen_self;
};

//! population type BCxSyDH (as BCxSy, and DH in last gen.)
class popBCxSyDH : public pop_base
{
public:
	popBCxSyDH(int x, int y) : pop_base(x+2*y+1), ngen_BC(x), ngen_self(y) {}
	virtual ~popBCxSyDH(){}

  ibd::OrdGeno gen_off(const std::vector<ibd::OrdGeno>& p, ibd::InhVector x) const;
private:
	int ngen_BC;
	int ngen_self;
};

//! three way cross, followed by selfing
class popC3Sx : public pop_base
{
public:
	popC3Sx(int x) : pop_base(2*x+1), ngen_self(x) {}
	virtual ~popC3Sx(){}

  ibd::OrdGeno gen_off(const std::vector<ibd::OrdGeno>& p, ibd::InhVector x) const;
private:
	int ngen_self;
};

//! three way cross, followed by selfing and DH generation
class popC3SxDH : public pop_base
{
public:
	popC3SxDH(int x) : pop_base(2*x+2), ngen_self(x) {}
	virtual ~popC3SxDH(){}

  ibd::OrdGeno gen_off(const std::vector<ibd::OrdGeno>& p, ibd::InhVector x) const;
private:
	int ngen_self;
};

//! AxB, followed by x gen. of RC with C, and y gen. of selfing
class popC3RCxSy : public pop_base
{
public:
	popC3RCxSy(int x, int y) : pop_base(x+2*y), ngen_RC(x), ngen_self(y) {}
	virtual ~popC3RCxSy(){}

  ibd::OrdGeno gen_off(const std::vector<ibd::OrdGeno>& p, ibd::InhVector x) const;
private:
	int ngen_RC;
	int ngen_self;
};

//! AxB, followed by x gen. of RC with C, and y gen. of selfing, and DH gen.
class popC3RCxSyDH : public pop_base
{
public:
	popC3RCxSyDH(int x, int y) : pop_base(x+2*y+1), ngen_RC(x), ngen_self(y) {}
	virtual ~popC3RCxSyDH(){}

  ibd::OrdGeno gen_off(const std::vector<ibd::OrdGeno>& p, ibd::InhVector x) const;
private:
	int ngen_RC;
	int ngen_self;
};

//! four way cross, followed by selfing
class popC4Sx : public pop_base
{
public:
	popC4Sx(int x) : pop_base(2*x+2), ngen_self(x) {}
	virtual ~popC4Sx(){}

  ibd::OrdGeno gen_off(const std::vector<ibd::OrdGeno>& p, ibd::InhVector x) const;
private:
	int ngen_self;
};

//! four way cross, followed by selfing and DH generation
class popC4SxDH : public pop_base
{
public:
	popC4SxDH(int x) : pop_base(2*x+3), ngen_self(x) {}
	virtual ~popC4SxDH(){}

  ibd::OrdGeno gen_off(const std::vector<ibd::OrdGeno>& p, ibd::InhVector x) const;
private:
	int ngen_self;
};

pop_base * init_pop(const std::string& poptype);

#endif

