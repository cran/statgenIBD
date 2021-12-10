#ifndef LOCUS_HEADER
#define LOCUS_HEADER

#include <string>
#include <vector>

#include "misc.h"

class Locus
{
public:
	Locus() {}
	Locus(std::string c, double p, const std::string& name = "")
		: chr(c), pos(p), loc_name(name) {}
	std::string GetChr() const { return chr; }
	double GetPosition() const { return pos; }
	std::string GetName() const { return loc_name; }
private:
  std::string chr;
	double pos;
	std::string loc_name;
	friend int compare(const Locus& , const Locus&);
	friend double recomb(const Locus& , const Locus&);
};

const std::string EVAL_POS = "__EVALPOS";
const std::string EXTR_POS = "__EXTRPOS";

int compare(const Locus&, const Locus&);
double recomb(const Locus& , const Locus&);

inline bool operator==(const Locus& locA, const Locus& locB)
{ return (compare(locA,locB) == 0); }

inline bool operator!=(const Locus& locA, const Locus& locB)
{ return (compare(locA,locB) != 0); }

inline bool operator<=(const Locus& locA, const Locus& locB)
{ return (compare(locA,locB) <= 0); }

inline bool operator>=(const Locus& locA, const Locus& locB)
{ return (compare(locA,locB) >= 0); }

inline bool operator<(const Locus& locA, const Locus& locB)
{ return (compare(locA,locB) < 0); }

inline bool operator>(const Locus& locA, const Locus& locB)
{ return (compare(locA,locB) > 0); }

typedef double (*MapFunction)(double cM);
inline double invhaldane(double cM) { return 0.5*(1.0-exp(-0.02*cM));}
inline double kosambi(double cM) { return 0.5*(exp(0.04*cM)-1.0)/(exp(0.04*cM)+1.0);}

extern MapFunction mapfunction;

typedef std::vector<Locus> LinkageMap;

std::vector<double> make_rec_map(const LinkageMap&);

int pos_qtl(const LinkageMap& Markermap, const Locus& QTLpos);

std::vector<ibd::Interval> make_intervals(const LinkageMap& markermap);

double total_length(const LinkageMap& markermap);

LinkageMap generate_extended_map(const LinkageMap& Markermap, double max_step_size);

LinkageMap generate_grid_map(const LinkageMap& Markermap, double max_step_size);

bool eval_pos(const Locus& loc);

#endif


