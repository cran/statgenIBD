#include "convert.h"
#include "Loc.h"

using namespace std;
using namespace ibd;

MapFunction mapfunction = invhaldane;

int compare(const Locus& locA,
            const Locus& locB)
{
	if (locA.chr < locB.chr) return -1;
	if (locA.chr > locB.chr) return  1;

	if (locA.pos < locB.pos) return -1;
	if (locA.pos > locB.pos) return  1;

	return 0;
}

double recomb(const Locus& locA,
              const Locus& locB)
{
	if (locA.chr != locB.chr)
		return 0.5;
	double distance = fabs(locA.pos-locB.pos);
	return mapfunction(distance);
}

vector<double> make_rec_map(const LinkageMap& linkage_map)
{
	const int nr_intervals = linkage_map.size() - 1; // number of intervals
	vector<double> r(nr_intervals);
	for (int i=0;i<nr_intervals;i++)
		r[i] = recomb(linkage_map[i],linkage_map[i+1]);
	return r;
}

int pos_qtl(const LinkageMap& Markermap,
            const Locus& QTLpos)
{
	int nloc = Markermap.size();
	for (int i=0;i<nloc-1;i++)
		if (QTLpos >= Markermap[i] && QTLpos <= Markermap[i+1])
			return i;
	throw ibd_error("Evaluation point not in interval!");
	return 0; // dummy
}

vector<ibd::Interval> make_intervals(const LinkageMap& markermap)
{
	vector<Interval> result;
	int nloc = markermap.size();
	int chr_nr = -1;
	double left = 0.0;
	double right;
	for (int i=0;i<nloc;i++)
	{
		if (chr_nr != markermap[i].GetChr())
		{
			chr_nr = markermap[i].GetChr();
			left = markermap[i].GetPosition();
		}
		if ((i==nloc-1) || markermap[i].GetChr() != markermap[i+1].GetChr())
		{
			right = markermap[i].GetPosition();
			if (right-left < 10.0) // short interval
				result.push_back(Interval(-5.0+0.5*(right+left),5.0+0.5*(right+left)));
			else
				result.push_back(Interval(left,right));
		}
	}
	return result;
}

double total_length(const LinkageMap& markermap)
{
	vector<Interval> intervals = make_intervals(markermap);
	double sum = 0.0;
	int Nchr = intervals.size();
	for (int i=0;i<Nchr;i++)
		sum += intervals[i].Length();
	return sum;
}

LinkageMap generate_extended_map(const LinkageMap& Markermap,
                                 double max_step_size)
{
	LinkageMap extended_map;
	int nloc = Markermap.size();
	int m;
	for (m=0;m<nloc-1;m++)
	{
		Locus Left = Markermap[m];
		Locus Right = Markermap[m+1];
		extended_map.push_back(Left);
		if (Left.GetChr() == Right.GetChr())
		{
			double left_pos = Left.GetPosition();
			double right_pos = Right.GetPosition();
			double dist = right_pos - left_pos;
			int N = (int) ceil(dist/max_step_size - 1.0e-5);
			for (int i=1;i<N;i++)
			{
				int chr = Left.GetChr();
				double pos = round(left_pos + (i/(1.0*N))*dist,2);
				string name = "EXT_" + stringify(chr) + "_" + stringify(pos);
				extended_map.push_back(Locus(chr,pos,name));
			}
		}
	}
	extended_map.push_back(Markermap[m]);

	return extended_map;
}

LinkageMap generate_grid_map(const LinkageMap& Markermap,
                             double max_step_size)
{
  LinkageMap grid_map;
  int nloc = Markermap.size();
  int chr_nr = -1;
  double left = 0.0;
  double right;
  for (int i=0;i<nloc;i++)
  {
    if (chr_nr != Markermap[i].GetChr())
    {
      chr_nr = Markermap[i].GetChr();
      left = Markermap[i].GetPosition();
    }
    if ((i==nloc-1) || Markermap[i].GetChr() != Markermap[i+1].GetChr())
    {
      right = Markermap[i].GetPosition();
      double dist = right - left;
      int N = (int) ceil(dist/max_step_size - 1.0e-5);
      for (int i=1;i<N;i++)
      {
        int chr = chr_nr;
        double pos = round(left + (i/(1.0*N))*dist,2);
        string name = "EXT_" + stringify(chr) + "_" + stringify(pos);
        grid_map.push_back(Locus(chr,pos,name));
      }
    }
  }
  return grid_map;
}

bool eval_pos(const Locus& loc)
{
	return (loc.GetName().find(EVAL_POS) == 0);
}
