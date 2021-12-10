#include <fstream>
#include <Rcpp.h>

#include "read_map.h"
#include "convert.h"

using namespace ibd;
using namespace std;

LinkageMap select_chr(const LinkageMap& markermap,
                      std::string sel_chr)
{
  const int M = markermap.size();
  LinkageMap sel_markermap;
  for (int m=0;m<M;m++)
  {
    Locus loc = markermap[m];
    if (loc.GetChr() == sel_chr)
      sel_markermap.push_back(loc);
  }
  return sel_markermap;
}

LinkageMap adjust_markermap(const LinkageMap& markermap)
{
  const int M = markermap.size();
  const double eps = 1.0e-6;
  const double delta = 0.001;
  LinkageMap new_markermap;
  for (int m=0;m<M;m++)
  {
    Locus loc = markermap[m];
    if (!new_markermap.empty())
    {
      if (loc.GetChr() == new_markermap.back().GetChr())
      {
        double dist = loc.GetPosition() - new_markermap.back().GetPosition();
        if (dist < eps)
        {
          string n = loc.GetName();
          double p = new_markermap.back().GetPosition();
          std::string c = loc.GetChr();
          loc = Locus(c,p+delta,n);
        }
      }
    }
    new_markermap.push_back(loc);
  }
  return new_markermap;
}

LinkageMap reduce_markermap(const LinkageMap& markermap,
                            const vector<string>& markers)
{
  const int M = markermap.size();
  LinkageMap new_markermap;
  for (int m=0;m<M;m++)
  {
    Locus loc = markermap[m];
    if (((find(markers.begin(),markers.end(),loc.GetName()) != markers.end())
           || (loc.GetName() == EVAL_POS) || loc.GetName() == EXTR_POS))
    {
      new_markermap.push_back(loc);
    }
  }
  return new_markermap;
}

void print_marker_warnings(const LinkageMap& markermap,
                           const vector<string>& markers)
{
  const int M = markermap.size();
  vector<string> ignored_markers;
  for (vector<string>::const_iterator it = markers.begin();it!=markers.end();it++)
  {
    string name = *it;
    bool marker_found = false;
    for (int m=0;m<M;m++)
    {
      if (markermap[m].GetName() == name)
      {
        marker_found = true;
        break;
      }
    }
    if (!marker_found)
      ignored_markers.push_back(name);
  }
  if (!ignored_markers.empty())
  {
    Rcpp::Rcout << endl << "WARNING: " << endl
               << "The following markers will be ignored (not defined in the map file): " << endl;
    for (vector<string>::const_iterator it=ignored_markers.begin();it!=ignored_markers.end();it++)
      Rcpp::Rcout << *it << endl;
    Rcpp::Rcout << endl << endl;
  }
}

LinkageMap read_map_file(const string& filename)
{
  LinkageMap markermap;
  ifstream inp(filename.c_str());
  if (!inp)
    throw ibd_error("Cannot read file " + filename);
  string str;
  while (inp)
  {
    double pos;
    string name;
    std::string chr;
    inp >> name >> chr >> pos;
    Locus loc(chr,pos,name);
    markermap.push_back(loc);
    eatcomment(inp);
  }
  return markermap;
}

LinkageMap read_eval_pos_df(const Rcpp::DataFrame& evalposdf)
{
  Rcpp::NumericVector pos = evalposdf["pos"];
  Rcpp::CharacterVector chr = evalposdf["chr"];
  LinkageMap positions;
  for (int m=0;m<evalposdf.nrows();m++)
  {
    std::string chrM = Rcpp::as<std::string>(chr[m]);    
    std::string name = "EVAL_" + chrM + "_" + stringify(pos[m]);
    double posM = pos[m];
	Locus loc(chrM,posM,name);
	positions.push_back(loc);    
  }

  sort(positions.begin(),positions.end());
  
  return(positions);
}


