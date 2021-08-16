#include "mainR.h"
#include "read_map.h"
#include "crosses.h"
#include "popt.h"
#include "util_genetics.h"

using namespace std;
using namespace ibd;
using namespace Rcpp;

vector<IndProp> make_ped_file(const string& poptype,
                              const vector<string>& ID)
{
  vector<IndProp> result;
  int npar;
  if (poptype.find("C4") == 0)
    npar = 4;
  else if (poptype.find("C3") == 0)
    npar = 3;
  else
    npar = 2;
  // header, parents
  string par1, par2;
  for (int i=0; i<npar;i++)
  {
    IndProp tmp = IndProp(ID[i],"*","INBPAR","0","0");
    result.push_back(tmp);
  }
  if (npar == 4) {
    IndProp tmp1 = IndProp("H1","*","HYBRID",ID[0],ID[1]);
    IndProp tmp2 = IndProp("H2","*","HYBRID",ID[2],ID[3]);
    result.push_back(tmp1);
    result.push_back(tmp2);
    par1 = "H1";
    par2 = "H2";
  } else if (npar == 3)
  {
    IndProp tmp = IndProp("H","*","HYBRID",ID[0],ID[1]);
    result.push_back(tmp);
    par1 = "H";
    par2 = ID[2];
  } else
  {
    par1 = ID[0];
    par2 = ID[1];
  }
  for (unsigned int i=npar;i<ID.size();i++)
  {
    IndProp tmp = IndProp(ID[i],"ID_FAM",poptype,par1,par2);
    result.push_back(tmp);
  }
  return result;
}

void marker_selection(LinkageMap& markermap,
                      LinkageMap& eval_pos,
                      int sel_chr,
                      bool analysis_fam,
                      bool pos_option)
{
  markermap = adjust_markermap(markermap);
  if (sel_chr != -1)
  {
    eval_pos = select_chr(eval_pos,sel_chr);
    markermap = select_chr(markermap,sel_chr);
  }
  const int M = eval_pos.size();
  if (M < 1)
    throw ibd_error("no evaluation points!");
  if (analysis_fam)
  {
    if (markermap.empty())
    {
      markermap.push_back(eval_pos.front());
      markermap.push_back(eval_pos.back());
    }
    else if (eval_pos.front() < markermap.front())
    {
      markermap.push_back(eval_pos.front());
      sort(markermap.begin(),markermap.end());
    }
    if (eval_pos.back() > markermap.back())
      markermap.push_back(eval_pos.back());
  }
  else
  {
    if (pos_option)
    {
      const int M = eval_pos.size();
      for (int m=0;m<M;m++)
        markermap.push_back(eval_pos[m]);
      sort(markermap.begin(),markermap.end());
    }
  }
  if (markermap.size() == 1)
  {
    Locus loc = markermap[0];
    markermap.push_back(Locus(loc.GetChr()+1,loc.GetPosition(),EXTR_POS));
  }
}

int linking_data(matrix<score>& geno,
                 const LinkageMap& markermap,
                 const vector<IndProp>& pop,
                 const matrix<score>& geno_locfile,
                 const vector<string>& ID,
                 const vector<string>& marker_names)
{
  // put all data in right format
  const int M = markermap.size();
  const int N = pop.size();
  map<string,int> ID_ndx = make_index(ID);

  map<string,int> marker_names_ndx = make_index(marker_names);

  geno = matrix<score>(N,M,Uscore);
  for (int i=0;i<N;i++)
  {
    map<string,int>::const_iterator it=ID_ndx.find(pop[i].GetID());
    if (it != ID_ndx.end())
    {
      if (pop[i].IsHybrid())
      {
        Rcout << "!Warning: Genotypic data for " << pop[i].GetID()
              << " will be ignored" << endl;
      }
      else
      {
        int ndx_ind = it->second;
        for (int m=0;m<M;m++)
        {
          string locname = markermap[m].GetName();
          map<string,int>::const_iterator it = marker_names_ndx.find(locname);
          if (it != marker_names_ndx.end())
          {
            int ndx_M = it->second;
            geno[i][m] = geno_locfile[ndx_ind][ndx_M];
          }
        }
      }
    }
  }
  return 0;
}

int count_scores(const vector<score>& geno)
{
  int M = geno.size();
  int cnt = 0;
  for (int m=0;m<M;m++)
  {
    if (geno[m] != Uscore)
      cnt++;
  }
  return cnt;
}

int main_pedigreeR(arma::cube& Z,
                   vector<string>& parents,
                   vector<string>& offspring,
                   LinkageMap& positions,
                   const string& poptype,
                   const string& locfile,
                   const string& mapfile,
                   const DataFrame& eval_pos_df,
                   const double& max_step_size,
                   const bool& grid = true,
                   const bool& verbose = false)
{
  // read all the data
  matrix<score> geno;
  LinkageMap markermap = read_map_file(mapfile);
  bool analysis_family = true;
  int sel_chr = -1;

  matrix<score> geno_locfile;
  vector<string> ID, marker_names;
  if (verbose)
  {
    Rcout << "reading data .............." << endl;
  }
  read_flapjackfile(ID, marker_names, geno_locfile, locfile);

  vector<IndProp> pop = make_ped_file(poptype, ID);

  markermap = reduce_markermap(markermap, marker_names);
  print_marker_warnings(markermap, marker_names);
  LinkageMap eval_pos;
  bool pos_option = eval_pos_df.length() > 0;
  if (pos_option)
  {
    eval_pos = read_eval_pos_df(eval_pos_df);
  }
  else if (max_step_size > 0) 
  {
    if (grid)
    {
      eval_pos = generate_grid_map(markermap, max_step_size);
    }
    else {
      eval_pos = generate_extended_map(markermap, max_step_size);
    }
  }
  else
    eval_pos = markermap;
  // Count number of inbred founders:
  int Nfnd = 0;
  const int Npop = pop.size();
  for (int i = 0; i < Npop; i++)
  {
    if (pop[i].IsFounder())
      Nfnd++;
  }
  marker_selection(markermap, eval_pos, sel_chr, analysis_family, pos_option);
  linking_data(geno, markermap, pop, geno_locfile, ID, marker_names);
  // start of analysis.
  Z = analysis_cross(pop, geno, markermap, eval_pos, verbose);
  const int npar = count_parents(pop);
  for (int i=0; i < npar; i++)
    parents.push_back(ID[i]);
  for (unsigned int i = npar; i < ID.size(); i++)
    offspring.push_back(ID[i]);

  positions = eval_pos;

  return 0;
}



