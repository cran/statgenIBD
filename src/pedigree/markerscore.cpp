#include <iostream>
#include <fstream>
#include <algorithm>

#include "misc.h"
#include "ibdexcept.h"
#include "convert.h"
#include "markerscore.h"

using namespace std;
using namespace ibd;

score::score(int a, int b) : pair<int,int>(a,b)
{
	if (a < b)
		std::swap(this->first,this->second);
}

bool score::homozygous() const
{
	return !(this->first == U_haplo_sc || this->second == U_haplo_sc ||
			 this->first != this->second);
}

string score::print_string() const
{
	string f="*";
	string s="*";
	if (this->first != U_haplo_sc)
		f = stringify(this->first);
	if (this->second != U_haplo_sc)
		s = stringify(this->second);
	return "(" + f + "," + s + ")";
}

std::ostream& operator<<(std::ostream& outp,
                         const score& score)
{
	string str = score.print_string();
	int nr_whitespaces = outp.width() - str.length();
	outp.width(0);
	if (nr_whitespaces > 0)
		outp << string(nr_whitespaces,' ');
	outp << str;
	return outp;
}

// reads a character from input and checks whether this character is equal to c
void check_char(istream& s,
                char c)
{
	char x;
	s >> x;
	if (x != c)
		throw ibd_error("missing " + string(1,c));
}

int read_allele(istream& s)
{
	char c;
	s >> c;
	if (c == '-')
		return U_haplo_sc;
	s.putback(c);

	char x;
	if (s >> x)
	{
		return int(x);
	}

	throw ibd_error("error while reading allele");
	return 0; // dummy
}

score read_score(istream& line_stream, char delimit)
{
  string field;
  getline(line_stream,field,delimit);
  field.erase(remove_if(field.begin(), field.end(),::isspace), field.end());

  istringstream s(field);
  if (field == "-") return Uscore;

  string::size_type pos = field.find('/');
	if (pos == field.npos)     // score is of the form "a"
	{
     int a = read_allele(s);
     return score(a,a);
	}

  int a = read_allele(s);
	check_char(s,'/');
	int b = read_allele(s);

	return score(a,b);
}

bool check_score(const OrdGeno& g, const score& sc)
{
	if (sc == Uscore)
		return true;
	if (sc.second == U_haplo_sc)
		return (sc.first == g.first || sc.first == g.second);
	return (((g.first == sc.first)&&(g.second == sc.second))||
			((g.first == sc.second)&&(g.second == sc.first)));
}

string trim(const string& str)
{
    const char* whitespace = " \t\v\r\n";
    std::size_t start = str.find_first_not_of(whitespace);
    std::size_t end = str.find_last_not_of(whitespace);
    return start == end ? std::string() : str.substr(start, end - start + 1);
}


vector<string> read_tab_delimited_line(istream& inp)
{
  vector<string> result;
	string line,str;
	getline(inp,line);
  istringstream line_stream(line);
	while (getline(line_stream,str,'\t')) 
	{
		result.push_back(trim(str));
	}
	return result;
}

int read_flapjackfile(vector<string>& geno,
                      vector<string>& markers,
                      matrix<score>& scores,
                      const string filename)
{
  ifstream inp;
  OpenFile(inp,filename);

  string str;
  getline(inp,str,'\t');
	markers = read_tab_delimited_line(inp);

  const int M = markers.size();
	string line;
	while (getline(inp,line))
	{
		if (line.empty()) break;
    istringstream line_stream(line);
    string g;
    getline(line_stream,g,'\t');
	
	// Check if g not already in geno.
	if (std::find(geno.begin(), geno.end(), g) != geno.end())
		throw ibd_error("Genotypes in genofile should be unique");
	
    geno.push_back(g);

    vector<score> ind_scores;
    for (int m=0;m<M;m++)
    {
       score sc = read_score(line_stream,'\t');
       ind_scores.push_back(sc);
    }
    scores.push_back(ind_scores);
	}
  return 0;
}

// npar: number of parents
map<score,int> ndx_score(int npar)
{
	int k=0;
	map<score,int> ndx;
	for (int i=0;i<npar;i++)
		ndx[score(i,i)] = k++;
	if (npar == 2)
	{
		ndx[score(0,1)] = k++;
	}
	else if (npar == 3)
	{
		for (int i=0;i<2;i++)
			ndx[score(i,2)] = k++;
	}
	else if (npar == 4)
	{
		for (int i=0;i<2;i++)
			for (int j=2;j<4;j++)
				ndx[score(i,j)] = k++;
	}
	else
		throw ibd_error("npar > 4 !");
	return ndx;
}

