#ifndef UTIL_GENETICS_HEADER
#define UTIL_GENETICS_HEADER

#include <string>
#include <vector>

//! Properties of an individual
/*!
This class contains the name (identity) of an individual, the names of the parents, and the type.

\warning At the moment 'Type' is defined as a string, it is better to define it as
enum for example.
*/
class IndProp
{
public:
	//! Default constructor
	IndProp() {}
	//! The constructor
	IndProp(std::string id, std::string fam, std::string type, std::string par1, std::string par2);

	//! Returns true if the individual is member of a family
	/*!
	  \warning At the moment we assume that DH indicates that the individual is member of a family.
	*/
	bool IsMemberFamily() const { return FAM != "*"; }

	//! Returns true if the individual is a (inbred) founder
	bool IsFounder() const { return TYPE == "INBFND"; }

	//! Returns true if the individual is a (inbred) parent
	bool IsInbredParent() const { return TYPE == "INBPAR"; }

	//! Returns true if the individual is a RIL (recombinant inbred line)
	bool IsRIL() const { return TYPE == "RIL"; }

	//! Returns true if the individual is a hybrid
	bool IsHybrid() const { return TYPE == "HYBRID"; }

	std::string GetID() const { return ID; }
	std::string GetFam() const { return FAM; }
	std::string GetP1() const { return P1; }
	std::string GetP2() const { return P2; }
	std::string GetType() const { return TYPE; }

private:
	std::string ID, FAM, TYPE, P1, P2;
};

int ndxID(const std::vector<IndProp>& pop, const std::string& ID);

bool match(int& x, const std::string& str, const char * pat);
bool match(int& x, int& y, const std::string& str, const char * pat);

#endif
