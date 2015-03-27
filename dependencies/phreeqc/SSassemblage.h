#if !defined(SSASSEMBLAGE_H_INCLUDED)
#define SSASSEMBLAGE_H_INCLUDED


#include <cassert>				// assert
#include <map>					// std::map
#include <string>				// std::string
#include <list>					// std::list
#include <vector>				// std::vector

#include "NumKeyword.h"
#include "NameDouble.h"
#include "SS.h"

class cxxSS;

//#include "cxxMix.h"
class cxxMix;

class cxxSSassemblage:public cxxNumKeyword
{

public:
	cxxSSassemblage(PHRQ_io * io = NULL);
	cxxSSassemblage(const std::map < int, cxxSSassemblage > &entity_map,
		cxxMix & mx, int n_user, PHRQ_io * io = NULL);
	~cxxSSassemblage();

	//void dump_xml(std::ostream& os, unsigned int indent = 0)const;

	void dump_raw(std::ostream & s_oss, unsigned int indent, int *n_out=NULL) const;

	void read_raw(CParser & parser, bool check = true);

	void totalize(Phreeqc * phreeqc_ptr);

	const cxxNameDouble & Get_totals() const {return this->totals;}
	std::map < std::string, cxxSS > & Get_SSs(void) {return SSs;}
	const std::map < std::string, cxxSS > & Get_SSs(void)const {return SSs;}
	void Set_SSs(std::map < std::string, cxxSS > & ss) {SSs = ss;}
	bool Get_new_def(void) const {return new_def;}
	void Set_new_def(bool tf) {new_def = tf;}
	std::vector<cxxSS *> Vectorize(void);
	void add(const cxxSSassemblage & addee, LDBLE extensive);
	cxxSS *Find(const std::string &s);

protected:
	// SOLID_SOLUTION_MODIFY candidate
	std::map < std::string, cxxSS > SSs;
	// SOLID_SOLUTION keyword data
	bool new_def;
	// internal variables
	cxxNameDouble totals;
	const static std::vector < std::string > vopts;
};

#endif // !defined(SSASSEMBLAGE_H_INCLUDED)
