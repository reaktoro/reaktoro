#if !defined(GASCOMP_H_INCLUDED)
#define GASCOMP_H_INCLUDED

#include <cassert>				// assert
#include <map>					// std::map
#include <string>				// std::string
#include <list>					// std::list
#include <vector>				// std::vector

#include "NameDouble.h"

class cxxGasComp: public PHRQ_base
{

  public:
	cxxGasComp(PHRQ_io *io=NULL);

	virtual ~cxxGasComp();


	void dump_raw(std::ostream & s_oss, unsigned int indent) const;

	bool read_raw(CParser & parser, bool check=true);

	std::string Get_phase_name(void) const {return this->phase_name;}
	void Set_phase_name(std::string s) {this->phase_name = s;}
	LDBLE Get_p_read() const {return this->p_read;}
	void Set_p_read(LDBLE t) {this->p_read = t;}
	LDBLE Get_moles() const {return this->moles;}
	void Set_moles(LDBLE t) {this->moles = t;}
	LDBLE Get_initial_moles() const {return this->initial_moles;}
	void Set_initial_moles(LDBLE t) {this->initial_moles = t;}

	void add(const cxxGasComp & addee, LDBLE extensive);
	void multiply(LDBLE extensive);
	void Serialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles);
	void Deserialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles, int &ii, int &dd);
	
  protected:
	std::string phase_name;
	// GAS_PHASE_MODIFY candidates
	LDBLE moles;
	// GAS_PHASE_MODIFY candidates with new_def=true
	LDBLE p_read;
	// internal workspace
	LDBLE initial_moles;
	const static std::vector < std::string > vopts;
};

#endif // !defined(GASCOMP_H_INCLUDED)
