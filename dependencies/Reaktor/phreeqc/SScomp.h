#if !defined(SSCOMP_H_INCLUDED)
#define SSCOMP_H_INCLUDED

#include <cassert>				// assert
#include <map>					// std::map
#include <string>				// std::string
#include <list>					// std::list
#include <vector>				// std::vector

#include "NameDouble.h"

class cxxSScomp: public PHRQ_base
{

  public:
	cxxSScomp(PHRQ_io *io=NULL);
	virtual ~cxxSScomp();

	void dump_raw(std::ostream & s_oss, unsigned int indent) const;
	void read_raw(CParser & parser, bool check = true);

	const std::string &Get_name() const	{return this->name;}
	void Set_name(std::string s) {this->name = s;}

	LDBLE Get_initial_moles() const {return this->initial_moles;}
	void Set_initial_moles(LDBLE t) {this->initial_moles = t;}
	LDBLE Get_moles() const {return this->moles;}
	void Set_moles(LDBLE t) {this->moles = t;}
	LDBLE Get_init_moles() const {return this->init_moles;}
	void Set_init_moles(LDBLE t) {this->init_moles = t;}
	LDBLE Get_delta() const {return this->delta;}
	void Set_delta(LDBLE t) {this->delta = t;}
	LDBLE Get_fraction_x() const {return this->fraction_x;}
	void Set_fraction_x(LDBLE t) {this->fraction_x = t;}
	LDBLE Get_log10_lambda() const {return this->log10_lambda;}
	void Set_log10_lambda(LDBLE t) {this->log10_lambda = t;}
	LDBLE Get_log10_fraction_x() const {return this->log10_fraction_x;}
	void Set_log10_fraction_x(LDBLE t) {this->log10_fraction_x = t;}
	LDBLE Get_dn() const {return this->dn;}
	void Set_dn(LDBLE t) {this->dn = t;}
	LDBLE Get_dnc() const {return this->dnc;}
	void Set_dnc(LDBLE t) {this->dnc = t;}
	LDBLE Get_dnb() const {return this->dnb;}
	void Set_dnb(LDBLE t) {this->dnb = t;}

	void multiply(LDBLE extensive);

protected:
	std::string name;
	// SOLID_SOLUTION_MODIFY candidate identifier
	LDBLE moles;

	// Solid solution workspace variables
	LDBLE initial_moles;
	LDBLE init_moles;
	LDBLE delta;
	LDBLE fraction_x;
	LDBLE log10_lambda;
	LDBLE log10_fraction_x;
	LDBLE dn, dnc, dnb;
	const static std::vector < std::string > vopts;
};

#endif // !defined(SSCOMP_H_INCLUDED)
