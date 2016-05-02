#if !defined(PPASSEMBLAGECOMP_H_INCLUDED)
#define PPASSEMBLAGECOMP_H_INCLUDED

#include <cassert>				// assert
#include <map>					// std::map
#include <string>				// std::string
#include <list>					// std::list
#include <vector>				// std::vector

#include "NameDouble.h"

class cxxPPassemblageComp: public PHRQ_base
{

  public:
	cxxPPassemblageComp(PHRQ_io *io=NULL);
	virtual ~cxxPPassemblageComp();

	void dump_xml(std::ostream & os, unsigned int indent = 0) const;
	void dump_raw(std::ostream & s_oss, unsigned int indent) const;
	void read_raw(CParser & parser, bool check = true);

	const std::string &Get_name() const	{return this->name;}
	void Set_name(const char * s)
	{
		if(s != NULL)
			this->name = std::string(s);
		else
			this->name.clear();
	}
	const std::string &Get_add_formula() const {return this->add_formula;}
	void Set_add_formula(const char * s)
	{
		if(s != NULL)
			this->add_formula = std::string(s);
		else
			this->add_formula.clear();
	}

	void totalize(Phreeqc * phreeqc_ptr);
	const cxxNameDouble & Get_totals() const {return (this->totals);}
	void Get_totals(cxxNameDouble & nd) {this->totals = nd;}
	LDBLE Get_si() const {return this->si;}
	void Set_si(LDBLE t) {this->si = t;}
	LDBLE Get_si_org() const {return this->si_org;}
	void Set_si_org(LDBLE t) {this->si_org = t;}
	LDBLE Get_moles() const {return this->moles;}
	void Set_moles(LDBLE t) {this->moles = t;}
	LDBLE Get_delta() const {return this->delta;}
	void Set_delta(LDBLE t) {this->delta = t;}
	LDBLE Get_initial_moles() const {return this->initial_moles;}
	void Set_initial_moles(LDBLE t) {this->initial_moles = t;}

	bool Get_force_equality() const {return this->force_equality;}
	void Set_force_equality(bool tf) {this->force_equality = tf;}
	bool Get_dissolve_only() const {return this->dissolve_only;}
	void Set_dissolve_only(bool tf) {this->dissolve_only = tf;}
	bool Get_precipitate_only() const {return this->precipitate_only;}
	void Set_precipitate_only(bool tf) {this->precipitate_only = tf;}

	void add(const cxxPPassemblageComp & comp, LDBLE extensive);
	void multiply(LDBLE extensive);
	void Serialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles);
	void Deserialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles, int &ii, int &dd);

protected:
	std::string name;
	std::string add_formula;
	LDBLE si;
	LDBLE si_org;
	LDBLE moles;
	LDBLE delta;
	LDBLE initial_moles;
	bool force_equality;
	bool dissolve_only;
	bool precipitate_only;
	cxxNameDouble totals;
	const static std::vector < std::string > vopts;
public:

};

#endif // !defined(PPASSEMBLAGECOMP_H_INCLUDED)
