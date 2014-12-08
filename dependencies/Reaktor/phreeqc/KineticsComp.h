#if !defined(KINETICSCOMP_H_INCLUDED)
#define KINETICSCOMP_H_INCLUDED

#include <cassert>				// assert
#include <map>					// std::map
#include <string>				// std::string
#include <list>					// std::list
#include <vector>				// std::vector

#include "NameDouble.h"

class cxxKineticsComp: public PHRQ_base
{

public:
	cxxKineticsComp(PHRQ_io *io=NULL);
	cxxKineticsComp(struct kinetics_comp *, PHRQ_io *io=NULL);
	virtual ~cxxKineticsComp();

	void dump_xml(std::ostream & os, unsigned int indent = 0) const;

	void dump_raw(std::ostream & s_oss, unsigned int indent) const;

	void read_raw(CParser & parser, bool check = true);

	const std::string &Get_rate_name() const
	{
		return this->rate_name;
	}
	void Set_rate_name(const char * s)
	{
		if (s != NULL)
			this->rate_name = std::string(s);
		else
			this->rate_name.clear();
	}

	cxxNameDouble &Get_namecoef(void) {return namecoef;}
	const cxxNameDouble &Get_namecoef(void)const {return namecoef;}
	void Set_namecoef(const cxxNameDouble nd) {namecoef = nd;}
	LDBLE Get_tol(void) const {return tol;}	
	void Set_tol(LDBLE t) {tol = t;}	
	LDBLE Get_m(void) const {return m;}	
	void Set_m(LDBLE t) {m = t;}	
	LDBLE Get_m0(void) const {return m0;}
	void Set_m0(LDBLE t) {m0 = t;}	
	LDBLE Get_moles(void) const {return moles;}
	void Set_moles(LDBLE t) {moles = t;}
	LDBLE Get_initial_moles(void) const {return initial_moles;}
	void Set_initial_moles(LDBLE t) {initial_moles = t;}	
	std::vector < LDBLE > &Get_d_params(void) {return d_params;};
	const std::vector < LDBLE > &Get_d_params(void)const {return d_params;};
	std::vector < std::string > &Get_c_params(void) {return c_params;};
	cxxNameDouble &Get_moles_of_reaction(void) {return moles_of_reaction;}
	const cxxNameDouble &Get_moles_of_reaction(void)const {return moles_of_reaction;}
	void Set_moles_of_reaction(const cxxNameDouble nd) {moles_of_reaction = nd;}

	void add(const cxxKineticsComp & comp, LDBLE extensive);
	void multiply(LDBLE extensive);

  protected:
	  std::string rate_name;
	  // KINETICS_MODIFY candidates
	  cxxNameDouble namecoef;		//stoichiometry of reaction
	  LDBLE tol;
	  LDBLE m;
	  LDBLE m0;
	  std::vector< LDBLE > d_params;
	  std::vector< std::string > c_params;  // Not used
	  // kinetics workspace variables
	  LDBLE moles;
	  LDBLE initial_moles;
	  cxxNameDouble moles_of_reaction;
	  const static std::vector < std::string > vopts;
  public:

};

#endif // !defined(KINETICSCOMP_H_INCLUDED)
