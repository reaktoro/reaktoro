#if !defined(SURFACECOMP_H_INCLUDED)
#define SURFACECOMP_H_INCLUDED

#include <cassert>				// assert
#include <map>					// std::map
#include <string>				// std::string
#include <list>					// std::list
#include <vector>				// std::vector
#include "NameDouble.h"

class cxxSurfaceComp: public PHRQ_base
{

public:

	cxxSurfaceComp(PHRQ_io *io=NULL);
	virtual ~cxxSurfaceComp();

	void dump_xml(std::ostream & os, unsigned int indent = 0) const;
	void dump_raw(std::ostream & s_oss, unsigned int indent) const;
	void read_raw(CParser & parser, bool check = true);
	void add(const cxxSurfaceComp & comp, LDBLE extensive);
	void multiply(LDBLE extensive);

	const std::string &Get_formula() const {return this->formula;}
	void Set_formula(const char * f) {this->formula = f ? f : "";}
	LDBLE Get_formula_z(void) const {return formula_z;};
	void Set_formula_z(LDBLE t) {this->formula_z = t;}
	LDBLE Get_moles(void) const {return moles;}
	void Set_moles(LDBLE t) {this->moles = t;}
	cxxNameDouble & Get_totals() {return (this->totals);}
	void Set_totals(cxxNameDouble & nd) {this->totals = nd;}
	LDBLE Get_la(void) const {return la;};
	void Set_la(LDBLE t) {this->la = t;}
	LDBLE Get_charge_balance() const {return this->charge_balance;}
	void Set_charge_balance(LDBLE d) {this->charge_balance = d;}
	const std::string &Get_charge_name() const {return this->charge_name;}
	void Set_charge_name(const char * f) {this->charge_name = f ? f : "";}
	const std::string &Get_phase_name() const {return this->phase_name;}
	void Set_phase_name(const char * f) {this->phase_name = f ? f : "";}
	LDBLE Get_phase_proportion(void) const {return phase_proportion;}
	void Set_phase_proportion(LDBLE d) {this->phase_proportion = d;}
	const std::string &Get_rate_name() const {return this->rate_name;}
	void Set_rate_name(const char * f) {this->rate_name = f ? f : "";}
	LDBLE Get_Dw(void) const {return Dw;}
	void Set_Dw(LDBLE d) {this->Dw = d;}
	const std::string &Get_master_element() const {return this->master_element;}
	void Set_master_element(const char * f) {this->master_element = f ? f : "";}

protected:
	std::string formula;
	LDBLE formula_z;
	LDBLE moles;
	cxxNameDouble totals;
	LDBLE la;
	std::string charge_name;
	LDBLE charge_balance;
	std::string phase_name;
	LDBLE phase_proportion;
	std::string rate_name;
	LDBLE Dw;
	std::string master_element;
	const static std::vector < std::string > vopts;
public:

};

#endif // !defined(SURFACECOMP_H_INCLUDED)
