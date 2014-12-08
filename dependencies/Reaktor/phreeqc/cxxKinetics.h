#if !defined(CXXKINETICS_H_INCLUDED)
#define CXXKINETICS_H_INCLUDED

#include <cassert>				// assert
#include <map>					// std::map
#include <string>				// std::string
#include <list>					// std::list
#include <vector>				// std::vector
#include "phrqtype.h"
#include "NumKeyword.h"
#include "KineticsComp.h"
#include "PHRQ_base.h"
class cxxMix;

class cxxKinetics:public cxxNumKeyword
{

  public:
	cxxKinetics(PHRQ_io *io=NULL);
	cxxKinetics(const std::map < int, cxxKinetics > &entity_map, cxxMix & mx,
				int n_user, PHRQ_io *io=NULL);
	~cxxKinetics();

	//void dump_xml(std::ostream& os, unsigned int indent = 0)const;

	void dump_raw(std::ostream & s_oss, unsigned int indent, int * n_out=NULL) const;

	void read_raw(CParser & parser, bool check = true);

	std::vector < LDBLE > &Get_steps(void) {return steps;}
	const std::vector < LDBLE > &Get_steps(void)const {return steps;}
	LDBLE Get_step_divide(void) const {return step_divide;}
	void Set_step_divide(LDBLE t) {step_divide = t;}
	int Get_rk(void) const {return rk;};
	void Set_rk(int t) {rk = t;}
	int Get_bad_step_max(void) const {return bad_step_max;}
	void Set_bad_step_max(int t) {bad_step_max = t;}
	bool Get_use_cvode(void) const {return use_cvode;}
	void Set_use_cvode(bool tf) {use_cvode = tf;}
	int Get_cvode_steps(void) const {return cvode_steps;}
	void Set_cvode_steps(int t) {cvode_steps = t;}
	int Get_cvode_order(void) const {return cvode_order;}
	void Set_cvode_order(int t) {cvode_order = t;}
	std::vector < cxxKineticsComp > &Get_kinetics_comps(void) {return kinetics_comps;}
	const std::vector < cxxKineticsComp > &Get_kinetics_comps(void)const {return kinetics_comps;}
	cxxNameDouble & Get_totals(void) {return this->totals;}
	void Set_totals(const cxxNameDouble &nd) {totals = nd;}
	bool Get_equalIncrements(void) const {return equalIncrements;}
	void Set_equalIncrements(bool tf) {equalIncrements = tf;}
	int Get_count(void) const {return count;}
	void Set_count(int i) {count = i;}
	int Get_reaction_steps(void) const;
	cxxKineticsComp * Find(const std::string &str);
	LDBLE Current_step(const bool incremental_reactions, const int reaction_step) const;

  protected:
	void add(const cxxKinetics & addee, LDBLE extensive);

  protected:
	// KINETICS_MODIFY candidates
	std::vector < cxxKineticsComp > kinetics_comps;
	std::vector < LDBLE >steps;
	int count;
	bool equalIncrements;
	LDBLE step_divide;
	int rk;
	int bad_step_max;
	bool use_cvode;
	int cvode_steps;
	int cvode_order;
	// internal variables
	cxxNameDouble totals;
	const static std::vector < std::string > vopts;
};

#endif // !defined(CXXKINETICS_H_INCLUDED)
