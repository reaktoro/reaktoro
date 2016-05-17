#if !defined(EXCHCOMP_H_INCLUDED)
#define EXCHCOMP_H_INCLUDED

#include <cassert>				// assert
#include <map>					// std::map
#include <string>				// std::string
#include <list>					// std::list
#include <vector>				// std::vector

#include "NameDouble.h"

class cxxExchComp: public PHRQ_base
{

  public:
	cxxExchComp(PHRQ_io *io=NULL);
	virtual ~cxxExchComp();

	void dump_xml(std::ostream & os, unsigned int indent = 0) const;

	void dump_raw(std::ostream & s_oss, unsigned int indent) const;

	void read_raw(CParser & parser, bool check=true);

	const std::string &Get_formula() const
	{
		return this->formula;
	}
	void Set_formula(const char *cstring)
	{
		if (cstring != NULL)
			this->formula = std::string(cstring);
		else
			this->formula.clear();
	}
	LDBLE Get_la() const
	{
		return this->la;
	}
	void Set_la(LDBLE d)
	{
		this->la = d;
	}
	LDBLE Get_charge_balance() const
	{
		return this->charge_balance;
	}
	void Set_charge_balance(LDBLE d)
	{
		this->charge_balance = d;
	}
	const std::string &Get_phase_name() const
	{
		return this->phase_name;
	}
	void Set_phase_name(const char *cstring)
	{
		if (cstring != NULL)
			this->phase_name = std::string(cstring);
		else
			this->phase_name.clear();
	}
	LDBLE Get_phase_proportion() const
	{
		return this->phase_proportion;
	}
	void Set_phase_proportion(LDBLE d)
	{
		this->phase_proportion = d;
	}
	const std::string &Get_rate_name() const
	{
		return this->rate_name;
	}
	void Set_rate_name(const char *cstring)
	{
		if (cstring != NULL)
			this->rate_name = std::string(cstring);
		else
			this->rate_name.clear();
	}
	LDBLE Get_formula_z() const
	{
		return this->formula_z;
	}
	void Set_formula_z(LDBLE d)
	{
		this->formula_z = d;
	}
	void Set_totals(struct elt_list *e_l, int count)
	{
		this->totals = cxxNameDouble(e_l, count);
	}
	void Set_totals(struct elt_list *e_l)
	{
		this->totals = cxxNameDouble(e_l);
	}
	void Set_totals(cxxNameDouble nd)
	{
		this->totals = nd;
	}
	cxxNameDouble & Get_totals() {return (this->totals);}
	const cxxNameDouble & Get_totals()const {return (this->totals);}

	void add(const cxxExchComp & comp, LDBLE extensive);
	void multiply(LDBLE extensive);
	void Serialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles);
	void Deserialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles, int &ii, int &dd);

protected:
	std::string formula;
	// EXCHANGE_MODIFY candidates
	cxxNameDouble totals;
	LDBLE la;
	LDBLE charge_balance;
	std::string phase_name;
	LDBLE phase_proportion;
	std::string rate_name;
	LDBLE formula_z;			// charge on formula
	const static std::vector < std::string > vopts;
  public:

};

#endif // !defined(EXCHCOMP_H_INCLUDED)
