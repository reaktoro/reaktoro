#if !defined(PPASSEMBLAGE_H_INCLUDED)
#define PPASSEMBLAGE_H_INCLUDED

#include <cassert>				// assert
#include <map>					// std::map
#include <string>				// std::string
#include <list>					// std::list
#include <vector>				// std::vector

#include "NumKeyword.h"
#include "PPassemblageComp.h"
class cxxMix;

class cxxPPassemblage:public cxxNumKeyword
{

  public:
	cxxPPassemblage(PHRQ_io * io=NULL);
	cxxPPassemblage(const std::map < int, cxxPPassemblage > &entity_map,
					  cxxMix & mx, int n_user, PHRQ_io * io=NULL);
	~cxxPPassemblage();

	void dump_raw(std::ostream & s_oss, unsigned int indent, int *n_out=NULL) const;

	void read_raw(CParser & parser, bool check = true);

	const cxxNameDouble & Get_assemblage_totals() const
	{
		return this->assemblage_totals;
	};
	const cxxNameDouble & Get_eltList() const
	{
		return this->eltList;
	};
	void Set_eltList(cxxNameDouble & nd) {this->eltList = nd;}
	std::map <std::string, cxxPPassemblageComp > & Get_pp_assemblage_comps() 
	{
		return this->pp_assemblage_comps;
	};
	const std::map <std::string, cxxPPassemblageComp > & Get_pp_assemblage_comps() const
	{
		return this->pp_assemblage_comps;
	};
	void  Set_pp_assemblage_comps(std::map <std::string, cxxPPassemblageComp > & c) 
	{
		this->pp_assemblage_comps = c;
	};
	bool Get_new_def(void) const {return this->new_def;}
	void Set_new_def(bool tf) {this->new_def = tf;}

	cxxPPassemblageComp *Find(const std::string name);

	void totalize(Phreeqc * phreeqc_ptr);

protected:
	void add(const cxxPPassemblage & addee, LDBLE extensive);
	// not written
	void dump_xml(std::ostream & os, unsigned int indent = 0) const;

protected:
	bool new_def;
	std::map <std::string, cxxPPassemblageComp > pp_assemblage_comps;
	cxxNameDouble eltList;              // list of elements in phases (and alternate reactions)
	cxxNameDouble assemblage_totals;    // after totalize, total moles of elements in the PPassemblage
	const static std::vector < std::string > vopts;
};

#endif // !defined(PPASSEMBLAGE_H_INCLUDED)
