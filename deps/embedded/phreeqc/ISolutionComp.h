#if !defined(ISOLUTIONCOMP_H_INCLUDED)
#define ISOLUTIONCOMP_H_INCLUDED

#include <string.h>             // ::strcmp
#include <string>
#include <map>                  // std::map
#include <vector>
#include <set>
#include "phrqtype.h"           // LDBLE
#include "Parser.h"             // CParser
#include "PHRQ_base.h"          // PHRQ_base

// forward declarations
class cxxSolution;
class cxxISolution;             // reqd for read and dump_xml
class PHRQ_io;

class cxxISolutionComp: public PHRQ_base
{
  public:
	cxxISolutionComp(PHRQ_io *io=NULL);
	virtual ~cxxISolutionComp(void);

  public:

	CParser::STATUS_TYPE read(const char *line, cxxSolution *solution_ptr);

	void dump_xml(std::ostream & os, unsigned int indent = 0) const;

	const std::string &Get_description() const
	{
		return this->description;
	}
	void Set_description(const char *l_description)
	{
		if (l_description != NULL)
			this->description = std::string(l_description);
		else
			this->description.clear();
	}

	LDBLE Get_moles() const
	{
		return this->moles;
	}
	void Set_moles(LDBLE l_moles)
	{
		this->moles = l_moles;
	}

	LDBLE Get_input_conc() const
	{
		return this->input_conc;
	}
	void Set_input_conc(LDBLE l_input_conc)
	{
		this->input_conc = l_input_conc;
	}

	std::string Get_units()const
	{
		return this->units;
	}
	void Set_units(const char *l_units)
	{
		if (l_units != NULL)
			this->units = std::string(l_units);
		else
			this->units.clear();
	}

	const std::string &Get_equation_name() const
	{
		return this->equation_name;
	}
	void Set_equation_name(const char *l_equation_name)
	{
		if (l_equation_name != NULL)
			this->equation_name = std::string(l_equation_name);
		else
			this->equation_name.clear();

	}

	LDBLE Get_phase_si() const
	{
		return this->phase_si;
	}
	void Set_phase_si(int l_phase_si)
	{
		this->phase_si = l_phase_si;
	}
	std::string Get_pe_reaction() const
	{
		return this->pe_reaction;
	}
	void Set_pe_reaction(const std::string & pe_r)
	{
		this->pe_reaction = pe_r;
	}
	const std::string &Get_as() const
	{
		return this->as;
	}
	void Set_as(const char *l_as)
	{
		if (l_as != NULL)
			this->as = std::string(l_as);
		else
			this->as.clear();
	}
	LDBLE Get_gfw() const
	{
		return this->gfw;
	};
	void Set_gfw(LDBLE l_gfw)
	{
		this->gfw = l_gfw;
	}

	bool operator<(const cxxISolutionComp & conc) const
	{
		return ::strcmp(this->description.c_str(), conc.description.c_str()) < 0;
	}

  protected:
	  std::string description;
	  LDBLE moles;
	  LDBLE input_conc;
	  std::string units;
	  std::string equation_name;
	  LDBLE phase_si;
	  std::string pe_reaction;
	  std::string as;
	  LDBLE gfw;
};

#endif // ISOLUTIONCOMP_H_INCLUDED
