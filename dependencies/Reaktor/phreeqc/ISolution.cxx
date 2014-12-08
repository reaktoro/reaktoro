// ISolution.cxx: implementation of the cxxSolutionxx class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <cassert>				// assert
#include <algorithm>			// std::sort
#include <sstream>

#include "Utils.h" // define first
#include "Phreeqc.h"
#include "ISolution.h"
#include "phqalloc.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxISolution::cxxISolution(PHRQ_io *io)
:
units("mMol/kgw")
{
	default_pe = "pe";
	cxxChemRxn temp_pe_reactions;
	pe_reactions[default_pe] = temp_pe_reactions;

}
cxxISolution::~cxxISolution()
{
}

#ifdef SKIP_OR_MOVE_TO_STRUCTURES
void
cxxISolution::ConvertUnits(Phreeqc * phreeqc_ptr)
  //
  // Converts from input units to moles per kilogram water
  //
{
	LDBLE sum_solutes = 0;
	// foreach conc
	std::map < std::string, cxxISolutionComp >::iterator iter =
		this->comps.begin();
	for (; iter != this->comps.end(); ++iter)
	{
		struct master *master_ptr = phreeqc_ptr-> master_bsearch(iter->first.c_str());
		if (master_ptr != NULL && (master_ptr->minor_isotope == TRUE))
			continue;
		//if (iter->second.Get_description() == "H(1)" || iter->second.Get_description() == "E") continue;
		if (strcmp(iter->second.Get_description().c_str(), "H(1)") == 0
			|| strcmp(iter->second.Get_description().c_str(), "E"))
			continue;
		if (iter->second.get_input_conc() <= 0.0)
			continue;
/*
*   Convert liters to kg solution
*/
		LDBLE moles = iter->second.get_input_conc();
		if (this->units.find("/l") != std::string::npos)
		{
			moles /= this->density;
		}
/*
* Convert to moles
*/
		//set gfw for element
		iter->second.set_gfw(phreeqc_ptr);
		// convert to moles
		if (iter->second.get_units().find("g/") != std::string::npos)
		{
			if (iter->second.get_gfw() != 0)
			{
				moles /= iter->second.get_gfw();
			}
			else
			{
				std::ostringstream oss;
				oss << "Could not find gfw, " << iter->second.
					Get_description();
				error_msg(oss.str().c_str(), CONTINUE);
			}
		}
/*
*   Convert milli or micro
*/
		char c = iter->second.get_units().c_str()[0];
		if (c == 'm')
		{
			moles *= 1e-3;
		}
		else if (c == 'u')
		{
			moles *= 1e-6;
		}
		iter->second.set_moles(moles);
/*
*   Sum grams of solute, convert from moles necessary
*/
		sum_solutes += moles * (iter->second.get_gfw());
	}
/*
 *   Convert /kgs to /kgw
 */
	LDBLE l_mass_water;
	if ((this->units.find("kgs") != std::string::npos) ||
		(this->units.find("/l") != std::string::npos))
	{
		l_mass_water = 1.0 - 1e-3 * sum_solutes;
		for (; iter != this->comps.end(); ++iter)
		{
			iter->second.set_moles(iter->second.get_moles() / l_mass_water);
		}
	}
/*
 *   Scale by mass of water in solution
 */
	l_mass_water = this->mass_water;
	for (; iter != this->comps.end(); ++iter)
	{
		iter->second.set_moles(iter->second.get_moles() * l_mass_water);
	}
}
#endif

#ifdef SKIP
void
cxxISolution::dump_xml(std::ostream & os, unsigned int indent) const const
{
	unsigned int i;

	for (i = 0; i < indent; ++i)
		os << Utilities::INDENT;
	os << "<solution>\n";

	cxxNumKeyword::dump_xml(os, indent);

	for (i = 0; i < indent + 1; ++i)
		os << Utilities::INDENT;
	os << "<temp>" << this->get_tc() << "</temp>" << "\n";

	for (i = 0; i < indent + 1; ++i)
		os << Utilities::INDENT;
	os << "<pH>" << this->get_ph() << "</pH>" << "\n";

	for (i = 0; i < indent + 1; ++i)
		os << Utilities::INDENT;
	os << "<pe>" << this->get_solution_pe() << "</pe>" << "\n";

	assert(this->pe.size() > 0);
	assert(this->default_pe >= 0);
	assert(this->pe.size() > (unsigned int) this->default_pe);

	for (i = 0; i < indent + 1; ++i)
		os << Utilities::INDENT;
	os << "<units>" << this->get_units() << "</units>" << "\n";

	for (i = 0; i < indent + 1; ++i)
		os << Utilities::INDENT;
	os << "<density>" << this->get_density() << "</density>" << "\n";

	// foreach conc
	if (!this->totals.empty())
	{
		for (i = 0; i < indent + 1; ++i)
			os << Utilities::INDENT;
		os << "<totals>\n";

		std::vector < cxxISolutionComp >::const_iterator iter =
			this->totals.begin();
		for (; iter != this->totals.end(); ++iter)
		{
			(*iter).dump_xml(*this, os, indent + 2);
		}

		for (i = 0; i < indent + 1; ++i)
			os << Utilities::INDENT;
		os << "</totals>\n";
	}

	// foreach isotope
	if (!this->isotopes.empty())
	{
		for (i = 0; i < indent + 1; ++i)
			os << Utilities::INDENT;
		os << "<isotopes>\n";

		std::list < cxxIsotope >::const_iterator iter =
			this->isotopes.begin();
		for (; iter != this->isotopes.end(); ++iter)
		{
			(*iter).dump_xml(os, indent + 2);
		}

		for (i = 0; i < indent + 1; ++i)
			os << Utilities::INDENT;
		os << "</isotopes>\n";
	}

	for (i = 0; i < indent + 1; ++i)
		os << Utilities::INDENT;
	os << "<water>" << this->get_mass_water() << "</water>" << "\n";

	for (i = 0; i < indent; ++i)
		os << Utilities::INDENT;
	os << "</solution>" << "\n";
}
#endif
