// SScomp.cxx: implementation of the cxxSScomp class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <cassert>				// assert
#include <algorithm>			// std::sort

#include "Utils.h"				// define first
#include "Phreeqc.h"
#include "SScomp.h"
#include "phqalloc.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxSScomp::cxxSScomp(PHRQ_io *io)
:
PHRQ_base(io)
	//
	// default constructor for cxxSScomp 
	//
{
	name = "";
	initial_moles = 0;
	moles = 0;
	init_moles = 0;
	delta = 0;
	fraction_x = 0;
	log10_lambda = 0;
	log10_fraction_x = 0;
	dn = 0;
	dnc = 0; 
	dnb = 0;
}
#ifdef SKIP
cxxSScomp::cxxSScomp(struct pure_phase * pure_phase_ptr, PHRQ_io *io)
:
PHRQ_base(io)
	//
	// constructor for cxxSScomp from struct pure_phase
	//
{
	this->Set_name(pure_phase_ptr->name);
	this->Set_add_formula(pure_phase_ptr->add_formula);
	si = pure_phase_ptr->si;
	si_org = pure_phase_ptr->si_org;
	moles = pure_phase_ptr->moles;
	delta = pure_phase_ptr->delta;
	initial_moles = pure_phase_ptr->initial_moles;
	force_equality = (pure_phase_ptr->force_equality == TRUE);
	dissolve_only = (pure_phase_ptr->dissolve_only == TRUE);
	precipitate_only = (pure_phase_ptr->precipitate_only == TRUE);
}
#endif
cxxSScomp::~cxxSScomp()
{
}

void
cxxSScomp::dump_raw(std::ostream & s_oss, unsigned int indent) const
{
	unsigned int i;
	s_oss.precision(DBL_DIG - 1);
	std::string indent0(""), indent1(""), indent2("");
	for (i = 0; i < indent; ++i)
		indent0.append(Utilities::INDENT);
	for (i = 0; i < indent + 1; ++i)
		indent1.append(Utilities::INDENT);
	for (i = 0; i < indent + 2; ++i)
		indent2.append(Utilities::INDENT);

	s_oss << indent1 << "# SOLID_SOLUTION_MODIFY candidate identifiers #\n";
	s_oss << indent1 << "-moles               " << this->moles << "\n";

	s_oss << indent1 << "# Solid solution workspace variables #\n";
	s_oss << indent1 << "-initial_moles       " << this->initial_moles << "\n";
	s_oss << indent1 << "-init_moles          " << this->init_moles << "\n";
	s_oss << indent1 << "-delta               " << this->delta << "\n";
	s_oss << indent1 << "-fraction_x          " << this->fraction_x << "\n";
	s_oss << indent1 << "-log10_lambda        " << this->log10_lambda << "\n";
	s_oss << indent1 << "-log10_fraction_x    " << this->log10_fraction_x << "\n";
	s_oss << indent1 << "-dn                  " << this->dn << "\n";
	s_oss << indent1 << "-dnc                 " << this->dnc << "\n";
	s_oss << indent1 << "-dnb                 " << this->dnb << "\n";
}

void
cxxSScomp::read_raw(CParser & parser, bool check)
{
	std::string str;

	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	int opt_save;

	opt_save = CParser::OPT_ERROR;
	bool initial_moles_defined(false);
	bool moles_defined(false);

	for (;;)
	{
		int opt = parser.get_option(vopts, next_char);
		if (opt == CParser::OPT_DEFAULT)
		{
			opt = opt_save;
		}

		switch (opt)
		{
		case CParser::OPT_EOF:
			break;
		case CParser::OPT_KEYWORD:
			break;
		case CParser::OPT_DEFAULT:
		case CParser::OPT_ERROR:
			opt = CParser::OPT_KEYWORD;
			// Allow return to Exchange for more processing
			break;

		case 0:				// name
			parser.error_msg("-Name ignored. Define with -component.");
			break;

		case 1:				// initial_moles
			if (!(parser.get_iss() >> this->initial_moles))
			{
				this->initial_moles = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for initial_moles.",
								 PHRQ_io::OT_CONTINUE);
			}
			initial_moles_defined = true;
			break;

		case 2:				// moles
			if (!(parser.get_iss() >> this->moles))
			{
				this->moles = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for moles.",
								 PHRQ_io::OT_CONTINUE);
			}
			moles_defined = true;
			break;

		case 3:				// init_moles
			if (!(parser.get_iss() >> this->init_moles))
			{
				this->init_moles = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for init_moles.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;

		case 4:				// delta
			if (!(parser.get_iss() >> this->delta))
			{
				this->delta = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for delta.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;

		case 5:				// fraction_x
			if (!(parser.get_iss() >> this->fraction_x))
			{
				this->fraction_x = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for fraction_x.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;


		case 6:				// log10_lambda
			if (!(parser.get_iss() >> this->log10_lambda))
			{
				this->log10_lambda = 0;
				parser.incr_input_error();
				parser.error_msg("Expected boolean value for log10_lambda.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;

		case 7:				// log10_fraction_x
			if (!(parser.get_iss() >> this->log10_fraction_x))
			{
				this->log10_fraction_x = 0;
				parser.incr_input_error();
				parser.error_msg("Expected boolean value for log10_fraction_x.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;

		case 8:				// dn
			if (!(parser.get_iss() >> this->dn))
			{
				this->dn = 0;
				parser.incr_input_error();
				parser.error_msg("Expected boolean value for dn.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;
		case 9:				// dnc
			if (!(parser.get_iss() >> this->dnc))
			{
				this->dnc = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for dnc.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;
		case 10:				// dnb
			if (!(parser.get_iss() >> this->dnb))
			{
				this->dnb = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for dnb.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;
		}
		if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD)
			break;
	}
	// members that must be defined
	if (check)
	{
		if (moles_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Moles not defined for PPassemblageComp input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (initial_moles_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Initial_moles not defined for PPassemblageComp input.",
				PHRQ_io::OT_CONTINUE);
		}
	}
}

#ifdef SKIP
void
cxxSScomp::totalize(Phreeqc * phreeqc_ptr)
{
	this->totals.clear();
	// component structures
	if (this->add_formula.size() != 0)
		return;
	struct phase *phase_ptr;
	int l;
	phase_ptr = phreeqc_ptr-> phase_bsearch(this->name.c_str(), &l, FALSE);
	if (phase_ptr != NULL)
	{
		cxxNameDouble phase_formula(phase_ptr->next_elt);
		this->totals.add_extensive(phase_formula, this->moles);
	}
	else
	{
		assert(false);
	}
	return;
}
#endif
void
cxxSScomp::multiply(LDBLE extensive)
{
	this->moles *= extensive;
	this->delta *= extensive;
	this->initial_moles *= extensive;
}
const std::vector< std::string >::value_type temp_vopts[] = {
	std::vector< std::string >::value_type("name"),	            // 0                 
	std::vector< std::string >::value_type("initial_moles"),	// 1
	std::vector< std::string >::value_type("moles"),	        // 2
	std::vector< std::string >::value_type("init_moles"),	    // 3
	std::vector< std::string >::value_type("delta"),	        // 4
	std::vector< std::string >::value_type("fraction_x"),	    // 5     
	std::vector< std::string >::value_type("log10_lambda"),	    // 6
	std::vector< std::string >::value_type("log10_fraction_x"),	// 7
	std::vector< std::string >::value_type("dn"),	            // 8
	std::vector< std::string >::value_type("dnc"),	            // 9
	std::vector< std::string >::value_type("dnb") 	            // 10
};									   
const std::vector< std::string > cxxSScomp::vopts(temp_vopts, temp_vopts + sizeof temp_vopts / sizeof temp_vopts[0]);	