// KineticsComp.cxx: implementation of the cxxKineticsComp class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <cassert>				// assert
#include <algorithm>			// std::sort

#include "Utils.h"				// define first
#include "Phreeqc.h"
#include "KineticsComp.h"
//#include "Dictionary.h"
#include "phqalloc.h"
#include "Dictionary.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxKineticsComp::cxxKineticsComp(PHRQ_io *io)
:
PHRQ_base(io)
	//
	// default constructor for cxxKineticsComp 
	//
{
	tol = 1e-8;
	m = -1;
	m0 = -1;
	moles = 0.0;
	initial_moles = 0;
	namecoef.type = cxxNameDouble::ND_NAME_COEF;
}
cxxKineticsComp::~cxxKineticsComp()
{
}

#ifdef SKIP
void
cxxKineticsComp::dump_xml(std::ostream & s_oss, unsigned int indent) const const
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

	// Kinetics_Comp element and attributes

	s_oss << indent0 << "formula=\"" << this->formula << "\"" << "\n";
	s_oss << indent0 << "moles=\"" << this->moles << "\"" << "\n";
	s_oss << indent0 << "la=\"" << this->la << "\"" << "\n";
	s_oss << indent0 << "charge_balance=\"" << this->
		charge_balance << "\"" << "\n";
	if (this->phase_name != NULL)
	{
		s_oss << indent0 << "phase_name=\"" << this->
			phase_name << "\"" << "\n";
	}
	if (this->rate_name != NULL)
	{
		s_oss << indent0 << "rate_name=\"" << this->
			rate_name << "\"" << "\n";
	}
	s_oss << indent0 << "phase_proportion=\"" << this->
		phase_proportion << "\"" << "\n";

	// totals
	s_oss << indent0;
	s_oss << "<totals " << "\n";
	this->totals.dump_xml(s_oss, indent + 1);

	// formula_totals
	s_oss << indent0;
	s_oss << "<formula_totals " << "\n";
	this->formula_totals.dump_xml(s_oss, indent + 1);
}
#endif
void
cxxKineticsComp::dump_raw(std::ostream & s_oss, unsigned int indent) const
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

	// Kinetics_Comp element and attributes
	s_oss << indent1 << "# KINETICS_MODIFY candidate identifiers #\n";
	s_oss << indent1 << "-tol                   " << this->tol << "\n";
	s_oss << indent1 << "-m                     " << this->m << "\n";
	s_oss << indent1 << "-m0                    " << this->m0 << "\n";

	// namecoef
	s_oss << indent1;
	s_oss << "-namecoef" << "\n";
	this->namecoef.dump_raw(s_oss, indent + 2);

	// d_params
	s_oss << indent1;
	s_oss << "-d_params" << "\n";
	{
		int i = 0;
		s_oss << indent2;
		for (std::vector < LDBLE >::const_iterator it = d_params.begin();
			 it != d_params.end(); it++)
		{
			if (i++ == 5)
			{
				s_oss << "\n";
				s_oss << indent2;
				i = 0;
			}
			s_oss << *it << " ";
		}
		s_oss << "\n";
	}

	s_oss << indent1 << "# KineticsComp workspace variables #\n";
	s_oss << indent1 << "-moles                 " << this->moles << "\n";
	s_oss << indent1 << "-initial_moles         " << this->initial_moles << "\n";
}

void
cxxKineticsComp::read_raw(CParser & parser, bool check)
{
	std::string str;


	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	int opt_save;

	std::vector < LDBLE > temp_d_params;
	opt_save = CParser::OPT_ERROR;
	bool tol_defined(false);
	bool m_defined(false);
	bool m0_defined(false);
	bool d_params_defined(false);

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
			// Allow return to Kinetics for more processing
			break;

		case 0:				// rate_name not used
			parser.warning_msg("Rate_name ignored. Define in -comp.");
			break;

		case 1:				// tol
			if (!(parser.get_iss() >> this->tol))
			{
				this->tol = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for tol.",
								 PHRQ_io::OT_CONTINUE);
			}
			tol_defined = true;
			break;

		case 2:				// m
			if (!(parser.get_iss() >> this->m))
			{
				this->m = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for m.",
								 PHRQ_io::OT_CONTINUE);
			}
			m_defined = true;
			break;

		case 3:				// m0
			if (!(parser.get_iss() >> this->m0))
			{
				this->m0 = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for m0.",
								 PHRQ_io::OT_CONTINUE);
			}
			m0_defined = true;
			break;


		case 4:				// moles
			if (!(parser.get_iss() >> this->moles))
			{
				this->moles = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for moles.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;

		case 5:				// namecoef
			if (this->namecoef.read_raw(parser, next_char) !=
				CParser::PARSER_OK)
			{
				parser.incr_input_error();
				parser.
					error_msg
					("Expected element name and molality for namecoef.",
					 PHRQ_io::OT_CONTINUE);
			}
			opt_save = 5;
			break;

		case 6:				// d_params
			while (parser.copy_token(token, next_char) == CParser::TT_DIGIT)
			{
				double dd;
				sscanf(token.c_str(), "%lf", &dd);
				temp_d_params.push_back((LDBLE) dd);
				d_params_defined = true;
			}
			opt_save = 6;
			break;
		case 7:				// initial_moles
			if (!(parser.get_iss() >> this->initial_moles))
			{
				this->moles = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for initial_moles.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;
		}
		if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD)
			break;
	}

	if (d_params_defined)
	{
		this->d_params = temp_d_params;
	}
	if (check)
	{
		// members that must be defined
		if (tol_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Tol not defined for KineticsComp input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (m_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("M not defined for KineticsComp input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (m0_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("M0 not defined for KineticsComp input.",
				PHRQ_io::OT_CONTINUE);
		}
	}
}

void
cxxKineticsComp::add(const cxxKineticsComp & addee, LDBLE extensive)
{
	if (extensive == 0.0)
		return;
	if (addee.rate_name.size() == 0)
		return;
	// this and addee must have same name
	// otherwise generate a new KineticsComp with multiply
	if (this->rate_name.size() == 0 && addee.rate_name.size() == 0)
	{
		return;
	}
	assert(this->rate_name == addee.rate_name);
	this->m += addee.m * extensive;
	this->m0 += addee.m0 * extensive;
	this->moles += addee.moles * extensive;
}

void
cxxKineticsComp::multiply(LDBLE extensive)
{
	this->m *= extensive;
	this->m0 *= extensive;
	this->moles *= extensive;
}
void
cxxKineticsComp::Serialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles)
{
	ints.push_back(dictionary.Find(this->rate_name));
	this->namecoef.Serialize(dictionary, ints, doubles);
	doubles.push_back(this->tol);
	doubles.push_back(this->m);
	doubles.push_back(this->m0);
	ints.push_back((int) this->d_params.size());
	for (size_t i = 0; i < this->d_params.size(); i++)
	{
		doubles.push_back(d_params[i]);
	}
	doubles.push_back(this->moles);
	doubles.push_back(this->initial_moles);
	this->moles_of_reaction.Serialize(dictionary, ints, doubles);
}
void
cxxKineticsComp::Deserialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles, int &ii, int &dd)
{
	this->rate_name = dictionary.GetWords()[ints[ii++]];
	this->namecoef.Deserialize(dictionary, ints, doubles, ii, dd);
	this->tol = doubles[dd++];
	this->m = doubles[dd++];
	this->m0 = doubles[dd++];
	int n = ints[ii++];
	this->d_params.clear();
	for (int j = 0; j < n; j++)
	{
		this->d_params.push_back(doubles[dd++]);
	}
	this->moles = doubles[dd++];
	this->initial_moles = doubles[dd++];
	this->moles_of_reaction.Deserialize(dictionary, ints, doubles, ii, dd);
}
const std::vector< std::string >::value_type temp_vopts[] = {
	std::vector< std::string >::value_type("rate_name_not_used"),	// 0
	std::vector< std::string >::value_type("tol"),	                // 1
	std::vector< std::string >::value_type("m"),	                // 2
	std::vector< std::string >::value_type("m0"),	                // 3
	std::vector< std::string >::value_type("moles"),	            // 4
	std::vector< std::string >::value_type("namecoef"),	            // 5
	std::vector< std::string >::value_type("d_params"),	            // 6
	std::vector< std::string >::value_type("initial_moles") 	    // 7
};									   
const std::vector< std::string > cxxKineticsComp::vopts(temp_vopts, temp_vopts + sizeof temp_vopts / sizeof temp_vopts[0]);
