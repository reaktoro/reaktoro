// SurfaceComp.cxx: implementation of the cxxSurfaceComp class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <cassert>				// assert
#include <algorithm>			// std::sort

#include "Utils.h"				// define first
#include "Phreeqc.h"
#include "SurfaceComp.h"
#include "phqalloc.h"
#include "Dictionary.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxSurfaceComp::cxxSurfaceComp(PHRQ_io *io)
:
PHRQ_base(io)
//
// default constructor for cxxSurfaceComp 
//
{
	formula_z = 0.0;
	moles = 0.0;
	totals.type = cxxNameDouble::ND_ELT_MOLES;
	la = 0.0;
	charge_balance = 0;
	phase_proportion = 0.0;
	Dw = 0.0;
}
cxxSurfaceComp::~cxxSurfaceComp()
{
}

void
cxxSurfaceComp::dump_xml(std::ostream & s_oss, unsigned int indent) const
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

	// Surf_Comp element and attributes

	//s_oss << indent0 << "formula=\"" << this->formula << "\"" << "\n";
	s_oss << indent0 << "formula_z=\"" << this->formula_z << "\"" << "\n";
	s_oss << indent0 << "moles=\"" << this->moles << "\"" << "\n";
	s_oss << indent0 << "la=\"" << this->la << "\"" << "\n";
	s_oss << indent0 << "charge_balance=\"" << this->charge_balance << "\"" << "\n";
	s_oss << indent0 << "phase_proportion=\"" << this->phase_proportion << "\"" << "\n";
	s_oss << indent0 << "Dw=\"" << this->Dw << "\"" << "\n";
	s_oss << indent0 << "charge_name=\"" << this->charge_name << "\"" << "\n";
	if (this->phase_name.size() != 0)
	{
		s_oss << indent0 << "phase_name=\"" << this->phase_name << "\"" << "\n";
	}
	if (this->rate_name.size() != 0)
	{
		s_oss << indent0 << "rate_name=\"" << this->rate_name << "\"" << "\n";
	}

	// totals
	s_oss << indent0;
	s_oss << "<totals " << "\n";
	this->totals.dump_xml(s_oss, indent + 1);
}

void
cxxSurfaceComp::dump_raw(std::ostream & s_oss, unsigned int indent) const
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

	s_oss << indent0 << "# SURFACE_MODIFY candidate identifiers #\n";
	s_oss << indent0 << "-formula_z               " << this->formula_z << "\n";
	s_oss << indent0 << "-moles                   " << this->moles << "\n";
	s_oss << indent0 << "-la                      " << this->la << "\n";
	s_oss << indent0 << "-charge_balance          " << this->charge_balance << "\n";
	if (this->phase_name.size() != 0)
	{
		s_oss << indent0 << "-phase_name              " << this->phase_name << "\n";
	}
	if (this->rate_name.size() != 0)
	{
		s_oss << indent0 << "-rate_name               " << this->rate_name << "\n";
	}
	s_oss << indent0 << "-phase_proportion        " << this->phase_proportion << "\n";
	s_oss << indent0 << "-Dw                      " << this->Dw << "\n";
	s_oss << indent0 << "-charge_name             " << this->charge_name << "\n";
	s_oss << indent0 << "-master_element          " << this->master_element << "\n";
	// totals
	s_oss << indent0;
	s_oss << "-totals" << "\n";
	this->totals.dump_raw(s_oss, indent + 1);
}

void
cxxSurfaceComp::read_raw(CParser & parser, bool check)
{
	std::string str;


	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	int opt_save;

	opt_save = CParser::OPT_ERROR;

	bool master_element_defined(false);
	bool charge_name_defined(false);
	bool moles_defined(false);
	bool la_defined(false);
	bool charge_balance_defined(false);
	bool formula_z_defined(false);
	bool Dw_defined(false);
	bool totals_defined(false);

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
			// Allow return to Surface for more processing
			break;

		case 0:				// formula
			output_msg("-formula is obsolete in surface comp raw");
			break;

		case 1:				// moles
			if (!(parser.get_iss() >> this->moles))
			{
				this->moles = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for moles.",
					PHRQ_io::OT_CONTINUE);
			}
			moles_defined = true;
			break;

		case 2:				// la
			if (!(parser.get_iss() >> this->la))
			{
				this->la = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for la.",
					PHRQ_io::OT_CONTINUE);
			}
			la_defined = true;
			break;
		case 3:				// charge_number not used
			parser.warning_msg("-charge_number identifier is obsolete.");
			break;

		case 4:				// charge_balance
			if (!(parser.get_iss() >> this->charge_balance))
			{
				this->charge_balance = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for charge_balance.",
					PHRQ_io::OT_CONTINUE);
			}
			charge_balance_defined = true;
			break;

		case 5:				// phase_name
			if (!(parser.get_iss() >> str))
			{
				this->phase_name.clear();
				parser.incr_input_error();
				parser.error_msg("Expected string value for phase_name.",
					PHRQ_io::OT_CONTINUE);
			}
			else
			{
				this->phase_name = str;
			}
			break;

		case 6:				// rate_name
			if (!(parser.get_iss() >> str))
			{
				this->rate_name.clear();
				parser.incr_input_error();
				parser.error_msg("Expected string value for rate_name.",
					PHRQ_io::OT_CONTINUE);
			}
			else
			{
				this->rate_name = str;
			}
			break;

		case 7:				// phase_proportion
			if (!(parser.get_iss() >> this->phase_proportion))
			{
				this->phase_proportion = 0;
				parser.incr_input_error();
				parser.
					error_msg("Expected numeric value for phase_proportion.",
					PHRQ_io::OT_CONTINUE);
			}
			break;

		case 8:				// totals
			if (this->totals.read_raw(parser, next_char) !=
				CParser::PARSER_OK)
			{
				parser.incr_input_error();
				parser.
					error_msg
					("Expected element name and molality for SurfaceComp totals.",
					PHRQ_io::OT_CONTINUE);
			}
			totals_defined = true;
			opt_save = 8;
			break;

		case 9:				// formula_z
			if (!(parser.get_iss() >> this->formula_z))
			{
				this->formula_z = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for formula_z.",
					PHRQ_io::OT_CONTINUE);
			}
			formula_z_defined = true;
			break;

		case 10:				// formula_totals
			parser.warning_msg("-formula_totals is an obsolete identifier.");
			break;

		case 11:				// Dw
			if (!(parser.get_iss() >> this->Dw))
			{
				this->Dw = 0.0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for Dw.",
					PHRQ_io::OT_CONTINUE);
			}
			Dw_defined = true;
			break;
		case 12:				// charge_name
			if (!(parser.get_iss() >> str))
			{
				this->charge_name.clear();
				parser.incr_input_error();
				parser.error_msg("Expected string value for charge_name.",
					PHRQ_io::OT_CONTINUE);
			}
			else
			{
				this->charge_name = str;
			}
			charge_name_defined = true;
			break;
		case 13:				// master_element
			if (!(parser.get_iss() >> str))
			{
				this->master_element.clear();
				parser.incr_input_error();
				parser.error_msg("Expected string value for master_element.",
					PHRQ_io::OT_CONTINUE);
			}
			else
			{
				this->master_element = str;
			}
			master_element_defined = true;
			break;
		}
		if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD)
			break;
	}
	if (check)
	{
		// members that must be defined
		if (charge_name_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Charge_name not defined for SurfaceComp input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (formula_z_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Formula_z not defined for ExchComp input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (moles_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Moles not defined for SurfaceComp input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (la_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("La not defined for SurfaceComp input.",
				PHRQ_io::OT_CONTINUE);
		}

		if (charge_balance_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Charge_balance not defined for SurfaceComp input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (Dw_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Dw not defined for SurfaceComp input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (master_element_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Master_element name not defined for SurfaceComp input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (totals_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Totals not defined for SurfaceComp input.",
				PHRQ_io::OT_CONTINUE);
		}
	}
}

void
cxxSurfaceComp::add(const cxxSurfaceComp & addee, LDBLE extensive)
{
	if (extensive == 0.0)
		return;
	if (addee.formula.size() == 0)
		return;

	if (this->formula.size() == 0 && addee.formula.size() == 0)
	{
		return;
	}
	assert(this->formula == addee.formula);
	assert(this->formula_z == addee.formula_z);
	if (this->formula.size() == 0 && addee.formula.size() != 0)
	{
		this->formula = addee.formula;
	}

	// this and addee must have same formula
	// otherwise generate a new exchcomp with multiply
	LDBLE ext1, ext2, f1, f2;
	ext1 = this->moles;
	ext2 = addee.moles * extensive;
	if (ext1 + ext2 != 0)
	{
		f1 = ext1 / (ext1 + ext2);
		f2 = ext2 / (ext1 + ext2);
	}
	else
	{
		f1 = 0.5;
		f2 = 0.5;
	}
	this->moles += addee.moles * extensive;
	this->totals.add_extensive(addee.totals, extensive);
	this->la = f1 * this->la + f2 * addee.la;
	this->charge_balance += addee.charge_balance * extensive;
	if (this->phase_name != addee.phase_name)
	{
		std::ostringstream oss;
		oss <<
			"Cannot mix two Surface components with same formula and different related phases, "
			<< this->formula;
		error_msg(oss.str().c_str(), CONTINUE);
		return;
	}
	else if (this->phase_name.size() != 0)
	{
		this->phase_proportion =
			this->phase_proportion * f1 + addee.phase_proportion * f2;
	}

	if (this->rate_name != addee.rate_name)
	{
		std::ostringstream oss;
		oss <<
			"Cannot mix two exchange components with same formula and different related kinetics, "
			<< this->formula;
		error_msg(oss.str().c_str(), CONTINUE);
		return;
	}
	else if (this->rate_name.size() != 0)
	{
		//LDBLE phase_proportion;
		this->phase_proportion =
			this->phase_proportion * f1 + addee.phase_proportion * f2;
	}
	if ((this->rate_name.size() != 0 && addee.phase_name.size() != 0) ||
		(this->phase_name.size() != 0 && addee.rate_name.size() != 0))
	{
		std::ostringstream oss;
		oss <<
			"Cannot mix exchange components related to phase with exchange components related to kinetics, "
			<< this->formula;
		error_msg(oss.str().c_str(), CONTINUE);
		return;
	}

}
void
cxxSurfaceComp::multiply(LDBLE extensive)
{
	this->moles *= extensive;

	this->totals.multiply(extensive);

	this->charge_balance *= extensive;
}
void
cxxSurfaceComp::Serialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles)
{
	ints.push_back(dictionary.Find(this->formula));
	doubles.push_back(this->formula_z);
	doubles.push_back(this->moles);
	this->totals.Serialize(dictionary, ints, doubles);
	doubles.push_back(this->la);
	ints.push_back(dictionary.Find(this->charge_name));
	doubles.push_back(this->charge_balance);
	ints.push_back(dictionary.Find(this->phase_name));
	doubles.push_back(this->phase_proportion);
	ints.push_back(dictionary.Find(this->rate_name));
	doubles.push_back(this->Dw);
	ints.push_back(dictionary.Find(this->master_element));
}

void
cxxSurfaceComp::Deserialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles, int &ii, int &dd)
{
	this->formula = dictionary.GetWords()[ints[ii++]];
	this->formula_z = doubles[dd++];
	this->moles = doubles[dd++];
	this->totals.Deserialize(dictionary, ints, doubles, ii, dd);
	this->la = doubles[dd++];
	this->charge_name = dictionary.GetWords()[ints[ii++]];
	this->charge_balance = doubles[dd++];
	this->phase_name = dictionary.GetWords()[ints[ii++]];
	this->phase_proportion = doubles[dd++];
	this->rate_name = dictionary.GetWords()[ints[ii++]];
	this->Dw = doubles[dd++];
	this->master_element = dictionary.GetWords()[ints[ii++]];
}

const std::vector< std::string >::value_type temp_vopts[] = {
	std::vector< std::string >::value_type("formula"),	        // 0 
	std::vector< std::string >::value_type("moles"),	        // 1
	std::vector< std::string >::value_type("la"),	            // 2 
	std::vector< std::string >::value_type("charge_number"),	// 3 
	std::vector< std::string >::value_type("charge_balance"),	// 4
	std::vector< std::string >::value_type("phase_name"),	    // 5 
	std::vector< std::string >::value_type("rate_name"),	    // 6 
	std::vector< std::string >::value_type("phase_proportion"),	// 7 
	std::vector< std::string >::value_type("totals"),	        // 8
	std::vector< std::string >::value_type("formula_z"),	    // 9
	std::vector< std::string >::value_type("formula_totals"),	// 10
	std::vector< std::string >::value_type("dw"),	            // 11
	std::vector< std::string >::value_type("charge_name"),	    // 12
	std::vector< std::string >::value_type("master_element") 	// 13
};									   
const std::vector< std::string > cxxSurfaceComp::vopts(temp_vopts, temp_vopts + sizeof temp_vopts / sizeof temp_vopts[0]);	
