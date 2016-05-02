// GasComp.cxx: implementation of the cxxGasComp class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <iostream>				// std::cout std::cerr
#include <cassert>				// assert
#include <algorithm>			// std::sort
#include <float.h>

#include "Utils.h"				// define first
#include "Phreeqc.h"
#include "GasComp.h"
#include "phqalloc.h"
#include "Dictionary.h"



//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxGasComp::cxxGasComp(PHRQ_io *io)
	//
	// default constructor for cxxExchComp 
	//
	: PHRQ_base(io)
{
	p_read = 0.0;
	moles = 0.0;
	initial_moles = 0.0;
}

cxxGasComp::~cxxGasComp()
{
}

void
cxxGasComp::dump_raw(std::ostream & s_oss, unsigned int indent) const
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
	s_oss << indent0 << "# GAS_PHASE_MODIFY candidate identifiers #\n";
	s_oss << indent0 << "-moles                   " << this->moles << "\n";

	s_oss << indent0 << "# GAS_PHASE_MODIFY candidate identifiers with new_def=true #\n";
	s_oss << indent0 << "-p_read                  " << this->p_read << "\n";

	s_oss << indent0 << "# GasComp workspace variables #\n";
	s_oss << indent0 << "-initial_moles           " << this->initial_moles << "\n";
}

bool
cxxGasComp::read_raw(CParser & parser, bool check)
{
	std::string str;

	int errors = parser.get_input_error();

	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;

	bool moles_defined(false);
	int opt;
	for (;;)
	{
		opt = parser.get_option(vopts, next_char);

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

		case 0:				// phase_name
			output_msg("-phase_name is obsolete. Define with -component\n");
			break;

		case 1:				// name
			output_msg("-name is obsolete. Define with -component\n");
			break;

		case 2:				// p_read
			if (!(parser.get_iss() >> this->p_read))
			{
				this->p_read = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for initial partial pressure.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;

		case 3:				// moles
			if (!(parser.get_iss() >> this->moles))
			{
				this->moles = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for moles.",
								 PHRQ_io::OT_CONTINUE);
			}
			moles_defined = true;
			break;

		case 4:				// initial_moles
			if (!(parser.get_iss() >> this->initial_moles))
			{
				this->initial_moles = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for initial_moles.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;
		}
		if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD)
			break;
	}
	if (check)
	{
		// members that must be defined
		if (moles_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Moles not defined for GasComp input.",
				PHRQ_io::OT_CONTINUE);
		}
	}
	return (parser.get_input_error() == errors);
}
void
cxxGasComp::add(const cxxGasComp & addee, LDBLE extensive)
{
	//LDBLE ext1, ext2;
	if (extensive == 0.0)
		return;
	if (addee.phase_name.size() == 0)
		return;

	/*
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
	*/

	assert(this->phase_name == addee.phase_name);

	this->p_read += addee.p_read * extensive;
	this->moles += addee.moles * extensive;
	this->initial_moles += addee.initial_moles * extensive;
}

void
cxxGasComp::multiply(LDBLE extensive)
{
	this->p_read *= extensive;
	this->moles *= extensive;
	this->initial_moles *= extensive;
}
void
cxxGasComp::Serialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles)
{
	ints.push_back(dictionary.Find(this->phase_name));
	doubles.push_back(this->moles);
	doubles.push_back(this->p_read);
	doubles.push_back(this->initial_moles);
}

void
cxxGasComp::Deserialize(Dictionary & dictionary, std::vector < int >&ints, 
	std::vector < double >&doubles, int &ii, int &dd)
{
	this->phase_name = dictionary.GetWords()[ints[ii++]];
	this->moles = doubles[dd++];
	this->p_read = doubles[dd++];
	this->initial_moles = doubles[dd++];
}

const std::vector< std::string >::value_type temp_vopts[] = {
	std::vector< std::string >::value_type("phase_name"),	// 0 
	std::vector< std::string >::value_type("name"),	        // 1
	std::vector< std::string >::value_type("p_read"),	    // 2 
	std::vector< std::string >::value_type("moles"),	    // 3 
	std::vector< std::string >::value_type("initial_moles")	// 4 
};									   
const std::vector< std::string > cxxGasComp::vopts(temp_vopts, temp_vopts + sizeof temp_vopts / sizeof temp_vopts[0]);	
