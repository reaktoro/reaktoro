// Reaction.cxx: implementation of the cxxReaction class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <cassert>				// assert
#include <algorithm>			// std::sort

#include "Utils.h"				// define first
#include "Phreeqc.h"
#include "Reaction.h"
#include "phqalloc.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxReaction::cxxReaction(PHRQ_io *io)
	//
	// default constructor for cxxReaction 
	//
:	cxxNumKeyword(io)
{
	this->Set_units("Mol");
	countSteps = 0;
	equalIncrements = false;
	reactantList.type = cxxNameDouble::ND_NAME_COEF;
	elementList.type = cxxNameDouble::ND_ELT_MOLES;
}
cxxReaction::~cxxReaction()
{
}

#ifdef SKIP
void
cxxReaction::dump_xml(std::ostream & s_oss, unsigned int indent) const const
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

	// Reaction element and attributes
	s_oss << indent0;
	s_oss << "<irrev " << "\n";

	s_oss << indent1;
	s_oss << "pitzer_irrev_gammas=\"" << this->
		pitzer_irrev_gammas << "\"" << "\n";

	// components
	s_oss << indent1;
	s_oss << "<component " << "\n";
	for (std::list < cxxReactionComp >::const_iterator it =
		 irrevComps.begin(); it != irrevComps.end(); ++it)
	{
		it->dump_xml(s_oss, indent + 2);
	}

	return;
}
#endif

void
cxxReaction::dump_raw(std::ostream & s_oss, unsigned int indent, int *n_out) const
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

	// Reaction element and attributes
	s_oss << indent0;
	int n_user_local = (n_out != NULL) ? *n_out : this->n_user;
	s_oss << "REACTION_RAW                 " << n_user_local << " " << this->description << "\n";

	s_oss << indent1;
	s_oss << "-reactant_list" << "\n";
	this->reactantList.dump_raw(s_oss, indent + 2);

	s_oss << indent1;
	s_oss << "-steps" << "\n";
	{
		int i = 0;
		s_oss << indent2;
		for (std::vector < LDBLE >::const_iterator it = this->steps.begin();
			 it != this->steps.end(); it++)
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

	s_oss << indent1;
	s_oss << "-count_steps               " << this->countSteps << "\n";

	s_oss << indent1;
	s_oss << "-equal_increments          " << this->equalIncrements << "\n";

	s_oss << indent1;
	s_oss << "-units                     " << this->units << "\n";

	s_oss << indent1 << "# REACTION workspace variables #\n";
	s_oss << indent1;
	s_oss << "-element_list" << "\n";
	this->elementList.dump_raw(s_oss, indent + 2);
}

void
cxxReaction::read_raw(CParser & parser, const bool check)
{

	int j;
	LDBLE d;
	CParser::TOKEN_TYPE k;

	// clear steps for modify operation, if steps are read
	bool cleared_once = false;

	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	int opt_save;
	bool useLastLine(false);

	// Read irrev number and description
	this->read_number_description(parser);

	opt_save = CParser::OPT_ERROR;
	bool units_defined(false);
	bool equalIncrements_defined(false);
	bool countSteps_defined(false);

	for (;;)
	{
		int opt;
		if (useLastLine == false)
		{
			opt = parser.get_option(vopts, next_char);
		}
		else
		{
			opt = parser.getOptionFromLastLine(vopts, next_char, true);
		}
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
			opt = CParser::OPT_EOF;
			parser.error_msg("Unknown input in IRREV_COMP_RAW keyword.",
							 PHRQ_io::OT_CONTINUE);
			parser.error_msg(parser.line().c_str(), PHRQ_io::OT_CONTINUE);
			useLastLine = false;
			break;

		case 0:				// units
			j = parser.copy_token(token, next_char);
			if (j == CParser::TT_EMPTY)
				break;
			this->Set_units(token.c_str());
			opt_save = CParser::OPT_DEFAULT;
			useLastLine = false;
			units_defined = true;
			break;

		case 1:				// reactant_list
			if (this->reactantList.read_raw(parser, next_char) !=
				CParser::PARSER_OK)
			{
				parser.incr_input_error();
				parser.error_msg("Expected reactant formula and coefficient.",
								 PHRQ_io::OT_CONTINUE);
			}
			opt_save = 1;
			useLastLine = false;
			break;

		case 2:				// element_list
			if (this->elementList.read_raw(parser, next_char) !=
				CParser::PARSER_OK)
			{
				parser.incr_input_error();
				parser.error_msg("Expected element formula and coefficient.",
								 PHRQ_io::OT_CONTINUE);
			}
			opt_save = 2;
			useLastLine = false;
			break;

		case 3:				// steps
			if (!cleared_once) 
			{
				this->steps.clear();
				cleared_once = true;
			}
			while ((k =
					parser.copy_token(token, next_char)) == CParser::TT_DIGIT)
			{
				std::istringstream iss(token);
				if (!(iss >> d))
				{
					parser.incr_input_error();
					parser.error_msg("Expected numeric value for steps.",
									 PHRQ_io::OT_CONTINUE);
				}
				else
				{
					this->steps.push_back(d);
				}
			}
			opt_save = 3;
			useLastLine = false;
			break;

		case 4:				// equal_increments
			if (!(parser.get_iss() >> this->equalIncrements))
			{
				this->equalIncrements = 0;
				parser.incr_input_error();
				parser.
					error_msg("Expected boolean value for equalIncrements.",
							  PHRQ_io::OT_CONTINUE);
			}
			opt_save = CParser::OPT_DEFAULT;
			useLastLine = false;
			equalIncrements_defined = true;
			break;

		case 5:				// countSteps
			if (!(parser.get_iss() >> this->countSteps))
			{
				this->countSteps = 0;
				parser.incr_input_error();
				parser.error_msg("Expected integer value for countSteps.",
								 PHRQ_io::OT_CONTINUE);
			}
			opt_save = CParser::OPT_DEFAULT;
			useLastLine = false;
			countSteps_defined = true;
			break;
		}
		if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD)
			break;
	}
	if (check)
	{
		// members that must be defined
		if (units_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Units not defined for REACTION_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (equalIncrements_defined == false)
		{
			parser.incr_input_error();
			parser.
				error_msg("Equal_increments not defined for REACTION_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (countSteps_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Count_steps not defined for REACTION_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
	}
}
int cxxReaction::
Get_reaction_steps(void) const
{
	if (equalIncrements) 
	{
		return this->countSteps;
	}
	return (int) this->steps.size();
}
const std::vector< std::string >::value_type temp_vopts[] = {
	std::vector< std::string >::value_type("units"),	        //0
	std::vector< std::string >::value_type("reactant_list"),	//1
	std::vector< std::string >::value_type("element_list"),	    //2
	std::vector< std::string >::value_type("steps"),	        //3
	std::vector< std::string >::value_type("equal_increments"),	//4
	std::vector< std::string >::value_type("count_steps") 	    //5
};									   
const std::vector< std::string > cxxReaction::vopts(temp_vopts, temp_vopts + sizeof temp_vopts / sizeof temp_vopts[0]);	