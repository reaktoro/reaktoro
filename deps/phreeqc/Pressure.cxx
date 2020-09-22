// Pressure.cxx: implementation of the cxxPressure class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <cassert>				// assert
#include <algorithm>			// std::sort

#include "Utils.h"				// define first
#include "Parser.h"
#include "Phreeqc.h"
#include "Pressure.h"
#include "phqalloc.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxPressure::cxxPressure(PHRQ_io *io)
	//
	// default constructor for cxxPressure 
	//
:	cxxNumKeyword(io)
{
	count = 0;
	equalIncrements = false;
}

cxxPressure::~cxxPressure()
{
}
int
cxxPressure::read(CParser & parser)
{
/*
 *      Reads pressure data for reaction steps
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
	// Number and description set in read_reaction_pressure

	PHRQ_io::LINE_TYPE lt;
	bool done = false;
	for (;;)
	{
		std::istream::pos_type ptr;
		std::istream::pos_type next_char = 0;
		std::string token, str;
		lt = parser.check_line(str, false, true, true, true);

		if (lt == PHRQ_io::LT_EMPTY || 
			lt == PHRQ_io::LT_KEYWORD ||
			lt == PHRQ_io::LT_EOF)
		{
			break;
		}
		if (lt == PHRQ_io::LT_OPTION)
		{
			this->error_msg("Expected numeric value for pressures.", PHRQ_io::OT_CONTINUE);
			break;
		}

		if (done)
		{
			this->error_msg("Unknown input following equal increment definition.", PHRQ_io::OT_CONTINUE);
			continue;
		}

		// LT_OK

		for (;;)
		{
			if (done) break;
			// new token
			std::string token;
			CParser::TOKEN_TYPE k =	parser.copy_token(token, next_char);

			// need new line
			if (k == CParser::TT_EMPTY)
			{
				break;
			}

			// read a pressure
			if (k == CParser::TT_DIGIT)
			{
				std::istringstream iss(token);
				LDBLE d;
				if (!(iss >> d))
				{
					this->error_msg("Expected numeric value for pressures.",
									 PHRQ_io::OT_CONTINUE);
				}
				else
				{
					this->pressures.push_back(d);
				}
				continue;
			}

			// non digit, must be "in"
			if (k == CParser::TT_UPPER || k == CParser::TT_LOWER)
			{
				if (this->pressures.size() != 2)
				{
					this->error_msg("To define equal increments, exactly two pressures should be defined.", CONTINUE);
				}
				else
				{
					int i = parser.copy_token(token, next_char);
					if (i == EMPTY)
					{
						error_msg("To define equal increments, define 'in n steps'.", CONTINUE);
					}
					else
					{
						std::istringstream iss(token);
						if ((iss >> i) && i > 0)
						{
							this->equalIncrements = true;
							this->count = i;
						}
						else
						{
							error_msg("Unknown input for pressure steps.", CONTINUE);
						}
					}
					done = true;
				}
				if (k == CParser::TT_UNKNOWN)
				{
					error_msg("Unknown input for pressure steps.", CONTINUE);
				}
			}
		} // tokens
	} // lines
	return lt;
}

void
cxxPressure::dump_raw(std::ostream & s_oss, unsigned int indent, int *n_out) const
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

	s_oss << indent0;
	int n_user_local = (n_out != NULL) ? *n_out : this->n_user;
	s_oss << "REACTION_PRESSURE_RAW        " << n_user_local << " " << this->description << "\n";

	s_oss << indent1;
	s_oss << "-count                     " << this->count << "\n";

	s_oss << indent1;
	s_oss << "-equal_increments          " << this->equalIncrements << "\n";

	// Temperature element and attributes

	s_oss << indent1;
	s_oss << "-pressures" << "\n";
	{
		int i = 0;
		s_oss << indent2;
		for (std::vector < LDBLE >::const_iterator it = this->pressures.begin();
			 it != this->pressures.end(); it++)
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
}

void
cxxPressure::read_raw(CParser & parser, bool check)
{
	// clear steps for modify operation, if pressures are read
	bool cleared_once = false;
	LDBLE d;
	CParser::TOKEN_TYPE k;

	std::istream::pos_type ptr;
	std::istream::pos_type next_char = 0;
	std::string token;
	int opt_save;
	bool useLastLine(false);

	// Number and description set in read_reaction_pressure_raw
	this->read_number_description(parser);

	opt_save = CParser::OPT_ERROR;
	bool equalIncrements_defined(false);
	bool count_defined(false);

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
			parser.error_msg("Unknown input in REACTION_PRESSURE_RAW keyword.",
							 PHRQ_io::OT_CONTINUE);
			parser.error_msg(parser.line().c_str(), PHRQ_io::OT_CONTINUE);
			useLastLine = false;
			break;

		case 0:				// pressures
			if (!cleared_once) 
			{
				this->pressures.clear();
				cleared_once = true;
			}
			while ((k =	parser.copy_token(token, next_char)) == CParser::TT_DIGIT)
			{
				std::istringstream iss(token);
				if (!(iss >> d))
				{
					parser.incr_input_error();
					parser.error_msg("Expected numeric value for pressures.",
									 PHRQ_io::OT_CONTINUE);
				}
				else
				{
					this->pressures.push_back(d);
				}
			}
			opt_save = 0;
			useLastLine = false;
			break;

		case 1:				// equal_increments
			if (!(parser.get_iss() >> this->equalIncrements))
			{
				this->equalIncrements = 0;
				parser.incr_input_error();
				parser.error_msg("Expected boolean value for equalIncrements.", PHRQ_io::OT_CONTINUE);
			}
			opt_save = CParser::OPT_DEFAULT;
			useLastLine = false;
			equalIncrements_defined = true;
			break;

		case 2:				// count
			if (!(parser.get_iss() >> this->count))
			{
				this->count = 0;
				parser.incr_input_error();
				parser.error_msg("Expected integer value for count.", PHRQ_io::OT_CONTINUE);
			}
			opt_save = CParser::OPT_DEFAULT;
			useLastLine = false;
			count_defined = true;
			break;
		}
		if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD)
			break;
	}
	// members that must be defined
	if (check)
	{
		if (equalIncrements_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Equal_increments not defined for REACTION_PRESSURE_RAW input.", 
				PHRQ_io::OT_CONTINUE);
		}
		if (count_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Count_temps not defined for REACTION_PRESSURE_RAW input.",
				 PHRQ_io::OT_CONTINUE);
		}
	}
}
#ifdef SKIP
void
cxxPressure::dump_xml(std::ostream & s_oss, unsigned int indent) const const
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

	// Temperature element and attributes
	s_oss << indent0;
	s_oss << "<temperature " << "\n";

	s_oss << indent1;
	s_oss << "pitzer_temperature_gammas=\"" << this->
		pitzer_temperature_gammas << "\"" << "\n";

	// components
	s_oss << indent1;
	s_oss << "<component " << "\n";
	for (std::list < cxxPressureComp >::const_iterator it =
		 temperatureComps.begin(); it != temperatureComps.end(); ++it)
	{
		it->dump_xml(s_oss, indent + 2);
	}

	return;
}
#endif
/* ---------------------------------------------------------------------- */
LDBLE cxxPressure::
Pressure_for_step(int step_number)
/* ---------------------------------------------------------------------- */
{
/*
 *   Determine pressure of reaction step
 */
	LDBLE p_temp;
	if (this->pressures.size() == 0) 
	{
		p_temp = 1;
	}
	else if (this->equalIncrements)
	{
		if (this->pressures.size() != 2)
		{
			error_msg("Number of pressures not equal to 2 for equal increments.", 0);
		}
		if (step_number > this->count)
		{
			p_temp = this->pressures[1];
		}
		else
		{
			LDBLE denom;
			denom = (this->count <= 1) ? 1 : (LDBLE) (this->count - 1);
			p_temp =  this->pressures[0] + ( this->pressures[1] - this->pressures[0]) *
				((LDBLE) (step_number - 1)) / (denom);
		}
	}
	else 
	{
		if (step_number > (int) this->pressures.size())
		{
			p_temp = this->pressures[this->pressures.size() - 1];
		}
		else
		{
			p_temp = this->pressures[step_number - 1];
		}

	}

	return (p_temp);
}
int cxxPressure::
Get_count(void) const
{
	if (equalIncrements) 
	{
		return this->count;
	}
	return (int) this->pressures.size();
}
void
cxxPressure::Serialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles)
{
	ints.push_back(this->n_user);
	{
		ints.push_back((int) this->pressures.size());
		for (size_t i = 0; i < this->pressures.size(); i++)
		{
			doubles.push_back(pressures[i]);
		}
	}
	ints.push_back(this->count);
	ints.push_back(this->equalIncrements ? 1 : 0);

}

void
cxxPressure::Deserialize(Dictionary & dictionary, std::vector < int >&ints, 
	std::vector < double >&doubles, int &ii, int &dd)
{
	this->n_user = ints[ii++];
	this->n_user_end = this->n_user;
	this->description = " ";

	{
		int count = ints[ii++];
		this->pressures.clear();
		for (int i = 0; i < count; i++)
		{
			this->pressures.push_back(doubles[dd++]);
		}
	}
	this->count = ints[ii++];
	this->equalIncrements = (ints[ii++] != 0);
}

const std::vector< std::string >::value_type temp_vopts[] = {
	std::vector< std::string >::value_type("pressures"),	        //0
	std::vector< std::string >::value_type("equal_increments"),	    //1
	std::vector< std::string >::value_type("count") 	            //2
};									   
const std::vector< std::string > cxxPressure::vopts(temp_vopts, temp_vopts + sizeof temp_vopts / sizeof temp_vopts[0]);	