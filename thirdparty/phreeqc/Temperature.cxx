// Temperature.cxx: implementation of the cxxTemperature class.
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
#include "Temperature.h"
#include "phqalloc.h"



//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxTemperature::cxxTemperature(PHRQ_io *io)
	//
	// default constructor for cxxTemperature 
	//
:	cxxNumKeyword(io)
{
	countTemps = 0;
	equalIncrements = false;
}
cxxTemperature::~cxxTemperature()
{
}

#ifdef SKIP
void
cxxTemperature::dump_xml(std::ostream & s_oss, unsigned int indent) const const
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
	for (std::list < cxxTemperatureComp >::const_iterator it =
		 temperatureComps.begin(); it != temperatureComps.end(); ++it)
	{
		it->dump_xml(s_oss, indent + 2);
	}

	return;
}
#endif

int
cxxTemperature::read(CParser & parser)
{
/*
 *      Reads temperature data for reaction steps
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
	// Number and description set in read_reaction_temperature

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
			this->error_msg("Expected numeric value for temperatures.", PHRQ_io::OT_CONTINUE);
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
					this->error_msg("Expected numeric value for temperatures.",
									 PHRQ_io::OT_CONTINUE);
				}
				else
				{
					this->temps.push_back(d);
				}
				continue;
			}

			// non digit, must be "in"
			if (k == CParser::TT_UPPER || k == CParser::TT_LOWER)
			{
				if (this->temps.size() != 2)
				{
					this->error_msg("To define equal increments, exactly two temperatures should be defined.", CONTINUE);
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
							this->countTemps = i;
						}
						else
						{
							error_msg("Unknown input for temperature steps.", CONTINUE);
						}
					}
					done = true;
				}
				if (k == CParser::TT_UNKNOWN)
				{
					error_msg("Unknown input for temperature steps.", CONTINUE);
				}
			}
		} // tokens
	} // lines
	return lt;
}

void
cxxTemperature::dump_raw(std::ostream & s_oss, unsigned int indent, int *n_out) const
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
	s_oss << "REACTION_TEMPERATURE_RAW     " << n_user_local << " " << this->description << "\n";

	s_oss << indent1;
	s_oss << "-count_temps               " << this->Get_countTemps() << "\n";

	s_oss << indent1;
	s_oss << "-equal_increments          " << this->equalIncrements << "\n";

	// Temperature element and attributes

	s_oss << indent1;
	s_oss << "-temps                     " << "\n";
	{
		int i = 0;
		s_oss << indent2;
		for (std::vector < LDBLE >::const_iterator it = this->temps.begin();
			 it != this->temps.end(); it++)
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
cxxTemperature::read_raw(CParser & parser, bool check)
{
	LDBLE d;
	CParser::TOKEN_TYPE k;
	// clear steps for modify operation, if pressures are read
	bool cleared_once = false;

	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	int opt_save;
	bool useLastLine(false);

	// Read temperature number and description
	this->read_number_description(parser);

	opt_save = CParser::OPT_ERROR;
	bool equalIncrements_defined(false);
	bool countTemps_defined(false);

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
			parser.error_msg("Unknown input in TEMPERATURE_COMP_RAW keyword.",
							 PHRQ_io::OT_CONTINUE);
			parser.error_msg(parser.line().c_str(), PHRQ_io::OT_CONTINUE);
			useLastLine = false;
			break;

		case 0:				// temps
			if (!cleared_once) 
			{
				this->temps.clear();
				cleared_once = true;
			}
			while ((k =
					parser.copy_token(token, next_char)) == CParser::TT_DIGIT)
			{
				std::istringstream iss(token);
				if (!(iss >> d))
				{
					parser.incr_input_error();
					parser.error_msg("Expected numeric value for temps.",
									 PHRQ_io::OT_CONTINUE);
				}
				else
				{
					this->temps.push_back(d);
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
				parser.
					error_msg("Expected boolean value for equalIncrements.",
							  PHRQ_io::OT_CONTINUE);
			}
			opt_save = CParser::OPT_DEFAULT;
			useLastLine = false;
			equalIncrements_defined = true;
			break;

		case 2:				// countTemps
			if (!(parser.get_iss() >> this->countTemps))
			{
				this->countTemps = 0;
				parser.incr_input_error();
				parser.error_msg("Expected integer value for countTemps.",
								 PHRQ_io::OT_CONTINUE);
			}
			opt_save = CParser::OPT_DEFAULT;
			useLastLine = false;
			countTemps_defined = true;
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
			parser.
				error_msg
				("Equal_increments not defined for REACTION_TEMPERATURE_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (countTemps_defined == false)
		{
			parser.incr_input_error();
			parser.
				error_msg
				("Count_temps not defined for REACTION_TEMPERATURE_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
	}
}
/* ---------------------------------------------------------------------- */
LDBLE cxxTemperature::
Temperature_for_step(int step_number)
/* ---------------------------------------------------------------------- */
{
/*
 *   Determine pressure of reaction step
 */
	LDBLE t_temp;
	if (this->temps.size() == 0) 
	{
		t_temp = 1;
	}
	else if (this->equalIncrements)
	{
		if (this->temps.size() != 2)
		{
			error_msg("Number of temperatures not equal to 2 for equal increments.", 0);
		}
		if (step_number > this->countTemps)
		{
			t_temp = this->temps[1];
		}
		else
		{
			LDBLE denom;
			denom = (this->countTemps <= 1) ? 1 : (LDBLE) (this->countTemps - 1);
			t_temp =  this->temps[0] + ( this->temps[1] - this->temps[0]) *
				((LDBLE) (step_number - 1)) / (denom);
		}
	}
	else 
	{
		if (step_number > (int) this->temps.size())
		{
			t_temp = this->temps[this->temps.size() - 1];
		}
		else
		{
			t_temp = this->temps[step_number - 1];
		}

	}

	return (t_temp);
}
int cxxTemperature::
Get_countTemps(void) const
{
	if (equalIncrements) 
	{
		return this->countTemps;
	}
	return (int) this->temps.size();
}
void
cxxTemperature::Serialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles)
{
	ints.push_back(this->n_user);
	{
		ints.push_back((int) this->temps.size());
		for (size_t i = 0; i < this->temps.size(); i++)
		{
			doubles.push_back(temps[i]);
		}
	}
	ints.push_back(this->countTemps);
	ints.push_back(this->equalIncrements ? 1 : 0);
}

void
cxxTemperature::Deserialize(Dictionary & dictionary, std::vector < int >&ints, 
	std::vector < double >&doubles, int &ii, int &dd)
{
	this->n_user = ints[ii++];
	this->n_user_end = this->n_user;
	this->description = " ";

	{
		int count = ints[ii++];
		this->temps.clear();
		for (int i = 0; i < count; i++)
		{
			this->temps.push_back(doubles[dd++]);
		}
	}
	this->countTemps = ints[ii++];
	this->equalIncrements = (ints[ii++] != 0);
}
const std::vector< std::string >::value_type temp_vopts[] = {
	std::vector< std::string >::value_type("temps"),	        //0
	std::vector< std::string >::value_type("equal_increments"),	//1
	std::vector< std::string >::value_type("count_temps") 	    //2
};									   
const std::vector< std::string > cxxTemperature::vopts(temp_vopts, temp_vopts + sizeof temp_vopts / sizeof temp_vopts[0]);	