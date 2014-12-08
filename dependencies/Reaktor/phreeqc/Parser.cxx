// Parser.cpp: implementation of the CParser class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif

#include <algorithm>			// std::transform
#include <map>					// std::map
#include <cassert>				// assert
#include <stdlib.h>
#include <iostream>				// std::cout std::cerr
#include "Utils.h"
#include <stdio.h>
#include "Parser.h"
#include "PHRQ_io.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
CParser::CParser(PHRQ_io *io):
PHRQ_base(io),
m_input_stream(std::cin), 
m_input_error(0),
m_next_keyword(Keywords::KEY_NONE)
{
	if (!io)
	{
		error_msg("This parser constructor requires non-null phrq_io", PHRQ_io::OT_STOP);
	}
	else
	{
		m_line_save = io->Get_m_line();
		m_line = io->Get_m_line();
		m_line_type = io->Get_m_line_type();
		m_line_iss.str(m_line);
		m_line_iss.seekg(0, std::ios_base::beg);
		m_line_iss.clear();
		echo_file = EO_ALL;
		echo_stream = EO_NONE;
		accumulate = false;
		phrq_io_only = true;
	}
}

CParser::CParser(std::istream & input, PHRQ_io *io):
PHRQ_base(io),
m_input_stream(input), 
m_input_error(0),
m_next_keyword(Keywords::KEY_NONE)
{
	m_line_save.reserve(80);
	m_line.reserve(80);
	echo_file = EO_ALL;
	echo_stream = EO_NONE;
	accumulate = false;
	phrq_io_only = false;
}

CParser::~CParser()
{
}

PHRQ_io::LINE_TYPE CParser::check_line(const std::string & str,
									   bool allow_empty, bool allow_eof,
									   bool allow_keyword, bool print)
{
	PHRQ_io::LINE_TYPE
		i;

	// Get line
	do
	{
		i = get_line();
		// reset iss
		m_line_iss.str(m_line);
		m_line_iss.seekg(0, std::ios_base::beg);
		m_line_iss.clear();

		// output for stream
		switch (this->echo_stream)
		{
		case EO_NONE:
			break;
		case EO_ALL:
			if (i != PHRQ_io::LT_EOF)
			{
				std::ostringstream msg;
				msg << "\t" << m_line_save << "\n";
				//get_output() << msg;
				io->output_msg(msg.str().c_str());
			}
			break;
		case EO_KEYWORDS:
			if (i == PHRQ_io::LT_KEYWORD)
			{
				std::ostringstream msg;
				msg << "\t" << m_line_save << "\n";
				//get_output() << msg;
				io->output_msg(msg.str().c_str());
			}
			break;
		case EO_NOKEYWORDS:
			if (i != PHRQ_io::LT_KEYWORD && i != PHRQ_io::LT_EOF)
			{
				std::ostringstream msg;
				msg << "\t" << m_line_save << "\n";
				//get_output() << msg;
				io->output_msg(msg.str().c_str());
			}
			break;
		}
		// output for file
		switch (this->echo_file)
		{
		case EO_NONE:
			break;
		case EO_ALL:
			if (i != PHRQ_io::LT_EOF)
			{
				std::ostringstream msg;
				msg << "\t" << m_line_save << "\n";
				this->echo_msg(msg.str());
			}
			break;
		case EO_KEYWORDS:
			if (i == PHRQ_io::LT_KEYWORD)
			{
				std::ostringstream msg;
				msg << "\t" << m_line_save << "\n";
				this->echo_msg(msg.str());
			}
			break;

		case EO_NOKEYWORDS:
			if (i != PHRQ_io::LT_KEYWORD && i != PHRQ_io::LT_EOF)
			{
				std::ostringstream msg;
				msg << "\t" << m_line_save << "\n";
				this->echo_msg(msg.str().c_str());
			}
			break;
		}
	}
	while (i == PHRQ_io::LT_EMPTY && allow_empty == false);

	// Check eof
	if (i == PHRQ_io::LT_EOF && allow_eof == false)
	{
		std::ostringstream msg;
		msg << "Unexpected eof while reading " << str <<
			"\nExecution terminated.\n";
		error_msg(msg.str().c_str(), PHRQ_io::OT_STOP);
	}

	// Check keyword
	if (i == PHRQ_io::LT_KEYWORD && allow_keyword == false)
	{
		std::ostringstream msg;
		msg << "Expected data for " << str <<
			", but got a keyword ending data block.";
		error_msg(msg.str().c_str(), PHRQ_io::OT_CONTINUE);
		incr_input_error();
	}
	m_line_type = i;
	return i;
}
PHRQ_io::LINE_TYPE CParser::get_line_phrq_io()
{
	m_line_type = io->get_line();
	m_line_save = io->Get_m_line_save();
	m_line = io->Get_m_line();
	m_next_keyword = io->Get_m_next_keyword();
	if (accumulate)
	{
		this->accumulated.append(m_line_save);
		this->accumulated.append("\n");
	}
	return m_line_type;
}

PHRQ_io::LINE_TYPE CParser::get_line()
{
	if (this->phrq_io_only)
	{
		return get_line_phrq_io();
	}
	PHRQ_io::LINE_TYPE return_value = PHRQ_io::LT_EMPTY;
	while (return_value == PHRQ_io::LT_EMPTY)
	{
		//
		// Eliminate all characters after # sign as a comment
		//

		//
		// Get line, check for eof
		//
		if (get_logical_line() == PHRQ_io::LT_EOF)
		{
			if (!m_input_stream.eof())
			{
				error_msg("Reading input file.", PHRQ_io::OT_CONTINUE);
				error_msg("istream::get() returned an error.", PHRQ_io::OT_STOP);
			}
			else
			{
				//{{MOD
				m_line.erase(m_line.begin(), m_line.end());	// m_line.clear();
				//}}MOD
				m_next_keyword = Keywords::KEY_END;
				return PHRQ_io::LT_EOF;
			}
		}

		//
		// Get long lines
		//
		bool
			empty = true;
		m_line = m_line_save.substr(0, m_line_save.find_first_of('#'));
		for (unsigned int i = 0; i < m_line.size(); ++i)
		{
			if (!::isspace(m_line[i]))
			{
				empty = false;
				break;
			}
		}

		if (this->accumulate)
		{
			this->accumulated.append(m_line_save);
			this->accumulated.append("\n");
		}
		//
		// New line character encountered 
		//
		return_value = (empty ? PHRQ_io::LT_EMPTY : PHRQ_io::LT_OK);
	}

	//
	// Determine return_value
	// 
	if (return_value == PHRQ_io::LT_OK)
	{
		if (check_key(m_line.begin(), m_line.end()))
		{
			return_value = PHRQ_io::LT_KEYWORD;
		}
		else
		{
			std::string::iterator beg = m_line.begin();
			std::string::iterator end = m_line.end();
			std::string token;
			copy_token(token, beg, end);

			if (token.size() > 1 && token[0] == '-' &&::isalpha(token[1]))
			{
				return_value = PHRQ_io::LT_OPTION;
			}
		}
	}
	return return_value;
}

/**
        Reads input stream until end of line, ";", or eof
        stores characters in line_save

        returns:
                EOF on empty line on end of file or
                OK otherwise
*/
PHRQ_io::LINE_TYPE CParser::get_logical_line()
{
	int
		j;
	unsigned int
		pos;
	char
		c;

	m_line_save.erase(m_line_save.begin(), m_line_save.end());	// m_line_save.clear();

	while ((j = m_input_stream.get()) != std::char_traits < char >::eof())
	{
		c = (char) j;
		if (c == '#')
		{
			// ignore all chars after # until newline
			do
			{
				c = (char) j;
				if (c == '\n')
				{
					break;
				}
				m_line_save += c;
			}
			while ((j =
					m_input_stream.get()) != std::char_traits <
				   char >::eof());
		}
		if (c == ';')
			break;
		if (c == '\n')
		{
			break;
		}
		if (c == '\\')
		{
			pos = (int) m_line_save.size();
			m_line_save += c;
			while ((j =
					m_input_stream.get()) != std::char_traits < char >::eof())
			{
				c = (char) j;
				if (c == '\\')
				{
					pos = (int) m_line_save.size();
					m_line_save += c;
					continue;
				}
				if (c == '\n')
				{
					// remove '\\'
					for (; pos < m_line_save.size(); pos++)
					{
						m_line_save[pos] = m_line_save[pos + 1];
					}
					m_line_save.erase(m_line_save.size() - 1, 1);
					break;
				}
				m_line_save += c;
				if (!::isspace(j))
					break;
			}
		}
		else
		{
			m_line_save += c;
		}
	}
	if (j == std::char_traits < char >::eof() && m_line_save.size() == 0)
	{
		return (PHRQ_io::LT_EOF);
	}
	return (PHRQ_io::LT_OK);
}

bool
CParser::check_key(std::string::iterator begin, std::string::iterator end)
{
	std::string lowercase;
	copy_token(lowercase, begin, end);
	std::transform(lowercase.begin(), lowercase.end(), lowercase.begin(),
				   tolower);
	m_next_keyword = Keywords::Keyword_search(lowercase);
	if (m_next_keyword == Keywords::KEY_NONE)
	{
		return false;
	}
	return true;
}
CParser::STATUS_TYPE CParser::check_units(std::string & tot_units,
										  bool alkalinity,
										  bool check_compatibility,
										  const std::string & default_units,
										  bool print)
{
/*
 *   Check if legitimate units
 *   Input:
 *           tot_units           character string to check,
 *           alkalinity          true if alkalinity, false if any other total,
 *           check_compatibility true check alk and default units, false otherwise
 *           default_units       character string of default units (check /L, /kg, etc)
 *           print               true print warning messages
 *   Output:
 *           tot_units           standard form for unit
 */
	using
		Utilities::str_tolower;
	using
		Utilities::replace;
	using
		Utilities::squeeze_white;

	static const char *
		units[] = {
		"Mol/l",				/* 0 */
		"mMol/l",				/* 1 */
		"uMol/l",				/* 2 */
		"g/l",					/* 3 */
		"mg/l",					/* 4 */
		"ug/l",					/* 5 */
		"Mol/kgs",				/* 6 */
		"mMol/kgs",				/* 7 */
		"uMol/kgs",				/* 8 */
		"g/kgs",				/* 9 = ppt */
		"mg/kgs",				/* 10 = ppm */
		"ug/kgs",				/* 11 = ppb */
		"Mol/kgw",				/* 12 = mol/kg H2O */
		"mMol/kgw",				/* 13 = mmol/kg H2O */
		"uMol/kgw",				/* 14 = umol/kg H2O */
		"g/kgw",				/* 15 = mol/kg H2O */
		"mg/kgw",				/* 16 = mmol/kg H2O */
		"ug/kgw",				/* 17 = umol/kg H2O */
		"eq/l",					/* 18 */
		"meq/l",				/* 19 */
		"ueq/l",				/* 20 */
		"eq/kgs",				/* 21 */
		"meq/kgs",				/* 22 */
		"ueq/kgs",				/* 23 */
		"eq/kgw",				/* 24 */
		"meq/kgw",				/* 25 */
		"ueq/kgw",				/* 26 */
	};

	squeeze_white(tot_units);
	str_tolower(tot_units);
	replace("milli", "m", tot_units);
	replace("micro", "u", tot_units);
	replace("grams", "g", tot_units);
	replace("gram", "g", tot_units);
	replace("moles", "Mol", tot_units);
	replace("mole", "Mol", tot_units);
	replace("mol", "Mol", tot_units);
	replace("liter", "l", tot_units);
	replace("kgh", "kgw", tot_units);
	replace("ppt", "g/kgs", tot_units);
	replace("ppm", "mg/kgs", tot_units);
	replace("ppb", "ug/kgs", tot_units);
	replace("equivalents", "eq", tot_units);
	replace("equivalent", "eq", tot_units);
	replace("equiv", "eq", tot_units);

	std::string::size_type end;
	if ((end = tot_units.find("/l")) != std::string::npos)
	{
		tot_units.resize(end + 2);
	}
	if ((end = tot_units.find("/kgs")) != std::string::npos)
	{
		tot_units.resize(end + 4);
	}
	if ((end = tot_units.find("/kgw")) != std::string::npos)
	{
		tot_units.resize(end + 4);
	}

	//
	//  Check if unit in list
	//
	bool
		found = false;
	for (unsigned int i = 0; i < sizeof(units) / sizeof(char *); ++i)
	{
		if (tot_units.compare(units[i]) == 0)
		{
			found = true;
			break;
		}
	}
	if (!found)
	{
		if (print)
		{
			std::ostringstream err;
			err << "Unknown unit, " << tot_units;
			this->error_msg(err.str().c_str(), PHRQ_io::OT_CONTINUE);
		}
		return PARSER_ERROR;
	}

	//
	//   Check if units are compatible with default_units
	//
	if (check_compatibility == false)
		return PARSER_OK;

	//
	//   Special cases for alkalinity
	//
	if (alkalinity == true && tot_units.find("Mol") != std::string::npos)
	{
		if (print)
		{
			this->warning_msg
				("Alkalinity given in moles, assumed to be equivalents.");
		}
		replace("Mol", "eq", tot_units);
	}
	if (alkalinity == false && tot_units.find("eq") != std::string::npos)
	{
		if (print)
		{
			this->error_msg("Only alkalinity can be entered in equivalents.",
					  PHRQ_io::OT_CONTINUE);
		}
		return PARSER_ERROR;
	}

	//
	//  See if default_units are compatible with tot_units
	//
	if (default_units.find("/l") != std::string::npos
		&& tot_units.find("/l") != std::string::npos)
		return PARSER_OK;
	if (default_units.find("/kgs") != std::string::npos
		&& tot_units.find("/kgs") != std::string::npos)
		return PARSER_OK;
	if (default_units.find("/kgw") != std::string::npos
		&& tot_units.find("/kgw") != std::string::npos)
		return PARSER_OK;

	std::string str = default_units;
	replace("kgs", "kg solution", str);
	replace("kgs", "kg solution", tot_units);
	replace("kgw", "kg water", str);
	replace("kgw", "kg water", tot_units);
	replace("/l", "/L", str);
	replace("Mol", "mol", str);
	replace("/l", "/L", tot_units);
	replace("Mol", "mol", tot_units);

	if (print)
	{
		std::ostringstream err;
		err << "Units for master species, " << tot_units <<
			", are not compatible with default units, " << str << ".";
		this->error_msg(err.str().c_str(), PHRQ_io::OT_CONTINUE);
	}
	return PARSER_ERROR;
}
CParser::TOKEN_TYPE CParser::token_type(const std::string & token)
{
	if (!token.empty())
	{
		if (::isupper(token[0]))
		{
			return CParser::TT_UPPER;
		}
		else if (::islower(token[0]))
		{
			return CParser::TT_LOWER;
		}
		else if (::isdigit(token[0]) || token[0] == '.' || token[0] == '-')
		{
			return CParser::TT_DIGIT;
		}
		else
		{
			assert(!::isspace(token[0]));
			return CParser::TT_UNKNOWN;
		}
	}
	else
	{
		return CParser::TT_EMPTY;
	}
}

CParser::TOKEN_TYPE CParser::peek_token()
{
	std::istringstream::pos_type pos = m_line_iss.tellg();
	std::string token;
	m_line_iss >> token;
	m_line_iss.seekg(pos);
	return token_type(token);
}

CParser::TOKEN_TYPE CParser::copy_token(std::string & token,
										std::string::iterator & begin,
										std::string::iterator & end)
{
	if (begin != end)
	{
		std::string::iterator b = begin;
		for (; b < end &&::isspace(*b); ++b);

		begin = b;
		for (; begin < end && !::isspace(*begin); ++begin);

		token.assign(b, begin);
	}
	else
	{
		token.resize(0);
	}

	return token_type(token);
}

CParser::TOKEN_TYPE CParser::copy_token(std::string & token,
										std::istream & is)
{
	is >> token;
	return token_type(token);
}

CParser::TOKEN_TYPE CParser::copy_token(std::string & token,
										std::istream::pos_type & pos)
{
	m_line_iss.seekg(pos);
	if (!(m_line_iss >> token))
	{
		token.erase(token.begin(), token.end());	// token.clear();
	}
	pos = m_line_iss.tellg();
	return token_type(token);
}

CParser::FIND_TYPE CParser::find_option(const std::string & item, int *n,
										const std::vector < std::string >
										&list, bool exact)
{
	std::string token(item);
	std::transform(token.begin(), token.end(), token.begin(), tolower);
	for (unsigned int i = 0; i < list.size(); i++)
	{
		if (exact == true)
		{
			if (list[i].compare(token) == 0)
			{
				*n = i;
				return FT_OK;
			}
		}
		else
		{
			if (list[i].find(token) == 0)
			{
				*n = i;
				return FT_OK;
			}
		}
	}

	*n = -1;
	return FT_ERROR;
}

int
CParser::get_option(const std::vector < std::string > &opt_list,
					std::string::iterator & next_char)
{
	//
	// Read a line and check for options
	//
	int j;
	int /*  opt_l,  */ opt;
	std::string::iterator opt_ptr;
	std::string option;

#if !defined(R_SO)
	fprintf(stderr, "Did not think this get_option was called\n");
#endif
	//
	// Read line
	//
	PHRQ_io::LINE_TYPE lt = check_line("get_option", false, true, true, true);
	if (lt == PHRQ_io::LT_EOF)
	{
		j = OPT_EOF;
	}
	else if (lt == PHRQ_io::LT_KEYWORD)
	{
		j = OPT_KEYWORD;
	}
	else if (lt == PHRQ_io::LT_OPTION)
	{
		opt_ptr = m_line.begin();
		std::string::iterator end = m_line.end();
		copy_token(option, opt_ptr, end);
		if (find_option(option, &opt, opt_list, false) == CParser::FT_OK)
		{
			j = opt;
			m_line_save.replace(m_line_save.find(option), option.size(),
								opt_list[opt]);
			m_line.replace(m_line.find(option), option.size(), opt_list[opt]);
			opt_ptr = m_line.begin();
			std::string::iterator end = m_line.end();
			copy_token(option, opt_ptr, end);
			next_char = opt_ptr;
			if (true)			// pr.echo_input == TRUE
			{
				if (true)		// database_file == NULL
				{
					//get_output() << "\t" << m_line_save << "\n";
					std::ostringstream msg;
					msg << "\t" << m_line_save << "\n";
					io->output_msg(msg.str().c_str());
				}
			}
		}
		else
		{
			if (true)			// (database_file == NULL)
			{
				//get_output() << "\t" << m_line_save << "\n";
				std::ostringstream msg;
				msg << "\t" << m_line_save << "\n";
				io->output_msg(msg.str().c_str());
			}

			std::ostringstream err;
			err << "Unknown option." << "\n";
			err << m_line_save << "\n";
			error_msg(err.str().c_str());

			j = OPT_ERROR;
			next_char = m_line.begin();
		}
	}
	else
	{
		opt_ptr = m_line.begin();
		std::string::iterator end = m_line.end();
		copy_token(option, opt_ptr, end);
		if (find_option(option, &opt, opt_list, true) == FT_OK)
		{
			j = opt;
			next_char = opt_ptr;
		}
		else
		{
			j = OPT_DEFAULT;
			next_char = m_line.begin();
		}
		if (true)				// pr.echo_input == TRUE
		{
			if (true)			// database_file == NULL
			{
#if !defined(R_SO)
				std::cout << "\t" << m_line_save << "\n";
#endif
			}
		}
	}
	return (j);
}

int
CParser::get_option(const std::vector < std::string > &opt_list,
					std::istream::pos_type & next_pos)
{
	//
	// Read a line and check for options
	//
	int j;
	int opt;
	std::istream::pos_type pos_ptr;
	std::string option;

	//
	// Read line
	//
	PHRQ_io::LINE_TYPE lt = check_line("get_option", false, true, true, true);
	if (lt == PHRQ_io::LT_EOF)
	{
		j = OPT_EOF;
	}
	else if (lt == PHRQ_io::LT_KEYWORD)
	{
		j = OPT_KEYWORD;
	}
	else if (lt == PHRQ_io::LT_OPTION)
	{
		std::string::iterator opt_ptr = m_line.begin();
		std::string::iterator end = m_line.end();
		copy_token(option, opt_ptr, end);
		if (find_option(option.substr(1), &opt, opt_list, false) == FT_OK)
		{
			// replace -option with option
			j = opt;
			m_line_save.replace(m_line_save.find(option), option.size(),
								opt_list[opt]);
			m_line.replace(m_line.find(option), option.size(), opt_list[opt]);

			// reset iss
			m_line_iss.str(m_line);
			m_line_iss.seekg(0, std::ios_base::beg);
			m_line_iss.clear();

			pos_ptr = 0;
			copy_token(option, pos_ptr);
			next_pos = pos_ptr;
		}
		else
		{
			j = OPT_ERROR;
			next_pos = pos_ptr;
		}
	}
	else
	{
		pos_ptr = m_line_iss.tellg();
		m_line_iss >> option;
		if (find_option(option, &opt, opt_list, true) == FT_OK)
		{
			j = opt;
			next_pos = m_line_iss.tellg();
		}
		else
		{
			j = OPT_DEFAULT;
			m_line_iss.seekg(pos_ptr);
			m_line_iss.clear();
			next_pos = pos_ptr;
		}
	}
	return (j);
}
CParser::STATUS_TYPE CParser::get_elt(std::string::iterator & begin,
									  const std::string::iterator end,
									  std::string & element)
{
	element.erase(element.begin(), element.end());	// element.clear();

	if (begin == end)
	{
		error_msg("Empty string in get_elt.  Expected an element name.",
				  PHRQ_io::OT_CONTINUE);
		return PARSER_ERROR;
	}

	//
	// Load name into char array element
	//
	char
		c = *begin;
	++begin;
	element.insert(element.end(), c);	// element.push_back(c);
	if (c == '[')
	{
		while ((c = *begin) != ']')
		{
			element.insert(element.end(), c);	// element.push_back(c);
			++begin;
			if ((c = *begin) == ']')
			{
				element.insert(element.end(), c);	// element.push_back(c);
				++begin;
				break;
			}
			else if (begin == end)
			{
				error_msg("No ending bracket (]) for element name",
						  PHRQ_io::OT_CONTINUE);
				incr_input_error();
				return PARSER_ERROR;
			}
		}
		while (::islower(c = *begin) || c == '_')
		{
			element.insert(element.end(), c);	// element.push_back(c);
			++begin;
			if (begin == end)
				break;
		}
	}
	else
	{
		while (::islower(c = *begin) || c == '_')
		{
			element.insert(element.end(), c);	// element.push_back(c);
			++begin;
			if (begin == end)
				break;
		}
	}
	return PARSER_OK;
}

CParser::STATUS_TYPE CParser::parse_couple(std::string & token)
{
	// Parse couple puts redox couples in standard form
	// "+" is removed and couples are rewritten in sort
	// order.

	if (Utilities::strcmp_nocase_arg1(token.c_str(), "pe") == 0)
	{
		Utilities::str_tolower(token);
		return PARSER_OK;
	}

	while (Utilities::replace("+", "", token));

	std::string::iterator ptr = token.begin();
	std::string elt1;
	get_elt(ptr, token.end(), elt1);

	if (*ptr != '(')
	{
		std::ostringstream err_msg;
		err_msg << "Element name must be followed by " <<
			"parentheses in redox couple, " << token << ".";
		error_msg(err_msg.str().c_str(), PHRQ_io::OT_CONTINUE);
		incr_input_error();
		return PARSER_ERROR;
	}

	int
		paren_count = 1;
	std::string paren1 = "(";
	while (ptr != token.end())
	{
		++ptr;
		if (*ptr == '/' || ptr == token.end())
		{
			std::ostringstream err_msg;
			err_msg << "End of line or  " "/"
				" encountered before end of parentheses, " << token << ".";
			error_msg(err_msg.str().c_str(), PHRQ_io::OT_CONTINUE);
			return PARSER_ERROR;
		}
		paren1.insert(paren1.end(), *ptr);	// element.push_back(c);
		if (*ptr == '(')
			++paren_count;
		if (*ptr == ')')
			--paren_count;
		if (paren_count == 0)
			break;
	}

	++ptr;
	if (ptr == token.end() || *ptr != '/')
	{
		std::ostringstream err_msg;
		err_msg << " " "/" " must follow parentheses " <<
			"ending first half of redox couple, " << token << ".";
		error_msg(err_msg.str().c_str(), PHRQ_io::OT_CONTINUE);
		return PARSER_ERROR;
	}
	++ptr;
	std::string elt2;
	get_elt(ptr, token.end(), elt2);
	if (elt1.compare(elt2) != 0)
	{
		std::ostringstream err_msg;
		err_msg << "Redox couple must be two redox states " <<
			"of the same element, " << token << ".";
		error_msg(err_msg.str().c_str(), PHRQ_io::OT_CONTINUE);
		return PARSER_ERROR;
	}
	if (*ptr != '(')
	{
		std::ostringstream err_msg;
		err_msg << "Element name must be followed by "
			"parentheses in redox couple, " << token << ".";
		error_msg(err_msg.str().c_str(), PHRQ_io::OT_CONTINUE);
		incr_input_error();
		return PARSER_ERROR;
	}
	std::string paren2 = "(";
	paren_count = 1;

	while (ptr != token.end())
	{
		++ptr;
		if (*ptr == '/' || ptr == token.end())
		{
			std::ostringstream err_msg;
			err_msg << "End of line or  " "/"
				" encountered before end of parentheses, " << token << ".";
			error_msg(err_msg.str().c_str(), PHRQ_io::OT_CONTINUE);
			return PARSER_ERROR;
		}
		paren2.insert(paren2.end(), *ptr);	// element.push_back(c);
		if (*ptr == '(')
			++paren_count;
		if (*ptr == ')')
			--paren_count;
		if (paren_count == 0)
			break;
	}
	if (paren1.compare(paren2) < 0)
	{
		token = elt1 + paren1 + std::string("/") + elt2 + paren2;
	}
	else if (paren1.compare(paren2) > 0)
	{
		token = elt2 + paren2 + std::string("/") + elt1 + paren1;
	}
	else
	{
		std::ostringstream err_msg;
		err_msg << "Both parts of redox couple are the same, " <<
			token << ".";
		error_msg(err_msg.str().c_str(), PHRQ_io::OT_CONTINUE);
		return PARSER_ERROR;
	}
	return PARSER_OK;
}

template <class T>
CParser::STATUS_TYPE CParser::addPair(std::map < std::string, T >&totals,
									  std::istream::pos_type & pos)
{
	std::string token;
	T d;
	CParser::TOKEN_TYPE j;

	m_line_iss.seekg(pos);

	j = copy_token(token, pos);

	if (j == CParser::TT_EMPTY)
		return PARSER_OK;

	if (!(m_line_iss >> d))
	{
		return PARSER_ERROR;
	}
	totals[token] = d;
	return PARSER_OK;
}

int
CParser::getOptionFromLastLine(const std::vector < std::string > &opt_list,
							   std::string::iterator & next_char, bool flag_error)
{
	//
	// Read a line and check for options
	//
	int j;
	int /*  opt_l,  */ opt;
	std::string::iterator opt_ptr;
	std::string option;

	//
	// Read line
	//
	PHRQ_io::LINE_TYPE lt = m_line_type;
	if (lt == PHRQ_io::LT_EOF)
	{
		j = OPT_EOF;
	}
	else if (lt == PHRQ_io::LT_KEYWORD)
	{
		j = OPT_KEYWORD;
	}
	else if (lt == PHRQ_io::LT_OPTION)
	{
		opt_ptr = m_line.begin();
		std::string::iterator end = m_line.end();
		copy_token(option, opt_ptr, end);
		if (find_option(option, &opt, opt_list, false) == CParser::FT_OK)
		{
			j = opt;
			m_line_save.replace(m_line_save.find(option), option.size(),
								opt_list[opt]);
			m_line.replace(m_line.find(option), option.size(), opt_list[opt]);
			opt_ptr = m_line.begin();
			std::string::iterator end = m_line.end();
			copy_token(option, opt_ptr, end);
			next_char = opt_ptr;
			if (true)			// pr.echo_input == TRUE
			{
				if (true)		// database_file == NULL
				{
					std::ostringstream msg;
					msg << "\t" << m_line_save << "\n";
					io->output_msg(msg.str().c_str());
				}
			}
		}
		else
		{
			if (flag_error)
			{
				if (true)			// (database_file == NULL)
				{
					std::ostringstream msg;
					msg << "\t" << m_line_save << "\n";
					io->output_msg(msg.str().c_str());
				}
				std::ostringstream err;
				err << "Unknown option." << "\n";
				err << m_line_save << "\n";
				error_msg(err.str().c_str());
			}
			j = OPT_ERROR;
			next_char = m_line.begin();
		}
	}
	else
	{
		opt_ptr = m_line.begin();
		std::string::iterator end = m_line.end();
		copy_token(option, opt_ptr, end);
		if (find_option(option, &opt, opt_list, true) == FT_OK)
		{
			j = opt;
			next_char = opt_ptr;
		}
		else
		{
			j = OPT_DEFAULT;
			next_char = m_line.begin();
		}
		if (true)				// pr.echo_input == TRUE
		{
			if (true)			// database_file == NULL
			{
#if !defined(R_SO)
				std::cout << "\t" << m_line_save << "\n";
#endif
			}
		}
	}
	return (j);
}

int
CParser::getOptionFromLastLine(const std::vector < std::string > &opt_list,
							   std::istream::pos_type & next_pos, bool flag_error)
{
	//
	// Read a line and check for options
	//
	int j;
	int opt;
	std::istream::pos_type pos_ptr;
	std::string option;

	//
	// Read line
	//
	PHRQ_io::LINE_TYPE lt = m_line_type;
	if (lt == PHRQ_io::LT_EOF)
	{
		j = OPT_EOF;
	}
	else if (lt == PHRQ_io::LT_KEYWORD)
	{
		j = OPT_KEYWORD;
	}
	else if (lt == PHRQ_io::LT_OPTION)
	{
		std::string::iterator opt_ptr = m_line.begin();
		std::string::iterator end = m_line.end();
		copy_token(option, opt_ptr, end);
		if (find_option(option.substr(1), &opt, opt_list, false) == FT_OK)
		{
			// replace -option with option
			j = opt;
			m_line_save.replace(m_line_save.find(option), option.size(),
								opt_list[opt]);
			m_line.replace(m_line.find(option), option.size(), opt_list[opt]);

			// reset iss
			m_line_iss.str(m_line);
			m_line_iss.seekg(0, std::ios_base::beg);
			m_line_iss.clear();

			pos_ptr = 0;
			copy_token(option, pos_ptr);
			next_pos = pos_ptr;
			//if (this->echo_file)			// pr.echo_input == TRUE
			if (false)
			{
				if (true)		// database_file == NULL
				{
					std::ostringstream msg;
					msg << "\t" << m_line_save << "\n";
					io->output_msg(msg.str().c_str());
				}
			}
		}
		else
		{
			if (flag_error) // don`t throw error
			{
				if (true)			// (database_file == NULL)
				{
					std::ostringstream msg;
					msg << "\t" << m_line_save << "\n";
					io->output_msg(msg.str().c_str());
				}
				error_msg("Unknown option.", PHRQ_io::OT_CONTINUE);
				error_msg(m_line_save.c_str(), PHRQ_io::OT_CONTINUE);
				incr_input_error();
			}
			j = OPT_ERROR;
			next_pos = pos_ptr;
		}
	}
	else
	{
		pos_ptr = 0;
		copy_token(option, pos_ptr);
		if (find_option(option, &opt, opt_list, true) == FT_OK)
		{
			j = opt;
			next_pos = pos_ptr;
		}
		else
		{
			j = OPT_DEFAULT;
			next_pos = 0;
		}
		if (true)				// pr.echo_input == TRUE
		{
			if (true)			// database_file == NULL
			{
				std::ostringstream msg;
				msg << "\t" << m_line_save << "\n";
				io->output_msg(msg.str().c_str());
			}
		}
	}
	return (j);
}
int CParser::
incr_input_error()
{
	return ++m_input_error;
}

CParser::TOKEN_TYPE CParser::copy_title(std::string & token,
										std::string::iterator & begin,
										std::string::iterator & end)
{
	if (begin != end)
	{
		std::string::iterator b = begin;
		std::string::iterator e = end;
		for (; b < end && (::isspace(*b) || (*b == ',')); ++b);
		begin = b;
		if (*begin == '"')
		{
			begin = ++b;
			for (; begin != end && !(*begin == '"'); ++begin);
			e = begin;
			if (begin != end && *begin == '"')
			{
				e = begin++;
			}
		}
		else if (*begin == '\'')
		{
			begin = ++b;
			for (; begin != end && !(*begin == '\''); ++begin);
			e = begin;
			if (begin != end && *begin == '\'')
			{
				e = begin++;
			}
		}
		else
		{
			for (; begin < end && !(*begin == ',') && !(::isspace(*begin)); ++begin);
			e = begin;
		}
		token.assign(b, e);
	}
	else
	{
		token.resize(0);
	}
	token = trim(token);
	return token_type(token);
}
bool CParser::get_true_false(std::istream::pos_type & pos, bool def)
{
	std::string token;
	this->copy_token(token, pos);
	std::string::iterator b = token.begin();
	for (; b != token.end() && (::isspace(*b)); ++b);
	if (b != token.end())
	{
		if (*b == 'f' || *b == 'F')
		{
			return false;
		}
		else if (*b == 't' || *b == 'T')
		{
			return true;
		}
	}
	return def;
}
CParser::TOKEN_TYPE CParser::get_rest_of_line(std::string &token)
{
	token.clear();
	int j;
	while ((j = m_line_iss.get()) != std::char_traits < char >::eof())
	{
		char c = (char) j;
		token += c;
	}
	token = trim(token);
	return token_type(token);
}
CParser::TOKEN_TYPE CParser::parse_delimited(std::string & source, std::string & result,
										const std::string& t = " \t")
{

	size_t pos = source.find_first_of(t);
	std::string temp;
	if (pos != std::string::npos)
	{
		result = source.substr(0, pos);
		temp = source.substr(pos+1);
		source = temp;
	}
	else
	{
		result = source;
		source.clear();
	}
	std::string str = result;

	return token_type(trim_left(str));
}