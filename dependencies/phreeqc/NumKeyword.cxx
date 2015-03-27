#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
// NumKeyword.cxx: implementation of the cxxNumKeyword class.
//
//////////////////////////////////////////////////////////////////////

#include <stdio.h>           // ::sscanf
#include "NumKeyword.h"
#include "Parser.h"
#include "Utils.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxNumKeyword::cxxNumKeyword(PHRQ_io *io)
{
	this->n_user = 1;
	this->n_user_end = 1;
	this->Set_io(io);
}

cxxNumKeyword::~cxxNumKeyword()
{
}

void
cxxNumKeyword::dump_xml(std::ostream & os, unsigned int indent) const
{
	unsigned int i;

	for (i = 0; i < indent + 1; ++i)
		os << "  ";
	os << "<n_user>" << this->n_user << "</n_user>" << "\n";

	for (i = 0; i < indent + 1; ++i)
		os << "  ";
	os << "<n_user_end>" << this->n_user_end << "</n_user_end>" << "\n";

	for (i = 0; i < indent + 1; ++i)
		os << "  ";
	os << "<Description>" << this->description << "</Description>" << "\n";
}

void
cxxNumKeyword::read_number_description(CParser & parser)
{
	std::string keyword;
	std::istream::pos_type ptr;

	// skip keyword
	parser.copy_token(keyword, ptr);

	// skip whitespace
	while (::isspace(parser.get_iss().peek()))
		parser.get_iss().ignore();

	// read number
	if (::isdigit(parser.get_iss().peek()) || parser.get_iss().peek() == '-')
	{
		parser.get_iss() >> this->n_user;
		char ch = (char)parser.get_iss().peek();
		if (ch == '-')
		{
			parser.get_iss() >> ch;			// eat '-'
			parser.get_iss() >> this->n_user_end;
			if (this->n_user_end < this->n_user)
			{
				this->n_user_end = this->n_user;
			}
		}
		else
		{
			this->n_user_end = this->n_user;
		}
	}
	else
	{
		this->n_user = this->n_user_end = 1;
	}

	// skip whitespace
	while (::isspace(parser.get_iss().peek()))
		parser.get_iss().ignore();

	// copy description
	std::getline(parser.get_iss(), this->description);
}


void
cxxNumKeyword::read_number_description(std::istream & is)
{
	// KEYWORD [[1[-20]] [This is the description]]

	// eat keyword
	std::string token;
	is >> token;

	// skip whitespace
	while (::isspace(is.peek()))
		is.ignore();

	if (::isdigit(is.peek()))
	{
		is >> this->n_user;
		char ch = (char)is.peek();
		if (ch == '-')
		{
			is >> ch;			// eat '-'
			is >> this->n_user_end;
		}
		else
		{
			this->n_user_end = this->n_user;
		}
	}
	else
	{
		this->n_user = this->n_user_end = 1;
	}

	while (::isspace(is.peek()))
		is.ignore();

	std::getline(is, this->description);
}
void
cxxNumKeyword::read_number_description(const std::string & line_in)
{
	std::string keyword, token;
	//std::istream::pos_type ptr;

	std::string line = line_in;
	std::string::iterator b = line.begin();
	std::string::iterator e = line.end();
	// skip keyword
	CParser::copy_token(keyword, b, e);

	// read number
	if (CParser::copy_token(token, b, e) == CParser::TT_DIGIT)
	{
		if (token[0] == '-') 
		{
			token = token.substr(1);
			Utilities::replace("-", " ", token);
			token = "-" + token;
		}
		else
		{
			Utilities::replace("-", " ", token);
		}
		int j = sscanf(token.c_str(), "%d%d", &this->n_user, &this->n_user_end);
		if (j == 0)
		{
			this->n_user = this->n_user_end = 1;
		}
		else if (j == 1)
		{
			this->n_user_end = this->n_user;
		}
		if (this->n_user_end < this->n_user)
		{
			this->n_user_end = this->n_user;
		}
	}
	else
	{
		this->n_user = this->n_user_end = 1;
	}

	// skip whitespace
	std::string::iterator ic;
	this->description.clear();
	for (ic = b; ic != e; ic++)
	{
		this->description += *ic;
	}
	trim_left(this->description);
}
