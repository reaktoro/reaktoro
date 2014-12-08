// Mix.cxx: implementation of the cxxMix class.
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
#include "cxxMix.h"
#include "phqalloc.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxMix::cxxMix(PHRQ_io *io)
	//
	// default constructor for cxxMix 
	//
:	cxxNumKeyword(io)
{
}

cxxMix::~cxxMix()
{
}

#ifdef SKIP
void
cxxMix::dump_xml(std::ostream & s_oss, unsigned int indent) const const
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

	// Mix element and attributes
	s_oss << indent0;
	s_oss << "<mix " << "\n";

	s_oss << indent1;
	s_oss << "pitzer_mix_gammas=\"" << this->
		pitzer_mix_gammas << "\"" << "\n";

	// components
	s_oss << indent1;
	s_oss << "<component " << "\n";
	for (std::list < cxxMixComp >::const_iterator it = mixComps.begin();
		 it != mixComps.end(); ++it)
	{
		it->dump_xml(s_oss, indent + 2);
	}

	return;
}
#endif

void
cxxMix::dump_raw(std::ostream & s_oss, unsigned int indent, int *n_out) const
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

	// Mix element and attributes
	s_oss << indent0;
	int n_user_local = (n_out != NULL) ? *n_out : this->n_user;
	s_oss << "MIX_RAW                      " << n_user_local << " " << this->description << "\n";

	for (std::map < int, LDBLE >::const_iterator it = this->mixComps.begin();
		 it != this->mixComps.end(); it++)
	{
		s_oss << indent1 << it->first << "     " << it->second << "\n";
	}
}

void
cxxMix::read_raw(CParser & parser)
{

	int i;
	LDBLE d;

	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	int opt_save;
	bool useLastLine(false);

	// Read mix number and description
	this->read_number_description(parser);

	opt_save = CParser::OPT_DEFAULT;

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
		case CParser::OPT_ERROR:
			opt = CParser::OPT_EOF;
			parser.error_msg("Unknown input in MIX_COMP_RAW keyword.",
							 PHRQ_io::OT_CONTINUE);
			parser.error_msg(parser.line().c_str(), PHRQ_io::OT_CONTINUE);
			useLastLine = false;
			break;

		case CParser::OPT_DEFAULT:	// solution number, mix fraction
			if (parser.copy_token(token, next_char) != CParser::TT_EMPTY)
			{
				std::istringstream iss(token);
				if (!(iss >> i))
				{
					parser.incr_input_error();
					parser.
						error_msg
						("Expected integer value for solution number.",
						 PHRQ_io::OT_CONTINUE);
					break;
				}
				if (!(parser.get_iss() >> d))
				{
					parser.incr_input_error();
					parser.
						error_msg
						("Expected numeric value for solution fraction.",
						 PHRQ_io::OT_CONTINUE);
					break;
				}
				this->mixComps[i] = d;
			}
			opt_save = CParser::OPT_DEFAULT;
			break;
		}
		if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD)
			break;
	}
}
void cxxMix::Vectorize(std::vector<int> &n, std::vector<LDBLE> &f)
{
	n.clear();
	f.clear();
	for (std::map < int, LDBLE >::const_iterator it = this->mixComps.begin();
		 it != this->mixComps.end(); it++)
	{
		n.push_back(it->first);
		f.push_back(it->second);
	}
}
const std::vector< std::string > cxxMix::vopts;