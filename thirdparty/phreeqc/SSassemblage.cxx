// SSassemblage.cxx: implementation of the cxxSSassemblage class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <float.h>
#include <cassert>				// assert
#include <algorithm>			// std::sort

#include "Utils.h"				// define first
#include "Phreeqc.h"
#include "SSassemblage.h"
#include "SS.h"
#include "cxxMix.h"
#include "phqalloc.h"
#include "Dictionary.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxSSassemblage::cxxSSassemblage(PHRQ_io * io)
	//
	// default constructor for cxxSSassemblage 
	//
:	cxxNumKeyword(io)
{
	new_def = false;
}

cxxSSassemblage::cxxSSassemblage(const std::map < int,
								 cxxSSassemblage > &entities, cxxMix & mix,
								 int l_n_user, PHRQ_io * io):
cxxNumKeyword(io)
{
	this->n_user = this->n_user_end = l_n_user;
//
//   Mix
//
	const std::map < int, LDBLE >&mixcomps = mix.Get_mixComps();
	std::map < int, LDBLE >::const_iterator it;
	for (it = mixcomps.begin(); it != mixcomps.end(); it++)
	{
		if (entities.find(it->first) != entities.end())
		{
			const cxxSSassemblage *entity_ptr =
				&(entities.find(it->first)->second);
			this->add(*entity_ptr, it->second);
		}
	}
	this->new_def = false;
}

cxxSSassemblage::~cxxSSassemblage()
{
}

#ifdef SKIP
void
cxxSSassemblage::dump_xml(std::ostream & s_oss, unsigned int indent) const const
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

	// SSassemblage element and attributes
	s_oss << indent0;
	s_oss << "<EQUILIBRIUM_PHASES " << "\n";

	// eltList
	this->eltList.dump_xml(s_oss, indent + 1);

	// SSs
	s_oss << indent1;
	s_oss << "<pure_phases " << "\n";
	for (std::list < cxxSS >::const_iterator it =
		 SSs.begin(); it != SSs.end(); ++it)
	{
		it->dump_xml(s_oss, indent + 2);
	}
}
#endif
void
cxxSSassemblage::dump_raw(std::ostream & s_oss, unsigned int indent, int *n_out) const
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

	// SSassemblage element and attributes
	s_oss << indent0;
	int n_user_local = (n_out != NULL) ? *n_out : this->n_user;
	s_oss << "SOLID_SOLUTIONS_RAW          " << n_user_local << " " << this->description << "\n";

	s_oss << indent1 << "# SOLID_SOLUTION_MODIFY candidate identifiers #\n";
	// SSs
	for (std::map < std::string, cxxSS >::const_iterator it =
		 SSs.begin(); it != SSs.end(); ++it)
	{
		s_oss << indent1;
		s_oss << "-solid_solution            " << it->first << "\n";
		(*it).second.dump_raw(s_oss, indent + 2);
	}

	s_oss << indent1 << "# SOLID_SOLUTION candidate identifiers with new_def=true #\n";
	s_oss << indent1;
	s_oss << "-new_def                   " << (this->new_def ? 1 : 0) << "\n";

	s_oss << indent1 << "# solid solution workspace variables #\n";
	s_oss << indent1;
	s_oss << "-SSassemblage_totals       " << "\n";
	this->totals.dump_raw(s_oss, indent + 2);
}
void
cxxSSassemblage::read_raw(CParser & parser, bool check)
{

	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	int opt_save;
	bool useLastLine(false);
	cxxSS *ss_ptr = NULL;

	// Read SSassemblage number and description
	this->read_number_description(parser);
	this->Set_new_def(false);

	opt_save = CParser::OPT_ERROR;

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
		useLastLine = false;
		switch (opt)
		{
		case CParser::OPT_EOF:
			break;
		case CParser::OPT_KEYWORD:
			break;
		case CParser::OPT_DEFAULT:
		case CParser::OPT_ERROR:
			opt = CParser::OPT_EOF;
			parser.
				error_msg("Unknown input in SOLID_SOLUTIONS_RAW or SOLID_SOLUTIONS_MODIFY keyword.",
						  PHRQ_io::OT_CONTINUE);
			parser.error_msg(parser.line().c_str(), PHRQ_io::OT_CONTINUE);
			break;

		case 0:				// solid_solution
			{
				std::string str;
				if (!(parser.get_iss() >> str))
				{
						parser.incr_input_error();
					parser.error_msg("Expected string value for solid solution name.",
									 PHRQ_io::OT_CONTINUE);
				}
				cxxSS temp_ss(this->Get_io());
				temp_ss.Set_name(str);
				ss_ptr = this->Find(str);
				if (ss_ptr)
				{ 
					temp_ss = *ss_ptr;
				}
				temp_ss.read_raw(parser, false);
				this->SSs[str] = temp_ss;
			}
			useLastLine = true;
			break;
		case 1:				// totals
			if (this->totals.read_raw(parser, next_char) !=
				CParser::PARSER_OK)
			{
				parser.incr_input_error();
				parser.
					error_msg
					("Expected element name and molality for SSassemblage totals.",
					 PHRQ_io::OT_CONTINUE);
			}
			opt_save = 1;
			break;
		case 2:				// new_def
			{
				int i;
				if (!(parser.get_iss() >> i))
				{
					parser.incr_input_error();
					parser.error_msg("Expected 0/1 for new_def.", PHRQ_io::OT_CONTINUE);
				}
				else
				{
					this->new_def = (i == 0) ? false : true;
				}
			}
			break;
		}
		if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD)
			break;
	}
}

void
cxxSSassemblage::totalize(Phreeqc * phreeqc_ptr)
{
	this->totals.clear();
	// component structures
	for (std::map < std::string, cxxSS >::iterator it =
		 SSs.begin(); it != SSs.end(); ++it)
	{
		(*it).second.totalize(phreeqc_ptr);
		this->totals.add_extensive((*it).second.Get_totals(), 1.0);
	}
	return;
}
void
cxxSSassemblage::add(const cxxSSassemblage & addee, LDBLE extensive)
		//
		// Add to existing ssassemblage to "this" ssassemblage
		//
{
	if (extensive == 0.0)
		return;

	for (std::map < std::string, cxxSS >::const_iterator itadd =
		 addee.SSs.begin(); itadd != addee.SSs.end();
		 ++itadd)
	{
		std::map < std::string, cxxSS >::iterator it =
			this->SSs.find((*itadd).first);
		if (it != this->SSs.end())
		{
			(*it).second.add((*itadd).second, extensive);
		}
		else
		{
			cxxSS entity = (*itadd).second;
			entity.multiply(extensive);
			std::string str(entity.Get_name());
			this->SSs[str] = entity;
		}
	}
}

std::vector<cxxSS *> cxxSSassemblage::
Vectorize(void)
{
	std::vector<cxxSS *> ss_v;
	std::map<std::string, cxxSS>::iterator it;
	for (it = this->SSs.begin(); it != this->SSs.end(); it++)
	{
		ss_v.push_back(&(it->second));
	}
	return ss_v;
}
cxxSS * cxxSSassemblage::
Find(const std::string &s)
{
	std::map<std::string, cxxSS>::iterator it;
	it = this->SSs.find(s);
	if (it != this->SSs.end())
		return &(it->second);
	return NULL;
}
void
cxxSSassemblage::Serialize(Dictionary & dictionary, std::vector < int >&ints, 
	std::vector < double >&doubles)
{
	/* int n_user; */
	ints.push_back(this->n_user);
	{
		ints.push_back((int) this->SSs.size());
		std::map < std::string, cxxSS >::iterator it;
		for (it = this->SSs.begin(); it != this->SSs.end();	it++)
		{
			(*it).second.Serialize(dictionary, ints, doubles);
		}
	}
	ints.push_back(this->new_def ? 1 : 0);
	this->totals.Serialize(dictionary, ints, doubles);
}

void
cxxSSassemblage::Deserialize(Dictionary & dictionary, std::vector < int >&ints, 
	std::vector < double >&doubles, int &ii, int &dd)
{

	this->n_user = ints[ii++];
	this->n_user_end = this->n_user;
	this->description = " ";
	{
		int count = ints[ii++];
		this->SSs.clear();
		for (int n = 0; n < count; n++)
		{
			cxxSS ssc;
			ssc.Deserialize(dictionary, ints, doubles, ii, dd);
			std::string str(ssc.Get_name());
			this->SSs[str] = ssc;
		}
	}
	this->new_def = (ints[ii++] != 0);
	this->totals.Deserialize(dictionary, ints, doubles, ii, dd);

}




const std::vector< std::string >::value_type temp_vopts[] = {
	std::vector< std::string >::value_type("solid_solution"),	    // 0
	std::vector< std::string >::value_type("ssassemblage_totals"),	// 1
	std::vector< std::string >::value_type("new_def") 	            // 2
};									   
const std::vector< std::string > cxxSSassemblage::vopts(temp_vopts, temp_vopts + sizeof temp_vopts / sizeof temp_vopts[0]);	
