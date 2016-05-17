// Exchange.cxx: implementation of the cxxExchange class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <iostream>				// std::cout std::cerr
#include <cassert>				// assert
#include <algorithm>			// std::sort

#include "Utils.h" // define first

#include "Phreeqc.h"
#include "cxxMix.h"
#include "Exchange.h"
#include "phqalloc.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxExchange::cxxExchange(PHRQ_io *io)
	//
	// default constructor for cxxExchange
	//
:	cxxNumKeyword(io)
{
	new_def = false;
	solution_equilibria = false;
	n_solution = -999;
	pitzer_exchange_gammas = true;
}
cxxExchange::cxxExchange(const std::map < int, cxxExchange > &entities,
						 cxxMix & mix, int l_n_user, PHRQ_io *io):
cxxNumKeyword(io)
{
	this->n_user = this->n_user_end = l_n_user;
	this->pitzer_exchange_gammas = true;
	this->new_def = false;
	this->n_solution = -999;
//
//   Mix exchangers
//
	const std::map < int, LDBLE >&mixcomps = mix.Get_mixComps();
	std::map < int, LDBLE >::const_iterator it;
	for (it = mixcomps.begin(); it != mixcomps.end(); it++)
	{
		if (entities.find(it->first) != entities.end())
		{
			const cxxExchange *entity_ptr =
				&(entities.find(it->first)->second);
			this->add(*entity_ptr, it->second);
			this->pitzer_exchange_gammas = entity_ptr->pitzer_exchange_gammas;
		}
	}
}


cxxExchange::~cxxExchange()
{
}
bool
cxxExchange::Get_related_phases() const
{
	for (size_t i = 0; i < this->exchange_comps.size(); i++)
	{
		if (this->exchange_comps[i].Get_phase_name().size() > 0)
		{
			return (true);
		}
	}
	return (false);
}

bool
cxxExchange::Get_related_rate() const
{
	for (size_t i = 0; i < this->exchange_comps.size(); i++)
	{
		if (this->exchange_comps[i].Get_rate_name().size() > 0)
		{
			return (true);
		}
	}
	return (false);
}

void
cxxExchange::dump_xml(std::ostream & s_oss, unsigned int indent) const
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

	// Exchange element and attributes
	s_oss << indent0;
	s_oss << "<exchange " << "\n";

	s_oss << indent1;
	s_oss << "pitzer_exchange_gammas=\"" << this->
		pitzer_exchange_gammas << "\"" << "\n";

	// components
	s_oss << indent1;
	s_oss << "<component " << "\n";
	for (size_t j = 0; j < this->exchange_comps.size(); j++)
	{
		this->exchange_comps[j].dump_xml(s_oss, indent + 2);
	}

	return;
}
#ifdef SKIP
void
cxxExchange::dump_xml(std::ostream & s_oss, unsigned int indent) const
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

	// Exchange element and attributes
	s_oss << indent0;
	s_oss << "<exchange " << "\n";

	s_oss << indent1;
	s_oss << "pitzer_exchange_gammas=\"" << this->
		pitzer_exchange_gammas << "\"" << "\n";

	// components
	s_oss << indent1;
	s_oss << "<component " << "\n";
	for (std::map < std::string, cxxExchComp >::const_iterator it = exchComps.begin();
		 it != exchComps.end(); ++it)
	{
		(*it).second.dump_xml(s_oss, indent + 2);
	}

	return;
}
#endif
void
cxxExchange::dump_raw(std::ostream & s_oss, unsigned int indent, int *n_out) const
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

	// Exchange element and attributes
	s_oss << indent0;
	int n_user_local = (n_out != NULL) ? *n_out : this->n_user;
	s_oss << "EXCHANGE_RAW                 " << n_user_local << " " << this->description << "\n";

	s_oss << indent1 << "# EXCHANGE_MODIFY candidate identifiers #\n";
	s_oss << indent1;
	s_oss << "-exchange_gammas           " << (this->pitzer_exchange_gammas ? 1 : 0) << "\n";
	// exchComps 
	for (size_t i = 0; i < this->exchange_comps.size(); i++)
	{
		s_oss << indent1;
		s_oss << "-component                 " << this->exchange_comps[i].Get_formula() << "\n";
		this->exchange_comps[i].dump_raw(s_oss, indent + 2);
	}

	s_oss << indent1 << "# EXCHANGE_MODIFY candidates with new_def=true #\n";
	s_oss << indent1;
	s_oss << "-new_def                   " << 0 << "\n";
	s_oss << indent1;
	s_oss << "-solution_equilibria       " << 0 << "\n";
	s_oss << indent1;
	s_oss << "-n_solution                " << this->n_solution << "\n";

	s_oss << indent1 << "# Exchange workspace variables #\n";
	s_oss << indent1;
	s_oss << "-totals" << "\n";
	this->totals.dump_raw(s_oss, indent + 1);

	return;
}
void
cxxExchange::read_raw(CParser & parser, bool check)
{
	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	bool useLastLine(false);

	// Read exchange number and description
	this->read_number_description(parser);
	this->Set_new_def(false);
	
	bool pitzer_exchange_gammas_defined(false);

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
			parser.error_msg("Unknown input in EXCH_COMP_RAW keyword.",
							 PHRQ_io::OT_CONTINUE);
			parser.error_msg(parser.line().c_str(), PHRQ_io::OT_CONTINUE);
			break;

		case 0:				// pitzer_exchange_gammas
		case 2:				// exchange_gammas
			if (!(parser.get_iss() >> this->pitzer_exchange_gammas))
			{
				this->pitzer_exchange_gammas = false;
				parser.incr_input_error();
				parser.
					error_msg
					("Expected boolean value for pitzer_exchange_gammas.",
					 PHRQ_io::OT_CONTINUE);
			}
			pitzer_exchange_gammas_defined = true;
			break;
		case 1:				// component
			{
				std::string str;
				if (!(parser.get_iss() >> str))
				{
					parser.incr_input_error();
					parser.error_msg("Expected string value for component name.",
									 PHRQ_io::OT_CONTINUE);
				}
				else
				{
					cxxExchComp temp_comp(this->io);
					temp_comp.Set_formula(str.c_str());
					cxxExchComp *comp_ptr = this->Find_comp(str);
					if (comp_ptr)
					{
						temp_comp = *comp_ptr;
					}
					temp_comp.read_raw(parser, check);
					if (comp_ptr)
					{
						for (size_t j = 0; j < this->exchange_comps.size(); j++)
						{
							if (Utilities::strcmp_nocase(this->exchange_comps[j].Get_formula().c_str(), str.c_str()) == 0)
							{
								this->exchange_comps[j] = temp_comp;
							}
						}
					}
					else
					{
						this->exchange_comps.push_back(temp_comp);
					}
					useLastLine = true;
				}
			}
			break;
		case 3:				// new_def
			if (!(parser.get_iss() >> this->new_def))
			{
				this->new_def = false;
				parser.incr_input_error();
				parser.
					error_msg
					("Expected boolean value for new_def.",
					 PHRQ_io::OT_CONTINUE);
			}
			break;
		case 4:				// solution_equilibria
			if (!(parser.get_iss() >> this->solution_equilibria))
			{
				this->solution_equilibria = false;
				parser.incr_input_error();
				parser.
					error_msg
					("Expected boolean value for solution_equilibria.",
					 PHRQ_io::OT_CONTINUE);
			}
			break;
		case 5:				// n_solution
			if (!(parser.get_iss() >> this->n_solution))
			{
				this->n_solution = -999;
				parser.incr_input_error();
				parser.
					error_msg
					("Expected integer value for n_solution.",
					 PHRQ_io::OT_CONTINUE);
			}
			break;
		case 6:				// totals
			if (this->totals.read_raw(parser, next_char) !=
				CParser::PARSER_OK)
			{
				parser.incr_input_error();
				parser.
					error_msg
					("Expected element name and molality for Exchange totals.",
					PHRQ_io::OT_CONTINUE);
			}
			break;
		}
		if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD)
			break;
	}
	// members that must be defined
	if (check)
	{
		if (pitzer_exchange_gammas_defined == false)
		{
			parser.incr_input_error();
			parser.
				error_msg
				("Pitzer_exchange_gammsa not defined for EXCHANGE_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
	}
	this->Sort_comps();
}
void
cxxExchange::add(const cxxExchange & addee, LDBLE extensive)
		//
		// Add existing exchange to "this" exchange
		//
{
	// exchComps
	if (extensive == 0.0)
		return;
	for (size_t i = 0; i < addee.exchange_comps.size(); i++)
	{
		size_t j;
		for (j = 0; j < this->Get_exchange_comps().size(); j++)
		if (addee.exchange_comps[i].Get_formula() == this->exchange_comps[j].Get_formula())
		{
			this->exchange_comps[j].add(addee.exchange_comps[i], extensive);
			break;
		}
		if (j == this->exchange_comps.size())
		{
			cxxExchComp exc = addee.exchange_comps[i];
			exc.multiply(extensive);
			this->exchange_comps.push_back(exc);
		}
	}
	this->pitzer_exchange_gammas = addee.pitzer_exchange_gammas;
}
void
cxxExchange::totalize()
{
	this->totals.clear();
	// component structures
	for (size_t i = 0; i < this->exchange_comps.size(); i++)
	{
		this->totals.add_extensive(this->exchange_comps[i].Get_totals(), 1.0);
		this->totals.add("Charge", this->exchange_comps[i].Get_charge_balance());
	}
	return;
}
bool 
cxxExchange::Get_pitzer_exchange_gammas() const
{
	return this->pitzer_exchange_gammas;
}
void 
cxxExchange::Set_pitzer_exchange_gammas(bool b)
{
	this->pitzer_exchange_gammas = b;
}
const cxxNameDouble & 
cxxExchange::Get_totals() const
{
	return totals;
}
cxxExchComp *cxxExchange::Find_comp(std::string s)
{
	for (size_t i = 0; i < this->exchange_comps.size(); i++)
	{
		cxxNameDouble nd(this->exchange_comps[i].Get_totals());
		cxxNameDouble::iterator nd_it;
		for (nd_it = nd.begin(); nd_it != nd.end(); nd_it++)
		{
			if(nd_it->first == s)
			{
				return (&(this->exchange_comps[i]));
			}
		}
	}
	return NULL;
}
void cxxExchange::
Sort_comps(void)
{
	// sort comps
	{
		std::map<std::string, cxxExchComp> comp_map;
		for (size_t i = 0; i < this->exchange_comps.size(); i++)
		{
			comp_map[this->exchange_comps[i].Get_formula()] = this->exchange_comps[i];
		}
		this->exchange_comps.clear();
		std::map<std::string, cxxExchComp>::iterator it;
		for (it = comp_map.begin(); it != comp_map.end(); it++)
		{
			this->exchange_comps.push_back(it->second);
		}
	}
}
/* ---------------------------------------------------------------------- */
void
cxxExchange::Serialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles)
/* ---------------------------------------------------------------------- */
{
	ints.push_back(this->n_user);
	ints.push_back((int) this->exchange_comps.size());
	for (size_t i = 0; i < this->exchange_comps.size(); i++)
	{
		exchange_comps[i].Serialize(dictionary, ints, doubles);
	}
	ints.push_back(this->pitzer_exchange_gammas ? 1 : 0);
	ints.push_back(this->new_def ? 1 : 0);
	ints.push_back(this->solution_equilibria ? 1 : 0);
	ints.push_back(this->n_solution);
	this->totals.Serialize(dictionary, ints, doubles);

}

/* ---------------------------------------------------------------------- */
void
cxxExchange::Deserialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles, int &ii, int &dd)
/* ---------------------------------------------------------------------- */
{
	this->n_user = ints[ii++];
	this->n_user_end = this->n_user;
	this->description = " ";

	int count = ints[ii++];
	this->exchange_comps.clear();
	for (int n = 0; n < count; n++)
	{
		cxxExchComp ec;
		ec.Deserialize(dictionary, ints, doubles, ii, dd);
		this->exchange_comps.push_back(ec);
	}
	this->pitzer_exchange_gammas = (ints[ii++] != 0);
	this->new_def = (ints[ii++] != 0);
	this->solution_equilibria = (ints[ii++] != 0);
	this->n_solution = ints[ii++];
	this->totals.Deserialize(dictionary, ints, doubles, ii, dd);

}


const std::vector< std::string >::value_type temp_vopts[] = {
	std::vector< std::string >::value_type("pitzer_exchange_gammas"),   // 0
	std::vector< std::string >::value_type("component"),		        // 1
	std::vector< std::string >::value_type("exchange_gammas"),	        // 2
	std::vector< std::string >::value_type("new_def"),		            // 3
	std::vector< std::string >::value_type("solution_equilibria"),	    // 4
	std::vector< std::string >::value_type("n_solution"),		        // 5
	std::vector< std::string >::value_type("totals")		            // 6
};									   
const std::vector< std::string > cxxExchange::vopts(temp_vopts, temp_vopts + sizeof temp_vopts / sizeof temp_vopts[0]);
