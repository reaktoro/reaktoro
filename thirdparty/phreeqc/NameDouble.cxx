// NameDouble.cxx: implementation of the cxxNameDouble class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <cassert>				// assert
#include <algorithm>			// std::sort
#include <map>					// std::sort
#include <iostream>				// std::cout std::cerr

#include "Utils.h"				// define first
#include "Phreeqc.h"
#include "NameDouble.h"
#include "Dictionary.h"
//#include "Dictionary.h"
#include "phqalloc.h"
#include "ISolutionComp.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxNameDouble::cxxNameDouble()
	//
	// default constructor for cxxNameDouble 
	//
{
	this->type = ND_ELT_MOLES;
}

cxxNameDouble::cxxNameDouble(struct elt_list *elt_list_ptr)
		//
		// constructor for cxxNameDouble from list of elt_list
		//
{
	int i;
	if (elt_list_ptr != NULL)
	{
		for (i = 0; elt_list_ptr[i].elt != NULL; i++)
		{
			(*this)[elt_list_ptr[i].elt->name] = elt_list_ptr[i].coef;
		}
	}
	this->type = ND_ELT_MOLES;
}

cxxNameDouble::cxxNameDouble(struct elt_list *elt_list_ptr, int count)
		//
		// constructor for cxxNameDouble from list of elt_list with known count
		//
{
	int i;
	if (elt_list_ptr != NULL)
	{
		for (i = 0; i < count; i++)
		{
			(*this)[elt_list_ptr[i].elt->name] = elt_list_ptr[i].coef;
		}
	}
	this->type = ND_ELT_MOLES;
}

cxxNameDouble::cxxNameDouble(const cxxNameDouble & old, LDBLE factor)
		//
		// constructor for cxxNameDouble from list of elt_list
		//
{
	for (cxxNameDouble::const_iterator it = old.begin(); it != old.end();
		 it++)
	{
		if (old.type == ND_ELT_MOLES)
		{
			if (it->second * factor > 0)
			{
				(*this)[(it->first)] = it->second * factor;
			}
		}
		else
		{
			(*this)[(it->first)] = it->second * factor;
		}
	}
	this->type = old.type;
}
cxxNameDouble::cxxNameDouble(std::map < std::string, cxxISolutionComp > &comps)
		//
		// constructor for cxxNameDouble from map of cxxISolutionComp
		//
{
	std::map < std::string, cxxISolutionComp >::iterator it;
	for (it = comps.begin(); it != comps.end(); it++)
	{
		(*this)[it->first] = it->second.Get_moles();
	}
	this->type = ND_ELT_MOLES;
}
#ifdef SKIP
cxxNameDouble::cxxNameDouble(struct master_activity *ma, int count,
							 cxxNameDouble::ND_TYPE l_type)
		//
		// constructor for cxxNameDouble from list of elt_list
		//
{
	int i;
	for (i = 0; i < count; i++)
	{
		if (ma[i].description == NULL)
			continue;
		(*this)[ma[i].description] = ma[i].la;
	}
	this->type = l_type;
}
#endif
cxxNameDouble::cxxNameDouble(struct name_coef *nc, int count)
		//
		// constructor for cxxNameDouble from list of elt_list
		//
{
	int i;
	for (i = 0; i < count; i++)
	{
		if (nc[i].name == NULL)
			continue;

		if ((*this).find(nc[i].name) == (*this).end())
		{
			(*this)[nc[i].name] = nc[i].coef;
		}
		else
		{
			(*this)[nc[i].name] = (*this).find(nc[i].name)->second + nc[i].coef;
		}

	}
	this->type = ND_NAME_COEF;
}

cxxNameDouble::~cxxNameDouble()
{
}

void
cxxNameDouble::dump_xml(std::ostream & s_oss, unsigned int indent) const
{
	unsigned int i;
	s_oss.precision(DBL_DIG - 1);
	std::string indent0(""), indent1("");
	for (i = 0; i < indent; ++i)
		indent0.append(Utilities::INDENT);
	for (i = 0; i < indent + 1; ++i)
		indent1.append(Utilities::INDENT);
	std::string xmlElement, xmlAtt1, xmlAtt2;

	switch ((*this).type)
	{
	case cxxNameDouble::ND_SPECIES_LA:
		xmlElement = "<soln_m_a ";
		xmlAtt1 = " m_a_desc=\"";
		xmlAtt1 = " m_a_la=\"";
		break;
	case cxxNameDouble::ND_SPECIES_GAMMA:
		xmlElement = "<soln_s_g ";
		xmlAtt1 = " m_a_desc=\"";
		xmlAtt1 = " m_a_la=\"";
		break;
	case cxxNameDouble::ND_ELT_MOLES:
		xmlElement = "<soln_total ";
		xmlAtt1 = " conc_desc=\"";
		xmlAtt1 = " conc_moles=\"";
		break;
	case cxxNameDouble::ND_NAME_COEF:
		xmlElement = "<NameCoef ";
		xmlAtt1 = " name=\"";
		xmlAtt1 = " coef=\"";
		break;
	}

	for (const_iterator it = (*this).begin(); it != (*this).end(); it++)
	{
		s_oss << indent0;
		s_oss << xmlElement << xmlAtt1 << it->first << xmlAtt2 << it->
			second << "/>" << "\n";
	}
}

void
cxxNameDouble::dump_raw(std::ostream & s_oss, unsigned int indent) const
{
	unsigned int i;
	s_oss.precision(DBL_DIG - 1);
	std::string indent0("");
	for (i = 0; i < indent; ++i)
		indent0.append(Utilities::INDENT);

	for (const_iterator it = (*this).begin(); it != (*this).end(); it++)
	{
		s_oss << indent0;
		if (it->first.size() < 29 - indent0.size())
		{
			s_oss << Utilities::pad_right(it->first, 29 - indent0.size())  << it->second << "\n";
		}
		else
		{
			s_oss << Utilities::pad_right(it->first, it->first.size()  + indent0.size())  << " " << it->second << "\n";
		}
	}
}

CParser::STATUS_TYPE cxxNameDouble::read_raw(CParser & parser,
											 std::istream::pos_type & pos)
{
	std::string token;
	LDBLE
		d;

	CParser::TOKEN_TYPE j;
	j = parser.copy_token(token, pos);

	if (j == CParser::TT_EMPTY)
		return CParser::PARSER_OK;

	if (!(parser.get_iss() >> d))
	{
		return CParser::PARSER_ERROR;
	}
	(*this)[token.c_str()] = d;
	return CParser::PARSER_OK;
}

void
cxxNameDouble::add_extensive(const cxxNameDouble & addee, LDBLE factor)
//
// Sums two name doubles, this + factor*nd2
//
{
	if (factor == 0)
		return;
	for (cxxNameDouble::const_iterator it = addee.begin(); it != addee.end();
		 it++)
	{
		cxxNameDouble::iterator current = (*this).find(it->first);
		if (current != (*this).end())
		{
			(*this)[it->first] = current->second + it->second * factor;
		}
		else
		{
			(*this)[it->first] = it->second * factor;
		}
	}
}
void
cxxNameDouble::add_intensive(const cxxNameDouble & addee, LDBLE f1,
							 LDBLE f2)
//
// Sums two name doubles, this*f1 + f2*nd2
//
{
	assert(f1 >= 0 && f2 >= 0);
	for (cxxNameDouble::const_iterator it = addee.begin(); it != addee.end();
		 it++)
	{
		cxxNameDouble::iterator current = (*this).find(it->first);
		if (current != (*this).end())
		{
			(*this)[it->first] = f1 * current->second + f2 * it->second;
		}
		else
		{
			(*this)[it->first] = f2 * it->second;
		}
	}
}
void
cxxNameDouble::add_log_activities(const cxxNameDouble & addee, LDBLE f1,
								  LDBLE f2)
//
// Sums two name doubles, this*f1 + f2*nd2, assuming log values
//
{
	assert(f1 >= 0 && f2 >= 0);
	for (cxxNameDouble::const_iterator it = addee.begin(); it != addee.end();
		 it++)
	{
		cxxNameDouble::iterator current = (*this).find(it->first);
		if (current != (*this).end())
		{
			LDBLE a1 = pow((LDBLE) 10., current->second);
			LDBLE a2 = pow((LDBLE) 10., it->second);
			(*this)[it->first] = log10(f1 * a1 + f2 * a2);
		}
		else
		{
			(*this)[it->first] = it->second + log10(f2);
		}
	}
}
cxxNameDouble 
cxxNameDouble::Simplify_redox(void) const
{
	cxxNameDouble const &nd = *this;
	std::basic_string < char >::size_type indexCh;
	cxxNameDouble new_totals;
	new_totals.type = cxxNameDouble::ND_ELT_MOLES;
	{
		std::string current_ename;
		std::string const *ename_ptr;
		cxxNameDouble::const_iterator it;

		// make list of elements in new_totals
		for (it = nd.begin(); it != nd.end(); ++it)
		{
			current_ename = it->first;
			if (it->first.size() < 4)
			{
				ename_ptr = &(it->first);
			}
			else
			{
				indexCh = it->first.find("(");
				if (indexCh != std::string::npos)
				{
					current_ename = it->first.substr(0, indexCh);
					ename_ptr = &(current_ename);
				}
				else
				{
					ename_ptr = &(it->first);
				}
			}
			if (current_ename == "H" || current_ename == "O" || current_ename == "Charge")
				continue;
			new_totals[*ename_ptr] = 0;
		}
	}

	// sum totals for elements
	{
		cxxNameDouble::const_iterator old_it = nd.begin();
		cxxNameDouble::iterator new_it = new_totals.begin();
		std::string old_ename;
		std::string const *old_ename_ptr;
		while (old_it != nd.end() && new_it != new_totals.end())
		{
			if (old_it->first.size() < 4)
			{
				old_ename_ptr = &old_it->first;
			}
			else
			{
				indexCh = old_it->first.find("(");
				if (indexCh != std::string::npos)
				{
					old_ename = old_it->first.substr(0, indexCh);
					old_ename_ptr = &old_ename;
				}
				else
				{
					old_ename_ptr = &old_it->first;
				}
			}
			int j = strcmp(new_it->first.c_str(), old_ename_ptr->c_str());
			if (j < 0)
			{
				new_it++;
			}
			else if (j == 0)
			{
				new_it->second += old_it->second;
				old_it++;
			}
			else 
			{
				old_it++;
			}
		}
	}
	return new_totals;
}
#ifdef SKIP
cxxNameDouble 
cxxNameDouble::Simplify_redox(void)
{
	// remove individual redox states from totals
	cxxNameDouble &nd = *this;
	std::set<std::string> list_of_elements;
	cxxNameDouble::const_iterator it;
	for (it = nd.begin(); it != nd.end(); ++it)
	{
		std::string current_ename(it->first);
		std::basic_string < char >::size_type indexCh;
		indexCh = current_ename.find("(");
		if (indexCh != std::string::npos)
		{
			current_ename = current_ename.substr(0, indexCh);
		}
		if (current_ename == "H" || current_ename == "O" || current_ename == "Charge")
			continue;
		list_of_elements.insert(current_ename);
	}

	cxxNameDouble new_totals;
	new_totals.type = cxxNameDouble::ND_ELT_MOLES;
	std::set<std::string>::iterator nt_it = list_of_elements.begin();
	for( ; nt_it != list_of_elements.end(); nt_it++)
	{
		new_totals[(*nt_it).c_str()] = nd.Get_total_element((*nt_it).c_str());
	}
	return new_totals;
}
#endif
void 
cxxNameDouble::Multiply_activities_redox(std::string str, LDBLE f)
{
	// update original master_activities using just computed factors
	cxxNameDouble::iterator it;
	LDBLE lg_f = log10(f);
	std::string redox_name = str;
	redox_name.append("(");

	for (it = this->begin(); it != this->end(); it++)
	{
		if (str[0] > it->first[0]) continue;
		if (it->first == str)
		{
			// Found exact match
			it->second += lg_f;
		}
		else 
		{
			// no exact match, current is element name, need to find all valences
			if (strstr(it->first.c_str(), redox_name.c_str()) == it->first.c_str())
			{
				it->second += lg_f;
			}
		}
		if (str[0] < it->first[0]) break;
	}
}
#ifdef SKIP
void 
cxxNameDouble::Multiply_activities_redox(std::string str, LDBLE f)
{
	// update original master_activities using just computed factors
	cxxNameDouble::iterator it;
	LDBLE lg_f = log10(f);
	std::string redox_name = str;
	redox_name.append("(");

	for (it = this->begin(); it != this->end(); it++)
	{
		if (it->first == str)
		{
			// Found exact match
			it->second += lg_f;
		}
		else 
		{
			// no exact match, current is element name, need to find all valences
			if (strstr(it->first.c_str(), redox_name.c_str()) == it->first.c_str())
			{
				it->second += lg_f;
			}
		}
	}
}
#endif
LDBLE
cxxNameDouble::Get_total_element(const char *string) const
{
	cxxNameDouble::const_iterator it;
	LDBLE d = 0.0;
	for (it = this->begin(); it != this->end(); ++it)
	{
		// C++ way to do it
		std::string ename(string);
		std::string current_ename(it->first);
		std::basic_string < char >::size_type indexCh;
		indexCh = current_ename.find("(");
		if (indexCh != std::string::npos)
		{
			current_ename = current_ename.substr(0, indexCh);
		}
		if (current_ename == ename)
		{
			d += it->second;
		}
	}
	return (d);
}
void
cxxNameDouble::add(const char *token, LDBLE total)
//
// add to total for a specified element
//
{
	char key[MAX_LENGTH];
	strcpy(key, token);

	cxxNameDouble::iterator current = (*this).find(key);
	if (current != (*this).end())
	{
		(*this)[key] = current->second + total;
	}
	else
	{
		(*this)[key] = total;
	}
}
void
cxxNameDouble::multiply(LDBLE extensive)
{
//
// Multiplies by extensive
//

	for (cxxNameDouble::iterator it = this->begin(); it != this->end(); it++)
	{
		it->second *= extensive;
	}
}
void
cxxNameDouble::merge_redox(const cxxNameDouble & source)
//
// Merges source into this
// Accounts for possible conflicts between redox state and
// totals
//
{
	for (cxxNameDouble::const_iterator sit = source.begin(); sit != source.end(); sit++)
	{
		
		std::string redox_name = sit->first;
		std::string elt_name;
		size_t pos = redox_name.find("(");

		bool redox;
		if (pos != std::string::npos) 
		{
			redox = true;
			elt_name = redox_name.substr(0, pos - 1);
		}
		else
		{
			redox = false;
			elt_name = redox_name;
		}
		if (redox)
		{
			// Remove elt_name, if present
			if ((*this).find(elt_name) != (*this).end())
			{
				(*this).erase((*this).find(elt_name));
			}
			// Put in redox name
			(*this)[redox_name] = sit->second;
		}
		else
		{
			std::string substring;
			substring.append(elt_name);
			substring.append("(");

			// Remove all redox
			bool deleted = true;
			while (deleted)
			{
				deleted = false;
				cxxNameDouble::iterator current = (*this).begin();
				for ( ; current != (*this).end(); current++)
				{
					if (current->first.find(substring) == 0)
					{
						(*this).erase(current);
						deleted = true;
						break;
					}
				}
			}
			// Put in elt name
			(*this)[elt_name] = sit->second;
		}
	}
}
struct DblCmp {     
	bool operator()(const std::pair<std::string, LDBLE> &lhs, const std::pair<std::string, LDBLE> &rhs) 
	{         
		return lhs.second > rhs.second;     
	} 
}; 
std::vector< std::pair<std::string, LDBLE> > 
cxxNameDouble::sort_second(void)
{
	std::vector< std::pair<std::string, LDBLE> > myvec(this->begin(), this->end()); 
	std::sort(myvec.begin(), myvec.end(), DblCmp());

	return myvec;
}
void
cxxNameDouble::Serialize(Dictionary &dictionary, std::vector < int >&ints,
						std::vector < double >&doubles)
{
	ints.push_back((int) (*this).size());
	for (const_iterator it = (*this).begin(); it != (*this).end(); it++)
	{
		int n = dictionary.Find(it->first);
		ints.push_back(n);
		doubles.push_back(it->second);
	}
}
void
cxxNameDouble::Deserialize(Dictionary &dictionary, std::vector<int> &ints, std::vector<double> &doubles, int &ii, int &dd)
{
	this->clear();
	int count = ints[ii++];
	for (int j = 0; j < count; j++)
	{
		int n = ints[ii++];
		assert(n >= 0);
		std::string str = dictionary.GetWords()[n];
		if (str.size() != 0)
		{
			(*this)[str] = doubles[dd++];
		}
	}
}
