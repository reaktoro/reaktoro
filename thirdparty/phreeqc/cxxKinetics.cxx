// cxxKinetics.cxx: implementation of the cxxKinetics class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <cassert>				// assert
#include <algorithm>			// std::sort

#include "Utils.h"				
#include "Phreeqc.h"
#include "cxxKinetics.h"
#include "cxxMix.h"
#include "phqalloc.h"
#include "PHRQ_io.h"
#include "Dictionary.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxKinetics::cxxKinetics(PHRQ_io *io)
	//
	// default constructor for cxxKinetics 
	//
:	cxxNumKeyword(io)
{
	step_divide = 1.0;
	rk = 3;
	bad_step_max = 500;
	use_cvode = false;
	cvode_steps = 100;
	cvode_order = 5;
	totals.type = cxxNameDouble::ND_ELT_MOLES;
	equalIncrements = false;
	count = 0;
}
cxxKinetics::cxxKinetics(const std::map < int, cxxKinetics > &entities,
						 cxxMix & mix, int l_n_user, PHRQ_io *io):
cxxNumKeyword(io)
{
	this->n_user = this->n_user_end = l_n_user;
	step_divide = 1.0;
	rk = 3;
	bad_step_max = 500;
	use_cvode = false;
	cvode_steps = 100;
	cvode_order = 5;
	totals.type = cxxNameDouble::ND_ELT_MOLES;
	equalIncrements = false;
	count = 0;
//
//   Mix
//
	const std::map < int, LDBLE >&mixcomps = mix.Get_mixComps();
	std::map < int, LDBLE >::const_iterator it;
	for (it = mixcomps.begin(); it != mixcomps.end(); it++)
	{
		if (entities.find(it->first) != entities.end())
		{
			const cxxKinetics *entity_ptr =
				&(entities.find(it->first)->second);
			this->add(*entity_ptr, it->second);
		}
	}
}

cxxKinetics::~cxxKinetics()
{
}

#ifdef SKIP
void
cxxKinetics::dump_xml(std::ostream & s_oss, unsigned int indent) const const
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

	// Kinetics element and attributes
	s_oss << indent0;
	s_oss << "<kinetics " << "\n";

	s_oss << indent1;
	s_oss << "pitzer_kinetics_gammas=\"" << this->
		pitzer_kinetics_gammas << "\"" << "\n";

	// components
	s_oss << indent1;
	s_oss << "<component " << "\n";
	for (std::list < cxxKineticsComp >::const_iterator it =
		 kineticsComps.begin(); it != kineticsComps.end(); ++it)
	{
		it->dump_xml(s_oss, indent + 2);
	}

	return;
}
#endif

void
cxxKinetics::dump_raw(std::ostream & s_oss, unsigned int indent, int * n_out) const
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

	// Kinetics element and attributes
	s_oss << indent0;
	int n_user_local = (n_out != NULL) ? *n_out : this->n_user;
	s_oss << "KINETICS_RAW                 " << n_user_local << " " << this->description << "\n";

	s_oss << indent1 << "# KINETICS_MODIFY candidate identifiers #\n";

	s_oss << indent1;
	s_oss << "-step_divide               " << this->step_divide << "\n";

	s_oss << indent1;
	s_oss << "-rk                        " << this->rk << "\n";

	s_oss << indent1;
	s_oss << "-bad_step_max              " << this->bad_step_max << "\n";

	s_oss << indent1;
	s_oss << "-use_cvode                 " << this->use_cvode << "\n";

	s_oss << indent1;
	s_oss << "-cvode_steps               " << this->cvode_steps << "\n";

	s_oss << indent1;
	s_oss << "-cvode_order               " << this->cvode_order << "\n";

	// kineticsComps structures
	for (size_t k = 0; k < this->kinetics_comps.size(); k++)
	{
		s_oss << indent1;
		s_oss << "-component                 " << this->kinetics_comps[k].Get_rate_name() << "\n";
		this->kinetics_comps[k].dump_raw(s_oss, indent + 2);
	}

	// equal_steps
	s_oss << indent1;
	s_oss << "-equal_increments           " << this->equalIncrements << "\n";

	// equal_steps
	s_oss << indent1;
	s_oss << "-count                     " << this->count << "\n";

	// steps
	s_oss << indent1;
	s_oss << "-steps             " << "\n";
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

	s_oss << indent1 << "# KINETICS workspace variables #\n";
	// totals
	s_oss << indent1;
	s_oss << "-totals                    " << "\n";
	this->totals.dump_raw(s_oss, indent + 2);
	return;
}
void
cxxKinetics::read_raw(CParser & parser, bool check)
{

	LDBLE d;


	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	int opt_save;
	bool useLastLine(false);
	std::vector < LDBLE > temp_steps;

	// Read kinetics number and description
	this->read_number_description(parser);

	opt_save = CParser::OPT_ERROR;
	bool step_divide_defined(false);
	bool rk_defined(false);
	bool bad_step_max_defined(false);
	bool use_cvode_defined(false);
	bool cvode_steps_defined(false);
	bool cvode_order_defined(false);
	bool steps_defined(false);

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
			parser.error_msg("Unknown input in KINETICS_COMP_RAW keyword.",
							 PHRQ_io::OT_CONTINUE);
			parser.error_msg(parser.line().c_str(), PHRQ_io::OT_CONTINUE);
			break;

		case 0:				// step_divide
			if (!(parser.get_iss() >> this->step_divide))
			{
				this->step_divide = 1.0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for step_divide.",
								 PHRQ_io::OT_CONTINUE);
			}
			step_divide_defined = true;
			break;

		case 1:				// rk
			if (!(parser.get_iss() >> this->rk))
			{
				this->rk = 3;
				parser.incr_input_error();
				parser.error_msg("Expected integer value for rk.",
								 PHRQ_io::OT_CONTINUE);
			}
			rk_defined = true;
			break;

		case 2:				// bad_step_max
			if (!(parser.get_iss() >> this->bad_step_max))
			{
				this->bad_step_max = 500;
				parser.incr_input_error();
				parser.error_msg("Expected integer value for bad_step_max.",
								 PHRQ_io::OT_CONTINUE);
			}
			bad_step_max_defined = true;
			break;

		case 3:				// use_cvode
			if (!(parser.get_iss() >> this->use_cvode))
			{
				this->use_cvode = false;
				parser.incr_input_error();
				parser.error_msg("Expected boolean value for use_cvode.",
								 PHRQ_io::OT_CONTINUE);
			}
			use_cvode_defined = true;
			break;

		case 4:				// component
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
					cxxKineticsComp temp_comp(this->io);
					temp_comp.Set_rate_name(str.c_str());
					cxxKineticsComp *comp_ptr = this->Find(str);
					if (comp_ptr)
					{
						temp_comp = *comp_ptr;
					}
					temp_comp.read_raw(parser, false);
					if (comp_ptr)
					{
						for (size_t j = 0; j < this->kinetics_comps.size(); j++)
						{
							if (Utilities::strcmp_nocase(this->kinetics_comps[j].Get_rate_name().c_str(), str.c_str()) == 0)
							{
								this->kinetics_comps[j] = temp_comp;
							}
						}
					}
					else
					{
						this->kinetics_comps.push_back(temp_comp);
					}
					useLastLine = true;
				}
			}
			break;

		case 5:				// totals
			if (this->totals.read_raw(parser, next_char) !=
				CParser::PARSER_OK)
			{
				parser.incr_input_error();
				parser.
					error_msg
					("Expected element name and molality for KineticsComp totals.",
					 PHRQ_io::OT_CONTINUE);
			}
			opt_save = 5;
			break;

		case 6:				// steps
			while (parser.copy_token(token, next_char) == CParser::TT_DIGIT)
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
					temp_steps.push_back(d);
					steps_defined = true;
				}
			}
			opt_save = 6;
			break;

		case 7:				// cvode_steps
			if (!(parser.get_iss() >> this->cvode_steps))
			{
				this->cvode_steps = 100;
				parser.incr_input_error();
				parser.error_msg("Expected integer value for cvode_steps.",
								 PHRQ_io::OT_CONTINUE);
			}
			cvode_steps_defined = true;
			break;

		case 8:				// cvode_order
			if (!(parser.get_iss() >> this->cvode_order))
			{
				this->cvode_order = 5;
				parser.incr_input_error();
				parser.error_msg("Expected integer value for cvode_order.",
								 PHRQ_io::OT_CONTINUE);
			}
			cvode_order_defined = true;

			break;
		case 9:				// equalIncrements
		case 11:			// equal_increments
			if (!(parser.get_iss() >> this->equalIncrements))
			{
				this->use_cvode = false;
				parser.incr_input_error();
				parser.error_msg("Expected boolean value for equalIncrements.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;
		case 10:				// count
			if (!(parser.get_iss() >> this->count))
			{
				this->count = 0;
				parser.incr_input_error();
				parser.error_msg("Expected integer value for count.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;

		}
		if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD)
			break;
	}
	if (steps_defined)
	{
		this->steps = temp_steps;
	}
	if (check)
	{
		// members that must be defined
		if (step_divide_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Step_divide not defined for KINETICS_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (rk_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Rk not defined for KINETICS_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (bad_step_max_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Bad_step_max not defined for KINETICS_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (use_cvode_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Use_cvode not defined for KINETICS_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (cvode_steps_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Cvode_steps not defined for KINETICS_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (cvode_order_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Cvode_order not defined for KINETICS_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
	}
}
void
cxxKinetics::add(const cxxKinetics & addee, LDBLE extensive)
		//
		// Add to existing ppassemblage to "this" ppassemblage
		//
{
	if (extensive == 0.0)
		return;
	for (size_t i_add = 0; i_add < addee.kinetics_comps.size(); i_add++)
	{
		bool found(false);
		size_t i;
		for (i = 0; i < this->kinetics_comps.size(); i++)
		{
			if (this->kinetics_comps[i].Get_rate_name() == addee.kinetics_comps[i_add].Get_rate_name())
			{
				found = true;
				break;
			}
		}
		if (found)
		{
			this->kinetics_comps[i].add(addee.kinetics_comps[i_add], extensive);
		}
		else
		{
			cxxKineticsComp entity = addee.kinetics_comps[i_add];
			entity.multiply(extensive);
			this->kinetics_comps.push_back(entity);
		}
	}
	this->steps = addee.steps;
	this->step_divide = addee.step_divide;
	this->rk = addee.rk;
	this->bad_step_max = addee.bad_step_max;
	this->use_cvode = addee.use_cvode;
	this->cvode_steps = addee.cvode_steps;
	this->cvode_order = addee.cvode_order;
	this->equalIncrements = addee.equalIncrements;
	this->count = addee.count;
}
cxxKineticsComp * cxxKinetics::
Find(const std::string &s)
{
	for (size_t i = 0; i < this->kinetics_comps.size(); i++)
	{
		if (Utilities::strcmp_nocase(this->kinetics_comps[i].Get_rate_name().c_str(), s.c_str()) == 0)
		{
			return & (this->kinetics_comps[i]);
		}
	}
	return NULL;
}
int cxxKinetics::
Get_reaction_steps(void) const
{
	if(equalIncrements)
	{
		return count;
	}
	return ((int) steps.size());
}
LDBLE cxxKinetics::
Current_step(bool incremental_reactions, int reaction_step) const
{
	if (this->steps.size() == 0)
		return 1;
	LDBLE kin_time = 1;
	if (!incremental_reactions)
	{
		if (!this->equalIncrements)
		{
			if (reaction_step > (int) this->steps.size())
			{
				kin_time = this->steps[this->steps.size() - 1];
			}
			else
			{
				kin_time = this->steps[reaction_step - 1];
			}
		}
		else 
		{
			if (reaction_step > this->count)
			{
				kin_time = this->steps[0];
			}
			else
			{
				kin_time = (LDBLE) reaction_step * this->steps[0] /	((LDBLE) (this->count));
			}
		}
	}
	else
	{
		/* incremental reactions */
		if (!this->equalIncrements)
		{
			if (reaction_step > (int) this->steps.size())
			{
				kin_time = this->steps[this->steps.size() - 1];
			}
			else
			{
				kin_time = this->steps[reaction_step - 1];
			}
		}
		else 
		{
			if (reaction_step > this->count)
			{
				kin_time = 0;
			}
			else
			{
				kin_time = this->steps[0] / ((LDBLE) (this->count));
			}
		}
	}
	return kin_time;
}
void
cxxKinetics::Serialize(Dictionary & dictionary, std::vector < int >&ints, 
	std::vector < double >&doubles)
{
	ints.push_back(this->n_user);
	ints.push_back((int) this->kinetics_comps.size());
	for (size_t i = 0; i < this->kinetics_comps.size(); i++)
	{
		this->kinetics_comps[i].Serialize(dictionary, ints, doubles);
	}
	ints.push_back((int) this->steps.size());
	for (size_t i = 0; i < this->steps.size(); i++)
	{
		doubles.push_back(this->steps[i]);
	}
	ints.push_back(this->count);
	ints.push_back(this->equalIncrements ? 1 : 0);
	doubles.push_back(this->step_divide);
	ints.push_back(this->rk);
	ints.push_back(this->bad_step_max);
	ints.push_back(this->use_cvode ? 1 : 0);
	ints.push_back(this->cvode_steps);
	ints.push_back(this->cvode_order);
	this->totals.Serialize(dictionary, ints, doubles);
}

void
cxxKinetics::Deserialize(Dictionary & dictionary, std::vector < int >&ints, 
	std::vector < double >&doubles, int &ii, int &dd)
{
	this->n_user = ints[ii++];
	this->n_user_end = this->n_user;
	this->description = " ";

	int n = ints[ii++];
	this->kinetics_comps.clear();
	for (int i = 0; i < n; i++)
	{
		cxxKineticsComp kc;
		kc.Deserialize(dictionary, ints, doubles, ii, dd);
		this->kinetics_comps.push_back(kc);
	}
	n = ints[ii++];
	this->steps.clear();
	for (int i = 0; i < n; i++)
	{
		this->steps.push_back(doubles[dd++]);
	}
	this->count = ints[ii++];
	this->equalIncrements = (ints[ii++] != 0);
	this->step_divide = doubles[dd++];
	this->rk = ints[ii++];
	this->bad_step_max = ints[ii++];
	this->use_cvode = (ints[ii++] != 0);
	this->cvode_steps = ints[ii++];
	this->cvode_order = ints[ii++];
	this->totals.Deserialize(dictionary, ints, doubles, ii, dd);
}


const std::vector< std::string >::value_type temp_vopts[] = {
	std::vector< std::string >::value_type("step_divide"),             // 0 
	std::vector< std::string >::value_type("rk"),                      // 1 
	std::vector< std::string >::value_type("bad_step_max"),            // 2 
	std::vector< std::string >::value_type("use_cvode"),               // 3 
	std::vector< std::string >::value_type("component"),               // 4 
	std::vector< std::string >::value_type("totals"),                  // 5 
	std::vector< std::string >::value_type("steps"),                   // 6 
	std::vector< std::string >::value_type("cvode_steps"),             // 7 
	std::vector< std::string >::value_type("cvode_order"),             // 8 
	std::vector< std::string >::value_type("equalincrements"),         // 9 
	std::vector< std::string >::value_type("count"),                   // 10
	std::vector< std::string >::value_type("equal_increments")         // 11
};
const std::vector< std::string > cxxKinetics::vopts(temp_vopts, temp_vopts + sizeof temp_vopts / sizeof temp_vopts[0]);
