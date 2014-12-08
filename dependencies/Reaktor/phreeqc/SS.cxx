// SSassemblageSS.cxx: implementation of the cxxSS class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <cassert>				// assert
#include <algorithm>			// std::sort

#include "Phreeqc.h"
#include "Utils.h"				// define first
#include "SS.h"
//#include "Dictionary.h"
#include "phqalloc.h"



//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxSS::cxxSS(PHRQ_io *io)
:
PHRQ_base(io)
	//
	// default constructor for cxxSS 
	//
{

	total_moles = 0;
	dn = 0;
	a0 = 0;
	a1 = 0;
	ag0 = 0;
	ag1 = 0;
	ss_in = false;
	miscibility = false;
	spinodal = false;
	tk = 298.15;
	xb1 = 0;
	xb2 = 0;
	input_case = SS_PARM_NONE;
	for (int i = 0; i < 4; i++)
	{
		p.push_back(0);
	}
}
cxxSS::~cxxSS()
{
}

#ifdef SKIP
void
cxxSS::dump_xml(std::ostream & s_oss, unsigned int indent) const const
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

	// S_S element and attributes

	s_oss << indent0 << "name=\"" << this->name << "\"" << "\n";
	s_oss << indent0 << "add_formula=\"" << this->
		add_formula << "\"" << "\n";
	s_oss << indent0 << "si=\"" << this->si << "\"" << "\n";
	s_oss << indent0 << "moles=\"" << this->moles << "\"" << "\n";
	s_oss << indent0 << "delta=\"" << this->delta << "\"" << "\n";
	s_oss << indent0 << "initial_moles=\"" << this->
		initial_moles << "\"" << "\n";
	s_oss << indent0 << "dissolve_only=\"" << this->
		dissolve_only << "\"" << "\n";

}
#endif
void
cxxSS::dump_raw(std::ostream & s_oss, unsigned int indent) const
{
	unsigned int i;
	s_oss.precision(DBL_DIG - 1);
	std::string indent0(""), indent1("");
	for (i = 0; i < indent; ++i)
		indent0.append(Utilities::INDENT);
	for (i = 0; i < indent + 1; ++i)
		indent1.append(Utilities::INDENT);

	s_oss << indent0 << "# SOLID_SOLUTION_MODIFY candidate identifiers #\n";
	for (size_t i = 0; i < this->ss_comps.size(); i++)
	{
		s_oss << indent0 << "-component               " << this->ss_comps[i].Get_name() << "\n";
		this->ss_comps[i].dump_raw(s_oss, indent + 1);
	}

	s_oss << indent0 << "# SOLID_SOLUTION_MODIFY candidate identifiers with new_def=true #\n";
	s_oss << indent0 << "-tk                      " << this->tk << "\n";
	s_oss << indent0 << "-input_case              " << this->input_case << "\n";
	s_oss << indent0 << "-p			              " << 
		p[0] << "\t"  << 
		p[1] << "\t"  << 
		p[2] << "\t"  << 
		p[3] << "\n";

	s_oss << indent0 << "# solid solution workspace variables #\n";
	s_oss << indent0 << "-ag0                     " << this->ag0 << "\n";
	s_oss << indent0 << "-ag1                     " << this->ag1 << "\n";
	s_oss << indent0 << "-a0                      " << this->a0 << "\n";
	s_oss << indent0 << "-a1                      " << this->a1 << "\n";
	s_oss << indent0 << "-xb1                     " << this->xb1 << "\n";
	s_oss << indent0 << "-xb2                     " << this->xb2 << "\n";
	s_oss << indent0 << "-miscibility             " << this->miscibility << "\n";
	s_oss << indent0 << "-spinodal                " << this->spinodal << "\n";
	s_oss << indent0 << "-ss_in                   " << this->ss_in << "\n";
	s_oss << indent0 << "-total_moles             " << this->total_moles    << "\n";
	s_oss << indent0 << "-dn                      " << this->dn             << "\n";
	s_oss << indent0 << "-totals                  " << "\n";
	this->totals.dump_raw(s_oss, indent + 1);
}

void
cxxSS::read_raw(CParser & parser, bool check)
{
	std::string str;

	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	int opt_save;

	opt_save = CParser::OPT_ERROR;
	bool a0_defined(false);
	bool a1_defined(false);
	bool ag0_defined(false);
	bool ag1_defined(false);
	bool miscibility_defined(false);
	bool xb1_defined(false);
	bool xb2_defined(false);
	bool useLastLine = false;

	for (size_t i = this->Get_p().size(); i < 4; i++)
	{
		this->p.push_back(0.0);
	}

	for (;;)
	{
		int opt;
		if (useLastLine == false)
		{
			opt = parser.get_option(vopts, next_char);
		}
		else
		{
			opt = parser.getOptionFromLastLine(vopts, next_char, false);
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
			opt = CParser::OPT_KEYWORD;
			// Allow return to Exchange for more processing
			break;

		case 0:				// name not used
			parser.warning_msg("-ss_name not used, defined in -solid_solution");
			opt_save = CParser::OPT_DEFAULT;
			break;

		case 1:				// total_moles
			if (!(parser.get_iss() >> this->total_moles))
			{
				this->total_moles = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for total_moles.", PHRQ_io::OT_CONTINUE);
			}
			opt_save = CParser::OPT_DEFAULT;
			break;

		case 2:				// a0
			if (!(parser.get_iss() >> this->a0))
			{
				this->a0 = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for a0.",
								 PHRQ_io::OT_CONTINUE);
			}
			a0_defined = true;
			opt_save = CParser::OPT_DEFAULT;
			break;

		case 3:				// a1
			if (!(parser.get_iss() >> this->a1))
			{
				this->a1 = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for a1.",
								 PHRQ_io::OT_CONTINUE);
			}
			a1_defined = true;
			opt_save = CParser::OPT_DEFAULT;
			break;

		case 4:				// components
		case 12:			// component
			{
				if (!(parser.get_iss() >> str))
				{
					this->name.clear();
					parser.incr_input_error();
					parser.error_msg("Expected string value for component name.",
									 PHRQ_io::OT_CONTINUE);
				}
				else
				{
					cxxSScomp temp_comp(io);
					temp_comp.Set_name(str);
					cxxSScomp * comp_ptr = this->Find(str.c_str());
					if (comp_ptr)
					{
						temp_comp = *comp_ptr;	
					}
					temp_comp.read_raw(parser, false);
					if (comp_ptr)
					{
						for (size_t j = 0; j < this->ss_comps.size(); j++)
						{
							if (Utilities::strcmp_nocase(this->ss_comps[j].Get_name().c_str(), str.c_str()) == 0)
							{
								this->ss_comps[j] = temp_comp;
							}
						}
					}
					else
					{
						this->ss_comps.push_back(temp_comp);
					}
					useLastLine = true;
				}
			}
			
			opt_save = CParser::OPT_DEFAULT;
			break;

		case 5:				// miscibility
			if (!(parser.get_iss() >> this->miscibility))
			{
				this->miscibility = 0;
				parser.incr_input_error();
				parser.error_msg("Expected boolean value for miscibility.",
								 PHRQ_io::OT_CONTINUE);
			}
			miscibility_defined = true;
			opt_save = CParser::OPT_DEFAULT;
			break;

		case 6:				// spinodal
			if (!(parser.get_iss() >> this->spinodal))
			{
				this->spinodal = 0;
				parser.incr_input_error();
				parser.error_msg("Expected boolean value for spinodal.", PHRQ_io::OT_CONTINUE);
			}
			opt_save = CParser::OPT_DEFAULT;
			break;

		case 7:				// tk
			if (!(parser.get_iss() >> this->tk))
			{
				this->tk = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for tk.", PHRQ_io::OT_CONTINUE);
			}
			opt_save = CParser::OPT_DEFAULT;
			break;

		case 8:				// xb1
			if (!(parser.get_iss() >> this->xb1))
			{
				this->xb1 = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for xb1.",
								 PHRQ_io::OT_CONTINUE);
			}
			xb1_defined = true;
			opt_save = CParser::OPT_DEFAULT;
			break;

		case 9:				// xb2
			if (!(parser.get_iss() >> this->xb2))
			{
				this->xb2 = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for xb2.",
								 PHRQ_io::OT_CONTINUE);
			}
			xb2_defined = true;
			opt_save = CParser::OPT_DEFAULT;
			break;

		case 10:				// ag0
			if (!(parser.get_iss() >> this->ag0))
			{
				this->ag0 = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for ag0.",
								 PHRQ_io::OT_CONTINUE);
			}
			ag0_defined = true;
			opt_save = CParser::OPT_DEFAULT;
			break;

		case 11:				// ag1
			if (!(parser.get_iss() >> this->ag1))
			{
				this->ag1 = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for ag1.",
								 PHRQ_io::OT_CONTINUE);
			}
			ag1_defined = true;
			opt_save = CParser::OPT_DEFAULT;
			break;
		case 13:				// input_case
			{
				int i;
				if (!(parser.get_iss() >> i))
				{
					this->input_case = cxxSS::SS_PARM_NONE;
					parser.incr_input_error();
					parser.error_msg("Expected integer value for parameter type.", PHRQ_io::OT_CONTINUE);
				}
				else
				{
					this->input_case = (cxxSS::SS_PARAMETER_TYPE) i;
				}
			}
			opt_save = CParser::OPT_DEFAULT;
			break;
		case 14:				// p
			{
				this->p.clear();
				this->p.assign(4, 0.0);
				for (int i = 0; i < 4; i++)
				{
					if (!(parser.get_iss() >> this->p[i]))
					parser.error_msg("Expected 4 parameters.");
				}
			}
			break;

		case 15:				// ss_in
			{
				if (!(parser.get_iss() >> this->ss_in))
				{
					this->ss_in = false;
					parser.incr_input_error();
					parser.error_msg("Expected boolean value for ss_in.",
						PHRQ_io::OT_CONTINUE);
				}
			}
			break;
		case 16:				// totals
			if (this->totals.read_raw(parser, next_char) !=	CParser::PARSER_OK)
			{
				parser.incr_input_error();
				parser.error_msg
					("Expected element name and molality for cxxSS totals.",
					 PHRQ_io::OT_CONTINUE);
			}
			opt_save = 16;
			break;
		case 17:				// dn
			if (!(parser.get_iss() >> this->dn))
			{
				this->dn = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for dn.", PHRQ_io::OT_CONTINUE);
			}
			opt_save = CParser::OPT_DEFAULT;
			break;
		}
		if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD)
			break;
	}
	if (check)
	{
		// members that must be defined
		if (a0_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("A0 not defined for SSassemblageSS input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (a1_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("A1 not defined for SSassemblageSS input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (ag0_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Ag0 not defined for SSassemblageSS input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (ag1_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Ag1 not defined for SSassemblageSS input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (miscibility_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Miscibility not defined for SSassemblageSS input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (xb1_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Xb1 not defined for SSassemblageSS input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (xb2_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Xb2 not defined for SSassemblageSS input.",
				PHRQ_io::OT_CONTINUE);
		}
	}
}

void
cxxSS::totalize(Phreeqc * phreeqc_ptr)
{
	this->totals.clear();
	// component structures
	for (size_t i = 0; i < this->ss_comps.size(); i++)
	{
		struct phase *phase_ptr;
		int l;
		phase_ptr = phreeqc_ptr-> phase_bsearch(ss_comps[i].Get_name().c_str(), &l, FALSE);
		if (phase_ptr != NULL)
		{
			cxxNameDouble phase_formula(phase_ptr->next_elt);
			this->totals.add_extensive(phase_formula, ss_comps[i].Get_moles());
		}
		else
		{
			assert(false);
		}
	}
	return;
}
void
cxxSS::add(const cxxSS & addee_in, LDBLE extensive)
{
	if (extensive == 0.0)
		return;
	if (addee_in.name.size() == 0)
		return;
	cxxSS addee = addee_in;
	for (size_t j = 0; j < addee.Get_ss_comps().size(); j++)
	{
		size_t i; 
		for (i = 0; i < this->ss_comps.size(); i++)
		{
			//look for component in this
			if (Utilities::strcmp_nocase(this->ss_comps[i].Get_name().c_str(), 
				addee.Get_ss_comps()[j].Get_name().c_str()) == 0)
			{
				LDBLE d = this->ss_comps[i].Get_initial_moles() + 
					addee.Get_ss_comps()[j].Get_initial_moles() * extensive;
				this->ss_comps[i].Set_initial_moles(d);
				d = this->ss_comps[i].Get_moles() + 
					addee.Get_ss_comps()[j].Get_moles() * extensive;
				this->ss_comps[i].Set_moles(d);
				d = this->ss_comps[i].Get_init_moles() + 
					addee.Get_ss_comps()[j].Get_init_moles() * extensive;
				this->ss_comps[i].Set_init_moles(d);
				d = this->ss_comps[i].Get_delta() + 
					addee.Get_ss_comps()[j].Get_delta() * extensive;
				this->ss_comps[i].Set_delta(d);
				break;
			}
		}
		if (i == this->ss_comps.size())
		{
			cxxSScomp comp = addee.Get_ss_comps()[j];
			comp.multiply(extensive);
			this->Get_ss_comps().push_back(comp);
		}
	}
}
void
cxxSS::multiply(LDBLE extensive)
{
	size_t i; 
	for (i = 0; i < this->ss_comps.size(); i++)
	{
		LDBLE d = this->ss_comps[i].Get_initial_moles() * extensive;
		this->ss_comps[i].Set_initial_moles(d);
		d = this->ss_comps[i].Get_moles() * extensive;
		this->ss_comps[i].Set_moles(d);
		d = this->ss_comps[i].Get_init_moles() * extensive;
		this->ss_comps[i].Set_init_moles(d);
		d = this->ss_comps[i].Get_delta() * extensive;
		this->ss_comps[i].Set_delta(d);
	}
}
cxxSScomp *
cxxSS::Find(const char * comp_name)
{

	for (size_t i = 0; i < this->ss_comps.size(); i++)
	{
		if (ss_comps[i].Get_name() == comp_name)
			return &(ss_comps[i]);
	}
	return NULL;
}
const std::vector< std::string >::value_type temp_vopts[] = {
	std::vector< std::string >::value_type("ss_name"),	    // 0                                   
	std::vector< std::string >::value_type("total_moles"),	// 1   
	std::vector< std::string >::value_type("a0"),           // 2   
	std::vector< std::string >::value_type("a1"),	        // 3
	std::vector< std::string >::value_type("components"),	// 4
	std::vector< std::string >::value_type("miscibility"),	// 5
	std::vector< std::string >::value_type("spinodal"),	    // 6
	std::vector< std::string >::value_type("tk"),	        // 7
	std::vector< std::string >::value_type("xb1"),	        // 8
	std::vector< std::string >::value_type("xb2"),	        // 9
	std::vector< std::string >::value_type("ag0"),	        // 10
	std::vector< std::string >::value_type("ag1"),	        // 11
	std::vector< std::string >::value_type("component"),	// 12
	std::vector< std::string >::value_type("input_case"),   //13
	std::vector< std::string >::value_type("p"),            //14
	std::vector< std::string >::value_type("ss_in"),        //15
	std::vector< std::string >::value_type("totals"),       //16
	std::vector< std::string >::value_type("dn")            //17
};									   
const std::vector< std::string > cxxSS::vopts(temp_vopts, temp_vopts + sizeof temp_vopts / sizeof temp_vopts[0]);	