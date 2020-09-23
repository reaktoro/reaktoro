#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <cassert>
#include <sstream>				// std::ostrstream

#include "Utils.h"
#include "Phreeqc.h"
#include "SolutionIsotope.h"
#include "phqalloc.h"
#include "Dictionary.h"

cxxSolutionIsotope::cxxSolutionIsotope(PHRQ_io *io)
:
PHRQ_base(io),
isotope_number(0.0)
{
	isotope_number = 0;
	elt_name.clear();
	isotope_name.clear();
	total = 0;
	ratio = -9999.9;
	ratio_uncertainty = 1;
	ratio_uncertainty_defined = false;
	x_ratio_uncertainty = 0;
	coef = 0;
}
cxxSolutionIsotope::~cxxSolutionIsotope(void)
{
}

void
cxxSolutionIsotope::dump_xml(std::ostream & s_oss, unsigned int indent) const
{
	unsigned int i;

	std::string indent0(""), indent1("");
	for (i = 0; i < indent; ++i)
		indent0.append(Utilities::INDENT);
	for (i = 0; i < indent + 1; ++i)
		indent1.append(Utilities::INDENT);

	s_oss << indent0;
	s_oss << "<soln_isotope=\"" << "\n";

	s_oss << indent1;
	s_oss << "iso_isotope_number=\"" << this->
		isotope_number << "\"" << "\n";

	s_oss << indent1;
	s_oss << "iso_elt_name=\"" << this->elt_name << "\"" << "\n";

	s_oss << indent1;
	s_oss << "iso_isotope_name=\"" << this->isotope_name << "\"" << "\n";

	s_oss << indent1;
	s_oss << "iso_total=\"" << this->total << "\"" << "\n";

	s_oss << indent1;
	s_oss << "iso_ratio=\"" << this->ratio << "\"" << "\n";

	if (this->ratio_uncertainty != NAN)
	{
		s_oss << indent1;
		s_oss << "iso_ratio_uncertainty=\"" << this->
			ratio_uncertainty << "\"" << "\n";
	}
	s_oss << indent0;
	s_oss << "\">" << "\n";
}

void
cxxSolutionIsotope::dump_raw(std::ostream & s_oss, unsigned int indent) const
{
	unsigned int i;

	std::string indent0("");
	for (i = 0; i < indent; ++i)
		indent0.append(Utilities::INDENT);
	std::string indent1 = indent0;
	indent1.append(Utilities::INDENT);

	s_oss << indent0;
	s_oss << indent0 << "-isotope_number                    " << this->isotope_number << "\n";
	s_oss << indent0 << "-elt_name                          " << this->elt_name << "\n";
	s_oss << indent0 << "-total                             " << this->total << "\n";
	s_oss << indent0 << "-ratio                             " << this->ratio << "\n";
	s_oss << indent0 << "-ratio_uncertainty_defined         " << this->ratio_uncertainty_defined << "\n";
	if (this->ratio_uncertainty_defined)
	{
		s_oss << indent0 << "-ratio_uncertainty                 " << this->ratio_uncertainty << "\n";
	}
	s_oss << indent0 << "-x_ratio_uncertainty               " << this->x_ratio_uncertainty << "\n";
	s_oss << indent0 << "-coef                              " << this->coef << "\n";
}
void cxxSolutionIsotope::read_raw(CParser & parser, bool check )
{

	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	int opt_save;

	opt_save = CParser::OPT_ERROR;
	bool isotope_number_defined(false);
	bool elt_name_defined(false);
	bool total_defined(false);
	bool ratio_defined(false);

	for (;;)
	{
		int opt = parser.get_option(vopts, next_char);
		if (opt == CParser::OPT_DEFAULT)
		{
			opt = opt_save;
		}
		
		opt_save = CParser::OPT_DEFAULT;
		switch (opt)
		{
		case CParser::OPT_EOF:
			break;
		case CParser::OPT_KEYWORD:
			break;
		case CParser::OPT_DEFAULT:
		case CParser::OPT_ERROR:
			opt = CParser::OPT_EOF;
			parser.error_msg("Unknown input in isotopes of SOLUTION_RAW keyword.",
							 PHRQ_io::OT_CONTINUE);
			parser.error_msg(parser.line().c_str(), PHRQ_io::OT_CONTINUE);
			continue;

		case 0:				// isotope_number                  
			if (!(parser.get_iss() >> this->isotope_number))
			{
				this->isotope_number = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for isotope_number.",
								 PHRQ_io::OT_CONTINUE);
			}
			isotope_number_defined = true;
			break;

		case 1:				// elt_name                
			if (!(parser.get_iss() >> this->elt_name))
			{
				this->elt_name.clear();
				parser.incr_input_error();
				parser.error_msg("Expected character value for elt_name.",
								 PHRQ_io::OT_CONTINUE);
			}
			elt_name_defined = true;
			break;

		case 2:				// total
			if (!(parser.get_iss() >> this->total))
			{
				this->total = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for total.",
								 PHRQ_io::OT_CONTINUE);
			}
			total_defined = true;
			break;

		case 3:				// ratio
			if (!(parser.get_iss() >> this->ratio))
			{
				this->ratio = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for ratio.",
								 PHRQ_io::OT_CONTINUE);
			}
			total_defined = true;
			break;

		case 4:				// ratio_uncertainty_defined
			if (!(parser.get_iss() >> this->ratio_uncertainty_defined))
			{
				this->ratio_uncertainty_defined = 0;
				parser.incr_input_error();
				parser.error_msg("Expected boolean value for ratio_uncertainty_defined.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;

		case 5:				// ratio_uncertainty
			if (!(parser.get_iss() >> this->ratio_uncertainty))
			{
				this->ratio_uncertainty = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for ratio_uncertainty.",
								 PHRQ_io::OT_CONTINUE);
			}
			ratio_uncertainty_defined = true;
			break;

		case 6:				// x_ratio_uncertainty
			if (!(parser.get_iss() >> this->x_ratio_uncertainty))
			{
				this->x_ratio_uncertainty = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for x_ratio_uncertainty.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;
		case 7:				// coef
			if (!(parser.get_iss() >> this->coef))
			{
				this->coef = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for coef.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;
		}

		if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD)
			break;
	}
	if (check)
	{
		// all members must be defined
		if (isotope_number_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("isotope_number not defined for isotopes SOLUTION_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (elt_name_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("elt_name not defined for isotopes SOLUTION_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (total_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("total not defined for isotopes SOLUTION_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (ratio_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("ratio not defined for isotopes SOLUTION_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (ratio_uncertainty_defined == false)
		{
			parser.incr_input_error();
			parser.
				error_msg("ratio_uncertainty not defined for isotopes SOLUTION_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
	}
	return;
}
bool
cxxSolutionIsotope::operator<(const cxxSolutionIsotope & isotope) const
{
	int i = Utilities::strcmp_nocase(this->elt_name.c_str(), isotope.elt_name.c_str());
	if (i != 0)
		return (i < 0);
	return (this->isotope_number < isotope.isotope_number);
}

void
cxxSolutionIsotope::add(const cxxSolutionIsotope & isotope_ptr,
						LDBLE intensive, LDBLE extensive)
{
	if ((this->isotope_number == isotope_ptr.isotope_number) &&
		(this->elt_name == isotope_ptr.elt_name) &&
		(this->isotope_name == isotope_ptr.isotope_name))
	{
		this->total += isotope_ptr.total * extensive;
		this->ratio += isotope_ptr.ratio * intensive;
		this->ratio_uncertainty += isotope_ptr.ratio_uncertainty * intensive;
		this->ratio_uncertainty_defined = (this->ratio_uncertainty_defined
										   || isotope_ptr.
										   ratio_uncertainty_defined);
	}
}
void
cxxSolutionIsotope::multiply(LDBLE extensive)
{
	this->total *= extensive;
}
void 
cxxSolutionIsotope::Serialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles)
{
	doubles.push_back(this->isotope_number);
	ints.push_back(dictionary.Find(this->elt_name));
	ints.push_back(dictionary.Find(this->isotope_name));
	doubles.push_back(this->total);
	doubles.push_back(this->ratio);
	doubles.push_back(this->ratio_uncertainty);
	ints.push_back(this->ratio_uncertainty_defined ? 1 : 0);
	doubles.push_back(this->x_ratio_uncertainty);
	doubles.push_back(this->coef);
}
void 
cxxSolutionIsotope::Deserialize(Dictionary & dictionary, std::vector < int >&ints, std::vector < double >&doubles, int &ii, int &dd)
{
	this->isotope_number = doubles[dd++];
	this->elt_name = dictionary.GetWords()[ints[ii++]];
	this->isotope_name = dictionary.GetWords()[ints[ii++]];
	this->total = doubles[dd++];
	this->ratio = doubles[dd++];
	this->ratio_uncertainty = doubles[dd++];
	this->ratio_uncertainty_defined = (ints[ii++] != 0);
	this->x_ratio_uncertainty = doubles[dd++];
	this->coef = doubles[dd++];
}

const std::vector< std::string >::value_type temp_vopts[] = {
	std::vector< std::string >::value_type("isotope_number"),	            // 0 
	std::vector< std::string >::value_type("elt_name"),	                    // 1 
	std::vector< std::string >::value_type("total"),	                    // 2 
	std::vector< std::string >::value_type("ratio"),	                    // 3 
	std::vector< std::string >::value_type("ratio_uncertainty_defined"),	// 4 
	std::vector< std::string >::value_type("ratio_uncertainty"),	        // 5 
	std::vector< std::string >::value_type("x_ratio_uncertainty"),	        // 6 
	std::vector< std::string >::value_type("coef") 	                        // 7 
};									   
const std::vector< std::string > cxxSolutionIsotope::vopts(temp_vopts, temp_vopts + sizeof temp_vopts / sizeof temp_vopts[0]);
