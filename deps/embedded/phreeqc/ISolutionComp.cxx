#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <cassert>

#include "Utils.h"
#include "ISolutionComp.h"
#include "Parser.h"
#include "Solution.h"
#include "phqalloc.h"

cxxISolutionComp::cxxISolutionComp(PHRQ_io *io):
PHRQ_base(io),
moles(0.0), 
input_conc(0.0), 
phase_si(0.0),
gfw(0.0)
{
}
cxxISolutionComp::~cxxISolutionComp(void)
{
}

#ifdef SKIP_OR_MOVE_TO_STRUCTURES
struct conc *
cxxISolutionComp::cxxISolutionComp2conc(Phreeqc * phreeqc_ptr, const std::map < std::string,
										cxxISolutionComp > &totals)
		// for ISolutions
		// takes a std::vector cxxISolutionComp structures
		// returns list of conc structures
{
	struct conc *c;
	c = (struct conc *)
		phreeqc_ptr-> PHRQ_malloc((size_t) ((totals.size() + 1) * sizeof(struct conc)));
	if (c == NULL)
		phreeqc_ptr-> malloc_error();
	int i = 0;
	for (std::map < std::string, cxxISolutionComp >::const_iterator it = totals.begin();
		 it != totals.end(); ++it)
	{
		c[i].description = phreeqc_ptr-> string_duplicate(it->second.description.c_str());
		c[i].moles = it->second.moles;
		c[i].input_conc = it->second.input_conc;
		if (it->second.units.size() == 0)
			c[i].units = NULL;
		else
			c[i].units = phreeqc_ptr-> string_hsave(it->second.units.c_str());
		if (it->second.equation_name.size() == 0)
			c[i].equation_name = NULL;
		else
			c[i].equation_name = phreeqc_ptr-> string_hsave(it->second.equation_name.c_str());
		c[i].phase_si = it->second.phase_si;
		c[i].n_pe = it->second.n_pe;
		c[i].as = phreeqc_ptr-> string_hsave(it->second.as.c_str());
		c[i].gfw = it->second.gfw;
		//c[i].skip                = 0;
		c[i].phase = NULL;
		i++;
	}
	c[i].description = NULL;
	return (c);
}
#endif

#ifdef SKIP_OR_MOVE_TO_STRUCTURES
void
cxxISolutionComp::set_gfw(Phreeqc * phreeqc_ptr)
{
// return gfw
	if (this->gfw > 0.0)
		return;
// calculate gfw from as or from master species gfw
	if (this->as.size() != 0)
	{
		/* use given chemical formula to calculate gfw */
		LDBLE l_gfw;
		if (phreeqc_ptr-> compute_gfw(this->as.c_str(), &l_gfw) == ERROR)
		{
			std::ostringstream oss;
			oss << "Could not compute gfw, " << this->as;
			error_msg(oss.str().c_str(), CONTINUE);
			return;
		}
		//if (this->description == "Alkalinity" && this->as == "CaCO3") 
		if (strcmp(this->description.c_str(), "Alkalinity") == 0
			&& strcmp(this->as.c_str(), "CaCO3"))
		{
			l_gfw /= 2.;
		}
		this->gfw = l_gfw;
		return;
	}
	/* use gfw of master species */
	std::string str(this->description);
	struct master *master_ptr = phreeqc_ptr-> master_bsearch(str.c_str());
	if (master_ptr != NULL)
	{
		/* use gfw for element redox state */
		this->gfw = master_ptr->gfw;
		return;
	}
	std::ostringstream oss;
	oss << "Could not find gfw, " << this->description;
	error_msg(oss.str().c_str(), CONTINUE);
	return;
}
#endif
#ifdef SKIP
cxxISolutionComp::STATUS_TYPE cxxISolutionComp::read(CParser & parser,
													 cxxISolution & solution)
{
	// std::string& str = parser.line(); 
	std::string str = parser.line();

	// defaults set in ctor

	// Remove space between "kg" and "solution" or "water" in units
	Utilities::replace("Kg", "kg", str);
	Utilities::replace("KG", "kg", str);
	while (Utilities::replace("kg ", "kg", str));

	std::istream::pos_type ptr = 0;

	//
	// Read master species list for mass balance equation
	//
	std::string token;
	std::string token1;
	int
		count_redox_states = 0;
	CParser::TOKEN_TYPE j;
	while (((j = parser.copy_token(token, ptr)) == CParser::TT_UPPER) ||
		   (token[0] == '[') ||
		   (Utilities::strcmp_nocase_arg1(token.c_str(), "ph") == 0) ||
		   (Utilities::strcmp_nocase_arg1(token.c_str(), "pe") == 0))
	{
		++count_redox_states;
		Utilities::replace("(+", "(", token);
		if (count_redox_states > 1)
			token1 += " ";
		token1 += token;
	}
	if (count_redox_states == 0)
	{
		parser.incr_input_error();
		parser.
			error_msg
			("No element or master species given for concentration input.",
			 PHRQ_io::OT_CONTINUE);
		return cxxISolutionComp::ERROR;
	}
	description = token1;

	// Determine if reading alkalinity, allow equivalents for units
	Utilities::str_tolower(token1);
	bool
		alk = false;
	if (token1.find("alk") == 0)
	{
		alk = true;
	}

	// Read concentration
	if (!(std::istringstream(token) >> this->input_conc))
	{
		std::ostringstream err;
		err << "Concentration data error for " << token1 <<
			" in solution input.";
		parser.error_msg(err, PHRQ_io::OT_CONTINUE);
		return cxxISolutionComp::ERROR;
	}
	if ((j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY)
		return cxxISolutionComp::OK;

	// Read optional data
	token1 = token;

	// Check for units info
	if (parser.check_units(token1, alk, false, solution.get_units(), false) ==
		CParser::OK)
	{
		if (parser.
			check_units(token1, alk, false, solution.get_units(),
						true) == CParser::OK)
		{
			this->units = token1;
			if ((j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY)
				return cxxISolutionComp::OK;
		}
		else
		{
			return cxxISolutionComp::ERROR;
		}
	}

	// Check for "as" followed by formula to be used for gfw
	token1 = token;
	Utilities::str_tolower(token1);
	if (token1.compare("as") == 0)
	{
		parser.copy_token(token, ptr);
		this->as = token;
		if ((j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY)
			return cxxISolutionComp::OK;
	}
	// Check for "gfw" followed by gram formula weight
	else if (token1.compare("gfw") == 0)
	{
		if (parser.copy_token(token, ptr) != CParser::TT_DIGIT)
		{
			parser.error_msg("Expecting gram formula weight.",
							 PHRQ_io::OT_CONTINUE);
			return cxxISolutionComp::ERROR;
		}
		else
		{
			parser.get_iss() >> this->gfw;
			if ((j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY)
				return cxxISolutionComp::OK;
		}
	}

	// Check for redox couple for pe
	if (Utilities::strcmp_nocase_arg1(token.c_str(), "pe") == 0)
	{
		this->n_pe = cxxPe_Data::store(solution.pe, token);
		if ((j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY)
			return cxxISolutionComp::OK;
	}
	else if (token.find("/") != std::string::npos)
	{
		if (parser.parse_couple(token) == CParser::OK)
		{
			this->n_pe = cxxPe_Data::store(solution.pe, token);
			if ((j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY)
				return cxxISolutionComp::OK;
		}
		else
		{
			return cxxISolutionComp::ERROR;
		}
	}

	// Must have phase
	this->equation_name = token;
	if ((j = parser.copy_token(token, ptr)) == CParser::TT_EMPTY)
		return cxxISolutionComp::OK;

	// Check for saturation index
	if (!(std::istringstream(token) >> this->phase_si))
	{
		parser.error_msg("Expected saturation index.", PHRQ_io::OT_CONTINUE);
		return cxxISolutionComp::ERROR;
	}
	return cxxISolutionComp::OK;
}
#endif
#ifdef SKIP
void
cxxISolutionComp::dump_xml(std::ostream & s_oss, unsigned int indent) const const
{
	unsigned int i;
	std::string indent0("");
	for (i = 0; i < indent; ++i)
		indent0.append(Utilities::INDENT);

	s_oss << indent0;
	s_oss << "<soln_total";

	s_oss << " conc_desc=\"" << this->description << "\"";

	s_oss << " conc_moles=\"" << this->moles << "\"";

	s_oss << "\">" << "\n";
}
#endif
/* ---------------------------------------------------------------------- */
CParser::STATUS_TYPE cxxISolutionComp::
read(const char *line_in, cxxSolution *solution_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Remove space between "kg" and "solution" or "water" in units
 */
	std::string line = line_in;
	Utilities::replace("Kg", "kg", line);
	Utilities::replace("KG", "kg", line);
	while (Utilities::replace("kg ", "kg", line) == TRUE);
/*
 *   Read master species list for mass balance equation
 */
	std::string master_list;
	std::string token;
	std::string::iterator b = line.begin(); 
	std::string::iterator e = line.end(); 
	{
		int j;
		while (((j = CParser::copy_token(token, b, e)) == CParser::TT_UPPER) ||
			(token[0] == '[') ||
			(Utilities::strcmp_nocase(token.c_str(), "ph") == 0) ||
			(Utilities::strcmp_nocase(token.c_str(), "pe") == 0))
		{
			Utilities::replace("(+", "(", token);
			if (master_list.size() > 0)
				master_list.append(" ");
			master_list.append(token);
		}
	}
	if (master_list.size() == 0)
	{
		error_msg
			("No element or master species given for concentration input.",
			  PHRQ_io::OT_CONTINUE);
		return (CParser::PARSER_ERROR);
	}
	this->Set_description(master_list.c_str());
/*
 *   Determine if reading alkalinity, allow equivalents for units
 */
	bool alk;
	Utilities::str_tolower(master_list);
	if (strstr(master_list.c_str(), "alk") == master_list.c_str())
	{
		alk = true;
	}
	else
	{
		alk = false;
	}
/*
 *   Read concentration
 */
	{
		LDBLE dummy;
		int j = sscanf(token.c_str(), SCANFORMAT, &dummy);
		if (j == 0)
		{
			std::ostringstream errstr;
			errstr << "Concentration data error for " << master_list << " in solution input.";
			error_msg(errstr.str().c_str(),  PHRQ_io::OT_CONTINUE);
			return (CParser::PARSER_ERROR);
		}
		else
		{
			this->Set_input_conc(dummy);
		}
		if ((j = CParser::copy_token(token, b, e)) == CParser::TT_EMPTY)
			return (CParser::PARSER_OK);
	}
/*
 *   Read optional data
 */
	std::string token1 = token;

/*
 *   Check for units info
 */
	CParser parser(this->io);
	if (solution_ptr->Get_initial_data() == NULL)
	{
		error_msg("Initial_data instance not defined in cxxISolutionComp::read", 1);
	}
	if (parser.check_units(token1, alk, false, solution_ptr->Get_initial_data()->Get_units().c_str(), false) == CParser::PARSER_OK)
	{
		if (parser.check_units(token1, alk, false, solution_ptr->Get_initial_data()->Get_units().c_str(), true) == CParser::PARSER_OK)
		{
			this->units = token1;
			if ((CParser::copy_token(token, b, e)) == CParser::TT_EMPTY)
				return (CParser::PARSER_OK);
		}
		else
		{
			return (CParser::PARSER_ERROR);
		}
	}
/*
 *   Check for "as" followed by formula to be used for gfw
 */
	token1 = token;
	Utilities::str_tolower(token1);
	if (strcmp(token1.c_str(), "as") == 0)
	{
		CParser::copy_token(token, b, e);
		this->as = token;
		if ((CParser::copy_token(token, b, e)) == CParser::TT_EMPTY)
			return (CParser::PARSER_OK);
/*
 *   Check for "gfw" followed by gram formula weight
 */
	}
	else if (strcmp(token1.c_str(), "gfw") == 0 || strcmp(token1.c_str(), "gfm") == 0)
	{
		if (CParser::copy_token(token, b, e) != DIGIT)
		{
			error_msg("Expecting gram formula weight.",  PHRQ_io::OT_CONTINUE);
			return (CParser::PARSER_ERROR);
		}
		else
		{
			sscanf(token.c_str(), SCANFORMAT, &this->gfw);
			if ((CParser::copy_token(token, b, e)) == CParser::TT_EMPTY)
				return (CParser::PARSER_OK);
		}
	}
/*
 *   Check for redox couple for pe
 */
	if (Utilities::strcmp_nocase(token.c_str(), "pe") == 0)
	{
		this->pe_reaction = token;
		if ((CParser::copy_token(token, b, e)) == CParser::TT_EMPTY)
			return (CParser::PARSER_OK);
	}
	else if (strstr(token.c_str(), "/") != NULL)
	{
		if (parser.parse_couple(token) == CParser::PARSER_OK)
		{
			this->pe_reaction = token;
			if ((CParser::copy_token(token, b, e)) == CParser::TT_EMPTY)
				return (CParser::PARSER_OK);
		}
		else
		{
			return (CParser::PARSER_ERROR);
		}
	}
/*
 *   Must have phase
 */
	this->equation_name = token;
	if (CParser::copy_token(token, b, e) == CParser::TT_EMPTY)
		return (CParser::PARSER_OK);
/*
 *   Check for saturation index
 */
	{
		int j = sscanf(token.c_str(), SCANFORMAT,
			&(this->phase_si));
		if (j != 1)
		{
			error_msg("Expected saturation index.",  PHRQ_io::OT_CONTINUE);
			return (CParser::PARSER_ERROR);
		}
	}
	return (CParser::PARSER_OK);

}