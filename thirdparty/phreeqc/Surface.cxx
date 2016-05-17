// Surface.cxx: implementation of the cxxSurface class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <cassert>				// assert
#include <algorithm>			// std::sort
#include "Utils.h"				// define first
#include "Phreeqc.h"
#include "Surface.h"
#include "cxxMix.h"
#include "phqalloc.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxSurface::cxxSurface(PHRQ_io *io)
	//
	// default constructor for cxxSurface 
	//
:	cxxNumKeyword(io)
{
	new_def = false;
	type = DDL;
	dl_type = NO_DL;
	sites_units = SITES_ABSOLUTE;
	only_counter_ions = false;
	thickness = 1e-8;
	debye_lengths = 0.0;
	DDL_viscosity = 1.0;
	DDL_limit = 0.8;
	transport = false;
	solution_equilibria = false;
	n_solution = -999;
}
cxxSurface::cxxSurface(std::map < int, cxxSurface > &entities,
					   cxxMix & mix, int l_n_user, PHRQ_io *io):
cxxNumKeyword(io)
{
	this->n_user = this->n_user_end = l_n_user;
	this->new_def = false;
	type = DDL;
	dl_type = NO_DL;
	sites_units = SITES_ABSOLUTE;
	only_counter_ions = false;
	thickness = 1e-8;
	debye_lengths = 0.0;
	DDL_viscosity = 1.0;
	DDL_limit = 0.8;
	transport = false;
	solution_equilibria = false;
	n_solution = -999;
//
//   Mix surfaces
//
	const std::map < int, LDBLE >&mixcomps = mix.Get_mixComps();
	std::map < int, LDBLE >::const_iterator it;
	for (it = mixcomps.begin(); it != mixcomps.end(); it++)
	{
		if (entities.find(it->first) != entities.end())
		{
			this->add(entities.find(it->first)->second, it->second);
		}
	}
}

cxxSurface::~cxxSurface()
{
}

bool
cxxSurface::Get_related_phases() const
{
	for (size_t i = 0; i != this->surface_comps.size(); i++)
	{
		if (this->surface_comps[i].Get_phase_name().size() == 0)
			continue;
		return (true);
	}
	return (false);
}

bool
cxxSurface::Get_related_rate() const
{
	for (size_t i = 0; i != this->surface_comps.size(); i++)
	{
		if (this->surface_comps[i].Get_rate_name().size() == 0)
			continue;
		return (true);
	}
	return (false);
}
#ifdef SKIP
void
cxxSurface::dump_xml(std::ostream & s_oss, unsigned int indent) const
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

	// Surface element and attributes
	s_oss << indent0;
	s_oss << "<surface " << "\n";

	s_oss << indent1;
	s_oss << "surface_type=\"" << this->type << "\"" << "\n";

	s_oss << indent1;
	s_oss << "dl_type=\"" << this->dl_type << "\"" << "\n";

	s_oss << indent1;
	s_oss << "sites_units=\"" << this->sites_units << "\"" << "\n";

	s_oss << indent1;
	s_oss << "only_counter_ions=\"" << this->
		only_counter_ions << "\"" << "\n";

	s_oss << indent1;
	s_oss << "thickness=\"" << this->thickness << "\"" << "\n";

	s_oss << indent1;
	s_oss << "debye_lengths=\"" << this->debye_lengths << "\"" << "\n";

	s_oss << indent1;
	s_oss << "DDL_viscosity=\"" << this->DDL_viscosity << "\"" << "\n";

	s_oss << indent1;
	s_oss << "DDL_limit=\"" << this->DDL_limit << "\"" << "\n";

	s_oss << indent1;
	s_oss << "transport=\"" << this->transport << "\"" << "\n";

	// surface component structures
	s_oss << indent1;
	s_oss << "<component " << "\n";
	{
		for (std::map < std::string, cxxSurfaceComp >::const_iterator it =
			 this->surfaceComps.begin(); it != this->surfaceComps.end(); ++it)
		{
			(*it).second.dump_xml(s_oss, indent + 2);
		}
	}
	// surface charge structures
	s_oss << indent1;
	s_oss << "<charge_component " << "\n";
	for (std::map < std::string, cxxSurfaceCharge >::const_iterator it =
		 surface_charges.begin(); it != surface_charges.end(); ++it)
	{
		(*it).second.dump_xml(s_oss, indent + 2);
	}

	return;
}
#endif
void
cxxSurface::dump_raw(std::ostream & s_oss, unsigned int indent, int *n_out) const
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

	// Surface element and attributes
	s_oss << indent0;
	int n_user_local = (n_out != NULL) ? *n_out : this->n_user;
	s_oss << "SURFACE_RAW                  " << n_user_local << " " << this->description << "\n";
	s_oss << indent1 << "# SURFACE_MODIFY candidate identifiers #\n";
	s_oss << indent1;
	s_oss << "-type                      " << this->type << "\n";
	s_oss << indent1;
	s_oss << "-dl_type                   " << this->dl_type << "\n";
	s_oss << indent1;
	s_oss << "-only_counter_ions         " << this->only_counter_ions << "\n";
	s_oss << indent1;
	s_oss << "-thickness                 " << this->thickness << "\n";
	s_oss << indent1;
	s_oss << "-debye_lengths             " << this->debye_lengths << "\n";
	s_oss << indent1;
	s_oss << "-DDL_viscosity             " << this->DDL_viscosity << "\n";
	s_oss << indent1;
	s_oss << "-DDL_limit                 " << this->DDL_limit << "\n";
	// surfaceComps 
	for (size_t i = 0; i != this->surface_comps.size(); i++)
	{
		const cxxSurfaceComp * comp_ptr = &(this->surface_comps[i]);
		s_oss << indent1;
		s_oss << "-component                 " << comp_ptr->Get_formula() << "\n";
		comp_ptr->dump_raw(s_oss, indent + 2);
	}
	// surface charge 
	for (size_t i = 0; i != this->surface_charges.size(); i++)
	{
		const cxxSurfaceCharge * charge_ptr = &(this->surface_charges[i]);
		s_oss << indent1;
		s_oss << "-charge_component          " << charge_ptr->Get_name() << "\n";
		charge_ptr->dump_raw(s_oss, indent + 2);
	}

	s_oss << indent1 << "# SURFACE_MODIFY candidates with new_def=true #\n";
	s_oss << indent1;
	s_oss << "-new_def                   " << this->new_def << "\n";
	s_oss << indent1;
	s_oss << "-sites_units               " << this->sites_units << "\n";
	s_oss << indent1;
	s_oss << "-solution_equilibria       " << this->solution_equilibria << "\n";
	s_oss << indent1;
	s_oss << "-n_solution                " << this->n_solution << "\n";

	s_oss << indent1 << "# Surface workspace variables #\n";
	s_oss << indent1;
	s_oss << "-transport                 " << this->transport << "\n";
	s_oss << indent1;
	s_oss << "-totals                    " << "\n";
	this->totals.dump_raw(s_oss, indent + 2);

	return;
}

void
cxxSurface::read_raw(CParser & parser, bool check)
{
	int i = 0;
	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	bool useLastLine(false);

	// Read surface number and description
	this->read_number_description(parser);
	this->Set_new_def(false);

	bool only_counter_ions_defined(false);
	bool thickness_defined(false);
	bool type_defined(false);
	bool dl_type_defined(false);
	bool sites_units_defined(false);
	bool debye_lengths_defined(false);
	bool DDL_viscosity_defined(false);
	bool DDL_limit_defined(false);
	bool transport_defined(false);

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
			parser.error_msg("Unknown input in SURFACE keyword.",
							 PHRQ_io::OT_CONTINUE);
			parser.error_msg(parser.line().c_str(), PHRQ_io::OT_CONTINUE);
			break;

		case 0:				// diffuse_layer
			parser.incr_input_error();
			parser.error_msg("Diffuse layer is obsolete, use -type.",
							 PHRQ_io::OT_CONTINUE);
			break;

		case 1:				// edl
			parser.incr_input_error();
			parser.error_msg("-edl is obsolete, use -type.",
							 PHRQ_io::OT_CONTINUE);
			break;

		case 2:				// only_counter_ions
			if (!(parser.get_iss() >> this->only_counter_ions))
			{
				this->only_counter_ions = false;
				parser.incr_input_error();
				parser.
					error_msg("Expected boolean value for only_counter_ions.",
							  PHRQ_io::OT_CONTINUE);
			}
			only_counter_ions_defined = true;
			break;

		case 3:				// donnan
			parser.incr_input_error();
			parser.error_msg("-Donnan is obsolete, use -dl_type.",
							 PHRQ_io::OT_CONTINUE);
			break;

		case 4:				// thickness
			if (!(parser.get_iss() >> this->thickness))
			{
				this->thickness = 0.0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for thickness.",
								 PHRQ_io::OT_CONTINUE);
			}
			thickness_defined = true;
			break;
		case 5:				// component
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
					cxxSurfaceComp temp_comp(this->io);
					temp_comp.Set_formula(str.c_str());
					cxxSurfaceComp *comp_ptr = this->Find_comp(str);
					if (comp_ptr)
					{
						temp_comp = *comp_ptr;
					}
					temp_comp.read_raw(parser, check);
					if (comp_ptr)
					{
						for (size_t j = 0; j < this->surface_comps.size(); j++)
						{
							if (Utilities::strcmp_nocase(this->surface_comps[j].Get_formula().c_str(), str.c_str()) == 0)
							{
								this->surface_comps[j] = temp_comp;
							}
						}
					}
					else
					{
						this->surface_comps.push_back(temp_comp);
					}
					useLastLine = true;
				}
			}
			break;
		case 6:				// charge_component
			{
				std::string str;
				if (!(parser.get_iss() >> str))
				{
					parser.incr_input_error();
					parser.error_msg("Expected string value for charge name.",
									 PHRQ_io::OT_CONTINUE);
				}
				else
				{
					cxxSurfaceCharge temp_charge(this->io);
					temp_charge.Set_name(str.c_str());
					cxxSurfaceCharge *charge_ptr = this->Find_charge(str);
					if (charge_ptr)
					{
						temp_charge = *charge_ptr;
					}
					temp_charge.read_raw(parser, check);
					if (charge_ptr)
					{
						for (size_t j = 0; j < this->surface_charges.size(); j++)
						{
							if (Utilities::strcmp_nocase(this->surface_charges[j].Get_name().c_str(), str.c_str()) == 0)
							{
								this->surface_charges[j] = temp_charge;
							}
						}
					}
					else
					{
						this->surface_charges.push_back(temp_charge);
					}
					useLastLine = true;
				}
			}	
			useLastLine = true;
			break;
		case 7:				// type
			i = 0;
			if (!(parser.get_iss() >> i))
			{
				this->type = NO_EDL;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for type.",
								 PHRQ_io::OT_CONTINUE);
			}
			this->type = (SURFACE_TYPE) i;
			type_defined = true;
			break;
		case 8:				// dl_type
			i = 0;
			if (!(parser.get_iss() >> i))
			{
				this->dl_type = NO_DL;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for dl_type.",
								 PHRQ_io::OT_CONTINUE);
			}
			this->dl_type = (DIFFUSE_LAYER_TYPE) i;
			dl_type_defined = true;
			break;
		case 9:				// sites_units
			i = 0;
			if (!(parser.get_iss() >> i))
			{
				this->sites_units = SITES_ABSOLUTE;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for sites_units.",
								 PHRQ_io::OT_CONTINUE);
			}
			this->sites_units = (SITES_UNITS) i;
			sites_units_defined = true;
			break;

		case 10:				// debye_lengths
			if (!(parser.get_iss() >> this->debye_lengths))
			{
				this->debye_lengths = 0.0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for debye_lengths.",
								 PHRQ_io::OT_CONTINUE);
			}
			debye_lengths_defined = true;
			break;

		case 11:				// DDL_viscosity
			if (!(parser.get_iss() >> this->DDL_viscosity))
			{
				this->DDL_viscosity = 0.0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for DDL_viscosity.",
								 PHRQ_io::OT_CONTINUE);
			}
			DDL_viscosity_defined = true;
			break;

		case 12:				// DDL_limit
			if (!(parser.get_iss() >> this->DDL_limit))
			{
				this->DDL_limit = 0.0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for DDL_limit.",
								 PHRQ_io::OT_CONTINUE);
			}
			DDL_limit_defined = true;
			break;

		case 13:				// transport
			if (!(parser.get_iss() >> this->transport))
			{
				this->transport = false;
				parser.incr_input_error();
				parser.error_msg("Expected boolean value for transport.",
								 PHRQ_io::OT_CONTINUE);
			}
			transport_defined = true;
			break;

		case 14:				// new_def
			if (!(parser.get_iss() >> this->new_def))
			{
				this->new_def = false;
				parser.incr_input_error();
				parser.error_msg("Expected boolean value for new_def.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;

		case 15:				// new_def
			if (!(parser.get_iss() >> this->solution_equilibria))
			{
				this->solution_equilibria = false;
				parser.incr_input_error();
				parser.error_msg("Expected boolean value for solution_equilibria.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;

		case 16:				// n_solution
			if (!(parser.get_iss() >> this->n_solution))
			{
				this->n_solution = -999;
				parser.incr_input_error();
				parser.error_msg("Expected integer value for n_solution.",
								 PHRQ_io::OT_CONTINUE);
			}
			break;
		case 17:				// totals
			if (this->totals.read_raw(parser, next_char) !=	CParser::PARSER_OK)
			{
				parser.incr_input_error();
				parser.
					error_msg
					("Expected element name and molality for Surface totals.",
					 PHRQ_io::OT_CONTINUE);
			}
			break;
		}
		if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD)
			break;
	}
	if (check)
	{
		// members that must be defined
		if (only_counter_ions_defined == false)
		{
			parser.incr_input_error();
			parser.
				error_msg("Only_counter_ions not defined for SURFACE_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (thickness_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Thickness not defined for SURFACE_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (type_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Surface type not defined for SURFACE_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (dl_type_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Dl_type not defined for SURFACE_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (sites_units_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Sites_units not defined for SURFACE_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (debye_lengths_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Debye_lengths not defined for SURFACE_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (DDL_viscosity_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("DDL_viscosity not defined for SURFACE_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (DDL_limit_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("DDL_limit not defined for SURFACE_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (transport_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Transport not defined for SURFACE_RAW input.",
				PHRQ_io::OT_CONTINUE);
		}
	}
	this->Sort_comps();
}

void
cxxSurface::totalize()
{
	this->totals.clear();
	// component structures
	for (size_t i = 0; i != this->surface_comps.size(); i++)
	{
		cxxSurfaceComp * comp_ptr = &(this->surface_comps[i]);
		this->totals.add_extensive(comp_ptr->Get_totals(), 1.0);
		this->totals.add("Charge", comp_ptr->Get_charge_balance());
	}
	return;
}
void
cxxSurface::add(const cxxSurface & addee_in, LDBLE extensive)
		//
		// Add surface to "this" surface
		//
{
	cxxSurface addee = addee_in;
	if (extensive == 0.0)
		return;
	if (this->surface_comps.size() == 0)
	{
		this->only_counter_ions = addee.only_counter_ions;
		this->dl_type = addee.dl_type;
		this->type = addee.type;
		this->sites_units = addee.sites_units;
		this->thickness = addee.thickness;
		this->debye_lengths = addee.debye_lengths;
		this->DDL_viscosity = addee.DDL_viscosity;
		this->DDL_limit = addee.DDL_limit;
		this->solution_equilibria = addee.solution_equilibria;
		this->n_solution = addee.n_solution;
		this->transport = addee.transport;
	}

	for (size_t i_add = 0; i_add < addee.Get_surface_comps().size(); i_add++)
	{
		const cxxSurfaceComp & comp_add_ptr = addee.Get_surface_comps()[i_add];
		size_t i_this;
		for (i_this = 0; i_this < this->surface_comps.size(); i_this++)
		{
			cxxSurfaceComp & comp_this_ptr = this->surface_comps[i_this];
			if(comp_add_ptr.Get_formula() == comp_this_ptr.Get_formula())
			{
				comp_this_ptr.add(addee.surface_comps[i_add], extensive);
				break;
			}
		}
		if (i_this == this->surface_comps.size())
		{
			
			cxxSurfaceComp entity = comp_add_ptr;
			entity.multiply(extensive);
			this->surface_comps.push_back(entity);
		}
	}
	for (size_t i_add = 0; i_add < addee.Get_surface_charges().size(); i_add++)
	{
		const cxxSurfaceCharge & charge_add_ptr = addee.Get_surface_charges()[i_add];

		size_t i_this;
		for (i_this = 0; i_this < this->surface_charges.size(); i_this++)
		{
			cxxSurfaceCharge & charge_this_ptr = this->Get_surface_charges()[i_this];
			if(charge_add_ptr.Get_name() == charge_this_ptr.Get_name())
			{
				charge_this_ptr.add(addee.surface_charges[i_add], extensive);
				break;
			}
		}
		if (i_this == this->surface_charges.size())
		{
			
			cxxSurfaceCharge entity = charge_add_ptr;
			entity.multiply(extensive);
			this->surface_charges.push_back(entity);
		}
	}
}

void
cxxSurface::multiply(LDBLE extensive)
		//
		// Add surface to "this" surface
		//
{

	for (size_t i = 0; i < this->surface_comps.size(); i++)
	{
		cxxSurfaceComp *comp_ptr = &(this->surface_comps[i]);
		comp_ptr->multiply(extensive);
	}
	for (size_t i = 0; i < this->surface_charges.size(); i++)
	{
		cxxSurfaceCharge *charge_ptr = &(this->surface_charges[i]);
		charge_ptr->multiply(extensive);
	}
}
cxxSurfaceComp * cxxSurface::
Find_comp(std::string str)
{
		for (size_t i = 0; i < this->surface_comps.size(); i++)
		{
			if (Utilities::strcmp_nocase(str.c_str(), this->surface_comps[i].Get_formula().c_str()) == 0)
				return &(this->surface_comps[i]);
		}
		return NULL;
}

cxxSurfaceCharge * cxxSurface::
Find_charge(std::string str)
{
	for (size_t i = 0; i < this->surface_charges.size(); i++)
	{
		if (Utilities::strcmp_nocase(str.c_str(), this->surface_charges[i].Get_name().c_str()) == 0)
			return &(this->surface_charges[i]);
	}
	return NULL;
}
const cxxSurfaceCharge * cxxSurface::
Find_charge(std::string str)const
{
	for (size_t i = 0; i < this->surface_charges.size(); i++)
	{
		if (Utilities::strcmp_nocase(str.c_str(), this->surface_charges[i].Get_name().c_str()) == 0)
			return &(this->surface_charges[i]);
	}
	return NULL;
}
void cxxSurface::
Sort_comps(void)
{
	// sort comps
	{
		std::map<std::string, cxxSurfaceComp> comp_map;
		for (size_t i = 0; i < this->surface_comps.size(); i++)
		{
			comp_map[this->surface_comps[i].Get_formula()] = this->surface_comps[i];
		}
		this->surface_comps.clear();
		std::map<std::string, cxxSurfaceComp>::iterator it;
		for (it = comp_map.begin(); it != comp_map.end(); it++)
		{
			this->surface_comps.push_back(it->second);
		}
	}

	// sort charge too
	{
		std::map<std::string, cxxSurfaceCharge> charge_map;
		for (size_t i = 0; i < this->surface_charges.size(); i++)
		{
			charge_map[this->surface_charges[i].Get_name()] = this->surface_charges[i];
		}
		this->surface_charges.clear();
		std::map<std::string, cxxSurfaceCharge>::iterator it;
		for (it = charge_map.begin(); it != charge_map.end(); it++)
		{
			this->surface_charges.push_back(it->second);
		}
	}
}
/* ---------------------------------------------------------------------- */
void
cxxSurface::Serialize(Dictionary & dictionary, std::vector < int >&ints, 
	std::vector < double >&doubles)
/* ---------------------------------------------------------------------- */
{

	ints.push_back(this->n_user);
	ints.push_back((int) this->surface_comps.size());
	{
		for (size_t i = 0; i < this->surface_comps.size(); i++)
		{
			surface_comps[i].Serialize(dictionary, ints, doubles);
		}	
	}
	ints.push_back((int) this->surface_charges.size());
	{
		for (size_t i = 0; i < 	this->surface_charges.size(); i++)
		{
			surface_charges[i].Serialize(dictionary, ints, doubles);
		}
	}
	ints.push_back(this->new_def ? 1 : 0);
	ints.push_back((int) this->type);
	ints.push_back((int) this->dl_type);
	ints.push_back((int) this->sites_units);
	ints.push_back(this->only_counter_ions ? 1 : 0);
	doubles.push_back(this->thickness);
	doubles.push_back(this->debye_lengths);
	doubles.push_back(this->DDL_viscosity);
	doubles.push_back(this->DDL_limit);
	ints.push_back(this->transport ? 1 : 0);
	this->totals.Serialize(dictionary, ints, doubles);
	ints.push_back(this->solution_equilibria ? 1 : 0);
	ints.push_back((int) this->n_solution);

}

/* ---------------------------------------------------------------------- */
void
cxxSurface::Deserialize(Dictionary & dictionary, std::vector < int >&ints, 
	std::vector < double >&doubles, int &ii, int &dd)
/* ---------------------------------------------------------------------- */
{
	this->n_user = ints[ii++];
	this->n_user_end = this->n_user;
	this->description = " ";
	{
		int count = ints[ii++];
		this->surface_comps.clear();
		for (int n = 0; n < count; n++)
		{
			cxxSurfaceComp sc;
			sc.Deserialize(dictionary, ints, doubles, ii, dd);
			this->surface_comps.push_back(sc);
		}
	}
	{
		int count = ints[ii++];
		this->surface_charges.clear();
		for (int n = 0; n < count; n++)
		{
			cxxSurfaceCharge sc;
			sc.Deserialize(dictionary, ints, doubles, ii, dd);
			this->surface_charges.push_back(sc);
		}
	}
	this->new_def = (ints[ii++] != 0);
	this->type = (SURFACE_TYPE) ints[ii++];
	this->dl_type = (DIFFUSE_LAYER_TYPE) ints[ii++];
	this->sites_units = (SITES_UNITS) ints[ii++];
	this->only_counter_ions = (ints[ii++] != 0);
	this->thickness = doubles[dd++];
	this->debye_lengths = doubles[dd++];
	this->DDL_viscosity = doubles[dd++];
	this->DDL_limit = doubles[dd++];
	this->transport = (ints[ii++] != 0);
	this->totals.Deserialize(dictionary, ints, doubles, ii, dd);
	this->solution_equilibria = (ints[ii++] != 0);
	this->n_solution = ints[ii++];

}


const std::vector< std::string >::value_type temp_vopts[] = {
	std::vector< std::string >::value_type("diffuse_layer"),	    // 0 
	std::vector< std::string >::value_type("edl"),	                // 1 
	std::vector< std::string >::value_type("only_counter_ions"),	// 2 
	std::vector< std::string >::value_type("donnan"),	            // 3 
	std::vector< std::string >::value_type("thickness"),	        // 4 
	std::vector< std::string >::value_type("component"),	        // 5
	std::vector< std::string >::value_type("charge_component"),	    // 6 
	std::vector< std::string >::value_type("type "),	            // 7
	std::vector< std::string >::value_type("dl_type"),	            // 8 
	std::vector< std::string >::value_type("sites_units"),	        // 9 
	std::vector< std::string >::value_type("debye_lengths"),	    // 10
	std::vector< std::string >::value_type("ddl_viscosity"),	    // 11
	std::vector< std::string >::value_type("ddl_limit"),	        // 12
	std::vector< std::string >::value_type("transport"),	        // 13
	std::vector< std::string >::value_type("new_def"),	            // 14
	std::vector< std::string >::value_type("solution_equilibria"),	// 15
	std::vector< std::string >::value_type("n_solution"),	        // 16
	std::vector< std::string >::value_type("totals") 	            // 17
};									   
const std::vector< std::string > cxxSurface::vopts(temp_vopts, temp_vopts + sizeof temp_vopts / sizeof temp_vopts[0]);	
