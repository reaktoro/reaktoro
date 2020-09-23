// SurfaceCharge.cxx: implementation of the cxxSurfaceCharge class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <cassert>				// assert
#include <algorithm>			// std::sort

#include "Utils.h"				// define first
#include "Phreeqc.h"
#include "SurfaceCharge.h"
#include "phqalloc.h"
#include "Dictionary.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

cxxSurfaceCharge::cxxSurfaceCharge(PHRQ_io *io)
:
PHRQ_base(io)
//
// default constructor for cxxSurfaceCharge 
//
{
	specific_area = 0.0;
	grams = 0.0;
	charge_balance = 0.0;
	mass_water = 0.0;
	la_psi = 0.0;
	capacitance[0] = 1.0;
	capacitance[1] = 5.0;
	sigma0 = sigma1 = sigma2 = sigmaddl = 0;
	diffuse_layer_totals.type = cxxNameDouble::ND_ELT_MOLES;
}
cxxSurfaceCharge::~cxxSurfaceCharge()
{
}

void
cxxSurfaceCharge::dump_xml(std::ostream & s_oss, unsigned int indent) const
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

	// Surf_Charge element and attributes
	s_oss << indent0 << "name=\"" << this->name << "\"" << "\n";
	s_oss << indent0 << "specific_area=\"" << this->
		specific_area << "\"" << "\n";
	s_oss << indent0 << "grams=\"" << this->grams << "\"" << "\n";
	s_oss << indent0 << "charge_balance=\"" << this->
		charge_balance << "\"" << "\n";
	s_oss << indent0 << "mass_water=\"" << this->
		mass_water << "\"" << "\n";
	s_oss << indent0 << "la_psi=\"" << this->la_psi << "\"" << "\n";
	s_oss << indent0 << "capacitance=\"" << this->
		capacitance[0] << " " << this->capacitance[0] << "\"" << "\n";

	// totals
	s_oss << indent0;
	s_oss << "<diffuse_layer_totals " << "\n";
	this->diffuse_layer_totals.dump_xml(s_oss, indent + 1);

}

void
cxxSurfaceCharge::dump_raw(std::ostream & s_oss, unsigned int indent) const
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

	// Surf_Charge element and attributes
	s_oss << indent0 << "# SURFACE_MODIFY candidate identifiers #\n";
	s_oss << indent0 << "-specific_area           " << this->specific_area << "\n";
	s_oss << indent0 << "-grams                   " << this->grams << "\n";
	s_oss << indent0 << "-charge_balance          " << this->charge_balance << "\n";
	s_oss << indent0 << "-mass_water              " << this->mass_water << "\n";
	s_oss << indent0 << "-la_psi                  " << this->la_psi << "\n";
	s_oss << indent0 << "-capacitance0            " << this->capacitance[0] << "\n";
	s_oss << indent0 << "-capacitance1            " << this->capacitance[1] << "\n";
	// totals
	s_oss << indent0;
	s_oss << "-diffuse_layer_totals" << "\n";
	this->diffuse_layer_totals.dump_raw(s_oss, indent + 1);

	// DL species
	//s_oss << indent0;
	//s_oss << "-diffuse_layer_species" << "\n";
	if (dl_species_map.size() > 0)
	{
		s_oss << indent0;
		s_oss << "-diffuse_layer_species" << "\n";
		std::map<int, double>::const_iterator it = this->dl_species_map.begin();
		for ( ; it != dl_species_map.end(); it++)
		{
			s_oss << indent1;
			s_oss << it->first << " " << it->second << "\n";
		}
	}
	s_oss << indent0 << "# Surface workspace variables #\n";
	s_oss << indent0 << "-sigma0                  " << this->sigma0 << "\n";
	s_oss << indent0 << "-sigma1                  " << this->sigma1 << "\n";
	s_oss << indent0 << "-sigma2                  " << this->sigma2 << "\n";
	s_oss << indent0 << "-sigmaddl                " << this->sigmaddl << "\n";
	std::map<LDBLE, cxxSurfDL>::const_iterator git;
	for (git = this->g_map.begin(); git != g_map.end(); git++)
	{
		s_oss << indent0 << "-g_map                   " << git->first << "\t";
		s_oss << git->second.Get_g() << "\t";
		s_oss << git->second.Get_dg() << "\t";
		s_oss << git->second.Get_psi_to_z() << "\n";
	}
}

void
cxxSurfaceCharge::read_raw(CParser & parser, bool check)
{
	std::string str;


	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	int opt_save;

	opt_save = CParser::OPT_ERROR;
	bool specific_area_defined(false);
	bool grams_defined(false);
	bool charge_balance_defined(false);
	bool mass_water_defined(false);
	bool la_psi_defined(false);
	bool capacitance0_defined(false);
	bool capacitance1_defined(false);
	bool g_map_first(true);

	for (;;)
	{
		int opt = parser.get_option(vopts, next_char);
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
		case CParser::OPT_DEFAULT:
		case CParser::OPT_ERROR:
			opt = CParser::OPT_KEYWORD;
			// Allow return to Surface for more processing
			break;

		case 0:				// name
			warning_msg("-name ignored. Defined with -charge_component.");
			break;

		case 1:				// specific_area
			if (!(parser.get_iss() >> this->specific_area))
			{
				this->specific_area = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for specific_area.",
					PHRQ_io::OT_CONTINUE);
			}
			specific_area_defined = true;
			break;

		case 2:				// grams
			if (!(parser.get_iss() >> this->grams))
			{
				this->grams = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for grams.",
					PHRQ_io::OT_CONTINUE);
			}
			grams_defined = true;
			break;


		case 3:				// charge_balance
			if (!(parser.get_iss() >> this->charge_balance))
			{
				this->charge_balance = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for charge_balance.",
					PHRQ_io::OT_CONTINUE);
			}
			charge_balance_defined = true;
			break;

		case 4:				// mass_water
			if (!(parser.get_iss() >> this->mass_water))
			{
				this->mass_water = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for mass_water.",
					PHRQ_io::OT_CONTINUE);
			}
			mass_water_defined = true;
			break;


		case 5:				// la_psi
			if (!(parser.get_iss() >> this->la_psi))
			{
				this->la_psi = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for la_psi.",
					PHRQ_io::OT_CONTINUE);
			}
			la_psi_defined = true;
			break;


		case 6:				// diffuse_layer_totals
			if (this->diffuse_layer_totals.read_raw(parser, next_char) !=
				CParser::PARSER_OK)
			{
				parser.incr_input_error();
				parser.
					error_msg
					("Expected element name and molality for SurfaceCharge diffuse_layer_totals.",
					PHRQ_io::OT_CONTINUE);
			}
			opt_save = 6;
			break;

		case 7:				// la_psi1
			parser.warning_msg("-la_psi1 identifier not used");
			break;

		case 8:				// la_psi2
			parser.warning_msg("-la_psi2 identifier not used");
			break;

		case 9:				// capacitance0
			if (!(parser.get_iss() >> this->capacitance[0]))
			{
				this->capacitance[0] = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for capacitance0.",
					PHRQ_io::OT_CONTINUE);
			}
			capacitance0_defined = true;
			break;

		case 10:				// capacitance1
			if (!(parser.get_iss() >> this->capacitance[1]))
			{
				this->capacitance[1] = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for capacitance1.",
					PHRQ_io::OT_CONTINUE);
			}
			capacitance1_defined = true;
			break;
		case 11:				// sigma0
			if (!(parser.get_iss() >> this->sigma0))
			{
				this->sigma0 = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for sigma0.",
					PHRQ_io::OT_CONTINUE);
			}
			break;
		case 12:				// sigma1
			if (!(parser.get_iss() >> this->sigma1))
			{
				this->sigma1 = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for sigma1.",
					PHRQ_io::OT_CONTINUE);
			}
			break;
		case 13:				// sigma2
			if (!(parser.get_iss() >> this->sigma2))
			{
				this->sigma2 = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for sigma2.",
					PHRQ_io::OT_CONTINUE);
			}
			break;
		case 14:				// sigmaddl
			if (!(parser.get_iss() >> this->sigmaddl))
			{
				this->sigmaddl = 0;
				parser.incr_input_error();
				parser.error_msg("Expected numeric value for sigmaddl.",
					PHRQ_io::OT_CONTINUE);
			}
			break;
		case 15:				// g_map
			{
				if (g_map_first)
				{
					this->g_map.clear();
					g_map_first = false;
				}
				LDBLE z, dummy;
				if (!(parser.get_iss() >> z))
					break;
				cxxSurfDL temp_surf_dl;
				this->g_map[z] = temp_surf_dl;
				std::map<LDBLE, cxxSurfDL>::iterator git = g_map.find(z);
				if (!(parser.get_iss() >> dummy))
					break;
				else
					git->second.Set_g(dummy);
				if (!(parser.get_iss() >> dummy))
					break;
				else
					git->second.Set_dg(dummy);
				if (!(parser.get_iss() >> dummy))
					break;
				else
					git->second.Set_psi_to_z(dummy);
			}
			break;
		case 16:				// dl_species_map
			int s_num;
			if (parser.peek_token() != CParser::TT_EMPTY)
			{
				if (!(parser.get_iss() >> s_num))
				{
					parser.incr_input_error();
					parser.error_msg("Expected integer for species number.",
						PHRQ_io::OT_CONTINUE);
				}
				else
				{
					double d; 
					if (!(parser.get_iss() >> d))
					{
						parser.incr_input_error();
						parser.error_msg("Expected double for species concentration.",
							PHRQ_io::OT_CONTINUE);
					}
					this->dl_species_map[s_num] = d;
				}
			}
			opt_save = 16;



			break;
		}
		if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD)
			break;
	}
	if (check)
	{
		// members that must be defined
		if (specific_area_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Specific_area not defined for SurfaceCharge input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (grams_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Grams not defined for SurfaceCharge input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (charge_balance_defined == false)
		{
			parser.incr_input_error();
			parser.
				error_msg("Charge_balance not defined for SurfaceCharge input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (mass_water_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Mass_water not defined for SurfaceCharge input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (la_psi_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("La_psi not defined for SurfaceCharge input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (capacitance0_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Capacitance0 not defined for SurfaceCharge input.",
				PHRQ_io::OT_CONTINUE);
		}
		if (capacitance1_defined == false)
		{
			parser.incr_input_error();
			parser.error_msg("Capacitance1 not defined for SurfaceCharge input.",
				PHRQ_io::OT_CONTINUE);
		}
	}
}

void
cxxSurfaceCharge::add(const cxxSurfaceCharge & addee, LDBLE extensive)
{
	if (extensive == 0.0)
		return;

	if (this->name.size() == 0 && addee.name.size() == 0)
	{
		return;
	}
	assert(this->name == addee.name);

	LDBLE ext1, ext2, f1, f2;
	ext1 = this->specific_area * this->grams;
	ext2 = addee.specific_area * addee.grams * extensive;
	if (ext1 + ext2 != 0)
	{
		f1 = ext1 / (ext1 + ext2);
		f2 = ext2 / (ext1 + ext2);
	}
	else
	{
		f1 = 0.5;
		f2 = 0.5;
	}
	this->specific_area = f1 * this->specific_area + f2 * addee.specific_area;
	this->grams += addee.grams * extensive;
	this->charge_balance += addee.charge_balance * extensive;
	this->mass_water += addee.mass_water * extensive;
	this->la_psi = this->la_psi * f1 + addee.la_psi * f2;
	this->capacitance[0] =
		this->capacitance[0] * f1 + this->capacitance[0] * f2;
	this->capacitance[1] =
		this->capacitance[1] * f1 + this->capacitance[1] * f2;
	this->diffuse_layer_totals.add_extensive(addee.diffuse_layer_totals, extensive);
}

void
cxxSurfaceCharge::multiply(LDBLE extensive)
{
	this->grams *= extensive;
	this->charge_balance *= extensive;
	this->mass_water *= extensive;
	this->diffuse_layer_totals.multiply(extensive);
}
void
cxxSurfaceCharge::Serialize(Dictionary & dictionary, std::vector < int >&ints, 
	std::vector < double >&doubles)
{

	ints.push_back(dictionary.Find(this->name));
	doubles.push_back(this->specific_area);
	doubles.push_back(this->grams);
	doubles.push_back(this->charge_balance);
	doubles.push_back(this->mass_water);
	doubles.push_back(this->la_psi);
	doubles.push_back(this->capacitance[0]);
	doubles.push_back(this->capacitance[1]);
	this->diffuse_layer_totals.Serialize(dictionary, ints, doubles);
	doubles.push_back(this->sigma0);
	doubles.push_back(this->sigma1);
	doubles.push_back(this->sigma2);
	doubles.push_back(this->sigmaddl);
	ints.push_back((int) this->g_map.size());
	{
		std::map<LDBLE, cxxSurfDL>::iterator it;
		for (it = this->g_map.begin(); it != this->g_map.end(); it++) 
		{
			doubles.push_back(it->first);
			it->second.Serialize(dictionary, ints, doubles);
		}
	}
	ints.push_back((int) this->dl_species_map.size());
	{
		std::map<int, double>::iterator it;
		for (it = this->dl_species_map.begin(); it != this->dl_species_map.end(); it++) 
		{
			ints.push_back(it->first);
			doubles.push_back(it->second);
		}
	}
}

void
cxxSurfaceCharge::Deserialize(Dictionary & dictionary, std::vector < int >&ints, 
	std::vector < double >&doubles, int &ii, int &dd)
{
	this->name = dictionary.GetWords()[ints[ii++]];
	this->specific_area = doubles[dd++];
	this->grams = doubles[dd++];
	this->charge_balance = doubles[dd++];
	this->mass_water = doubles[dd++];
	this->la_psi = doubles[dd++];
	this->capacitance[0] = doubles[dd++];
	this->capacitance[1] = doubles[dd++];
	this->diffuse_layer_totals.Deserialize(dictionary, ints, doubles, ii, dd);
	this->sigma0 = doubles[dd++];
	this->sigma1 = doubles[dd++];
	this->sigma2 = doubles[dd++];
	this->sigmaddl = doubles[dd++];
	{
		this->g_map.clear();
		int count = ints[ii++];
		for (int i = 0; i < count; i++)
		{
			double d = doubles[dd++];
			cxxSurfDL sdl;
			sdl.Deserialize(dictionary, ints, doubles, ii, dd);
			this->g_map[d] = sdl;
		}
	}
	{
		this->dl_species_map.clear();
		int count = ints[ii++];
		for (int i = 0; i < count; i++)
		{
			int j = ints[ii++];
			double d = doubles[dd++];
			dl_species_map[j] = d;
		}
	}


}
void
cxxSurfDL::Serialize(Dictionary & dictionary, std::vector < int >&ints, 
	std::vector < double >&doubles)
{
	doubles.push_back(this->g);
	doubles.push_back(this->dg);
	doubles.push_back(this->psi_to_z);
}

void
cxxSurfDL::Deserialize(Dictionary & dictionary, std::vector < int >&ints, 
	std::vector < double >&doubles, int &ii, int &dd)
{
	this->g = doubles[dd++];
	this->dg = doubles[dd++];
	this->psi_to_z = doubles[dd++];
}
const std::vector< std::string >::value_type temp_vopts[] = {
	std::vector< std::string >::value_type("name"),	                // 0 
	std::vector< std::string >::value_type("specific_area"),	    // 1 
	std::vector< std::string >::value_type("grams"),	            // 2 
	std::vector< std::string >::value_type("charge_balance"),	    // 3 
	std::vector< std::string >::value_type("mass_water"),	        // 4 
	std::vector< std::string >::value_type("la_psi"),	            // 5 
	std::vector< std::string >::value_type("diffuse_layer_totals"),	// 6 
	std::vector< std::string >::value_type("la_psi1"),	            // 7 
	std::vector< std::string >::value_type("la_psi2"),	            // 8 
	std::vector< std::string >::value_type("capacitance0"),	        // 9 
	std::vector< std::string >::value_type("capacitance1"),	        // 10 
	std::vector< std::string >::value_type("sigma0"),	            // 11 
	std::vector< std::string >::value_type("sigma1"),	            // 12 
	std::vector< std::string >::value_type("sigma2"),	            // 13 
	std::vector< std::string >::value_type("sigmaddl"),	            // 14
	std::vector< std::string >::value_type("g_map"),	            // 15
	std::vector< std::string >::value_type("diffuse_layer_species") // 16
};									   
const std::vector< std::string > cxxSurfaceCharge::vopts(temp_vopts, temp_vopts + sizeof temp_vopts / sizeof temp_vopts[0]);	