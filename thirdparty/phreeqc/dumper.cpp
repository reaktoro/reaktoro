#include <algorithm>			// std::replace

#include "dumper.h"
#include "Parser.h"
#include "PHRQ_io.h"

dumper::dumper(PHRQ_io *io)
:
PHRQ_base(io)
{
	this->file_name = "dump.out";
	this->append = false;
	this->on = false;
}
dumper::dumper(CParser & parser, PHRQ_io *io)
:
PHRQ_base(io)
{
	this->file_name = "dump.out";
	this->append = false;
	this->Read(parser);
}

dumper::~dumper(void)
{
}

void dumper::SetAll(bool tf)
{
	this->binList.SetAll(tf);

}
bool dumper::Read(CParser & parser)
{
	bool return_value(true);

	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	int opt_save;

	opt_save = CParser::OPT_DEFAULT;
	bool cleared_once = false;
	this->on = true;

	for (;;)
	{
		int opt;
		StorageBinListItem cells;
		opt = parser.get_option(vopts, next_char);
		if (opt == CParser::OPT_DEFAULT)
		{
			opt = opt_save;
		}
		else
		{
			opt_save = opt;
		}
		if (opt > 1 && !cleared_once)
		{
			binList.SetAll(false);
			cleared_once = true;
		}
		// Select StorageBinListItem
		StorageBinListItem *item = NULL;
		switch (opt)
		{
		case 3:					// cell
		case 4:					// cells
			item = &cells;
			break;
		case 5:
		case 6:
			item = &(this->binList.Get_solution());
			break;
		case 7:
		case 8:
		case 9:
		case 10:
			item = &(this->binList.Get_pp_assemblage());
			break;
		case 11:
			item = &(this->binList.Get_exchange());
			break;
		case 12:
			item = &(this->binList.Get_surface());
			break;
		case 13:
		case 14:
		case 15:
			item = &(this->binList.Get_ss_assemblage());
			break;
		case 16:
		case 17:
			item = &(this->binList.Get_gas_phase());
			break;
		case 18:
			item = &(this->binList.Get_kinetics());
			break;

		case 19:	// mix
			item = &(this->binList.Get_mix());
			break;
		case 20:	// reaction
		case 21:	// reactions
			item = &(this->binList.Get_reaction());
			break;
		case 22:	// temperature
		case 23:	// reaction_temperature
		case 24:    // reaction_temperatures
			item = &(this->binList.Get_temperature());
			break;
		case 25:	// pressure
		case 26:	// reaction_pressure
		case 27:    // reaction_pressures
			item = &(this->binList.Get_pressure());
			break;
		default:
			break;
		}

		// Read dump entity list of numbers or number ranges for line, store in item
		if ((opt > 2))
		{
			for (;;)
			{ 
				CParser::TOKEN_TYPE j = parser.copy_token(token, next_char);
				if (item && j == CParser::TT_DIGIT)
				{
					item->Augment(token);
				}
				else if (item && j == CParser::TT_EMPTY)
				{
					item->Augment(token);
					break;
				}
				else
				{
					parser.error_msg("Expected single number or range of numbers.",
						PHRQ_io::OT_CONTINUE);
				}
			}
		}

		if (opt == 3 || opt == 4)
		{
			this->binList.TransferAll(cells);
		}
		// Process other identifiers
		std::set < int >::iterator it;
		switch (opt)
		{
		case CParser::OPT_EOF:
			break;
		case CParser::OPT_KEYWORD:
			break;
		case 3:
		case 4:
		case 5:
		case 6:
		case 7:
		case 8:
		case 9:
		case 10:
		case 11:
		case 12:
		case 13:
		case 14:
		case 15:
		case 16:
		case 17:
		case 18:
		case 19:
		case 20:
		case 21:
		case 22:
		case 23:
		case 24:
		case 25:
		case 26:
		case 27:
			break;
		case 0:				//file
			std::getline(parser.get_iss(), this->file_name);
			trim(this->file_name);
			if (this->file_name.size() == 0)
			{
				this->file_name = "dump.out";
			}

			break;
		case 1:				//append
			{
				parser.copy_token(token, next_char);
				//if (!(parser.get_iss() >> this->append))
				this->append = true;
				if (token.c_str()[0] == 'f' || token.c_str()[0] == 'F')
				{
					this->append = false;
				}
			}
			break;
		case 2:			//all
			this->SetAll(true);
			break;
		default:
		case CParser::OPT_DEFAULT:
		case CParser::OPT_ERROR:
			opt = CParser::OPT_EOF;
			parser.error_msg("Unknown input reading DUMP definition.",
							 PHRQ_io::OT_CONTINUE);
			parser.error_msg(parser.line().c_str(), PHRQ_io::OT_CONTINUE);
			return_value = false;
			break;
		}
		
		if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD)
			break;
	}
	return(return_value);
}

bool dumper::Get_bool_any(void)
{
	return (
		Get_bool_solution()			||
		Get_bool_pp_assemblage()	||
		Get_bool_exchange()			||
		Get_bool_surface()			||
		Get_bool_ss_assemblage()	||
		Get_bool_gas_phase()		||
		Get_bool_kinetics()			||
		Get_bool_mix()				||
		Get_bool_reaction()			||
		Get_bool_temperature()		||
		Get_bool_pressure()
		);
}
const std::vector< std::string >::value_type temp_vopts[] = {
	std::vector< std::string >::value_type("file"),		                // 0 
	std::vector< std::string >::value_type("append"),		            // 1 
	std::vector< std::string >::value_type("all"),			            // 2 
	std::vector< std::string >::value_type("cell"),			            // 3 
	std::vector< std::string >::value_type("cells"),		            // 4 
	std::vector< std::string >::value_type("solution"),		            // 5 
	std::vector< std::string >::value_type("solutions"),		        // 6 
	std::vector< std::string >::value_type("pp_assemblage"),	        // 7 
	std::vector< std::string >::value_type("pp_assemblages"),	        // 8 
	std::vector< std::string >::value_type("equilibrium_phase"),	    // 9 
	std::vector< std::string >::value_type("equilibrium_phases"),	    // 10
	std::vector< std::string >::value_type("exchange"),		            // 11
	std::vector< std::string >::value_type("surface"),		            // 12
	std::vector< std::string >::value_type("ss_assemblage"),	        // 13
	std::vector< std::string >::value_type("solid_solution"),	        // 14
	std::vector< std::string >::value_type("solid_solutions"),	        // 15
	std::vector< std::string >::value_type("gas_phase"),		        // 16
	std::vector< std::string >::value_type("gas_phases"),		        // 17
	std::vector< std::string >::value_type("kinetics"),		            // 18
	std::vector< std::string >::value_type("mix"),			            // 19
	std::vector< std::string >::value_type("reaction"),		            // 20
	std::vector< std::string >::value_type("reactions"),		        // 21
	std::vector< std::string >::value_type("temperature"),		        // 22
	std::vector< std::string >::value_type("reaction_temperature"),	    // 23
	std::vector< std::string >::value_type("reaction_temperatures"),    // 24
	std::vector< std::string >::value_type("pressure"),		            // 25
	std::vector< std::string >::value_type("reaction_pressure"),	    // 26
	std::vector< std::string >::value_type("reaction_pressures")	    // 27
};
const std::vector< std::string > dumper::vopts(temp_vopts, temp_vopts + sizeof temp_vopts / sizeof temp_vopts[0]);
