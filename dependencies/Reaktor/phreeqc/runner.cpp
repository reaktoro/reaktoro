#include "runner.h"
#include "Parser.h"
#include "NA.h"
#include "Utils.h"
runner::runner(PHRQ_io *io)
:
PHRQ_base(io)
{
	this->time_step = NA;
	this->start_time = NA;
	this->run_cells = false;

}
runner::runner(CParser & parser, PHRQ_io *io)
:
PHRQ_base(io)
{
	this->time_step = NA;
	this->start_time = NA;
	this->run_cells = false;
	this->Read(parser);
}

runner::~runner(void)
{
}
bool runner::Read(CParser & parser)
{

	bool return_value(true);

	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	int opt_save;

	this->Get_cells().Set_defined(true);
	opt_save = CParser::OPT_DEFAULT;

	StorageBinListItem item;
	for (;;)
	{
		int opt;
		opt = parser.get_option(vopts, next_char);
		if (opt == CParser::OPT_DEFAULT)
		{
			opt = opt_save;
		}
		else
		{
			opt_save = opt;
		}

		// Process other identifiers
		std::set < int >::iterator it;
		switch (opt)
		{
		case CParser::OPT_EOF:
			break;
		case CParser::OPT_KEYWORD:
			break;

		case 0:
		case 1:
			for (;;)
			{ 
				CParser::TOKEN_TYPE j = parser.copy_token(token, next_char);
				if (j == CParser::TT_DIGIT)
				{
					item.Augment(token);
				}
				else if (j == CParser::TT_EMPTY)
				{
					item.Augment(token);
					break;
				}
				else
				{
					parser.error_msg("Expected single number or range of numbers.",
						PHRQ_io::OT_CONTINUE);
				}
			}			
			break;
		case 2:				//start_time
			if (!(parser.get_iss() >> this->start_time))
			{
				parser.error_msg("Expected start_time for RUN_CELLS.",
					PHRQ_io::OT_CONTINUE);
				parser.error_msg(parser.line().c_str(), PHRQ_io::OT_CONTINUE);
				break;
			}
			{
				std::string token;
				if (parser.get_iss() >> token)
				{
					token = trim(token);
					this->start_time = Utilities::convert_time(this->start_time, token, "s");
				}
			}
			break;
		case 3:				//time_step
		case 4:				//time_steps
		case 5:				//step
		case 6:				//steps
			if (!(parser.get_iss() >> this->time_step))
			{
				parser.error_msg("Expected time_step for RUN_CELLS.",
					PHRQ_io::OT_CONTINUE);
				parser.error_msg(parser.line().c_str(), PHRQ_io::OT_CONTINUE);
				break;
			}
			{
				std::string token;
				if (parser.get_iss() >> token)
				{
					token = trim(token);
					this->time_step = Utilities::convert_time(this->time_step, token, "s");
				}
			}
			break;
		default:
		case CParser::OPT_DEFAULT:
		case CParser::OPT_ERROR:
			opt = CParser::OPT_EOF;
			parser.error_msg("Unknown input reading RUN_CELLS definition.",
							 PHRQ_io::OT_CONTINUE);
			parser.error_msg(parser.line().c_str(), PHRQ_io::OT_CONTINUE);
			return_value = false;
			break;
		}
		if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD)
			break;
	}
	if (item.Get_numbers().size() > 0)
	{
		this->cells = item;
	}
	return(return_value);
}
const std::vector< std::string >::value_type temp_vopts[] = {
		std::vector< std::string >::value_type("cell"),		 	 // 0 
		std::vector< std::string >::value_type("cells"),		 // 1 
		std::vector< std::string >::value_type("start_time"),	 // 2 
		std::vector< std::string >::value_type("time_step"),	 // 3 
		std::vector< std::string >::value_type("time_steps"),	 // 4 
		std::vector< std::string >::value_type("step"),		 	 // 5 
		std::vector< std::string >::value_type("steps")			 // 6 
};									   
const std::vector< std::string > runner::vopts(temp_vopts, temp_vopts + sizeof temp_vopts / sizeof temp_vopts[0]);	
