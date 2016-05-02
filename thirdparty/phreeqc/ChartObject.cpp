// ChartObject.cpp: implementation of the ChartObject class.
//
//////////////////////////////////////////////////////////////////////
#ifdef MULTICHART
#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include "Utils.h"
#include <iostream>
#include "ChartObject.h"
#include "Parser.h"
#include <fstream>
#include <math.h>
#include <iomanip>
#include "Phreeqc.h"
#include "phqalloc.h"

#include "Form1.h"
using namespace zdg_ui2;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ChartObject::ChartObject(PHRQ_io *io)
	//
	// default constructor for cxxExchComp
	//
	:
cxxNumKeyword(io)
{
	new_ug = false;
	FirstCallToUSER_GRAPH = true;
	
	update_time_chart = 150;
	PanelHeight = 510;
	PanelWidth = 640;

	Symbol_map["Square"] = 0;
	Symbol_map["Diamond"] = 1;
	Symbol_map["Triangle"] = 2;
	Symbol_map["Circle"] = 3;
	Symbol_map["XCross"] = 4;
	Symbol_map["Plus"] = 5;
	Symbol_map["Star"] = 6;
	Symbol_map["TriangleDown"] = 7;
	Symbol_map["HDash"] = 8;
	Symbol_map["VDash"] = 9;
	Symbol_map["None"] = 10;

	Color_vector.push_back("Red");
	Color_vector.push_back("Green");
	Color_vector.push_back("Blue");
	Color_vector.push_back("Orange");
	Color_vector.push_back("Magenta");
//	Color_vector.push_back("Yellow");
	Color_vector.push_back("Black");
	Color_vector.push_back("Cyan");
	Color_vector.push_back("Brown");
	Color_vector.push_back("Lime");
	Color_vector.push_back("Gray");

	chart_title.clear();
	axis_titles.clear();
	
	int i;
	for (i = 0; i < 5; i++)
	{
		axis_scale_x[i] = NA;
		axis_scale_y[i] = NA;
		axis_scale_y2[i] = NA;
	}

	chart_type = 0;
	graph_initial_solutions = false;
	shifts_as_points = false;
	connect_simulations = true;
	colnr = 0;
	ColumnOffset = 0;
	prev_advection_step = 0;
	prev_transport_step = 0;
	prev_sim_no = 0;
	end_timer = false;
	done = false;

	curve_added = false;
	point_added = false;

	user_graph = new rate;
	user_graph->commands = NULL;
	user_graph->name = NULL;
	user_graph->new_def = 0;
	user_graph->linebase = user_graph->loopbase = user_graph->varbase = NULL;

	default_symbol = 0;
	default_symbol_csv = 0;
	default_color = 0;
	default_color_csv = 0;

	graph_x = NA;
	graph_y.clear();
	secondary_y.clear();

	usingResource = 0;
	form_started = false;
	active = true;
	detach = false;

	batch_background = true;
	batch_grid = true;
	batch = ChO_NO_BATCH;
}

ChartObject::~ChartObject()
{
	while (0 != System::Threading::Interlocked::CompareExchange(this->usingResource, 7, 0))
	{
		System::Threading::Thread::Sleep(5);
	}
	this->Rate_free();
	delete this->user_graph;

	std::vector<CurveObject *>::iterator it;
	for (it = this->CurvesCSV.begin(); it != CurvesCSV.end(); it++)
	{
		delete *it;
	}

	for (it = this->Curves.begin(); it != Curves.end(); it++)
	{
		delete *it;
	}
}

bool 
ChartObject::Set_axis_scale(CParser & parser)
{
	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	LDBLE *scale_ptr = NULL;
	std::vector<std::string> string_vector;
	size_t j = 0;

	// rereads option
	parser.copy_token(token, next_char);

	// reads parameters
	while ((j < 6) && (parser.copy_token(token, next_char)) != CParser::TT_EMPTY)
	{
		string_vector.push_back(token);
		j++;
	}
	if (string_vector.size() == 0)
	{
		error_msg("No axis defined for scales", PHRQ_io::OT_CONTINUE);
		return false;
	}
	std::string str = string_vector[0];
	std::string type;
	Utilities::str_tolower(str);
	if (str.substr(0,2) == "y2" || str.substr(0,1) == "s")
	{
		type = "sy";
		scale_ptr = this->axis_scale_y2;
	}
	else if (str.substr(0,1) == "y")
	{
		type = "y";
		scale_ptr = this->axis_scale_y;
	}
	else if (str.substr(0,1) == "x")
	{
		type = "x";
		scale_ptr = this->axis_scale_x;
	}
	else
	{
		std::ostringstream estream;
		estream << "Found " << str;
		estream << ", but expect axis type \'x\', \'y\', \'y2\'or \'sy\'.";
		estream << "\n";
		error_msg(estream.str().c_str(), PHRQ_io::OT_CONTINUE);
		return false;
	}

	for (j = 1; j < string_vector.size() && j < 5; j++)
	{
		std::string s = string_vector[j];
		if (s[0] == 'a' || s[0] == 'A')
		{
			scale_ptr[j - 1] = NA;
		}
		else if (CParser::token_type(s) == CParser::TT_DIGIT)
		{
			scale_ptr[j - 1] = atof(s.c_str());
		}
		else
		{
			std::ostringstream estream;
			estream << "Found " << s;
			estream << ", but expect number or 'a(uto)'.";
			estream << "\n";
			error_msg(estream.str().c_str(), CONTINUE);
			return false;
		}
	}
	if (string_vector.size() == 6)
	{
		std::string s = string_vector[5];
		if (s[0] == 't' || s[0] == 'T' || s[0] == 'l' || s[0] == 'L')
		{
			scale_ptr[4] = 10.0;
			if (((fabs(scale_ptr[0] - NA) > 1e-3) && scale_ptr[0] <=0) ||
				((fabs(scale_ptr[1] - NA) > 1e-3) && scale_ptr[1] <=0))
			{
				std::ostringstream estream;
				estream << "MIN and MAX must be > 0 for log " << type << "-scale.";
				estream << "\n";
				error_msg(estream.str().c_str(), CONTINUE);
				return false;
			}

		}
	}
	if ((fabs(scale_ptr[0] - NA) > 1e-3) && (fabs(scale_ptr[1] - NA) > 1e-3))
	{
		if (scale_ptr[0] > scale_ptr[1])
		{
			std::ostringstream estream;
			estream << "Maximum must be larger than minimum of axis_scale " << type << "-scale." << "\n";
			estream << "Switching values for MIN and MAX. " << "\n";
			warning_msg(estream.str().c_str());
			LDBLE t;
			t = scale_ptr[0];
			scale_ptr[0] = scale_ptr[1];
			scale_ptr[1] = scale_ptr[0];
			return false;
		}
	}
	return true;
}

bool
ChartObject::Read(CParser & parser)
{

	std::istream::pos_type ptr;
	std::istream::pos_type next_char;
	std::string token;
	int opt_save;
	bool useLastLine(false);
	bool new_command_lines(false);
	// if no definitions in USER_GRAPH data block, deactivate plot
	bool no_def(true);
	if (this->FirstCallToUSER_GRAPH)
	{
		this->new_ug = true;
	}
	else
	{
		this->ColumnOffset = (int) this->Curves.size();
		this->new_ug = true;
	}
	this->new_plotxy_curves.clear();

	// Read number and description
	{
		this->read_number_description(parser);
	}

	opt_save = CParser::OPT_DEFAULT;

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
		switch (opt)
		{
		case CParser::OPT_EOF:
			break;
		case CParser::OPT_KEYWORD:
			break;
		case CParser::OPT_ERROR:
			opt = CParser::OPT_EOF;
			error_msg("Unknown input in USER_GRAPH keyword.", CONTINUE);
			error_msg(parser.line().c_str(), CONTINUE);
			useLastLine = false;
			break;
		case 0:				/* start */
			opt_save = CParser::OPT_DEFAULT;
			break;
		case 1:				/* end */
			opt_save = CParser::OPT_DEFAULT;
			break;
		case 2:				/* headings */
		case 3:				/* heading */
			this->new_headings.clear();
			while (parser.copy_token(token, next_char) != CParser::TT_EMPTY)
			{
				this->new_headings.push_back(token);
			}
			this->headings_original = this->new_headings;
			break;
		case 4:				/* chart title */
			{
				std::string tok;
				parser.get_rest_of_line(tok);
				std::string::iterator b = tok.begin();
				std::string::iterator e = tok.end();
				CParser::copy_title(this->chart_title, b, e);
			}
			break;
		case 5:	/* axis titles */
			{
				this->axis_titles.clear();
				std::string l; 
				parser.get_rest_of_line(l);
				std::string tok;
				std::string::iterator b = l.begin();
				std::string::iterator e = l.end();
				this->axis_titles.clear();
				while(parser.copy_title(tok,  b, e) != CParser::TT_EMPTY)
				{
					this->axis_titles.push_back(tok);
				}
			}
			while (this->axis_titles.size() < 3) this->axis_titles.push_back("");
			break;
		case 6:	/* axis scales */
			{ 
				this->Set_axis_scale(parser);
			}
			break;
		case 7: /* initial_solutions */
			this->graph_initial_solutions = parser.get_true_false(next_char, FALSE);
			break;
		case 8: /* plot_concentration_vs */
			parser.copy_token(token, next_char);
			Utilities::str_tolower(token);
			if (token[0] == 'x' || token[0] == 'd')
				chart_type = 0;
			else if (token[0] == 't')
				chart_type = 1;
			else
			{
				std::ostringstream estream;
				estream << "Found " << token << ", but expect plot type: (\'x\' or \'dist\') for distance, (\'t\') for time.";
				error_msg(estream.str().c_str(), CONTINUE);			
			}
			break;
		case 9:  /* shifts_as_points */
			this->shifts_as_points = parser.get_true_false(next_char, true);
			if (this->shifts_as_points)
				this->chart_type = 0;
			else
				this->chart_type = 1;
			break;
		case 10: /* grid_offset */
#ifdef PHREEQ98
			/*
			i = copy_token(token, &next_char, &l);
			str_tolower(token);
			if (i == DIGIT)
				sscanf(token, "%d", &RowOffset);
			i = copy_token(token, &next_char, &l);
			str_tolower(token);
			if (i == DIGIT)
				sscanf(token, "%d", &ColumnOffset);
			*/
#endif
			break;
		case 11: /* connect_simulations */
			this->connect_simulations = parser.get_true_false(next_char, true);
			break;
		case 12: /* plot_csv_file */
		case 17: /* plot_tsv_file */
			{
				std::string file_name;
				parser.get_rest_of_line(file_name);
				file_name = trim(file_name);
				this->OpenCSVFile(file_name);
				this->csv_file_names.push_back(file_name);
			}
			break;
		case 13: /* clear */
		case 14: /* detach */
			this->detach = true;
			break;
			/* End of modifications */
		case 15: /* active */
			this->active = parser.get_true_false(next_char, true);
			if (this->active && this->Get_connect_simulations())
			{
				this->ColumnOffset = 0;//this->Curves.size();
				this->new_ug = false;
			}
			break;
			/* End of modifications */
		case 16: /* batch */
			{
				this->batch = ChartObject::ChO_BATCH_ONLY;
				std::string rest_of_line, lc_rest_of_line;
				parser.get_rest_of_line(rest_of_line);
				lc_rest_of_line = rest_of_line;
				Utilities::str_tolower(lc_rest_of_line);
	
				std::vector<std::string> file_types;
				file_types.push_back(".emf");
				file_types.push_back(".png");
				file_types.push_back(".jpg");
				file_types.push_back(".gif");
				file_types.push_back(".tiff");
				file_types.push_back(".bmp");
				file_types.push_back(".jpeg");

				size_t first, last;
				size_t i;
				for (i = 0; i < file_types.size(); i++)
				{
					first = lc_rest_of_line.rfind(file_types[i].c_str());
					if (first != std::string::npos)
					{
						break;
					}
				}
				if (i >=  file_types.size())
				{
					std::ostringstream estream;
					estream << "Batch file name must have suffix emf, png, jpg, jpeg, gif, tiff, or bmp.";
					error_msg(estream.str().c_str(), CONTINUE);	
					break;
				}
				switch (i)
				{
				case 0:
					this->batch = ChartObject::ChO_EMF;
					last = first + 4;
					break;
				case 1:
					this->batch = ChartObject::ChO_PNG;
					last = first + 4;
					break;
				case 2:
					this->batch = ChartObject::ChO_JPG;
					last = first + 4;
					break;
				case 3:
					this->batch = ChartObject::ChO_GIF;
					last = first + 4;
					break;
				case 4:
					this->batch = ChartObject::ChO_TIFF;
					last = first + 5;
					break;
				case 5:
					this->batch = ChartObject::ChO_BMP;
					last = first + 4;
					break;
				case 6:
					this->batch = ChartObject::ChO_JPG;
					last = first + 5;
					break;
				}

				this->batch_fn = rest_of_line.substr(0, last);
				if (last+1 < rest_of_line.size())
				{
					token = rest_of_line.substr(last);
					token = trim(token);
					std::string tf;
					std::string::iterator b = token.begin();
					std::string::iterator e = token.end();

					CParser::copy_token(tf, b, e);
					if (tf.size() > 0)
					{
						Utilities::str_tolower(tf);
						if (tf[0] == 'f') 
						{
							this->batch_background = false;
						}
					}
					CParser::copy_token(tf, b, e);
					if (tf.size() > 0)
					{
						Utilities::str_tolower(tf);
						if (tf[0] == 'f') 
						{
							this->batch_grid = false;
						}
					}
				}
			}
			break;
		case CParser::OPT_DEFAULT:	// Read Basic commands
			{
				if (!new_command_lines)
				{
					this->rate_command_list.clear();
					this->rate_command_list_original.clear();
					new_command_lines = true;
				}
				this->rate_new_def = true;
				/* read command */
				std::string cmd(parser.line());
				this->rate_command_list_original.push_back(cmd);
				std::string cmd_lower = cmd;
				Utilities::str_tolower(cmd_lower);
				if ((cmd_lower.find("graph_y") != std::string::npos) || 
					(cmd_lower.find("graph_sy") != std::string::npos))
				{
					//Number of curves not known here
					//Curves are created in Basic cmds
				}
				if (cmd_lower.find("plot_xy") != std::string::npos)
				{
					//Curves are created in Basic cmds
					CurveObject new_curve = ExtractCurveInfo(cmd); // truncates cmd
					// Add to list of new plotxy curves
					this->new_plotxy_curves.push_back(new_curve);
				}
				this->rate_command_list.push_back(cmd);
			}
			opt_save = CParser::OPT_DEFAULT;
			break;
		}
		if (opt == CParser::OPT_EOF || opt == CParser::OPT_KEYWORD)
			break;
		no_def = false;
	}

	// install new plotxy commands
	// disable this user_graph if USER_GRAPH block is empty
	//if (new_command_lines || no_def) this->Set_rate_struct();
	if (new_command_lines) this->Set_rate_struct();
	if (no_def) this->Set_active(false);
	return true;
}
bool
ChartObject::OpenCSVFile(std::string file_name)
{
	std::string token;

	std::ifstream f_csv(file_name.c_str(), std::ifstream::in);

	if (!f_csv.is_open())
	{
		std::ostringstream estream;
		estream << "Could not open tsv file for USER_GRAPH " << file_name << "\n       Please, give the full path + filename.";
		error_msg(estream.str().c_str(), CONTINUE);
		return false;
	}

	CParser parser(f_csv, this->Get_io());
	parser.set_echo_file(CParser::EO_NONE);	

	/* Get lines */
	int linenr = 0; 

	// temporary storage for new CSV curves
	std::vector<std::string> headings;
	std::vector<CurveObject *> csv_curves;

	// Read file line by line into temporary curves
	while (parser.check_line("cvs file", false, true, true, false) != PHRQ_io::LT_EOF)
	{
		// Headings line
		if (linenr == 0)
		{
			// Skip x in 1st column
			parser.get_iss() >> token;

			// Read rest of headings
			while (parser.get_iss() >> token)
			{
				headings.push_back(token);
			}
			// add curves to temporary csv_curves
			size_t i; 
			for (i = 0; i < headings.size(); i++)
			{
				CurveObject *c = new CurveObject;
				c->Set_id(headings[i]);
				c->Set_line_w(0);
				std::string sym = "";
				this->Get_legal_symbol_csv(sym);
				c->Set_symbol(sym);
				csv_curves.push_back(c);
			}

			linenr++;
			continue;
		}

		token = parser.line();
		std::string tok1;
		CParser::TOKEN_TYPE tt = CParser::parse_delimited(token, tok1, "\t");

		// Curve properties lines
		if (linenr < 6 && tt != CParser::TT_DIGIT)
		{
			Utilities::str_tolower(tok1);
			std::string tok2;
			size_t i=0;
			while (token.size() != 0 && i < csv_curves.size())
			{
				CParser::parse_delimited(token, tok2, "\t");
				tok2 = trim(tok2);
				if (!strncmp(tok1.c_str(), "color", 5))
				{
					csv_curves[i]->Set_color(tok2);
				}
				else if (!strncmp(tok1.c_str(), "symbol", 5) && (strstr(tok1.c_str(), "_si") == NULL)
					&& (strstr(tok1.c_str(), "-si") == NULL))
				{
					csv_curves[i]->Set_symbol(tok2);
				}
				else if (!strncmp(tok1.c_str(), "symbol_size", 8) || !strncmp(tok1.c_str(), "symbol-size", 8))
				{
					if (tok2.size() > 0)
					{
						csv_curves[i]->Set_symbol_size(atof(tok2.c_str()));
					}
				}
				else if (!strncmp(tok1.c_str(), "line_w", 5) || !strncmp(tok1.c_str(), "line-w", 5))
				{
					csv_curves[i]->Set_line_w(atof(tok2.c_str()));
				}
				else if (!strncmp(tok1.c_str(), "y_axis", 5) || !strncmp(tok1.c_str(), "y-axis", 5))
				{
					if (tok2.size() > 0)
					{
						csv_curves[i]->Set_y_axis(atoi(tok2.c_str()));
					}
				}
				i++;
			}
			linenr++;
			continue;
		}

		// Curve data
		if (linenr < 6)	linenr = 6;
		if (tt != CParser::TT_DIGIT) 
		{
			linenr++;	
			continue;
		} 

		// x value for all curves
		LDBLE x_value = atof(tok1.c_str());

		// y values for curves
		std::string tok2;
		size_t i=0;
		while (token.size() != 0 && i < csv_curves.size())
		{
			CParser::parse_delimited(token, tok2, "\t");
			Utilities::squeeze_white(tok2);
			if (tok2.size() == 0)
			{
				//csv_curves[i].Get_x().push_back(NA);
				//csv_curves[i].Get_y().push_back(NA);
			}
			else
			{
				csv_curves[i]->Get_x().push_back(x_value);
				csv_curves[i]->Get_y().push_back(atof(tok2.c_str()));
			}
			i++;
		}
		linenr++;
	}

	// Append new curves

	std::vector<CurveObject *>::iterator it = csv_curves.begin();
	for (; it != csv_curves.end(); it++)
	{
		if ((*it)->Get_x().size() > 0)
		{
			std::string col = (*it)->Get_color();
			this->Get_color_string_csv(col);
			(*it)->Set_color(col);
			this->CurvesCSV.push_back(*it);
			this->Set_curve_added(true);
		}
	}
	return true;
}

CurveObject 
ChartObject::ExtractCurveInfo(std::string & cmd_line)
{
	/* plot_xy x, tot("Cl"), color = Red, symbol = Circle, symbol_size = 0.0, line_w = 1.0, y_axis = 2 */

	// Make copy of cmd_line
	//int curvenr = (int) this->Curves.size();
	std::string str_line = cmd_line;

	// Massage line
	while (Utilities::replace(" ,",",",str_line));
	while (Utilities::replace("\t,",",",str_line));
	while (Utilities::replace(",","#",str_line));
	while (Utilities::replace("#",", ", str_line));

	int sel;
	std::string token, tok1, tok2;
	std::string revised_line;
	std::string::iterator b = str_line.begin();
	std::string::iterator e = str_line.end();
	// new curve
	CurveObject new_curve;

	while (CParser::copy_token(token, b, e) != CParser::TT_EMPTY)
	{
		sel = -1;
		std::string token_save = token;
		Utilities::str_tolower(token);
		tok1 = token;
		if (!strncmp(tok1.c_str(), "color", 5))
		{
			sel = 0;
		}
		else if (!strncmp(tok1.c_str(), "symbol", 5) && (strstr(tok1.c_str(), "_si") == NULL)
			&& (strstr(tok1.c_str(), "-si") == NULL))
		{
			sel = 1;
		}
		else if (!strncmp(tok1.c_str(), "symbol_size", 8) || !strncmp(tok1.c_str(), "symbol-size", 8))
		{
			sel = 2;
		}
		else if (!strncmp(tok1.c_str(), "line_w", 5) || !strncmp(tok1.c_str(), "line-w", 5))
		{
			sel = 3;
		}
		else if (!strncmp(tok1.c_str(), "y_axis", 5) || !strncmp(tok1.c_str(), "y-axis", 5))
		{
			sel = 4;
		}

		// Append to command line
		if (sel < 0)
		{
			revised_line.append(token_save);
			revised_line.append(" ");
		}
		// Parse plot_xy pair
		else
		{
			// remove preceding comma
			std::string comma = revised_line.substr(revised_line.size() - 2);
			if (comma == ", ")
			{
				revised_line = revised_line.substr(0, revised_line.size() - 2);
				revised_line.append(" ");
			}

			token = token_save;
			size_t p1 = token.find("=");
			std::string tok2;
			// "=" found
			if(p1 != std::string::npos)
			{
				if (p1 != token.size() - 1)
				{
					// color=Red
					tok2 = token.substr(p1 + 1);
				}
				else
				{
					// color= Red
					CParser::copy_token(tok2, b, e);
				}
			}
			else
			{
				// no "=" found
				CParser::copy_token(tok2, b, e);
				p1 = tok2.find("=");
				if (tok2 == "=")
				{
					// color = Red
					CParser::copy_token(tok2, b, e);
				}
				else if (p1 != std::string::npos)
				{
					// color =Red
					tok2 = tok2.substr(p1 + 1);
				}
				else
				{
					// color Red
					tok2 = tok2;
				}
			}
			// remove any commas
			while(Utilities::replace(",","",tok1));
			while(Utilities::replace(",","",tok2));

			// tok1 is name, tok2 is value

			switch (sel)
			{
			case 0:	
				new_curve.Set_color(tok2);
				break;
			case 1:	
				new_curve.Set_symbol(tok2);
				break;
			case 2:
				new_curve.Set_symbol_size(atof(tok2.c_str()));
				break;
			case 3:
				new_curve.Set_line_w(atof(tok2.c_str()));
				break;
			case 4:
				new_curve.Set_y_axis(atoi(tok2.c_str()));
				break;
			}
		}
	}
	
	cmd_line = revised_line;
	return new_curve;
}
void 
ChartObject::Set_rate_struct(void)
{

	if (rate_command_list.size() == 0) return;
	std::list<std::string>::iterator it = rate_command_list.begin();
	std::ostringstream oss;
	for (; it != rate_command_list.end(); it++)
	{
		oss << *it << "\n";
	}
	this->Rate_free();
	if (this->phreeqc_ptr)
	{
		this->user_graph->commands = (char *) phreeqc_ptr-> PHRQ_malloc((oss.str().size()) + 100 * sizeof(char));
	}
	::strcpy(this->user_graph->commands, oss.str().c_str());
	this->user_graph->new_def = this->rate_new_def;
	this->user_graph->loopbase = NULL;
	this->user_graph->varbase = NULL;
	this->user_graph->linebase = NULL;
	this->user_graph->name = NULL;
}
void
ChartObject::Get_legal_symbol(std::string &sym)
{
	std::map<std::string, int>::iterator it;
	if ((it = this->Symbol_map.find(sym)) == this->Symbol_map.end())
	{
		sym = "Default";
		for (it = this->Symbol_map.begin(); it != this->Symbol_map.end(); it++)
		{
			if (default_symbol == it->second)
			{
				sym = it->first;
				break;
			}
		}
		default_symbol++;
		default_symbol = default_symbol % this->Symbol_map.size();
	}
	return;
}
void
ChartObject::Get_legal_symbol_csv(std::string &sym)
{
	std::map<std::string, int>::iterator it;
	if ((it = this->Symbol_map.find(sym)) == this->Symbol_map.end())
	{
		sym = "Default";
		for (it = this->Symbol_map.begin(); it != this->Symbol_map.end(); it++)
		{
			if (this->default_symbol_csv == it->second)
			{
				sym = it->first;
				break;
			}
		}
		default_symbol_csv++;
		default_symbol_csv = default_symbol_csv % this->Symbol_map.size();
	}
	return;
}
void
ChartObject::Get_color_string(std::string &color)
{
	if (color.size() == 0)
	{
		color = Color_vector[this->default_color];
		default_color++;
		default_color = default_color % this->Color_vector.size();
	}
	return;
}
void
ChartObject::Get_color_string_csv(std::string &color)
{
	if (color.size() == 0)
	{
		color = Color_vector[this->default_color_csv];
		default_color_csv++;
		default_color_csv = default_color_csv % this->Color_vector.size();
	}
	return;
}
ZedGraph::SymbolType 
ChartObject::Return_SymbolType(const std::string sym)
{
	int i;
	std::map<std::string, int>::iterator it;
	if ((it = this->Symbol_map.find(sym)) != this->Symbol_map.end())
	{
		i = it->second;
	}
	else
	{
		i = default_symbol++;
		default_symbol = default_symbol % this->Symbol_map.size();
	}

	switch (i)
	{
		case 0:
			return SymbolType::Square;
			break;
		case 1:
			return SymbolType::Diamond;
			break;
		case 2:
			return SymbolType::Triangle;
			break;
		case 3:
			return SymbolType::Circle;
			break;
		case 4:
			return SymbolType::XCross;
			break;
		case 5:
			return SymbolType::Plus;
			break;
		case 6:
			return SymbolType::Star;
			break;
		case 7:
			return SymbolType::TriangleDown;
			break;
		case 8:
			return SymbolType::HDash;
			break;
		case 9:
			return SymbolType::VDash;
			break;
		case 10:
			return SymbolType::None;
			break;
		default:
			return SymbolType::Default;
			break;
	}
}
void 
ChartObject::SaveCurvesToFile(std::string &file_name)
{
	// reimplemented in Form
	// This version not currently used
	std::ofstream f_out(file_name.c_str(), std::ifstream::out);

	if (!f_out.is_open())
	{
		std::ostringstream estream;
		estream << "Could not open csv file for USER_GRAPH " << file_name;
		error_msg(estream.str().c_str(), CONTINUE);
		return;
	}

	// list of curves
	std::vector<CurveObject *> all_curves; 
	size_t i;
	for (i = 0; i < this->CurvesCSV.size(); i++)
	{
		all_curves.push_back(this->CurvesCSV[i]);
	}
	for (i = 0; i < this->Curves.size(); i++)
	{
		all_curves.push_back(Curves[i]);
	}
	// write headings
	size_t max_points = 0; 
	std::vector<CurveObject *>::iterator it = all_curves.begin();
	f_out.precision(4);
	i = 0;
	for ( ; it != all_curves.end(); it++)
	{
		f_out.width(12);
		f_out << "x" << "\t";
		f_out.width(12);
		if ((*it)->Get_id().size() > 0) 
		{
			f_out << (*it)->Get_id() << "\t";
		}
		else
		{
			f_out << "y" << "\t";;
		}
		if ((*it)->Get_x().size() > max_points)
			max_points = (*it)->Get_x().size();
	}
	f_out << "\n";

	// write data
	size_t i2 = 0;
	f_out << std::scientific;
	f_out.precision(4);
	while (i2 < max_points)
	{
		for (it = all_curves.begin(); it != all_curves.end(); it++)
		{
			if (i2 < (*it)->Get_x().size())
			{
				f_out.width(12);
				f_out << (*it)->Get_x()[i2] << "\t";
				f_out.width(12);
				f_out << (*it)->Get_y()[i2] << "\t";
			}
			else if (i2 < max_points)
			{
				f_out.width(13);
				f_out << "\t";
				f_out.width(13);
				f_out << "\t";
			}
		}
		f_out << "\n";
		i2++;
	}
	f_out.close();
	return;
}
// file only used with MULTICHART
bool
ChartObject::start_chart(void)
{
	Application::EnableVisualStyles();
	Application::SetCompatibleTextRenderingDefault(true); 

	// needed to send ChartObject pointer to thread
	Thread ^t = gcnew Thread(
                gcnew ParameterizedThreadStart(Form1::ThreadForm));

	t->SetApartmentState(ApartmentState::STA);
	t->IsBackground = false;
	t->Priority = ThreadPriority::Normal;

	ChartObj ^p = gcnew ChartObj(this);
	t->Start(p);
	this->form_started = true;

	//Thread::Sleep( 1 ); /* this when debugging... */
	//_beginthread(void (Form1::ThreadForm), 0, NULL);
	return true;
}

void 
ChartObject::Rate_free(void)
{

	
	if (this->phreeqc_ptr)
	{
		user_graph->commands = (char *) phreeqc_ptr-> free_check_null(user_graph->commands);
	}
	if (user_graph->linebase != NULL)
	{
		char cmd[] = "new; quit";
		if (this->phreeqc_ptr)
		{
			phreeqc_ptr-> basic_run(cmd, user_graph->linebase, user_graph->varbase, user_graph->loopbase);
		}
		user_graph->linebase = NULL;
		user_graph->varbase = NULL;
		user_graph->loopbase = NULL;
	}
}
void ChartObject::Initialize_graph_pts(void)
{
	graph_x = NA;
	graph_y.clear();
	secondary_y.clear();
}

void 
ChartObject::Finalize_graph_pts(void)
{
	if (graph_x != NA)
	{
		std::map<int, LDBLE>::iterator it;
		for (it = graph_y.begin(); it != graph_y.end(); it++)
		{
			Curves[it->first]->Get_x().push_back(graph_x);
			Curves[it->first]->Get_y().push_back(it->second);
			// Mark added curve for first point, might have been invisible in DefineCurves
			if (Curves[it->first]->Get_x().size() == 1)
				this->Set_curve_added(true);
		}
	}
	if (graph_x != NA)
	{
		std::map<int, bool>::iterator it;
		for (it  =  secondary_y.begin(); it != secondary_y.end(); it++)
		{
			Curves[it->first]->Set_y_axis(2);
		}
	}
}
void 
ChartObject::Add_new_series(void)
{
	std::vector<CurveObject *> Curves;
	size_t i;
	size_t e = this->Curves.size();
	for (i = this->ColumnOffset; i < e; i++)
	{
		this->Add_curve(false, 
			this->Curves[i]->Get_id(),
			this->Curves[i]->Get_line_w(),
			"",
			this->Curves[i]->Get_symbol_size(),
			this->Curves[i]->Get_y_axis(),
			"");
	}
	this->ColumnOffset = (int) e;
	this->curve_added = true;
}
void 
ChartObject::Add_curve(bool plotxy, std::string id, 
					   LDBLE line_width, 
					   std::string symbol,
					   LDBLE symbol_size, 
					   int y_axis,
					   std::string color)				
{
	CurveObject *c = new CurveObject;
	c->Set_id(id);
	c->Set_line_w(line_width);
	this->Get_legal_symbol(symbol);
	c->Set_symbol(symbol);
	c->Set_symbol_size(symbol_size);
	if (this->CurvesCSV.size() > this->Curves.size() && !plotxy)
	{
		c->Set_symbol_size(0.0);
	}
	c->Set_y_axis(y_axis);
	this->Get_color_string(color);
	c->Set_color(color);

	this->Curves.push_back(c);
}
void 
ChartObject::Set_rate_new_def(bool tf)
{
	this->rate_new_def = tf;
	if (this->user_graph != NULL)
	{
		if (tf)
		{
			this->user_graph->new_def = 1;
		}
		else
		{
			this->user_graph->new_def = 0;
		}
	}	
}
void 
ChartObject::dump(std::ostream & oss, unsigned int indent)
{
	size_t i;
	oss.precision(DBL_DIG - 1);
	std::string indent0(""), indent1("");
	for (i = 0; i < indent; ++i)
		indent0.append(Utilities::INDENT);
	for (i = 0; i < indent + 1; ++i)
		indent1.append(Utilities::INDENT);
	oss << indent0 << "USER_GRAPH " << this->n_user << " " << this->description << "\n";

	// chart title
	oss << indent1 << "-chart_title \"" << this->chart_title << "\"" << "\n";

	// axis titles
	if (this->axis_titles.size() > 0)
	{
		oss << indent1 << "-axis_titles "; 

		for (i = 0; i < this->axis_titles.size(); i++)
		{
			oss <<  "\"" << axis_titles[i] << "\" ";
		}
		oss << "\n";
	}

	// axis_scale_x
	LDBLE *scale_ptr = this->axis_scale_x;
	{
		oss << indent1 << "-axis_scale x_axis ";
		for (i = 0; i < 4; i++)
		{
			if (scale_ptr[i] == NA)
			{
				oss << " auto";
			}
			else
			{
				oss << " " << scale_ptr[i];
			}
		}
		if (scale_ptr[4] == 10.0)
		{
			oss << " log";
		}
		oss << "\n";
	}
	// axis_scale_y
	scale_ptr = this->axis_scale_y;
	{
		oss << indent1 << "-axis_scale y_axis ";
		for (i = 0; i < 4; i++)
		{
			if (scale_ptr[i] == NA)
			{
				oss << " auto";
			}
			else
			{
				oss << " " << scale_ptr[i];
			}
		}
		if (scale_ptr[4] == 10.0)
		{
			oss << " log";
		}
		oss << "\n";
	}
	// axis_scale_sy
	scale_ptr = this->axis_scale_y2;
	{
		oss << indent1 << "-axis_scale sy_axis ";
		for (i = 0; i < 4; i++)
		{
			if (scale_ptr[i] == NA)
			{
				oss << " auto";
			}
			else
			{
				oss << " " << scale_ptr[i];
			}
		}
		if (scale_ptr[4] == 10.0)
		{
			oss << " log";
		}
		oss << "\n";
	}
	// chart type
	if (this->chart_type == 0)
	{
		oss << indent1 << "-plot_concentration_vs x" << "\n";
	}
	else
	{
		oss << indent1 << "-plot_concentration_vs t" << "\n";
	}
	// graph_initial_solutions
	if (this->graph_initial_solutions)
	{
		oss << indent1 << "-initial_solutions true" << "\n";
	}
	else
	{
		oss << indent1 << "-initial_solutions false" << "\n";
	}
	// connect_simulations
	if (this->connect_simulations)
	{
		oss << indent1 << "-connect_simulations true" << "\n";
	}
	else
	{
		oss << indent1 << "-connect_simulations false" << "\n";
	}
	// csv files
	for (i = 0; i < this->csv_file_names.size(); i++)
	{
		oss << indent1 << "-plot_tsv_file " << this->csv_file_names[i] << "\n";
	}

	// headings
	if (this->headings_original.size() > 0)
	{
		oss << indent1 << "-headings "; 
		for (i = 0; i < this->headings_original.size(); i++)
		{
			oss << this->headings_original[i] << " ";
		}
		oss << "\n";
	}

	// commands
	oss << indent1 << "-start" << "\n";
	std::list<std::string>::iterator it = rate_command_list_original.begin();
	for (; it != rate_command_list_original.end(); it++)
	{
		oss << *it << "\n";
	}
	oss << indent1 << "-end" << "\n";

	/*
	struct rate *user_graph;
	// C++ for rate struct
	std::string rate_name;
	std::list<std::string> rate_command_list;
	bool rate_new_def;

	int default_symbol;
	int default_symbol_csv;
	int default_color;
	int default_color_csv;

	// temporary storage before stored graph_x/y/sy data are stored in curves
	// Initialize_graph_pts and Finalize_graph_pts use this storage.
	LDBLE graph_x;
	std::map<int, LDBLE> graph_y;
	std::map<int, bool> secondary_y;

	// temporary plotxy curve definitions before stored in curves
	// a plotxy curve is copied to Curves when cmdplotxy is encountered
	// this keeps order correct between plotxy and graph_x/y/sy
	std::vector<CurveObject> new_plotxy_curves;

	// temporary headings until stored during basic_run
	std::vector<std::string> new_headings;
	bool active;
	bool detach;
	bool form_started;
	*/
}

const std::vector< std::string >::value_type temp_vopts[] = {
	std::vector< std::string >::value_type("start"),                   // 0 
	std::vector< std::string >::value_type("end"),	                   // 1 
	std::vector< std::string >::value_type("heading"),	               // 2 
	std::vector< std::string >::value_type("headings"),                // 3 
	std::vector< std::string >::value_type("chart_title"),	           // 4 
	std::vector< std::string >::value_type("axis_titles"),	           // 5 
	std::vector< std::string >::value_type("axis_scale"),	           // 6 
	std::vector< std::string >::value_type("initial_solutions"),       // 7 
	std::vector< std::string >::value_type("plot_concentration_vs"),   // 8 
	std::vector< std::string >::value_type("shifts_as_points"),        // 9 
	std::vector< std::string >::value_type("grid_offset"),	           // 10
	std::vector< std::string >::value_type("connect_simulations"),	   // 11
	std::vector< std::string >::value_type("plot_csv_file"),	       // 12	
	std::vector< std::string >::value_type("clear"),	               // 13
	std::vector< std::string >::value_type("detach"),	               // 14
	std::vector< std::string >::value_type("active"),	               // 15
	std::vector< std::string >::value_type("batch"),                   // 16
	std::vector< std::string >::value_type("plot_tsv_file")            // 17
};
const std::vector< std::string > ChartObject::vopts(temp_vopts, temp_vopts + sizeof temp_vopts / sizeof temp_vopts[0]);


#endif // MULTICHART