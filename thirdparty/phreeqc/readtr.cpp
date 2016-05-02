#include <iostream>				/* std::cout std::cerr */
#include <sstream>
#include <fstream>
#include "StorageBin.h"
#include "SS.h"

typedef unsigned char boolean;
#include "Phreeqc.h"
#include "phqalloc.h"
#include "Utils.h"


#define OPTION_EOF -1
#define OPTION_KEYWORD -2
#define OPTION_ERROR -3
#define OPTION_DEFAULT -4
#define OPTION_DEFAULT2 -5

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_transport(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads advection and column information
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
	char *ptr;
	int i, j, l;
	int count_length, count_disp, count_punch, count_print, count_por;
	int count_length_alloc, count_disp_alloc, count_por_alloc;
	char token[MAX_LENGTH];
	char *description;
	int n_user, n_user_end;
	LDBLE *length, *disp, *pors;
	int *punch_temp, *print_temp;
	int return_value, opt, opt_save;
	char *next_char, *next_char_save;
	char file_name[MAX_LENGTH];

	const char *opt_list[] = {
		"cells",				/* 0 */
		"shifts",				/* 1 */
		"print",				/* 2 */
		"selected_output",		/* 3 */
		"bcond",				/* 4 */
		"timest",				/* 5 */
		"diffc",				/* 6 */
		"tempr",				/* 7 */
		"length",				/* 8 */
		"disp",					/* 9 */
		"punch",				/* 10 */
		"stagnant",				/* 11 */
		"bc",					/* 12 */
		"boundary_conditions",	/* 13 */
		"time_step",			/* 14 */
		"temp_retardation_factor",	/* 15 */
		"diffusion_coefficient",	/* 16 */
		"dispersivity",			/* 17 */
		"direction",			/* 18 */
		"temperature_retardation_factor",	/* 19 */
		"print_cells",			/* 20 */
		"selected_cells",		/* 21 */
		"flow_direction",		/* 22 */
		"flow",					/* 23 */
		"lengths",				/* 24 */
		"dispersivities",		/* 25 */
		"dump",					/* 26 */
		"output",				/* 27 */
		"output_frequency",		/* 28 */
		"selected_output_frequency",	/* 29 */
		"punch_cells",			/* 30 */
		"dump_frequency",		/* 31 */
		"dump_restart",			/* 32 */
		"punch_frequency",		/* 33 */
		"print_frequency",		/* 34 */
		"correct_disp",			/* 35 */
		"initial_time",			/* 36 */
		"warning",				/* 37 */
		"warnings",				/* 38 */
		"thermal_diffusion",	/* 39 */
		"multi_d",				/* 40 */
		"interlayer_d",			/* 41 */
		"porosities",			/* 42 */
		"porosity"				/* 43 */
	};
	int count_opt_list = 44;

	strcpy(file_name, "phreeqc.dmp");
/*
 *   Initialize
 */
	simul_tr++;
	if (simul_tr == 1)
	{
		correct_disp = FALSE;
		old_cells = 0;
		max_cells = 0;
		all_cells = 0;
	}
	else
		old_cells = count_cells;
	count_length = count_disp = count_punch = count_print = count_por = 0;

	length = (LDBLE *) PHRQ_malloc(sizeof(LDBLE));
	if (length == NULL)
		malloc_error();

	disp = (LDBLE *)PHRQ_malloc(sizeof(LDBLE));
	if (disp == NULL)
		malloc_error();

	pors = (LDBLE *)PHRQ_malloc(sizeof(LDBLE));
	if (pors == NULL)
		malloc_error();

	punch_temp = (int *)PHRQ_malloc(sizeof(int));
	if (punch_temp == NULL)
		malloc_error();

	print_temp = (int *) PHRQ_malloc(sizeof(int));
	if (print_temp == NULL)
		malloc_error();

	count_length_alloc = count_disp_alloc = count_por_alloc = 1;
	transport_start = 1;
/*
 *   Read transport number (not currently used)
 */
	ptr = line;
	read_number_description(ptr, &n_user, &n_user_end, &description);
	description = (char *) free_check_null(description);
/*
 *   Set use data to last read
 */
	use.Set_trans_in(true);
/*
 *   Read lines
 */
	opt_save = OPTION_DEFAULT;
	return_value = UNKNOWN;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT)
			opt = opt_save;
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
		case OPTION_DEFAULT:
			input_error++;
			error_msg("Unknown input in TRANSPORT keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* cells */
			sscanf(next_char, "%d", &count_cells);
			opt_save = OPTION_DEFAULT;
			break;
		case 1:				/* shifts */
			if (copy_token(token, &next_char, &l) == DIGIT)
				sscanf(token, "%d", &count_shifts);
			else
			{
				warning_msg
					("Expected the number of shifts. One shift is assumed.");
				count_shifts = 1;
			}
			j = copy_token(token, &next_char, &l);
			if (j != EMPTY)
			{
				if (j == DIGIT)
					sscanf(token, "%d", &ishift);
				else
				{
					input_error++;
					error_msg
						("Expected shift direction, -1, 0, 1. Use -direction instead.",
						 CONTINUE);
					ishift = 1;
				}
			}
			opt_save = OPTION_DEFAULT;
			break;
		case 2:				/* print */
		case 20:				/* print_cells */
			print_temp =
				read_list_ints_range(&next_char, &count_print, TRUE,
									 print_temp);
			opt_save = 2;
			break;
		case 3:				/* selected_output */
		case 29:				/* selected_output_frequency */
		case 33:				/* punch_frequency */
			sscanf(next_char, "%d", &punch_modulus);
			opt_save = OPTION_DEFAULT;
			if (punch_modulus <= 0)
			{
				error_string = sformatf(
						"Punch frequency must be greater than 0. Frequency set to 1000.");
				warning_msg(error_string);
				punch_modulus = 1000;
			}

			break;
		case 4:				/* bcond */
		case 12:				/* bc */
		case 13:				/* boundary_conditions */
			/* first cell boundary condition */
			i = copy_token(token, &next_char, &l);
			str_tolower(token);
			if (i == DIGIT)
			{
				sscanf(token, "%d", &bcon_first);
				if (bcon_first < 1 || bcon_first > 3)
				{
					input_error++;
					error_msg
						("Expected boundary condition to be 'constant' (1), 'closed' (2) , or 'flux' (3).",
						 CONTINUE);
				}
			}
			else if (i == EMPTY)
				bcon_first = 3;
			else if (strstr(token, "co") == token)
				bcon_first = 1;
			else if (strstr(token, "cl") == token)
				bcon_first = 2;
			else if (strstr(token, "f") == token)
				bcon_first = 3;
			else
			{
				input_error++;
				error_msg
					("Expected boundary condition to be 'constant', 'closed', or 'flux'.",
					 CONTINUE);
			}

			/* last cell boundary condition */
			i = copy_token(token, &next_char, &l);
			str_tolower(token);
			if (i == DIGIT)
			{
				sscanf(token, "%d", &bcon_last);
				if (bcon_last < 1 || bcon_last > 3)
				{
					input_error++;
					error_msg
						("Expected boundary condition to be 'constant' (1), 'closed' (2) , or 'flux' (3).",
						 CONTINUE);
				}
			}
			else if (i == EMPTY)
				bcon_last = 3;
			else if (strstr(token, "co") == token)
				bcon_last = 1;
			else if (strstr(token, "cl") == token)
				bcon_last = 2;
			else if (strstr(token, "f") == token)
				bcon_last = 3;
			else
			{
				input_error++;
				error_msg
					("Expected boundary condition to be 'constant', 'closed', or 'flux'.",
					 CONTINUE);
			}
			opt_save = OPTION_DEFAULT;
			break;
		case 5:					/* timest */
		case 14:				/* time_step */
			if (copy_token(token, &next_char, &l) == DIGIT)
				sscanf(token, SCANFORMAT, &timest);
			{
				std::string stdtoken;
				j = copy_token(stdtoken, &next_char);
				if (j == UPPER || j == LOWER)
				{
					timest = Utilities::convert_time(timest, stdtoken, "s");
					j = copy_token(stdtoken, &next_char);
				}
				if (j == DIGIT)
				{
					sscanf(stdtoken.c_str(), SCANFORMAT, &mcd_substeps);
				}
			}
			//if (copy_token(token, &next_char, &l) == DIGIT)
			//	sscanf(token, SCANFORMAT, &mcd_substeps);
			if (mcd_substeps < 1)
			{
				mcd_substeps = 1.0;
				warning_msg("Substep factor in MCD must be >= 1.0\n"
							"mcd_substeps = 1.0 assumed.");
			}
			opt_save = OPTION_DEFAULT;
			break;
		case 6:				/* diffc */
		case 16:				/* diffusion_coefficient */
			sscanf(next_char, SCANFORMAT, &diffc);
			opt_save = OPTION_DEFAULT;
			break;
		case 7:				/* tempr */
		case 15:				/* temp_retardation_factor */
		case 19:				/* temperature_retardation_factor */
		case 39:				/* thermal_diffusion */
			if (copy_token(token, &next_char, &l) == DIGIT)
				sscanf(token, SCANFORMAT, &tempr);
			if (tempr < 1)
			{
				tempr = 1;
				warning_msg
					("Temperature retardation factor < 1 is not possible.\n"
					 "Temperature retardation factor = 1 assumed.");
			}
			j = copy_token(token, &next_char, &l);
			if (j == DIGIT)
				sscanf(token, SCANFORMAT, &heat_diffc);
			opt_save = OPTION_DEFAULT;
			break;
		case 8:				/* length */
		case 24:				/* lengths */
			if (read_line_LDBLEs
				(next_char, &length, &count_length,
				 &count_length_alloc) == ERROR)
			{
				input_error++;
				error_msg("Reading lengths in TRANSPORT keyword.\n",
						  CONTINUE);
			}
			opt_save = 8;
			break;
		case 9:				/* disp */
		case 17:				/* dispersivity */
		case 25:				/* dispersivities */
			if (read_line_LDBLEs
				(next_char, &disp, &count_disp, &count_disp_alloc) == ERROR)
			{
				input_error++;
				error_msg("Reading dispersivities in TRANSPORT keyword.\n",
						  CONTINUE);
			}
			opt_save = 9;
			break;
		case 10:				/* punch */
		case 21:				/* selected_cells */
		case 30:				/* punch_cells */
			punch_temp =
				read_list_ints_range(&next_char, &count_punch, TRUE,
									 punch_temp);
			opt_save = 10;
			break;
		case 11:				/* stagnant */
			if (copy_token(token, &next_char, &l) != EMPTY)
			{
				/* exchange factor */
				if (sscanf(token, "%d", &(stag_data->count_stag)) != 1)
				{
					input_error++;
					error_string = sformatf(
							"Expecting number of stagnant layers.");
					error_msg(error_string, CONTINUE);
					break;
				}

				/* exchange factor */
				j = copy_token(token, &next_char, &l);
				if (j != EMPTY)
				{
					if (sscanf(token, SCANFORMAT, &(stag_data->exch_f)) != 1)
					{
						input_error++;
						error_string = sformatf(
								"Expecting exchange factor for stagnant layers.");
						error_msg(error_string, CONTINUE);
						break;
					}
					copy_token(token, &next_char, &l);
					if (sscanf(token, SCANFORMAT, &(stag_data->th_m)) != 1)
					{
						input_error++;
						error_string = sformatf(
								"Expecting porosity in the mobile zone.");
						error_msg(error_string, CONTINUE);
						break;
					}
					copy_token(token, &next_char, &l);
					if (sscanf(token, SCANFORMAT, &(stag_data->th_im)) != 1)
					{
						input_error++;
						error_string = sformatf(
								"Expecting porosity in the immobile zone.");
						error_msg(error_string, CONTINUE);
						break;
					}
				}
			}
			opt_save = OPTION_DEFAULT;
			break;
		case 18:				/* direction */
		case 22:				/* flow_direction */
		case 23:				/* flow */
			copy_token(token, &next_char, &l);
			str_tolower(token);
			if (strstr(token, "f") == token)
				ishift = 1;
			else if (strstr(token, "b") == token)
				ishift = -1;
			else if (strstr(token, "d") == token)
				ishift = 0;
			else if (strstr(token, "n") == token)
				ishift = 0;
			else
			{
				input_error++;
				error_msg
					("Expected flow direction to be 'forward', 'back', or 'no_flow'.",
					 CONTINUE);
			}
			opt_save = OPTION_DEFAULT;
			break;
		case 26:				/* dump */
			dump_in = TRUE;
			next_char_save = next_char;
			if (copy_token(file_name, &next_char, &l) == EMPTY)
				strcpy(file_name, "phreeqc.dmp");
			else
			{
				string_trim(next_char_save);
				strcpy(file_name, next_char_save);
			}
			opt_save = OPTION_DEFAULT;
			break;
		case 27:				/* output */
		case 28:				/* output_frequency */
		case 34:				/* print_frequency */
			sscanf(next_char, "%d", &print_modulus);
			opt_save = OPTION_DEFAULT;
			if (print_modulus <= 0)
			{
				error_string = sformatf(
						"Print frequency must be greater than 0. Frequency set to 1000.");
				warning_msg(error_string);
				print_modulus = 1000;
			}
			break;
		case 31:				/* dump_frequency */
			dump_in = TRUE;
			if (copy_token(token, &next_char, &l) == DIGIT)
				sscanf(token, "%d", &dump_modulus);
			else
			{
				warning_msg("Expected integer value for dump_frequency.");
				dump_modulus = 0;
			}
			opt_save = OPTION_DEFAULT;
			break;
		case 32:				/* dump_restart */
			dump_in = TRUE;
			if (copy_token(token, &next_char, &l) == DIGIT)
				sscanf(token, "%d", &transport_start);
			else
			{
				warning_msg
					("Expected shift number to start calculations, 1 will be used.");
				transport_start = 1;
			}
			opt_save = OPTION_DEFAULT;
			break;
		case 35:				/* correct_dispersion */
			correct_disp = get_true_false(next_char, TRUE);
			opt_save = OPTION_DEFAULT;
			break;
		case 36:				/* initial_time */
			if (copy_token(token, &next_char, &l) == DIGIT)
				sscanf(token, SCANFORMAT, &initial_total_time);
			{
				std::string stdtoken;
				j = copy_token(stdtoken, &next_char);
				if (j == UPPER || j == LOWER)
				{
					initial_total_time = Utilities::convert_time(initial_total_time, stdtoken, "s");
				}
			}
			opt_save = OPTION_DEFAULT;
			break;
		case 37:				/* warning */
		case 38:				/* warnings */
			transport_warnings = get_true_false(next_char, TRUE);
			break;
		case 40:				/* multicomponent diffusion */
			copy_token(token, &next_char, &l);
			str_tolower(token);
			if (strstr(token, "f") == token)
				multi_Dflag = 0;
			else if (strstr(token, "t") == token)
				multi_Dflag = 1;
			else
			{
				input_error++;
				error_msg
					("Expected multicomponent diffusion flag: 'true' or 'false'.",
					 CONTINUE);
			}
			default_Dw = 1e-9;
			multi_Dpor = 0.3;
			multi_Dpor_lim = 0.0;
			multi_Dn = 1.0;
			if (copy_token(token, &next_char, &l) == EMPTY)
				break;
			else
			{
				/* default species diffusion coeff */
				if (sscanf(token, SCANFORMAT, &default_Dw) != 1)
				{
					input_error++;
					error_string = sformatf(
							"Expected default species diffusion coefficient in water at 25oC, m2/s.");
					error_msg(error_string, CONTINUE);
					break;
				}
			}
			if (copy_token(token, &next_char, &l) == EMPTY)
				break;
			else
			{
				/* porosity */
				if (sscanf(token, SCANFORMAT, &multi_Dpor) != 1)
				{
					input_error++;
					error_string = sformatf(
							"Expected porosity to calculate diffusion coefficient.");
					error_msg(error_string, CONTINUE);
					break;
				}
			}
			if (copy_token(token, &next_char, &l) == EMPTY)
				break;
			else
			{
				/* porosity */
				if (sscanf(token, SCANFORMAT, &multi_Dpor_lim) != 1)
				{
					input_error++;
					error_string = sformatf(
							"Expected porosity limit for diffusive transport.");
					error_msg(error_string, CONTINUE);
					break;
				}
			}
			if (copy_token(token, &next_char, &l) == EMPTY)
				break;
			else
			{
				if (sscanf(token, SCANFORMAT, &multi_Dn) != 1)
				{
					input_error++;
					error_string = sformatf(
							"Expected exponent for porosity reduction of diffusion coefficient (Dp = Dw * (por)^n).");
					error_msg(error_string, CONTINUE);
					break;
				}
			}
			opt_save = OPTION_DEFAULT;
			break;
		case 41:				/* interlayer diffusion */
			copy_token(token, &next_char, &l);
			str_tolower(token);
			if (strstr(token, "f") == token)
				interlayer_Dflag = 0;
			else if (strstr(token, "t") == token)
				interlayer_Dflag = 1;
			else
			{
				input_error++;
				error_msg
					("Expected interlayer diffusion flag: 'true' or 'false'.",
					 CONTINUE);
			}
			interlayer_Dpor = 0.1;
			interlayer_Dpor_lim = 0.0;
			interlayer_tortf = 100.0;
			if (copy_token(token, &next_char, &l) == EMPTY)
				break;
			else
			{
				/* porosity */
				if (sscanf(token, SCANFORMAT, &interlayer_Dpor) != 1)
				{
					input_error++;
					error_string = sformatf( "Expected interlayer porosity.");
					error_msg(error_string, CONTINUE);
					break;
				}
			}
			if (copy_token(token, &next_char, &l) == EMPTY)
				break;
			else
			{
				/* porosity limit */
				if (sscanf(token, SCANFORMAT, &interlayer_Dpor_lim) != 1)
				{
					input_error++;
					error_string = sformatf(
							"Expected interlayer porosity limit for diffusive transport.");
					error_msg(error_string, CONTINUE);
					break;
				}
			}
			if (copy_token(token, &next_char, &l) == EMPTY)
				break;
			else
			{
				if (sscanf(token, SCANFORMAT, &interlayer_tortf) != 1)
				{
					input_error++;
					error_string = sformatf(
							"Expected interlayer tortuosity factor (Dp = Dw /t_f).");
					error_msg(error_string, CONTINUE);
					break;
				}
			}
			opt_save = OPTION_DEFAULT;
			break;
		case 42:				/* porosities */
		case 43:				/* porosity */
			if (read_line_LDBLEs
				(next_char, &pors, &count_por,
				&count_por_alloc) == ERROR)
			{
				input_error++;
				error_msg("Reading porosities in TRANSPORT keyword.\n",
					CONTINUE);
			}
			opt_save = 42;
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
/*
 *   Determine number of cells
 */
	max_cells = count_cells;
	if (count_length > max_cells)
		max_cells = count_length;
	if (count_disp > max_cells)
		max_cells = count_disp;
	//if (count_por > max_cells)
	//	max_cells = count_por;
	if (max_cells > count_cells)
	{
		if (max_cells == count_length)
		{
			sprintf(token,
					"Number of cells is increased to number of 'lengths' %d.",
					count_length);
			warning_msg(token);
		}
		else if (max_cells == count_disp)
		{
			sprintf(token,
				"Number of cells is increased to number of dispersivities %d.",
				count_disp);
			warning_msg(token);
		}
		//else
		//{
		//	sprintf(token,
		//		"Number of cells is increased to number of porosities %d.",
		//		count_por);
		//	warning_msg(token);
		//}
	}
/*
 *   Allocate space for cell_data
 */
	cell_data = (struct cell_data *) PHRQ_realloc(cell_data,
		(size_t) (max_cells *	(1 + stag_data->count_stag) + 1) * sizeof(struct cell_data));
	if (cell_data == NULL)
		malloc_error();

	// initialize new cells
	int all_cells_now = max_cells * (1 + stag_data->count_stag) + 1;
	if (all_cells_now > all_cells)
	{
		for (int i = all_cells; i < all_cells_now; i++)
		{
			cell_data[i].length = 1.0;
			cell_data[i].mid_cell_x = 1.0;
			cell_data[i].disp = 1.0;
			cell_data[i].temp = 25.0;
			cell_data[i].por = 0.3;   
			cell_data[i].por_il = 0.01;
			cell_data[i].punch = FALSE;
			cell_data[i].print = FALSE;
		}
		all_cells = all_cells_now;
	}

/*
 *   Fill in data for lengths
 */
	if (count_length == 0)
	{
		if (old_cells < max_cells)
		{
			error_string = sformatf(
					"No cell-lengths were read; length = 1 m assumed.");
			warning_msg(error_string);
			for (i = 0; i < max_cells; i++)
				cell_data[i].length = 1.0;
		}
	}
	else
	{
		for (i = 0; i < count_length; i++)
		{
			cell_data[i].length = length[i];
		}
		if (max_cells > count_length)
		{
			error_string = sformatf(
					"Cell-lengths were read for %d cells. Last value is used till cell %d.",
					count_length, max_cells);
			warning_msg(error_string);
			for (i = count_length - 1; i < max_cells; i++)
				cell_data[i].length = length[count_length - 1];
		}
	}
	cell_data[0].mid_cell_x = cell_data[0].length / 2;
	for (i = 1; i < max_cells; i++)
	{
		cell_data[i].mid_cell_x = cell_data[i - 1].mid_cell_x +
			(cell_data[i - 1].length + cell_data[i].length) / 2;
	}
	cell_data[max_cells].mid_cell_x =
		cell_data[max_cells - 1].mid_cell_x + cell_data[max_cells - 1].length;
/*
 *   Fill in data for dispersivities
 */
	if (count_disp == 0)
	{
		if (old_cells < max_cells)
		{
			error_string = sformatf(
				"No dispersivities were read; disp = 0 assumed.");
			warning_msg(error_string);
			for (i = 0; i < max_cells; i++)
				cell_data[i].disp = 0.0;
		}
	}
	else
	{
		for (i = 0; i < count_disp; i++)
			cell_data[i].disp = disp[i];
		if (max_cells > count_disp)
		{
			error_string = sformatf(
				"Dispersivities were read for %d cells. Last value is used till cell %d.",
				count_disp, max_cells);
			warning_msg(error_string);
			for (i = count_disp - 1; i < max_cells; i++)
				cell_data[i].disp = disp[count_disp - 1];
		}
	}
/*
 *   Fill in data for porosities
 */
	if (count_por == 0)
	{
		if (old_cells < all_cells && multi_Dflag)
		{
			multi_Dpor = (multi_Dpor < 1e-10 ? 1e-10 : multi_Dpor);
			if (multi_Dpor > 1e-10)
				error_string = sformatf(
				"No porosities were read; used the value %8.2e from -multi_D.", multi_Dpor);
			else
				error_string = sformatf(
				"No porosities were read; used the minimal value %8.2e from -multi_D.", multi_Dpor);
			warning_msg(error_string);
			//for (i = old_cells + 1; i < all_cells; i++)
			for (i = old_cells; i < all_cells; i++)
				cell_data[i].por = multi_Dpor;
		}
	}
	else
	{
		if ((stag_data->exch_f > 0) && (stag_data->count_stag == 1))
		{
			error_string = sformatf(
				"Mobile porosities were read, but mobile/immobile porosity was also defined in -stagnant. Using the values from -stagnant for mobile/immobile exchange and tortuosity factors.");
			warning_msg(error_string);
			for (i = 1; i <= max_cells; i++)
				cell_data[i].por = stag_data->th_m;
			for (i++; i <= 2 * max_cells + 1; i++)
				cell_data[i].por = stag_data->th_im;
		}
		else
		{
			for (i = 0; i < count_por; i++)
				cell_data[i].por = pors[i];
			if (max_cells > count_por)
			{
				error_string = sformatf(
					"Porosities were read for %d cells. Last value is used till cell %d.",
					count_por, all_cells - 1);
				warning_msg(error_string);
				for (i = count_por - 1; i < all_cells; i++)
					cell_data[i].por = pors[count_por - 1];
			}
		}
	}
	if (interlayer_Dflag && !multi_Dflag)
	{
		input_error++;
		error_string = sformatf(
			"-multi_D must be defined, when -interlayer_D true.");
		error_msg(error_string, CONTINUE);

	}
	for (i = 0; i < all_cells; i++)
	{
		interlayer_Dpor = (interlayer_Dpor < 1e-10 ? 1e-10 : interlayer_Dpor);
		cell_data[i].por_il = interlayer_Dpor;
	}
	count_cells = max_cells;
/*
 *  Account for stagnant cells
 */
	if (stag_data->count_stag > 0)
	{
		max_cells = count_cells * (1 + stag_data->count_stag) + 1;
		for (i = 0; i < count_cells; i++)
		{
			for (l = 1; l <= stag_data->count_stag; l++)
				cell_data[i + 1 + l * count_cells].mid_cell_x =
					cell_data[i].mid_cell_x;
		}
	}
/*
 *   Fill in data for punch
 */
	if (count_punch != 0)
	{
		for (i = 0; i < max_cells; i++)
			cell_data[i].punch = FALSE;
		for (i = 0; i < count_punch; i++)
		{
			if (punch_temp[i] > max_cells || punch_temp[i] < 1)
			{
				error_string = sformatf(
						"Cell number for punch is out of range, %d. Request ignored.",
						punch_temp[i]);
				warning_msg(error_string);
			}
			else
				cell_data[punch_temp[i] - 1].punch = TRUE;
		}
	}
	else if (simul_tr == 1)
		for (i = 0; i < max_cells; i++)
			cell_data[i].punch = TRUE;
/*
 *   Fill in data for print
 */
	if (count_print != 0)
	{
		for (i = 0; i < max_cells; i++)
			cell_data[i].print = FALSE;
		for (i = 0; i < count_print; i++)
		{
			if (print_temp[i] > max_cells || print_temp[i] < 1)
			{
				error_string = sformatf(
						"Cell number for print is out of range, %d. Request ignored.",
						print_temp[i]);
				warning_msg(error_string);
			}
			else
				cell_data[print_temp[i] - 1].print = TRUE;
		}
	}
	else if (simul_tr == 1)
		for (i = 0; i < max_cells; i++)
			cell_data[i].print = TRUE;
//#define OLD_POROSITY
#if defined(OLD_POROSITY)
/*
 *   Fill in porosities
 */
	if (interlayer_Dflag && !multi_Dflag)
	{
		input_error++;
		error_string = sformatf(
				"-multi_D must be defined, when -interlayer_D true.");
		error_msg(error_string, CONTINUE);

	}
	for (i = 0; i < max_cells; i++)
	{
		multi_Dpor = (multi_Dpor < 1e-10 ? 1e-10 : multi_Dpor);    //Fix for Jenkins !!!!!!!!!!!!
		//if (cell_data[i].por < 0)
		{
			cell_data[i].por = multi_Dpor;                         //Fix for Jenkins !!!!!!!!!!!!
		}
		interlayer_Dpor = (interlayer_Dpor < 1e-10 ? 1e-10 : interlayer_Dpor);
		cell_data[i].por_il = interlayer_Dpor;
	}
#endif
	//{
	//	for (int i = 0; i < all_cells; i++)
	//	{
	//		std::cerr << i << "  " << cell_data[i].por << std::endl;
	//	}
	//}
	
/*
 *   Calculate dump_modulus
 */
	if (dump_in == TRUE)
	{
		if (dump_modulus == 0)
		{
			warning_msg
				("Expected dump_modulus. Value of 'shifts/2' will be used.");
			dump_modulus = count_shifts / 2;
			if (dump_modulus == 0)
				dump_modulus = 1;
		}
		if (transport_start > count_shifts)
		{
			input_error++;
			error_string = sformatf(
					"Starting shift for transport, %d, is greater than number of shifts, %d.",
					transport_start, count_shifts);
			error_msg(error_string, CONTINUE);
		}
	}
/*
 *  Check boundary conditions
 */
	if ((ishift != 0) && ((bcon_first == 2) || (bcon_last == 2)))
	{
		warning_msg
			("Boundary condition = 'closed' not possible with advective transport.\n\t Boundary condition = 'flux' assumed.");
		if (bcon_first == 2)
			bcon_first = 3;
		if (bcon_last == 2)
			bcon_last = 3;
	}
/*
 *  Retain data from previous run
 */
	if (simul_tr > 1)
	{
		if ((count_length == 0) && (count_disp == 0) && (count_por == 0))
			dup_print("Column data retained from former run", TRUE);
	}
/*
 *  Check heat_diffc
 */
	if (heat_diffc < 0)
		heat_diffc = diffc;
	else if (stag_data->count_stag == 1)
	{
		if (stag_data->exch_f > 0)
		{
			if (diffc <= 0 && heat_diffc > 0)
			{
				input_error++;
				error_string = sformatf(
						"Must enter diffusion coefficient (-diffc) when modeling thermal diffusion.");
				error_msg(error_string, CONTINUE);
			}
			else if (heat_diffc > diffc)
			{
				error_string = sformatf(
						"Thermal diffusion is calculated assuming exchange factor was for\n\t effective (non-thermal) diffusion coefficient = %e.",
						(double) diffc);
				warning_msg(error_string);
			}
		}
		else
		{
			if (heat_diffc > diffc)
			{
				input_error++;
				error_string = sformatf(
						"Must enter value for mobile/stagnant exchange factor when modeling thermal diffusion.");
				error_msg(error_string, CONTINUE);
			}
		}
	}
	else if (stag_data->count_stag > 1 && heat_diffc > diffc)
	{
		input_error++;
		error_string = sformatf(
				"Only one stagnant layer permitted (-stag) when modeling thermal diffusion.");
		error_msg(error_string, CONTINUE);
	}
/*
 *   free storage for length, disp, punch
 */
	length = (LDBLE *) free_check_null(length);
	disp = (LDBLE *) free_check_null(disp);
	pors = (LDBLE *) free_check_null(pors);
	punch_temp = (int *) free_check_null(punch_temp);
	print_temp = (int *) free_check_null(print_temp);

	if (dump_in == TRUE)
	{
		dump_file_name_cpp.clear();
		dump_file_name_cpp.append(file_name);
	}
	return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_line_LDBLEs(char *next_char, LDBLE ** d, int *count_d, int *count_alloc)
/* ---------------------------------------------------------------------- */
{
	int i, j, l, n;
	LDBLE value;
	char token[MAX_LENGTH];

	for (;;)
	{
		j = copy_token(token, &next_char, &l);
		if (j == EMPTY)
			break;
		if (j != DIGIT)
			return (ERROR);
		if (replace("*", " ", token) == TRUE)
		{
			if (sscanf(token, "%d" SCANFORMAT, &n, &value) != 2)
				return (ERROR);
		}
		else
		{
			sscanf(token, SCANFORMAT, &value);
			n = 1;
		}
		for (;;)
		{
			if ((*count_d) + n > (*count_alloc))
			{
				*count_alloc *= 2;
				*d = (LDBLE *) PHRQ_realloc(*d,
											(size_t) (*count_alloc) *
											sizeof(LDBLE));
				if (*d == NULL)
					malloc_error();
			}
			else
				break;
		}
		for (i = 0; i < n; i++)
			(*d)[(*count_d) + i] = value;
		*count_d += n;
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
dump_cpp(void)
/* ---------------------------------------------------------------------- */
{
/*
 * dumps solution compositions to file
 */

	int j, l;

	if (dump_in == FALSE || pr.dump == FALSE)
		return (OK);

	cxxStorageBin phreeqcBin(phrq_io);
	phreeqc2cxxStorageBin(phreeqcBin);

	std::ofstream fs(dump_file_name_cpp.c_str());
	if (!fs.is_open())
	{
		error_string = sformatf( "Can`t open file, %s.", dump_file_name_cpp.c_str());
		input_error++;
		error_msg(error_string, CONTINUE);
		return (OK);
	}
	
	fs << "# Dumpfile" << "\n" << "# Transport simulation " << simul_tr << "  Shift " << transport_step << "\n" << "#" << "\n";
	phreeqcBin.dump_raw(fs, 0);
	fs << "END" << "\n";

	char token[MAX_LENGTH];
	sprintf(token, "KNOBS\n");
	fs << token; 
	sprintf(token, "\t-iter%15d\n", itmax);
	fs << token; 
	sprintf(token, "\t-tol %15.3e\n", (double) ineq_tol);
	fs << token; 
	sprintf(token, "\t-step%15.3e\n", (double) step_size);
	fs << token;
	sprintf(token, "\t-pe_s%15.3e\n", (double) pe_step_size);
	fs << token;
	sprintf(token, "\t-diag      ");
	fs << token;
	if (diagonal_scale == TRUE)
	{
		sprintf(token, "true\n");
		fs << token;
	}
	else
	{
		sprintf(token, "false\n");
		fs << token;
	}
	std::map < int, SelectedOutput >::iterator so_it = SelectedOutput_map.begin();
	for ( ; so_it != SelectedOutput_map.end(); so_it++)
	{
		current_selected_output = &(so_it->second);

		sprintf(token, "SELECTED_OUTPUT %d\n", current_selected_output->Get_n_user());
		fs << token ;
		//sprintf(token, "\t-file  %-15s\n", "sel_o$$$.prn");
		//fs << token;
		fs << "\t-file  " << "sel_o$$$" << current_selected_output->Get_n_user() << ".prn\n";
		//if (punch.count_totals != 0)
		if (current_selected_output->Get_totals().size() > 0)
		{
			sprintf(token, "\t-tot ");
			fs << token;
			for (size_t i = 0; i < current_selected_output->Get_totals().size(); i++)
			{
				sprintf(token, "  %s", current_selected_output->Get_totals()[i].first.c_str());
				fs << token;
			}
			sprintf(token, "\n");
			fs << token;
		}
		if (current_selected_output->Get_molalities().size() > 0)
		{
			sprintf(token, "\t-mol ");
			fs << token;
			for (size_t i = 0; i < current_selected_output->Get_molalities().size(); i++)
			{
				sprintf(token, "  %s", current_selected_output->Get_molalities()[i].first.c_str());
				fs << token;
			}
			sprintf(token, "\n");
			fs << token;
		}
		if (current_selected_output->Get_activities().size() > 0)
		{
			sprintf(token, "\t-act ");
			fs << token;
			for (size_t i = 0; i < current_selected_output->Get_activities().size(); i++)
			{
				sprintf(token, "  %s", current_selected_output->Get_activities()[i].first.c_str());
				fs << token;
			}
			sprintf(token, "\n");
			fs << token;
		}
		if (current_selected_output->Get_pure_phases().size() > 0)
		{
			sprintf(token, "\t-equ ");
			fs << token;
			for (size_t i = 0; i < current_selected_output->Get_pure_phases().size(); i++)
			{
				sprintf(token, "  %s", current_selected_output->Get_pure_phases()[i].first.c_str());
				fs << token;
			}
			sprintf(token, "\n");
			fs << token;
		}
		if (current_selected_output->Get_si().size() > 0)
		{
			sprintf(token, "\t-si ");
			fs << token;
			for (size_t i = 0; i < current_selected_output->Get_si().size(); i++)
			{
				sprintf(token, "  %s", current_selected_output->Get_si()[i].first.c_str());
				fs << token;
			}
			sprintf(token, "\n");
			fs << token;
		}
		if (current_selected_output->Get_gases().size() > 0)
		{
			sprintf(token, "\t-gas ");
			fs << token;
			for (size_t i = 0; i < current_selected_output->Get_gases().size(); i++)
			{
				sprintf(token, "  %s", current_selected_output->Get_gases()[i].first.c_str());
				fs << token;
			}
			sprintf(token, "\n");
			fs << token;
		}
		if (current_selected_output->Get_s_s().size() > 0)
		{
			sprintf(token, "\t-solid_solutions ");
			fs << token;
			for (size_t i = 0; i < current_selected_output->Get_s_s().size(); i++)
			{
				sprintf(token, "  %s", current_selected_output->Get_s_s()[i].first.c_str());
				fs << token;
			}
			sprintf(token, "\n");
			fs << token;
		}
		if (current_selected_output->Get_kinetics().size() > 0)
		{
			sprintf(token, "\t-kin ");
			fs << token;
			for (size_t i = 0; i < current_selected_output->Get_kinetics().size(); i++)
			{
				sprintf(token, "  %s", current_selected_output->Get_kinetics()[i].first.c_str());
				fs << token;
			}
			sprintf(token, "\n");
			fs << token;
		}
	}
	sprintf(token, "TRANSPORT\n");
	fs << token;
	sprintf(token, "\t-cells %6d\n", count_cells);
	fs << token;
	sprintf(token, "\t-shifts%6d%6d\n", count_shifts, ishift);
	fs << token;
	sprintf(token, "\t-output_frequency %6d\n", print_modulus);
	fs << token;
	sprintf(token, "\t-selected_output_frequency %6d\n",
		punch_modulus);
	fs << token;
	sprintf(token, "\t-bcon  %6d%6d\n", bcon_first, bcon_last);
	fs << token;
	sprintf(token, "\t-timest %13.5e\n", (double) timest);
	fs << token;
	if (!high_precision)
	{
		sprintf(token, "\t-diffc  %13.5e\n", (double) diffc);
		fs << token;
	}
	else
	{
		sprintf(token, "\t-diffc  %20.12e\n", (double) diffc);
		fs << token;
	}
	sprintf(token, "\t-tempr  %13.5e\n", (double) tempr);
	fs << token;
	if (correct_disp == TRUE)
	{
		sprintf(token, "\t-correct_disp %s\n", "True");
		fs << token;
	}
	else
	{
		sprintf(token, "\t-correct_disp %s\n", "False");
		fs << token;
	}
	sprintf(token, "\t-length\n");
	fs << token;
	for (int i = 0; i < count_cells; i++)
	{
		sprintf(token, "%12.3e", (double) cell_data[i].length);
		fs << token;
		if (i > 0 && (i % 8) == 0)
		{
			sprintf(token, "\n");
			fs << token;
		}
	}
	sprintf(token, "\n");
	fs << token;
	sprintf(token, "\t-disp\n");
	fs << token;
	for (int i = 0; i < count_cells; i++)
	{
		if (!high_precision)
		{
			sprintf(token, "%12.3e", (double) cell_data[i].disp);
			fs << token;
		}
		else
		{
			sprintf(token, "%20.12e", (double) cell_data[i].disp);
			fs << token;
		}
		if (i > 0 && (i % 8) == 0)
		{
			sprintf(token, "\n");
			fs << token;
		}
	}
	sprintf(token, "\n");
	fs << token;
	sprintf(token, "\t-punch_cells");
	fs << token;
	if (stag_data->count_stag > 0)
		j = 1 + (1 + stag_data->count_stag) * count_cells;
	else
		j = count_cells;
	l = 0;
	for (int i = 0; i < j; i++)
	{
		if (cell_data[i].punch != TRUE)
			continue;
		sprintf(token, "  %d", i + 1);
		fs << token;
		l++;
		if ((l % 20) == 0)
		{
			sprintf(token, "\n");
			fs << token;
		}
	}
	sprintf(token, "\n");
	fs << token;
	sprintf(token, "\t-print_cells");
	fs << token;
	if (stag_data->count_stag > 0)
		j = 1 + (1 + stag_data->count_stag) * count_cells;
	else
		j = count_cells;
	l = 0;
	for (int i = 0; i < j; i++)
	{
		if (cell_data[i].print != TRUE)
			continue;
		sprintf(token, "  %d", i + 1);
		fs << token;
		l++;
		if ((l % 20) == 0)
		{
			sprintf(token, "\n");
			fs << token;
		}
	}
	sprintf(token, "\n");
	fs << token;
	sprintf(token, "\t-dump            $$$.dmp\n");
	fs << token;
	sprintf(token, "\t-dump_frequency  %d\n", dump_modulus);
	fs << token;
	sprintf(token, "\t-dump_restart    %d\n", transport_step + 1);
	fs << token;

#if defined MULTICHART
	// user graphs
	chart_handler.dump(fs, 0);
#endif

	sprintf(token, "END\n");
	fs << token;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
dump(void)
/* ---------------------------------------------------------------------- */
{
/*
 * dumps solution compositions to file
 */
	if (dump_in == FALSE || pr.dump == FALSE)
		return (OK);

	dump_cpp();
	return OK;

}
