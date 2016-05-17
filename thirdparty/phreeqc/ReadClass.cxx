#ifdef _DEBUG
#pragma warning(disable : 4786)	// disable truncation warning (Only used by debugger)
#endif
#include <iostream>
#include <fstream>
#include <sstream>

#include "Utils.h"
#include "Phreeqc.h"
#include "Parser.h"
#include "Solution.h"
#include "Exchange.h"
#include "Surface.h"
#include "PPassemblage.h"
#include "cxxKinetics.h"
#include "SSassemblage.h"
#include "GasPhase.h"
#include "Reaction.h"
#include "cxxMix.h"
#include "Temperature.h"
#include "dumper.h"
#include "runner.h"
#include "cxxMix.h"
#include "Surface.h"
#include "phqalloc.h"

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_dump(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads DUMP data block
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
	int return_value;
	/*
	 *  Make parser
	 */
	std::istringstream iss_in;
	return_value = streamify_to_next_keyword(iss_in);
	CParser parser(iss_in, phrq_io);
	assert(!reading_database());

	//For testing, need to read line to get started
	parser.set_echo_file(CParser::EO_NONE);
	std::vector < std::string > vopts;
	std::istream::pos_type next_char;
	parser.get_option(vopts, next_char);

	if (pr.echo_input == FALSE)
	{
		parser.set_echo_file(CParser::EO_NONE);
	}
	else
	{
		parser.set_echo_file(CParser::EO_NOKEYWORDS);
	}

	dump_info.Read(parser);

	// Need to output the next keyword
	if (return_value == OPTION_KEYWORD) echo_msg(sformatf( "\t%s\n", line));
	return (return_value);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_delete(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads DELETE data block
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
	int return_value;
	/*
	 *  Make parser
	 */
	std::istringstream iss_in;
	return_value = streamify_to_next_keyword(iss_in);
	CParser parser(iss_in, phrq_io);
	//assert(!reading_database());

	//For testing, need to read line to get started
	parser.set_echo_file(CParser::EO_NONE);
	std::vector < std::string > vopts;
	std::istream::pos_type next_char;
	parser.get_option(vopts, next_char);

	if (pr.echo_input == FALSE)
	{
		parser.set_echo_file(CParser::EO_NONE);
	}
	else
	{
		parser.set_echo_file(CParser::EO_NOKEYWORDS);
	}

	delete_info.Read(parser);


	// Need to output the next keyword
	if (return_value == OPTION_KEYWORD) echo_msg(sformatf( "\t%s\n", line));
	return (return_value);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_run_cells(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads DELETE data block
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
	int return_value;
	/*
	 *  Make parser
	 */
	std::istringstream iss_in;
	return_value = streamify_to_next_keyword(iss_in);
	CParser parser(iss_in, phrq_io);
	assert(!reading_database());

	//For testing, need to read line to get started
	parser.set_echo_file(CParser::EO_NONE);
	std::vector < std::string > vopts;
	std::istream::pos_type next_char;
	parser.get_option(vopts, next_char);

	if (pr.echo_input == FALSE)
	{
		parser.set_echo_file(CParser::EO_NONE);
	}
	else
	{
		parser.set_echo_file(CParser::EO_NOKEYWORDS);
	}

	runner r(parser, phrq_io);
	run_info = r;


	// Need to output the next keyword
	if (return_value == OPTION_KEYWORD) echo_msg(sformatf( "\t%s\n", line));
	return (return_value);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
streamify_to_next_keyword(std::istringstream & lines)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads to next keyword or eof
 *
 *   Returns:
 *       OPTION_KEYWORD
 *       OPTION_EOF
 *       
 */
	// Handle echo
	int save_echo_input = pr.echo_input;
	pr.echo_input = FALSE;

	std::string accumulate(line);
	accumulate.append("\n");
	int j;
	for (;;)
	{
		j = check_line("Streamify", FALSE, TRUE, TRUE, FALSE);
		if (j == EOF)
		{						/* end of file */
			break;
		}
		else if (j == KEYWORD)
		{						/* keyword */
			break;
		}
		else
		{
			accumulate.append(line);
			accumulate.append("\n");
		}
	}

	lines.str(accumulate);
	pr.echo_input = save_echo_input;
	if (j == EOF) return (OPTION_EOF);
	if (j == KEYWORD) return (OPTION_KEYWORD);


	return (OPTION_ERROR);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
dump_entities(void)
/* ---------------------------------------------------------------------- */
{
	if (!dump_info.Get_on() || pr.dump == FALSE)
	{
		return(OK);
	}
	dump_info.Set_on(false);
	if (!dump_info.Get_bool_any())
	{
		return(OK);
	}

	if (this->phrq_io)
	{
		std::ios_base::openmode mode = std::ios_base::out;
		if (dump_info.Get_append())
		{
			mode = std::ios_base::app;
		}
		if (this->phrq_io->dump_open(dump_info.Get_file_name().c_str(), mode))
		{
			dump_ostream(*this->phrq_io->Get_dump_ostream());
			this->phrq_io->dump_close();
		}
		else
		{
			error_string = sformatf( "Unable to open dump file \"%s\"", dump_info.Get_file_name().c_str());
			error_msg(error_string, STOP);
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
delete_entities(void)
/* ---------------------------------------------------------------------- */
{
	if (!delete_info.Get_solution().Get_defined() &&
		!delete_info.Get_pp_assemblage().Get_defined() &&
		!delete_info.Get_exchange().Get_defined() &&
		!delete_info.Get_surface().Get_defined() &&
		!delete_info.Get_ss_assemblage().Get_defined() &&
		!delete_info.Get_gas_phase().Get_defined() &&
		!delete_info.Get_kinetics().Get_defined() &&
		!delete_info.Get_mix().Get_defined() &&
		!delete_info.Get_reaction().Get_defined() &&
		!delete_info.Get_temperature().Get_defined() &&
		!delete_info.Get_pressure().Get_defined())
	{
		return(OK);
	}

	// solutions
	if (delete_info.Get_solution().Get_defined())
	{
		if (delete_info.Get_solution().Get_numbers().size() == 0)
		{
			Rxn_solution_map.clear();
		}
		else
		{
			std::set < int >::iterator it;
			for (it = delete_info.Get_solution().Get_numbers().begin(); it != delete_info.Get_solution().Get_numbers().end(); it++)
			{
				Rxn_solution_map.erase(*it);
			}
		}
	}

	// pp_assemblages
	if (delete_info.Get_pp_assemblage().Get_defined())
	{
		if (delete_info.Get_pp_assemblage().Get_numbers().size() == 0)
		{
			Rxn_pp_assemblage_map.clear();
		}
		else
		{
			std::set < int >::iterator it;
			for (it = delete_info.Get_pp_assemblage().Get_numbers().begin(); it != delete_info.Get_pp_assemblage().Get_numbers().end(); it++)
			{
				Rxn_pp_assemblage_map.erase(*it);
			}
		}
	}
	// exchangers
	if (delete_info.Get_exchange().Get_defined())
	{
		if (delete_info.Get_exchange().Get_numbers().size() == 0)
		{
			Rxn_exchange_map.clear();
		}
		else
		{
			std::set < int >::iterator it;
			for (it = delete_info.Get_exchange().Get_numbers().begin(); it != delete_info.Get_exchange().Get_numbers().end(); it++)
			{
				Rxn_exchange_map.erase(*it);
			}
		}
	}
	// surfaces
	if (delete_info.Get_surface().Get_defined())
	{
		if (delete_info.Get_surface().Get_numbers().size() == 0)
		{
			Rxn_surface_map.clear();
		}
		else
		{
			std::set < int >::iterator it;
			for (it = delete_info.Get_surface().Get_numbers().begin(); it != delete_info.Get_surface().Get_numbers().end(); it++)
			{
				Rxn_surface_map.erase(*it);
			}
		}
	}
	// ss_assemblages
	if (delete_info.Get_ss_assemblage().Get_defined())
	{
		if (delete_info.Get_ss_assemblage().Get_numbers().size() == 0)
		{
			Rxn_ss_assemblage_map.clear();
		}
		else
		{
			std::set < int >::iterator it;
			for (it = delete_info.Get_ss_assemblage().Get_numbers().begin(); it != delete_info.Get_ss_assemblage().Get_numbers().end(); it++)
			{
				Rxn_ss_assemblage_map.erase(*it);
			}
		}
	}
	// gas_phases
	if (delete_info.Get_gas_phase().Get_defined())
	{
		if (delete_info.Get_gas_phase().Get_numbers().size() == 0)
		{
			Rxn_gas_phase_map.clear();
		}
		else
		{
			std::set < int >::iterator it;
			for (it = delete_info.Get_gas_phase().Get_numbers().begin(); it != delete_info.Get_gas_phase().Get_numbers().end(); it++)
			{
				Rxn_gas_phase_map.erase(*it);
			}
		}
	}
	// kinetics
	if (delete_info.Get_kinetics().Get_defined())
	{
		if (delete_info.Get_kinetics().Get_numbers().size() == 0)
		{
			Rxn_kinetics_map.clear();
		}
		else
		{
			std::set < int >::iterator it;
			for (it = delete_info.Get_kinetics().Get_numbers().begin(); it != delete_info.Get_kinetics().Get_numbers().end(); it++)
			{
				Rxn_kinetics_map.erase(*it);
			}
		}
	}
	// mixes
	if (delete_info.Get_mix().Get_defined())
	{
		if (delete_info.Get_mix().Get_numbers().size() == 0)
		{
			Rxn_mix_map.clear();
		}
		else
		{
			std::set < int >::iterator it;
			for (it = delete_info.Get_mix().Get_numbers().begin(); it != delete_info.Get_mix().Get_numbers().end(); it++)
			{
				Rxn_mix_map.erase(*it);
			}
		}
	}
	// reactions
	if (delete_info.Get_reaction().Get_defined())
	{
		if (delete_info.Get_reaction().Get_numbers().size() == 0)
		{
			Rxn_reaction_map.clear();
		}
		else
		{
			std::set < int >::iterator it;
			for (it = delete_info.Get_reaction().Get_numbers().begin(); it != delete_info.Get_reaction().Get_numbers().end(); it++)
			{
				Rxn_reaction_map.erase(*it);
			}
		}
	}
	// temperatures
	if (delete_info.Get_temperature().Get_defined())
	{
		if (delete_info.Get_temperature().Get_numbers().size() == 0)
		{
			Rxn_temperature_map.clear();
		}
		else
		{
			std::set < int >::iterator it;
			for (it = delete_info.Get_temperature().Get_numbers().begin(); it != delete_info.Get_temperature().Get_numbers().end(); it++)
			{
				Rxn_temperature_map.erase(*it);
			}
		}
	}
	// pressures
	if (delete_info.Get_pressure().Get_defined())
	{
		if (delete_info.Get_pressure().Get_numbers().size() == 0)
		{
			Rxn_pressure_map.clear();
		}
		else
		{
			std::set < int >::iterator it;
			for (it = delete_info.Get_pressure().Get_numbers().begin(); it != delete_info.Get_pressure().Get_numbers().end(); it++)
			{
				Rxn_pressure_map.erase(*it);
			}
		}
	}
	// Turn off delete until next read
	delete_info.SetAll(false);
	return (OK);
}
#ifdef USE_OPTIMIZED_BUT_NOT_MUCH
/* ---------------------------------------------------------------------- */
int Phreeqc::
run_as_cells(void)
/* ---------------------------------------------------------------------- */
{
	struct save save_data;
	LDBLE kin_time;
	int count_steps, use_mix;
	char token[2 * MAX_LENGTH];

	state = REACTION;
	if (run_info.Get_cells().Get_numbers().size() == 0 ||
		!(run_info.Get_cells().Get_defined())) return(OK);

	// running cells
	run_info.Set_run_cells(true);

	dup_print("Beginning of run as cells.", TRUE);
	LDBLE initial_total_time_save;
	if (run_info.Get_start_time() != NA)
	{
		initial_total_time_save = run_info.Get_start_time();
	}
	else
	{
		initial_total_time_save = initial_total_time;
	}

	std::set < int >::iterator it = run_info.Get_cells().Get_numbers().begin();

	for ( ; it != run_info.Get_cells().Get_numbers().end(); it++)
	{
		int i = *it;
		if (i < 0) continue;
		initial_total_time = initial_total_time_save;
		cxxKinetics *kinetics_ptr = NULL;
/*
 *   Run reaction step
 */
		/*
		*   Find maximum number of steps
		*/
		dup_print("Beginning of batch-reaction calculations.", TRUE);
		count_steps = 1;
		if (cxxReaction *rxn_ptr = Utilities::Rxn_find(Rxn_reaction_map, i))
		{
			int count = rxn_ptr->Get_reaction_steps();
			if (count > count_steps)
				count_steps = count;
		}
		if (cxxKinetics *rxn_ptr = Utilities::Rxn_find(Rxn_kinetics_map, i))
		{
			kinetics_ptr = rxn_ptr;
			if (rxn_ptr->Get_reaction_steps() > count_steps)
				count_steps = rxn_ptr->Get_reaction_steps();
		}
		if (cxxTemperature *rxn_ptr = Utilities::Rxn_find(Rxn_temperature_map, i))
		{
			int count = rxn_ptr->Get_countTemps();
			if (count > count_steps)
			{
				count_steps = count;
			}
		}
		if (cxxPressure *rxn_ptr = Utilities::Rxn_find(Rxn_pressure_map, i))
		{
			int count = rxn_ptr->Get_count();
			if (count > count_steps)
			{
				count_steps = count;
			}
		}
		count_total_steps = count_steps;
		if (count_steps > 1)
		{
			state = ADVECTION;
			set_advection(i, TRUE, TRUE, i);
			/*
			*  save data for saving solutions
			*/
			memcpy(&save_data, &save, sizeof(struct save));
			/* 
			*Copy everything to -2
			*/
			copy_use(-2);
			rate_sim_time_start = 0;
			rate_sim_time = 0;
			for (reaction_step = 1; reaction_step <= count_steps; reaction_step++)
			{
				sprintf(token, "Reaction step %d.", reaction_step);
				if (reaction_step > 1 && incremental_reactions == FALSE)
				{
					copy_use(-2);
				}
				set_initial_moles(-2);
				dup_print(token, FALSE);
				/*
				*  Determine time step for kinetics
				*/
				kin_time = 0.0;
				if (use.Get_kinetics_in() == TRUE)
				{
					// runner kin_time
					// equivalent to kin_time in count_steps
					if (run_info.Get_time_step() != NA)
					{
						if (incremental_reactions == FALSE)
						{
							/* not incremental reactions */
							kin_time = reaction_step * run_info.Get_time_step() / ((LDBLE) count_steps);
						}
						else
						{
							/* incremental reactions */
							kin_time = run_info.Get_time_step() / ((LDBLE) count_steps);
						}
					}
					// runner kin_time not defined
					else
					{
						cxxKinetics *kinetics_ptr = Utilities::Rxn_find(Rxn_kinetics_map, -2);
						kin_time = kinetics_ptr->Current_step((incremental_reactions==TRUE), reaction_step);
					}
				}
				if (incremental_reactions == FALSE ||
					(incremental_reactions == TRUE && reaction_step == 1))
				{
					use_mix = TRUE;
				}
				else
				{
					use_mix = FALSE;
				}
				/*
				*   Run reaction step
				*/
				run_reactions(-2, kin_time, use_mix, 1.0);
				if (incremental_reactions == TRUE)
				{
					rate_sim_time_start += kin_time;
					rate_sim_time = rate_sim_time_start;
				}
				else
				{
					rate_sim_time = kin_time;
				}
				punch_all();
				print_all();
				/* saves back into -2 */
				if (reaction_step < count_steps)
				{
					saver();
				}
			}
			/*
			*   save end of reaction
			*/
			memcpy(&save, &save_data, sizeof(struct save));
			if (use.Get_kinetics_in() == TRUE)
			{
				Utilities::Rxn_copy(Rxn_kinetics_map, -2, use.Get_n_kinetics_user());
			}
			saver();
		}
		else
			// only 1 step, no worries about incremental reactions
		{
			state = TRANSPORT;

			rate_sim_time_start = 0;
			rate_sim_time = 0;
			reaction_step = 1;

			sprintf(token, "Reaction step %d.", reaction_step);

			dup_print(token, FALSE);
			/*
			*  Determine time step for kinetics
			*/
			kin_time = 0.0;
			if (kinetics_ptr)
			{
				// runner kin_time
				// equivalent to kin_time in count_steps
				if (run_info.Get_time_step() != NA)
				{
					kin_time = run_info.Get_time_step();
				}
				// runner kin_time not defined
				else
				{
					kin_time = kinetics_ptr->Get_steps()[0];
				}
			}
			/*
			*   Run reaction step
			*/
			use_mix = TRUE;
			run_reactions(i, kin_time, use_mix, 1.0);
			rate_sim_time = kin_time;
			saver();
		}
	}
	initial_total_time += rate_sim_time;
	run_info.Get_cells().Set_defined(false);
	// not running cells
	run_info.Set_run_cells(false);
	return (OK);
}
#else
/* ---------------------------------------------------------------------- */
int Phreeqc::
run_as_cells(void)
/* ---------------------------------------------------------------------- */
{
	struct save save_data;
	LDBLE kin_time;
	int count_steps, use_mix;
	char token[2 * MAX_LENGTH];

	state = REACTION;
	if (run_info.Get_cells().Get_numbers().size() == 0 ||
		!(run_info.Get_cells().Get_defined())) return(OK);

	// running cells
	run_info.Set_run_cells(true);

	dup_print("Beginning of run as cells.", TRUE);
	LDBLE initial_total_time_save;
	if (run_info.Get_start_time() != NA)
	{
		initial_total_time_save = run_info.Get_start_time();
	}
	else
	{
		initial_total_time_save = initial_total_time;
	}

	std::set < int >::iterator it = run_info.Get_cells().Get_numbers().begin();

	for ( ; it != run_info.Get_cells().Get_numbers().end(); it++)
	{
		int i = *it;
		if (i < 0) continue;
		if (Utilities::Rxn_find(Rxn_solution_map, i) == NULL
			&& Utilities::Rxn_find(Rxn_mix_map, i) == NULL)
			continue;
		initial_total_time = initial_total_time_save;
		set_advection(i, TRUE, TRUE, i);
/*
 *   Run reaction step
 */
		/*
		*   Find maximum number of steps
		*/
		dup_print("Beginning of batch-reaction calculations.", TRUE);
		count_steps = 1;
		if (!this->run_cells_one_step)
		{
			if (use.Get_reaction_in() == TRUE && use.Get_reaction_ptr() != NULL)
			{
				int count = use.Get_reaction_ptr()->Get_reaction_steps();
				if (count > count_steps)
					count_steps = count;
			}
			if (use.Get_kinetics_in() == TRUE && use.Get_kinetics_ptr() != NULL)
			{
				if (use.Get_kinetics_ptr()->Get_reaction_steps() > count_steps)
					count_steps = use.Get_kinetics_ptr()->Get_reaction_steps();
			}
			if (use.Get_temperature_in() == TRUE && use.Get_temperature_ptr() != NULL)
			{
				int count = use.Get_temperature_ptr()->Get_countTemps();
				if (count > count_steps)
				{
					count_steps = count;
				}
			}
			if (use.Get_pressure_in() == TRUE && use.Get_pressure_ptr() != NULL)
			{
				int count = use.Get_pressure_ptr()->Get_count();
				if (count > count_steps)
				{
					count_steps = count;
				}
			}
		}
		count_total_steps = count_steps;
		/*
		*  save data for saving solutions
		*/
		memcpy(&save_data, &save, sizeof(struct save));
		/* 
		*Copy everything to -2
		*/
		copy_use(-2);
		rate_sim_time_start = 0;
		rate_sim_time = 0;
		for (reaction_step = 1; reaction_step <= count_steps; reaction_step++)
		{
			sprintf(token, "Reaction step %d.", reaction_step);
			if (reaction_step > 1 && incremental_reactions == FALSE)
			{
				copy_use(-2);
			}
			set_initial_moles(-2);
			dup_print(token, FALSE);
			/*
			*  Determine time step for kinetics
			*/
			kin_time = 0.0;
			if (use.Get_kinetics_in() == TRUE)
			{
				// runner kin_time
				// equivalent to kin_time in count_steps
				if (run_info.Get_time_step() != NA)
				{
					if (incremental_reactions == FALSE)
					{
						/* not incremental reactions */
						kin_time = reaction_step * run_info.Get_time_step() / ((LDBLE) count_steps);
					}
					else
					{
						/* incremental reactions */
						kin_time = run_info.Get_time_step() / ((LDBLE) count_steps);
					}
				}
				// runner kin_time not defined
				else
				{
					cxxKinetics *kinetics_ptr = Utilities::Rxn_find(Rxn_kinetics_map, -2);
					kin_time = kinetics_ptr->Current_step((incremental_reactions==TRUE), reaction_step);
				}
			}
			if (incremental_reactions == FALSE ||
				(incremental_reactions == TRUE && reaction_step == 1))
			{
				use_mix = TRUE;
			}
			else
			{
				use_mix = FALSE;
			}
			/*
			*   Run reaction step
			*/
			run_reactions(-2, kin_time, use_mix, 1.0);
			if (incremental_reactions == TRUE)
			{
				rate_sim_time_start += kin_time;
				rate_sim_time = rate_sim_time_start;
			}
			else
			{
				rate_sim_time = kin_time;
			}
			if (state != ADVECTION)
			{
				punch_all();
				print_all();
			}
			/* saves back into -2 */
			if (reaction_step < count_steps)
			{
				saver();
			}
		}
		/*
		*   save end of reaction
		*/
		memcpy(&save, &save_data, sizeof(struct save));
		if (use.Get_kinetics_in() == TRUE)
		{
			Utilities::Rxn_copy(Rxn_kinetics_map, -2, use.Get_n_kinetics_user());
		}
		saver();
	}
	initial_total_time += rate_sim_time;
	run_info.Get_cells().Set_defined(false);
	// not running cells
	run_info.Set_run_cells(false);
	return (OK);
}
#endif
/* ---------------------------------------------------------------------- */
void Phreeqc::
dump_ostream(std::ostream& os)
/* ---------------------------------------------------------------------- */
{
	// solutions
	if (dump_info.Get_bool_solution())
	{
		if (dump_info.Get_solution().size() == 0)
		{
			Utilities::Rxn_dump_raw(Rxn_solution_map, os, 0);
		}
		else
		{
			std::set < int >::iterator it;
			for (it = dump_info.Get_solution().begin(); it != dump_info.Get_solution().end(); it++)
			{
				cxxSolution *p = Utilities::Rxn_find(Rxn_solution_map, *it);
				if (p != NULL && p->Get_n_user() >= 0)
				{
					p->dump_raw(os, 0);
				}
			}
		}
	}

	// pp_assemblages
	if (dump_info.Get_bool_pp_assemblage())
	{
		if (dump_info.Get_pp_assemblage().size() == 0)
		{
			Utilities::Rxn_dump_raw(Rxn_pp_assemblage_map, os, 0);
		}
		else
		{
			std::set < int >::iterator it;
			for (it = dump_info.Get_pp_assemblage().begin(); it != dump_info.Get_pp_assemblage().end(); it++)
			{
				cxxPPassemblage *p = Utilities::Rxn_find(Rxn_pp_assemblage_map, *it);

				if (p != NULL && p->Get_n_user() >= 0)
				{
					p->dump_raw(os, 0);
				}
			}
		}
	}
	// exchanges
	if (dump_info.Get_bool_exchange())
	{
		if (dump_info.Get_exchange().size() == 0)
		{
			Utilities::Rxn_dump_raw(Rxn_exchange_map, os, 0);
		}
		else
		{
			std::set < int >::iterator it;
			for (it = dump_info.Get_exchange().begin(); it != dump_info.Get_exchange().end(); it++)
			{
				cxxExchange *p = Utilities::Rxn_find(Rxn_exchange_map, *it);

				if (p != NULL && p->Get_n_user() >= 0)
				{
					p->dump_raw(os, 0);
				}
			}
		}
	}

	// surfaces
	if (dump_info.Get_bool_surface())
	{
		if (dump_info.Get_surface().size() == 0)
		{
			Utilities::Rxn_dump_raw(Rxn_surface_map, os, 0);
		}
		else
		{
			std::set < int >::iterator it;
			for (it = dump_info.Get_surface().begin(); it != dump_info.Get_surface().end(); it++)
			{
				cxxSurface *p = Utilities::Rxn_find(Rxn_surface_map, *it);

				if (p != NULL && p->Get_n_user() >= 0)
				{
					p->dump_raw(os, 0);
				}
			}
		}
	}
	// ss_assemblages
	if (dump_info.Get_bool_ss_assemblage())
	{
		if (dump_info.Get_ss_assemblage().size() == 0)
		{
			Utilities::Rxn_dump_raw(Rxn_ss_assemblage_map, os, 0);
		}
		else
		{
			std::set < int >::iterator it;
			for (it = dump_info.Get_ss_assemblage().begin(); it != dump_info.Get_ss_assemblage().end(); it++)
			{
				cxxSSassemblage *p = Utilities::Rxn_find(Rxn_ss_assemblage_map, *it);

				if (p != NULL && p->Get_n_user() >= 0)
				{
					p->dump_raw(os, 0);
				}
			}
		}
	}
	// gas_phases
	if (dump_info.Get_bool_gas_phase())
	{
		if (dump_info.Get_gas_phase().size() == 0)
		{
			Utilities::Rxn_dump_raw(Rxn_gas_phase_map, os, 0);
		}
		else
		{
			std::set < int >::iterator it;
			for (it = dump_info.Get_gas_phase().begin(); it != dump_info.Get_gas_phase().end(); it++)
			{
				cxxGasPhase *p = Utilities::Rxn_find(Rxn_gas_phase_map, *it);

				if (p != NULL && p->Get_n_user() >= 0)
				{
					p->dump_raw(os, 0);
				}
			}
		}
	}

	// kinetics
	if (dump_info.Get_bool_kinetics())
	{
		if (dump_info.Get_kinetics().size() == 0)
		{
			Utilities::Rxn_dump_raw(Rxn_kinetics_map, os, 0);
		}
		else
		{
			std::set < int >::iterator it;
			for (it = dump_info.Get_kinetics().begin(); it != dump_info.Get_kinetics().end(); it++)
			{
				cxxKinetics *p = Utilities::Rxn_find(Rxn_kinetics_map, *it);

				if (p != NULL && p->Get_n_user() >= 0)
				{
					p->dump_raw(os, 0);
				}
			}
		}
	}
	// mix
	if (dump_info.Get_bool_mix())
	{
		if (dump_info.Get_mix().size() == 0)
		{
			Utilities::Rxn_dump_raw(Rxn_mix_map, os, 0);
		}
		else
		{
			std::set < int >::iterator it;
			for (it = dump_info.Get_mix().begin(); it != dump_info.Get_mix().end(); it++)
			{
				cxxMix *p = Utilities::Rxn_find(Rxn_mix_map, *it);

				if (p != NULL && p->Get_n_user() >= 0)
				{
					p->dump_raw(os, 0);
				}
			}
		}
	}

	// reaction
	if (dump_info.Get_bool_reaction())
	{
		if (dump_info.Get_reaction().size() == 0)
		{
			Utilities::Rxn_dump_raw(Rxn_reaction_map, os, 0);
		}
		else
		{
			std::set < int >::iterator it;
			for (it = dump_info.Get_reaction().begin(); it != dump_info.Get_reaction().end(); it++)
			{
				cxxReaction *p = Utilities::Rxn_find(Rxn_reaction_map, *it);

				if (p != NULL && p->Get_n_user() >= 0)
				{
					p->dump_raw(os, 0);
				}
			}
		}
	}

	// temperature
	if (dump_info.Get_bool_temperature())
	{
		if (dump_info.Get_temperature().size() == 0)
		{
			Utilities::Rxn_dump_raw(Rxn_temperature_map, os, 0);
		}
		else
		{
			std::set < int >::iterator it;
			for (it = dump_info.Get_temperature().begin(); it != dump_info.Get_temperature().end(); it++)
			{
				cxxTemperature *p = Utilities::Rxn_find(Rxn_temperature_map, *it);

				if (p != NULL && p->Get_n_user() >= 0)
				{
					p->dump_raw(os, 0);
				}
			}
		}
	}
	// pressure
	if (dump_info.Get_bool_pressure())
	{
		if (dump_info.Get_pressure().size() == 0)
		{
			Utilities::Rxn_dump_raw(Rxn_pressure_map, os, 0);
		}
		else
		{
			std::set < int >::iterator it;
			for (it = dump_info.Get_pressure().begin(); it != dump_info.Get_pressure().end(); it++)
			{
				cxxPressure *p = Utilities::Rxn_find(Rxn_pressure_map, *it);

				if (p != NULL && p->Get_n_user() >= 0)
				{
					p->dump_raw(os, 0);
				}
			}
		}
	}
	// Turn off any reaction calculation
	os << "USE mix none" << "\n";
	os << "USE reaction none" << "\n";
	os << "USE reaction_temperature none" << "\n";
	os << "USE reaction_pressure none" << "\n";

	// Turn off dump until next read
	dump_info.SetAll(false);
}
#if defined MULTICHART
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_user_graph_handler(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads USER_GRAPH_DATA_BLOCK data block
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
	int return_value;

	/*
	 *  Make parser
	 */
	std::istringstream iss_in;
	return_value = streamify_to_next_keyword(iss_in);
	CParser parser(iss_in, phrq_io);

	//For testing, need to read line to get started
	std::vector < std::string > vopts;
	std::istream::pos_type next_char;

	if (pr.echo_input == FALSE)
	{
		parser.set_echo_file(CParser::EO_NONE);
	}
	else
	{
		parser.set_echo_file(CParser::EO_NOKEYWORDS);
	}

	assert(!reading_database());

	bool success = chart_handler.Read(this, parser);

	// Need to output the next keyword
	if (return_value == OPTION_KEYWORD) echo_msg(sformatf( "\t%s\n", line));
	return (return_value);
}
#endif
