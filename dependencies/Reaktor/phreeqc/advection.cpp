#include "Utils.h"
#include "Phreeqc.h"
#include "phqalloc.h"
#include "cxxKinetics.h"
#include "Solution.h"


/* ---------------------------------------------------------------------- */
int Phreeqc::
advection(void)
/* ---------------------------------------------------------------------- */
{
	int i;
	LDBLE kin_time;
/*
 *   Calculate advection
 */
	state = ADVECTION;
/*	mass_water_switch = TRUE; */
/*
 *   Check existence of all solutions
 */
	for (i = 0; i <= count_ad_cells; i++)
	{
		if (Utilities::Rxn_find(Rxn_solution_map, i) == NULL)
		//if (solution_bsearch(i, &n, TRUE) == NULL)
		{
			input_error++;
			error_string = sformatf(
					"Solution %d is needed for advection, but is not defined.",
					i);
			error_msg(error_string, CONTINUE);
		}
	}
/*
 *   Check kinetics logic
 */
	kin_time = advection_kin_time;
	if (kin_time <= 0.0)
	{
		for (i = 1; i <= count_ad_cells; i++)
		{
			if (Utilities::Rxn_find(Rxn_kinetics_map, i) != NULL)
			{
				input_error++;
				error_string = sformatf(
						"KINETIC reaction(s) defined, but time_step is not defined in ADVECTION keyword.");
				error_msg(error_string, CONTINUE);
				break;
			}
		}
	}
/*
 *   Quit on error
 */
	if (get_input_errors() > 0)
	{
		error_msg("Program terminating due to input errors.", STOP);
	}
/*
 *   Equilibrate solutions with phases, exchangers, surfaces
 */
	last_model.force_prep = TRUE;
	rate_sim_time_start = 0;
	for (advection_step = 1; advection_step <= count_ad_shifts;
		 advection_step++)
	{
		log_msg(sformatf(
				   "\nBeginning of advection time step %d, cumulative pore volumes %f.\n",
				   advection_step,
				   (double) (((LDBLE) advection_step) /
							 ((LDBLE) count_ad_cells))));
		if (pr.use == TRUE && pr.all == TRUE)
		{
			output_msg(sformatf(
					   "Beginning of advection time step %d, cumulative pore volumes %f.\n",
					   advection_step,
					   (double) (((LDBLE) advection_step) /
								 ((LDBLE) count_ad_cells))));
		}
/*
 *  Advect
 */
		for (i = count_ad_cells; i > 0; i--)
		{
			//solution_duplicate(i - 1, i);
			Utilities::Rxn_copy(Rxn_solution_map, i -1, i);
		}
/*
 *  Equilibrate and (or) mix
 */
		for (i = 1; i <= count_ad_cells; i++)
		{
			set_initial_moles(i);
			cell_no = i;
			set_advection(i, TRUE, TRUE, i);
			run_reactions(i, kin_time, TRUE, 1.0);
			if (advection_kin_time_defined == TRUE)
			{
				rate_sim_time = rate_sim_time_start + kin_time;
			}
			log_msg(sformatf( "\nCell %d.\n\n", i));
			if (pr.use == TRUE && pr.all == TRUE &&
				advection_step % print_ad_modulus == 0 &&
				advection_print[i - 1] == TRUE)
			{
				output_msg(sformatf( "\nCell %d.\n\n", i));
			}
			if (advection_step % punch_ad_modulus == 0 &&
				advection_punch[i - 1] == TRUE)
			{
				punch_all();
			}
			if (advection_step % print_ad_modulus == 0 &&
				advection_print[i - 1] == TRUE)
			{
				print_all();
			}
			if (i > 1)
				Utilities::Rxn_copy(Rxn_solution_map, -2, i - 1);
				//solution_duplicate(-2, i - 1);
			saver();
		}
		Utilities::Rxn_copy(Rxn_solution_map, -2, count_ad_cells);
		//solution_duplicate(-2, count_ad_cells);
		rate_sim_time_start += kin_time;
	}
	initial_total_time += rate_sim_time_start;
	/* free_model_allocs(); */
	mass_water_switch = FALSE;
	return (OK);
}
