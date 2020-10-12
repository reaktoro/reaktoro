#include "Utils.h"
#include "Phreeqc.h"
#include "phqalloc.h"
#include "Exchange.h"
#include "GasPhase.h"
#include "PPassemblage.h"
#include "SSassemblage.h"
#include "cxxKinetics.h"
#include "Solution.h"
/* ---------------------------------------------------------------------- */
int Phreeqc::
transport(void)
/* ---------------------------------------------------------------------- */
{
	int i, j, k, n;
	int j_imm;
	LDBLE b, f, mix_f_m, mix_f_imm;
	LDBLE water_m, water_imm;
	int first_c, last_c, b_c;
	int max_iter;
	char token[MAX_LENGTH];
	LDBLE kin_time, stagkin_time, kin_time_save;

	int punch_boolean;
	LDBLE step_fraction;

	state = TRANSPORT;
	diffc_tr = diffc;
	diffc_max = 0.0;
	transp_surf = warn_fixed_Surf = warn_MCD_X = 0;
/*	mass_water_switch = TRUE; */
/*
 *   Check existence of solutions
 */
	j = -1;
	/* check column solutions */
	for (i = 1; i <= count_cells; i++)
	{
		use.Set_solution_ptr(Utilities::Rxn_find(Rxn_solution_map, i));
		if (use.Get_solution_ptr() == NULL)
		{
			input_error++;
			error_string = sformatf(
					"Solution %d is needed for transport, but is not defined.",
					i);
			error_msg(error_string, CONTINUE);
		}
		else
			cell_data[i - 1].temp = use.Get_solution_ptr()->Get_tc();
	}

	if (multi_Dflag)
	{
		sol_D = (struct sol_D *) PHRQ_malloc((size_t)
			 (count_cells + 2 + stag_data->count_stag * count_cells) *
			 sizeof(struct sol_D));
		if (sol_D == NULL)
			malloc_error();
		sol_D_dbg = sol_D;
		for (i = 0; i < count_cells + 2 + stag_data->count_stag * count_cells;
			 i++)
		{
			sol_D[i].count_spec = 0;
			sol_D[i].count_exch_spec = 0;
			sol_D[i].exch_total = 0;
			sol_D[i].x_max = 0;
			sol_D[i].spec = NULL;
		}
	}
	/* check solution 0 */
	if (Utilities::Rxn_find(Rxn_solution_map, 0) == NULL)
	{
		if (ishift == 1)
		{
			input_error++;
			error_string = sformatf(
					"Solution 0 is needed for transport, but is not defined.");
			error_msg(error_string, CONTINUE);
		}
		else
			Utilities::Rxn_copy(Rxn_solution_map, 1, 0);
	}

	/* check solution count_cells */
	if (Utilities::Rxn_find(Rxn_solution_map, count_cells + 1) == NULL)
	{
		if (ishift == -1)
		{
			input_error++;
			error_string = sformatf(
					"Solution %d is needed for transport, but is not defined.",
					count_cells + 1);
			error_msg(error_string, CONTINUE);
		}
		else
			Utilities::Rxn_copy(Rxn_solution_map, count_cells, count_cells + 1);
	}
/*
 *   Initialize temperature in stagnant cells ...
 */
	for (n = 1; n <= stag_data->count_stag; n++)
	{
		for (i = 1; i <= count_cells; i++)
		{
			k = i + 1 + n * count_cells;
			use.Set_solution_ptr(Utilities::Rxn_find(Rxn_solution_map, k));
			if (use.Get_solution_ptr() != NULL)
				cell_data[k - 1].temp = use.Get_solution_ptr()->Get_tc();
		}
	}
/*
 * First equilibrate solutions
 */
	dup_print("Equilibrating initial solutions", TRUE);
	transport_step = 0;
	for (i = 0; i <= count_cells + 1; i++)
	{
		if ((bcon_first == 2 && i == 0) ||
			(bcon_last == 2 && i == count_cells + 1))
			continue;
		set_initial_moles(i);
		cell_no = i;
		set_and_run_wrapper(i, NOMIX, FALSE, i, 0.0);
		if (use.Get_surface_ptr() != NULL && use.Get_surface_ptr()->Get_transport())
			transp_surf = TRUE;
		if (transp_surf && !multi_Dflag)
		{
			error_string = sformatf(
					"-multi_d must be defined for surface transport");
			error_msg(error_string, CONTINUE);
		}
		if (multi_Dflag == TRUE)
		{
			fill_spec(cell_no);
		}
		if (cell_no > 0 && cell_no <= count_cells)
		{
			if ((cell_data[i - 1].punch == TRUE)
				&& (cell_no != count_cells + 1))
				punch_all();
			if ((cell_data[i - 1].print == TRUE)
				&& (cell_no != count_cells + 1))
				print_all();
		}
/*    if (i > 0 && i <= count_cells)*/
		saver();
	}
/*
 * Also stagnant cells
 */
	for (n = 1; n <= stag_data->count_stag; n++)
	{
		for (i = 1; i <= count_cells; i++)
		{
			k = i + 1 + n * count_cells;
			cell_no = k;
			if (Utilities::Rxn_find(Rxn_solution_map, k) != 0)
			{
				set_initial_moles(k);
				set_and_run_wrapper(k, NOMIX, FALSE, k, 0.0);
				if (multi_Dflag == TRUE)
				{
					fill_spec(cell_no);
				}
				if (cell_data[k - 1].punch == TRUE)
					punch_all();
				if ((cell_data[k - 1].print == TRUE)
					&& (transport_step % print_modulus == 0))
					print_all();
				saver();
			}
		}
	}
/*
 *  Initialize mixing factors, define kinetics times
 *  for multicomponent diffusion, limit mixing by diffc_max (usually from H+)
 */
	if (multi_Dflag == TRUE)
		diffc_tr = diffc_max;
	if ((stag_data->exch_f > 0) && (stag_data->count_stag == 1))
	{
		/* multi_D calc's are OK if all cells have the same amount of water */
/* if (multi_Dflag == TRUE)
   {
     sprintf(token, "multi_D calc's and stagnant: define MIXing factors explicitly, or \n\t give all cells the same amount of water.");
    	warning_msg(token);
   }
 */
		
		Rxn_mix_map.clear();

/*
 * stagnant mix factors go in mix[0 .. count_cells]
 */

	}
/*
 * mix[] is extended in init_mix(), to accommodate column mix factors
 */
	nmix = init_mix();
	heat_nmix = init_heat_mix(nmix);
	if (nmix < 2)
		stagkin_time = timest;
	else
		stagkin_time = timest / nmix;
	if (ishift != 0)
		kin_time = timest / (1 + nmix);
	else
		kin_time = stagkin_time;
	kin_time_save = kin_time;

/* Reaction defined for a shift... */
	//if (ishift == 0 && nmix == 1)
	//	step_fraction = 1.0;
	//else
	//	step_fraction = 1.0 / (1.0 + nmix);
	if (!ishift)
	{
		if (nmix < 2)
			step_fraction = 1.0;
		else
			step_fraction = 1.0 / nmix;
	}
	else
		step_fraction = 1.0 / (1.0 + nmix);
/*
 *   Set boundary conditions, transport direction
 */
	last_model.force_prep = TRUE;
	if ((ishift == 0) || (bcon_first == 1) || (bcon_last == 1))
		b_c = 1;
	else
		b_c = 0;
	if (ishift >= 0)
	{
		last_c = count_cells;
		first_c = 1;
	}
	else
	{
		last_c = 1;
		first_c = count_cells;
	}
/*
 * Define stagnant/mobile mix structure, if not read explicitly.
 *
 * With count_stag = 1, mix factors are calculated from exchange factor à
 * (= exch_f), mobile é_m (= th_m) and immobile é_im (= th_im) porosity.
 * These variables are read under keyword TRANSPORT, after stagnant, in
 * structure stag_data.
 * MIX 'cell_no' in input file can be an alternative for the calculation here.
 */
	if ((stag_data->exch_f > 0) && (stag_data->count_stag == 1))
	{
		b = stag_data->th_m / (stag_data->th_m + stag_data->th_im);
		f = exp(-stag_data->exch_f * stagkin_time / (b * stag_data->th_im));
		mix_f_imm = b - b * f;
		mix_f_m = mix_f_imm * stag_data->th_im / stag_data->th_m;
		n = 0;
		for (i = 1; i <= count_cells; i++)
		{
			j = i;
			j_imm = j + (1 + count_cells);
			if (Utilities::Rxn_find(Rxn_solution_map, j) == NULL)
				error_msg
					("Could not find mobile cell solution in TRANSPORT.",
					 STOP);
			if (Utilities::Rxn_find(Rxn_solution_map, j_imm) == NULL)
				error_msg
					("Could not find immobile cell solution in TRANSPORT.",
					 STOP);
			water_m = Utilities::Rxn_find(Rxn_solution_map, j)->Get_mass_water();
			water_imm = Utilities::Rxn_find(Rxn_solution_map, j_imm)->Get_mass_water();
/*
 * Define C_m = (1 - mix_f_m) * C_m0  +  mix_f_m) * C_im0
 */
			{
				cxxMix temp_mix;
				temp_mix.Set_n_user(j);
				temp_mix.Set_n_user_end(j);
				temp_mix.Add(j, 1 - mix_f_m);
				temp_mix.Add(j_imm, mix_f_m * water_m / water_imm);
				Rxn_mix_map[j] = temp_mix;
				n++;
			}
/*
 * Define C_im = mix_f_imm * C_m0  +  (1 - mix_f_imm) * C_im0,  or...
 */
			{
				cxxMix temp_mix;
				temp_mix.Set_n_user(j_imm);
				temp_mix.Set_n_user_end(j_imm);
				temp_mix.Add(j_imm, 1 - mix_f_imm);
				temp_mix.Add(j, mix_f_imm * water_imm / water_m);
				Rxn_mix_map[j_imm] = temp_mix;
				n++;
			}
		}

		if (heat_nmix > 0)
		{
/*
 * Assumption: D_e used for calculating exch_f in input file equals diffc
 */
			f = stag_data->exch_f * (heat_diffc - diffc) / diffc / tempr;
			f = exp(-f * stagkin_time / (b * stag_data->th_im));
			heat_mix_f_imm = b - b * f;
			heat_mix_f_m =
				heat_mix_f_imm * stag_data->th_im / stag_data->th_m;
		}
	}
/*
 *   Stop if error
 */
	if (get_input_errors() > 0)
	{
		error_msg("Program terminating due to input errors.", STOP);
	}
/*
 * Now transport
 */
	max_iter = 0;
	for (transport_step = transport_start; transport_step <= count_shifts;
		 transport_step++)
	{
		/*
		 *  Set initial moles of phases
		 */
		for (i = 1; i <= count_cells; i++)
				set_initial_moles(i);
		/*
		 * Also stagnant cells
		 */
		for (n = 1; n <= stag_data->count_stag; n++)
		{
			for (i = 1; i <= count_cells; i++)
			{
				k = i + 1 + n * count_cells;
				cell_no = k;
				if (Utilities::Rxn_find(Rxn_solution_map, k) != 0)
				{
					set_initial_moles(k);
				}
			}
		}
/*
 * Start diffusing if boundary cond = 1, (fixed c, or closed)
 */
		if (b_c == 1)
		{
			/* For half of mixing steps */
			for (j = 1; j <= floor((LDBLE) nmix / 2); j++)
			{
				rate_sim_time_start =
					(transport_step - 1) * timest + (j - 1) * kin_time;
				rate_sim_time = rate_sim_time_start + kin_time;

				mixrun = j;
				if (multi_Dflag)
					sprintf(token,
							"Transport step %3d. Multicomponent diffusion run %3d.",
							transport_step, j);
				else
					sprintf(token, "Transport step %3d. Mixrun %3d.",
							transport_step, j);
				dup_print(token, FALSE);

				if (heat_nmix > 0)
				{
					heat_mix(heat_nmix);
					/* equilibrate again ... */
					for (i = 1; i <= count_cells; i++)
					{
						cell_no = i;
						set_and_run_wrapper(i, NOMIX, FALSE, i, 0.0);
						if (multi_Dflag)
							fill_spec(i);
						saver();
					}
				}
				/* Go through cells */
				if (transp_surf)
				{
					if (disp_surf(stagkin_time) == ERROR)
						error_msg("Error in surface transport, stopping.",
								  STOP);
				}
				if (multi_Dflag)
					multi_D(stagkin_time, 1, FALSE);

				for (i = 1; i <= count_cells; i++)
				{

					if (iterations > max_iter)
						max_iter = iterations;
					cell_no = i;
					mixrun = j;
					if (multi_Dflag)
						sprintf(token,
								"Transport step %3d. MCDrun %3d. Cell %3d. (Max. iter %3d)",
								transport_step, j, i, max_iter);
					else
						sprintf(token,
								"Transport step %3d. Mixrun %3d. Cell %3d. (Max. iter %3d)",
								transport_step, j, i, max_iter);
					status(0, token);

					run_reactions(i, kin_time, DISP, step_fraction);
					if (multi_Dflag)
						fill_spec(i);

					/* punch and output file */
					if ((ishift == 0) && (j == nmix)
						&& ((stag_data->count_stag == 0)
							|| Utilities::Rxn_find(Rxn_solution_map, i + 1 + count_cells) == 0))
					{
						if ((cell_data[i - 1].punch == TRUE)
							&& (transport_step % punch_modulus == 0))
							punch_all();
						if ((cell_data[i - 1].print == TRUE)
							&& (transport_step % print_modulus == 0))
							print_all();
					}
					if (i > 1)
						Utilities::Rxn_copy(Rxn_solution_map, -2, i - 1);
					saver();

					/* maybe sorb a surface component... */
					if ((ishift == 0) && (j == nmix)
						&& ((stag_data->count_stag == 0)
							|| Utilities::Rxn_find(Rxn_solution_map, i + 1 + count_cells) == 0))
					{
						if (change_surf_count > 0)
						{
							for (k = 0; k < change_surf_count; k++)
							{
								if (change_surf[k].cell_no != i)
									break;
								reformat_surf(change_surf[k].comp_name,
											  change_surf[k].fraction,
											  change_surf[k].new_comp_name,
											  change_surf[k].new_Dw,
											  change_surf[k].cell_no);
								change_surf[k].cell_no = -99;
							}
							change_surf_count = 0;
						}
					}
				}
				Utilities::Rxn_copy(Rxn_solution_map, -2, count_cells);

				/* Stagnant zone mixing after completion of each
				   diffusive/dispersive step ...  */
				rate_sim_time_start =
					(transport_step - 1) * timest + (j - 1) * stagkin_time;
				rate_sim_time = rate_sim_time_start + stagkin_time;

				if (stag_data->count_stag > 0)
				{
					if ((ishift == 0) && (j == nmix))
						punch_boolean = TRUE;
					else
						punch_boolean = FALSE;
					for (i = 1; i <= count_cells; i++)
						mix_stag(i, stagkin_time, punch_boolean,
								 step_fraction);
				}
			}
		}
/*
 * Advective transport
 */
		if (ishift != 0)
		{
			sprintf(token, "Transport step %3d.", transport_step);
			dup_print(token, FALSE);
			if (b_c == 1)
				rate_sim_time_start =
					(transport_step - 1) * timest + (j - 1) * kin_time;
			else
				rate_sim_time_start = (transport_step - 1) * timest;
			rate_sim_time = rate_sim_time_start + kin_time;

			/* halftime kinetics for resident water in first cell ... */
			if (Utilities::Rxn_find(Rxn_kinetics_map, first_c) != NULL && count_cells > 1)
			{
				cell_no = first_c;
				kin_time = kin_time_save / 2;
				run_reactions(first_c, kin_time, NOMIX, 0.0);
				saver();
				kin_time = kin_time_save;
			}

			/* for each cell in column */
			for (i = last_c; i != (first_c - ishift); i -= ishift)
				Utilities::Rxn_copy(Rxn_solution_map, i - ishift, i);

/* if boundary_solutions must be flushed by the flux from the column...
      if (ishift == 1 && bcon_last == 3)
	solution_duplicate (last_c, last_c + 1);
      else if (ishift == -1 && bcon_first == 3)
	solution_duplicate (last_c, last_c - 1);
 */
			if (transp_surf)
			{
				for (i = last_c + ishift; i != (first_c - ishift);
					 i -= ishift)
				{
					if ((ishift == 1 && i == last_c + 1 && bcon_last != 3) ||
						(ishift == -1 && i == last_c - 1 && bcon_first != 3))
						continue;
					cxxSurface * surface_ptr = Utilities::Rxn_find(Rxn_surface_map, i - ishift);
					if (surface_ptr == NULL)
					{
						if ((Utilities::Rxn_find(Rxn_surface_map,i) != NULL) &&
							((i == 0 && bcon_first == 3)
							|| (i == count_cells + 1 && bcon_last == 3)))
						{
							Rxn_surface_map.erase(i);
						}
						continue;
					}
					if (surface_ptr->Get_transport())
					{
						cxxSurface * surface_ptr1 = Utilities::Rxn_find(Rxn_surface_map, i);
						if (surface_ptr1 == NULL)
						{
							cxxSurface surf;
							surf.Set_n_user(i);
							surf.Set_n_user_end(i);
							Rxn_surface_map[i] = surf;
						}
						if (i == first_c)
						{
							Rxn_surface_map[i] = mobile_surface_copy(surface_ptr, i, false);
						}
						else
						{
							Rxn_surface_map[i] = mobile_surface_copy(surface_ptr, i, true);
						}
					}
				}
			}

/*
 * thermal diffusion when nmix = 0...
 */
			if ((nmix == 0) && (heat_nmix > 0))
			{
				heat_mix(heat_nmix);
				/* equilibrate again ... */
				for (i = 1; i <= count_cells; i++)
				{
					cell_no = i;
					set_and_run_wrapper(i, NOMIX, FALSE, i, 0.0);
					if (multi_Dflag)
						fill_spec(i);
					saver();
				}
			}

			for (i = 1; i <= count_cells; i++)
			{
				if (i == first_c && count_cells > 1)
					kin_time /= 2;
				cell_no = i;
				mixrun = 0;
				if (multi_Dflag)
					sprintf(token,
							"Transport step %3d. MCDrun %3d. Cell %3d. (Max. iter %3d)",
							transport_step, 0, i, max_iter);
				else
					sprintf(token,
							"Transport step %3d. Mixrun %3d. Cell %3d. (Max. iter %3d)",
							transport_step, 0, i, max_iter);
				status(0, token);
				run_reactions(i, kin_time, NOMIX, step_fraction);
				if (multi_Dflag == TRUE)
					fill_spec(i);
				if (iterations > max_iter)
					max_iter = iterations;
				if ((nmix == 0) && ((stag_data->count_stag == 0) ||
									(Utilities::Rxn_find(Rxn_solution_map, i + 1 + count_cells) == 0)))
				{
					if ((cell_data[i - 1].punch == TRUE)
						&& (transport_step % punch_modulus == 0))
						punch_all();
					if ((cell_data[i - 1].print == TRUE)
						&& (transport_step % print_modulus == 0))
						print_all();
				}
				if (i == first_c && count_cells > 1)
					kin_time = kin_time_save;
				saver();

				/* maybe sorb a surface component... */
				if ((nmix == 0) && ((stag_data->count_stag == 0) ||
									(Utilities::Rxn_find(Rxn_solution_map, i + 1 + count_cells) == 0)))
				{
					if (change_surf_count > 0)
					{
						for (k = 0; k < change_surf_count; k++)
						{
							if (change_surf[k].cell_no != i)
								break;
							reformat_surf(change_surf[k].comp_name,
										  change_surf[k].fraction,
										  change_surf[k].new_comp_name,
										  change_surf[k].new_Dw,
										  change_surf[k].cell_no);
							change_surf[k].cell_no = -99;
						}
						change_surf_count = 0;
					}
				}

				/* If nmix is zero, stagnant zone mixing after
				   advective step ... */
				if ((nmix == 0) && (stag_data->count_stag > 0))
				{
					mix_stag(i, stagkin_time, TRUE, step_fraction);
				}
			}
		}
/*
 * Further dispersive and diffusive transport
 */
		if (b_c != 1)
			j = 1;
		for (j = j; j <= nmix; j++)
		{
			if (multi_Dflag)
				sprintf(token,
						"Transport step %3d. Multicomponent diffusion run %3d.",
						transport_step, j);
			else
				sprintf(token, "Transport step %3d. Mixrun %3d.",
						transport_step, j);
			dup_print(token, FALSE);
			rate_sim_time_start =
				(transport_step - 1) * timest + (j - 1) * kin_time;
			if (ishift != 0)
				rate_sim_time_start += kin_time;
			rate_sim_time = rate_sim_time_start + kin_time;

			if (heat_nmix > 0)
			{
				heat_mix(heat_nmix);
				/* equilibrate again ... */
				for (i = 1; i <= count_cells; i++)
				{
					cell_no = i;
					set_and_run_wrapper(i, NOMIX, FALSE, i, 0.0);
					if (multi_Dflag)
						fill_spec(i);
					saver();
				}
			}
			if (transp_surf)
			{
				if (disp_surf(stagkin_time) == ERROR)
					error_msg("Error in surface transport, stopping.", STOP);
			}

			if (multi_Dflag == TRUE)
				multi_D(stagkin_time, 1, FALSE);

			/* for each cell in column */
			for (i = 1; i <= count_cells; i++)
			{
				if (iterations > max_iter)
					max_iter = iterations;
				cell_no = i;
				mixrun = j;
				if (multi_Dflag)
					sprintf(token,
							"Transport step %3d. MCDrun %3d. Cell %3d. (Max. iter %3d)",
							transport_step, j, i, max_iter);
				else
					sprintf(token,
							"Transport step %3d. Mixrun %3d. Cell %3d. (Max. iter %3d)",
							transport_step, j, i, max_iter);
				status(0, token);

				run_reactions(i, kin_time, DISP, step_fraction);
				if (multi_Dflag == TRUE)
					fill_spec(i);
				if ((j == nmix) && ((stag_data->count_stag == 0) ||
									(Utilities::Rxn_find(Rxn_solution_map, i + 1 + count_cells) == 0)))
				{
					if ((cell_data[i - 1].punch == TRUE)
						&& (transport_step % punch_modulus == 0))
						punch_all();
					if ((cell_data[i - 1].print == TRUE)
						&& (transport_step % print_modulus == 0))
						print_all();
				}
				if (i > 1)
					Utilities::Rxn_copy(Rxn_solution_map, -2, i - 1);
				saver();

				/* maybe sorb a surface component... */
				if ((j == nmix) && ((stag_data->count_stag == 0) ||
									(Utilities::Rxn_find(Rxn_solution_map, i + 1 + count_cells) == 0)))
				{
					if (change_surf_count > 0)
					{
						for (k = 0; k < change_surf_count; k++)
						{
							if (change_surf[k].cell_no != i)
								break;
							reformat_surf(change_surf[k].comp_name,
										  change_surf[k].fraction,
										  change_surf[k].new_comp_name,
										  change_surf[k].new_Dw,
										  change_surf[k].cell_no);
							change_surf[k].cell_no = -99;
						}
						change_surf_count = 0;
					}
				}
			}
			Utilities::Rxn_copy(Rxn_solution_map, -2, count_cells);

			/* Stagnant zone mixing after completion of each
			   diffusive/dispersive step ... */
			rate_sim_time_start =
				(transport_step - 1) * timest + (j - 1) * stagkin_time;
			rate_sim_time = rate_sim_time_start + stagkin_time;

			if (stag_data->count_stag > 0)
			{
				if (j == nmix)
					punch_boolean = TRUE;
				else
					punch_boolean = FALSE;
				for (i = 1; i <= count_cells; i++)
					mix_stag(i, stagkin_time, punch_boolean, step_fraction);
			}
		}
		if (dump_modulus != 0 && (transport_step % dump_modulus) == 0)
			dump();
	}
	screen_msg("\n");

	/* free_model_allocs(); */
/*
 * free mix structures
 */
	Dispersion_mix_map.clear();
	if ((stag_data->exch_f > 0) && (stag_data->count_stag == 1))
	{
		Rxn_mix_map.clear();
	}

	if (heat_nmix > 0)
	{
		heat_mix_array = (LDBLE *) free_check_null(heat_mix_array);
		temp1 = (LDBLE *) free_check_null(temp1);
		temp2 = (LDBLE *) free_check_null(temp2);
	}
	if (multi_Dflag)
	{
		for (i = 0; i < count_cells + 2 + stag_data->count_stag * count_cells;
			 i++)
			sol_D[i].spec = (struct spec *) free_check_null(sol_D[i].spec);
		sol_D = (struct sol_D *) free_check_null(sol_D);
	}

	initial_total_time += rate_sim_time;
	rate_sim_time = 0;
	mass_water_switch = FALSE;
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
init_mix(void)
/* ---------------------------------------------------------------------- */
{
	LDBLE dav, lav, mixf, mf12, maxmix, corr_disp, diffc_here, mD;
	int i, l_nmix;

	std::vector<LDBLE> m, m1;
	for(i = 0; i < count_cells + 1; i++)
	{
		m.push_back(0);
		m1.push_back(0);
	}
	corr_disp = 1.;
	if (correct_disp == TRUE && ishift != 0)
	{
		if (bcon_first == 3)
			corr_disp += 1. / count_cells;
		if (bcon_last == 3)
			corr_disp += 1. / count_cells;
	}
	maxmix = 0.0;
	if (multi_Dflag)
	{
		for (i = 1; i < count_cells; i++)
		{
			lav = (cell_data[i - 1].length + cell_data[i].length) / 2;
			if (ishift != 0)
				dav = (cell_data[i - 1].disp + cell_data[i].disp) / 2;
			else
				dav = 0;
			mixf = dav * corr_disp / cell_data[i].length;
			if (mixf > maxmix)
				maxmix = mixf;
			m[i] = mixf;			/* m[i] has mixf with lower cell */
			mD = diffc_max * timest / (lav * lav);
			if (mD > maxmix)
				maxmix = mD;
		}
/*
 * Also for boundary cells
 */
		if (bcon_first == 1)
		{
			lav = cell_data[0].length;
			if (ishift != 0)
				dav = cell_data[0].disp;
			else
				dav = 0;

			mixf = dav / lav;
			if (mixf > maxmix)
				maxmix = mixf;
			m[0] = 2 * mixf;
			mD = diffc_max * timest / (lav * lav);
			if (mD > maxmix)
				maxmix = mD;
		}
		else
			m[0] = 0;

		if (bcon_last == 1)
		{
			lav = cell_data[count_cells - 1].length;
			if (ishift != 0)
				dav = cell_data[count_cells - 1].disp;
			else
				dav = 0;

			mixf = dav / lav;
			if (mixf > maxmix)
				maxmix = mixf;
			m[count_cells] = 2 * mixf;
			mD = diffc_max * timest / (lav * lav);
			if (mD > maxmix)
				maxmix = mD;
		}
		else
			m[count_cells] = 0;

/*
 * Find number of mixes
 */
		if (maxmix == 0)
		{
			l_nmix = 0;
			if (mcd_substeps > 1 && stag_data->count_stag > 0)
				l_nmix = (int) ceil(mcd_substeps);
		}
		else
		{
			if ((bcon_first == 1) || (bcon_last == 1))
				l_nmix = 1 + (int) floor(4.5 * maxmix);
			else
				l_nmix = 1 + (int) floor(3.0 * maxmix);

			if ((ishift != 0) && ((bcon_first == 1) || (bcon_last == 1)))
			{
				if (l_nmix < 2)
					l_nmix = 2;
			}
			if (mcd_substeps > 1)
				l_nmix = (int) ceil(l_nmix * mcd_substeps);

			for (i = 0; i <= count_cells; i++)
				m[i] /= l_nmix;
		}
/*
 * Fill mix structure
 */
	
		if (l_nmix != 0)
		{
			for (i = 1; i <= count_cells; i++)
			{
				cxxMix temp_mix;
				temp_mix.Set_n_user(i);
				temp_mix.Set_n_user_end(i);

				temp_mix.Add(i - 1, m[i - 1]);
				temp_mix.Add(i + 1, m[i]);
				temp_mix.Add(i, 1.0 - m[i - 1] - m[i]);
				Dispersion_mix_map[i] = temp_mix;
			}
		}
		return (l_nmix);
	}
	else // multi_D false
	{
		diffc_here = diffc_tr;
/*
 * Define mixing factors among inner cells
 */
		for (i = 1; i < count_cells; i++)
		{
// find mix with lower numbered cell...
			lav = (cell_data[i - 1].length + cell_data[i].length) / 2;
			if (ishift != 0)
				dav = (cell_data[i - 1].disp + cell_data[i].disp) / 2;
			else
				dav = 0;

			mixf = (diffc_here * timest / lav + dav) * corr_disp / cell_data[i].length;
			m[i] = mixf;			/* m[i] has mixf with lower cell */

// and with higher numbered cell...
			mixf = (diffc_here * timest / lav + dav) * corr_disp / cell_data[i - 1].length;
			mf12 = m[i] + mixf;
			if (mf12 > maxmix)
				maxmix = mf12;
			m1[i] = mixf;			/* m1[i] has mixf with higher cell */
		}
/*
 * Also for boundary cells
 */
		if (bcon_first == 1)
		{
			lav = cell_data[0].length;
			if (ishift != 0)
				dav = cell_data[0].disp;
			else
				dav = 0;

			mixf = (diffc_here * timest / lav + dav) / lav;
			mf12 = m1[1] + 2 * mixf;
			if (mf12 > maxmix)
				maxmix = mf12;
			m[0] = 2 * mixf;
		}
		else
			m[0] = 0;

		if (bcon_last == 1)
		{
			lav = cell_data[count_cells - 1].length;
			if (ishift != 0)
				dav = cell_data[count_cells - 1].disp;
			else
				dav = 0;

			mixf = (diffc_here * timest / lav + dav) / lav;
			mf12 = m[count_cells - 1] + 2 * mixf;
			if (mf12 > maxmix && !multi_Dflag)
				maxmix = mf12;
			m1[count_cells] = 2 * mixf;
		}
		else
			m[count_cells] = 0;

/*
 * Find number of mixes
 */
		if (maxmix == 0)
		{
			l_nmix = 0;
		}
		else
		{
			l_nmix = 1 + (int) floor(1.5 * maxmix);

			if ((ishift != 0) && ((bcon_first == 1) || (bcon_last == 1)))
			{
				if (l_nmix < 2)
					l_nmix = 2;
			}

			for (i = 0; i <= count_cells; i++)
			{
				m[i] /= l_nmix;
				m1[i] /= l_nmix;
			}
		}
		/*
		 * Fill mix structure
		 */
		
		if (l_nmix != 0)
		{
			for (i = 1; i <= count_cells; i++)
			{
				cxxMix temp_mix;
				temp_mix.Set_n_user(i);
				temp_mix.Set_n_user_end(i);

				temp_mix.Add(i - 1, m[i - 1]);
				temp_mix.Add(i + 1, m1[i]);
				temp_mix.Add(i, 1.0 - m[i - 1] - m1[i]);
				Dispersion_mix_map[i] = temp_mix;
			}
		}
		return (l_nmix);
	}
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
mix_stag(int i, LDBLE kin_time, int l_punch, LDBLE step_fraction)
/* ---------------------------------------------------------------------- */
{
	int j, n, k;
	LDBLE t_imm;
	cxxSolution *ptr_imm, *ptr_m;
/*
 * Kinetics in transport cell is done while transporting
 */
	for (n = 1; n <= stag_data->count_stag; n++)
	{
		k = i + 1 + n * count_cells;
		if ((ptr_imm = Utilities::Rxn_find(Rxn_solution_map, k)) != NULL)
		{
			if (n == 1)
			{
				if (heat_nmix > 0)
				{
					ptr_m = Utilities::Rxn_find(Rxn_solution_map, i);
					t_imm =
						heat_mix_f_imm * ptr_m->Get_tc() + (1 -
													  heat_mix_f_imm) *
						ptr_imm->Get_tc();
					ptr_m->Set_tc(
						heat_mix_f_m * ptr_imm->Get_tc() + (1 -
													  heat_mix_f_m) *
						ptr_m->Get_tc());
					cell_data[i - 1].temp = ptr_m->Get_tc();
					cell_data[k - 1].temp= t_imm = ptr_imm->Get_tc();
					/* equilibrate again ... */
					cell_no = i;
					set_and_run_wrapper(i, NOMIX, FALSE, i, 0.0);
					if (multi_Dflag == TRUE)
						fill_spec(cell_no);
					saver();
					cell_no = k;
					set_and_run_wrapper(k, NOMIX, FALSE, k, 0.0);
					if (multi_Dflag == TRUE)
						fill_spec(cell_no);
					saver();
				}
/*
 * Mobile cell, kinetics already done ...
 */
				cell_no = i;
				if (transp_surf)
				{
					if (diff_stag_surf(i) == ERROR)
						error_msg("Error in surface transport, stopping.",
								  STOP);
				}
				if (multi_Dflag == TRUE)
					multi_D(1.0, i, TRUE);
				set_and_run_wrapper(i, STAG, FALSE, -2, 0.0);
				if (multi_Dflag == TRUE)
					fill_spec(cell_no);
				use.Set_kinetics_ptr(Utilities::Rxn_find(Rxn_kinetics_map, i));
				if (use.Get_kinetics_ptr() != NULL)
				{
					use.Set_n_kinetics_user(i);
					use.Set_kinetics_in(true);
				}

				if (l_punch && (cell_data[i - 1].print == TRUE) &&
					(transport_step % print_modulus == 0))
					print_all();
				if (l_punch && (cell_data[i - 1].punch == TRUE) &&
					(transport_step % punch_modulus == 0))
					punch_all();
				saver();

				/* maybe sorb a surface component... */
				if (l_punch && change_surf_count)
				{
					for (j = 0; j < change_surf_count; j++)
					{
						if (change_surf[j].cell_no != i)
							break;
						reformat_surf(change_surf[j].comp_name,
									  change_surf[j].fraction,
									  change_surf[j].new_comp_name,
									  change_surf[j].new_Dw,
									  change_surf[j].cell_no);
						change_surf[j].cell_no = -99;
					}
					change_surf_count = 0;
				}
			}

			cell_no = k;
			run_reactions(k, kin_time, STAG, step_fraction);
			if (multi_Dflag == TRUE)
				fill_spec(cell_no);

			if ((cell_data[k - 1].print == TRUE) && (l_punch == TRUE) &&
				(transport_step % print_modulus == 0))
				print_all();
			if ((cell_data[k - 1].punch == TRUE) && (l_punch == TRUE) &&
				(transport_step % punch_modulus == 0))
				punch_all();
			saver();

			/* maybe sorb a surface component... */
			if (l_punch && change_surf_count)
			{
				for (j = 0; j < change_surf_count; j++)
				{
					if (change_surf[j].cell_no != k)
						break;
					reformat_surf(change_surf[j].comp_name,
								  change_surf[j].fraction,
								  change_surf[j].new_comp_name,
								  change_surf[j].new_Dw,
								  change_surf[j].cell_no);
					change_surf[j].cell_no = -99;
				}
				change_surf_count = 0;
			}
		}
	}
	for (n = 1; n <= stag_data->count_stag; n++)
	{
		k = i + 1 + n * count_cells;
		if (Utilities::Rxn_find(Rxn_solution_map, k) != 0)
		{
			Utilities::Rxn_copy(Rxn_solution_map, -2 - k, k);
			if (n == 1)
				Utilities::Rxn_copy(Rxn_solution_map, -2, i);
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
init_heat_mix(int l_nmix)
/* ---------------------------------------------------------------------- */
{
	LDBLE lav, mixf, maxmix, corr_disp;
	int i, k, n;
	int l_heat_nmix;
	LDBLE t0;
/*
 * Check for need to model thermal diffusion...
 */
	if (heat_diffc <= diffc)
		return (0);
	if (count_cells < 2)
		return (0);

	l_heat_nmix = 0;
	t0 = Utilities::Rxn_find(Rxn_solution_map, 0)->Get_tc();
	for (i = 0; i < count_cells; i++)
	{
		if (fabs(cell_data[i].temp - t0) > 1.0)
		{
			l_heat_nmix = 1;
			break;
		}
	}
	if (l_heat_nmix == 0)
	{
		if (fabs(Utilities::Rxn_find(Rxn_solution_map, count_cells + 1)->Get_tc() - t0) > 1.0)
			l_heat_nmix = 1;
		for (n = 1; n <= stag_data->count_stag; n++)
		{
			for (i = 1; i < count_cells; i++)
			{
				k = i + 1 + n * count_cells;
				if (Utilities::Rxn_find(Rxn_solution_map, k) != 0)
				{
					if (fabs(cell_data[k - 1].temp - t0) > 1.0)
					{
						l_heat_nmix = 1;
						break;
					}
				}
			}
		}
	}
	if (l_heat_nmix == 0)
		return (0);
/*
 * Initialize arrays...
 */
	heat_mix_array = (LDBLE *) PHRQ_malloc((count_cells + 2) * sizeof(LDBLE));
	if (heat_mix_array == NULL)
		malloc_error();

	temp1 = (LDBLE *) PHRQ_malloc((count_cells + 2) * sizeof(LDBLE));
	if (temp1 == NULL)
		malloc_error();

	temp2 = (LDBLE *) PHRQ_malloc((count_cells + 2) * sizeof(LDBLE));
	if (temp2 == NULL)
		malloc_error();
/*
 * Define mixing factors among inner cells...
 */
	corr_disp = 1.;
	if (correct_disp == TRUE && ishift != 0)
	{
		if (bcon_first == 3)
			corr_disp += 1. / count_cells;
		if (bcon_last == 3)
			corr_disp += 1. / count_cells;
	}
	if (l_nmix > 0)
		corr_disp /= l_nmix;
	maxmix = 0.0;
	for (i = 1; i < count_cells; i++)
	{
		lav = (cell_data[i - 1].length + cell_data[i].length) / 2;
		mixf =
			(heat_diffc -
			 diffc_tr) * timest * corr_disp / tempr / (lav * lav);
		if (mixf > maxmix)
			maxmix = mixf;
		heat_mix_array[i + 1] = mixf;	/* m[i] has mixf with lower cell */
	}
/*
 * Also for boundary cells
 */
	if (bcon_first == 1)
	{
		lav = cell_data[0].length;
		mixf =
			(heat_diffc -
			 diffc_tr) * timest * corr_disp / tempr / (lav * lav);
		if (2 * mixf > maxmix)
			maxmix = 2 * mixf;
		heat_mix_array[1] = 2 * mixf;
	}
	else
		heat_mix_array[1] = 0;

	if (bcon_last == 1)
	{
		lav = cell_data[count_cells - 1].length;
		mixf =
			(heat_diffc -
			 diffc_tr) * timest * corr_disp / tempr / (lav * lav);
		if (2 * mixf > maxmix)
			maxmix = 2 * mixf;
		heat_mix_array[count_cells + 1] = 2 * mixf;
	}
	else
		heat_mix_array[count_cells + 1] = 0;
/*
 * Find number of mixes
 */
	if (maxmix == 0)
		l_heat_nmix = 0;
	else
	{
		l_heat_nmix = 1 + (int) floor(3.0 * maxmix);
		for (i = 1; i <= count_cells + 1; i++)
			heat_mix_array[i] /= l_heat_nmix;
	}

	return (l_heat_nmix);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
heat_mix(int l_heat_nmix)
/* ---------------------------------------------------------------------- */
{
	int i, j;

	for (i = 1; i <= count_cells; i++)
		temp1[i] = Utilities::Rxn_find(Rxn_solution_map, i)->Get_tc();
	temp1[0] = Utilities::Rxn_find(Rxn_solution_map, 0)->Get_tc();
	temp1[count_cells + 1] =
		Utilities::Rxn_find(Rxn_solution_map, (count_cells + 1))->Get_tc();

	for (i = 1; i <= l_heat_nmix; i++)
	{
		for (j = 1; j <= count_cells; j++)
			temp2[j] =
				heat_mix_array[j] * temp1[j - 1] + heat_mix_array[j +
																  1] *
				temp1[j + 1] + (1 - heat_mix_array[j] -
								heat_mix_array[j + 1]) * temp1[j];
		for (j = 1; j <= count_cells; j++)
			temp1[j] = temp2[j];
	}

	for (i = 1; i <= count_cells; i++)
	{
		cell_data[i - 1].temp = temp1[i];
		Utilities::Rxn_find(Rxn_solution_map, i)->Set_tc(temp1[i]);
	}

	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
set_initial_moles(int i)
/* ---------------------------------------------------------------------- */
{
	cxxKinetics *kinetics_ptr;
	char token[MAX_LENGTH], token1[MAX_LENGTH], *ptr;
	int j, k, l;
	/*
	 *   Pure phase assemblage
	 */
	{
		cxxPPassemblage * pp_assemblage_ptr = Utilities::Rxn_find(Rxn_pp_assemblage_map, i);
		if (pp_assemblage_ptr != NULL)
		{
			std::map<std::string, cxxPPassemblageComp>::iterator it;
			it =  pp_assemblage_ptr->Get_pp_assemblage_comps().begin();
			for ( ; it != pp_assemblage_ptr->Get_pp_assemblage_comps().end(); it++)
			{
				it->second.Set_initial_moles(it->second.Get_moles());
				if (it->second.Get_initial_moles() < 0)
					it->second.Set_initial_moles(0.0);
			}
		}
	}
	/*
	 *   Gas phase
	 */
	{
		cxxGasPhase * gas_phase_ptr = Utilities::Rxn_find(Rxn_gas_phase_map, i);
		if (gas_phase_ptr != NULL)
		{
			std::vector<cxxGasComp> gc = gas_phase_ptr->Get_gas_comps();
			for (size_t l = 0; l < gc.size(); l++)
			{
				gc[l].Set_initial_moles(gc[l].Get_moles());
			}
			gas_phase_ptr->Set_gas_comps(gc);
		}
	}
	/*
	 *   Kinetics
	 */
	kinetics_ptr = Utilities::Rxn_find(Rxn_kinetics_map, i);
	if (kinetics_ptr != NULL)
	{
		for (j = 0; j < (int) kinetics_ptr->Get_kinetics_comps().size(); j++)
		{
			cxxKineticsComp *kinetics_comp_ptr = &(kinetics_ptr->Get_kinetics_comps()[j]);
			kinetics_comp_ptr->Set_initial_moles(kinetics_comp_ptr->Get_m());
		}
	}
	/*
	 *   Solid solutions
	 */
	{
		cxxSSassemblage *ss_assemblage_ptr = Utilities::Rxn_find(Rxn_ss_assemblage_map, i);
		if (ss_assemblage_ptr != NULL)
		{
			std::vector<cxxSS *> ss_ptrs = ss_assemblage_ptr->Vectorize();
			for (k = 0; k < (int) ss_ptrs.size(); k++)
			{
				cxxSS * ss_ptr = ss_ptrs[k];
				for (j = 0; j < (int) ss_ptr->Get_ss_comps().size(); j++)
				{
					cxxSScomp *comp_ptr = &(ss_ptr->Get_ss_comps()[j]);
					comp_ptr->Set_init_moles(comp_ptr->Get_moles());
				}
			}
		}
	}
	/*
	 *   For interlayer diffusion: add tiny bit of exchanger if absent
	 */
	cxxExchange * exchange_ptr = Utilities::Rxn_find(Rxn_exchange_map, i);
	if (interlayer_Dflag && exchange_ptr == NULL)
	{
		cxxExchange temp_exchange;
		temp_exchange.Set_n_user_both(i);
		temp_exchange.Set_description("Interlayer diffusion: added 2e-10 moles X-");
		use.Set_exchange_in(true);
		use.Set_n_exchange_user(i);

		temp_exchange.Set_new_def(true);
		temp_exchange.Set_solution_equilibria(true);
		temp_exchange.Set_n_solution(i);

		cxxExchComp comp;
		count_elts = 0;
		paren_count = 0;
		strcpy(token, "X");
		ptr = token;
		get_elts_in_species(&ptr, 2e-10);
		ptr = token;
		LDBLE z;
		get_token(&ptr, token1, &z, &l);
		comp.Set_formula(token1);
		comp.Set_formula_z(z);
		comp.Set_totals(elt_list_NameDouble());
		comp.Set_charge_balance(0.0);
		temp_exchange.Get_exchange_comps().push_back(comp);
		Rxn_exchange_map[i] = temp_exchange;

		state = INITIAL_EXCHANGE;
		initial_exchangers(TRUE);
		state = TRANSPORT;
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
fill_spec(int l_cell_no)
/* ---------------------------------------------------------------------- */
{
/* copy species activities into sol_D.spec... */

	int i, i2, count_spec, count_exch_spec;
	char token[MAX_LENGTH];
	const char * name;
	struct species *s_ptr, *s_ptr2;
	struct master *master_ptr;
	LDBLE dum, dum2;
	LDBLE lm;
	LDBLE por, por_il, temp_factor, temp_il_factor, viscos;
	bool x_max_done = false;

	s_ptr2 = NULL;

	sol_D[l_cell_no].spec =
		(struct spec *) free_check_null(sol_D[l_cell_no].spec);
	sol_D[l_cell_no].spec =
		(struct spec *) PHRQ_malloc((size_t) count_species_list *
									sizeof(struct spec));
	if (sol_D[l_cell_no].spec == NULL)
		malloc_error();

	temp_factor = temp_il_factor = 1.0;
	if (l_cell_no == 0)
	{
		por = cell_data[0].por;
		por_il = cell_data[0].por_il;
	}
	else if (l_cell_no == count_cells + 1)
	{
		por = cell_data[count_cells - 1].por;
		por_il = cell_data[count_cells - 1].por_il;
	}
	else
	{
		por = cell_data[l_cell_no - 1].por;
		por_il = cell_data[l_cell_no - 1].por_il;
	}
	if (por < multi_Dpor_lim)
		por = temp_factor = 0.0;

	if (por_il < interlayer_Dpor_lim)
		por_il = temp_il_factor = 0.0;
/*
 * correct diffusion coefficient for temperature and viscosity, D_T = D_298 * Tk * viscos_298 / (298 * viscos)
 */
	viscos = viscosity();
/*
 * put temperature factor in por_factor which corrects for porous medium...
 */
	temp_factor *= tk_x * 0.88862 / (298.15 * viscos);
	temp_il_factor *= tk_x * 0.88862 / (298.15 * viscos);

	count_spec = count_exch_spec = 0;
/*
 * sort species by name...
 */
	if (count_species_list > 0)
		qsort(&species_list[0], (size_t) count_species_list,
			  (size_t) sizeof(struct species_list), sort_species_name);

	for (i = 0; i < count_species_list; i++)
	{
/*
 *   copy species data
 */
		s_ptr = species_list[i].s;

		if (s_ptr->type == EX && !interlayer_Dflag)
			continue;
		if (s_ptr->type == SURF)
			continue;
		if (i > 0 && strcmp(s_ptr->name, species_list[i - 1].s->name) == 0)
			continue;
		if (s_ptr == s_h2o)
			continue;

		if (s_ptr->type == EX)
		{
			if (s_ptr->moles > 1e-30)
			{
				/* find exchanger's name, use only master exchanger 'X' */
				if (species_list[i].master_s->secondary != NULL)
					master_ptr = species_list[i].master_s->secondary;
				else
					master_ptr = species_list[i].master_s->primary;
				if (s_ptr->equiv != 0.0)
					dum = fabs(s_ptr->equiv) / master_ptr->total;
				else
				{
					if (species_list[i].master_s->z == 0)
						dum = 1 / master_ptr->total;
					else
						dum = 1;
				}
				name = master_ptr->elt->name;
				if (strcmp(name, "X") != 0)
				{
					if (!warn_MCD_X)
					{
						sprintf(token,
							"MCD found more than 1 exchanger, uses X for interlayer diffusion.");
						warning_msg(token);
						warn_MCD_X = 1;
					}
					continue;
				}
				dum2 = s_ptr->moles * dum;	/* equivalent fraction */
				sol_D[l_cell_no].spec[count_spec].name = s_ptr->name;
				  //string_hsave(s_ptr->name);
				sol_D[l_cell_no].spec[count_spec].type = EX;
				sol_D[l_cell_no].spec[count_spec].c = dum2;
				sol_D[l_cell_no].spec[count_spec].lg = s_ptr->lg - log10(dum);
				sol_D[l_cell_no].spec[count_spec].a = dum2 * pow(10,
															   sol_D
															   [l_cell_no].
															   spec
															   [count_spec].
															   lg);
				sol_D[l_cell_no].exch_total = master_ptr->total;
				if (transport_step == 0 && !x_max_done)
				{
					x_max_done = true;
					dum = master_ptr->total / Utilities::Rxn_find(Rxn_solution_map, l_cell_no)->Get_mass_water();
					if (dum > sol_D[1].x_max)
						sol_D[1].x_max = dum;
				}

				/* find the aqueous species in the exchange reaction... */
				for (i2 = 0; (s_ptr->rxn->token[i2].s != NULL); i2++)
				{
					if ((s_ptr2 = s_ptr->rxn->token[i2].s)->type == AQ)
						break;
				}
				/* copy its name and Dw and charge... */
				sol_D[l_cell_no].spec[count_spec].aq_name = s_ptr2->name;
				  //string_hsave(s_ptr2->name);
				sol_D[l_cell_no].spec[count_spec].z = s_ptr2->z;
				if (s_ptr2->dw == 0)
					sol_D[l_cell_no].spec[count_spec].Dwt =
						default_Dw * temp_il_factor;
				else
					sol_D[l_cell_no].spec[count_spec].Dwt =
						s_ptr2->dw * temp_il_factor;
				count_exch_spec++;
				count_spec++;
			}
			continue;
		}

		lm = s_ptr->lm;
		if (lm > MIN_LM)
		{
			//sol_D[l_cell_no].spec[count_spec].name = string_hsave(s_ptr->name);
			sol_D[l_cell_no].spec[count_spec].name = s_ptr->name;
			sol_D[l_cell_no].spec[count_spec].type = AQ;
			sol_D[l_cell_no].spec[count_spec].c =
				s_ptr->moles / mass_water_aq_x;
			sol_D[l_cell_no].spec[count_spec].a = under(lm + s_ptr->lg);
			sol_D[l_cell_no].spec[count_spec].lm = lm;
			sol_D[l_cell_no].spec[count_spec].lg = s_ptr->lg;
			sol_D[l_cell_no].spec[count_spec].z = s_ptr->z;
			if (s_ptr->dw == 0)
				sol_D[l_cell_no].spec[count_spec].Dwt =
					default_Dw * temp_factor;
			else
				sol_D[l_cell_no].spec[count_spec].Dwt = s_ptr->dw * temp_factor;
			if (sol_D[l_cell_no].spec[count_spec].Dwt * pow(por, multi_Dn) >
				diffc_max)
				diffc_max =
					sol_D[l_cell_no].spec[count_spec].Dwt * pow(por, multi_Dn);
			sol_D[l_cell_no].spec[count_spec].erm_ddl = s_ptr->erm_ddl;

			count_spec++;
		}
	}
	sol_D[l_cell_no].spec =
		(struct spec *) PHRQ_realloc(sol_D[l_cell_no].spec,
									 (size_t) count_spec *
									 sizeof(struct spec));
	if (sol_D[l_cell_no].spec == NULL)
		malloc_error();

	sol_D[l_cell_no].count_spec = count_spec;
	sol_D[l_cell_no].count_exch_spec = count_exch_spec;

	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
sort_species_name(const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
	const struct species_list *nptr1, *nptr2;

	nptr1 = (const struct species_list *) ptr1;
	nptr2 = (const struct species_list *) ptr2;

	return (strcmp(nptr1->s->name, nptr2->s->name));
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
multi_D(LDBLE DDt, int mobile_cell, int stagnant)
/* ---------------------------------------------------------------------- */
{
	/*
	 * 1. determine mole transfer (mol/s) of solute species for the interface between 2 cells.
	 * 2. sum up as mole transfer of master_species
	 * 3. add moles of master_species to the 2 cells
	 *      NOTE. Define the water content of stagnant cells relative to the
	 *      mobile cell (with, for example, 1 kg water)
	 *      Define properties of each interface only 1 time with MIX.
	 */
	int icell, jcell, i, l, n, length, length2, il_calcs;
	int i1;
	int first_c, last_c;
	char token[MAX_LENGTH];
	LDBLE mixf, temp;

	for (n = 0; n < (stagnant ? stag_data->count_stag : 1); n++)
	{
		icell = mobile_cell + 1 + n * count_cells;
		if (stagnant)
		{
			if (n == 0)
				icell -= 1;
			/*
			 *    find the mix ptr for icell and go along the cells that mix with it
			 */
			use.Set_mix_ptr(Utilities::Rxn_find(Rxn_mix_map, icell));
			if (use.Get_mix_ptr() == NULL)
				continue;
			first_c = 0;
			last_c = (int) (use.Get_mix_ptr()->Get_mixComps().size() - 1);
		}
		else
		{						/* regular column... */
			if (bcon_first == 1)
				first_c = 0;
			else
				first_c = 1;
			if (bcon_last == 1)
				last_c = count_cells;
			else
				last_c = count_cells - 1;
		}

		for (i = first_c; i <= last_c; i++)
		{
			if (stagnant)
			{
				std::vector<int> n_solution;
				std::vector<LDBLE> fraction;
				(use.Get_mix_ptr())->Vectorize(n_solution, fraction);

				if ((jcell = n_solution[i]) <= icell)
					continue;

				mixf = fraction[i];
				if (mcd_substeps > 1)
					mixf /= nmix;
			}
			else
			{					/* regular column... */
				icell = i;
				jcell = i + 1;
				mixf = 1.0;
			}
			/*
			 * 1. obtain J_ij...
			 */
			il_calcs = find_J(icell, jcell, mixf, DDt, stagnant);
			/*
			 * 2. sum up the primary or secondary master_species
			 */
			if (!il_calcs)
			{
				tot1_h = tot1_o = tot2_h = tot2_o = 0.0;
				m_s = (struct M_S *) free_check_null(m_s);
				m_s = (struct M_S *) PHRQ_malloc((size_t) count_elements *
												 sizeof(struct M_S));
				if (m_s == NULL)
					malloc_error();
				for (i1 = 0; i1 < count_elements; i1++)
				{
					m_s[i1].name = NULL;
					m_s[i1].tot1 = 0;
					m_s[i1].tot2 = 0;
				}
				count_m_s = 0;
			}
			fill_m_s(J_ij, J_ij_count_spec);

			/*
			 * 3. find the solutions, add or subtract the moles...
			 */
			if (i > 0 || stagnant)
			{
				use.Set_solution_ptr(Utilities::Rxn_find(Rxn_solution_map, icell));
				use.Get_solution_ptr()->Set_total_h(use.Get_solution_ptr()->Get_total_h() - tot1_h);
				use.Get_solution_ptr()->Set_total_o(use.Get_solution_ptr()->Get_total_o() - tot1_o);
				use.Get_solution_ptr()->Set_cb(use.Get_solution_ptr()->Get_cb() - J_ij_sum);
				for (l = 0; l < count_m_s; l++)
				{
					temp = 0.0;
					length = (int) strlen(m_s[l].name);
					cxxNameDouble::iterator it;
					for (it = use.Get_solution_ptr()->Get_totals().begin(); 
						it != use.Get_solution_ptr()->Get_totals().end(); it++)
					{
						LDBLE moles = it->second;
						length2 =
							(int) (size_t) strcspn(it->first.c_str(), "(");
						if (strncmp
							(m_s[l].name,
							 it->first.c_str(),
							 length) == 0 && length == length2)
						{
							if (moles <	m_s[l].tot1)
							{
								temp = moles;
								it->second = 0.0;
								/* see if other redox states have more moles... */
								cxxNameDouble::iterator kit = it;
								kit++;
								for ( ; kit != use.Get_solution_ptr()->Get_totals().end(); kit++)
								{
									length2 = (int) (size_t) strcspn(
											kit->first.c_str(), "(");
									if (strncmp(m_s[l].name,
										kit->first.c_str(), length) == 0
										&& length == length2)
									{
										temp += kit->second;
										if (temp < m_s[l].tot1)
										{
											kit->second = 0;
										}
										else
										{
											kit->second = temp - m_s[l].tot1;
											temp = 0.0;
											break;
										}
									}
								}
								if (temp != 0.0 && m_s[l].tot1 - temp > 1e-12)
								{
									sprintf(token,
											"Negative concentration in MCD: added %.1e moles %s in cell %d.",
											(double) (m_s[l].tot1 - temp),
											m_s[l].name, icell);
									warning_msg(token);
								}
							}
							else
								it->second -= m_s[l].tot1;
							break;
						}
					}
					if (it == use.Get_solution_ptr()->Get_totals().end())
					{
						use.Get_solution_ptr()->Get_totals()[m_s[l].name] = -m_s[l].tot1;
						if (-m_s[l].tot1 < 0)
						{
							if (-m_s[l].tot1 < -1e-12)
							{
								sprintf(token,
										"Negative concentration in MCD: added %.2e moles %s in cell %d",
										(double) m_s[l].tot1, m_s[l].name, icell);
								warning_msg(token);
							}
							use.Get_solution_ptr()->Get_totals()[m_s[l].name] = 0;
						}
					}
				}
			}
			if (i < count_cells || stagnant)
			{
				use.Set_solution_ptr(Utilities::Rxn_find(Rxn_solution_map, jcell));
				dummy = use.Get_solution_ptr()->Get_total_h();
				use.Get_solution_ptr()->Set_total_h(dummy + tot2_h);
				dummy = use.Get_solution_ptr()->Get_total_o();
				use.Get_solution_ptr()->Set_total_o(dummy + tot2_o);
				dummy = use.Get_solution_ptr()->Get_cb();
				use.Get_solution_ptr()->Set_cb(dummy + J_ij_sum);
				for (l = 0; l < count_m_s; l++)
				{
					temp = 0.0;
					length = (int) strlen(m_s[l].name);
					cxxNameDouble::iterator it;
					for (it = use.Get_solution_ptr()->Get_totals().begin(); 
						it != use.Get_solution_ptr()->Get_totals().end(); it++)
					{
						length2 = (int) (size_t) strcspn(
							it->first.c_str(), "(");
						if (strncmp(m_s[l].name,
							it->first.c_str(), length) == 0 
							&& length == length2)
						{
							if (it->second < -m_s[l].tot2)
							{
								temp = it->second;
								it->second = 0;
								/* see if other redox states have more moles... */
								cxxNameDouble::iterator kit = it;
								kit++;
								for ( ; kit != use.Get_solution_ptr()->Get_totals().end(); kit++)
								{
									length2 = (int) (size_t) strcspn(
										kit->first.c_str(), "(");
									if (strncmp
										(m_s[l].name,
										kit->first.c_str(), length) == 0
										&& length == length2)
									{
										temp += kit->second;
										if (temp < -m_s[l].tot2)
										{
											kit->second = 0;
										}
										else
										{
											kit->second = temp + m_s[l].tot2;
											temp = 0.0;
											break;
										}
									}
								}
								if (temp != 0.0
									&& -m_s[l].tot2 - temp > 1e-12)
								{
									sprintf(token,
											"Negative concentration in MCD: added %.3e moles %s in cell %d",
											(double) (-m_s[l].tot2 - temp),
											m_s[l].name, jcell);
									warning_msg(token);
								}
							}
							else
								it->second += m_s[l].tot2;
							break;
						}
					}
					if (it == use.Get_solution_ptr()->Get_totals().end())
					{
						use.Get_solution_ptr()->Get_totals()[m_s[l].name] = m_s[l].tot2;
						if (m_s[l].tot2 < 0)
						{
							if (m_s[l].tot2 < -1e-12)
							{
								sprintf(token,
										"Negative concentration in MCD: added %.4e moles %s in cell %d",
										(double) -m_s[l].tot2, m_s[l].name, jcell);
								warning_msg(token);
							}
							use.Get_solution_ptr()->Get_totals()[m_s[l].name] = 0;
						}
					}
				}
			}
		}
	}
	m_s = (struct M_S *) free_check_null(m_s);
	J_ij = (struct J_ij *) free_check_null(J_ij);
	J_ij_il = (struct J_ij *) free_check_null(J_ij_il);
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
fill_m_s(struct J_ij *l_J_ij, int l_J_ij_count_spec)
/* ---------------------------------------------------------------------- */
{
/*  sum up the primary or secondary master_species from solute species
 *      H and O go in tot1&2_h and tot1&2_o
 */
	int j, k, l;
	char *ptr;

	for (j = 0; j < l_J_ij_count_spec; j++)
	{
		{
			char * temp_name = string_duplicate(l_J_ij[j].name);
			ptr = temp_name;
			count_elts = 0;
			get_elts_in_species(&ptr, 1);
			free_check_null(temp_name);
		}
		for (k = 0; k < count_elts; k++)
		{
			if (strcmp(elt_list[k].elt->name, "X") == 0)
				continue;
			if (strcmp(elt_list[k].elt->name, "H") == 0)
			{
				tot1_h += elt_list[k].coef * l_J_ij[j].tot1;
				tot2_h += elt_list[k].coef * l_J_ij[j].tot2;
			}
			else if (strcmp(elt_list[k].elt->name, "O") == 0)
			{
				tot1_o += elt_list[k].coef * l_J_ij[j].tot1;
				tot2_o += elt_list[k].coef * l_J_ij[j].tot2;
			}
			else
			{
				for (l = 0; l < count_m_s; l++)
				{
					if (strcmp(m_s[l].name, elt_list[k].elt->name) == 0)
					{
						m_s[l].tot1 += elt_list[k].coef * l_J_ij[j].tot1;
						m_s[l].tot2 += elt_list[k].coef * l_J_ij[j].tot2;
						break;
					}
				}
				if (l == count_m_s)
				{
					//m_s[l].name = string_hsave(elt_list[k].elt->name);
					m_s[l].name = elt_list[k].elt->name;
					m_s[l].tot1 = elt_list[k].coef * l_J_ij[j].tot1;
					m_s[l].tot2 = elt_list[k].coef * l_J_ij[j].tot2;
					count_m_s++;
				}
			}
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
find_J(int icell, int jcell, LDBLE mixf, LDBLE DDt, int stagnant)
/* ---------------------------------------------------------------------- */
{
	/* mole transfer of the individual master_species:
	 * Eqn 1:
	 * J_ij = DDt * (A_ij / lav) * (-D_i*grad(c) + D_i*z_i*c_i * SUM(D_i*z_i*grad(c)) / SUM(D_i*(z_i)^2*c_i))
	 * regular column, stagnant FALSE:
	 *    D_i = temperature-corrected Dw
	 *    A_ij = A_icell * A_jcell
	 *    A_icell = (L porewater in i_cell / length_icell) / tort_f_icell /
	 *       (length_icell / 2)
	 *    lav = A_icell + A_jcell
	 *    grad(c) is concentration difference in icell and jcell, 
		   for activity corrections see Appelo & Wersin, 2007.
	 *  stagnant TRUE:
	 * J_ij = mixf_ij * (-D_i*grad(c) + D_i*z_i*c_i * SUM(D_i*z_i*grad(c)) / SUM(D_i*(z_i)^2*c_i))
	 *    mixf_ij = mixf / (Dw / init_tort_f) / new_tort_f * new_por / init_por
	 *    mixf is defined in MIX; Dw is default multicomponent diffusion coefficient;
	 *    init_tort_f equals multi_Dpor^(-multi_Dn); new_pf = new tortuosity factor.
	 * Interlayer diffusion (IL) takes the gradient in the equivalent concentrations on X-.
		surface area A for IL:
		stagnant: mixf_il is mixf * por_il / por.
					por_il = interlayer porosity, from -interlayer_D true 'por_il'.
		            por = total porosity, from -multi_D true 'multi_Dpor'.
		            **nov. 12, 2011**: 
							mixf is corrected, * (1 - por_il / por).
							new_pf = (por - por_il)^(-multi_Dn).
		in regular column, A is calc'd from (free + DL porewater) and cell-length.
		for IL: A * por_il / (por - por_il).

		por_il is entered as a single value. It is limited to 0.999 * por.
		por_il in a cell is reduced by conc of X- / (max conc of X- of all cells)

		IL-water = (free + DL porewater) * por_il / (por - por_il).
	 */
	int i, i_max, j, j_max, k, k_il, only_counter, il_calcs;
	int i1;
	LDBLE lav, A1, A2, A_ij, A_ij_il, ddlm, aq1, aq2, mixf_il;
	LDBLE dl_s, dl_aq1, dl_aq2, c_dl, visc1, visc2, dum, dum2, tort1, tort2;
	LDBLE por_il1, por_il2, por_il12;
	LDBLE c, Dz2c, Dz2c_dl, Dz2c_il, aq_il1, aq_il2;
	LDBLE cec1, cec2, cec12, rc1, rc2;
	struct V_M
	{
		LDBLE grad, D, z, Dz, Dzc, Dzc_dl, g_dl;
		int o_c;
	} *V_M, *V_M_il;
	cxxSurface *s_ptr1, *s_ptr2;
	cxxSurfaceCharge *s_charge_ptr, *s_charge_ptr1, *s_charge_ptr2;
	char token[MAX_LENGTH], token1[MAX_LENGTH];

	V_M = V_M_il = NULL;
	/* check for immediate return and interlayer diffusion calcs... */
	if (interlayer_Dflag)
	{
		il_calcs = 1;
		if (icell == 0 && cell_data[0].por_il < interlayer_Dpor_lim)
			il_calcs = 0;
		else if (icell == count_cells &&
				 cell_data[count_cells - 1].por_il < interlayer_Dpor_lim)
			il_calcs = 0;
		else if (icell > 0
				 && (cell_data[icell - 1].por_il < interlayer_Dpor_lim
					 || cell_data[jcell - 1].por_il < interlayer_Dpor_lim))
			il_calcs = 0;
	}
	else
		il_calcs = 0;

	if (stagnant)
	{
		if (!il_calcs && (cell_data[icell - 1].por < multi_Dpor_lim
						  || cell_data[jcell - 1].por < multi_Dpor_lim))
			return (OK);
	}
	else
	{							/* regular column... */
		if (icell == 0)
		{
			if (!il_calcs && cell_data[0].por < multi_Dpor_lim)
				return (OK);
		}
		else if (icell == count_cells)
		{
			if (!il_calcs && cell_data[count_cells - 1].por < multi_Dpor_lim)
				return (OK);
		}
		else
		{
			if (!il_calcs && (cell_data[icell - 1].por < multi_Dpor_lim
							  || cell_data[jcell - 1].por < multi_Dpor_lim))
				return (OK);
		}
	}

	/* do the calcs */
	aq1 = Utilities::Rxn_find(Rxn_solution_map, icell)->Get_mass_water();
	aq2 = Utilities::Rxn_find(Rxn_solution_map, jcell)->Get_mass_water();
	/*
	 * check if DL calculations must be made, find amounts of water...
	 */
	s_charge_ptr1 = s_charge_ptr2 = NULL;
	s_ptr1 = s_ptr2 = NULL;
	dl_s = dl_aq1 = dl_aq2 = 0.0;
	visc1 = visc2 = 1.0;
	only_counter = FALSE;

	s_ptr1 = Utilities::Rxn_find(Rxn_surface_map, icell);
	if (s_ptr1 != NULL)
	{
		if (s_ptr1->Get_dl_type() != cxxSurface::NO_DL)
		{
			if (s_ptr1->Get_only_counter_ions())
				only_counter = TRUE;
			/* find the one (and only one...) immobile surface comp with DL... */
			for (i = 0; i < (int) s_ptr1->Get_surface_comps().size(); i++)
			{
				cxxSurfaceComp * comp_i_ptr = &(s_ptr1->Get_surface_comps()[i]); 
				if (comp_i_ptr->Get_Dw() == 0)
				{
					s_charge_ptr1 = s_ptr1->Find_charge(comp_i_ptr->Get_charge_name());
					dl_aq1 = s_charge_ptr1->Get_mass_water();
					visc1 = s_ptr1->Get_DDL_viscosity();
					/* check for more comps with Dw = 0 */
					for (j = i + 1; j < (int) s_ptr1->Get_surface_comps().size(); j++)
					{
						cxxSurfaceComp * comp_j_ptr = &(s_ptr1->Get_surface_comps()[j]);
						if (comp_j_ptr->Get_Dw() == 0
							&& (comp_j_ptr->Get_charge_name() != 
							comp_i_ptr->Get_charge_name()))
						{
							if (!warn_fixed_Surf)
							{
								k = (int) strcspn(comp_i_ptr->Get_formula().c_str(), "_");
								strncpy(token1, comp_i_ptr->Get_formula().c_str(), k);
								token1[k] = '\0';
								sprintf(token,
									"MCD found more than 1 fixed surface with a DDL,\n\t uses the 1st in alphabetical order: %s.",
								token1);
								warning_msg(token);
								warn_fixed_Surf = 1;
							}
							break;
						}
					}
					break;
				}
			}
		}
	}
	s_ptr2 = Utilities::Rxn_find(Rxn_surface_map, jcell);
	if (s_ptr2 != NULL)
	{
		if (s_ptr2->Get_dl_type() != cxxSurface::NO_DL)
		{
			if (s_ptr2->Get_only_counter_ions())
				only_counter = TRUE;
			for (i = 0; i < (int) s_ptr2->Get_surface_comps().size(); i++)
			{
				cxxSurfaceComp * comp_i_ptr = &(s_ptr2->Get_surface_comps()[i]);
				if (comp_i_ptr->Get_Dw() == 0)
				{
					s_charge_ptr2 = s_ptr2->Find_charge(comp_i_ptr->Get_charge_name());
					dl_aq2 = s_charge_ptr2->Get_mass_water();
					visc2 = s_ptr2->Get_DDL_viscosity();
					/* check for more comps with Dw = 0 */
					for (j = i + 1; j < (int) s_ptr2->Get_surface_comps().size(); j++)
					{
						cxxSurfaceComp * comp_j_ptr = &(s_ptr2->Get_surface_comps()[j]);
						if (comp_j_ptr->Get_Dw() == 0
							&& (comp_j_ptr->Get_charge_name() != 
							comp_i_ptr->Get_charge_name()))
						{
							if (!warn_fixed_Surf)
							{
								k = (int) strcspn(comp_i_ptr->Get_formula().c_str(), "_");
								strncpy(token1, comp_i_ptr->Get_formula().c_str(), k);
								token1[k] = '\0';
								sprintf(token,
									"MCD found more than 1 fixed surface with a DDL,\n\t uses the 1st in alphabetical order: %s.",
									token1);
								warning_msg(token);
								warn_fixed_Surf = 1;
							}
							break;
						}
					}
					break;
				}
			}
		}
	}
	if (!stagnant)
	{
		if (icell == 0)
			visc1 = visc2;
		else if (icell == count_cells)
			visc2 = visc1;
	}
	/* in each cell: DL surface = mass_water_DL / (cell_length)
	   free pore surface = mass_water_free / (cell_length)
	   determine DL surface as a fraction of the total pore surface... */
	if (dl_aq1 > 0)
		dl_s = dl_aq1 / (dl_aq1 + aq1);
	if (dl_aq2 > 0)
	{
		dum = dl_aq2 / (dl_aq2 + aq2);
		if (dl_aq1 > 0)
		/* average the 2... */
			dl_s = (dl_s + dum) / 2;
		else
		/* there is one DL surface... */
			dl_s = dum;
	}

	por_il1 = por_il2 = por_il12 = 0.0;
	cec1 = cec2 = cec12 = rc1 = rc2 = 0.0;
	if (il_calcs)
	{
		/* find interlayer porosity por_il, 
		   make it relative to exchange capacity (mol X/L), highest X in sol_D[1].x_max (mol X / L).
		   Find amounts of IL water and cec.
		   Must do this separately, since por and por_il are in cell_data structure. */
		if (icell == 0)
		{
			por_il1 = sol_D[0].exch_total / aq1 / sol_D[1].x_max *
				cell_data[0].por_il;
			por_il2 = sol_D[1].exch_total / aq2 / sol_D[1].x_max *
				cell_data[0].por_il;
			if (sol_D[0].exch_total > 3e-10 && sol_D[1].exch_total > 3e-10)
				/* take the average... */
				por_il12 = (por_il1 + por_il2) / 2;
			else
				/* at column ends, take the clay... */
				por_il12 = (por_il1 >= por_il2 ? por_il1 : por_il2);
			if (por_il12 > 0.999 * cell_data[0].por)
				por_il12 = 0.999 * cell_data[0].por;

			if (por_il2 > 0.999 * cell_data[0].por)
				por_il2 = 0.999 * cell_data[0].por;
			aq_il2 = (aq2 + dl_aq2) * por_il2 /
                (cell_data[0].por - por_il2);
			/* Assume interlayer water is proportional with CEC... */
			aq_il1 = aq_il2 * sol_D[0].exch_total / sol_D[1].exch_total;
		}
		else if (icell == count_cells)
		{
			por_il1 = sol_D[count_cells].exch_total / aq1 / sol_D[1].x_max *
				cell_data[count_cells - 1].por_il;
			por_il2 = sol_D[count_cells + 1].exch_total / aq2 / sol_D[1].x_max *
				cell_data[count_cells - 1].por_il;
			if (sol_D[count_cells].exch_total > 3e-10 && sol_D[count_cells + 1].exch_total > 3e-10)
				por_il12 = (por_il1 + por_il2) / 2;
			else
				por_il12 = (por_il1 >= por_il2 ? por_il1 : por_il2);
			if (por_il12 > 0.999 * cell_data[count_cells - 1].por)
				por_il12 = 0.999 * cell_data[count_cells - 1].por;

			if (por_il1 > 0.999 * cell_data[count_cells - 1].por)
				por_il1 = 0.999 * cell_data[count_cells - 1].por;
			aq_il1 = (aq1 + dl_aq1) * por_il1 /
                (cell_data[count_cells - 1].por - por_il1);
			aq_il2 = aq_il1 * sol_D[count_cells + 1].exch_total /
				sol_D[count_cells].exch_total;
		}
		else
		{
			por_il1 = sol_D[icell].exch_total / aq1 / sol_D[1].x_max *
				cell_data[icell - 1].por_il;
			por_il2 = sol_D[jcell].exch_total / aq2 / sol_D[1].x_max *
				cell_data[jcell - 1].por_il;
		
			if (sol_D[icell].exch_total > 3e-10 && sol_D[jcell].exch_total > 3e-10)
				por_il12 = (por_il1 + por_il2) / 2;
			else
				por_il12 = (por_il1 >= por_il2 ? por_il1 : por_il2);
			if (por_il12 > 0.999 * cell_data[icell - 1].por || por_il12 > 0.999 * cell_data[jcell - 1].por)
				por_il12 = (cell_data[icell - 1].por >= cell_data[jcell - 1].por ? 
					0.999 * cell_data[jcell - 1].por : 
					0.999 * cell_data[icell - 1].por);
			aq_il1 = (aq1 + dl_aq1) * por_il1 /
                (cell_data[icell - 1].por - por_il1);
			aq_il2 = (aq2 + dl_aq2) * por_il2 / 
                (cell_data[jcell - 1].por - por_il2);
		}
		if (por_il12 == 0)
			il_calcs = 0;
		else
		{
			dum = sol_D[icell].exch_total;
			dum2 = sol_D[jcell].exch_total;
			rc1 = (dum2 > dum ? dum / dum2 : 1);
			rc2 = (dum > dum2 ? dum2 / dum : 1);
			if (sol_D[icell].exch_total > 3e-10)
				cec1 = dum / aq_il1;
			else
				cec1 = 2e-10;
			if (sol_D[jcell].exch_total > 3e-10)
				cec2 = dum2 / aq_il2;
			else
				cec2 = 2e-10;
			/* and the largest for calculating the mass transfer... */
			cec12 = (cec1 > cec2 ? cec1 : cec2);
		}
	}

	/* In stagnant calc's, find mixf_il for IL diffusion, correct mixf.
	   In regular column, find surface areas A and A_il */
	tort1 = tort2 = lav = 1.0;
	A_ij = A_ij_il = mixf_il = 0.0;
	if (stagnant)
	{
		mixf /= (default_Dw * pow(multi_Dpor, multi_Dn) * multi_Dpor);
		dum = (cell_data[icell - 1].por <= cell_data[jcell - 1].por ?
				cell_data[icell - 1].por : cell_data[jcell - 1].por);
		if (il_calcs)
		{
			mixf_il = mixf * por_il12 / interlayer_tortf;
			dum -= por_il12;
		}
		mixf *= (dum * pow(dum, multi_Dn));
	}
	else
	{							/* regular column... */
		if (icell == 0)
		{
			tort1 = tort2 = pow(cell_data[0].por, -multi_Dn);
			lav = cell_data[0].length;
			A_ij = (aq2 + dl_aq2) / (lav * 0.5 * lav);
			if (il_calcs)
				A_ij_il =
					A_ij * por_il12 / ((cell_data[0].por - por_il12) *
											  interlayer_tortf);
			A_ij /= tort1;
		}
		else if (icell == count_cells)
		{
			tort1 = tort2 = pow(cell_data[count_cells - 1].por, -multi_Dn);
			lav = cell_data[count_cells - 1].length;
			A_ij = (aq1 + dl_aq1) / (lav * 0.5 * lav);
			if (il_calcs)
				A_ij_il = A_ij * por_il12 /
					((cell_data[count_cells - 1].por - por_il12) * interlayer_tortf);
			A_ij /= tort2;
		}
		else
		{
			tort1 = pow(cell_data[icell - 1].por, -multi_Dn);
			tort2 = pow(cell_data[jcell - 1].por, -multi_Dn);
			A1 = (aq1 + dl_aq1) / (cell_data[icell - 1].length *
						0.5 * cell_data[icell - 1].length);
			A2 = (aq2 + dl_aq2) / (cell_data[jcell - 1].length *
						0.5 * cell_data[jcell - 1].length);
			if (il_calcs)
			{
				dum = A1 * por_il12 /
					((cell_data[icell - 1].por - por_il12) * interlayer_tortf);
				dum2 = A2 * por_il12 /
					((cell_data[jcell - 1].por - por_il12) * interlayer_tortf);
				A_ij_il = dum * dum2 / (dum + dum2);
			}
			A1 /= tort1;
			A2 /= tort2;
			A_ij = A1 * A2 / (A1 + A2);
		}
	}
	/* diffuse... */
	J_ij_count_spec = 0;
	J_ij_sum = 0.0;
	/*
	 * malloc sufficient space...
	 */
	k = sol_D[icell].count_spec + sol_D[jcell].count_spec;

	J_ij = (struct J_ij *) free_check_null(J_ij);
	J_ij = (struct J_ij *) PHRQ_malloc((size_t) k * sizeof(struct J_ij));
	if (J_ij == NULL)
		malloc_error();

	V_M = (struct V_M *) PHRQ_malloc((size_t) k * sizeof(struct V_M));
	if (V_M == NULL)
		malloc_error();

	for (i = 0; i < k; i++)
	{
		J_ij[i].tot1 = 0.0;
		V_M[i].grad = 0.0;
		V_M[i].D = 0.0;
		V_M[i].Dz = 0.0;
		V_M[i].Dzc = 0.0;
		V_M[i].Dzc_dl = 0.0;
		V_M[i].g_dl = 1.0;
		V_M[i].o_c = 1;
	}
	Dz2c = Dz2c_dl = Dz2c_il = 0.0;

	if (il_calcs)
	{
		/* also for interlayer cations */
		k = sol_D[icell].count_exch_spec + sol_D[jcell].count_exch_spec;

		J_ij_il = (struct J_ij *) free_check_null(J_ij_il);
		J_ij_il = (struct J_ij *) PHRQ_malloc((size_t) k * sizeof(struct J_ij));
		if (J_ij_il == NULL)
			malloc_error();

		V_M_il = (struct V_M *) PHRQ_malloc((size_t) k * sizeof(struct V_M));
		if (V_M_il == NULL)
			malloc_error();
		for (i = 0; i < k; i++)
		{
			J_ij_il[i].tot1 = 0.0;
			V_M_il[i].grad = 0.0;
			V_M_il[i].D = 0.0;
			V_M_il[i].Dz = 0.0;
			V_M_il[i].Dzc = 0.0;
			V_M_il[i].Dzc_dl = 0.0;
			V_M_il[i].g_dl = 1.0;
			V_M_il[i].o_c = 1;
		}
	}
	/*
	 * coefficients in Eqn (1)...
	 */
	i = j = k = k_il = 0;
	i_max = sol_D[icell].count_spec;
	j_max = sol_D[jcell].count_spec;

	while (i < i_max || j < j_max)
	{
		if (j == j_max
			|| (i < i_max
				&& strcmp(sol_D[icell].spec[i].name,
						  sol_D[jcell].spec[j].name) < 0))
		{
			/* species 'name' is only in icell */
			if (il_calcs && sol_D[icell].spec[i].type == EX)
			{
				//J_ij_il[k_il].name = string_hsave(sol_D[icell].spec[i].name);
				J_ij_il[k_il].name = sol_D[icell].spec[i].name;
				V_M_il[k_il].D = sol_D[icell].spec[i].Dwt;
				V_M_il[k_il].z = sol_D[icell].spec[i].z;
				V_M_il[k_il].Dz = V_M_il[k_il].D * V_M_il[k_il].z;
				V_M_il[k_il].Dzc =
					V_M_il[k_il].Dz * sol_D[icell].spec[i].c * cec12 / (2 *
																		V_M_il
																		[k_il].
																		z);
				Dz2c_il += V_M_il[k_il].Dzc * V_M_il[k_il].z;
				V_M_il[k_il].grad = -sol_D[icell].spec[i].c * cec12 / V_M_il[k_il].z;	/* use equivalent fraction */
				k_il++;
			}
			else
			{
				//J_ij[k].name = string_hsave(sol_D[icell].spec[i].name);
				J_ij[k].name = sol_D[icell].spec[i].name;
				V_M[k].D = sol_D[icell].spec[i].Dwt;
				V_M[k].z = sol_D[icell].spec[i].z;
				V_M[k].Dz = V_M[k].D * V_M[k].z;
				V_M[k].Dzc = V_M[k].Dz * sol_D[icell].spec[i].c / 2;
				if (dl_s > 0)
				{
					s_charge_ptr = (dl_aq1 > 0) ? s_charge_ptr1 : s_charge_ptr2;
					LDBLE g = s_charge_ptr->Get_g_map()[V_M[k].z].Get_g();
					{
						{
							if (only_counter)
							{
								if ((s_charge_ptr->Get_la_psi() < 0 && V_M[k].z < 0)
									|| (s_charge_ptr->Get_la_psi() > 0
										&& V_M[k].z > 0))
								{
									V_M[k].o_c = 0;
									V_M[k].Dzc_dl = 0;
								}
								else	/* assume for counter ions in the DDL the free pore space conc's... */
								{
									V_M[k].g_dl = 1.0;
									V_M[k].Dzc_dl =
										V_M[k].Dz * sol_D[icell].spec[i].c /
										2;
								}
							}
							else
							{
								if (dl_aq1 > 0)
								{
									V_M[k].g_dl =
										(1 +
										g * aq1 / dl_aq1) *
										sol_D[icell].spec[i].erm_ddl;
									V_M[k].Dzc_dl =
										V_M[k].Dz * sol_D[icell].spec[i].c /
										2 * V_M[k].g_dl;
								}
								else
									V_M[k].Dzc_dl =
										V_M[k].Dz * sol_D[icell].spec[i].c / 2;
							}
							//break;
						}
					}
					Dz2c_dl += V_M[k].Dzc_dl * V_M[k].z;
				}
				Dz2c += V_M[k].Dzc * V_M[k].z;
				V_M[k].grad = -sol_D[icell].spec[i].c; /* assume d log(gamma) / d log(c) = 0 */
				k++;
			}
			if (i < i_max)
				i++;
		}
		else if (i == i_max
				 || (j < j_max
					 && strcmp(sol_D[icell].spec[i].name,
							   sol_D[jcell].spec[j].name) > 0))
		{
			/* species 'name' is only in jcell */
			if (il_calcs && sol_D[jcell].spec[j].type == EX)
			{
				//J_ij_il[k_il].name = string_hsave(sol_D[jcell].spec[j].name);
				J_ij_il[k_il].name = sol_D[jcell].spec[j].name;
				V_M_il[k_il].D = sol_D[jcell].spec[j].Dwt;
				V_M_il[k_il].z = sol_D[jcell].spec[j].z;
				V_M_il[k_il].Dz = V_M_il[k_il].D * V_M_il[k_il].z;
				V_M_il[k_il].Dzc =
					V_M_il[k_il].Dz * sol_D[jcell].spec[j].c * cec12 / (2 *
																		V_M_il
																		[k_il].
																		z);
				Dz2c_il += V_M_il[k_il].Dzc * V_M_il[k_il].z;
				V_M_il[k_il].grad = sol_D[jcell].spec[j].c * cec12 / V_M_il[k_il].z;	/* use equivalent fraction */
				k_il++;
			}
			else
			{
				//J_ij[k].name = string_hsave(sol_D[jcell].spec[j].name);
				J_ij[k].name = sol_D[jcell].spec[j].name;
				V_M[k].D = sol_D[jcell].spec[j].Dwt;
				V_M[k].z = sol_D[jcell].spec[j].z;
				V_M[k].Dz = V_M[k].D * V_M[k].z;
				V_M[k].Dzc = V_M[k].Dz * sol_D[jcell].spec[j].c / 2;
				if (dl_s > 0)
				{
					s_charge_ptr = (dl_aq2 > 0) ? s_charge_ptr2 : s_charge_ptr1;
					LDBLE g = s_charge_ptr->Get_g_map()[V_M[k].z].Get_g();
					{
						{
							if (only_counter)
							{
								if ((s_charge_ptr->Get_la_psi() < 0 && V_M[k].z < 0)
									|| (s_charge_ptr->Get_la_psi() > 0
										&& V_M[k].z > 0))
								{
									V_M[k].o_c = 0;
									V_M[k].Dzc_dl = 0;
								}
								else	/* assume for counter ions in the DDL the free pore space conc's... */
								{
									V_M[k].g_dl = 1.0;
									V_M[k].Dzc_dl =
										V_M[k].Dz * sol_D[jcell].spec[j].c /
										2;
								}
							}
							else
							{
								if (dl_aq2 > 0)
								{
									V_M[k].g_dl =
										(1 +
										 g * aq2 /
										 dl_aq2) *
										sol_D[jcell].spec[j].erm_ddl;
									V_M[k].Dzc_dl =
										V_M[k].Dz * sol_D[jcell].spec[j].c /
										2 * V_M[k].g_dl;
								}
								else
									V_M[k].Dzc_dl =
										V_M[k].Dz * sol_D[jcell].spec[j].c /
										2;
							}
							//break;
						}
					}
					Dz2c_dl += V_M[k].Dzc_dl * V_M[k].z;
				}
				Dz2c += V_M[k].Dzc * V_M[k].z;
				V_M[k].grad = sol_D[jcell].spec[j].c;  /* assume d log(gamma) / d log(c) = 0 */
				k++;
			}
			if (j < j_max)
				j++;
		}
		else if (strcmp(sol_D[icell].spec[i].name, sol_D[jcell].spec[j].name)
				 == 0)
		{
			/* species 'name' is in both cells */
			if (il_calcs && sol_D[icell].spec[i].type == EX)
			{
				//J_ij_il[k_il].name = string_hsave(sol_D[icell].spec[i].name);
				J_ij_il[k_il].name = sol_D[icell].spec[i].name;
				if (sol_D[icell].spec[i].Dwt == 0
					|| sol_D[jcell].spec[j].Dwt == 0)
					V_M_il[k_il].D = 0.0;
				else
					V_M_il[k_il].D =
						(sol_D[icell].spec[i].Dwt +
						 sol_D[jcell].spec[j].Dwt) / 2;

				V_M_il[k_il].z = sol_D[icell].spec[i].z;
				V_M_il[k_il].Dz = V_M_il[k_il].D * V_M_il[k_il].z;
				V_M_il[k_il].Dzc =
					V_M_il[k_il].Dz * (sol_D[icell].spec[i].c * cec1 +
									   sol_D[jcell].spec[j].c * cec2) / (2 *
																		 V_M_il
																		 [k_il].
																		 z);
				Dz2c_il += V_M_il[k_il].Dzc * V_M_il[k_il].z;
				V_M_il[k_il].grad = (sol_D[jcell].spec[j].c - sol_D[icell].spec[i].c) * cec12 / V_M_il[k_il].z;	/* use equivalent fraction */
				k_il++;
			}
			else
			{
				//J_ij[k].name = string_hsave(sol_D[icell].spec[i].name);
				J_ij[k].name = sol_D[icell].spec[i].name;
				if (sol_D[icell].spec[i].Dwt == 0
					|| sol_D[jcell].spec[j].Dwt == 0)
					V_M[k].D = 0.0;
				else
					V_M[k].D =
						(sol_D[icell].spec[i].Dwt +
						 sol_D[jcell].spec[j].Dwt) / 2;

				V_M[k].z = sol_D[icell].spec[i].z;
				V_M[k].Dz = V_M[k].D * V_M[k].z;
				V_M[k].Dzc =
					V_M[k].Dz * (sol_D[icell].spec[i].c +
								 sol_D[jcell].spec[j].c) / 2;
				/*    Dzc[k] = Dz[k] * (sol_D[icell].spec[i].c > sol_D[jcell].spec[j].c ? sol_D[icell].spec[i].c : sol_D[jcell].spec[j].c);
				 */
				if (dl_s > 0)
				{
					c_dl = 0.0;
					if (dl_aq1 > 0)
					{
						LDBLE g = s_charge_ptr1->Get_g_map()[V_M[k].z].Get_g();
						{
							{
								if (only_counter)
								{
									if ((s_charge_ptr1->Get_la_psi() < 0
										 && V_M[k].z < 0)
										|| (s_charge_ptr1->Get_la_psi() > 0
											&& V_M[k].z > 0))
									{
										V_M[k].o_c = 0;
										V_M[k].Dzc_dl = 0;
									}
									else	/* assume for counter ions in the DDL the free pore space conc's... */
									{
										V_M[k].g_dl = 1.0;
										c_dl = sol_D[icell].spec[i].c / 2;
									}
								}
								else
								{
									V_M[k].g_dl =
										(1 +
										 g * aq1 /
										 dl_aq1) *
										sol_D[icell].spec[i].erm_ddl;
									c_dl =
										sol_D[icell].spec[i].c / 2 *
										V_M[k].g_dl;
								}
								//break;
							}
						}
					}
					else
						c_dl = sol_D[icell].spec[i].c / 2;

					if (dl_aq2 > 0)
					{
						LDBLE g = s_charge_ptr2->Get_g_map()[V_M[k].z].Get_g();
						{
							{
								if (only_counter)
								{
									if ((s_charge_ptr2->Get_la_psi() < 0
										 && V_M[k].z < 0)
										|| (s_charge_ptr2->Get_la_psi() > 0
											&& V_M[k].z > 0))
									{
										V_M[k].o_c = 0;
										V_M[k].Dzc_dl = 0;
									}
									else	/* assume for counter ions in the DDL the free pore space conc's... */
									{
										dum = 1.0;
										c_dl +=
											sol_D[jcell].spec[j].c / 2 * dum;
										V_M[k].g_dl =
											(V_M[k].g_dl + dum) / 2;
									}
								}
								else
								{
									dum =
										(1 +
										 g * aq2 /
										 dl_aq2) *
										sol_D[jcell].spec[j].erm_ddl;
									c_dl += sol_D[jcell].spec[j].c / 2 * dum;
									V_M[k].g_dl = (V_M[k].g_dl + dum) / 2;
								}
								//break;
							}
						}
					}
					else if (V_M[k].o_c == 1)
						c_dl += sol_D[jcell].spec[j].c / 2;

					V_M[k].Dzc_dl = V_M[k].Dz * c_dl;
					Dz2c_dl += V_M[k].Dzc_dl * V_M[k].z;
				}
				Dz2c += V_M[k].Dzc * V_M[k].z;
				V_M[k].grad =
					(sol_D[jcell].spec[j].c - sol_D[icell].spec[i].c);
				ddlm = sol_D[jcell].spec[j].lm - sol_D[icell].spec[i].lm;
				if (fabs(ddlm) > 1e-10)
					V_M[k].grad *=
						(1 +
						 (sol_D[jcell].spec[j].lg -
						  sol_D[icell].spec[i].lg) / ddlm);
				k++;
			}
			if (i < i_max)
				i++;
			if (j < j_max)
				j++;
		}
	}
	/*
	 * fill in J_ij...
	 */
	if (Dz2c == 0)
		k = 0;
	J_ij_count_spec = i_max = k;
	J_ij_sum = 0;
	c = c_dl = 0.0;
	for (i = 0; i < i_max; i++)
	{
		c += V_M[i].Dz * V_M[i].grad;
		c_dl += V_M[i].o_c * V_M[i].Dz * V_M[i].g_dl * V_M[i].grad;
	}
	for (i = 0; i < i_max; i++)
	{
		J_ij[i].tot1 = -V_M[i].D * V_M[i].grad + c * V_M[i].Dzc / Dz2c;
		J_ij[i].tot1 *= (1 - dl_s);
		if (Dz2c_dl > 0)
		{
			dum =
				(-V_M[i].D * V_M[i].g_dl * V_M[i].grad +
				 c_dl * V_M[i].Dzc_dl / Dz2c_dl) * (2 / (visc1 + visc2));
			if ((J_ij[i].tot1 <= 0 && dum <= 0)
				|| (J_ij[i].tot1 > 0 && dum > 0))
			{
				J_ij[i].tot1 += V_M[i].o_c * dum * dl_s;
			}
		}
		/*
		 * multiply with timestep...
		 * for stagnant, DDt = 1, the timestep is in mixf.
		 * NOTE (for stagnant). The timestep calculated in init_mix for MCD (by PHREEQC) must be equal
		 *  or smaller than the timestep taken (by the user) for calculating mixf in MIX.
		 *  Make this timestep small enough, consider the largest Dw in phreeqd.dat (usually H+).
		 *  Dw used for calculating mixf must be given as default_Dw in the input file.
		 */
		if (stagnant)
			J_ij[i].tot1 *= mixf;
		else
			J_ij[i].tot1 *= A_ij * DDt;
		J_ij_sum += V_M[i].z * J_ij[i].tot1;
		J_ij[i].tot2 = J_ij[i].tot1;
	}
	/*
	 * calculate interlayer mass transfer...
	 */
	if (il_calcs && Dz2c_il != 0 && k_il > 0)
	{
		
		cxxExchange *ex_ptr1 = Utilities::Rxn_find(Rxn_exchange_map, icell);
		cxxExchange *ex_ptr2 = Utilities::Rxn_find(Rxn_exchange_map, jcell);
		c = 0.0;
		i_max = k_il;
		for (i = 0; i < i_max; i++)
			c += V_M_il[i].Dz * V_M_il[i].grad;
		for (i = 0; i < i_max; i++)
		{
			J_ij_il[i].tot1 = -V_M_il[i].D * V_M_il[i].grad +
				c * V_M_il[i].Dzc / Dz2c_il;
			if (stagnant)
				J_ij_il[i].tot1 *= mixf_il;
			else
				J_ij_il[i].tot1 *= A_ij_il * DDt;
			J_ij_sum += V_M_il[i].z * J_ij_il[i].tot1;
			J_ij_il[i].tot2 = J_ij_il[i].tot1;
		}

		/* express the transfer in elemental moles... */
		tot1_h = tot1_o = tot2_h = tot2_o = 0.0;
		m_s = (struct M_S *) free_check_null(m_s);
		m_s = (struct M_S *) PHRQ_malloc((size_t) count_elements *
										 sizeof(struct M_S));
		if (m_s == NULL)
			malloc_error();
		for (i1 = 0; i1 < count_elements; i1++)
		{
			m_s[i1].name = NULL;
			m_s[i1].tot1 = 0;
			m_s[i1].tot2 = 0;
		}
		count_m_s = 0;
		fill_m_s(J_ij_il, k_il);

		/* do the mass transfer... */
		if (icell > 0 || stagnant)
		{
			size_t k;
			for (k = 0; k < ex_ptr1->Get_exchange_comps().size(); k++)
			{
				cxxNameDouble nd(ex_ptr1->Get_exchange_comps()[k].Get_totals());
				cxxNameDouble::iterator it = nd.begin();
				i_max = 0;
				for (; it != nd.end(); it++)
				{ 
					if (strcmp("X", it->first.c_str()) == 0)
						i_max = 1;
				}
				if (i_max)
					break;
			}

			if (k < ex_ptr1->Get_exchange_comps().size())
			{
				cxxExchComp &comp_ref = ex_ptr1->Get_exchange_comps()[k];
				cxxNameDouble nd(comp_ref.Get_totals());
				cxxNameDouble::iterator it = nd.begin();
				/* transfer O and H... */
				for (; it != nd.end(); it++)
				{				
					struct element *elt_ptr = element_store(it->first.c_str());
					LDBLE coef = it->second;
					if (strcmp("H", elt_ptr->name) == 0)
					{
						if (coef < rc1 * tot1_h)
						{
							tot1_h -= coef;
							comp_ref.Get_totals().insert("H", 0);
						}
						else
						{
							comp_ref.Get_totals().insert("H", coef - rc1 * tot1_h);
							tot1_h *= (1 - rc1);
						}
					}
					else if (strcmp("O", elt_ptr->name) == 0)
					{
						if (coef < rc1 * tot1_o)
						{
							tot1_o -= coef;
							comp_ref.Get_totals().insert("O", 0);
						}
						else
						{
							comp_ref.Get_totals().insert("O", coef - rc1 * tot1_o);
							tot1_o *= (1 - rc1);
						}
					}
				}
				/* transfer other elements... */
				j_max = 0;		/* if j_max turns true, reallocate the exchange structure */
				for (j = 0; j < count_m_s; j++)
				{
					// Make sure list includes each element
					comp_ref.Get_totals().add(m_s[j].name, 0);

					cxxNameDouble nd(comp_ref.Get_totals());
					cxxNameDouble::iterator it = nd.begin();
					for (; it != nd.end(); it++)
					{
						struct element *elt_ptr = element_store(it->first.c_str());
						LDBLE coef = it->second;
						if (strcmp(m_s[j].name, elt_ptr->name) != 0)
							continue;

						/* rc1 part goes to exchange species... */
						if (coef < rc1 * m_s[j].tot1)
						{
							m_s[j].tot1 -= coef;
							comp_ref.Get_totals().insert(m_s[j].name, 0);
						}
						else
						{
							comp_ref.Get_totals().insert(m_s[j].name, coef - rc1 * m_s[j].tot1);
							m_s[j].tot1 *= (1 - rc1);
						}
					}
				}
			}
		}
		if (icell < count_cells || stagnant)
		{
			size_t k;
			for (k = 0; k < ex_ptr2->Get_exchange_comps().size(); k++)
			{
				cxxExchComp &comp_ref = ex_ptr2->Get_exchange_comps()[k];
				cxxNameDouble nd(comp_ref.Get_totals());
				cxxNameDouble::iterator it = nd.begin();
				i_max = 0;
				for (; it != nd.end(); it++)
				{ 
					if (strcmp("X", it->first.c_str()) == 0)
						i_max = 1;
				}
				if (i_max)
					break;
			}
			if (k < ex_ptr2->Get_exchange_comps().size())
			{
				cxxExchComp &comp_ref = ex_ptr2->Get_exchange_comps()[k];
				cxxNameDouble nd(comp_ref.Get_totals());
				cxxNameDouble::iterator it = nd.begin();
				/* transfer O and H... */
				for (; it != nd.end(); it++)
				{
					struct element *elt_ptr = element_store(it->first.c_str());
					LDBLE coef = it->second;

					if (strcmp("H", elt_ptr->name) == 0)
					{
						if (coef < -rc2 * tot2_h)
						{
							tot2_h += coef;
							comp_ref.Get_totals().insert("H", 0);
						}
						else
						{
							comp_ref.Get_totals().insert("H", coef + rc2 * tot2_h);
							tot2_h *= (1 - rc2);
						}
					}
					else if (strcmp("O", elt_ptr->name) == 0)
					{
						if (coef < -rc2 * tot2_o)
						{
							tot2_o += coef;
							comp_ref.Get_totals().insert("O", 0);
						}
						else
						{
							comp_ref.Get_totals().insert("O", coef + rc2 * tot2_o);
							tot2_o *= (1 - rc2);
						}
					}
				}
				/* transfer other elements... */
				for (j = 0; j < count_m_s; j++)
				{
					// Make sure list includes each element
					comp_ref.Get_totals().add(m_s[j].name, 0);

					cxxNameDouble nd(comp_ref.Get_totals());
					cxxNameDouble::iterator it = nd.begin();
					for (; it != nd.end(); it++)
					{
						struct element *elt_ptr = element_store(it->first.c_str());
						LDBLE coef = it->second;
						if (strcmp(m_s[j].name, elt_ptr->name) != 0)
							continue;

						/* rc2 part goes to exchange species... */
						if (coef < -rc2 * m_s[j].tot2)
						{
							m_s[j].tot2 += coef;
							comp_ref.Get_totals().insert(m_s[j].name, 0);
						}
						else
						{
							comp_ref.Get_totals().insert(m_s[j].name, coef + rc2 * m_s[j].tot2);
							m_s[j].tot2 *= (1 - rc2);
						}
					}
				}
			}
		}
	}
	/* appt 3 July 07, improved convergence without transporting charge imbalance */
	J_ij_sum = 0;
	V_M = (struct V_M *) free_check_null(V_M);
	if (il_calcs)
		V_M_il = (struct V_M *) free_check_null(V_M_il);
	return (il_calcs);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
disp_surf(LDBLE DDt)
/* ---------------------------------------------------------------------- */
/*
 *  Disperse/diffuse surfaces:
 *       surface[n1] = SS(mixf * surface[n2]) + (1 - SS(mixf) * surface[i1])
 *  where SS is Sum of, f2 is a mixing factor.
 *  Implementation:
 *  Mobile surface comps and charges are mixed face by face and 1 by 1 in surface[n1]:
  Step (from cell 1 to count_cells + 1):
 *  1. Surface n1 is made a copy of cell[i1]'s surface, if it exists, or else
 *       b. a copy of the first encountered mobile surface[n2] from cell i2 = i1 - 1.
 *  2  a. Column surfaces are mixed by dispersion/diffusion
 *       b. Column surfaces are mixed by MCD (if multi_Dflag is true)
 *  3. Surfaces n1 and n2 are stored in a temporary surface, with n_user = i1 or i2.
 *  4. When done for all cells, new surfaces are copied into the cells.
 *       NOTE... The surfaces' diffuse_layer, edl and phases/kinetics relations must be identical,
 *		       but only mobile surface_comp's (Dw > 0) and their surface_charge's are transported.
 */
{
	int i, i1, i2, k, k1; 
	int charge_done, surf1, surf2;
	LDBLE f1, f2, mixf, mixf_store, mcd_mixf;
	LDBLE lav, A_ij, por, Dp1, Dp2;
	cxxMix * mix_ptr;
	cxxSurface *surface_ptr1, *surface_ptr2;
	LDBLE viscos_f;
/*
 * temperature and viscosity correction for MCD coefficient, D_T = D_298 * Tk * viscos_298 / (298 * viscos)
 */
	viscos_f = viscosity();
	viscos_f = tk_x * 0.88862 / (298.15 * viscos_f);

	cxxSurface surface_n1, surface_n2;
	cxxSurface *surface_n2_ptr;

	std::map<int, cxxSurface> Rxn_temp_surface_map;

	for (i1 = 1; i1 <= count_cells + 1; i1++)
	{
		if (i1 <= count_cells && cell_data[i1 - 1].por < multi_Dpor_lim)
			continue;

		if (i1 == 1 && bcon_first != 1)
			continue;
		if (i1 == count_cells + 1 && bcon_last != 1)
			continue;

		i2 = i1 - 1;
		if (i2 > 0 && cell_data[i2 - 1].por < multi_Dpor_lim)
			continue;
/*
 * step 1. define surface n1 from cell i1, if it exists...
 */
		surface_n1.Set_n_user(-99);
		surface_n2.Set_n_user(-99);

		surface_ptr1 = Utilities::Rxn_find(Rxn_surface_map, i1);
		if (surface_ptr1 != NULL)
		{
			surface_n1 = *surface_ptr1;
		}

		surface_ptr2 = Utilities::Rxn_find(Rxn_surface_map, i2);
		surf1 = surf2 = 0;
		if (surface_ptr2 != NULL)
		{
			if (surface_ptr2->Get_transport())
				surf2 = 1;
		}
		if (surface_ptr1 != NULL)
		{
			if (surface_ptr1->Get_transport())
				surf1 = 1;
		}

		mixf = mixf_store = 0;
		if (surf1 || surf2)
		{
/*
 * step 2a. Dispersive mixing of mobile surfaces from column cells i2 < i1...
 */
			if (i1 <= count_cells)
				mix_ptr = Utilities::Rxn_find(Dispersion_mix_map, i1);
			else
				mix_ptr = Utilities::Rxn_find(Dispersion_mix_map, count_cells);

			std::vector<int> num;
			std::vector<LDBLE> frac;
			mix_ptr->Vectorize(num, frac);
			for(size_t i3 = 0; i3 < num.size(); i3++)
			{
				if (i1 <= count_cells)
				{
					if (num[i3] == i2)
					{
						mixf = mixf_store = frac[i3];
						break;
					}
				}
				else if (num[i3] == i1)
				{
					mixf = mixf_store = frac[i3];
					break;
				}
			}

			/* step 1b. If cell i1 has no surface, get the mobile comps from surface i2... */
			if (surface_n1.Get_n_user() == -99)
			{
				surface_n1 = mobile_surface_copy(surface_ptr2, i1, false);
				/* limit charges to 1... */
				if (surface_n1.Get_surface_charges().size() > 1)
				{
					std::string charge_name = surface_n1.Get_surface_charges()[0].Get_name();
					surface_n1 = sum_surface_comp(&surface_n1, 0, 
						&surface_n1, charge_name, 1, surface_n1.Get_surface_comps()[0].Get_Dw());
				}
				f1 = 0;
			}
			else
				f1 = 1;
			/* find the (possibly modified) surface in the previous cell i2... */
			f2 = 1;
			if (i2 > 0 || bcon_first == 1)
			{
				surface_n2_ptr = Utilities::Rxn_find(Rxn_temp_surface_map, i2);
				if (surface_n2_ptr != NULL)
				{
					surface_n2 = *surface_n2_ptr;
				}

				/* if not found... */
				else
				{
					/* copy it from surface_ptr2... */
					if (surface_ptr2 != NULL)
					{
						surface_n2 = *surface_ptr2;
						surface_n2.Set_n_user_both(i2);
					}
					else
					{
						/* or make it a mobile copy of the surface in cell i1... */
						surface_n2 = mobile_surface_copy(surface_ptr1, i2, false);
						/* limit charges to 1... */
						if (surface_n2.Get_surface_charges().size() > 1)
						{
							std::string charge_name = surface_n2.Get_surface_charges()[0].Get_name();
							surface_n2 = sum_surface_comp(&surface_n2, 0, &surface_n2, charge_name, 1,
								surface_n2.Get_surface_comps()[0].Get_Dw());
						}
						f2 = 0;
					}
				}
			}

/*
 * For MCD, step 2b. Find physical dimensions and porosity of the cells...
 */
			if (i1 == 1)
			{
				por = cell_data[0].por;
				lav = cell_data[0].length / 2;
				A_ij = Utilities::Rxn_find(Rxn_solution_map, 1)->Get_mass_water() / cell_data[0].length;
			}
			else if (i1 == count_cells + 1)
			{
				por = cell_data[count_cells - 1].por;
				lav = cell_data[count_cells - 1].length / 2;
				A_ij =
					Utilities::Rxn_find(Rxn_solution_map, count_cells)->Get_mass_water() /
					cell_data[count_cells - 1].length;
			}
			else
			{
				por = cell_data[i2 - 1].por;
				lav =
					(cell_data[i1 - 1].length + cell_data[i2 - 1].length) / 2;
				A_ij =
					Utilities::Rxn_find(Rxn_solution_map, i1)->Get_mass_water() / 
					(cell_data[i1 - 1].length * cell_data[i1 - 1].por);
				A_ij +=
					Utilities::Rxn_find(Rxn_solution_map, i2)->Get_mass_water() / 
					(cell_data[i2 - 1].length * cell_data[i2 - 1].por);
				A_ij /= 2;
				A_ij *=
					(cell_data[i1 - 1].por <
					 cell_data[i2 - 1].por ? cell_data[i1 -
													   1].por : cell_data[i2 -
																		  1].
					 por);
			}

			/* mix in comps with the same charge structure... */
			if (surf2)
			{
				for (k = 0; k < (int) surface_ptr2->Get_surface_comps().size(); k++)
				{
					cxxSurfaceComp *comp_k_ptr = &(surface_ptr2->Get_surface_comps()[k]);
					std::string charge_name = comp_k_ptr->Get_charge_name();
					if (comp_k_ptr->Get_Dw() == 0)
						continue;
					
					charge_done = FALSE;
					for (k1 = 0; k1 < k; k1++)
					{
						cxxSurfaceComp *comp_k1_ptr = &(surface_ptr2->Get_surface_comps()[k1]);
						if (comp_k_ptr->Get_charge_name() ==
							comp_k1_ptr->Get_charge_name())
						{
							charge_done = TRUE;
							break;
						}
					}
					if (charge_done)
						continue;

					/* Define the MCD diffusion coefficient... */
					mcd_mixf = 0;
					if (multi_Dflag)
					{
						Dp2 = comp_k_ptr->Get_Dw() * pow(por, multi_Dn);
						Dp1 = 0;
						if (surface_ptr1 != NULL && surface_ptr1->Get_transport())
						{
							for (k1 = 0; k1 < (int) surface_ptr1->Get_surface_comps().size(); k1++)
							{
								cxxSurfaceComp *comp_k1_ptr = &(surface_ptr1->Get_surface_comps()[k1]);
								if (strcmp(comp_k1_ptr->Get_formula().c_str(),
									comp_k_ptr->Get_formula().c_str()) != 0)
									continue;
								Dp1 =
									comp_k1_ptr->Get_Dw() *
									pow(cell_data[i1 - 1].por, multi_Dn);
								break;
							}
						}
						if (Dp1 > 0)
							Dp2 = (Dp2 + Dp1) / 2;

						/* and the mixing factor... */
						mcd_mixf = Dp2 * viscos_f * A_ij * DDt / lav;
					}
					mixf = mixf_store + mcd_mixf;

					if (mixf < 1e-8)
						mixf = 0;
					if (mixf > 0.99999999)
						mixf = 0.99999999;
					if (i1 <= count_cells)
					{
						surface_n1 = sum_surface_comp(&surface_n1, f1, surface_ptr2, charge_name, mixf,
							surface_ptr2->Get_surface_comps()[k].Get_Dw());
						f1 = 1;
					}
					if (i2 > 0)
					{
						surface_n2 = sum_surface_comp(&surface_n2, f2, surface_ptr2, charge_name, -mixf,
							surface_ptr2->Get_surface_comps()[k].Get_Dw());
						f2 = 1;
					}
					surface_n1.Set_n_user_both(i1);
				}
			}
			if (surf1)
			{
				for (k = 0; k < (int) surface_ptr1->Get_surface_comps().size(); k++)
				{
					cxxSurfaceComp * comp_k_ptr = &(surface_ptr1->Get_surface_comps()[k]);
					std::string charge_name = comp_k_ptr->Get_charge_name();
					if (comp_k_ptr->Get_Dw() == 0)
						continue;
					charge_done = FALSE;
					for (k1 = 0; k1 < k; k1++)
					{
						cxxSurfaceComp * comp_k1_ptr = &(surface_ptr1->Get_surface_comps()[k1]);
						if (comp_k_ptr->Get_charge_name() ==
							comp_k1_ptr->Get_charge_name())
						{
							charge_done = TRUE;
							break;
						}
					}
					if (charge_done)
						continue;

					/* Define the MCD diffusion coefficient... */
					mcd_mixf = 0;
					if (multi_Dflag)
					{
						if (i1 <= count_cells)
						{
							Dp1 =
								comp_k_ptr->Get_Dw() *
								pow(cell_data[i1 - 1].por, multi_Dn);
						}
						else
						{
							Dp1 = comp_k_ptr->Get_Dw() * pow(por, multi_Dn);
						}
						Dp2 = 0;
						if (surface_ptr2 != NULL && surface_ptr2->Get_transport())
						{
							for (k1 = 0; k1 < (int) surface_ptr2->Get_surface_comps().size(); k1++)
							{
								cxxSurfaceComp * comp_k1_ptr = &(surface_ptr2->Get_surface_comps()[k1]);
								if (strcmp(comp_k1_ptr->Get_formula().c_str(),
									comp_k_ptr->Get_formula().c_str()) != 0)
									continue;
								Dp2 = comp_k1_ptr->Get_Dw() * pow(por, multi_Dn);
								break;
							}
						}
						if (Dp2 > 0)
							Dp1 = (Dp1 + Dp2) / 2;

						/* and the mixing factor... */
						mcd_mixf = Dp1 * viscos_f * A_ij * DDt / lav;
					}
					mixf = mixf_store + mcd_mixf;

					if (mixf < 1e-8)
						mixf = 0;
					if (mixf > 0.99999999)
						mixf = 0.99999999;
					if (i2 > 0)
					{
						surface_n2 = sum_surface_comp
							(&surface_n2, f2, surface_ptr1, charge_name, mixf,
							 surface_ptr1->Get_surface_comps()[k].Get_Dw());
						f2 = 1;
					}
					if (i1 <= count_cells)
					{
						surface_n1 = sum_surface_comp
							(&surface_n1, f1, surface_ptr1, charge_name, -mixf,
							 surface_ptr1->Get_surface_comps()[k].Get_Dw());
						f1 = 1;
					}
					surface_n2.Set_n_user_both(i2);
				}
			}
		}

/*
 *  Step 3. copy surface[n1] and [n2] in a new temporary surface...
 */
		if (surface_n1.Get_n_user() == -99)
			continue;

		Rxn_temp_surface_map[i1] = surface_n1;
		{
			cxxSurface t;
			surface_n1 = t;
		}
		surface_n1.Set_n_user_both(-99);

		if (surface_n2.Get_n_user() == i2)
		{
			surface_n2.Set_n_user_both(i2);
			Rxn_temp_surface_map[i2] = surface_n2;
			{
				cxxSurface t;
				surface_n2 = t;
			}
			surface_n2.Set_n_user_both(-99);
		}
	}
/*
 * Step 4. Dispersion/diffusion is done. New surfaces can be copied in the cell's surface...
 */
	std::map<int, cxxSurface>::iterator jit = Rxn_temp_surface_map.begin();
	for ( ; jit != Rxn_temp_surface_map.end(); jit++)
	{
		i = jit->first;
		assert (i == jit->second.Get_n_user());
		if ((i == 0 && bcon_first == 1)	|| (i == count_cells + 1 && bcon_last == 1))
		{
			continue;
		}
		if (i >= 0 && i <= 1 + count_cells * (1 + stag_data->count_stag))
		{
			surface_ptr1 = Utilities::Rxn_find(Rxn_surface_map, i);
			if (surface_ptr1 != NULL)
			{
				Rxn_surface_map[i] = jit->second;
			}
			else
			{
				//Add to map
				Rxn_surface_map[i] = jit->second;
			}
		}
	}

	return (OK);
}

/* ---------------------------------------------------------------------- */
cxxSurface Phreeqc::
sum_surface_comp(cxxSurface *source1, LDBLE f1,  cxxSurface *source2,
				 std::string charge_name, LDBLE f2, LDBLE new_Dw)
/* ---------------------------------------------------------------------- */
{
/*
 *   Takes fraction f1 of the 1st surface, adds fraction f2 of the 2nd surface's comps[k] and its charge.
 *   The result goes in target
 */
	int new_n_user;
	cxxSurface *surface_ptr1, *surface_ptr2;
	std::string token;
/*
 *   Find surfaces
 */
	surface_ptr1 = source1;
	if (surface_ptr1 == NULL)
	{
		error_string = sformatf( "Null pointer for surface 1 in sum_surface.");
		error_msg(error_string, STOP);
		input_error++;
		return (ERROR);
	}
	surface_ptr2 = source2;
/*
 *   Store data for structure surface
 */
	new_n_user = surface_ptr1->Get_n_user();
	cxxSurface temp_surface(*surface_ptr1);
	temp_surface.Set_n_user_both(new_n_user);
	temp_surface.Set_description("Copy");
	temp_surface.Set_solution_equilibria(false);
	temp_surface.Set_n_solution(-99);
/*
 *   Multiply component compositions by f1
 */
	temp_surface.multiply(f1);
/*
 *   Add in surface_ptr2
 */
	// Only components with same charge as component k
	cxxSurface addee(*surface_ptr2);
	addee.Get_surface_comps().clear();
	addee.Get_surface_charges().clear();

	for (std::vector<cxxSurfaceComp>::iterator it = surface_ptr2->Get_surface_comps().begin();
		it != surface_ptr2->Get_surface_comps().end(); it++)
	{
		if (it->Get_charge_name() == charge_name)
		{
			addee.Get_surface_comps().push_back(*it);
		}
	}
	for (std::vector<cxxSurfaceCharge>::iterator it = surface_ptr2->Get_surface_charges().begin();
		it != surface_ptr2->Get_surface_charges().end(); it++)
	{
		if (it->Get_name() == charge_name)
		{
			addee.Get_surface_charges().push_back(*it);
		}
	}

	if (f2 == 0)
		f2 = 1e-30;
	temp_surface.add(addee, f2);
	temp_surface.Set_transport(false);
	for (size_t i = 0; i < temp_surface.Get_surface_comps().size(); i++)
	{
		if (temp_surface.Get_surface_comps()[i].Get_charge_name() == charge_name)
		{
			temp_surface.Get_surface_comps()[i].Set_Dw(new_Dw);
		}
		if (temp_surface.Get_surface_comps()[i].Get_Dw() > 0)
		{
			temp_surface.Set_transport(true);
		}
	}
/*
 *   Finish up
 */
	temp_surface.Sort_comps();
	return (temp_surface);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
check_surfaces(cxxSurface *surface_ptr1, cxxSurface *surface_ptr2)
/* ---------------------------------------------------------------------- */
{
/*  checks if surfaces can be mixed...
 */
	int n_user1, n_user2, return_code;

	return_code = OK;
	n_user1 = surface_ptr1->Get_n_user();
	n_user2 = surface_ptr2->Get_n_user();

	if (surface_ptr1->Get_dl_type() != surface_ptr2->Get_dl_type())
	{
		error_string = sformatf(
				"Surfaces %d and %d differ in definition of diffuse layer. Cannot mix.",
				n_user1, n_user2);
		error_msg(error_string, STOP);
		return_code = ERROR;
		input_error++;
	}

	if (surface_ptr1->Get_type() != surface_ptr2->Get_type())
	{
		error_string = sformatf(
				"Surfaces %d and %d differ in use of electrical double layer. Cannot mix.",
				n_user1, n_user2);
		error_msg(error_string, STOP);
		return_code = ERROR;
		input_error++;
	}
	if (surface_ptr1->Get_only_counter_ions() != surface_ptr2->Get_only_counter_ions())
	{
		error_string = sformatf(
				"Surfaces %d and %d differ in use of only counter ions in the diffuse layer. Cannot mix.",
				n_user1, n_user2);
		error_msg(error_string, STOP);
		return_code = ERROR;
		input_error++;
	}
	if (surface_ptr1->Get_related_phases() != surface_ptr2->Get_related_phases())
	{
		error_string = sformatf(
				"Surfaces %d and %d differ in use of related phases (sites proportional to moles of an equilibrium phase). Cannot mix.",
				n_user1, n_user2);
		error_msg(error_string, STOP);
		return_code = ERROR;
		input_error++;
	}
	if (surface_ptr1->Get_related_rate() != surface_ptr2->Get_related_rate())
	{
		error_string = sformatf(
				"Surfaces %d and %d differ in use of related rate (sites proportional to moles of a kinetic reactant). Cannot mix.",
				n_user1, n_user2);
		error_msg(error_string, STOP);
		return_code = ERROR;
		input_error++;
	}

	return (return_code);
}
/* ---------------------------------------------------------------------- */
cxxSurface Phreeqc::
mobile_surface_copy(cxxSurface *surface_old_ptr,
					 int n_user_new, bool move_old)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies mobile comps from surface_old_ptr to surf_ptr1,
 *   comps and charges with Dw > 0 are moved if move_old == TRUE, else copied.
 *   NOTE... when all comps are moved, the old surface is deleted and surfaces are sorted again,
 *		 which will modify pointers and surface numbers.
 *   User number of new surface structure is n_user_new, structure is freed when n_user_new is already defined
 */
	cxxSurface temp_surface(*surface_old_ptr);
/*
 *   Store moving surface's properties in temp_surface
 */
	temp_surface.Set_n_user_both(n_user_new);
	std::ostringstream desc;
	desc << "Surface defined in simulation " << simulation << ".";
	temp_surface.Set_description(desc.str().c_str());
	temp_surface.Set_solution_equilibria(false);
	temp_surface.Set_transport(true);


	size_t count_comps = surface_old_ptr->Get_surface_comps().size();
	int i1, i2;
	i1 = i2 = 0;
	temp_surface.Get_surface_comps().clear();
	temp_surface.Get_surface_charges().clear();
	/* see if comps must be moved, Dw > 0 */
	for (size_t i = 0; i < count_comps; i++)
	{
		cxxSurfaceComp &comp_ref = surface_old_ptr->Get_surface_comps()[i];
		if (comp_ref.Get_Dw() > 0)
		{
			i1++;

			// copy comp
			temp_surface.Get_surface_comps().push_back(comp_ref);

			// copy charge, if needed
			cxxSurfaceCharge *charge_ptr = temp_surface.Find_charge(comp_ref.Get_charge_name());
			if (charge_ptr == NULL)
			{
				i2++;
				cxxSurfaceCharge *old_charge_ptr = surface_old_ptr->Find_charge(comp_ref.Get_charge_name());
				temp_surface.Get_surface_charges().push_back(*old_charge_ptr);
			}
		}
	}

	if (i1 > 0)
	{
		/* OK, store moved parts from old surface, but first:
		   get immobile surface comps from new surface... */
		cxxSurface *surf_ptr = Utilities::Rxn_find(Rxn_surface_map, n_user_new);
		if (surf_ptr != NULL)
		{
			for (size_t k = 0; k < surf_ptr->Get_surface_comps().size(); k++)
			{
				cxxSurfaceComp &comp_ref = surf_ptr->Get_surface_comps()[k];
				if (comp_ref.Get_Dw() > 0)
					continue;
				bool charge_done(false);
				for (size_t k1 = 0; k1 < k; k1++)
				{
					if (surf_ptr->Get_surface_comps()[k1].Get_charge_name() ==
						comp_ref.Get_charge_name())
					{
						charge_done = true;
						break;
					}
				}
				if (charge_done)
					continue;
				temp_surface = sum_surface_comp(&temp_surface, 1, surf_ptr, comp_ref.Get_charge_name(), 1, 0);
			}
		}
		// Return value is temp_surface
		temp_surface.Set_n_user_both(n_user_new);
	}
	
	/* delete moved parts from old surface */
	if (move_old && i1 > 0)
	{
		cxxSurface replace_old(temp_surface);
		int n_user_old = surface_old_ptr->Get_n_user();
		if ((size_t) i1 != count_comps)
		{
			/* redefine old surface with only unmovable comps (Dw = 0) */
			/* other temp_surface members were set above */
			replace_old.Set_n_user_both(n_user_old);
			replace_old.Set_transport(false);
			replace_old.Get_surface_comps().clear();
			replace_old.Get_surface_charges().clear();

			i1 = i2 = 0;
			for (size_t i = 0; i < count_comps; i++)
			{
				if (surface_old_ptr->Get_surface_comps()[i].Get_Dw() == 0)
				{
					i1++;
					// copy surface comp
					cxxSurfaceComp & comp_ref = surface_old_ptr->Get_surface_comps()[i];
					replace_old.Get_surface_comps().push_back(comp_ref);
					cxxSurfaceCharge *charge_ptr = replace_old.Find_charge(comp_ref.Get_charge_name());

					// copy surface charge if necessary
					if (charge_ptr == NULL)
					{
						i2++;
						cxxSurfaceCharge *old_charge_ptr = surface_old_ptr->Find_charge(comp_ref.Get_charge_name());
						replace_old.Get_surface_charges().push_back(*old_charge_ptr);
					}
				}
			}
			if (replace_old.Get_surface_comps().size() == 0)
			{
				Rxn_surface_map.erase(surface_old_ptr->Get_n_user());
			}
			else
			{
				replace_old.Sort_comps();
				Rxn_surface_map[surface_old_ptr->Get_n_user()] = replace_old;
			}
		}
	}
	temp_surface.Sort_comps();
	return (temp_surface);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
diff_stag_surf(int mobile_cell)
/* ---------------------------------------------------------------------- */
/*
 *  Diffuse stagnant and mobile surfaces, following the steps of disp_surf.
 *  First the mobile/stagnant surfaces are mixed, then the stagnant surfaces
 *  when not already done.
 *  If mixing factors among the cells are defined expicitly, it is assumed that
 *  mixing with a lower numbered cell was done when that cell was processed:
 *  for any cell in MCD, need only include the mixing factors for higher numbered cells.
 */
{
	int i, i1, i2, k, k1, ns;
	int charge_done, surf1, surf2;
	LDBLE f1, f2, mixf, mixf_store;
	LDBLE Dp1, Dp2;
	cxxMix *mix_ptr;
	cxxSurface *surface_ptr1, *surface_ptr2;
	LDBLE viscos_f;
/*
 * temperature and viscosity correction for MCD coefficient, D_T = D_298 * Tk * viscos_298 / (298 * viscos)
 */
	viscos_f = viscosity();
	viscos_f = tk_x * 0.88862 / (298.15 * viscos_f);

	cxxSurface surface_n1, surface_n2;
	cxxSurface *surface_n1_ptr = &surface_n1;
	cxxSurface *surface_n2_ptr;
	std::map<int, cxxSurface> Rxn_temp_surface_map;

	for (ns = 0; ns < stag_data->count_stag; ns++)
	{

		i1 = mobile_cell + 1 + ns * count_cells;
		if (ns == 0)
			i1--;

		if (cell_data[i1 - 1].por < multi_Dpor_lim)
			continue;
		surface_n1.Set_n_user_both(-99);
		surface_n2.Set_n_user_both(-99);
		surface_ptr1 = Utilities::Rxn_find(Rxn_surface_map, i1);
/*
 * step 2a. mix surfaces...
 */
		mix_ptr = Utilities::Rxn_find(Rxn_mix_map, i1);
		if (mix_ptr == NULL)
			continue;

		std::vector<int> num;
		std::vector<LDBLE> frac;
		mix_ptr->Vectorize(num, frac);
		for (size_t i3 = 0; i3 < num.size(); i3++)
		{
			if ((i2 = num[i3]) <= i1)
				continue;
			if (cell_data[i2 - 1].por < multi_Dpor_lim)
				continue;
			surface_ptr2 = Utilities::Rxn_find(Rxn_surface_map, i2);
			surf1 = surf2 = 0;
			if (surface_ptr2 != NULL)
			{
				if (surface_ptr2->Get_transport())
					surf2 = 1;
			}
			if (surface_ptr1 != NULL)
			{
				if (surface_ptr1->Get_transport())
					surf1 = 1;
			}
			if (!surf1 && !surf2)
				continue;
			mixf = mixf_store = frac[i3];;

			/* find the (possibly modified) surface in cell i1... */
			f1 = 1;
			surface_n1_ptr = Utilities::Rxn_find(Rxn_temp_surface_map, i1);
			if (surface_n1_ptr != NULL)
			{
				surface_n1 = *surface_n1_ptr;
			}
			/* if not found... */
			else
			{
				/* copy it from surface_ptr1... */
				if (surface_ptr1 != NULL)
				{
					surface_n1 = *surface_ptr1;
				}
				else
				{
					/* or make it a mobile copy of the surface in cell i2... */
					surface_n1 = mobile_surface_copy(surface_ptr2, i1, false);

					/* limit comps to 1... */
					if (surface_n1.Get_surface_charges().size() > 1)
					{
						std::string charge_name = surface_n1.Get_surface_charges()[0].Get_name();
						surface_n1 = sum_surface_comp(&surface_n1, 0, &surface_n1, charge_name, 1,
							surface_n1.Get_surface_comps()[0].Get_Dw());
					}
					f1 = 0;
				}
			}
			/* find the (possibly modified) surface in cell i2... */
			f2 = 1;
			surface_n2_ptr = Utilities::Rxn_find(Rxn_temp_surface_map, i2);

			if (surface_n2_ptr != NULL)
			{
				surface_n2 = *surface_n2_ptr;
			}
			/* if not found... */
			else
			{
				/* copy it from surface_ptr2... */
				if (surface_ptr2 != NULL)
				{
					surface_n2 = *surface_ptr2;
				}
				else
				{
					/* or make it a mobile copy of the surface in cell i1... */
					surface_n2 = mobile_surface_copy(surface_ptr1, i2, false);
					/* limit comps to 1... */
					if (surface_n2.Get_surface_charges().size() > 1)
					{
						std::string charge_name = surface_n2.Get_surface_charges()[0].Get_name();
						surface_n2 = sum_surface_comp(&surface_n2, 0, &surface_n2, charge_name, 1,
							surface_n2.Get_surface_comps()[0].Get_Dw());
					}
					f2 = 0;
				}
			}

			/* For MCD, step 2b. Adapt mixf to default values... */
			if (multi_Dflag)
			{
				mixf_store *=
					(cell_data[i1 - 1].por <=
					 cell_data[i2 - 1].por ? cell_data[i1 -
													   1].por : cell_data[i2 -
																		  1].
					 por);
				mixf_store /= (default_Dw * pow(multi_Dpor, multi_Dn) *
							   multi_Dpor);
			}

			/* mix in comps with the same charge structure... */
			if (surf2)
			{
				for (k = 0; k < (int) surface_ptr2->Get_surface_comps().size(); k++)
				{
					cxxSurfaceComp *comp_k_ptr = &(surface_ptr2->Get_surface_comps()[k]);
					std::string charge_name = comp_k_ptr->Get_charge_name();
					if (comp_k_ptr->Get_Dw() == 0)
						continue;
					charge_done = FALSE;
					for (k1 = 0; k1 < k; k1++)
					{
						cxxSurfaceComp *comp_k1_ptr = &(surface_ptr2->Get_surface_comps()[k1]);
						if (comp_k_ptr->Get_charge_name() ==
							comp_k1_ptr->Get_charge_name())
						{
							charge_done = TRUE;
							break;
						}
					}
					if (charge_done)
						continue;

					/* find diffusion coefficients of surfaces... */
					if (multi_Dflag)
					{
						Dp2 = comp_k_ptr->Get_Dw() *
							pow(cell_data[i2 - 1].por, multi_Dn) * viscos_f;
						Dp1 = 0;
						if (surf1)
						{
							for (k1 = 0; k1 < (int) surface_ptr1->Get_surface_comps().size(); k1++)
							{
								cxxSurfaceComp *comp_k1_ptr = &(surface_ptr1->Get_surface_comps()[k1]);
								if (strcmp(comp_k1_ptr->Get_formula().c_str(),
									comp_k_ptr->Get_formula().c_str()) != 0)
									continue;
								Dp1 =
									comp_k1_ptr->Get_Dw() *
									pow(cell_data[i1 - 1].por,
										multi_Dn) * viscos_f;
								break;
							}
						}
						if (Dp1 > 0)
							Dp2 = (Dp2 + Dp1) / 2;

						/* and adapt the mixing factor... */
						mixf = mixf_store * Dp2;
						mixf /= Utilities::Rxn_find(Rxn_solution_map, i2)->Get_mass_water();
					}

					if (mixf < 1e-8)
						mixf = 0;
					if (mixf > 0.99999999)
						mixf = 0.99999999;
					surface_n1 = sum_surface_comp(&surface_n1, f1, surface_ptr2, charge_name, mixf,
						surface_ptr2->Get_surface_comps()[k].Get_Dw());
					f1 = 1;

					surface_n2 = sum_surface_comp(&surface_n2, f2, surface_ptr2, charge_name, -mixf,
						surface_ptr2->Get_surface_comps()[k].Get_Dw());
				}
			}

			if (surf1)
			{
				for (k = 0; k < (int) surface_ptr1->Get_surface_comps().size(); k++)
				{
					cxxSurfaceComp *comp_k_ptr = &(surface_ptr1->Get_surface_comps()[k]);
					std::string charge_name = comp_k_ptr->Get_charge_name();
					if (comp_k_ptr->Get_Dw() == 0)
						continue;
					charge_done = FALSE;
					for (k1 = 0; k1 < k; k1++)
					{
						cxxSurfaceComp *comp_k1_ptr = &(surface_ptr1->Get_surface_comps()[k1]);
						if (comp_k_ptr->Get_charge_name() ==
							comp_k1_ptr->Get_charge_name())
						{
							charge_done = TRUE;
							break;
						}
					}
					if (charge_done)
						continue;

					/* find diffusion coefficients of surfaces... */
					if (multi_Dflag)
					{
						Dp1 =
							comp_k_ptr->Get_Dw() *
							pow(cell_data[i1 - 1].por, multi_Dn) * viscos_f;

						Dp2 = 0;
						if (surf2)
						{
							for (k1 = 0; k1 < (int) surface_ptr2->Get_surface_comps().size(); k1++)
							{
								cxxSurfaceComp *comp_k1_ptr = &(surface_ptr2->Get_surface_comps()[k1]);
								if (strcmp(comp_k1_ptr->Get_formula().c_str(),
									comp_k_ptr->Get_formula().c_str()) != 0)
									continue;
								Dp2 =
									comp_k1_ptr->Get_Dw() *
									pow(cell_data[i2 - 1].por,
										multi_Dn) * viscos_f;
								break;
							}
						}
						if (Dp2 > 0)
							Dp1 = (Dp1 + Dp2) / 2;

						/* and adapt the mixing factor... */
						mixf = mixf_store * Dp1;
						mixf /= Utilities::Rxn_find(Rxn_solution_map, i1)->Get_mass_water();
					}

					if (mixf < 1e-8)
						mixf = 0;
					if (mixf > 0.99999999)
						mixf = 0.99999999;
					surface_n2 = sum_surface_comp(&surface_n2, f2, surface_ptr1, charge_name, mixf,
						surface_ptr1->Get_surface_comps()[k].Get_Dw());
					f2 = 1;

					surface_n1 = sum_surface_comp(&surface_n1, f1, surface_ptr1, charge_name, -mixf,
						surface_ptr1->Get_surface_comps()[k].Get_Dw());
				}
			}

/*
 *  Step 3. copy surface[n1] and [n2] in a new temporary surface...
 */
			if (surface_n1.Get_n_user() == -99)
				continue;

			surface_n1.Set_n_user_both(i1);
			Rxn_temp_surface_map[i1] = surface_n1;

			assert(surface_n2.Get_n_user() != -99);
			assert(surface_n2.Get_n_user() == i2);
			surface_n2.Set_n_user_both(i2);
			Rxn_temp_surface_map[i2] = surface_n2;
		}
	}
/*
 * Step 4. Diffusion is done. New surfaces can be copied in the cells...
 */
	std::map<int, cxxSurface>::iterator jit = Rxn_temp_surface_map.begin();
	for ( ; jit != Rxn_temp_surface_map.end(); jit++)
	{
		i = jit->first;
		assert(i == jit->second.Get_n_user());
		if ((i == 0 && bcon_first == 1)	|| (i == count_cells + 1 && bcon_last == 1))
		{
			continue;
		}
		if (i >= 0 && i <= 1 + count_cells * (1 + stag_data->count_stag))
		{
			surface_ptr1 = Utilities::Rxn_find(Rxn_surface_map, i);
			if (surface_ptr1 != NULL)
			{
				Rxn_surface_map[i] = jit->second;
			}
			else
			{
				//Add to map
				Rxn_surface_map[i] = jit->second;
			}
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
reformat_surf(const char *comp_name, LDBLE fraction, const char *new_comp_name,
			  LDBLE new_Dw, int l_cell)
/* ---------------------------------------------------------------------- */
{
	cxxSurface *surface_ptr;

	if ((surface_ptr = Utilities::Rxn_find(Rxn_surface_map, l_cell)) == NULL)
		return (OK);
	if (surface_ptr->Find_charge(comp_name) == NULL)
		return (OK);
	// Assume comp_name is charge name
	std::string old_charge_name = comp_name;
	std::string new_charge_name = new_comp_name;
	if (fraction > 0.99999999)
		fraction = 0.99999999;

	cxxSurface temp_surface(*surface_ptr);

	cxxSurface change;

	for (size_t i = 0; i < temp_surface.Get_surface_comps().size(); i++)
	{
		cxxSurfaceComp *comp_ptr = &(temp_surface.Get_surface_comps()[i]);
		std::string charge_name = comp_ptr->Get_charge_name();
		if (charge_name == old_charge_name)
		{
			cxxSurfaceComp comp(*comp_ptr);
			comp.multiply(fraction);
			std::string std_comp_name = comp_ptr->Get_formula();
			Utilities::replace(comp_name, new_comp_name, std_comp_name);
			comp.Set_formula(std_comp_name.c_str());
			comp.Set_charge_name(new_comp_name);
			cxxNameDouble nd;
			cxxNameDouble::iterator it;
			for (it = comp.Get_totals().begin(); it != comp.Get_totals().end(); it++)
			{
				std::string tot_name(it->first);
				Utilities::replace(comp_name, new_comp_name, tot_name);
				nd[tot_name] = it->second;
			}
			comp.Set_totals(nd);
			change.Get_surface_comps().push_back(comp);
			comp_ptr->multiply(1 - fraction);
		}
	}
	for (size_t i = 0; i < temp_surface.Get_surface_charges().size(); i++)
	{
		cxxSurfaceCharge *charge_ptr = &(temp_surface.Get_surface_charges()[i]);
		std::string charge_name = charge_ptr->Get_name();
		if (charge_name == old_charge_name)
		{
			cxxSurfaceCharge charge(*charge_ptr);
			charge.multiply(fraction);
			std::string std_charge_name = charge_ptr->Get_name();
			Utilities::replace(comp_name, new_comp_name, std_charge_name);
			charge.Set_name(std_charge_name.c_str());
			change.Get_surface_charges().push_back(charge);
			charge_ptr->multiply(1 - fraction);
		}
	}
	temp_surface.add(change, 1.0);
	// change Dw
	for (size_t i = 0; i < temp_surface.Get_surface_comps().size(); i++)
	{
		cxxSurfaceComp *comp_ptr = &(temp_surface.Get_surface_comps()[i]);
		std::string charge_name = comp_ptr->Get_charge_name();
		if (charge_name == new_charge_name)
		{
			comp_ptr->Set_Dw(new_Dw);
		}
	}
	temp_surface.Set_transport(false);
	for (size_t i = 0; i < temp_surface.Get_surface_comps().size(); i++)
	{
		if (temp_surface.Get_surface_comps()[i].Get_Dw() > 0)
		{
			temp_surface.Set_transport(true);
			break;
		}
	}
	temp_surface.Sort_comps();
	Rxn_surface_map[l_cell] = temp_surface;
	return OK;
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
viscosity(void)
/* ---------------------------------------------------------------------- */
{
  LDBLE viscos;
/* from Atkins, 1994. Physical Chemistry, 5th ed. */
	viscos =
		pow((LDBLE) 10.,
			-(1.37023 * (tc_x - 20) +
			  0.000836 * (tc_x - 20) * (tc_x - 20)) / (109 + tc_x));
  return viscos;
}
