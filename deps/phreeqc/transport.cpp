#include "Utils.h"
#include "Phreeqc.h"
#include "phqalloc.h"
#include "Exchange.h"
#include "GasPhase.h"
#include "PPassemblage.h"
#include "SSassemblage.h"
#include "cxxKinetics.h"
#include "Solution.h"
#include <limits.h>

LDBLE F_Re3 = F_C_MOL / (R_KJ_DEG_MOL * 1e3);
LDBLE tk_x2; // average tk_x of icell and jcell
LDBLE dV_dcell; // difference in Volt among icell and jcell
int find_current;
char token[MAX_LENGTH];

// implicit...
std::set <std::string> dif_spec_names;
std::set <std::string> dif_els_names;
std::map<int, std::map<std::string, double> > neg_moles;
std::map<std::string, double> els;
double *Ct2, *l_tk_x2, **A, **LU, **mixf, **mixf_stag;
int mixf_comp_size = 0;

struct CURRENT_CELLS
{
	LDBLE dif, ele, R; // diffusive and electric components, relative cell resistance
} *current_cells;
LDBLE sum_R, sum_Rd; // sum of R, sum of (current_cells[0].dif - current_cells[i].dif) * R
struct V_M   // For calculating Vinograd and McBain's zero-charge, diffusive tranfer of individual solutes
{
	LDBLE grad, D, z, c, zc, Dz, Dzc;
	LDBLE b_ij; // harmonic mean of cell properties, with EDL enrichment
};
struct CT /* summed parts of V_M and mcd transfer in a timestep for all cells, for free + DL water */
{
	LDBLE kgw, dl_s, Dz2c, Dz2c_stag, visc1, visc2, J_ij_sum;
	LDBLE A_ij_il, Dz2c_il, mixf_il;
	int J_ij_count_spec, J_ij_il_count_spec;
	struct V_M *v_m, *v_m_il;
	struct J_ij *J_ij, *J_ij_il;
	int count_m_s;
	struct M_S *m_s;
	int v_m_size, J_ij_size, m_s_size;
} *ct = NULL;
struct MOLES_ADDED /* total moles added to balance negative conc's */
{
	char *name;
	LDBLE moles;
} *moles_added;
int count_moles_added;

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

	int punch_boolean = 0;
	LDBLE step_fraction;

	state = TRANSPORT;
	diffc_tr = diffc;
	diffc_max = 0.0;
	transp_surf = warn_fixed_Surf = warn_MCD_X = 0;
	dV_dcell = current_A = 0.0;
	current_cells = NULL;

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
				"Solution %d is needed for transport, but is not defined.", i);
			error_msg(error_string, CONTINUE);
		}
		else
			cell_data[i].temp = use.Get_solution_ptr()->Get_tc();
	}

	if (multi_Dflag)
	{
		sol_D = (struct sol_D *) PHRQ_malloc((size_t) (all_cells)* sizeof(struct sol_D));
		if (sol_D == NULL)
			malloc_error();
		for (i = 0; i < all_cells; i++)
		{
			sol_D[i].count_spec = 0;
			sol_D[i].count_exch_spec = 0;
			sol_D[i].exch_total = 0;
			sol_D[i].x_max = 0;
			sol_D[i].viscos_f = 1.0;
			sol_D[i].tk_x = 298.15;
			sol_D[i].spec = NULL;
			sol_D[i].spec_size = 0;
		}
		//sol_D_dbg = sol_D;

		ct = (struct CT *) PHRQ_malloc((size_t) (all_cells)* sizeof(struct CT));
		if (ct == NULL)
			malloc_error();
		for (i = 0; i < all_cells; i++)
		{
			ct[i].kgw = 1.0;
			ct[i].dl_s = 0.0;
			ct[i].Dz2c = ct[i].Dz2c_stag = 0.0;
			ct[i].visc1 = ct[i].visc2 = 0.0;
			ct[i].J_ij_sum = 0.0;
			ct[i].A_ij_il = 0.0;
			ct[i].Dz2c_il = 0.0;
			ct[i].mixf_il = 0.0;
			ct[i].J_ij_count_spec = 0;
			ct[i].J_ij_il_count_spec = 0;
			ct[i].v_m = NULL;
			ct[i].v_m_il = NULL;
			ct[i].J_ij = NULL;
			ct[i].J_ij_il = NULL;
			ct[i].m_s = NULL;
			ct[i].v_m_size = ct[i].J_ij_size = ct[i].m_s_size = 0;
		}
		count_moles_added = count_elements;
		moles_added = (struct MOLES_ADDED *) PHRQ_malloc((size_t) (count_moles_added)* sizeof(struct MOLES_ADDED));
		if (moles_added == NULL)
			malloc_error();

		for (i = 0; i < count_moles_added; i++)
		{
			moles_added[i].name = NULL;
			moles_added[i].moles = 0;
		}
	}
	try
	{
		/* check solution 0 */
		use.Set_solution_ptr(Utilities::Rxn_find(Rxn_solution_map, 0));
		if (!use.Get_solution_ptr())
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
		else
		{
			if ((cell_data[0].potV = use.Get_solution_ptr()->Get_potV()))
				dV_dcell = 1;
		}

		/* check solution count_cells */
		use.Set_solution_ptr(Utilities::Rxn_find(Rxn_solution_map, count_cells + 1));
		if (!use.Get_solution_ptr())
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
		else
		{
			if ((cell_data[count_cells + 1].potV = use.Get_solution_ptr()->Get_potV()))
				dV_dcell = 1;
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
					cell_data[k].temp = use.Get_solution_ptr()->Get_tc();
				if (n == 1 && implicit && use.Get_solution_ptr() == NULL)
				{
					input_error++;
					error_string = sformatf(
						"Stagnant solution %d not found for implicit diffusion with 1 stagnant layer.\n       Please define it, or set -implicit false.", k);
					error_msg(error_string, CONTINUE);
				}
			}
		}

		if (fix_current && !dV_dcell && !implicit)
		{
			warning_msg("fix_current (A) was defined, but potential in a boundary cell was not.\n\tUsing 1e-8 V in cell 0 for calculating the potentials.");
			cell_data[0].potV = 1e-8;
			dV_dcell = 1;
		}
		if (dV_dcell)
		{
			if (ishift)
			{
				input_error++;
				error_string = sformatf(
					"Electro-diffusion cannot be combined with advective transport.");
				error_msg(error_string, CONTINUE);
			}
			if (!multi_Dflag)
			{
				input_error++;
				error_string = sformatf(
					"Electrical Field (potential) was defined, but needs -multi_D.");
				error_msg(error_string, CONTINUE);
			}
			else
			{
				if (bcon_first != 1)
				{
					bcon_first = 1;
					error_string = sformatf(
						"Electrical Field (potential) was defined, assuming constant boundary condition for first cell.");
					warning_msg(error_string);
				}
				if (bcon_last != 1)
				{
					bcon_last = 1;
					error_string = sformatf(
						"Electrical Field (potential) was defined, assuming constant boundary condition for last cell.");
					warning_msg(error_string);
				}
				current_cells = (struct CURRENT_CELLS *) PHRQ_malloc((size_t)
					(count_cells + 1) * sizeof(struct CURRENT_CELLS));
				if (current_cells == NULL)
					malloc_error();
				for (int i = 0; i < count_cells + 1; i++)
				{
					current_cells[i].dif = 0.0;
					current_cells[i].ele = 0.0;
					current_cells[i].R = 0.0;
				}
			}
		}
		if (implicit && interlayer_Dflag)
		{
			input_error++;
			error_string = sformatf(
				"Interlayer diffusion cannot be calculated implicitly. Please remove -implicit true.");
			error_msg(error_string, CONTINUE);
		}
		if (implicit && !multi_Dflag)
		{
			input_error++;
			error_string = sformatf(
				"Implicit diffusion needs diffusion coefficients for individual species. Please add -multi_d true.");
			error_msg(error_string, CONTINUE);
		}
		if (implicit && !timest)
		{
			error_string = sformatf(
				"Time step not defined for diffusion, using -time_step 1 # seconds");
			warning_msg(error_string);
			timest = 1;
		}
		//if (implicit && stag_data->count_stag > 1)
		//{
		//	error_string = sformatf(
		//		"Sorry, implicit diffusion can handle only 1 stagnant layer for now. Please remove -implicit true, or set -stagnant 1.");
		//	//error_msg(error_string, CONTINUE);
		//}
		if (implicit && current_cells == NULL)
		{
			current_cells = (struct CURRENT_CELLS *) PHRQ_malloc((size_t)
				(count_cells + 1) * sizeof(struct CURRENT_CELLS));
			if (current_cells == NULL)
				malloc_error();
			for (int i = 0; i < count_cells + 1; i++)
			{
				current_cells[i].dif = 0.0;
				current_cells[i].ele = 0.0;
				current_cells[i].R = 0.0;
			}
		}

		/*
		* First equilibrate solutions
		*/
		dup_print("Equilibrating initial solutions", TRUE);
		transport_step = 0;
		for (i = 0; i <= count_cells + 1; i++)
		{
			if (!implicit && ((bcon_first == 2 && i == 0) ||
				(bcon_last == 2 && i == count_cells + 1)))
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
				fill_spec(cell_no, 0);
			}
			print_punch(i, true);

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
						fill_spec(cell_no, i);
					}
					print_punch(k, true);
					saver();
				}
			}
		}
		/*
		*  Initialize mixing factors, define kinetics times
		*  for multicomponent diffusion, limit mixing by diffc_max (usually from H+)
		*/
		if (multi_Dflag)
			diffc_tr = diffc_max;
		if ((stag_data->exch_f > 0) && (stag_data->count_stag == 1))
		{
			/* multi_D calc's are OK if all cells have the same amount of water */
			 if (multi_Dflag == TRUE)
			{
				sprintf(token, "multi_D calc's and stagnant: define MIXing factors explicitly, or \n\t give in -multi_D the Dw used for calculating the mobile-immobile exchange factor.");
				warning_msg(token);
			}
			
			Rxn_mix_map.clear();
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
		* With count_stag = 1, mix factors are calculated from exchange factor alpha
		* (= exch_f), mobile th_m and immobile th_im porosity.
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
			for (j = 1; j <= count_cells; j++)
			{
				j_imm = j + (1 + count_cells);
				if (Utilities::Rxn_find(Rxn_solution_map, j) == NULL)
					error_msg
					("Could not find mobile cell solution in TRANSPORT.",
						STOP);
				if (Utilities::Rxn_find(Rxn_solution_map, j_imm) == NULL)
					//error_msg
						//("Could not find immobile cell solution in TRANSPORT.",
							//STOP);
					continue;
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
		if (implicit)
			sprintf(token, "\nCalculating implicit transport: %d (mobile) cells, %d shifts, %d mixruns, max. mixf = %g.\n\n",
				count_cells, count_shifts - transport_start + 1, nmix, max_mixf);
		else
			sprintf(token, "\nCalculating transport: %d (mobile) cells, %d shifts, %d mixruns...\n\n",
				count_cells, count_shifts - transport_start + 1, nmix);
		screen_msg(token);
		max_iter = 0;
		for (transport_step = transport_start; transport_step <= count_shifts;
			transport_step++)
		{
			/*
			*  Set initial moles of phases
			*/
			for (i = 0; i <= count_cells + 1; i++)
			{
				if (!dV_dcell && !implicit && (i == 0 || i == count_cells + 1))
					continue;
				set_initial_moles(i);
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
						set_initial_moles(k);
				}
			}
			/*
			* Start diffusing if boundary cond = 1, (fixed c, or closed)
			*/
			if (b_c == 1)
			{
				/* For half of mixing steps */
				for (j = 1; j <= floor((LDBLE)nmix / 2); j++)
				{
					rate_sim_time_start =
						(transport_step - 1) * timest + (j - 1) * kin_time;
					rate_sim_time = rate_sim_time_start + kin_time;

					mixrun = j;
					if (multi_Dflag && j == floor((LDBLE)nmix / 2))
					{
						//sprintf(token,
						//		"Transport step %3d. Multicomponent diffusion run %3d.",
						//		transport_step, j);
						//dup_print(token, FALSE);
					}
					else if (!multi_Dflag)
					{
						sprintf(token, "Transport step %3d. Mixrun %3d.",
							transport_step, j);
						dup_print(token, FALSE);
					}

					if (heat_nmix > 0 && !implicit)
					{
						heat_mix(heat_nmix);
						/* equilibrate again ... */
						for (i = 1; i <= count_cells; i++)
						{
							cell_no = i;
							set_and_run_wrapper(i, NOMIX, FALSE, i, 0.0);
							if (multi_Dflag)
								fill_spec(i, i - 1);
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
					if (implicit)
						diffuse_implicit(stagkin_time, stag_data->count_stag);
					else if (multi_Dflag)
						multi_D(stagkin_time, 1, FALSE);

					for (i = 0; i <= count_cells + 1; i++)
					{
						if (!dV_dcell && (i == 0 || i == count_cells + 1) && !implicit)
						{
							if (j == nmix && stag_data->count_stag == 0 &&
								(cell_data[0].print || cell_data[0].punch ||
									cell_data[count_cells + 1].print || cell_data[count_cells + 1].punch))
								print_punch(i, false);
							continue;
						}
						if (overall_iterations > max_iter)
							max_iter = overall_iterations;
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

						if (i == 0 || i == count_cells + 1)
						{
							run_reactions(i, kin_time, NOMIX, step_fraction); // nsaver = i
						}
						else
						{
							run_reactions(i, kin_time, DISP, step_fraction); // n_saver = -2
						}
						if (multi_Dflag)
							fill_spec(i, 0);

						/* punch and output file */
						if (ishift == 0 && j == nmix && (stag_data->count_stag == 0 || (implicit && stag_data->count_stag == 1)))
							print_punch(i, true);
						if (i > 1)
							Utilities::Rxn_copy(Rxn_solution_map, -2, i - 1);
						saver();
					}

					if (!dV_dcell)
						Utilities::Rxn_copy(Rxn_solution_map, -2, count_cells);
					/* Stagnant zone mixing after completion of each
					diffusive/dispersive step ...  */
					rate_sim_time_start =
						(transport_step - 1) * timest + (j - 1) * stagkin_time;
					rate_sim_time = rate_sim_time_start + stagkin_time;

					if (stag_data->count_stag > 0)
					{
						if (ishift == 0 && j == nmix)
							punch_boolean = TRUE;
						else
							punch_boolean = FALSE;
						for (i = 0; i <= count_cells + 1; i++) // allow for stagnant cell mixing with boundary cells
						{
							mix_stag(i, stagkin_time, punch_boolean, step_fraction);
						}
					}
					if (ishift == 0 && j == nmix && stag_data->count_stag > 0)
					{
						for (n = 1; n <= stag_data->count_stag; n++)
						{
							for (i = 1; i <= count_cells; i++)
							{
								k = i + 1 + n * count_cells;
								if (Utilities::Rxn_find(Rxn_solution_map, k) == NULL)
									continue;
								print_punch(k, false);
							}
						}
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
							if ((Utilities::Rxn_find(Rxn_surface_map, i) != NULL) &&
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
				if (nmix == 0 && heat_nmix > 0 && !implicit)
				{
					heat_mix(heat_nmix);
					/* equilibrate again ... */
					for (i = 1; i <= count_cells; i++)
					{
						cell_no = i;
						set_and_run_wrapper(i, NOMIX, FALSE, i, 0.0);
						if (multi_Dflag)
							fill_spec(i, i - 1);
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
						fill_spec(i, i - 1);
					if (overall_iterations > max_iter)
						max_iter = overall_iterations;
					if (nmix == 0 && stag_data->count_stag == 0)
						print_punch(i, true);
					if (i == first_c && count_cells > 1)
						kin_time = kin_time_save;
					saver();

					/* If nmix is zero, stagnant zone mixing after advective step ... */
					if ((nmix == 0) && (stag_data->count_stag > 0))
					{
						mix_stag(i, stagkin_time, TRUE, step_fraction);
					}
				}
				if (nmix == 0 && stag_data->count_stag > 0)
				{
					for (n = 1; n <= stag_data->count_stag; n++)
					{
						for (i = 1; i <= count_cells; i++)
						{
							k = i + 1 + n * count_cells;
							if (Utilities::Rxn_find(Rxn_solution_map, k) == NULL)
								continue;
							print_punch(k, false);
						}
					}
				}
			}
			/*
			* Further dispersive and diffusive transport
			*/
			if (b_c != 1)
				j = 1;
			for (; j <= nmix; j++) // loop on j
			{
				if (multi_Dflag && j == nmix && (transport_step % print_modulus == 0))
				{
					sprintf(token,
						"Transport step %3d. Multicomponent diffusion run %3d.",
						transport_step, j);
					dup_print(token, FALSE);
				}
				else if (!multi_Dflag)
				{
					sprintf(token, "Transport step %3d. Mixrun %3d.",
						transport_step, j);
					dup_print(token, FALSE);
				}
				rate_sim_time_start =
					(transport_step - 1) * timest + (j - 1) * kin_time;
				if (ishift != 0)
					rate_sim_time_start += kin_time;
				rate_sim_time = rate_sim_time_start + kin_time;

				if (heat_nmix > 0 && !implicit)
				{
					heat_mix(heat_nmix);
					/* equilibrate again ... */
					for (i = 1; i <= count_cells; i++)
					{
						cell_no = i;
						set_and_run_wrapper(i, NOMIX, FALSE, i, 0.0);
						if (multi_Dflag)
							fill_spec(i, i - 1);
						saver();
					}
				}
				if (transp_surf)
				{
					if (disp_surf(stagkin_time) == ERROR)
						error_msg("Error in surface transport, stopping.", STOP);
				}
				if (implicit)
					diffuse_implicit(stagkin_time, stag_data->count_stag);
				else if (multi_Dflag)
					multi_D(stagkin_time, 1, FALSE);

				/* for each cell in column */
				for (i = 0; i <= count_cells + 1; i++)
				{
					if (!dV_dcell && (i == 0 || i == count_cells + 1) && !implicit)
					{
						if (j == nmix && stag_data->count_stag == 0 &&
							(cell_data[0].print || cell_data[0].punch ||
								cell_data[count_cells + 1].print || cell_data[count_cells + 1].punch))
							print_punch(i, false);
						continue;
					}
					if (overall_iterations > max_iter)
						max_iter = overall_iterations;
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

					if (i == 0 || i == count_cells + 1)
						run_reactions(i, kin_time, NOMIX, step_fraction);
					else
					{
						run_reactions(i, kin_time, DISP, step_fraction);
					}
					if (multi_Dflag == TRUE)
						fill_spec(i, 0);
					if (j == nmix && (stag_data->count_stag == 0 || (implicit && stag_data->count_stag == 1)))
						print_punch(i, true);
					if (i > 1)
						Utilities::Rxn_copy(Rxn_solution_map, -2, i - 1);
					saver();
				}
				if (!dV_dcell)
					Utilities::Rxn_copy(Rxn_solution_map, -2, count_cells);
				/* Stagnant zone mixing after completion of each diffusive/dispersive step ... */
				rate_sim_time_start =
					(transport_step - 1) * timest + (j - 1) * stagkin_time;
				rate_sim_time = rate_sim_time_start + stagkin_time;

				if (stag_data->count_stag > 0)
				{
					if (j == nmix)
						punch_boolean = TRUE;
					else
						punch_boolean = FALSE;
					for (i = 0; i <= count_cells + 1; i++) // allow for stagnant cell mixing with boundary cells
					{
						mix_stag(i, stagkin_time, punch_boolean, step_fraction);
					}
				}
				if (j == nmix && ((stag_data->count_stag > 0/* && !implicit) || (implicit && stag_data->count_stag == 1*/)))
				{
					for (n = 1; n <= stag_data->count_stag; n++)
					{
						for (i = 1; i <= count_cells; i++)
						{
							k = i + 1 + n * count_cells;
							if (Utilities::Rxn_find(Rxn_solution_map, k) == NULL)
								continue;
							cell_no = k;
							print_punch(k, false);
						}
					}
				}
			}
			if (dump_modulus != 0 && (transport_step % dump_modulus) == 0)
				dump();
		}
		screen_msg("\n");

		if (implicit && neg_moles.size())
		{
			std::map<int, std::map<std::string, double> > ::iterator it1;
			std::map<std::string, double> ::iterator it2;
			for (i = 0; i <= all_cells; i++)
			{
				if ((it1 = neg_moles.find(i)) != neg_moles.end() && (els = it1->second).size())
				{
					for (it2 = els.begin(); it2 != els.end(); it2++)
					{
						for (int i1 = 0; i1 < count_moles_added; i1++)
						{
							if (moles_added[i1].name && !strcmp(moles_added[i1].name, it2->first.c_str()))
							{
								moles_added[i1].moles -= it2->second;
								break;
							}
							else if (!moles_added[i1].moles)
							{
								moles_added[i1].name = string_duplicate(it2->first.c_str());
								moles_added[i1].moles -= it2->second;
								break;
							}
						}
					}
				}
			}
		}

		if (multi_Dflag && moles_added[0].moles > 0)
		{
			sprintf(token,
				"\nFor balancing negative concentrations in MCD, added in total to the system:");
			if (phrq_io)
				phrq_io->warning_msg(token);
			for (i = 0; i < count_moles_added; i++)
			{
				if (!moles_added[i].moles)
					break;
				sprintf(token,
					"\t %.4e moles %s.",
					(double)moles_added[i].moles, moles_added[i].name);
				if (phrq_io)
					phrq_io->warning_msg(token);
			}
		}
	}
	catch (...)
	{
		transport_cleanup();
	    throw;
	}
	transport_cleanup();
	initial_total_time += rate_sim_time;
	rate_sim_time = 0;
	mass_water_switch = FALSE;
	return (OK);
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
transport_cleanup(void)
/* ---------------------------------------------------------------------- */
{
	int i;
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
		heat_mix_array = (LDBLE *)free_check_null(heat_mix_array);
		temp1 = (LDBLE *)free_check_null(temp1);
		temp2 = (LDBLE *)free_check_null(temp2);
	}
	if (multi_Dflag)
	{
		for (i = 0; i < all_cells; i++)
		{
			sol_D[i].spec = (struct spec *) free_check_null(sol_D[i].spec);
		}
		sol_D = (struct sol_D *) free_check_null(sol_D);
		for (int i = 0; i < all_cells; i++)
		{
			ct[i].v_m = (struct V_M *) free_check_null(ct[i].v_m);
			ct[i].v_m_il = (struct V_M *) free_check_null(ct[i].v_m_il);
			ct[i].J_ij = (struct J_ij *) free_check_null(ct[i].J_ij);
			ct[i].J_ij_il = (struct J_ij *) free_check_null(ct[i].J_ij_il);
			ct[i].m_s = (struct M_S *) free_check_null(ct[i].m_s);
		}
		ct = (struct CT *) free_check_null(ct);
		for (int i = 0; i < count_moles_added; i++)
		{
			moles_added[i].name = (char *)free_check_null(moles_added[i].name);
		}
		moles_added = (struct MOLES_ADDED *) free_check_null(moles_added);
	}
	if (implicit)
	{
		int l_stag = (stag_data->count_stag < 2 ? stag_data->count_stag : 0);
		Ct2 = (LDBLE *)free_check_null(Ct2);
		l_tk_x2 = (LDBLE *)free_check_null(l_tk_x2);
		if (A)
		{
			for (i = 0; i < count_cells + 2 + l_stag * count_cells; i++)
			{
				A[i] = (LDBLE *)free_check_null(A[i]);
				LU[i] = (LDBLE *)free_check_null(LU[i]);
			}
		}
		if (mixf)
		{
			for (i = 0; i < count_cells + 2; i++)
			{
				mixf[i] = (LDBLE *)free_check_null(mixf[i]);
				if (l_stag)
					mixf_stag[i] = (LDBLE *)free_check_null(mixf_stag[i]);
				if (!dV_dcell)
				{
					cell_data[i].potV = 0.0;
					use.Set_solution_ptr(Utilities::Rxn_find(Rxn_solution_map, i));
					use.Get_solution_ptr()->Set_potV(0);
				}
			}
		}
		A = (LDBLE **)free_check_null(A);
		LU = (LDBLE **)free_check_null(LU);
		mixf = (LDBLE **)free_check_null(mixf);
		mixf_stag = (LDBLE **)free_check_null(mixf_stag);
		dif_spec_names.clear();
		mixf_comp_size = 0;
	}
	current_cells = (struct CURRENT_CELLS *) free_check_null(current_cells);
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
print_punch(int i, boolean active)
/* ---------------------------------------------------------------------- */
{
	if ((!(cell_data[i].punch && (transport_step % punch_modulus == 0)) &&
		!(cell_data[i].print && (transport_step % print_modulus == 0))) ||
		(bcon_first == 2 && i == 0) ||
		(bcon_last == 2 && i == count_cells + 1))
		return;
	if (!active)
		run_reactions(i, 0, NOMIX, 0);
	cell_no = i;
	use.Set_kinetics_ptr(Utilities::Rxn_find(Rxn_kinetics_map, i));
	if (use.Get_kinetics_ptr() != NULL)
	{
		use.Set_n_kinetics_user(i);
		use.Set_kinetics_in(true);
	}
	if (cell_data[i].punch && (transport_step % punch_modulus == 0))
		punch_all();
	if (cell_data[i].print && (transport_step % print_modulus == 0))
		print_all();

	/* maybe sorb a surface component... */
	if (change_surf_count > 0)
	{
		for (int k = 0; k < change_surf_count; k++)
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
		save.n_surface_user = save.n_solution_user;
		save.n_surface_user_end = save.n_solution_user_end;
	}
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
init_mix(void)
/* ---------------------------------------------------------------------- */
{
	LDBLE lav, dav = 0.0, mf12, maxmix, corr_disp, diffc_here, mD;
	bool warning = false;
	int i, l_nmix;
	LDBLE *m, *m1;
	m = (LDBLE *)PHRQ_malloc((count_cells + 1) * sizeof(LDBLE));
	if (m == NULL)
		malloc_error();
	m1 = (LDBLE *)PHRQ_malloc((count_cells + 1) * sizeof(LDBLE));
	if (m1 == NULL)
		malloc_error();
	for (i = 0; i < count_cells + 1; i++)
	{
		m[i] = 0;
		m1[i] = 0;
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
		if (dV_dcell)
			dV_dcell = (cell_data[count_cells + 1].potV - cell_data[0].potV) / count_cells;
		for (i = 1; i <= count_cells; i++)
		{
			// estimate the number of mixes
			if (i < count_cells)
			{
				lav = (cell_data[i + 1].length + cell_data[i].length) / 2;
				mD = diffc_max * timest / (lav * lav);
				if (mD > maxmix)
					maxmix = mD;
			}
			if (ishift != 0)
			{
				/* define dispersive mixing (diffusion is done in multi_D),
				with harmonic mean of cell properties (manual eqn 49):
				M1 = -a1 * (v1 * timest / 2) * A1 * (c12 - c1) / (h1 / 2)
				a1 is dispersivity in cell 1, A1 = V1 / h1, v1 * timest / 2 = h1 / 2.
				All cells have equal V in advective transport, delta(c) = M1 / V1.
				*/
				if (i < count_cells)
				{
					if (cell_data[i].disp)
						dav = cell_data[i].length / cell_data[i].disp;
					if (cell_data[i + 1].disp)
						dav += cell_data[i + 1].length / cell_data[i + 1].disp;
					if (dav)
						m1[i] = 2 * corr_disp / dav;		/* m1[i] has mixf with higher cell */
				}
				if (i > 1)
				{
					if (cell_data[i].disp)
						dav = cell_data[i].length / cell_data[i].disp;
					if (cell_data[i - 1].disp)
						dav += cell_data[i - 1].length / cell_data[i - 1].disp;
					if (dav)
						m[i] = 2 * corr_disp / dav;		/* m[i] has mixf with lower cell */
				}
				mf12 = m[i] + m1[i];
				if (mf12 > maxmix)
					maxmix = mf12;
			}
		}
		/*
		* Also for boundary cells
		*/
		if (bcon_first == 1)
		{
			mD = 2 * diffc_max * timest / (cell_data[1].length * cell_data[1].length);
			if (mD > maxmix)
				maxmix = mD;
			if (ishift != 0)
			{
				m[1] = 2 * cell_data[1].disp / cell_data[1].length * corr_disp;
				mf12 = m[1] + m1[1];
				if (mf12 > maxmix)
					maxmix = mf12;
			}
		}
		if (bcon_last == 1)
		{
			mD = 2 * diffc_max * timest / (cell_data[count_cells].length * cell_data[count_cells].length);
			if (mD > maxmix)
				maxmix = mD;
			if (ishift != 0)
			{
				m1[count_cells] = 2 * cell_data[count_cells].disp / cell_data[count_cells].length * corr_disp;
				mf12 = m[count_cells] + m1[count_cells];
				if (mf12 > maxmix)
					maxmix = mf12;
			}
		}
		/*
		* Find number of mixes
		*/
		if (maxmix == 0)
		{
			l_nmix = 0;
			if (mcd_substeps > 1 && stag_data->count_stag > 0)
				l_nmix = (int) ceil(mcd_substeps);
		}
		else if (implicit)
		{
			l_nmix = 1;
			if (maxmix > max_mixf)
				l_nmix = 1 + (int)floor(maxmix / max_mixf);
			if (ishift != 0 && (bcon_first == 1 || bcon_last == 1))
			{
				if (l_nmix < 2)
					l_nmix = 2;
			}
			if (mcd_substeps > 1)
				l_nmix = (int)ceil(l_nmix * mcd_substeps);
		}
		else
		{
			if (2.25 * maxmix + 1.0 > (double)INT_MAX)
			{
				m = (LDBLE *)free_check_null(m);
				m1 = (LDBLE *)free_check_null(m1);
				char token[MAX_LENGTH];
				sprintf(token, "Calculated number of mixes %g, is beyond program limit,\nERROR: please set implicit true, or decrease time_step, or increase cell-lengths.", 2.25 * maxmix);
				error_msg(token, STOP);
			}
			if (bcon_first == 1 || bcon_last == 1)
				l_nmix = 1 + (int)floor(2.25 * maxmix);
			else
				l_nmix = 1 + (int)floor(1.5 * maxmix);

			if (ishift != 0 && (bcon_first == 1 || bcon_last == 1))
			{
				if (l_nmix < 2)
					l_nmix = 2;
			}
			if (mcd_substeps > 1)
				l_nmix = (int)ceil(l_nmix * mcd_substeps);
		}

		for (i = 1; i <= count_cells; i++)
		{
			m[i] /= l_nmix;
			m1[i] /= l_nmix;
			/*
			* Fill mix structure
			*/
			cxxMix temp_mix;
			temp_mix.Set_n_user(i);
			temp_mix.Set_n_user_end(i);

			temp_mix.Add(i - 1, m[i]);
			temp_mix.Add(i + 1, m1[i]);
			temp_mix.Add(i, 1.0 - m[i] - m1[i]);
			Dispersion_mix_map[i] = temp_mix;
		}
		m = (LDBLE *)free_check_null(m);
		m1 = (LDBLE *)free_check_null(m1);
		return l_nmix;
	}
	else // multi_D false
	{
		diffc_here = 2 * diffc_tr * timest;
		/*
		* Define mixing factors among inner cells
		*/
		for (i = 1; i <= count_cells; i++)
		{
			if (i < count_cells)
			{
				// find mix with higher numbered cell...
				if (ishift != 0)
				{
					if (cell_data[i].disp)
						dav = cell_data[i].length / cell_data[i].disp;
					if (cell_data[i + 1].disp)
						dav += cell_data[i + 1].length / cell_data[i + 1].disp;
					if (dav)
						m1[i] = 2 / dav;
				}
				// add diffusive mixf...
				m1[i] += diffc_here /
					(cell_data[i].length * cell_data[i].length + cell_data[i].length * cell_data[i + 1].length);
				m1[i] *= corr_disp;		/* m1[i] has mixf with higher cell */
			}
			if (i > 1)
			{
				// and with lower numbered cell...
				if (ishift != 0)
				{
					if (cell_data[i].disp)
						dav = cell_data[i].length / cell_data[i].disp;
					if (cell_data[i - 1].disp)
						dav += cell_data[i - 1].length / cell_data[i - 1].disp;
					if (dav)
						m[i] = 2 / dav;
				}
				// add diffusive mixf...
				m[i] += diffc_here /
					(cell_data[i].length * cell_data[i].length + cell_data[i].length * cell_data[i - 1].length);
				m[i] *= corr_disp;		/* m[i] has mixf with lower numbered cell */
				if (m[i] != m1[i - 1] && !warning && (!dav || m[i] / (2 / dav) > 1.00001))
				{
					warning_msg("Unequal cell-lengths may give mass-balance error, consider using -multi_D");
					warning = true;
				}
			}
			mf12 = m[i] + m1[i];
			if (mf12 > maxmix)
				maxmix = mf12;
		}
		/*
		* Also for boundary cells
		*/
		if (bcon_first == 1)
		{
			m[1] = diffc_here / (cell_data[1].length * cell_data[1].length);
			if (ishift != 0)
				m[1] += cell_data[1].disp / cell_data[1].length;
			mf12 = m[1] + m1[1];
			if (mf12 > maxmix)
				maxmix = mf12;
		}
		if (bcon_last == 1)
		{
			m1[count_cells] = diffc_here / (cell_data[count_cells].length * cell_data[count_cells].length);
			if (ishift != 0)
				m1[count_cells] += cell_data[count_cells].disp / cell_data[count_cells].length;
			mf12 = m[count_cells] + m1[count_cells];
			if (mf12 > maxmix)
				maxmix = mf12;
		}
		/*
		* Find number of mixes
		*/
		if (maxmix == 0)
			l_nmix = 0;
		else
		{
			if (1.5 * maxmix > (double)INT_MAX)
			{
				m = (LDBLE *)free_check_null(m);
				m1 = (LDBLE *)free_check_null(m1);
				char token[MAX_LENGTH];
				sprintf(token, "Calculated number of mixes %g, is beyond program limit,\nERROR: please set implicit true, or decrease time_step, or increase cell-lengths.", 1.5 * maxmix);
				error_msg(token, STOP);
			}
			l_nmix = 1 + (int) floor(1.5 * maxmix);

			if ((ishift != 0) && ((bcon_first == 1) || (bcon_last == 1)))
			{
				if (l_nmix < 2)
					l_nmix = 2;
			}

			for (i = 1; i <= count_cells; i++)
			{
				m[i] /= l_nmix;
				m1[i] /= l_nmix;
				/*
				* Fill mix structure
				*/
				cxxMix temp_mix;
				temp_mix.Set_n_user(i);
				temp_mix.Set_n_user_end(i);

				temp_mix.Add(i - 1, m[i]);
				temp_mix.Add(i + 1, m1[i]);
				temp_mix.Add(i, 1.0 - m[i] - m1[i]);
				Dispersion_mix_map[i] = temp_mix;
			}
		}
		m = (LDBLE *)free_check_null(m);
		m1 = (LDBLE *)free_check_null(m1);
		return (l_nmix);
	}
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
mix_stag(int i, LDBLE kin_time, int l_punch, LDBLE step_fraction)
/* ---------------------------------------------------------------------- */
{
	int n, k;
	LDBLE t_imm;
	cxxSolution *ptr_imm, *ptr_m;
	k = -1000; // compiler says k may be undefined
	ptr_imm = NULL;
	boolean done_mixing = false;
	/*
	* Kinetics in transport cell is done while transporting
	*/
	for (n = 1; n <= stag_data->count_stag; n++)
	{
		if (i == 0 || i == count_cells + 1)
		{
			use.Set_mix_ptr(NULL);
			use.Set_mix_in(false);
			use.Set_mix_ptr(Utilities::Rxn_find(Rxn_mix_map, i));
			if (use.Get_mix_ptr())
			{
				for (std::map < int, LDBLE >::const_iterator it = use.Get_mix_ptr()->Get_mixComps().begin();
					it != use.Get_mix_ptr()->Get_mixComps().end(); it++)
				{
					if (it->first > i && it->first < all_cells)
					{
						k = it->first;
						ptr_imm = Utilities::Rxn_find(Rxn_solution_map, k);
						break;
					}
				}
			}
		}
		else
		{
			k = i + 1 + n * count_cells;
			if (k < all_cells)
				ptr_imm = Utilities::Rxn_find(Rxn_solution_map, k);
		}
		if (ptr_imm != NULL)
		{
			if (n == 1)
			{
				if (heat_nmix > 0 && (!implicit || (implicit && stag_data->count_stag > 1)))
				{
					ptr_m = Utilities::Rxn_find(Rxn_solution_map, i);
					t_imm =
						heat_mix_f_imm * ptr_m->Get_tc() + (1 - heat_mix_f_imm) * ptr_imm->Get_tc();
					ptr_m->Set_tc(heat_mix_f_m * ptr_imm->Get_tc() + (1 - heat_mix_f_m) * ptr_m->Get_tc());
					cell_data[i].temp = ptr_m->Get_tc();
					cell_data[k].temp = t_imm = ptr_imm->Get_tc();
					/* equilibrate again ... */
					cell_no = i;
					set_and_run_wrapper(i, NOMIX, FALSE, i, 0.0);
					if (multi_Dflag == TRUE)
						fill_spec(cell_no, 0);
					saver();
					cell_no = k;
					set_and_run_wrapper(k, NOMIX, FALSE, k, 0.0);
					if (multi_Dflag == TRUE)
						fill_spec(cell_no, i);
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
				if (!implicit || (implicit && stag_data->count_stag > 1))
				{
					if (multi_Dflag == TRUE)
						multi_D(1.0, i, 2);
					set_and_run_wrapper(i, STAG, FALSE, -2, 0.0);
					if (multi_Dflag == TRUE)
						fill_spec(cell_no, 0);
					if (l_punch)
						print_punch(i, true);
					saver(); // save solution i in -2, original can be used in other stagnant mixes
				}
			}

			cell_no = k;
			if (implicit)
				run_reactions(k, kin_time, NOMIX, step_fraction);
			else
				run_reactions(k, kin_time, STAG, step_fraction);
			if (multi_Dflag == TRUE)
				fill_spec(cell_no, i);
			saver(); // save solution k in -2 - k, original k can be used in other stagnant mixes

			done_mixing = true;
		}
		else if (n == 1 && l_punch && !implicit)
			print_punch(i, false);
	}

	if (done_mixing) // after all mixing is done, the temporal solution becomes the original for the next timestep
	{
		for (n = 1; n <= stag_data->count_stag; n++)
		{
			k = i + 1 + n * count_cells;
			if (Utilities::Rxn_find(Rxn_solution_map, k) != 0)
			{
				Utilities::Rxn_copy(Rxn_solution_map, -2 - k, k);
				if (n == 1 && !implicit)
					Utilities::Rxn_copy(Rxn_solution_map, -2, i);
			}
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
init_heat_mix(int l_nmix)
/* ---------------------------------------------------------------------- */
{
	LDBLE lav, mixf, maxmix, corr_disp, l_diffc;
	int i, k, n;
	int l_heat_nmix;
	LDBLE t0;
	/*
	* Check for need to model thermal diffusion...
	*/
	if (heat_diffc <= diffc && !implicit)
		return (0);
	if (count_cells < 2)
		return (0);

	l_heat_nmix = 0;
	if (implicit)
		l_diffc = heat_diffc;
	else
		l_diffc = heat_diffc - diffc_tr;
	t0 = Utilities::Rxn_find(Rxn_solution_map, 0)->Get_tc();
	for (i = 1; i <= count_cells; i++)
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
					if (fabs(cell_data[k].temp - t0) > 1.0)
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
	heat_mix_array = (LDBLE *)PHRQ_malloc((count_cells + 2) * sizeof(LDBLE));
	if (heat_mix_array == NULL)
		malloc_error();

	temp1 = (LDBLE *)PHRQ_malloc((count_cells + 2) * sizeof(LDBLE));
	if (temp1 == NULL)
		malloc_error();

	temp2 = (LDBLE *)PHRQ_malloc((count_cells + 2) * sizeof(LDBLE));
	if (temp2 == NULL)
		malloc_error();
	/*
	* Define mixing factors among inner cells...
	*/
	corr_disp = 1.;
	if (correct_disp == TRUE && ishift != 0)
	{
		mixf = (l_nmix > 1 ? l_nmix : 1);
		if (bcon_first == 3)
			corr_disp += 1. / count_cells / mixf;
		if (bcon_last == 3)
			corr_disp += 1. / count_cells / mixf;
	}
	maxmix = 0.0;
	for (i = 1; i < count_cells; i++)
	{
		lav = (cell_data[i + 1].length + cell_data[i].length) / 2;
		mixf = (l_diffc) * timest * corr_disp / tempr / (lav * lav);
		if (mixf > maxmix)
			maxmix = mixf;
		heat_mix_array[i + 1] = mixf;	/* m[i] has mixf with lower cell */
	}
	/*
	* Also for boundary cells
	*/
	if (bcon_first == 1)
	{
		lav = cell_data[1].length;
		mixf = (l_diffc) * timest * corr_disp / tempr / (lav * lav);
		if (2 * mixf > maxmix)
			maxmix = 2 * mixf;
		heat_mix_array[1] = 2 * mixf;
	}
	else
		heat_mix_array[1] = 0;

	if (bcon_last == 1)
	{
		lav = cell_data[count_cells].length;
		mixf = (l_diffc) * timest * corr_disp / tempr / (lav * lav);
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
		if (implicit)
		{
			LDBLE viscos_f;
			l_heat_nmix = l_nmix;
			for (i = 1; i <= count_cells + 1; i++)
			{
				heat_mix_array[i - 1] = heat_mix_array[i] / l_heat_nmix;	/* for implicit, m[i] has mixf with higher cell */
				viscos_f = sol_D[i - 1].viscos_f * exp(heat_diffc / sol_D[i - 1].tk_x - heat_diffc / 298.15);
				viscos_f += sol_D[i].viscos_f * exp(heat_diffc / sol_D[i].tk_x - heat_diffc / 298.15);
				heat_mix_array[i - 1] *= (viscos_f / 2);
			}
		}
		else
		{
			l_heat_nmix = 1 + (int)floor(3.0 * maxmix);
			for (i = 1; i <= count_cells + 1; i++)
				heat_mix_array[i] /= l_heat_nmix;
		}
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
			heat_mix_array[j] * temp1[j - 1] + heat_mix_array[j + 1] * temp1[j + 1] + 
				(1 - heat_mix_array[j] - heat_mix_array[j + 1]) * temp1[j];
		for (j = 1; j <= count_cells; j++)
			temp1[j] = temp2[j];
	}

	for (i = 1; i <= count_cells; i++)
	{
		cell_data[i].temp = temp1[i];
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
			it = pp_assemblage_ptr->Get_pp_assemblage_comps().begin();
			for (; it != pp_assemblage_ptr->Get_pp_assemblage_comps().end(); it++)
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
fill_spec(int l_cell_no, int ref_cell)
/* ---------------------------------------------------------------------- */
{
	/* copy species activities into sol_D.spec... */

	int i, i1, i2, i3, count_spec, count_exch_spec, size_xt;
	char token[MAX_LENGTH];
	const char * name;
	struct species *s_ptr, *s_ptr2;
	struct master *master_ptr;
	LDBLE dum, dum2;
	LDBLE lm;
	LDBLE por, por_il, viscos_f, viscos_il_f, viscos;
	bool x_max_done = false;
	std::set <std::string> loc_spec_names;

	s_ptr2 = NULL;

	// for implicit
	std::set <std::string> ::iterator it;
	std::pair <std::set<std::string> ::iterator, bool> name_ret;
	if (implicit && !l_cell_no)
		dif_spec_names.clear();
	size_xt = 5;

	//sol_D[l_cell_no].spec = (struct spec *) free_check_null(sol_D[l_cell_no].spec);
	if (sol_D[l_cell_no].spec == NULL)
	{
		sol_D[l_cell_no].spec = (struct spec *) PHRQ_malloc((size_t)(count_species_list + size_xt) * sizeof(struct spec));
		sol_D[l_cell_no].spec_size = count_species_list + size_xt;
	}
	else if (count_species_list + size_xt > sol_D[l_cell_no].spec_size)
	{
		sol_D[l_cell_no].spec = (struct spec *) PHRQ_realloc(sol_D[l_cell_no].spec, (size_t)(count_species_list + size_xt) * sizeof(struct spec));
		sol_D[l_cell_no].spec_size = count_species_list + size_xt;
	}
	if (sol_D[l_cell_no].spec == NULL)
		malloc_error();

	for (i = 0; i < sol_D[l_cell_no].spec_size; i++)
	{
		sol_D[l_cell_no].spec[i].name = NULL;
		sol_D[l_cell_no].spec[i].aq_name = NULL;
		sol_D[l_cell_no].spec[i].type = -1;
		sol_D[l_cell_no].spec[i].a = 0.0;
		sol_D[l_cell_no].spec[i].lm = min_dif_LM;
		sol_D[l_cell_no].spec[i].lg = -0.04;
		sol_D[l_cell_no].spec[i].c = 0.0;
		sol_D[l_cell_no].spec[i].z = 0.0;
		sol_D[l_cell_no].spec[i].Dwt = 0.0;
		sol_D[l_cell_no].spec[i].dw_t = 0.0;
		sol_D[l_cell_no].spec[i].erm_ddl = 0.0;
		sol_D[l_cell_no].count_exch_spec = sol_D[l_cell_no].count_spec = 0;
	}

	sol_D[l_cell_no].tk_x = tk_x;

	viscos_f = viscos_il_f = 1.0;
	if (l_cell_no == 0)
	{
		por = cell_data[1].por;
		por_il = cell_data[1].por_il;
	}
	else if (l_cell_no == count_cells + 1)
	{
		por = cell_data[count_cells].por;
		por_il = cell_data[count_cells].por_il;
	}
	else
	{
		por = cell_data[l_cell_no].por;
		por_il = cell_data[l_cell_no].por_il;
	}
	if (por < multi_Dpor_lim)
		por = viscos_f = 0.0;

	if (por_il < interlayer_Dpor_lim)
		por_il = viscos_il_f = 0.0;
	/*
	* correct diffusion coefficient for temperature and viscosity, D_T = D_298 * Tk * viscos_298 / (298 * viscos)
	*   modify viscosity effect: Dw(TK) = Dw(298.15) * exp(dw_t / TK - dw_t / 298.15), SC data from Robinson and Stokes, 1959
	*/
	viscos = viscos_0;
	/*
	* put temperature factor in por_factor which corrects for porous medium...
	*/
	viscos_f *= tk_x * viscos_0_25 / (298.15 * viscos);
	viscos_il_f *= tk_x * viscos_0_25 / (298.15 * viscos);
	sol_D[l_cell_no].viscos_f = tk_x * viscos_0_25 / (298.15 * viscos);

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
		if (s_ptr->type > HPLUS && 
			!(s_ptr->type == EX && interlayer_Dflag))
			continue;
		//if (s_ptr->type == EX && !interlayer_Dflag)
		//	continue;
		//if (s_ptr->type == SURF)
		//	continue;
		if (i > 0 && strcmp(s_ptr->name, species_list[i - 1].s->name) == 0)
			continue;
		//if (s_ptr == s_h2o)
		//	continue;

		if (s_ptr->type == EX)
		{
			if (s_ptr->lm > min_dif_LM)
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
				sol_D[l_cell_no].spec[count_spec].a = dum2 * pow(10, sol_D[l_cell_no].spec[count_spec].lg);
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
					sol_D[l_cell_no].spec[count_spec].Dwt = default_Dw * viscos_il_f;
				else
				{
					if (s_ptr2->dw_t)
					{
						sol_D[l_cell_no].spec[count_spec].Dwt = s_ptr2->dw *
							exp(s_ptr2->dw_t / 298.15 - s_ptr2->dw_t / tk_x) * viscos_il_f;
						sol_D[l_cell_no].spec[count_spec].dw_t = s_ptr2->dw_t;
					}
					else
						sol_D[l_cell_no].spec[count_spec].Dwt = s_ptr2->dw * viscos_il_f;
				}

				//if (implicit) // && name_ret.second && (l_cell_no > 1 || (l_cell_no == 1 && bcon_first != 2)))
				//{
				//	// name_ret = dif_spec_names.insert(s_ptr->name); // but not in implicit now...
				//	// must fill in the spec in previous cells, order the names, see below for aqueous species...
				//}
				count_exch_spec++;
				count_spec++;
			}
			continue;
		}

		lm = s_ptr->lm;

		if (lm > min_dif_LM)
		{
			if (implicit/* && l_cell_no < count_cells + 2 */)
			{
				name_ret = dif_spec_names.insert(s_ptr->name);
				loc_spec_names.insert(s_ptr->name);
				i2 = 0;
				for (it = dif_spec_names.begin(); it != dif_spec_names.end(); it++)
				{
					if (*it == s_ptr->name)
						break;
					i2++;
				}
				if (i2 > count_spec)
				{
					if (l_cell_no == 0)
					{
						for (i1 = i2; i1 > count_spec; i1--)
						{
							it = dif_spec_names.find(s_ptr->name);
							dif_spec_names.erase(--it);
						}
					}
					else
					{
						// there are species before s_ptr->name that must be included first
						i1 = 0;
						for (it = loc_spec_names.begin(); it != loc_spec_names.end(); it++)
						{
							if (*it == s_ptr->name)
								break;
							i1++;
						}
						i3 = i2 - i1;
						if (i3 + count_spec + 1 > sol_D[l_cell_no].spec_size)
						{
							sol_D[l_cell_no].spec = (struct spec *) PHRQ_realloc(sol_D[l_cell_no].spec, (size_t)(i3 + count_spec + 1 + size_xt) * sizeof(struct spec));
							if (sol_D[l_cell_no].spec == NULL)
								malloc_error();
							sol_D[l_cell_no].spec_size = i3 + count_spec + 1 + size_xt;
						}
						for (; i1 < i2; i1++) // i1 is loop variable 
						{
							memmove(&sol_D[l_cell_no].spec[i1], &sol_D[ref_cell].spec[i1], sizeof(struct spec));
							sol_D[l_cell_no].spec[i1].c = 0.0;
							sol_D[l_cell_no].spec[i1].a = 0.0;
							sol_D[l_cell_no].spec[i1].lm = min_dif_LM;
							sol_D[l_cell_no].spec[i1].lg = -0.04;

							loc_spec_names.insert(sol_D[l_cell_no].spec[i1].name);
							count_spec++;
						}
					}
				}
			}
			if (count_spec >= sol_D[l_cell_no].spec_size)
			{
				sol_D[l_cell_no].spec = (struct spec *) PHRQ_realloc(sol_D[l_cell_no].spec, (size_t)(count_spec + size_xt) * sizeof(struct spec));
				if (sol_D[l_cell_no].spec == NULL)
					malloc_error();
				sol_D[l_cell_no].spec_size = count_spec + size_xt;
			}
			sol_D[l_cell_no].spec[count_spec].name = s_ptr->name;
			sol_D[l_cell_no].spec[count_spec].type = AQ;
			sol_D[l_cell_no].spec[count_spec].c = s_ptr->moles / mass_water_aq_x;
			sol_D[l_cell_no].spec[count_spec].a = under(lm + s_ptr->lg);
			sol_D[l_cell_no].spec[count_spec].lm = lm;
			sol_D[l_cell_no].spec[count_spec].lg = s_ptr->lg;
			sol_D[l_cell_no].spec[count_spec].z = s_ptr->z;
			if (s_ptr->dw == 0)
				sol_D[l_cell_no].spec[count_spec].Dwt = default_Dw * viscos_f;
			else
			{
				if (s_ptr->dw_t)
				{
					sol_D[l_cell_no].spec[count_spec].Dwt = s_ptr->dw *
						exp(s_ptr->dw_t / tk_x - s_ptr->dw_t / 298.15) * viscos_f;
					sol_D[l_cell_no].spec[count_spec].dw_t = s_ptr->dw_t;
				}
				else
					sol_D[l_cell_no].spec[count_spec].Dwt = s_ptr->dw * viscos_f;
			}
			if (correct_Dw)
			{
				calc_SC(); // note that neutral species are corrected as if z = 1, but is viscosity-dependent
				sol_D[l_cell_no].spec[count_spec].Dwt = s_ptr->dw_corr * viscos_f;
			}
			if (l_cell_no <= count_cells + 1 && sol_D[l_cell_no].spec[count_spec].Dwt * pow(por, multi_Dn) > diffc_max)
				diffc_max = sol_D[l_cell_no].spec[count_spec].Dwt * pow(por, multi_Dn);
			sol_D[l_cell_no].spec[count_spec].erm_ddl = s_ptr->erm_ddl;

			if (implicit/* && l_cell_no < count_cells + 2 */)
			{
				if (name_ret.second && l_cell_no)
				{
					// must fill in the spec in previous cells, order the names...
					i3 = 0;
					for (it = dif_spec_names.begin(); it != dif_spec_names.end(); ++it)
					{
						if (it == name_ret.first)
							break;
						i3++;
					}
					for (i1 = 0; i1 < l_cell_no; i1++)
					{
						i2 = sol_D[i1].count_spec + 1;
						if (i2 > sol_D[i1].spec_size)
						{
							sol_D[i1].spec = (struct spec *) PHRQ_realloc(sol_D[i1].spec, (size_t)(i2 + size_xt) * sizeof(struct spec));
							if (sol_D[i1].spec == NULL)
								malloc_error();
							sol_D[i1].spec_size = i2 + size_xt;
						}
						i2--;
						for (; i2 > i3; i2--) // i2 is loop variable
							sol_D[i1].spec[i2] = sol_D[i1].spec[i2 - 1];

						memmove(&sol_D[i1].spec[i2], &sol_D[l_cell_no].spec[i2], sizeof(struct spec));
						sol_D[i1].spec[i2].a = 0.0;
						sol_D[i1].spec[i2].lm = min_dif_LM;
						sol_D[i1].spec[i2].lg = -0.04;
						sol_D[i1].spec[i2].c = 0.0;
						sol_D[i1].count_spec += 1;
					}
				}
			}
			count_spec++;
		}
	}
	sol_D[l_cell_no].count_spec = count_spec;
	sol_D[l_cell_no].count_exch_spec = count_exch_spec;
	if (implicit/* && l_cell_no < count_cells + 2 */ && loc_spec_names.size() < dif_spec_names.size())
	{
		i3 = (int) dif_spec_names.size();
		if (i3 > sol_D[l_cell_no].spec_size)
		{
			sol_D[l_cell_no].spec = (struct spec *) PHRQ_realloc(sol_D[l_cell_no].spec, (size_t)(i3 + size_xt) * sizeof(struct spec));
			if (sol_D[l_cell_no].spec == NULL)
				malloc_error();
			sol_D[l_cell_no].spec_size = i3 + size_xt;
		}
		for (i1 = count_spec; i1 < i3; i1++)
		{
			memmove(&sol_D[l_cell_no].spec[i1], &sol_D[ref_cell].spec[i1], sizeof(struct spec));
			sol_D[l_cell_no].spec[i1].c = 0.0;
			sol_D[l_cell_no].spec[i1].a = 0.0;
			sol_D[l_cell_no].spec[i1].lm = min_dif_LM;
			sol_D[l_cell_no].spec[i1].lg = -0.04;
			sol_D[l_cell_no].count_spec += 1;

			loc_spec_names.insert(sol_D[l_cell_no].spec[i1].name);
			count_spec++;
		}
		//if (loc_spec_names != dif_spec_names)
		//{
		//	error_string = sformatf(
		//		"Species in implicit diffusion in cell %d are %; different from previous cells with %d species.",
		//		l_cell_no, loc_spec_names.size(), dif_spec_names.size());
		//	error_msg(error_string, CONTINUE);
		//}
	}
	if (implicit && l_cell_no/* && l_cell_no < count_cells + 2 */)
	{
		for (i = 0; i < count_spec; i++)
		{
			//name = sol_D[l_cell_no].spec[i].name;
			if (strcmp(sol_D[l_cell_no].spec[i].name, sol_D[ref_cell].spec[i].name))
			{
				error_string = sformatf(
					"Implicit diffusion: species %s in cell %d differs from species %s in previous cells.",
					sol_D[l_cell_no].spec[i].name, l_cell_no, sol_D[ref_cell].spec[i].name);
				error_msg(error_string, CONTINUE);
			}
		}
	}
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
void Phreeqc::
diffuse_implicit(LDBLE DDt, int stagnant)
/* ---------------------------------------------------------------------- */
{
	// The structure comes from example 11.3 in A&P and Appelo 2017, CCR 101, 102: 
	// Ct2 is implicitly solved for the individual species, corrected for electro-neutrality, this always gives Ct2 > 0.
	// Version 3.5.2 includes diffusion in 1 stagnant layer in the implicit calc'n.
	// Transport of aqueous species is summarized into master species.
	// With electro-migration, transport of anions and cations is calculated in opposite directions since (sign) J = - z * dV.
	// Only available moles are transported, thus are > 0, but if concentrations oscillate, 
	// change max_mixf in input file: -implicit true 1 # max_mixf = 1 (default).
	int i, icell, cp, comp;
	// ifirst = (bcon_first == 2 ? 1 : 0); ilast = (bcon_last == 2 ? count_cells - 1 : count_cells);
	int ifirst, ilast;
	int i_1, i0, i1, i2 = 0;
	//double mfr, mfr1, max_b = 0, b, grad, dVc, j_0e, min_dif_M = pow(10, min_dif_LM);
	double mfr, mfr1, grad, dVc, j_0e, min_dif_M = pow(10, min_dif_LM);
	LDBLE dum1, dum2, dum_stag = 0.0, min_mol;

	LDBLE dum = 0;
	//cxxSurfaceCharge * charge_ptr = NULL;
	if (stagnant > 1)
		stagnant = 0;
	cxxSolution *sptr1, *sptr2, *sptr_stag;
	std::vector<int> n_solution;
	std::vector<LDBLE> fraction;
	std::map<int, std::map<std::string, double> > ::iterator it1;
	std::map<std::string, double> ::iterator it2;
	// c= count_cells, c1=(c+1)= end boundary-cell, c_1=(c-1), c2=(c+2)= first stagnant cell, cc1=(c+c+1)= last stagnant cell
	int c = count_cells, c1 = c + 1, c2 = c + 2, c_1 = c - 1, cc = c + stagnant * c, cc1 = cc + 1;

	comp = sol_D[1].count_spec - sol_D[1].count_exch_spec;
	if (heat_nmix)
		comp += 1;

	if (Ct2 == NULL)
		Ct2 = (LDBLE *) PHRQ_malloc((size_t) (count_cells + 2 + stagnant * count_cells) * sizeof(LDBLE));
	if (Ct2 == NULL) malloc_error();
	if (l_tk_x2 == NULL)
		l_tk_x2 = (LDBLE *) PHRQ_malloc((size_t) (count_cells + 2 + stagnant * count_cells) * sizeof(LDBLE));
	if (l_tk_x2 == NULL) malloc_error();

	if (A == NULL)
	{
		A = (LDBLE **)PHRQ_malloc((size_t)(count_cells + 2 + stagnant * count_cells) * sizeof(LDBLE *));
		if (A == NULL) malloc_error();
		for (i = 0; i < count_cells + 2 + stagnant * count_cells; i++)
		{
			if (stagnant)
				A[i] = (LDBLE *)PHRQ_calloc((size_t)(2 * count_cells + 2), sizeof(LDBLE));
			else 
				A[i] = (LDBLE *)PHRQ_malloc((size_t)3 * sizeof(LDBLE));
			if (A[i] == NULL) malloc_error();
		}
	}
	if (LU == NULL)
	{
		LU = (LDBLE **)PHRQ_malloc((size_t)(count_cells + 2 + stagnant * count_cells) * sizeof(LDBLE *));
		if (LU == NULL) malloc_error();
		for (i = 0; i < count_cells + 2 + stagnant * count_cells; i++)
		{
			if (stagnant)
				LU[i] = (LDBLE *)PHRQ_calloc((size_t)(2 * count_cells + 2), sizeof(LDBLE));
			else 
				LU[i] = (LDBLE *)PHRQ_malloc((size_t)3 * sizeof(LDBLE));
			if (LU[i] == NULL) malloc_error();
		}
	}
	if (mixf == NULL)
	{
		mixf = (LDBLE **)PHRQ_malloc((size_t)(count_cells + 2) * sizeof(LDBLE *));
		if (mixf == NULL) malloc_error();
		for (i = 0; i < count_cells + 2; i++)
		{
			mixf[i] = NULL;
		}
	}
	if (stagnant)
	{
		if (mixf_stag == NULL)
		{
			mixf_stag = (LDBLE **)PHRQ_malloc((size_t)(count_cells + 2) * sizeof(LDBLE *));
			if (mixf_stag == NULL) malloc_error();
			for (i = 0; i < count_cells + 2; i++)
			{
				mixf_stag[i] = NULL;
			}
		}
	}
	if (comp + 2 > mixf_comp_size)
	{
		for (i = 0; i < count_cells + 2; i++)
		{
			if (mixf[i] == NULL)
				mixf[i] = (LDBLE *)PHRQ_malloc((size_t)(comp + 2) * sizeof(LDBLE));
			else
				mixf[i] = (LDBLE *)PHRQ_realloc(mixf[i], (size_t)(comp + 2) * sizeof(LDBLE));
			if (mixf[i] == NULL) malloc_error();
			if (stagnant)
			{
				if (mixf_stag[i] == NULL)
					mixf_stag[i] = (LDBLE *)PHRQ_malloc((size_t)(comp + 2) * sizeof(LDBLE));
				else
					mixf_stag[i] = (LDBLE *)PHRQ_realloc(mixf_stag[i], (size_t)(comp + 2) * sizeof(LDBLE));
				if (mixf_stag[i] == NULL) malloc_error();
				for (cp = 0; cp < comp; cp++)
					mixf_stag[i][cp] = 0;
			}
		}
		mixf_comp_size = comp + 2;
	}

	if (dV_dcell)
		find_current = 1;
	// obtain b_ij...
	dummy = default_Dw * pow(multi_Dpor, multi_Dn) * multi_Dpor;
	for (i = 0; i <= count_cells + 1; i++)
	{
		if (!heat_nmix)
		{
			if (i <= count_cells)
				l_tk_x2[i] = (sol_D[i].tk_x + sol_D[i + 1].tk_x) / 2;
			//if (stagnant && i > 0 && i <= c)
			//	l_tk_x2[i + c1] = (sol_D[i].tk_x + sol_D[i + c1].tk_x) / 2;
		}
		if (stagnant)
		{
			use.Set_mix_ptr(Utilities::Rxn_find(Rxn_mix_map, i));
			if (use.Get_mix_ptr() == NULL)
			{
				i1 = 0;
				if (i == 0 && (sptr1 = Utilities::Rxn_find(Rxn_solution_map, c2)) != NULL)
				{
					// the 1st boundary solution diffuses into 1st stagnant cell (=c2) by the ratio immobile_kgw in cell c2 / mobile_kgw in cell 1
					find_J(0, 1, 1, 1, 0);
					ct[c2].kgw = sptr1->Get_mass_water();
					for (cp = 0; cp < comp; cp++)
					{
						if (!heat_nmix || cp < comp - 1)
							mixf_stag[i][cp] = DDt * ct[i].v_m[cp].b_ij * ct[c2].kgw / ct[1].kgw;
						else if (heat_nmix && cp == comp - 1)
							mixf_stag[i][cp] = mixf_stag[i][0] * heat_diffc / tempr / sol_D[i].spec[0].Dwt;
					}
					i1 = c2;
				}
				else if (i == c1 && (sptr1 = Utilities::Rxn_find(Rxn_solution_map, cc1)) != NULL)
				{
					// the end boundary solution diffuses into stagnant cell cc1
					find_J(c, c1, 1, 1, 0);
					ct[cc1].kgw = sptr1->Get_mass_water();
					for (cp = 0; cp < comp; cp++)
					{
						if (!heat_nmix || cp < comp - 1)
							mixf_stag[i][cp] = DDt * ct[c].v_m[cp].b_ij * ct[cc1].kgw / ct[c].kgw;
						else if (heat_nmix && cp == comp - 1)
							mixf_stag[i][cp] = mixf_stag[i][0] * heat_diffc / tempr / sol_D[c].spec[0].Dwt;
					}
					i1 = cc1;
				}
				if (i1)
				{
					find_J(i, i1, 1, 1, 1);
				}
			}
			else
			{
				use.Get_mix_ptr()->Vectorize(n_solution, fraction);
				for (i1 = 0; i1 < (int) use.Get_mix_ptr()->Get_mixComps().size(); i1++)
				{
					if (n_solution[i1] > count_cells && n_solution[i1] <= cc1)
					{
						mfr = fraction[i1] / nmix;
						find_J(i, n_solution[i1], mfr, 1, 1);
						mfr *= (2.0 / dummy);
						for (cp = 0; cp < comp; cp++)
						{
							if (!heat_nmix || cp < comp - 1)
								mixf_stag[i][cp] = ct[i].v_m[cp].b_ij * mfr;
							else if (heat_nmix && cp == comp - 1)
								mixf_stag[i][cp] = mixf_stag[i][0] * heat_diffc / tempr / sol_D[i].spec[0].Dwt;
						}
					}
				}
			}
		}
		if (i <= count_cells)
			find_J(i, i + 1, 1, 1, 0);
	}
	if (dV_dcell)
		find_current = 0;

	ifirst = (bcon_first == 2 ? 1 : 0);
	ilast = (bcon_last == 2 ? count_cells - 1 : count_cells);

	for (cp = 0; cp < comp; cp++)
	{
		for (i = 0; i <= cc1; i++)
		{
			if (heat_nmix && cp == comp - 1)
				Ct2[i] = sol_D[i].tk_x;
			else
				Ct2[i] = sol_D[i].spec[cp].c;
		}
		// fill coefficient matrix A ...
		// boundary cells ...
		if (heat_nmix && cp == comp - 1)
		{
			mfr = mixf[0][cp] = heat_mix_array[0];
			mfr1 = mixf[1][cp] = heat_mix_array[1];
		}
		else
		{
			mfr = mixf[0][cp] = DDt * ct[0].v_m[cp].b_ij;
			mfr1 = mixf[1][cp] = DDt * ct[1].v_m[cp].b_ij;
		}
		if (bcon_first == 2)
		{
			A[1][0] = 0; A[1][1] = 1 + mfr1; A[1][2] = -mfr1;
			if (stagnant)
			{
				A[0][0] = A[c2][c2] = 1;
				if (mixf_stag[1][cp])
				{
					A[1][1] += mixf_stag[1][cp];
					A[1][c2] = A[c2][1] = -mixf_stag[1][cp];
					A[c2][c2] += mixf_stag[1][cp];
				}
			}
			else
			{
				A[0][1] = 1; A[0][2] = 0;
			}
		}
		else
		{
			A[1][0] = -mfr; A[1][1] = 1 + mfr + mfr1; A[1][2] = -mfr1;
			if (stagnant)
			{
				if (dV_dcell)
				{
					A[0][0] = 1 + mfr; A[0][1] = -mfr;
				}
				else
					A[0][0] = 1;
				A[c2][c2] = 1;
				if (mixf_stag[0][cp])
				{
					if (dV_dcell)
					{
						A[0][0] += mixf_stag[0][cp]; A[0][c2] = -mixf_stag[0][cp];
					}
					A[c2][0] = -mixf_stag[0][cp];
					A[c2][c2] += mixf_stag[0][cp];
				}
				if (mixf_stag[1][cp])
				{
					A[1][1] += mixf_stag[1][cp];
					A[1][c2] = A[c2][1] = -mixf_stag[1][cp];
					A[c2][c2] += mixf_stag[1][cp];
				}
			}
			else
			{
				if (dV_dcell)
				{
					A[0][1] = 1 + mfr; A[0][2] = -mfr;
				}
				else
				{
					A[0][1] = 1; A[0][2] = 0;
				}
			}
		}
		if (heat_nmix && cp == comp - 1)
		{
			mfr = mixf[c_1][cp] = heat_mix_array[c_1];
			mfr1 = mixf[count_cells][cp] = heat_mix_array[count_cells];
		}
		else
		{
			mfr = mixf[c_1][cp] = DDt * ct[c_1].v_m[cp].b_ij;
			mfr1 = mixf[count_cells][cp] = DDt * ct[count_cells].v_m[cp].b_ij;
		}
		if (bcon_last == 2)
		{
			if (stagnant)
			{
				A[count_cells][c_1] = -mfr; A[count_cells][count_cells] = 1 + mfr;
				A[c1][c1] = A[cc1][cc1] = 1;
				if (mixf_stag[count_cells][cp])
				{
					A[count_cells][count_cells] += mixf_stag[count_cells][cp];
					A[count_cells][cc1] = A[cc1][count_cells] = -mixf_stag[count_cells][cp];
					A[cc1][cc1] += mixf_stag[count_cells][cp];
				}
			}
			else
			{
				A[count_cells][0] = -mfr; A[count_cells][1] = 1 + mfr; A[count_cells][2] = 0;
				A[c1][0] = 0; A[c1][1] = 1;
			}
		}
		else
		{
			if (stagnant)
			{
				A[count_cells][c_1] = -mfr; A[count_cells][count_cells] = 1 + mfr + mfr1;
				A[count_cells][c1] = -mfr1;
				if (dV_dcell)
				{
					A[c1][count_cells] = -mfr1; A[c1][c1] = 1 + mfr1;
				}
				else
					A[c1][c1] = 1;
				A[cc1][cc1] = 1;
				if (mixf_stag[count_cells][cp])
				{
					A[count_cells][count_cells] += mixf_stag[count_cells][cp];
					A[count_cells][cc1] = A[cc1][count_cells] = -mixf_stag[count_cells][cp];
					A[cc1][cc1] += mixf_stag[count_cells][cp];
				}
				if (mixf_stag[c1][cp])
				{
					if (dV_dcell)
					{
						A[c1][c1] += mixf_stag[c1][cp]; A[c1][cc1] = -mixf_stag[c1][cp];
					}
					A[cc1][c1] = -mixf_stag[c1][cp];
					A[cc1][cc1] += mixf_stag[c1][cp];
				}
			}
			else
			{
				A[count_cells][0] = -mfr; A[count_cells][1] = 1 + mfr + mfr1; A[count_cells][2] = -mfr1;
				if (dV_dcell)
				{
					A[c1][0] = -mfr1; A[c1][1] = 1 + mfr1;
				}
				else
				{
					A[c1][0] = 0; A[c1][1] = 1;
				}
			}
		}
		// inner cells ... 
		for (i = 2; i < count_cells; i++)
		{
			if (heat_nmix && cp == comp - 1)
			{
				mfr = mixf[i - 1][cp] = heat_mix_array[i - 1];
				mfr1 = mixf[i][cp] = heat_mix_array[i];
			}
			else
			{
				mfr = mixf[i - 1][cp] = DDt * ct[i - 1].v_m[cp].b_ij;
				mfr1 = mixf[i][cp] = DDt * ct[i].v_m[cp].b_ij;
			}
			if (stagnant && mixf_stag[i][cp])
			{
				A[i][i - 1] = -mfr; A[i][i] = 1 + mfr + mfr1 + mixf_stag[i][cp]; A[i][i + 1] = -mfr1;
				A[i + c1][i + c1] = 1 + mixf_stag[i][cp];
				A[i][i + c1] = A[i + c1][i] = -mixf_stag[i][cp];
			}
			else
			{
				A[i][0] = -mfr; A[i][1] = 1 + mfr + mfr1; A[i][2] = -mfr1;
			}
		}
		if (stagnant)
			// decompose full A in LU, bypass A[][] = 0...
		{
			LDBLE s;
			LU[0][0] = A[0][0]; LU[1][0] = A[1][0]; LU[c2][0] = A[c2][0];
			LU[0][1] = A[0][1] / LU[0][0]; LU[0][c2] = A[0][c2] / LU[0][0];
			for (i = 1; i <= cc1; i++)
			{
				i_1 = i - 1;
				if (i < c2)
				{
					LU[i][i] = A[i][i] - LU[i][i_1] * LU[i_1][i];
					LU[i + 1][i] = A[i + 1][i];
					LU[i][i + 1] = A[i][i + 1] / LU[i][i];
					if (i == 1)
					{
						LU[c2][1] = A[c2][1] - LU[c2][0] * LU[0][1];
						LU[1][c2] = (A[1][c2] - LU[1][0] * LU[0][c2]) / LU[1][1];
					}
					else
					{
						for (i0 = 0; i0 <= i_1; i0++)
						{
							if (i0 == i_1 && i0 < c)
							{
								LU[c2 + i0][i] = A[c2 + i0][i];
								LU[i][c2 + i0] = A[i][c2 + i0] / LU[i][i];
							}
							else if (i0 < c)
							{
								LU[c2 + i0][i] = -LU[c2 + i0][i_1] * LU[i_1][i];
								LU[i][c2 + i0] = -LU[i][i_1] * LU[i_1][c2 + i0] / LU[i][i];
							}
							if (i > c && i0 == c_1)
							{
								LU[c2 + i0][i] += A[c2 + i0][i];
								LU[i][c2 + i0] += A[i][c2 + i0] / LU[i][i];
							}
						}
					}
				}
				else
				{
					// the i = count_cells + 2 and higher rows are filled from i0 onwards, need to subtract the L*U terms for the ith column of L.
					for (i0 = i; i0 <= cc1; i0++)
					{
						s = 0;
						if (i0 == i)
						{
							i2 = (i == c2 ? 0 : i - c1);
							for (i1 = i2; i1 <= i_1; i1++)
								s += LU[i0][i1] * LU[i1][i];
							LU[i][i] = A[i][i] - s;
							i2 += (i2 == 0 ? 2 : 1);
						}
						else
						{
							for (i1 = i2; i1 <= i_1; i1++)
								s += LU[i0][i1] * LU[i1][i];
							LU[i0][i] = -s;
						}
					}
					// and the ith row of U...
					if (i == cc1)
						break;
					for (i0 = i + 1; i0 <= cc1; i0++)
					{
						s = 0;
						if (i0 == i + 1)
							i2 = i - c;
						for (i1 = i2; i1 <= i_1; i1++)
							s += LU[i][i1] * LU[i1][i0];
						LU[i][i0] = -s / LU[i][i];
						i2++;
					}
				}
			}
			Ct2[0] /= LU[0][0];
			for (i = 1; i <= cc1; i++)
			{
				i_1 = i - 1;
				if (i < c2)
					s = LU[i][i_1] * Ct2[i_1];
				else
				{
					i2 = 0;
					s = 0;
					for (i1 = i2; i1 <= i_1; i1++)
					{
						s += LU[i][i1] * Ct2[i1];
						i2 += (i2 == 0 ? 2 : 1);
					}
				}
				Ct2[i] = (Ct2[i] - s) / LU[i][i];
			}
			for (i = cc; i >= 0; i--)
			{
				s = 0;
				for (i1 = i + 1; i1 <= cc1; i1++)
					s += LU[i][i1] * Ct2[i1];
				Ct2[i] -= s;
			}
		}
		else
		{
			// decompose A in LU : store L in A[..][0..1] and U in A[..][2] ...
			for (i = 1; i <= count_cells + 1; i++)
			{
				A[i - 1][2] = A[i - 1][2] / A[i - 1][1];
				A[i][1] -= A[i][0] * A[i - 1][2];
			}
			// solve Ct2 in A.Ct2 = L.U.Ct2 = Ct1, Ct1 was put in Ct2 ...
			// First, find y in L.y = Ct1, put y in Ct2
			Ct2[0] /= A[0][1];
			for (i = 1; i <= count_cells + 1; i++)
				Ct2[i] = (Ct2[i] - A[i][0] * Ct2[i - 1]) / A[i][1];
			// Now obtain Ct2 in U.Ct2 = 'y' ...
			for (i = count_cells; i >= 0; i--)
				Ct2[i] -= A[i][2] * Ct2[i + 1];
		}
		// Moles transported by concentration gradient from cell [i] to [i + 1] go in tot1,
		//        moles by stagnant exchange from cell [i] to [i + c1] go in tot_stag
		// Correct for electro-neutrality...
		for (i = ifirst; i <= ilast + 1; i++)
		{
			if (heat_nmix && cp == comp - 1)
			{
				l_tk_x2[i] = (Ct2[i + 1] + Ct2[i]) / 2;
				cell_data[i].temp = Ct2[i] - 273.15;
				sptr1 = Utilities::Rxn_find(Rxn_solution_map, i);
				sptr1->Set_tc(Ct2[i] - 273.15);
				if (stagnant && i > 0 && i < c1 && mixf_stag[i][cp])
				{
					i1 = i + c1;
					cell_data[i1].temp = Ct2[i1] - 273.15;
					sptr2 = Utilities::Rxn_find(Rxn_solution_map, i1);
					sptr2->Set_tc(Ct2[i1] - 273.15);
				}
				continue;
			}
			else if (i == ilast + 1 && (!stagnant || !mixf_stag[i][cp]))
				continue;

			if (!cp)
			{
				ct[i].J_ij_sum = 0;
				ct[i].Dz2c_stag = 0;
				if (i < ilast + 1)
					current_cells[i].ele = current_cells[i].dif = 0;
			}
			for (i0 = stagnant; i0 >= 0; i0--)
			{
				i1 = (i0 == 0 ? i + 1 : i == 0 ? c2 : i == c1 ? cc1 : i + c1);
				if (Utilities::Rxn_find(Rxn_solution_map, i1) == NULL)
					continue;
				ct[i].J_ij[cp].name = sol_D[i].spec[cp].name;
				ct[i].J_ij[cp].charge = ct[i].v_m[cp].z;
				grad = (Ct2[i1] - Ct2[i])/* * ct[i].v_m[cp].D*/; // .D has the d{lg(gamma)} / d{lg(molal)} correction
				if (!i0)
					ct[i].J_ij[cp].tot1 = -mixf[i][cp] * grad;
				else
					ct[i].J_ij[cp].tot_stag = -mixf_stag[i][cp] * grad;
				if (ct[i].v_m[cp].z)
				{
					if (i == 0)
						ct[i].v_m[cp].zc = ct[i].v_m[cp].z * Ct2[i1];
					else if (i == ilast)
					{
						if (!i0)
							ct[i].v_m[cp].zc = ct[i].v_m[cp].z * Ct2[i];
						else
							ct[i].v_m[cp].zc = ct[i].v_m[cp].z * Ct2[i1];
					}
					else
						ct[i].v_m[cp].zc = ct[i].v_m[cp].z * (Ct2[i] + Ct2[i1]) / 2;
					if (i0)
					{
						mixf_stag[i][cp] *= ct[i].v_m[cp].zc;
						ct[i].J_ij_sum += ct[i].v_m[cp].z * ct[i].J_ij[cp].tot_stag;
						ct[i].Dz2c_stag -= ct[i].v_m[cp].z * mixf_stag[i][cp];
					}
					else if (i < ilast + 1)
					{
						mixf[i][cp] *= ct[i].v_m[cp].zc;
						current_cells[i].dif += ct[i].v_m[cp].z * ct[i].J_ij[cp].tot1;
						current_cells[i].ele -= ct[i].v_m[cp].z * mixf[i][cp];
					}
				}
			}
		}
		// define element list
		if (cp == 0)
			dif_els_names.clear();
		if (heat_nmix && cp == comp - 1)
		{
			comp -= 1;
			break;
		}
		char * temp_name = string_duplicate(ct[1].J_ij[cp].name);
		char * ptr = temp_name;
		count_elts = 0;
		get_elts_in_species(&ptr, 1);
		free_check_null(temp_name);
		for (int k = 0; k < count_elts; k++)
		{
			if (!strcmp(elt_list[k].elt->name, "X")) continue;
			dif_els_names.insert(elt_list[k].elt->name);
		}
	}
	count_m_s = (int) dif_els_names.size();
	sum_R = sum_Rd = 0;
	for (i = ifirst; i <= ilast; i++)
	{
		ct[i].J_ij_count_spec = comp;
		current_cells[i].R = l_tk_x2[i] / (current_cells[i].ele * F_Re3);
		if (dV_dcell && !fix_current)
		{
			sum_R += current_cells[i].R;
			if (i > 1)
				sum_Rd += (current_cells[0].dif - current_cells[i].dif) * current_cells[i].R;
		}
	}
	if (dV_dcell)
	{
		if (fix_current)
		{
			int sign = (current_x >= 0 ? 1 : -1);
			current_x = sign * fix_current * DDt / F_C_MOL;
			j_0e = current_x - current_cells[ifirst].dif;
		}
		else
		{
			j_0e = (dV_dcell * count_cells - sum_Rd) / sum_R;
			current_x = j_0e + current_cells[ifirst].dif;
		}
	}
	else
	{
		current_x = 1e-10 / F_C_MOL;
		cell_data[ifirst].potV = 0e-8;
		j_0e = current_x - current_cells[ifirst].dif;
	}
	dVc = j_0e * current_cells[ifirst].R;
	cell_data[ifirst + 1].potV = cell_data[ifirst].potV + dVc;
	for (i = ifirst + 1; i < ilast; i++)
	{
		dVc = current_cells[i].R * (current_x - current_cells[i].dif);
		//if (((dV_dcell && (dVc * j_0e > 0)) ||
		//	(dV_dcell > 0 && (cell_data[i].potV + dVc) > (cell_data[count_cells + 1].potV)) ||
		//	(dV_dcell < 0 && (cell_data[i].potV + dVc) < (cell_data[count_cells + 1].potV))))
		if ((dV_dcell && (dVc * j_0e > 0)) ||
			((dV_dcell > 0) && ((cell_data[i].potV + dVc) > cell_data[count_cells + 1].potV)) ||
			((dV_dcell < 0) && ((cell_data[i].potV + dVc) < cell_data[count_cells + 1].potV)))
		{
			dVc = (cell_data[count_cells + 1].potV - cell_data[i].potV) / (count_cells + 1 - i);
		}
		cell_data[i + 1].potV = cell_data[i].potV + dVc;
	}
	if (!dV_dcell || fix_current)
	{
		dVc = current_cells[i].R * (current_x - current_cells[i].dif);
		cell_data[i + 1].potV = cell_data[i].potV + dVc;
	}

	for (cp = 0; cp < comp; cp++)
	{
		if (!ct[ifirst].v_m[cp].z)
			continue;
		for (i = ifirst; i <= ilast + stagnant; i++)
		{
			if (i < ilast + 1)
			{
				dVc = (cell_data[i + 1].potV - cell_data[i].potV) * F_Re3 / l_tk_x2[i];
				ct[i].J_ij[cp].tot1 -= mixf[i][cp] * dVc;
			}
			if (stagnant && ct[i].Dz2c_stag)
			{
				ct[i].J_ij[cp].tot_stag += (mixf_stag[i][cp] * ct[i].J_ij_sum / ct[i].Dz2c_stag);
			}
		}
	}
	current_A = current_x / DDt * F_C_MOL;

	for (i = ifirst; i <= ilast + stagnant + (bcon_last == 2 ? 1 : 0); i++)
	{
		if (i <= ilast + 1)
		{
			// preserve the potentials...
			sptr1 = Utilities::Rxn_find(Rxn_solution_map, i);
			sptr1->Set_potV(cell_data[i].potV);
		}
		// Translate transport of the solute species into master species...
		ct[i].count_m_s = count_m_s;
		if (ct[i].m_s_size == 0 && ct[i].m_s != NULL)
			ct[i].m_s = (struct M_S *) free_check_null(ct[i].m_s);
		if (ct[i].m_s == NULL)
		{
			ct[i].m_s = (struct M_S *) PHRQ_malloc((size_t)(count_m_s + 5) * sizeof(struct M_S));
			ct[i].m_s_size = count_m_s + 5;
		}
		else if (count_m_s > ct[i].m_s_size)
		{
			ct[i].m_s = (struct M_S *) PHRQ_realloc(ct[i].m_s, (size_t)(count_m_s + 5) * sizeof(struct M_S));
			ct[i].m_s_size = count_m_s + 5;
		}
		if (ct[i].m_s == NULL)
			malloc_error();
		std::set <std::string> ::iterator it = dif_els_names.begin();
		for (i1 = 0; i1 < count_m_s; i1++)
		{
			ct[i].m_s[i1].tot1 = ct[i].m_s[i1].tot2 = ct[i].m_s[i1].tot_stag = 0.0;
			ct[i].m_s[i1].charge = 0.0;
			ct[i].m_s[i1].name = (*it).c_str();
			it++;
		}
		fill_m_s(ct[i].J_ij, ct[i].J_ij_count_spec, i, stagnant);
	}
	/*
	* 3. find the solutions in the column, add or subtract the moles...
	*/
	cxxNameDouble::iterator kit;
	int if1, il1, incr;
	for (cp = 0; cp < count_m_s; cp++)
	{
		if (!dV_dcell || dV_dcell * ct[1].m_s[cp].charge <= 0)
		{
			if1 = ifirst; il1 = ilast + 1; incr = 1;
		}
		else
		{
			if1 = ilast; il1 = ifirst - 1; incr = -1;
		}
		for (icell = if1; icell != il1; icell += incr)
		{
			min_mol = min_dif_M * ct[icell].kgw;
			dum1 = dum2 = 0;
			sptr1 = Utilities::Rxn_find(Rxn_solution_map, icell);
			sptr2 = Utilities::Rxn_find(Rxn_solution_map, icell + 1);
			if (stagnant)
			{
				i1 = (icell == 0 ? c1 + 1 : icell == c1 ? cc1 : icell + c1);
				sptr_stag = Utilities::Rxn_find(Rxn_solution_map, i1);
			}
			else
				sptr_stag = NULL;

			if (!strcmp(ct[icell].m_s[cp].name, "H"))
			{
				dummy = ct[icell].m_s[cp].tot1;
				if (dV_dcell || (icell > 0 && icell <= ilast))
					sptr1->Set_total_h(sptr1->Get_total_h() - dummy);
				if (dV_dcell || (icell >= 0 && icell < ilast) || (icell == ilast && bcon_last == 2))
					sptr2->Set_total_h(sptr2->Get_total_h() + dummy);
				if (sptr_stag)
				{
					dummy = ct[icell].m_s[cp].tot_stag;
					if (dV_dcell || (icell > 0 && icell <= ilast))
						sptr1->Set_total_h(sptr1->Get_total_h() - dummy);
					if (icell == c)
					{
						// mix the boundary solution with the stagnant column end-cell
						dummy += ct[icell + 1].m_s[cp].tot_stag;
						if (dV_dcell)
							sptr2->Set_total_h(sptr2->Get_total_h() - ct[icell + 1].m_s[cp].tot_stag);
					}
					sptr_stag->Set_total_h(sptr_stag->Get_total_h() + dummy);
				}
				continue;
			}
			if (!strcmp(ct[icell].m_s[cp].name, "O"))
			{
				dummy = ct[icell].m_s[cp].tot1;
				if (dV_dcell || (icell > 0 && icell <= ilast))
					sptr1->Set_total_o(sptr1->Get_total_o() - dummy);
				if (dV_dcell || (icell >= 0 && icell < ilast) || (icell == ilast && bcon_last == 2))
					sptr2->Set_total_o(sptr2->Get_total_o() + dummy);
				if (sptr_stag)
				{
					dummy = ct[icell].m_s[cp].tot_stag;
					if (dV_dcell || (icell > 0 && icell <= ilast))
						sptr1->Set_total_o(sptr1->Get_total_o() - dummy);
					if (icell == c)
					{
						dummy += ct[icell + 1].m_s[cp].tot_stag;
						if (dV_dcell)
							sptr2->Set_total_o(sptr2->Get_total_o() - ct[icell + 1].m_s[cp].tot_stag);
					}
					sptr_stag->Set_total_o(sptr_stag->Get_total_o() + dummy);
				}
				//if (cp == count_m_s - 1) // transport the charge imbalance (but doesn't improve the change of H2 or O2).
				//{
				//	sptr1->Set_cb(sptr1->Get_cb() + ct[icell].J_ij_sum);
				//	sptr2->Set_cb(sptr2->Get_cb() + ct[icell + 1].J_ij_sum);
				//	if (sptr_stag)
				//		sptr_stag->Set_cb(sptr_stag->Get_cb() + ct[i1].J_ij_sum);
				//}
				continue;
			}
			// see if icell - incr has negative moles, subtract it, so that it is canceled...
			if (icell > 0 && icell <= ilast &&
				(it1 = neg_moles.find(icell - incr)) != neg_moles.end()
				&& (it2 = (els = it1->second).find(ct[icell].m_s[cp].name)) != els.end())
			{
				ct[icell].m_s[cp].tot1 += it2->second;
				neg_moles.erase(it1);
				els.erase(it2);
				neg_moles.insert(std::make_pair(icell - incr, els));
			}
			dum1 = sptr1->Get_totals()[ct[icell].m_s[cp].name];
			if (!dum1)
			{
				dum1 = moles_from_redox_states(sptr1, ct[icell].m_s[cp].name);
				if (dum1)
					sptr1->Get_totals()[ct[icell].m_s[cp].name] = dum1;
			}
			dum2 = sptr2->Get_totals()[ct[icell].m_s[cp].name];
			if (!dum2)
			{
				dum2 = moles_from_redox_states(sptr2, ct[icell].m_s[cp].name);
				if (dum2)
					sptr2->Get_totals()[ct[icell].m_s[cp].name] = dum2;
			}
			if (sptr_stag)
			{
				dum_stag = sptr_stag->Get_totals()[ct[icell].m_s[cp].name];
				if (!dum_stag)
				{
					dum_stag = moles_from_redox_states(sptr_stag, ct[icell].m_s[cp].name);
					if (dum_stag)
						sptr_stag->Get_totals()[ct[icell].m_s[cp].name] = dum_stag;
				}
			}

			// check for negative moles, add moles from the donnan layer when necessary and available...
			if (ct[icell].m_s[cp].tot1 > dum1 &&
				(dV_dcell || (icell > 0 && icell <= ilast)))
			{
				if (!ct[icell].dl_s)
					ct[icell].m_s[cp].tot1 = dum1;
				else
				{
					cxxSurface * s_ptr = Utilities::Rxn_find(Rxn_surface_map, icell);
					if (s_ptr)
					{
						dum1 += moles_from_donnan_layer(s_ptr, ct[icell].m_s[cp].name, ct[icell].m_s[cp].tot1 - dum1 + min_mol);
					}
				}
				sptr1->Get_totals()[ct[icell].m_s[cp].name] = dum1;
			}
			// check for negative moles, in the other cell...
			ct[icell].m_s[cp].tot2 = ct[icell].m_s[cp].tot1;
			if (-ct[icell].m_s[cp].tot2 > dum2 &&
				(dV_dcell || (icell >= 0 && icell < ilast) || (icell == ilast && bcon_last == 2)))
			{
				if (!ct[icell + 1].dl_s)
					ct[icell].m_s[cp].tot2 = -dum2;
				else
				{
					cxxSurface * s_ptr = Utilities::Rxn_find(Rxn_surface_map, icell + 1);
					if (s_ptr)
					{
						dum2 += moles_from_donnan_layer(s_ptr, ct[icell].m_s[cp].name, -(ct[icell].m_s[cp].tot2 + dum2 - min_mol));
					}
					sptr2->Get_totals()[ct[icell].m_s[cp].name] = dum2;
					if (-ct[icell].m_s[cp].tot2 > dum2)
						ct[icell].m_s[cp].tot2 = -dum2;
				}
			}
			if (fabs(ct[icell].m_s[cp].tot2) < fabs(ct[icell].m_s[cp].tot1))
				ct[icell].m_s[cp].tot1 = ct[icell].m_s[cp].tot2;

			if (!cp)
			{
				ct[icell].J_ij_sum = ct[icell + 1].J_ij_sum = 0.0;
				if (sptr_stag)
					ct[i1].J_ij_sum = 0.0;
			}

			if (dV_dcell || (icell > 0 && icell <= ilast))
			{
				dum1 -= ct[icell].m_s[cp].tot1;
				if (sptr_stag)
					dum1 -= ct[icell].m_s[cp].tot_stag;
				sptr1->Get_totals()[ct[icell].m_s[cp].name] = (dum1 > 0 ? dum1 : 0e-16);
				if (dum1 < -min_mol && incr > 0)
				{
					els.clear();
					els.insert(std::make_pair(ct[icell].m_s[cp].name, dum1));
					neg_moles.erase(icell);
					neg_moles.insert(std::make_pair(icell, els));
				}
				else
					ct[icell].J_ij_sum += dum1 * ct[icell].m_s[cp].charge;
			}

			if (dV_dcell || (icell >= 0 && icell < ilast) || (icell == ilast && bcon_last == 2))
			{
				dum2 += ct[icell].m_s[cp].tot1;
				sptr2->Get_totals()[ct[icell].m_s[cp].name] = (dum2 > 0 ? dum2 : 0e-16);
				if (dum2 < -min_mol && incr < 0)
				{
					els.clear();
					els.insert(std::make_pair(ct[icell].m_s[cp].name, dum2));
					neg_moles.erase(icell + 1);
					neg_moles.insert(std::make_pair(icell + 1, els));
				}
				else
					ct[icell + 1].J_ij_sum += dum2 * ct[icell].m_s[cp].charge;
			}

			if (sptr_stag)
			{
				if ((it1 = neg_moles.find(i1)) != neg_moles.end()
					&& (it2 = (els = it1->second).find(ct[icell].m_s[cp].name)) != els.end())
				{
					ct[icell].m_s[cp].tot_stag += it2->second;
					// the negative moles are canceled...
					neg_moles.erase(it1);
					els.erase(it2);
					neg_moles.insert(std::make_pair(i1, els));
				}
				dummy = ct[icell].m_s[cp].tot_stag;
				if (icell == c)
				{
					if (dV_dcell)
					{
						if ((dum = sptr2->Get_totals()[ct[icell + 1].m_s[cp].name] - ct[icell + 1].m_s[cp].tot_stag) > 0)
						{
							sptr2->Get_totals()[ct[icell].m_s[cp].name] = dum;
							dummy += ct[icell + 1].m_s[cp].tot_stag;
						}
						else
						{
							dummy += sptr2->Get_totals()[ct[icell].m_s[cp].name];
							sptr2->Get_totals()[ct[icell].m_s[cp].name] = 0.0;
							if (dum < -min_mol)
							{
								els.clear();
								els.insert(std::make_pair(ct[icell].m_s[cp].name, dum));
								neg_moles.erase(icell);
								neg_moles.insert(std::make_pair(icell, els));
							}
						}
					}
					else
						dummy += ct[icell + 1].m_s[cp].tot_stag;
				}
				if (-dummy > dum_stag)
				{
					if (ct[i1].dl_s)
					{
						cxxSurface * s_ptr = Utilities::Rxn_find(Rxn_surface_map, i1);
						if (s_ptr)
						{
							dum_stag += moles_from_donnan_layer(s_ptr, ct[icell].m_s[cp].name, -(dummy + dum_stag - min_mol));
						}
						sptr_stag->Get_totals()[ct[icell].m_s[cp].name] = dum_stag;
					}
				}
				dum_stag += dummy;
				sptr_stag->Get_totals()[ct[icell].m_s[cp].name] = (dum_stag > 0 ? dum_stag : 0e-16);
				if (dum_stag < -min_mol)
				{
					els.clear();
					els.insert(std::make_pair(ct[icell].m_s[cp].name, dum_stag));
					neg_moles.erase(i1);
					neg_moles.insert(std::make_pair(i1, els));
				}
				else
					ct[i1].J_ij_sum += dum_stag * ct[icell].m_s[cp].charge;
			}
			//if (cp == count_m_s - 1)
			//{
			//	sptr1->Set_cb(sptr1->Get_cb() + ct[icell].J_ij_sum);
			//	sptr2->Set_cb(sptr2->Get_cb() + ct[icell + 1].J_ij_sum);
			//  if (sptr_stag)
			//		sptr_stag->Set_cb(sptr_stag->Get_cb() + ct[i1].J_ij_sum);
			//}
		}
	}
	return;
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
moles_from_redox_states(cxxSolution *sptr, const char *name)
/* ---------------------------------------------------------------------- */
{
	int length, length2;
	cxxNameDouble::iterator kit;
	LDBLE dum = 0;

	length = (int)strlen(name);
	for (kit = sptr->Get_totals().begin(); kit != sptr->Get_totals().end(); kit++)
	{
		length2 = (int)(size_t)strcspn(kit->first.c_str(), "(");
		if (length == length2 && !strncmp(name, kit->first.c_str(), length2))
		{
			dum += kit->second; kit->second = 0;
		}
	}
	return(dum);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
add_MCD_moles(LDBLE moles, LDBLE min_mol, int i, cxxSolution *sptr, const char *name)
/* ---------------------------------------------------------------------- */
{
	cxxNameDouble::iterator kit;
	std::map<int, std::map<std::string, double> > ::iterator it1;
	std::map<std::string, double> ::iterator it2;
	LDBLE dum = sptr->Get_totals()[name];

	if (!dum)
		dum = moles_from_redox_states(sptr, name);
	if ((it1 = neg_moles.find(i)) != neg_moles.end()
		&& (it2 = (els = it1->second).find(name)) != els.end())
	{
		dum += it2->second;
		neg_moles.erase(it1);
		els.erase(it2);
		neg_moles.insert(std::make_pair(i, els));
	}
	dum += moles;
	if (dum < -min_mol)
	{
		if (ct[i].dl_s)
		{
			cxxSurface *su_ptr = Utilities::Rxn_find(Rxn_surface_map, i);
			if (su_ptr != NULL)
				dum += moles_from_donnan_layer(su_ptr, name, -dum + min_mol);
		}
	}
	sptr->Get_totals()[name] = (dum > 0 ? dum : 0e-16);
	if (dum < -min_mol)
	{
		els.insert(std::make_pair(name, dum));
		if ((it1 = neg_moles.find(i)) != neg_moles.end())
			neg_moles.erase(it1);
		neg_moles.insert(std::make_pair(i, els));
	}
	return(dum);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
moles_from_donnan_layer(cxxSurface *sptr, const char *name, LDBLE moles_needed)
/* ---------------------------------------------------------------------- */
{
	cxxNameDouble::iterator kit;
	LDBLE dum = 0;
	cxxSurfaceCharge * charge_ptr = NULL;

	for (size_t j = 0; j < sptr->Get_surface_charges().size(); j++)
	{
		if (sptr->Get_dl_type() == cxxSurface::DONNAN_DL)
		{
			charge_ptr = &(sptr->Get_surface_charges()[j]);
			for (kit = charge_ptr->Get_diffuse_layer_totals().begin(); kit != charge_ptr->Get_diffuse_layer_totals().end(); kit++)
			{
				if (strcmp(kit->first.c_str(), "H") == 0 || strcmp(kit->first.c_str(), "O") == 0)
					continue;
				if (strcmp(kit->first.c_str(), name) == 0)
				{
					if (kit->second > moles_needed)
					{
						dum += moles_needed; kit->second -= moles_needed;
					}
					else
					{
						dum += kit->second;  kit->second = 0;
					}
				}
			}
		}
	}
	return(dum);
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
	* If an electrical field is applied (dV_dcell != 0), the currents j_i = current_cells[i].ele + dif (C * s)
	are calculated for all cells. Then the ele part from cell 0 -> 1 is calculated:
	current_x = (j_0d + j_0e) = j_0 = j_1 = ... = j_i
	j_0e * (R0 + R1 + ...) + (j_0d - j_1d) * R1 + ... + (j_0d - j_id) * Ri = Vtot
	or
	j_0e * Sum_R + Sum_Rd = Vtot.
	Ri = dV_dcell / j_ie, the relative cell resistance.
	Solve j_0e, find (V1 - V0) = j_0e * R0. j_1e = current_x - j_1d, find (V2 - V1) = j_1e * R1, etc.
	*/
	int icell, jcell, i, l, n, length, length2, il_calcs;
	int i1, loop_f_c;
	int first_c, last_c, last_c2 = 0;
	char token[MAX_LENGTH];
	LDBLE mixf, temp;
	LDBLE dVtemp = 0.0;
	if (dV_dcell && stagnant)
	{
		dVtemp = dV_dcell;
		dV_dcell = 0;
	}

	icell = jcell = -1;
	first_c = last_c = -1;
	il_calcs = -1;

	current_x = sum_R = sum_Rd = 0.0;
	if (dV_dcell)
		find_current = loop_f_c = 1; // calculate J_ij once for all cells, loop to dV_dcell2.
	else
		find_current = loop_f_c = 0;

	for (int f_c = 0; f_c <= loop_f_c; f_c++)
	{
		for (n = 0; n <= (stagnant ? stag_data->count_stag : 0); n++) // allow for stagnant cell mixing with higher cells in the layer
			{
			icell = mobile_cell + 1 + n * count_cells;
			if (stagnant)
			{
				if (n == 0)
					icell -= 1;
				else if (mobile_cell == 0 || mobile_cell == count_cells + 1)
					continue;
				/*
				*    find the mix ptr for icell and go along the cells that mix with it
				*/
				use.Set_mix_ptr(Utilities::Rxn_find(Rxn_mix_map, icell));
				if (use.Get_mix_ptr() == NULL)
				{
					if (first_c < 0)
						first_c = icell;
					if (last_c2 < icell)
						last_c2 = icell;
					continue;
				}
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
				last_c2 = (dV_dcell ? count_cells + 1: count_cells);
			}

			for (i = first_c; i <= last_c; i++)
			{
				if (stagnant)
				{
					std::vector<int> n_solution;
					std::vector<LDBLE> fraction;
					(use.Get_mix_ptr())->Vectorize(n_solution, fraction);
					jcell = n_solution[i];
					if (jcell <= icell)
						continue;
					if (jcell >= all_cells || jcell < 0)
					{
						error_string = sformatf(
							"Multi_D asked for diffusion from cell %d to %d, %d is beyond the number of cells", icell, jcell, jcell);
						error_msg(error_string, CONTINUE);
					}
					if (jcell > last_c2)
						last_c2 = jcell;
					mixf = (nmix ? fraction[i] / nmix : fraction[i]);
				}
				else
				{					/* regular column... */
					icell = i;
					jcell = i + 1;
					mixf = 1.0;
				}
				if (dV_dcell)
				{
					tk_x2 = (sol_D[icell].tk_x + sol_D[jcell].tk_x) / 2;
				}
				/*
				* 1. obtain J_ij...
				*/
				il_calcs = (int) find_J(icell, jcell, mixf, DDt, stagnant);

				if (find_current)
				{
					if (i < last_c)
						continue;
					else
					{
						LDBLE dVc, j_0e;
						// distribute dV_dcell according to relative resistance, calculate current_x if not fixed
						j_0e = (dV_dcell * count_cells - sum_Rd) / sum_R;
						current_x = j_0e + current_cells[0].dif;
						if (fix_current)
						{
							int sign = (current_x >= 0 ? 1 : -1);
							current_x = sign * fix_current / F_C_MOL;
							j_0e = current_x - current_cells[0].dif;
						}
						dVc = j_0e * current_cells[0].R;
						cell_data[1].potV = cell_data[0].potV + dVc;
						for (i1 = 1; i1 < count_cells; i1++)
						{
							dVc = current_cells[i1].R * (current_x - current_cells[i1].dif);
							cell_data[i1 + 1].potV = cell_data[i1].potV + dVc;
						}
						if (fix_current)
						{
							dVc = current_cells[i1].R * (current_x - current_cells[i1].dif);
							cell_data[i1 + 1].potV = cell_data[i1].potV + dVc;
						}
						find_current = 0;
						continue;
					}
				}
				if (!ct[icell].J_ij_count_spec)
					continue;

				/*
				* 2. sum up the primary or secondary master_species
				*/
				if (!il_calcs)
				{
					tot1_h = tot1_o = tot2_h = tot2_o = 0.0;
					m_s = (struct M_S *) free_check_null(m_s);
					count_m_s = (ct[icell].J_ij_count_spec < count_moles_added ?
						ct[icell].J_ij_count_spec : count_moles_added);
					m_s = (struct M_S *) PHRQ_malloc((size_t) count_m_s *
						sizeof(struct M_S));
					if (m_s == NULL)
						malloc_error();
					for (i1 = 0; i1 < count_m_s; i1++)
					{
						m_s[i1].name = NULL;
						m_s[i1].tot1 = 0;
						m_s[i1].tot2 = 0;
					}
					count_m_s = 0;
				}
				fill_m_s(ct[icell].J_ij, ct[icell].J_ij_count_spec, icell, stagnant);

				/*
				* 3. find the solutions, add or subtract the moles...
				*/
				if (dV_dcell || (icell != 0 && icell != count_cells + 1))
				{
					use.Set_solution_ptr(Utilities::Rxn_find(Rxn_solution_map, icell));
					use.Get_solution_ptr()->Set_total_h(use.Get_solution_ptr()->Get_total_h() - tot1_h);
					use.Get_solution_ptr()->Set_total_o(use.Get_solution_ptr()->Get_total_o() - tot1_o);
					if (dV_dcell && (icell > 0 || fix_current))
					{
						use.Get_solution_ptr()->Set_potV(cell_data[icell].potV);
					}
					for (l = 0; l < count_m_s; l++)
					{
						length = (int) strlen(m_s[l].name);
						cxxNameDouble::iterator it;
						for (it = use.Get_solution_ptr()->Get_totals().begin();
							it != use.Get_solution_ptr()->Get_totals().end(); it++)
						{
							length2 =
								(int) (size_t) strcspn(it->first.c_str(), "(");
							if (strncmp(m_s[l].name, it->first.c_str(), length) == 0 && length == length2)
							{
								it->second -= m_s[l].tot1;
								break;
							}
						}
						if (it == use.Get_solution_ptr()->Get_totals().end())
						{
							use.Get_solution_ptr()->Get_totals()[m_s[l].name] = -m_s[l].tot1;
						}
					}
				}
				if (dV_dcell || jcell != count_cells + 1)
				{
					use.Set_solution_ptr(Utilities::Rxn_find(Rxn_solution_map, jcell));
					dummy = use.Get_solution_ptr()->Get_total_h();
					use.Get_solution_ptr()->Set_total_h(dummy + tot2_h);
					dummy = use.Get_solution_ptr()->Get_total_o();
					use.Get_solution_ptr()->Set_total_o(dummy + tot2_o);
					if (icell == count_cells && fix_current && !stagnant)
					{
						use.Get_solution_ptr()->Set_potV(cell_data[jcell].potV);
					}
					for (l = 0; l < count_m_s; l++)
					{
						length = (int) strlen(m_s[l].name);
						cxxNameDouble::iterator it;
						for (it = use.Get_solution_ptr()->Get_totals().begin();
							it != use.Get_solution_ptr()->Get_totals().end(); it++)
						{
							length2 = (int) (size_t) strcspn(it->first.c_str(), "(");
							if (strncmp(m_s[l].name, it->first.c_str(), length) == 0 && length == length2)
							{
								it->second += m_s[l].tot2;
								break;
							}
						}
						if (it == use.Get_solution_ptr()->Get_totals().end())
						{
							use.Get_solution_ptr()->Get_totals()[m_s[l].name] = m_s[l].tot2;
						}
					}
				}
			}
		}
	}
	// check for negative conc's...
	//if (stagnant)
	//	first_c = mobile_cell; // allow for stagnant cell mixing with boundary cell 0
	for (i = first_c; i <= last_c2; i++)
	{
		//if (stagnant && i > first_c && i <= count_cells + first_c)
		//	continue;

		use.Set_solution_ptr(Utilities::Rxn_find(Rxn_solution_map, i));
		if (!use.Get_solution_ptr())
			continue;
		cxxNameDouble::iterator it;
		for (it = use.Get_solution_ptr()->Get_totals().begin();
			it != use.Get_solution_ptr()->Get_totals().end(); it++)
		{
			if (strcmp(it->first.c_str(), "H(0)") == 0)
				continue;
			if (strcmp(it->first.c_str(), "O(0)") == 0)
				continue;
			LDBLE moles = it->second;
			if (moles < 0 && ct[i].dl_s)
			{
				use.Set_surface_ptr(Utilities::Rxn_find(Rxn_surface_map, i));
				cxxSurface * s_ptr = use.Get_surface_ptr();
				cxxSurfaceCharge * charge_ptr = NULL;
				cxxNameDouble::iterator jit;
				for (size_t j = 0; j < s_ptr->Get_surface_charges().size(); j++)
				{
					if (s_ptr->Get_dl_type() == cxxSurface::DONNAN_DL)
					{
						charge_ptr = &(s_ptr->Get_surface_charges()[j]);
						for (jit = charge_ptr->Get_diffuse_layer_totals().begin(); jit != charge_ptr->Get_diffuse_layer_totals().end(); jit++)
						{
							if (strcmp(jit->first.c_str(), "H") == 0 || strcmp(jit->first.c_str(), "O") == 0)
								continue;
							if (strcmp(jit->first.c_str(), it->first.c_str()) == 0)
							{
								moles += jit->second;
								it->second += jit->second;
								jit->second = 0;
							}
						}
					}
				}
			}
			if (moles < 0)
			{
				temp = moles;
				it->second = 0;
				/* see if other redox states have more moles... */
				length = (int) strlen(it->first.c_str());
				cxxNameDouble::iterator kit;
				for (kit = use.Get_solution_ptr()->Get_totals().begin();
					kit != use.Get_solution_ptr()->Get_totals().end(); kit++)
				{
					length2 = (int) (size_t) strcspn(kit->first.c_str(), "(");
					if (!strncmp(it->first.c_str(), kit->first.c_str(), length2))
					{
						temp += kit->second;
						if (temp < 0)
						{
							kit->second = 0;
						}
						else
						{
							kit->second = temp;
							break;
						}
					}
				}
				if (temp < -1e-12)
				{
					sprintf(token,
						"Negative concentration in MCD: added %.4e moles %s in cell %d",
						(double)-temp, it->first.c_str(), i);
					warning_msg(token);
					for (i1 = 0; i1 < count_moles_added; i1++)
					{
						if (moles_added[i1].name && !strcmp(moles_added[i1].name, it->first.c_str()))
						{
							moles_added[i1].moles -= temp;
							break;
						}
						else if (!moles_added[i1].moles)
						{
							moles_added[i1].name = string_duplicate(it->first.c_str());
							moles_added[i1].moles -= temp;
							break;
						}
					}
				}
			}
		}
	}

	m_s = (struct M_S *) free_check_null(m_s);

	for (i = first_c; i < last_c2; i++)
	{
		if (stagnant && i > first_c && i <= count_cells + first_c)
			continue;
		ct[i].J_ij = (struct J_ij *) free_check_null(ct[i].J_ij);
		if (il_calcs)
			ct[i].J_ij_il = (struct J_ij *) free_check_null(ct[i].J_ij_il);
		ct[i].v_m = (struct V_M *) free_check_null(ct[i].v_m);
	}
	if (dVtemp && stagnant)
	{
		dV_dcell = dVtemp;
	}

	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
fill_m_s(struct J_ij *l_J_ij, int l_J_ij_count_spec, int icell, int stagnant)
/* ---------------------------------------------------------------------- */
{
	/*  sum up the primary or secondary master_species from solute species
	*      H and O go in tot1&2_h and tot1&2_o
	*/
	int j, k, l;
	LDBLE fraction;
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
		if (implicit && stagnant < 2)
		{
			for (k = 0; k < count_elts; k++)
			{
				for (l = 0; l < count_m_s; l++)
				{
					if (strcmp(ct[icell].m_s[l].name, elt_list[k].elt->name) == 0)
					{
						fraction = fabs((double)elt_list[k].coef * l_J_ij[j].tot1) + fabs(ct[icell].m_s[l].tot1);
						if (fraction)
							fraction = fabs((double)elt_list[k].coef * l_J_ij[j].tot1) / fraction;
						else
							fraction = 1;
						ct[icell].m_s[l].tot1 += elt_list[k].coef * l_J_ij[j].tot1;
						ct[icell].m_s[l].charge *= (1 - fraction);
						ct[icell].m_s[l].charge += fraction * l_J_ij[j].charge;
						if (stagnant)
							ct[icell].m_s[l].tot_stag += elt_list[k].coef * l_J_ij[j].tot_stag;
						break;
					}
				}
			}
		}
		else
		{
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
						m_s[l].name = elt_list[k].elt->name;
						m_s[l].tot1 = elt_list[k].coef * l_J_ij[j].tot1;
						m_s[l].tot2 = elt_list[k].coef * l_J_ij[j].tot2;
						count_m_s++;
					}
				}
			}
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
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
	*    grad(c) is concentration difference in icell and jcell (dx is in lav),
	for activity corrections see Appelo & Wersin, 2007.
	** dec. 28, 2015**
	included aq_dl in the harmonic mean:
	J_ij = - b_i * b_j / (b_i + b_j) * (c_j - c_i) in mol/s (see ex 21 in the manual 3).
	b_i = A1 / (G_i * h_i / 2) * Dw for a pore without EDL. A1 = aq1 / h_i (m^2).
	with EDL (no aq_il_i in A1, for now):
	t_aq1 = aq1 + aq_dl_i. A1 = t_aq1 / h_i. f_free_i = aq1 / t_aq1.
	b_i_cat = A1 / (G_i * h_i / 2) * Dw * {f_free + (1 - f_free) * Bm}. Bm  = Boltzmann enrichment in EDL = g_dl.
	b_i_ani = A1 / (G_i * h_i / 2) * Dw * {f_free + (1 - f_free) / Bm)}.
	22/2/18: now calculates diffusion through EDL's of multiple, differently charged surfaces
	*  stagnant TRUE:
	*    same eqn for J_ij, but multplies with 2 * mixf. (times 2, because mixf = A / (G_i * h_i))
	*    mixf_ij = mixf / (Dw / init_tort_f) / new_tort_f * new_por / init_por
	*    mixf is defined in MIX; Dw is default multicomponent diffusion coefficient;
	*    init_tort_f equals multi_Dpor^(-multi_Dn); new_pf = new tortuosity factor.
	* Interlayer diffusion (IL) takes the gradient in the equivalent concentrations on X-.
	surface area A for IL:
	stagnant: ct[icell].mixf_il is mixf * por_il / por.
	por_il = interlayer porosity, from -interlayer_D true 'por_il'.
	por = free + DL porewater porosity, from -multi_D true 'multi_Dpor'.
	in regular column, A is calc'd from (free + DL porewater) and cell-length.
	for IL: A * por_il / por.

	por_il should be entered for the cell with the maximal cec.
	IL water is related to X-, thus the cec (eq/L IL water) is the same for all cells if X is difined.
	IL-water = (free + DL porewater) * por_il / por.
	for IL: A * aq_il / t_aq.
	*/
	int i, i_max, j, j_max, k, k_il, only_counter, il_calcs;
	int i1;
	LDBLE A1 = 0.0, A2 = 0.0, ddlm, aq1, aq2, t_aq1, t_aq2, f_free_i, f_free_j;
	LDBLE dl_aq1, dl_aq2, dum, dum1, dum2, tort1, tort2, b_i, b_j;
	LDBLE Sum_zM, aq_il1, aq_il2;
	LDBLE por_il1, por_il2, por_il12 = 0.0;
	LDBLE cec1, cec2, cec12 = 0.0, rc1 = 0.0, rc2 = 0.0;
	LDBLE dV, c1, c2;
	LDBLE max_b = 0.0;

	il_calcs = (interlayer_Dflag ? 1 : 0);
	cxxSurface *s_ptr1, *s_ptr2;
	LDBLE g_i, g_j;

	std::vector<cxxSurfaceCharge> s_charge_p;
	std::vector<cxxSurfaceCharge> s_charge_p1;
	std::vector<cxxSurfaceCharge> s_charge_p2;
	std::vector<cxxSurfaceCharge>::iterator it_sc;
	std::vector<cxxSurfaceComp> s_com_p;

	ct[icell].dl_s = ct[jcell].dl_s = dl_aq1 = dl_aq2 = 0.0;

	if (dV_dcell && !find_current)
		goto dV_dcell2;

	/* check for immediate return and interlayer diffusion calcs... */
	//ct[icell].J_ij_sum = 0.0;
	ct[icell].J_ij_count_spec = 0;
	if (!il_calcs)
	{
		if (stagnant)
		{
			if (cell_data[icell].por < multi_Dpor_lim || cell_data[jcell].por < multi_Dpor_lim)
			{
				if (!implicit)
					return (il_calcs);
			}
		}
		else
		{							/* regular column... */
			if ((icell == 0 && cell_data[1].por < multi_Dpor_lim)
				|| (icell == count_cells && cell_data[count_cells].por < multi_Dpor_lim)
				|| (icell != 0 && icell != count_cells && (cell_data[icell].por < multi_Dpor_lim
				|| cell_data[jcell].por < multi_Dpor_lim)))
			{
				if (dV_dcell)
				{
					current_cells[icell].R = -1e15;
					current_cells[icell].ele = dV_dcell / current_cells[icell].R;
					current_cells[icell].dif = 0;
					sum_R += current_cells[icell].R;
					sum_Rd += current_cells[0].dif * current_cells[icell].R;
				}
				if (!implicit)
					return (il_calcs);
			}
		}
	}

	/* do the calcs */
	aq1 = ct[icell].kgw = Utilities::Rxn_find(Rxn_solution_map, icell)->Get_mass_water();
	aq2 = ct[jcell].kgw = Utilities::Rxn_find(Rxn_solution_map, jcell)->Get_mass_water();
	/*
	* check if DL calculations must be made, find amounts of water...
	*/
	s_ptr1 = s_ptr2 = NULL;
	ct[icell].visc1 = ct[icell].visc2 = 1.0;
	only_counter = FALSE;
	
	s_ptr1 = Utilities::Rxn_find(Rxn_surface_map, icell);
	if (s_ptr1 != NULL)
	{
		if (s_ptr1->Get_dl_type() != cxxSurface::NO_DL)
		{
			s_charge_p.assign(s_ptr1->Get_surface_charges().begin(), s_ptr1->Get_surface_charges().end());
			s_com_p.assign(s_ptr1->Get_surface_comps().begin(), s_ptr1->Get_surface_comps().end());

			if (s_ptr1->Get_only_counter_ions())
				only_counter = TRUE;

			ct[icell].visc1 = s_ptr1->Get_DDL_viscosity();
			/* find the immobile surface charges with DL... */
			for (i = 0; i < (int)s_charge_p.size(); i++)
			{
				for (i1 = 0; i1 < (int)s_com_p.size(); i1++)
				{
					if (!(s_charge_p[i].Get_name().compare(s_com_p[i1].Get_charge_name())) && !s_com_p[i1].Get_Dw())
					{
						dl_aq1 += s_charge_p[i].Get_mass_water();
						s_charge_p1.push_back(s_charge_p[i]);
						break;
					}
				}
			}
		}
	}
	s_ptr2 = Utilities::Rxn_find(Rxn_surface_map, jcell);
	if (s_ptr2 != NULL)
	{
		if (s_ptr2->Get_dl_type() != cxxSurface::NO_DL)
		{
			s_charge_p.assign(s_ptr2->Get_surface_charges().begin(), s_ptr2->Get_surface_charges().end());
			s_com_p.assign(s_ptr2->Get_surface_comps().begin(), s_ptr2->Get_surface_comps().end());

			if (s_ptr2->Get_only_counter_ions())
				only_counter = TRUE;

			ct[icell].visc2 = s_ptr2->Get_DDL_viscosity();

			for (i = 0; i < (int)s_charge_p.size(); i++)
			{
				for (i1 = 0; i1 < (int)s_com_p.size(); i1++)
				{
					if (!(s_charge_p[i].Get_name().compare(s_com_p[i1].Get_charge_name())) && !s_com_p[i1].Get_Dw())
					{
						dl_aq2 += s_charge_p[i].Get_mass_water();
						s_charge_p2.push_back(s_charge_p[i]);
						break;
					}
				}
			}
		}
	}
	if (!stagnant)
	{
		if (icell == 0)
			ct[icell].visc1 = ct[icell].visc2;
		else if (icell == count_cells)
			ct[icell].visc2 = ct[icell].visc1;
	}

	/* in each cell: DL surface = mass_water_DL / (cell_length)
	free pore surface = mass_water_free / (cell_length)
	determine DL surface as a fraction of the total pore surface... */
	t_aq1 = aq1 + dl_aq1;
	t_aq2 = aq2 + dl_aq2;
	f_free_i = aq1 / t_aq1;
	f_free_j = aq2 / t_aq2;
	if (dl_aq1 > 0)
		ct[icell].dl_s = dl_aq1 / t_aq1;
	if (dl_aq2 > 0)
		ct[icell].dl_s = ct[jcell].dl_s = dl_aq2 / t_aq2;

	if (il_calcs)
	{
		/* find interlayer porosity por_il,
		make it relative to exchange capacity (mol X/L), highest X in sol_D[1].x_max (mol X / L).
		Find amounts of IL water and cec. */
		por_il1 = por_il2 = por_il12 = 0.0;
		cec1 = cec2 = cec12 = rc1 = rc2 = 0.0;
		if (icell == 0)
		{
			por_il1 = sol_D[0].exch_total / aq1 / sol_D[1].x_max *
				cell_data[1].por_il;
			por_il2 = sol_D[1].exch_total / aq2 / sol_D[1].x_max *
				cell_data[1].por_il;
			if (sol_D[0].exch_total > 3e-10 && sol_D[1].exch_total > 3e-10)
				/* take the harmonic mean... */
				por_il12 = 2 * por_il1 * por_il2 / (por_il1 + por_il2);
			else
				/* at column ends, take the clay... */
				por_il12 = (por_il1 >= por_il2 ? por_il1 : por_il2);

			aq_il2 = t_aq2 * por_il12 / cell_data[1].por;
			aq_il1 = t_aq2 * por_il12 / cell_data[1].por;
		}
		else if (icell == count_cells)
		{
			por_il1 = sol_D[count_cells].exch_total / aq1 / sol_D[1].x_max *
				cell_data[count_cells].por_il;
			por_il2 = sol_D[count_cells + 1].exch_total / aq2 / sol_D[1].x_max *
				cell_data[count_cells].por_il;
			if (sol_D[count_cells].exch_total > 3e-10 && sol_D[count_cells + 1].exch_total > 3e-10)
				por_il12 = 2 * por_il1 * por_il2 / (por_il1 + por_il2);
			else
				por_il12 = (por_il1 >= por_il2 ? por_il1 : por_il2);

			aq_il1 = t_aq1 * por_il12 / cell_data[count_cells].por;
			aq_il2 = t_aq2 * por_il12 / cell_data[count_cells].por;
		}
		else
		{
			por_il1 = sol_D[icell].exch_total / aq1 / sol_D[1].x_max *
				cell_data[icell].por_il;
			por_il2 = sol_D[jcell].exch_total / aq2 / sol_D[1].x_max *
				cell_data[jcell].por_il;

			if (sol_D[icell].exch_total > 3e-10 && sol_D[jcell].exch_total > 3e-10)
				por_il12 = 2 * por_il1 * por_il2 / (por_il1 + por_il2);
			else
				por_il12 = (por_il1 >= por_il2 ? por_il1 : por_il2);

			aq_il1 = t_aq1 * por_il12 / cell_data[icell].por;
			aq_il2 = t_aq2 * por_il12 / cell_data[jcell].por;
		}
		if (por_il12 < interlayer_Dpor_lim)
			il_calcs = 0;
		else
		{
			dum = sol_D[icell].exch_total;
			dum2 = sol_D[jcell].exch_total;
			// the rc's are for distributing the mole transfer (later on) over X and solution.
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
			/* and the largest for calculating the mass transfer since it is usually at the boundary... */
			cec12 = (cec1 > cec2 ? cec1 : cec2);
		}
	}

	/* Find ct[icell].mixf_il for IL diffusion.
	In stagnant calc's, correct mixf by default values, A = por / tort.
	In regular column, A = surface area / (0.5 * x * tort)*/

	tort1 = tort2 = 1.0;
	ct[icell].A_ij_il = ct[icell].mixf_il = 0.0;
	if (stagnant)
	{
		mixf /= (default_Dw * pow(multi_Dpor, multi_Dn) * multi_Dpor);
		if (il_calcs)
			ct[icell].mixf_il = mixf * por_il12 / interlayer_tortf;
	}
	if (icell == 0)
	{
		tort1 = tort2 = pow(cell_data[1].por, -multi_Dn);
		if (stagnant)
			A2 = cell_data[1].por / tort2;
		else
		{
			A2 = t_aq2 / (cell_data[1].length * 0.5 * cell_data[1].length);
			if (il_calcs && !stagnant)
				ct[icell].A_ij_il = A2 * por_il12 / (cell_data[1].por * interlayer_tortf);
			A2 /= tort2;
		}
		A1 = A2;
	}
	else if (icell == count_cells)
	{
		tort1 = tort2 = pow(cell_data[count_cells].por, -multi_Dn);
		if (stagnant)
			A1 = cell_data[count_cells].por / tort1;
		else
		{
			A1 = t_aq1 / (cell_data[count_cells].length * 0.5 * cell_data[count_cells].length);
			if (il_calcs && !stagnant)
				ct[icell].A_ij_il = A1 * por_il12 / (cell_data[count_cells].por * interlayer_tortf);
			A1 /= tort1;
		}
		A2 = A1;
	}
	else
	{
		tort1 = pow(cell_data[icell].por, -multi_Dn);
		tort2 = pow(cell_data[jcell].por, -multi_Dn);
		if (stagnant)
		{
			A1 = cell_data[icell].por / tort1;
			A2 = cell_data[jcell].por / tort2;
		}
		else
		{
			A1 = t_aq1 / (cell_data[icell].length * 0.5 * cell_data[icell].length);
			A2 = t_aq2 / (cell_data[jcell].length * 0.5 * cell_data[jcell].length);
			if (il_calcs && !stagnant)
			{
				dum = A1 * por_il12 / (cell_data[icell].por * interlayer_tortf);
				dum2 = A2 * por_il12 / (cell_data[jcell].por * interlayer_tortf);
				ct[icell].A_ij_il = dum * dum2 / (dum + dum2);
			}
			A1 /= tort1;
			A2 /= tort2;
		}
	}
	/* diffuse... */
	/*
	* malloc sufficient space...
	*/
	k = sol_D[icell].count_spec + sol_D[jcell].count_spec;

	if (ct[icell].J_ij == NULL)
	{
		ct[icell].J_ij = (struct J_ij *) PHRQ_malloc((size_t)(k) * sizeof(struct J_ij));
		ct[icell].J_ij_size = k;
	}
	else if (k > ct[icell].J_ij_size)
	{
		ct[icell].J_ij = (struct J_ij *) PHRQ_realloc(ct[icell].J_ij, (size_t)(k) * sizeof(struct J_ij));
		ct[icell].J_ij_size = k;
	}
	if (ct[icell].J_ij == NULL)
		malloc_error();

	if (ct[icell].v_m == NULL)
	{
		ct[icell].v_m = (struct V_M *) PHRQ_malloc((size_t)(k) * sizeof(struct V_M));
		ct[icell].v_m_size = k;
	}
	else if (k > ct[icell].v_m_size)
	{
		ct[icell].v_m = (struct V_M *) PHRQ_realloc(ct[icell].v_m, (size_t)(k) * sizeof(struct V_M));
		ct[icell].v_m_size = k;
	}
	if (ct[icell].v_m == NULL)
		malloc_error();

	for (i = 0; i < k; i++)
	{
		ct[icell].J_ij[i].tot1 = ct[icell].J_ij[i].tot_stag = 0.0;
		ct[icell].v_m[i].grad = 0.0;
		ct[icell].v_m[i].D = 1.0; // used for gamma correction
		ct[icell].v_m[i].z = 0.0;
		ct[icell].v_m[i].c = 0.0;
		ct[icell].v_m[i].zc = 0.0;
		//ct[icell].v_m[i].Dz = 0.0;
		//ct[icell].v_m[i].Dzc = 0.0;
		ct[icell].v_m[i].b_ij = 0.0;
	}
	ct[icell].Dz2c = ct[icell].Dz2c_il = 0.0;

	if (il_calcs)
	{
		/* also for interlayer cations */
		k = sol_D[icell].count_exch_spec + sol_D[jcell].count_exch_spec;

		ct[icell].J_ij_il = (struct J_ij *) free_check_null(ct[icell].J_ij_il);
		ct[icell].J_ij_il = (struct J_ij *) PHRQ_malloc((size_t) k * sizeof(struct J_ij));
		if (ct[icell].J_ij_il == NULL)
			malloc_error();

		ct[icell].v_m_il = (struct V_M *) free_check_null(ct[icell].v_m_il);
		ct[icell].v_m_il = (struct V_M *) PHRQ_malloc((size_t) k * sizeof(struct V_M));
		if (ct[icell].v_m_il == NULL)
			malloc_error();

		for (i = 0; i < k; i++)
		{
			ct[icell].J_ij_il[i].tot1 = 0.0;
			ct[icell].v_m_il[i].grad = 0.0;
			ct[icell].v_m_il[i].D = 0.0;
			ct[icell].v_m_il[i].z = 0.0;
			ct[icell].v_m_il[i].c = 0.0;
			ct[icell].v_m_il[i].zc = 0.0;
			ct[icell].v_m_il[i].Dz = 0.0;
			ct[icell].v_m_il[i].Dzc = 0.0;
			ct[icell].v_m_il[i].b_ij = 0.0;
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
			|| (i < i_max && strcmp(sol_D[icell].spec[i].name, sol_D[jcell].spec[j].name) < 0))
		{
			/* species 'name' is only in icell */
			if (il_calcs && sol_D[icell].spec[i].type == EX)
			{
				ct[icell].J_ij_il[k_il].name = sol_D[icell].spec[i].name;
				ct[icell].v_m_il[k_il].D = sol_D[icell].spec[i].Dwt;
				ct[icell].v_m_il[k_il].z = sol_D[icell].spec[i].z;
				ct[icell].v_m_il[k_il].Dz = ct[icell].v_m_il[k_il].D * ct[icell].v_m_il[k_il].z;
				dum = sol_D[icell].spec[i].c * cec12 / (2 * ct[icell].v_m_il[k_il].z);
				ct[icell].v_m_il[k_il].Dzc = ct[icell].v_m_il[k_il].Dz * dum;
				ct[icell].Dz2c_il += ct[icell].v_m_il[k_il].Dzc * ct[icell].v_m_il[k_il].z;
				ct[icell].v_m_il[k_il].grad = -sol_D[icell].spec[i].c * cec12 / ct[icell].v_m_il[k_il].z;	/* use equivalent fraction */
				k_il++;
			}
			else
			{
				ct[icell].J_ij[k].name = sol_D[icell].spec[i].name;
				ct[icell].v_m[k].z = sol_D[icell].spec[i].z;
				ct[icell].v_m[k].grad = -sol_D[icell].spec[i].c; /* assume d log(gamma) / d log(c) = 0 */
				c1 = sol_D[icell].spec[i].c / 2;
				ct[icell].v_m[k].c = c1;
				if (ct[icell].v_m[k].z)
					ct[icell].v_m[k].zc = ct[icell].v_m[k].z * c1;

				if (dV_dcell && ct[icell].v_m[k].z && !fix_current)
				{
					// compare diffusive and electromotive forces
					dum = ct[icell].v_m[k].grad;
					if (icell == 0)
						dum2 = (cell_data[1].potV - cell_data[0].potV) / (cell_data[1].length / 2);
					else if (icell == count_cells)
						dum2 = (cell_data[count_cells + 1].potV - cell_data[count_cells].potV) / (cell_data[count_cells].length / 2);
					else
						dum2 = (cell_data[jcell].potV - cell_data[icell].potV) / ((cell_data[jcell].length + cell_data[icell].length) / 2);
					dum2 *= F_Re3 / tk_x2 * ct[icell].v_m[k].z * c1;
					if (dum + dum2 > 0)
					{
						// step out: no transport against the dV_dcell gradient if c = 0 in jcell...
						if (i < i_max)
							i++;
						continue;
					}
				}

				g_i = g_j = 0;
				if (ct[icell].dl_s > 0)
				{
					if (dl_aq1)
					{
						for (it_sc = s_charge_p1.begin(); it_sc != s_charge_p1.end(); it_sc++)
						{
							g_i += it_sc->Get_z_gMCD_map()[ct[icell].v_m[k].z];
						}
						g_i *= sol_D[icell].spec[i].erm_ddl;
					}
					if (dl_aq2)
					{
						for (it_sc = s_charge_p2.begin(); it_sc != s_charge_p2.end(); it_sc++)
						{
							if (ct[icell].v_m[k].z == 0 || only_counter)
							{
								g_j += it_sc->Get_z_gMCD_map()[ct[icell].v_m[k].z];
							}
							else
							{
								if (abs(ct[icell].v_m[k].z) == 1)
									// there is always H+ and OH-...
									g_j += it_sc->Get_z_gMCD_map()[ct[icell].v_m[k].z];
								else
								{
									dum1 = it_sc->Get_mass_water() / mass_water_bulk_x;
									dum2 = it_sc->Get_z_gMCD_map()[1] / dum1;
									g_j += pow(dum2, ct[icell].v_m[k].z) * dum1;
								}
							}
						}
						g_j *= sol_D[icell].spec[i].erm_ddl;
					}
				}

				b_i = A1 * sol_D[icell].spec[i].Dwt * (f_free_i + g_i / ct[icell].visc1);
				b_j = A2 * (f_free_j + g_j / ct[icell].visc2);
				if (icell == count_cells && !stagnant)
					ct[icell].v_m[k].b_ij = b_i;
				else if (icell == all_cells - 1 && stagnant)
					ct[icell].v_m[k].b_ij = b_i / 2; /* with the mixf *= 2 for this 'reservoir' cell in the input */
				else
				{
					if (sol_D[icell].tk_x == sol_D[jcell].tk_x)
						b_j *= sol_D[icell].spec[i].Dwt;
					else
					{
						dum2 = sol_D[icell].spec[i].Dwt / sol_D[icell].viscos_f;
						dum2 *= exp(sol_D[icell].spec[i].dw_t / sol_D[jcell].tk_x - sol_D[icell].spec[i].dw_t / sol_D[icell].tk_x);
						dum2 *= sol_D[jcell].viscos_f;
						b_j *= dum2;
					}
					ct[icell].v_m[k].b_ij = b_i * b_j / (b_i + b_j);
					if (icell == 0 && !stagnant)
							ct[icell].v_m[k].b_ij = b_j;
					else if (icell == 3 && stagnant && !g_i && g_j)
						ct[icell].v_m[k].b_ij = b_j / 2; /* with the mixf *= 2 for stagnant cell 3 in the input */
				}

				if (ct[icell].v_m[k].z)
					ct[icell].Dz2c += ct[icell].v_m[k].b_ij * ct[icell].v_m[k].zc * ct[icell].v_m[k].z;

				k++;
			}
			if (i < i_max)
				i++;
		}

		else if (i == i_max ||
			(j < j_max && strcmp(sol_D[icell].spec[i].name, sol_D[jcell].spec[j].name) > 0))
		{
			/* species 'name' is only in jcell */
			if (il_calcs && sol_D[jcell].spec[j].type == EX)
			{
				ct[icell].J_ij_il[k_il].name = sol_D[jcell].spec[j].name;
				ct[icell].v_m_il[k_il].D = sol_D[jcell].spec[j].Dwt;
				ct[icell].v_m_il[k_il].z = sol_D[jcell].spec[j].z;
				ct[icell].v_m_il[k_il].Dz = ct[icell].v_m_il[k_il].D * ct[icell].v_m_il[k_il].z;
				ct[icell].v_m_il[k_il].Dzc = ct[icell].v_m_il[k_il].Dz * sol_D[jcell].spec[j].c *
					cec12 / (2 * ct[icell].v_m_il[k_il].z);
				ct[icell].Dz2c_il += ct[icell].v_m_il[k_il].Dzc * ct[icell].v_m_il[k_il].z;
				ct[icell].v_m_il[k_il].grad = sol_D[jcell].spec[j].c * cec12 / ct[icell].v_m_il[k_il].z;	/* use equivalent fraction */
				k_il++;
			}
			else
			{
				ct[icell].J_ij[k].name = sol_D[jcell].spec[j].name;
				ct[icell].v_m[k].z = sol_D[jcell].spec[j].z;
				ct[icell].v_m[k].grad = sol_D[jcell].spec[j].c;  /* assume d log(gamma) / d log(c) = 0 */
				c2 = sol_D[jcell].spec[j].c / 2;
				ct[icell].v_m[k].c = c2;
				if (ct[icell].v_m[k].z)
					ct[icell].v_m[k].zc = ct[icell].v_m[k].z * c2;

				if (dV_dcell && ct[icell].v_m[k].z && !fix_current)
				{
					// compare diffuse and electromotive forces
					dum = ct[icell].v_m[k].grad;
					if (icell == 0)
						dum2 = (cell_data[1].potV - cell_data[0].potV) / (cell_data[1].length / 2);
					else if (icell == count_cells)
						dum2 = (cell_data[count_cells + 1].potV - cell_data[count_cells].potV) / (cell_data[count_cells].length / 2);
					else
						dum2 = (cell_data[jcell].potV - cell_data[icell].potV) / ((cell_data[jcell].length + cell_data[icell].length) / 2);
					dum2 *= F_Re3 / tk_x2 * ct[icell].v_m[k].z * c2;
					// don't transport unavailable moles against the gradient
					if (dum + dum2 < 0)
					{
						// step out...
						if (j < j_max)
							j++;
						continue;
					}
				}
				g_i = g_j = 0;
				if (ct[icell].dl_s > 0)
				{
					if (dl_aq1)
					{
						for (it_sc = s_charge_p1.begin(); it_sc != s_charge_p1.end(); it_sc++)
						{
							if (ct[icell].v_m[k].z == 0 || only_counter)
							{
								g_i += it_sc->Get_z_gMCD_map()[ct[icell].v_m[k].z];
							}
							else
							{
								if (abs(ct[icell].v_m[k].z) == 1)
									// there is always H+ and OH-...
									g_i += it_sc->Get_z_gMCD_map()[ct[icell].v_m[k].z];
								else
								{
									dum1 = it_sc->Get_mass_water() / mass_water_bulk_x;
									dum2 = it_sc->Get_z_gMCD_map()[1] / dum1;
									g_i += pow(dum2, ct[icell].v_m[k].z) * dum1;
								}
							}
						}
						g_i *= sol_D[jcell].spec[j].erm_ddl;
					}
					if (dl_aq2)
					{
						for (it_sc = s_charge_p2.begin(); it_sc != s_charge_p2.end(); it_sc++)
						{
							g_j += it_sc->Get_z_gMCD_map()[ct[icell].v_m[k].z];
						}
						g_j *= sol_D[jcell].spec[j].erm_ddl;
					}
				}
				b_i = A1 * (f_free_i + g_i / ct[icell].visc1);
				b_j = A2 * sol_D[jcell].spec[j].Dwt * (f_free_j + g_j / ct[icell].visc2);
				if (icell == 0 && !stagnant)
					ct[icell].v_m[k].b_ij = b_j;
				else if (icell == 3 && stagnant && g_j && !g_i)
					ct[icell].v_m[k].b_ij = b_j / 2; /* with the mixf *= 2 for 'reservoir' cell 3 in the input */
				else
				{
					if (sol_D[icell].tk_x == sol_D[jcell].tk_x)
						b_i *= sol_D[jcell].spec[j].Dwt;
					else
					{
						dum2 = sol_D[jcell].spec[j].Dwt / sol_D[jcell].viscos_f;
						dum2 *= exp(sol_D[jcell].spec[j].dw_t / sol_D[icell].tk_x - sol_D[jcell].spec[j].dw_t / sol_D[jcell].tk_x);
						dum2 *= sol_D[icell].viscos_f;
						b_i *= dum2;
					}
					ct[icell].v_m[k].b_ij = b_i * b_j / (b_i + b_j);
					if (icell == count_cells && !stagnant)
						ct[icell].v_m[k].b_ij = b_i;
					else if (jcell == all_cells - 1 && stagnant && !g_j && g_i)
						ct[icell].v_m[k].b_ij = b_i / 2; /* with the mixf * 2 for this 'reservoir' cell in the input */
				}
				if (ct[icell].v_m[k].z)
					ct[icell].Dz2c += ct[icell].v_m[k].b_ij * ct[icell].v_m[k].zc * ct[icell].v_m[k].z;

				k++;
			}
			if (j < j_max)
				j++;
		}
		else if (strcmp(sol_D[icell].spec[i].name, sol_D[jcell].spec[j].name) == 0)
		{
			/* species 'name' is in both cells */
			if (il_calcs && sol_D[icell].spec[i].type == EX)
			{
				ct[icell].J_ij_il[k_il].name = sol_D[icell].spec[i].name;
				if (sol_D[icell].spec[i].Dwt == 0 || sol_D[jcell].spec[j].Dwt == 0)
					ct[icell].v_m_il[k_il].D = 0.0;
				else
					ct[icell].v_m_il[k_il].D =
					(sol_D[icell].spec[i].Dwt + sol_D[jcell].spec[j].Dwt) / 2;

				ct[icell].v_m_il[k_il].z = sol_D[icell].spec[i].z;
				ct[icell].v_m_il[k_il].Dz = ct[icell].v_m_il[k_il].D * ct[icell].v_m_il[k_il].z;
				ct[icell].v_m_il[k_il].Dzc = ct[icell].v_m_il[k_il].Dz * (sol_D[icell].spec[i].c +
					sol_D[jcell].spec[j].c) * cec12 / ct[icell].v_m_il[k_il].z;
				ct[icell].Dz2c_il += ct[icell].v_m_il[k_il].Dzc * ct[icell].v_m_il[k_il].z;
				ct[icell].v_m_il[k_il].grad = (sol_D[jcell].spec[j].c - sol_D[icell].spec[i].c) *
					cec12 / ct[icell].v_m_il[k_il].z;	/* use equivalent fraction */
				k_il++;
			}
			else
			{
				ct[icell].J_ij[k].name = sol_D[icell].spec[i].name;
				ct[icell].v_m[k].z = sol_D[icell].spec[i].z;
				ct[icell].v_m[k].grad = (sol_D[jcell].spec[j].c - sol_D[icell].spec[i].c);
				c1 = sol_D[icell].spec[i].c / 2;
				c2 = sol_D[jcell].spec[j].c / 2;
				ct[icell].v_m[k].c = c1 + c2;
				if (ct[icell].v_m[k].z)
					ct[icell].v_m[k].zc = ct[icell].v_m[k].z * ct[icell].v_m[k].c;

				if (dV_dcell && ct[icell].v_m[k].z && !fix_current && !implicit)
				{
					// compare diffuse and electromotive forces
					dum = ct[icell].v_m[k].grad;
					if (icell == 0)
						dum2 = (cell_data[1].potV - cell_data[0].potV) / (cell_data[1].length / 2);
					else if (icell == count_cells)
						dum2 = (cell_data[count_cells + 1].potV - cell_data[count_cells].potV) / (cell_data[count_cells].length / 2);
					else
						dum2 = (cell_data[jcell].potV - cell_data[icell].potV) / ((cell_data[jcell].length + cell_data[icell].length) / 2);
					dum2 *= F_Re3 / tk_x2 * ct[icell].v_m[k].z * (c1 + c2);
					// don't transport unavailable moles against the gradient
					if (fabs(dum) < fabs(dum2) &&
						((dum2 >= 0 && sol_D[jcell].spec[j].c * aq2 < 1e-12) ||
						(dum2 <= 0 && sol_D[icell].spec[i].c * aq1 < 1e-12)))
					{
						// step out:
						if (i < i_max)
							i++;
						if (j < j_max)
							j++;
						continue;
					}
				}
				g_i = g_j = 0;
				if (ct[icell].dl_s > 0)
				{
					if (dl_aq1)
					{
						for (it_sc = s_charge_p1.begin(); it_sc != s_charge_p1.end(); it_sc++)
						{
							g_i += it_sc->Get_z_gMCD_map()[ct[icell].v_m[k].z];
						}
						g_i *= sol_D[icell].spec[i].erm_ddl;
					}
					if (dl_aq2)
					{
						for (it_sc = s_charge_p2.begin(); it_sc != s_charge_p2.end(); it_sc++)
						{
							g_j += it_sc->Get_z_gMCD_map()[ct[icell].v_m[k].z];
						}
						g_j *= sol_D[jcell].spec[j].erm_ddl;
					}
				}
				b_i = A1 * sol_D[icell].spec[i].Dwt * (f_free_i + g_i / ct[icell].visc1);
				b_j = A2 * sol_D[jcell].spec[j].Dwt * (f_free_j + g_j / ct[icell].visc2);
				ct[icell].v_m[k].b_ij = b_i * b_j / (b_i + b_j);
				// but for boundary cells...
				if (stagnant > 1)
				{ /* for a diffusion experiment with well-mixed reservoir in cell 3 and the last stagnant cell,
					   and with the mixf * 2 for the boundary cells in the input... */
					if (icell == 3 && !g_i && g_j)
						ct[icell].v_m[k].b_ij = b_j / 2;
					else if (jcell == all_cells - 1 && !g_j && g_i)
						ct[icell].v_m[k].b_ij = b_i / 2;
				}
				else
				{
					if (icell == 0 || (icell == count_cells + 1 && jcell == count_cells + count_cells + 1))
						ct[icell].v_m[k].b_ij = b_j;
					else if (icell == count_cells && jcell == count_cells + 1)
						ct[icell].v_m[k].b_ij = b_i;
				}

				if (ct[icell].v_m[k].z)
					ct[icell].Dz2c += ct[icell].v_m[k].b_ij * ct[icell].v_m[k].zc * ct[icell].v_m[k].z;

				//ddlm = sol_D[jcell].spec[j].lm - sol_D[icell].spec[i].lm; // appt: this could give an incorrect large factor for implicit
				//if (fabs(ddlm) > 1e-10)
				//	ct[icell].v_m[k].grad *= (1 + (sol_D[jcell].spec[j].lg - sol_D[icell].spec[i].lg) / ddlm);
				ddlm = sol_D[jcell].spec[j].lm - sol_D[icell].spec[i].lm;
				dum1 = sol_D[jcell].spec[j].lg - sol_D[icell].spec[i].lg;
				dum = 1;
				if (ddlm)
					dum += dum1 / ddlm;
				if (dum > 0.2 && dum < 2.25)
				{
					ct[icell].v_m[k].grad *= dum;
					if (implicit)
						ct[icell].v_m[k].D = dum;
				}
				//if (ct[icell].v_m[k].b_ij > max_b)
				//	max_b = ct[icell].v_m[k].b_ij;
				k++;
			}
			if (i < i_max)
				i++;
			if (j < j_max)
				j++;
		}
	}
	if (!dV_dcell && !ct[icell].Dz2c)
		k = 0;
	ct[icell].J_ij_count_spec = i_max = k;
	ct[icell].J_ij_il_count_spec = k_il;
	if (icell == count_cells)
		ct[jcell].J_ij_count_spec = ct[icell].J_ij_count_spec;

	if (implicit && stagnant < 2)
		return(max_b);
	/*
	* fill in J_ij...
	*/
	if (dV_dcell)
	{
		current_cells[icell].ele = current_cells[icell].dif = 0;
		for (i = 0; i < ct[icell].J_ij_count_spec; i++)
		{
			if (!ct[icell].v_m[i].z)
				continue;
			current_cells[icell].ele -= ct[icell].v_m[i].b_ij * ct[icell].v_m[i].z *
				ct[icell].v_m[i].zc;
			current_cells[icell].dif -= ct[icell].v_m[i].b_ij * ct[icell].v_m[i].z *
				ct[icell].v_m[i].grad;
		}
		current_cells[icell].R = tk_x2 / (current_cells[icell].ele * F_Re3);
		sum_R += current_cells[icell].R;
		sum_Rd += (current_cells[0].dif - current_cells[icell].dif) * current_cells[icell].R;
		return(il_calcs);
	}

	// dV_dcell, 2nd pass, current_x has been calculated, and
	// voltage was adapted to give equal current in the cells.
dV_dcell2:

	//ct[icell].J_ij_sum = 0;
	Sum_zM = 0.0;

	for (i = 0; i < ct[icell].J_ij_count_spec; i++)
	{
		if (ct[icell].v_m[i].z)
			Sum_zM += ct[icell].v_m[i].b_ij * ct[icell].v_m[i].z * ct[icell].v_m[i].grad;
	}
	for (i = 0; i < ct[icell].J_ij_count_spec; i++)
	{
		ct[icell].J_ij[i].tot1 = -ct[icell].v_m[i].grad;
		ct[icell].J_ij[i].charge = ct[icell].v_m[i].z;
		if (!dV_dcell && ct[icell].v_m[i].z && ct[icell].Dz2c > 0)
		{
			ct[icell].J_ij[i].tot1 += Sum_zM * ct[icell].v_m[i].zc / ct[icell].Dz2c;
		}
		if (stagnant)
			ct[icell].J_ij[i].tot1 *= ct[icell].v_m[i].b_ij * 2 * mixf;
		else
			ct[icell].J_ij[i].tot1 *= ct[icell].v_m[i].b_ij * DDt;
		ct[icell].J_ij[i].tot2 = ct[icell].J_ij[i].tot1;
		//ct[icell].J_ij_sum += ct[icell].v_m[i].z * ct[icell].J_ij[i].tot1;
	}
	// assure that icell and jcell have dl water when checking negative conc's in MCD
	ct[icell].dl_s = dl_aq1;
	ct[jcell].dl_s = dl_aq2;

	if (dV_dcell)
	{
		//ct[icell].J_ij_sum = 0;
		dV = cell_data[jcell].potV - cell_data[icell].potV;
		dum = dV * F_Re3 / tk_x2;
		// perhaps adapt dV for getting equal current...
		//current_cells[icell].ele = current_cells[icell].dif = 0;
		//for (i = 0; i < ct[icell].J_ij_count_spec; i++)
		//{
		//	if (!ct[icell].v_m[i].z)
		//		continue;
		//	current_cells[icell].ele -= ct[icell].v_m[i].b_ij * ct[icell].v_m[i].z *
		//		ct[icell].v_m[i].zc * dum;
		//	current_cells[icell].dif -= ct[icell].v_m[i].b_ij * ct[icell].v_m[i].z *
		//		ct[icell].v_m[i].grad;
		//}
		//dum1 = (current_x - current_cells[icell].dif) / current_cells[icell].ele;
		//if (isnan(dum1))
		//	dum1 = 1;
		//dum *= dum1;
		for (i = 0; i < ct[icell].J_ij_count_spec; i++)
		{
			if (!ct[icell].v_m[i].z)
				continue;
			ct[icell].J_ij[i].tot1 -= ct[icell].v_m[i].b_ij *
				ct[icell].v_m[i].zc * dum * DDt;
			ct[icell].J_ij[i].tot2 = ct[icell].J_ij[i].tot1;
			//ct[icell].J_ij_sum += ct[icell].v_m[i].z * ct[icell].J_ij[i].tot1;
		}
		current_A = current_x * F_C_MOL;
	}
	/*
	* calculate interlayer mass transfer...
	*/
	if (il_calcs && ct[icell].Dz2c_il != 0 && ct[icell].J_ij_il_count_spec > 0)
	{
		cxxExchange *ex_ptr1 = Utilities::Rxn_find(Rxn_exchange_map, icell);
		cxxExchange *ex_ptr2 = Utilities::Rxn_find(Rxn_exchange_map, jcell);
		Sum_zM = 0.0;
		i_max = k_il = ct[icell].J_ij_il_count_spec;
		for (i = 0; i < i_max; i++)
			Sum_zM += ct[icell].v_m_il[i].Dz * ct[icell].v_m_il[i].grad;
		for (i = 0; i < i_max; i++)
		{
			ct[icell].J_ij_il[i].tot1 = -ct[icell].v_m_il[i].D * ct[icell].v_m_il[i].grad +
				Sum_zM * ct[icell].v_m_il[i].Dzc / ct[icell].Dz2c_il;
			if (stagnant)
				ct[icell].J_ij_il[i].tot1 *= ct[icell].mixf_il;
			else
				ct[icell].J_ij_il[i].tot1 *= ct[icell].A_ij_il * DDt;
			//ct[icell].J_ij_sum += ct[icell].v_m_il[i].z * ct[icell].J_ij_il[i].tot1;
			ct[icell].J_ij_il[i].tot2 = ct[icell].J_ij_il[i].tot1;
			// ct[icell].J_ij_il[i].charge = ct[icell].v_m_il[i].z; // appt not used now
		}

		/* express the transfer in elemental moles... */
		tot1_h = tot1_o = tot2_h = tot2_o = 0.0;
		m_s = (struct M_S *) free_check_null(m_s);
		m_s = (struct M_S *) PHRQ_malloc((size_t) count_moles_added *
			sizeof(struct M_S));
		if (m_s == NULL)
			malloc_error();
		for (i1 = 0; i1 < count_moles_added; i1++)
		{
			m_s[i1].charge = 0;
			m_s[i1].name = NULL;
			m_s[i1].tot1 = 0;
			m_s[i1].tot2 = 0;
		}
		count_m_s = 0;
		fill_m_s(ct[icell].J_ij_il, k_il, icell, stagnant);

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
	/* do not transport charge imbalance */
	//ct[icell].J_ij_sum = 0;
	//V_M = (struct V_M *) free_check_null(V_M);
	if (il_calcs)
		ct[icell].v_m_il = (struct V_M *) free_check_null(ct[icell].v_m_il);
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
	int i, i1, i2, k, k1/*, n1, n2*/; 
	int charge_done, surf1, surf2;
	LDBLE f1, f2, mixf, mixf_store, mcd_mixf;
	LDBLE lav, A_ij, por, Dp1, Dp2;
	cxxMix * mix_ptr;
	cxxSurface *surface_ptr1, *surface_ptr2;
	LDBLE viscos_f;
	/*
	* temperature and viscosity correction for MCD coefficient, D_T = D_298 * Tk * viscos_298 / (298 * viscos)
	*/
	viscos_f = viscos_0;
	viscos_f = tk_x * viscos_0_25 / (298.15 * viscos_f);

	//n1 = 0;
	//n2 = n1 + 1;
	cxxSurface surface_n1, surface_n2;
	cxxSurface *surface_n2_ptr;

	std::map<int, cxxSurface> Rxn_temp_surface_map;

	for (i1 = 1; i1 <= count_cells + 1; i1++)
	{
		if (i1 <= count_cells && cell_data[i1].por < multi_Dpor_lim)
			continue;

		if (i1 == 1 && bcon_first != 1)
			continue;
		if (i1 == count_cells + 1 && bcon_last != 1)
			continue;

		i2 = i1 - 1;
		if (i2 > 0 && cell_data[i2].por < multi_Dpor_lim)
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
			for (size_t i3 = 0; i3 < num.size(); i3++)
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
				por = cell_data[1].por;
				lav = cell_data[1].length / 2;
				A_ij = Utilities::Rxn_find(Rxn_solution_map, 1)->Get_mass_water() / cell_data[1].length;
			}
			else if (i1 == count_cells + 1)
			{
				por = cell_data[count_cells].por;
				lav = cell_data[count_cells].length / 2;
				A_ij =
					Utilities::Rxn_find(Rxn_solution_map, count_cells)->Get_mass_water() /
					cell_data[count_cells].length;
			}
			else
			{
				por = cell_data[i2].por;
				lav =
					(cell_data[i1].length + cell_data[i2].length) / 2;
				A_ij =
					Utilities::Rxn_find(Rxn_solution_map, i1)->Get_mass_water() /
					(cell_data[i1].length * cell_data[i1].por);
				A_ij +=
					Utilities::Rxn_find(Rxn_solution_map, i2)->Get_mass_water() /
					(cell_data[i2].length * cell_data[i2].por);
				A_ij /= 2;
				A_ij *=
					(cell_data[i1].por <
					cell_data[i2].por ? cell_data[i1].por : cell_data[i2].por);
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
									pow(cell_data[i1].por, multi_Dn);
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
								pow(cell_data[i1].por, multi_Dn);
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
	//n2 = 0;
	std::map<int, cxxSurface>::iterator jit = Rxn_temp_surface_map.begin();
	for (; jit != Rxn_temp_surface_map.end(); jit++)
	{
		i = jit->first;
		assert(i == jit->second.Get_n_user());
		if ((i == 0 && bcon_first == 1) || (i == count_cells + 1 && bcon_last == 1))
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
sum_surface_comp(cxxSurface *source1, LDBLE f1, cxxSurface *source2,
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
		error_string = sformatf("Null pointer for surface 1 in sum_surface.");
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
	int i, i1, i2, k, k1, /*n1, n2,*/ ns;
	int charge_done, surf1, surf2;
	LDBLE f1, f2, mixf, mixf_store;
	LDBLE Dp1, Dp2;
	cxxMix *mix_ptr;
	cxxSurface *surface_ptr1, *surface_ptr2;
	LDBLE viscos_f;
	/*
	* temperature and viscosity correction for MCD coefficient, D_T = D_298 * Tk * viscos_298 / (298 * viscos)
	*/
	viscos_f = viscos_0;
	viscos_f = tk_x * viscos_0_25 / (298.15 * viscos_f);

	cxxSurface surface_n1, surface_n2;
	cxxSurface *surface_n1_ptr = &surface_n1;
	cxxSurface *surface_n2_ptr;
	std::map<int, cxxSurface> Rxn_temp_surface_map;

	for (ns = 0; ns < stag_data->count_stag; ns++)
	{

		i1 = mobile_cell + 1 + ns * count_cells;
		if (ns == 0)
			i1--;

		if (cell_data[i1].por < multi_Dpor_lim)
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
			if (cell_data[i2].por < multi_Dpor_lim)
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
				//n1 = 0;
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
				//n2 = 1;
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
					(cell_data[i1].por <= cell_data[i2].por ? cell_data[i1].por : cell_data[i2].por);
				mixf_store /= (default_Dw * pow(multi_Dpor, multi_Dn) * multi_Dpor);
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
						if (comp_k_ptr->Get_charge_name() == comp_k1_ptr->Get_charge_name())
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
							pow(cell_data[i2].por, multi_Dn) * viscos_f;
						Dp1 = 0;
						if (surf1)
						{
							for (k1 = 0; k1 < (int) surface_ptr1->Get_surface_comps().size(); k1++)
							{
								cxxSurfaceComp *comp_k1_ptr = &(surface_ptr1->Get_surface_comps()[k1]);
								if (strcmp(comp_k1_ptr->Get_formula().c_str(),
									comp_k_ptr->Get_formula().c_str()) != 0)
									continue;
								Dp1 = comp_k1_ptr->Get_Dw() * pow(cell_data[i1].por, multi_Dn) * viscos_f;
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
						if (comp_k_ptr->Get_charge_name() == comp_k1_ptr->Get_charge_name())
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
						Dp1 = comp_k_ptr->Get_Dw() * pow(cell_data[i1].por, multi_Dn) * viscos_f;

						Dp2 = 0;
						if (surf2)
						{
							for (k1 = 0; k1 < (int) surface_ptr2->Get_surface_comps().size(); k1++)
							{
								cxxSurfaceComp *comp_k1_ptr = &(surface_ptr2->Get_surface_comps()[k1]);
								if (strcmp(comp_k1_ptr->Get_formula().c_str(),
									comp_k_ptr->Get_formula().c_str()) != 0)
									continue;
								Dp2 = comp_k1_ptr->Get_Dw() * pow(cell_data[i2].por, multi_Dn) * viscos_f;
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
	//n2 = 0;
	std::map<int, cxxSurface>::iterator jit = Rxn_temp_surface_map.begin();
	for (; jit != Rxn_temp_surface_map.end(); jit++)
	{
		i = jit->first;
		assert(i == jit->second.Get_n_user());
		if ((i == 0 && bcon_first == 1) || (i == count_cells + 1 && bcon_last == 1))
		{
			continue;
		}
		if (i >= 0 && i <= 1 + count_cells * (1 + stag_data->count_stag))
		{
			surface_ptr1 = Utilities::Rxn_find(Rxn_surface_map, i);
			//if (surface_ptr1 != NULL)
			//{
			//	Rxn_surface_map[i] = jit->second;
			//}
			//else
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

	/* from Atkins, 1994. Physical Chemistry, 5th ed. */
	//viscos =
	//	pow((LDBLE) 10.,
	//		-(1.37023 * (tc_x - 20) +
	//		  0.000836 * (tc_x - 20) * (tc_x - 20)) / (109 + tc_x));
	/* Huber et al., 2009, J. Phys. Chem. Ref. Data, Vol. 38, 101-125 */
	LDBLE H[4] = { 1.67752, 2.20462, 0.6366564, -0.241605 };
	LDBLE Tb = tk_x / 647.096, denom = H[0], mu0;
	int i, j, i1;
	for (i = 1; i < 4; i++)
		denom += H[i] / pow(Tb, i);
	mu0 = 100.0 * sqrt(Tb) / denom;

	LDBLE H2[42] =
	{ 5.20094e-1, 2.22531e-1, -2.81378e-1, 1.61913e-1, -3.25372e-2, 0, 0,
	8.50895e-2, 9.99115e-1, -9.06851e-1, 2.57399e-1, 0, 0, 0,
	-1.08374, 1.88797, -7.72479e-1, 0, 0, 0, 0,
	-2.89555e-1, 1.26613, -4.89837e-1, 0, 6.98452e-2, 0, -4.35673e-3,
	0, 0, -2.57040e-1, 0, 0, 8.72102e-3, 0,
	0, 1.20573e-1, 0, 0, 0, 0, -5.93264e-4 };
	LDBLE Rb = rho_0 / 0.322, Tb_1, Rb_1, S1 = 0.0, S2, mu1;
	Tb_1 = 1.0 / Tb - 1.0; Rb_1 = Rb - 1.0;
	for (i = 0; i < 6; i++)
	{
		S2 = 0.0;
		for (j = 0; j < 7; j++)
		{
			i1 = 7 * i;
			if (!H2[i1 + j])
				continue;
			if (j)
				S2 += H2[i1 + j] * pow(Rb_1, j);
			else
				S2 += H2[i1];
		}
		if (i)
			S1 += S2 * pow(Tb_1, i);
		else
			S1 += S2;
	}
	mu1 = exp(Rb * S1);
	viscos_0 = viscos = mu0 * mu1 / 1e3;
	viscos_0_25 = 0.8900239182946;
	//#define OLD_VISCOSITY
#ifdef OLD_VISCOSITY
	/* from Atkins, 1994. Physical Chemistry, 5th ed. */
	viscos =
		pow((LDBLE) 10.,
		-(1.37023 * (tc_x - 20) +
		0.000836 * (tc_x - 20) * (tc_x - 20)) / (109 + tc_x));
	viscos_0_25 = 0.88862;
#endif
	//return viscos;
	if (!print_viscosity)
		return viscos;

	/* (modified) Jones-Dole eqn for viscosity:
	viscos / viscos_0 =
	1 + A * eq_tot^0.5 +
	f_an * (Sum(B_i * m_i) +
	Sum(D_i * m_i * ((1 + f_I) * mu_x^d3_i + (m_i * f_z)^d3_i) / (2 + f_I)))
	A calculated from Falkenhagen-Dole
	B_i = b0 + b1*exp(b2 * tc), b0..2 in Jones_Dole[0..2], read in SOLUTION_SPECIES
	D_i = d1 * exp(d2 * tc), d1, 2 in Jones_Dole[3, 4]
	d3_i in Jones_Dole[5]
	Jones_Dole[6] contains the anion factor, 1 for Cl-, variable for other anions
	f_z = (z * z + |z|) / 2, the contribution of the ion to mu_x, if z = 0: f_z = mu_x / m_i
	f_I = variable, depends on d3_i > 1, or d3_i < 1.
	tc is limited to 200 C.


	A from Falkenhagen-Dole for a salt:
	A = 4.3787e-14 * TK**1.5 / (eps_r**0.5)* (z1 + z2)**-0.5 / (D1 * D2) * psi
	psi = (D1*z2 + D2*z1)/4 - z1*z2 * (D1-D2)**2 / ((D1*z1 + D2*z2)**0.5 + ((D1 + D2) * (z1 + z2))**0.5)**2
	D1, z1 for the cation, D2, |z2| for the anion of the salt.
	We use the harmonic mean of the Dw's, and the arithmetic mean of the z's,
	both weighted by the equivalent concentration.
	*/
	LDBLE D1, D2, z1, z2, m_plus, m_min, eq_plus, eq_min, eq_dw_plus, eq_dw_min, t1, t2, ta;
	LDBLE A, psi, Bc = 0, Dc = 0, Dw = 0.0, l_z, f_z, lm, V_an, m_an, V_Cl, tc;

	m_plus = m_min = eq_plus = eq_min = eq_dw_plus = eq_dw_min = V_an = m_an = V_Cl = ta = 0;

	tc = (tc_x > 200) ? 200 : tc_x;

	for (i = 0; i < count_s_x; i++)
	{
		if (s_x[i]->type != AQ && s_x[i]->type > HPLUS)
			continue;
		if ((lm = s_x[i]->lm) < -9)
			continue;
		if (s_x[i]->Jones_Dole[0] || s_x[i]->Jones_Dole[1] || s_x[i]->Jones_Dole[3])
		{
			t1 = s_x[i]->moles / mass_water_aq_x;
			l_z = fabs(s_x[i]->z);
			if (l_z)
				f_z = (l_z * l_z + l_z) / 2;
			else
				f_z = mu_x / t1;
			//if data at tc's other than 25 are scarce, put the values found for 25 C in [7] and [8], optimize [1], [2], and [4]...
			if (s_x[i]->Jones_Dole[7] || s_x[i]->Jones_Dole[8])
			{
				s_x[i]->Jones_Dole[0] = s_x[i]->Jones_Dole[7] -
					s_x[i]->Jones_Dole[1] * exp(-s_x[i]->Jones_Dole[2] * 25.0);
				s_x[i]->Jones_Dole[3] =
					s_x[i]->Jones_Dole[8] / exp(-s_x[i]->Jones_Dole[4] * 25.0);
			}
			// find B * m and D * m * mu^d3
			Bc += (s_x[i]->Jones_Dole[0] +
				s_x[i]->Jones_Dole[1] * exp(-s_x[i]->Jones_Dole[2] * tc)) *
				t1;
			// define f_I from the exponent of the D * m^d3 term...
			if (s_x[i]->Jones_Dole[5] >= 1)
				t2 = mu_x / 3 / s_x[i]->Jones_Dole[5];
			else if (s_x[i]->Jones_Dole[5] > 0.4)
				t2 = -0.8 / s_x[i]->Jones_Dole[5];
			else
				t2 = -1;
			Dc += (s_x[i]->Jones_Dole[3] * exp(-s_x[i]->Jones_Dole[4] * tc)) *
				t1 * (pow(mu_x, s_x[i]->Jones_Dole[5])*(1 + t2) + pow(t1 * f_z, s_x[i]->Jones_Dole[5])) / (2 + t2);
			//output_msg(sformatf("\t%s\t%e\t%e\t%e\n", s_x[i]->name, t1, Bc, Dc ));
		}
		// parms for A...
		if ((l_z = s_x[i]->z) == 0)
			continue;
		Dw = s_x[i]->dw;
		if (Dw)
		{
			Dw *= (0.89 / viscos_0 * tk_x / 298.15);
			if (s_x[i]->dw_t)
				Dw *= exp(s_x[i]->dw_t / tk_x - s_x[i]->dw_t / 298.15);
		}
		if (l_z < 0)
		{
			if (!strcmp(s_x[i]->name, "Cl-"))
				// volumina for f_an...
			{
				V_Cl = s_x[i]->logk[vm_tc];
				V_an += V_Cl * s_x[i]->moles;
				ta += s_x[i]->moles;
				m_an += s_x[i]->moles;
			}
			else// if (s_x[i]->Jones_Dole[6])
			{
				V_an += s_x[i]->logk[vm_tc] * s_x[i]->Jones_Dole[6] * s_x[i]->moles;
				ta += s_x[i]->moles;
				m_an += s_x[i]->moles;
			}
			if (Dw)
			{
				// anions for A...
				m_min += s_x[i]->moles;
				t1 = s_x[i]->moles * l_z;
				eq_min -= t1;
				eq_dw_min -= t1 / Dw;
			}
		}
		else if (Dw)
		{
			// cations for A...
			m_plus += s_x[i]->moles;
			t1 = s_x[i]->moles * l_z;
			eq_plus += t1;
			eq_dw_plus += t1 / Dw;
		}
	}
	if (m_plus && m_min && eq_dw_plus && eq_dw_min)
	{
		z1 = eq_plus / m_plus;     z2 = eq_min / m_min;
		D1 = eq_plus / eq_dw_plus; D2 = eq_min / eq_dw_min;

		t1 = (D1 - D2) / (sqrt(D1 * z1 + D2 * z2) + sqrt((D1 + D2) * (z1 + z2)));
		psi = (D1 * z2 + D2 * z1) / 4.0 - z1 * z2 * t1 * t1;
		// Here A is A * viscos_0, avoids multiplication later on...
		A = 4.3787e-14 * pow(tk_x, 1.5) / (sqrt(eps_r * (z1 + z2) / ((z1 > z2) ? z1 : z2)) * (D1 * D2)) * psi;
	}
	else
		A = 0;
	viscos = viscos_0 + A * sqrt((eq_plus + eq_min) / 2 / mass_water_aq_x);
	if (m_an)
	{
		V_an /= m_an;
		ta /= m_an;
	}
	if (!V_Cl)
		V_Cl = calc_vm_Cl();
	if (V_an && V_Cl && ta)
		viscos += (viscos_0 * (2 - ta * V_an / V_Cl) * (Bc + Dc));
	else
		viscos += (viscos_0 * (Bc + Dc));
	if (viscos < 0)
	{
		viscos = 0;
		warning_msg("viscosity < 0, reset to 0.");
	}
	return viscos;

}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
calc_vm_Cl(void)
/* ---------------------------------------------------------------------- */
{
	/*
	*  Calculate molar volume of Cl- with a Redlich type eqn:
	Vm = Vm0(tc) + (Av / 2) * z^2 * I^0.5 + coef(tc) * I^(b4).
	*    Vm0(tc) is calc'd using supcrt parms, or from millero[0] + millero[1] * tc + millero[2] * tc^2
	*    for Av * z^2 * I^0.5, see Redlich and Meyer, Chem. Rev. 64, 221.
	Av is in (cm3/mol)(mol/kg)^-0.5, = DH_Av.
	If b_Av != 0, the extended DH formula is used: I^0.5 /(1 + b_Av * DH_B * I^0.5).
	DH_Av and DH_B are from calc_dielectrics(tc, pa).
	*	  coef(tc) = logk[vmi1] + logk[vmi2] / (TK - 228) + logk[vmi3] * (TK - 228).
	*    b4 = logk[vmi4], or
	*	  coef(tc) = millero[3] + millero[4] * tc + millero[5] * tc^2
	*/
	LDBLE V_Cl = 0;
	LDBLE pb_s = 2600. + patm_x * 1.01325, TK_s = tc_x + 45.15, sqrt_mu = sqrt(mu_x);
	struct species *s_ptr;

	s_ptr = s_search("Cl-");
	if (!s_ptr)
		return V_Cl;

	if (s_ptr->logk[vma1])
	{
		/* supcrt volume at I = 0... */
		V_Cl = s_ptr->logk[vma1] + s_ptr->logk[vma2] / pb_s +
			(s_ptr->logk[vma3] + s_ptr->logk[vma4] / pb_s) / TK_s -
			s_ptr->logk[wref] * QBrn;
		/* the ionic strength term * I^0.5... */
		if (s_ptr->logk[b_Av] < 1e-5)
			V_Cl += s_ptr->z * s_ptr->z * 0.5 * DH_Av * sqrt_mu;
		else
		{
			/* limit the Debye-Hueckel slope by b... */
			/* pitzer... */
			//s_ptr->rxn_x->logk[vm_tc] += s_ptr->z * s_ptr->z * 0.5 * DH_Av *
			//	log(1 + s_ptr->logk[b_Av] * sqrt(mu_x)) / s_ptr->logk[b_Av];
			/* extended DH... */
			V_Cl += s_ptr->z * s_ptr->z * 0.5 * DH_Av *
				sqrt_mu / (1 + s_ptr->logk[b_Av] * DH_B * sqrt_mu);
		}
		/* plus the volume terms * I... */
		if (s_ptr->logk[vmi1] != 0.0 || s_ptr->logk[vmi2] != 0.0 || s_ptr->logk[vmi3] != 0.0)
		{
			LDBLE bi = s_ptr->logk[vmi1] + s_ptr->logk[vmi2] / TK_s + s_ptr->logk[vmi3] * TK_s;
			if (s_ptr->logk[vmi4] == 1.0)
				V_Cl += bi * mu_x;
			else
				V_Cl += bi * pow(mu_x, s_ptr->logk[vmi4]);
		}
	}
	else if (s_ptr->millero[0])
	{
		/* Millero volume at I = 0... */
		V_Cl = s_ptr->millero[0] + tc_x * (s_ptr->millero[1] + tc_x * s_ptr->millero[2]);
		if (s_ptr->z)
		{
			/* the ionic strength terms... */
			V_Cl += s_ptr->z * s_ptr->z * 0.5 * DH_Av * sqrt_mu +
				(s_ptr->millero[3] + tc_x * (s_ptr->millero[4] + tc_x * s_ptr->millero[5])) * mu_x;
		}
	}
	return V_Cl;
}
