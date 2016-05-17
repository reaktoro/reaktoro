#include "Utils.h"
#include "Phreeqc.h"
#include "phqalloc.h"
#include <stdarg.h>
#include <assert.h>
#include "Temperature.h"
#include "cxxMix.h"
#include "Exchange.h"
#include "GasPhase.h"
#include "Reaction.h"
#include "PPassemblage.h"
#include "SSassemblage.h"
#include "cxxKinetics.h"
#include "Solution.h"
/* ---------------------------------------------------------------------- */
int Phreeqc::
array_print(LDBLE * array_l, int row_count, int column_count,
			int l_max_column_count)
/* ---------------------------------------------------------------------- */
{
	int i, j, k;

	for (i = 0; i < row_count; i++)
	{
		k = 0;
		output_msg(sformatf("%d\n", i));
		for (j = 0; j < column_count; j++)
		{
			if (k > 7)
			{
				output_msg(sformatf("\n"));
				k = 0;
			}
			output_msg(sformatf("%11.2e",
					   (double) array_l[i * l_max_column_count + j]));
			k++;
		}
		if (k != 0)
		{
			output_msg(sformatf("\n"));
		}
		output_msg(sformatf("\n"));
	}
	output_msg(sformatf("\n"));
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
set_pr_in_false(void)
/* ---------------------------------------------------------------------- */
{
	// set pr_in to false for subsequent steps...
	if (use.Get_pp_assemblage_in())
	{
		for (int i = 0; i < count_unknowns; i++)
		{
			if (x[i]->type == PP)
				x[i]->phase->pr_in = false;
		}
	}
	if (use.Get_gas_phase_ptr())
	{
		cxxGasPhase *gas_phase_ptr = use.Get_gas_phase_ptr();
		for (size_t i = 0; i < gas_phase_ptr->Get_gas_comps().size(); i++)
		{
			cxxGasComp *gc_ptr = &(gas_phase_ptr->Get_gas_comps()[i]);
			int k;
			struct phase *phase_ptr = phase_bsearch(gc_ptr->Get_phase_name().c_str(), &k, FALSE);
			if (phase_ptr)
				phase_ptr->pr_in = false;
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
print_all(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print each block of results in sequence to file output
 *   also print selected output to punch_file
 *   Each routine is controlled by a variable in structure print.
 *   print.all == FALSE will turn off all prints.
 */
	if (pr.all == FALSE)
	{
		set_pr_in_false();
		return (OK);
	}
/*
 *   Makes sorted list including all species for each valence state
 */
	if (pr.surface == TRUE || pr.exchange == TRUE || pr.species == TRUE)
	{
		species_list_sort();
	}
/*
 *   Print results
 */
	s_h2o->lm = s_h2o->la;
	print_using();
	print_mix();
	print_reaction();
	print_kinetics();
	print_user_print();
	print_gas_phase();
	print_pp_assemblage();
	print_ss_assemblage();
	print_surface();
	print_exchange();
	print_initial_solution_isotopes();
	print_isotope_ratios();
	print_isotope_alphas();
	print_totals();
	print_eh();
	print_species();
	print_alkalinity();
	print_saturation_indices();
	if (!pr.saturation_indices)
		set_pr_in_false();
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
punch_all(void)
/* ---------------------------------------------------------------------- */
{
//#ifndef PHREEQ98		/* if not PHREEQ98 use the standard declaration */
//	if (pr.hdf == FALSE && (punch.in == FALSE || pr.punch == FALSE) && user_graph->commands == NULL)
//		return (OK);
///*
// *   Punch results
// */
//	if (state == TRANSPORT || state == PHAST || state == ADVECTION)
//	{
//		use.Get_kinetics_ptr() = kinetics_bsearch(use.Get_n_kinetics_user(), &i);
//	}
//	else if (use.Get_kinetics_in() != FALSE)
//	{
//		use.Get_kinetics_ptr() = kinetics_bsearch(-2, &i);
//	}
//#else /* if PHREEQ98 execute punch_user_graph first, that is, before punch.in and pr.punch is checked */
	if (state == TRANSPORT || state == PHAST || state == ADVECTION)
	{
		use.Set_kinetics_ptr(Utilities::Rxn_find(Rxn_kinetics_map, use.Get_n_kinetics_user()));
	}
	else if (use.Get_kinetics_in() != FALSE)
	{
		use.Set_kinetics_ptr(Utilities::Rxn_find(Rxn_kinetics_map, -2));
	}
#if defined PHREEQ98
	if (pr.user_graph == TRUE)
	{
		punch_user_graph();
	}
#elif  defined MULTICHART
	if (pr.user_graph == TRUE)
	{
			chart_handler.Punch_user_graph(this);
	}
#endif
	//if (pr.hdf == FALSE && (punch.in == FALSE || pr.punch == FALSE))
	if (pr.hdf == FALSE && (SelectedOutput_map.size() == 0 || pr.punch == FALSE))
		return (OK);
	std::map < int, SelectedOutput >::iterator so_it = SelectedOutput_map.begin();
	for ( ; so_it != SelectedOutput_map.end(); so_it++)
	{
		current_selected_output = &(so_it->second);
		if (pr.punch == FALSE ||
			current_selected_output == NULL ||
			!current_selected_output->Get_active() /* ||
			current_selected_output->Get_punch_ostream() == NULL*/)
			continue;
		phrq_io->Set_punch_ostream(current_selected_output->Get_punch_ostream());

		// UserPunch
		std::map < int, UserPunch >::iterator up_it = UserPunch_map.find(current_selected_output->Get_n_user());
		current_user_punch = up_it == UserPunch_map.end() ? NULL : &(up_it->second);

		punch_identifiers();
		punch_totals();
		punch_molalities();
		punch_activities();
		punch_pp_assemblage();
		punch_saturation_indices();
		punch_gas_phase();
		punch_kinetics();
		punch_ss_assemblage();
		punch_isotopes();
		punch_calculate_values();
		punch_user_punch();
		/*
		*   new line for punch_file
		*/
		punch_msg("\n");

		/*
		*   signal end of row
		*/
		fpunchf_end_row("\n");
		punch_flush();
	}
	current_selected_output = NULL;
	current_user_punch = NULL;
	phrq_io->Set_punch_ostream(NULL);
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
print_diffuse_layer(cxxSurfaceCharge *charge_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints total moles of each element in diffuse layer
 *   Remove comment to print total moles of each species
 */
	LDBLE mass_water_surface, r, sum_surfs;
	LDBLE molality, moles_excess, moles_surface, d;

	if (use.Get_surface_ptr() == NULL)
		return (OK);
/*
 *   Find position of component in surface charge data
 */
	int j;
	for (j = 0; j < count_unknowns; j++)
	{
		if (x[j]->type != SURFACE_CB)
			continue;
		cxxSurfaceCharge * charge_ptr_search = use.Get_surface_ptr()->Find_charge(x[j]->surface_charge);
		if (charge_ptr->Get_name() == charge_ptr_search->Get_name())
		{
			break;
		}
	}
	if (j >= count_unknowns)
	{
		error_string = sformatf(
				"In print_diffuse_layer: component not found, %s.",
				charge_ptr->Get_name().c_str());
		error_msg(error_string, STOP);
	}
/*
 *   Loop through all surface components, calculate each H2O surface (diffuse layer),
 *   H2O aq, and H2O bulk (diffuse layers plus aqueous).
 */

	if (mass_water_surfaces_x != 0)
	{
		d = 100 * charge_ptr->Get_mass_water() / mass_water_surfaces_x;
	}
	else
	{
		d = 0.0;
	}
	output_msg(sformatf(
			   "\tWater in diffuse layer: %8.3e kg, %4.1f%% of total DDL-water.\n",
			   (double) charge_ptr->Get_mass_water(), (double) d));
	if (use.Get_surface_ptr()->Get_debye_lengths() > 0)
	{
		sum_surfs = 0.0;
		for (j = 0; j < count_unknowns; j++)
		{
			if (x[j]->type != SURFACE_CB)
				continue;
			cxxSurfaceCharge * charge_ptr_search = use.Get_surface_ptr()->Find_charge(x[j]->surface_charge);
			sum_surfs +=
				charge_ptr_search->Get_specific_area() *
				charge_ptr_search->Get_grams();
		}
		r = 0.002 * mass_water_bulk_x / sum_surfs;
		output_msg(sformatf(
				   "\tRadius of total pore:   %8.3e m; of free pore: %8.3e m.\n",
				   (double) r, (double) (r - use.Get_surface_ptr()->Get_thickness())));
	}

	if (debug_diffuse_layer == TRUE)
	{
		output_msg(sformatf(
				   "\n\t\tDistribution of species in diffuse layer\n\n"));
		output_msg(sformatf(
				   "\n\tSpecies     \t    Moles   \tMoles excess\t      g\n"));
	}
	mass_water_surface = charge_ptr->Get_mass_water();
	count_elts = 0;
	paren_count = 0;
	for (j = 0; j < count_s_x; j++)
	{
		if (s_x[j]->type > HPLUS)
			continue;
		molality = under(s_x[j]->lm);
		moles_excess = mass_water_aq_x * molality * (charge_ptr->Get_g_map()[s_x[j]->z].Get_g() *
			s_x[j]->erm_ddl +
			mass_water_surface /
			mass_water_aq_x * (s_x[j]->erm_ddl - 1));
		moles_surface = mass_water_surface * molality + moles_excess;
		if (debug_diffuse_layer == TRUE)
		{
			output_msg(sformatf("\t%-12s\t%12.3e\t%12.3e\t%12.3e\n",
					   s_x[j]->name, moles_surface, moles_excess,
					   charge_ptr->Get_g_map()[s_x[j]->z].Get_g()));
		}
/*
 *   Accumulate elements in diffuse layer
 */
		add_elt_list(s_x[j]->next_elt, moles_surface);
	}
/*
	strcpy(token, s_h2o->name);
	ptr = &(token[0]);
	get_elts_in_species (&ptr, mass_water_surface / gfw_water);
 */
	if (count_elts > 0)
	{
		qsort(elt_list, (size_t) count_elts,
			  (size_t) sizeof(struct elt_list), elt_list_compare);
		elt_list_combine();
	}
/*
 *   Print totals
 */
	if (use.Get_surface_ptr()->Get_dl_type() != cxxSurface::DONNAN_DL)
	{
		output_msg(sformatf(
				   "\n\tTotal moles in diffuse layer (excluding water)\n\n"));
	}
	else
	{
		LDBLE exp_g =  charge_ptr->Get_g_map()[1].Get_g() * mass_water_aq_x / mass_water_surface + 1;
		LDBLE psi_DL = -log(exp_g) * R_KJ_DEG_MOL * tk_x / F_KJ_V_EQ;
		output_msg(sformatf(
				   "\n\tTotal moles in diffuse layer (excluding water), Donnan calculation."));
		output_msg(sformatf(
				   "\n\tDonnan Layer potential, psi_DL = %10.3e V.\n\tBoltzmann factor, exp(-psi_DL * F / RT) = %9.3e (= c_DL / c_free if z is +1).\n\n",
					 psi_DL, exp_g));
	}
	output_msg(sformatf("\tElement       \t     Moles\n"));
	for (j = 0; j < count_elts; j++)
	{
		output_msg(sformatf("\t%-14s\t%12.4e\n",
				   elt_list[j].elt->name, (double) elt_list[j].coef));
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
print_eh(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints eh calculated from redox couples
 *   Only calculates eh if two redox states of an element have mass balance
 *   equations.
 */
	int i, j, k, first;
	LDBLE pe, eh;
	struct master *master_ptr0, *master_ptr1;
	char token[MAX_LENGTH];

	if (pr.eh == FALSE || pr.all == FALSE)
		return (OK);

	tk_x = tc_x + 273.15;

	first = TRUE;
	for (i = 0; i < count_master; i++)
	{
		if (master[i]->in != TRUE)
			continue;
		if (master[i]->primary == TRUE)
			continue;
/*
 *   Secondary master species has mass balance equation
 */
		master_ptr0 = master[i]->elt->primary;
		for (k = i + 1; k < count_master; k++)
		{
			if (master[k]->in != TRUE)
				continue;
			master_ptr1 = master[k]->elt->primary;
			if (master_ptr1 != master_ptr0)
				break;
/*
 *  Another secondary master species of same element has mass balance equation
 *  Rewrite equations to calculate pe
 */
			rewrite_master_to_secondary(master[k], master[i]);
			trxn_swap("e-");
/* debug
			trxn_print();
 */
/*
 *   Calculate pe, eh
 */
			pe = -k_calc(trxn.logk, tk_x, patm_x * PASCAL_PER_ATM);
			for (j = 1; j < count_trxn; j++)
			{
				pe -= trxn.token[j].s->la * trxn.token[j].coef;
			}
			eh = ((LOG_10 * R_KJ_DEG_MOL * tk_x) / F_KJ_V_EQ) * pe;
/*
 *   Print heading
 */
			if (first == TRUE)
			{
				print_centered("Redox couples");
				output_msg(sformatf("\t%-15s%12s%12s\n\n",
						   "Redox couple", "pe", "Eh (volts)"));
				first = FALSE;
			}
/*
 *   Print result
 */
			strcpy(token, master[i]->elt->name);
			strcat(token, "/");
			strcat(token, master[k]->elt->name);
			output_msg(sformatf("\t%-15s%12.4f%12.4f\n", token,
					   (double) pe, (double) eh));
		}
	}
	if (first == FALSE)
		output_msg(sformatf("\n"));
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
print_exchange(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print moles of each exchange species
 */
	int i;
	cxxExchange * exchange_ptr;
	const char *name, *name1;
	struct master *master_ptr;
	LDBLE dum, dum2;
/*
 *  Print exchange data
 */
	exchange_ptr = use.Get_exchange_ptr();
	if (exchange_ptr == NULL || pr.exchange == FALSE || pr.all == FALSE)
		return (OK);

	if (state >= REACTION)
	{
		print_centered("Exchange composition");
	}
/*
 *   Print list of species
 */

	s_h2o->lm = s_h2o->la;
	name = s_hplus->secondary->elt->name;
	for (i = 0; i < count_species_list; i++)
	{
/*
 *   Get name of master species
 */
		if (species_list[i].s->type != EX)
			continue;
		if (species_list[i].master_s->secondary != NULL)
		{
			master_ptr = species_list[i].master_s->secondary;
			name1 = species_list[i].master_s->secondary->elt->name;
		}
		else
		{
			master_ptr = species_list[i].master_s->primary;
			name1 = species_list[i].master_s->primary->elt->name;
		}
/*
 *   Check if new master species, print total molality
 */
		if (name1 != name)
		{
			name = name1;
			output_msg(sformatf("%-14s%12.3e mol", name,
					   (double) master_ptr->unknown->moles));
			cxxExchange *exchange_ptr = (cxxExchange *) (use.Get_exchange_ptr());
			const cxxExchComp *exchange_comp_ptr = exchange_ptr->Find_comp(master_ptr->unknown->exch_comp);
			assert(exchange_comp_ptr);
			if (exchange_comp_ptr->Get_phase_name().size() > 0)
			{
				output_msg(sformatf("\t[%g (mol %s)/(mol %s)]",
					(double) exchange_comp_ptr->Get_phase_proportion(),
					exchange_comp_ptr->Get_formula().c_str(),
					exchange_comp_ptr->Get_phase_name().c_str()));
			}
			else if (exchange_comp_ptr->Get_rate_name().size() > 0)
			{
				output_msg(sformatf(
						   "\t[%g (mol %s)/(mol kinetic reactant %s)]",
						   (double) exchange_comp_ptr->Get_phase_proportion(),
						   exchange_comp_ptr->Get_formula().c_str(),
						   exchange_comp_ptr->Get_rate_name().c_str()));
			}
			output_msg(sformatf("\n\n"));
			/* Heading for species */
			output_msg(sformatf("\t%-15s%12s%12s%12s%10s\n", " ", " ",
					   "Equiv-  ", "Equivalent", "Log "));
			output_msg(sformatf("\t%-15s%12s%12s%12s%10s\n\n",
					   "Species", "Moles  ", "alents  ", "Fraction", "Gamma"));
		}
/*
 *   Print species data
 */
/* !!!!! */
		if (master_ptr->total > 1.0e-10)
		{
			if (species_list[i].s->equiv != 0.0)
			{
				dum = fabs(species_list[i].s->equiv) / master_ptr->total;
			}
			else
			{
				if (species_list[i].master_s->z == 0)
				{
					dum = 1 / master_ptr->total;
				}
				else
				{
					dum = 1;
				}
			}
			if (species_list[i].master_s->z != 0.0)
			{
				dum2 = fabs(species_list[i].master_s->z);
			}
			else
			{
				dum2 = 1;
			}
			output_msg(sformatf("\t%-15s%12.3e%12.3e%12.3e%10.3f\n",
					   species_list[i].s->name,
					   (double) species_list[i].s->moles,
					   (double) (species_list[i].s->moles * dum2 *
								 species_list[i].s->equiv),
					   (double) (species_list[i].s->moles *
								 dum /* / dum2 */ ),
					   (double) (species_list[i].s->lg - log10(dum))));
		}
	}
	output_msg(sformatf("\n"));
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
print_gas_phase(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints gas phase composition if present
 */
	LDBLE lp, moles, initial_moles, delta_moles;
	struct rxn_token *rxn_ptr;
	char info[MAX_LENGTH];
	bool PR = false;

	if (pr.gas_phase == FALSE || pr.all == FALSE)
		return (OK);
	if (use.Get_gas_phase_ptr() == NULL)
		return (OK);

	cxxGasPhase *gas_phase_ptr = use.Get_gas_phase_ptr();
	if (gas_phase_ptr->Get_v_m() >= 0.01)
		PR = true;
	if (gas_phase_ptr->Get_type() == cxxGasPhase::GP_PRESSURE)
	{
		if (gas_unknown == NULL)
			return (OK);
		if (gas_unknown->moles < 1e-12)
		{
			sprintf(info, "Fixed-pressure gas phase %d dissolved completely",
				   use.Get_n_gas_phase_user());
			print_centered(info);
			return (OK);
		}
		gas_phase_ptr->Set_total_moles(gas_unknown->moles);
		gas_phase_ptr->Set_volume(gas_phase_ptr->Get_total_moles() * R_LITER_ATM * tk_x /
			gas_phase_ptr->Get_total_p());
		if (PR)
			gas_phase_ptr->Set_volume(gas_phase_ptr->Get_v_m() * gas_unknown->moles);
	}
/*
 *   Print heading
 */
	print_centered("Gas phase");
	output_msg(sformatf("Total pressure: %5.2f      atmospheres",
			   (double) gas_phase_ptr->Get_total_p()));
	if (gas_phase_ptr->Get_total_p() >= 1500)
		output_msg(" WARNING: Program limit.\n");
	else if (PR)
		output_msg("          (Peng-Robinson calculation)\n");
	else
		output_msg(" \n");
	output_msg(sformatf("    Gas volume: %10.2e liters\n",
			   (double) gas_phase_ptr->Get_volume()));
	if(gas_phase_ptr->Get_total_moles() > 0)
	{
		if (PR)
		{
			output_msg(sformatf("  Molar volume: %10.2e liters/mole",
				    (double) (gas_phase_ptr->Get_v_m())));
		}
		else
		{
			output_msg(sformatf("  Molar volume: %10.2e liters/mole",
				(double) (gas_phase_ptr->Get_volume() / gas_phase_ptr->Get_total_moles())));
		}
	}
	if (/*!numerical_fixed_volume && */((PR && gas_phase_ptr->Get_v_m() <= 0.016)))
		output_msg(" WARNING: Program limit for Peng-Robinson.\n");
	else
		output_msg("\n");
	if (PR)
	output_msg(sformatf( "   P * Vm / RT: %8.5f  (Compressibility Factor Z) \n",
			     (double) (gas_phase_ptr->Get_total_p() * gas_phase_ptr->Get_v_m() / (R_LITER_ATM * tk_x))));


	output_msg(sformatf("\n%68s\n%78s\n", "Moles in gas",
			   "----------------------------------"));
	if (PR)
		output_msg(sformatf( "%-11s%12s%12s%7s%12s%12s%12s\n\n", "Component",
			   "log P", "P", "phi", "Initial", "Final", "Delta"));
	else
		output_msg(sformatf("%-18s%12s%12s%12s%12s%12s\n\n", "Component",
			   "log P", "P", "Initial", "Final", "Delta"));

	for (size_t j = 0; j < gas_phase_ptr->Get_gas_comps().size(); j++)
	{
/*
 *   Calculate partial pressure
 */
		cxxGasComp *gc_ptr = &(gas_phase_ptr->Get_gas_comps()[j]);
		int k;
		struct phase *phase_ptr = phase_bsearch(gc_ptr->Get_phase_name().c_str(), &k, FALSE);
		if (phase_ptr->in == TRUE)
		{
			lp = -phase_ptr->lk;
			for (rxn_ptr =
				 phase_ptr->rxn_x->token + 1;
				 rxn_ptr->s != NULL; rxn_ptr++)
			{
				lp += rxn_ptr->s->la * rxn_ptr->coef;
			}
			lp -= phase_ptr->pr_si_f;
			moles = phase_ptr->moles_x;
		}
		else
		{
			lp = -99.99;
			moles = 0;
			phase_ptr->p_soln_x = 0;
		}
/*
 *   Print gas composition
 */
		if (state != TRANSPORT && state != PHAST)
		{
			initial_moles = gc_ptr->Get_moles();
			delta_moles = phase_ptr->moles_x - gc_ptr->Get_moles();
		}
		else
		{
			initial_moles = gc_ptr->Get_initial_moles();
			delta_moles = gc_ptr->Get_initial_moles() -
				gc_ptr->Get_moles();
		}
		if (moles <= MIN_TOTAL)
			moles = 0.0;
		if (fabs(delta_moles) <= MIN_TOTAL)
			delta_moles = 0.0;
		if (PR)
		{
			output_msg(sformatf("%-11s%12.2f%12.3e%7.3f%12.3e%12.3e%12.3e\n",
				   phase_ptr->name,
				   (double) lp,
				   (double) phase_ptr->p_soln_x,
				   (double) phase_ptr->pr_phi,
				   (double) initial_moles,
				   (double) moles,
				   (double) delta_moles));
		}
		else
			output_msg(sformatf("%-18s%12.2f%12.3e%12.3e%12.3e%12.3e\n",
				   phase_ptr->name,
				   (double) lp,
				   (double) phase_ptr->p_soln_x,
				   (double) initial_moles,
				   (double) moles,
				   (double) delta_moles));
		if (!strcmp(phase_ptr->name, "H2O(g)") && phase_ptr->p_soln_x == 90)
			output_msg("       WARNING: The pressure of H2O(g) is above the program limit: use the polynomial for log_k.\n");

	}
	output_msg("\n");
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
print_ss_assemblage(void)
/* ---------------------------------------------------------------------- */
{
	/*
	 *   Prints solid solution composition if present
	 */
	int i, j;
	LDBLE delta_moles;
	LDBLE nb, nc, xb, xb1, xb2, xb1moles, xb2moles;

	if (pr.ss_assemblage == FALSE || pr.all == FALSE)
		return (OK);
	if (use.Get_ss_assemblage_ptr() == NULL)
		return (OK);
	/*
	 *   Print heading
	 */
	print_centered("Solid solutions");
	output_msg(sformatf("\n"));
	output_msg(sformatf("%-15s  %22s  %11s  %11s  %11s\n\n",
			   "Solid solution", "Component", "Moles", "Delta moles",
			   "Mole fract"));
	/*
	 *   Print solid solutions
	 */
	std::vector<cxxSS *> ss_ptrs = use.Get_ss_assemblage_ptr()->Vectorize();
	for (j = 0; j < (int) ss_ptrs.size(); j++)
	{
		cxxSS * ss_ptr = ss_ptrs[j];
		if (ss_ptr->Get_ss_in())
		{
			/* solid solution name, moles */
			output_msg(sformatf("%-15s  %22s  %11.2e\n",
					   ss_ptr->Get_name().c_str(), "  ",
					   (double)  ss_ptr->Get_total_moles()));
			/* component name, moles, delta moles, mole fraction */

			for (i = 0; i < (int) ss_ptr->Get_ss_comps().size(); i++)
			{
				cxxSScomp *comp_ptr = &(ss_ptr->Get_ss_comps()[i]);
				if (state != TRANSPORT && state != PHAST)
				{
					delta_moles =
						comp_ptr->Get_moles() -
						comp_ptr->Get_initial_moles() -
						comp_ptr->Get_delta();
				}
				else
				{
					delta_moles =
						comp_ptr->Get_moles() -
						comp_ptr->Get_init_moles();
				}
				output_msg(sformatf(
						   "%15s  %22s  %11.2e  %11.2e  %11.2e\n", " ",
						   comp_ptr->Get_name().c_str(),
						   (double) comp_ptr->Get_moles(), (double) delta_moles,
						   (double) (comp_ptr->Get_moles() /
									 ss_ptr->Get_total_moles())));
			}
			if (ss_ptr->Get_miscibility())
			{
				cxxSScomp *comp0_ptr = &(ss_ptr->Get_ss_comps()[0]);
				cxxSScomp *comp1_ptr = &(ss_ptr->Get_ss_comps()[1]);
				nc = comp0_ptr->Get_moles();
				nb = comp1_ptr->Get_moles();
				xb = nb / (nb + nc);
				xb1 = ss_ptr->Get_xb1();
				xb2 = ss_ptr->Get_xb2();

				if (xb > xb1 && xb < xb2)
				{
					xb2moles = (xb1 - 1) / xb1 * nb + nc;
					xb2moles = xb2moles / ((xb1 - 1) / xb1 * xb2 + (1 - xb2));
					xb1moles = (nb - xb2moles * xb2) / xb1;
					output_msg(sformatf(
							   "\n%14s  Solid solution is in miscibility gap\n",
							   " "));
					output_msg(sformatf(
							   "%14s  End members in pct of %s\n\n", " ",
							   comp1_ptr->Get_name().c_str()));
					output_msg(sformatf("%22s  %11g pct  %11.2e\n",
							   " ", (double) xb1, (double) xb1moles));
					output_msg(sformatf("%22s  %11g pct  %11.2e\n",
							   " ", (double) xb2, (double) xb2moles));
				}
			}
		}
		else
		{
			/* solid solution name, moles */
			output_msg(sformatf("%-15s  %22s  %11.2e\n",
					   ss_ptr->Get_name().c_str(), "  ",
					   (double) 0.0));
			/* component name, moles, delta moles, mole fraction */
			for (i = 0; i < (int) ss_ptr->Get_ss_comps().size(); i++)
			{
				cxxSScomp *comp_ptr = &(ss_ptr->Get_ss_comps()[i]);
				if (state != TRANSPORT && state != PHAST)
				{
					delta_moles =
						comp_ptr->Get_moles() -
						comp_ptr->Get_initial_moles() -
						comp_ptr->Get_delta();
				}
				else
				{
					delta_moles =
						comp_ptr->Get_moles() -
						comp_ptr->Get_init_moles();
				}
				output_msg(sformatf(
						   "%15s  %22s  %11.2e  %11.2e  %11.2e\n", " ",
						   comp_ptr->Get_name().c_str(),
						   (double) 0, (double) delta_moles, (double) 0));
			}
		}
	}
	output_msg(sformatf("\n"));
	return (OK);
}
/* ---------------------------------------------------------------------- */
 int Phreeqc::
print_reaction(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   prints irreversible reaction as defined and as
 *   relative moles of each element and total amount
 *   of reaction
 */
	cxxReaction *reaction_ptr;
	if (pr.use == FALSE || pr.all == FALSE)
		return (OK);
	if (state < REACTION || use.Get_reaction_in() == FALSE)
		return (OK);
	if (state == TRANSPORT && transport_step == 0)
		return (OK);
	reaction_ptr = use.Get_reaction_ptr();
/*
 *  Print amount of reaction
 */
	output_msg(sformatf("Reaction %d.\t%s\n\n", use.Get_n_reaction_user(),
			   reaction_ptr->Get_description().c_str()));
	output_msg(sformatf(
			   "\t%11.3e moles of the following reaction have been added:\n\n",
			   (double) step_x));
/*
 *  Print reaction
 */
	output_msg(sformatf("\t%-15s%10s\n", " ", "Relative"));
	output_msg(sformatf("\t%-15s%10s\n\n", "Reactant", "moles"));
	cxxNameDouble::const_iterator cit = reaction_ptr->Get_reactantList().begin();
	for ( ; cit != reaction_ptr->Get_reactantList().end(); cit++)
	{
		output_msg(sformatf("\t%-15s%13.5f\n",
				   cit->first.c_str(), (double) cit->second));
	}
	output_msg(sformatf("\n"));
/*
 *   Debug
 */

	output_msg(sformatf("\t%-15s%10s\n", " ", "Relative"));
	output_msg(sformatf("\t%-15s%10s\n", "Element", "moles"));
	cit = reaction_ptr->Get_elementList().begin();
	for ( ; cit != reaction_ptr->Get_elementList().end(); cit++)
	{
		struct element * elt_ptr = element_store(cit->first.c_str());
		assert(elt_ptr);
		output_msg(sformatf("\t%-15s%13.5f\n",
				   elt_ptr->name,
				   (double) cit->second));
	}
	output_msg(sformatf("\n"));
	return (OK);
}
/* ---------------------------------------------------------------------- */
 int Phreeqc::
print_kinetics(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   prints kinetic reaction,
 *   should be called only on final kinetic step
 */
	LDBLE sim_time;
	cxxKinetics *kinetics_ptr;
	if (pr.kinetics == FALSE || pr.all == FALSE)
		return (OK);
	if (state < REACTION)
		return (OK);
	kinetics_ptr = NULL;
	if (use.Get_kinetics_in() == TRUE)
	{
		if (state == TRANSPORT || state == PHAST || state == ADVECTION)
		{
			kinetics_ptr = Utilities::Rxn_find(Rxn_kinetics_map, use.Get_n_kinetics_user());
		}
		else
		{
			kinetics_ptr = Utilities::Rxn_find(Rxn_kinetics_map, -2);
		}
	}
	if (kinetics_ptr == NULL)
		return (OK);
/*
 *   determine time step
 */
	if (state == TRANSPORT || state == PHAST)
	{
		kin_time_x = timest;
	}
	else if (state == ADVECTION)
	{
		kin_time_x = advection_kin_time;
	}
	sim_time = 0.;
	if (run_info.Get_run_cells())
	{
		sim_time = rate_sim_time;
	}
	else
	{
		if (incremental_reactions == TRUE)
		{
			if (!kinetics_ptr->Get_equalIncrements())
			{
				for (int i = 0; i < reaction_step; i++)
				{
					if (i < (int) kinetics_ptr->Get_steps().size())
					{
						sim_time += kinetics_ptr->Get_steps()[i];
					}
					else
					{
						sim_time +=	kinetics_ptr->Get_steps().back();
					}
				}
			}
			else
			{
				if (reaction_step > kinetics_ptr->Get_count())
				{
					sim_time = kinetics_ptr->Get_steps().front();
				}
				else
				{
					sim_time =
						reaction_step * kinetics_ptr->Get_steps().front() /
						((LDBLE) (kinetics_ptr->Get_count()));
				}
			}
		}
	}
/*
 *  Print amount of reaction
 */
	if (phast == FALSE)
	{
		output_msg(sformatf("Kinetics %d.\t%s\n\n",
				   use.Get_n_kinetics_user(), kinetics_ptr->Get_description().c_str()));
	}
	else
	{
		output_msg(sformatf("Kinetics.\n\n"));
	}
/*
 *  Print reaction
 */
	if (state == TRANSPORT)
	{
		output_msg(sformatf("\tTime:      %g seconds\n",
				   (double) (initial_total_time + transport_step * timest)));
		output_msg(sformatf("\tTime step: %g seconds\n\n",
				   (double) kin_time_x));
	}
	else if (state == ADVECTION)
	{
		output_msg(sformatf("\tTime:      %g seconds\n",
				   (double) (initial_total_time +
							 advection_step * advection_kin_time)));
		output_msg(sformatf("\tTime step: %g seconds\n\n",
				   (double) kin_time_x));
	}
	else if (state == PHAST)
	{
		output_msg(sformatf("\tTime:      %g seconds\n",
				   (double) rate_sim_time_end));
		output_msg(sformatf("\tTime step: %g seconds\n\n",
				   (double) kin_time_x));
	}
	else if (state == REACTION)
	{
		if (incremental_reactions == FALSE)
		{
			output_msg(sformatf("\tTime step: %g seconds\n\n",
					   (double) kin_time_x));
		}
		else
		{
			output_msg(sformatf(
					   "\tTime step: %g seconds  (Incremented time: %g seconds)\n\n",
					   (double) kin_time_x, (double) sim_time));
		}
	}
	output_msg(sformatf("\t%-15s%12s%12s   %-15s%12s\n\n",
			   "Rate name", "Delta Moles", "Total Moles", "Reactant",
			   "Coefficient"));
	for (size_t i = 0; i < kinetics_ptr->Get_kinetics_comps().size(); i++)
	{
		cxxKineticsComp *kinetics_comp_ptr = &(kinetics_ptr->Get_kinetics_comps()[i]);
		if (state != TRANSPORT && state != PHAST)
		{
			output_msg(sformatf("\t%-15s%12.3e%12.3e",
					   kinetics_comp_ptr->Get_rate_name().c_str(),
					   (double) -kinetics_comp_ptr->Get_moles(),
					   (double) kinetics_comp_ptr->Get_m()));
		}
		else
		{
			output_msg(sformatf("\t%-15s%12.3e%12.3e",
					   kinetics_comp_ptr->Get_rate_name().c_str(),
					   (double) (kinetics_comp_ptr->Get_m() -
								 kinetics_comp_ptr->Get_initial_moles()),
					   (double) kinetics_comp_ptr->Get_m()));
		}
		cxxNameDouble::iterator it = kinetics_comp_ptr->Get_namecoef().begin();
		for ( ; it != kinetics_comp_ptr->Get_namecoef().end(); it++)
		{
			std::string name = it->first;
			LDBLE coef = it->second;
			if (it == kinetics_comp_ptr->Get_namecoef().begin())
			{
				output_msg(sformatf("   %-15s%12g\n",
						   name.c_str(),
						   (double) coef));
			}
			else
			{
				output_msg(sformatf("\t%39s   %-15s%12g\n", " ",
						    name.c_str(),
						   (double) coef));
			}
		}
	}
	output_msg(sformatf("\n"));
	return (OK);
}
#ifdef SKIP
/* ---------------------------------------------------------------------- */
int Phreeqc::
print_master_reactions(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Debugging print routine to test primary and secondary reactions
 */
	int i;
	struct rxn_token *next_token;

	for (i = 0; i < count_master; i++)
	{
		output_msg(sformatf("%s\t%s\n\tPrimary reaction\n",
				   master[i]->elt->name, master[i]->s->name));
		next_token = master[i]->rxn_primary->token;
		for (; next_token->s != NULL; next_token++)
		{
			output_msg(sformatf("\t\t%s\t%f\n", next_token->s->name,
					   (double) next_token->coef));
		}
		output_msg(sformatf("\n\tSecondary reaction:\n"));
		if (master[i]->rxn_secondary != NULL)
		{
			next_token = master[i]->rxn_secondary->token;
			for (; next_token->s != NULL; next_token++)
			{
				output_msg(sformatf("\t\t%s\t%f\n",
						   next_token->s->name, (double) next_token->coef));
			}
		}
		output_msg(sformatf("\n\tRedox reaction:\n"));
		if (*(master[i]->pe_rxn) != NULL)
		{
			next_token = (*(master[i]->pe_rxn))->token;
			for (; next_token->s != NULL; next_token++)
			{
				output_msg(sformatf("\t\t%s\t%f\n",
						   next_token->s->name, (double) next_token->coef));
			}
		}
		output_msg(sformatf("\n"));
	}
	return (OK);
}
#endif
/* ---------------------------------------------------------------------- */
 int Phreeqc::
print_mix(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   prints definition of mixing, solution number and multiplier
 */
	cxxMix * mix_ptr;
	cxxSolution *solution_ptr;

	if (pr.use == FALSE || pr.all == FALSE)
		return (OK);
	if (use.Get_mix_in() == FALSE || state < REACTION)
		return (OK);
	if (state == TRANSPORT)
	{
		mix_ptr = Utilities::Rxn_find(Rxn_mix_map, use.Get_n_mix_user());
	}
	else
	{
		mix_ptr = Utilities::Rxn_find(Rxn_mix_map, use.Get_n_mix_user_orig());
	}
	if (mix_ptr == NULL)
	{
		mix_ptr = use.Get_mix_ptr();
	}
/*
 *  Print mixture data
 */
	if (mix_ptr == NULL)
	{
		return (OK);

	}
	if (state == TRANSPORT)
	{
		output_msg(sformatf("Mixture %d.\t%s\n\n", use.Get_n_mix_user(),
				   mix_ptr->Get_description().c_str()));
	}
	else
	{
		output_msg(sformatf("Mixture %d.\t%s\n\n", mix_ptr->Get_n_user(),
				   mix_ptr->Get_description().c_str()));
	}
	std::map<int, LDBLE>::const_iterator cit;
	for (cit = mix_ptr->Get_mixComps().begin(); cit != mix_ptr->Get_mixComps().end(); cit++)
	{
		solution_ptr = Utilities::Rxn_find(Rxn_solution_map, cit->first);
		if (solution_ptr == NULL)
		{
			input_error++;
			return (ERROR);
		}
		output_msg(sformatf("\t%11.3e Solution %d\t%-55s\n",
				   (double) cit->second,
				   cit->first, solution_ptr->Get_description().c_str()));
	}
	output_msg(sformatf("\n"));
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
print_reaction(struct reaction *rxn_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Debugging print of individual chemical reactions for
 *   species or phases
 */
	int j;
	struct rxn_token *next_token;

	if (pr.use == FALSE || pr.all == FALSE)
		return (OK);

	output_msg(sformatf("%s\t\n", rxn_ptr->token[0].s->name));
	output_msg(sformatf("\n\tlog k:\n"));
	for (j = 0; j < MAX_LOG_K_INDICES; j++)
	{
		output_msg(sformatf("\t%f", (double) rxn_ptr->logk[j]));
	}
	output_msg(sformatf("\n\nReaction:\n"));
	for (next_token = rxn_ptr->token; next_token->s != NULL; next_token++)
	{
		output_msg(sformatf("\t\t%s\t%f\n", next_token->s->name,
				   (double) next_token->coef));
	}
	output_msg(sformatf("\n"));
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
print_saturation_indices(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints saturation indices of all applicable pure_phases
 */
	int i;
	LDBLE si, iap;
	LDBLE lk;
	LDBLE la_eminus;
	struct rxn_token *rxn_ptr;
	struct reaction *reaction_ptr;
	bool gas = true;

	if (pr.saturation_indices == FALSE || pr.all == FALSE)
		return (OK);
	if (state == INITIAL_SOLUTION)
	{
		iap = 0;
		for (size_t tok = 1; tok < pe_x[default_pe_x].Get_tokens().size(); tok++)
		{
			iap += pe_x[default_pe_x].Get_tokens()[tok].coef * pe_x[default_pe_x].Get_tokens()[tok].s->la;
			/* fprintf(output,"\t%s\t%f\t%f\n", rxn_ptr->s->name, rxn_ptr->coef, rxn_ptr->s->la ); */
		}
		lk = k_calc(pe_x[default_pe_x].Get_logk(), tk_x, patm_x * PASCAL_PER_ATM);
		la_eminus = lk + iap;
		/* fprintf(output,"\t%s\t%f\n", "pe", si ); */
	}
	else
	{
		la_eminus = s_eminus->la;
	}
/* If a fixed pressure gas-phase disappeared, no PR for the SI's of gases... */
	if (use.Get_gas_phase_ptr() != NULL)
	{
		cxxGasPhase * gas_phase_ptr = use.Get_gas_phase_ptr();
		if (gas_phase_ptr->Get_type() == cxxGasPhase::GP_PRESSURE)
		{
			if (gas_unknown == NULL || gas_unknown->moles < 1e-12)
				gas = false;
		}
	}

/*
 *   Print heading
 */
	print_centered("Saturation indices");
	output_msg(sformatf("  %-15s%9s%8s%9s%3d%4s%3d%4s\n\n", "Phase", "SI**",
			   "log IAP", "log K(", int(tk_x), " K, ", int(floor(patm_x + 0.5)), " atm)"));

	for (i = 0; i < count_phases; i++)
	{
		if (phases[i]->in == FALSE || phases[i]->type != SOLID)
			continue;
		/* check for solids and gases in equation */
		if (phases[i]->replaced)
			reaction_ptr = phases[i]->rxn_s;
		else
			reaction_ptr = phases[i]->rxn;
/*
 *   Print saturation index
 */
		reaction_ptr->logk[delta_v] = calc_delta_v(reaction_ptr, true) -
			 phases[i]->logk[vm0];
		if (reaction_ptr->logk[delta_v])
				mu_terms_in_logk = true;
		lk = k_calc(reaction_ptr->logk, tk_x, patm_x * PASCAL_PER_ATM);
		iap = 0.0;
		for (rxn_ptr = reaction_ptr->token + 1; rxn_ptr->s != NULL;
			 rxn_ptr++)
		{
			if (rxn_ptr->s != s_eminus)
			{
				iap += (rxn_ptr->s->lm + rxn_ptr->s->lg) * rxn_ptr->coef;
			}
			else
			{
				iap += la_eminus * rxn_ptr->coef;
			}
		}
		si = -lk + iap;

		output_msg(sformatf("  %-15s%7.2f  %8.2f%8.2f  %s",
				   phases[i]->name, (double) si, (double) iap, (double) lk,
				   phases[i]->formula));
		if (gas && phases[i]->pr_in && phases[i]->pr_p)
		{
			if (phases[i]->moles_x || state == INITIAL_SOLUTION)
			{
				output_msg(sformatf("\t%s%5.1f%s%5.3f",
					    " Pressure ", (double) phases[i]->pr_p, " atm, phi ", (double) phases[i]->pr_phi));
			} else
			{
				for (int j = 0; j < count_unknowns; j++)
				{
					if (x[j]->type != PP)
						continue;
					if (!strcmp(x[j]->phase->name, phases[i]->name))
					{
						if (x[j]->moles)
							output_msg(sformatf("\t%s%5.1f%s%5.3f",
								" Pressure ", (double) phases[i]->pr_p, " atm, phi ", (double) phases[i]->pr_phi));
						break;
					}
				}
			}
		}
		phases[i]->pr_in = false;
		output_msg("\n");
	}
	output_msg(sformatf("\n%s\n%s",
		"**For a gas, SI = log10(fugacity). Fugacity = pressure * phi / 1 atm.",
		"  For ideal gases, phi = 1."));
	output_msg("\n\n");

	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
print_pp_assemblage(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints saturation indices and masses of pure_phases in pp_assemblage
 */
	int j, k;
	LDBLE si, iap, lk;
	char token[MAX_LENGTH];
	struct rxn_token *rxn_ptr;
	struct phase *phase_ptr;

	if (pr.pp_assemblage == FALSE || pr.all == FALSE)
		return (OK);
	if (pure_phase_unknown == NULL)
		return (OK);
/*
 *   Print heading
 */
	print_centered("Phase assemblage");
	output_msg(sformatf("%73s\n", "Moles in assemblage"));
	output_msg(sformatf("%-14s%8s%2s%7s  %11s", "Phase", "SI", "  ", "log IAP",
			   "log K(T, P)"));
	output_msg(sformatf("  %8s%12s%12s", " Initial", " Final",
			   " Delta"));
	output_msg("\n\n");

	for (j = 0; j < count_unknowns; j++)
	{
		if (x[j]->type != PP)
			continue;
		//cxxPPassemblageComp * comp_ptr = pp_assemblage_ptr->Find(x[j]->pp_assemblage_comp_name);
		cxxPPassemblageComp * comp_ptr = (cxxPPassemblageComp * ) x[j]->pp_assemblage_comp_ptr;
/*
 *   Print saturation index
 */
		iap = 0.0;
		phase_ptr = x[j]->phase;
		if (x[j]->phase->rxn_x == NULL || phase_ptr->in == FALSE)
		{
			output_msg(sformatf("%-18s%23s", x[j]->phase->name,
					   "Element not present."));
		}
		else
		{
			phase_ptr = x[j]->phase;
			phase_ptr->rxn->logk[delta_v] = calc_delta_v(phase_ptr->rxn, true) -
				phase_ptr->logk[vm0];
			if (phase_ptr->rxn->logk[delta_v])
				mu_terms_in_logk = true;
			lk = k_calc(phase_ptr->rxn->logk, tk_x, patm_x * PASCAL_PER_ATM);
			for (rxn_ptr = phase_ptr->rxn->token + 1; rxn_ptr->s != NULL;
				 rxn_ptr++)
			{
				if (rxn_ptr->s != s_eminus)
				{
					iap += (rxn_ptr->s->lm + rxn_ptr->s->lg) * rxn_ptr->coef;
				}
				else
				{
					iap += s_eminus->la * rxn_ptr->coef;
				}
			}
			si = -lk + iap;
			/*
			   for (rxn_ptr = x[j]->phase->rxn_x->token + 1; rxn_ptr->s != NULL; rxn_ptr++) {
			   iap += rxn_ptr->s->la * rxn_ptr->coef;
			   }
			   si = -x[j]->phase->lk + iap;
			   output_msg(OUTPUT_MESSAGE,"\t%-15s%7.2f%8.2f%8.2f", x[j]->phase->name, (double) si, (double) iap, (double) x[j]->phase->lk);
			 */
			output_msg(sformatf("%-14s%8.2f  %7.2f  %8.2f",
					   x[j]->phase->name, (double) si, (double) iap, (double) lk));
		}
/*
 *   Print pure phase assemblage data
 */
		if (x[j]->moles < 0.0)
			x[j]->moles = 0.0;
		if (state != TRANSPORT && state != PHAST)
		{
			sprintf(token, "  %11.3e %11.3e %11.3e",
					(double) (comp_ptr->Get_moles() +
							  comp_ptr->Get_delta()), (double) x[j]->moles,
					(double) (x[j]->moles - comp_ptr->Get_moles() -
							  comp_ptr->Get_delta()));
		}
		else
		{
			sprintf(token, " %11.3e %11.3e %11.3e",
					(double) comp_ptr->Get_initial_moles(),
					(double) x[j]->moles,
					(double) (x[j]->moles - comp_ptr->Get_initial_moles()));
		}
		if (x[j]->moles <= 0.0)
		{
			for (k = 0; k < 11; k++)
			{
				token[13 + k] = ' ';
			}
		}
		if (comp_ptr->Get_add_formula().size() == 0)
		{
			output_msg(sformatf("%37s\n", token));
		}
		else
		{
			output_msg(sformatf("\n	 %-18s%-15s%36s\n",
					   comp_ptr->Get_add_formula().c_str(), " is reactant", token));
		}
	}
	output_msg("\n");
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
print_species(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints description of solution, uses array species_list for
 *   order of aqueous species.
 */
	int i;
	const char *name, *name1;
	struct master *master_ptr;
	LDBLE min;
	LDBLE lm;

	if (pr.species == FALSE || pr.all == FALSE)
		return (OK);
	min = -1000;
	print_centered("Distribution of species");
/*
 *   Heading for species
 */
	if (pitzer_model == TRUE)
	{
		if (ICON == TRUE)
		{
			output_msg(sformatf("%60s%10s\n", "MacInnes", "MacInnes"));
			output_msg(sformatf("%40s%10s%10s%10s%10s\n",
					   "MacInnes", "Log", "Log", "Log", "mole V"));
		}
		else
		{
			output_msg(sformatf("%60s%10s\n", "Unscaled", "Unscaled"));
			output_msg(sformatf("%40s%10s%10s%10s%10s\n",
					   "Unscaled", "Log", "Log", "Log", "mole V"));
		}
	}
	else
	{
		output_msg(sformatf("%50s%10s%10s%10s\n", "Log", "Log", "Log", "mole V"));
	}
	output_msg(sformatf("   %-13s%12s%12s%10s%10s%10s%10s\n\n", "Species",
#ifdef NO_UTF8_ENCODING
			   "Molality", "Activity", "Molality", "Activity", "Gamma", "cm3/mol"));
#else
			   "Molality", "Activity", "Molality", "Activity", "Gamma", "cm³/mol"));
#endif
/*
 *   Print list of species
 */
	s_h2o->lm = s_h2o->la;
	name = s_hplus->secondary->elt->name;
	for (i = 0; i < count_species_list; i++)
	{
/*
 *   Get name of master species
 */
		if (species_list[i].s->type == EX)
			continue;
		if (species_list[i].s->type == SURF)
			continue;
		if (species_list[i].master_s->secondary != NULL)
		{
			master_ptr = species_list[i].master_s->secondary;
			name1 = species_list[i].master_s->secondary->elt->name;
		}
		else
		{
			master_ptr = species_list[i].master_s->primary;
			name1 = species_list[i].master_s->primary->elt->name;
		}
/*
 *   Check if new master species, print total molality
 */
		if (name1 != name)
		{
			name = name1;
			output_msg(sformatf("%-11s%12.3e\n", name,
					   (double) (master_ptr->total / mass_water_aq_x)));
			min = censor * master_ptr->total / mass_water_aq_x;
			if (min > 0)
			{
				min = log10(min);
			}
			else
			{
				min = -1000.;
			}
		}
/*
 *   Print species data
 */
		if (species_list[i].s->lm > min)
		{
			if (species_list[i].s == s_h2o)
			{
				lm = log10(s_h2o->moles / mass_water_aq_x);
			}
			else
			{
				lm = species_list[i].s->lm;
			}
			output_msg(sformatf(
					   "   %-13s%12.3e%12.3e%10.3f%10.3f%10.3f",
					   species_list[i].s->name,
					   (double) ((species_list[i].s->moles) /
								 mass_water_aq_x),
					   (double) under(species_list[i].s->lm +
									  species_list[i].s->lg), (double) lm,
					   (double) (species_list[i].s->lm +
								 species_list[i].s->lg),
					   (double) species_list[i].s->lg));
			//if (species_list[i].s->logk[vm_tc] || !strcmp(species_list[i].s->name, "H+"))
			if (species_list[i].s->logk[vm_tc] || species_list[i].s == s_hplus)
				output_msg(sformatf("%10.2f\n",
					   (double) species_list[i].s->logk[vm_tc]));
			else
				output_msg(sformatf("     (0)  \n"));
		}
	}
	output_msg(sformatf("\n"));
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
print_surface(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints description of surface, including charge and potential,
 *   grams and specific area, moles of each species on surface sites,
 *   and description of diffuse layer if applicable.
 */
	cxxSurface *surface_ptr;
	std::string name, token;
	struct master *master_ptr;
	LDBLE molfrac, charge;
/*
 *  Print surface speciation
 */
	surface_ptr = use.Get_surface_ptr();
	if (surface_ptr == NULL || pr.surface == FALSE || pr.all == FALSE)
		return (OK);
	if (surface_ptr->Get_type() == cxxSurface::CD_MUSIC)
		return (print_surface_cd_music());

	if (state >= REACTION)
	{
		print_centered("Surface composition");
	}
/*
 *   Print list of species
 */

	s_h2o->lm = s_h2o->la;
	if (use.Get_surface_ptr()->Get_type() == cxxSurface::DDL)
	{
		output_msg(sformatf("%-14s\n", "Diffuse Double Layer Surface-Complexation Model\n"));
	}
	else if (use.Get_surface_ptr()->Get_type() == cxxSurface::CCM)
	{
		output_msg(sformatf("%-14s\n", "Constant Capacitance Surface-Complexation Model\n"));
	}
	for (int j = 0; j < count_unknowns; j++)
	{
		/*if (use.Get_surface_ptr()->edl == TRUE) { */
		if (use.Get_surface_ptr()->Get_type() == cxxSurface::DDL || use.Get_surface_ptr()->Get_type() == cxxSurface::CCM)
		{
			if (x[j]->type != SURFACE_CB)
				continue;
			name = x[j]->master[0]->elt->name;
			Utilities::replace("_psi", "", name);
		}
		else
		{
			if (x[j]->type != SURFACE)
				continue;
			token = x[j]->master[0]->elt->name;
			Utilities::replace("_", " ", token);
			std::string::iterator b = token.begin();
			std::string::iterator e = token.end();
			CParser::copy_token(name, b, e);
		}
		output_msg(sformatf("%-14s\n", name.c_str()));
/*
 *   Description of surface
 */
		if (dl_type_x != cxxSurface::NO_DL)
		{
			output_msg(sformatf(
					   "\t%11.3e  Surface + diffuse layer charge, eq\n",
					   (double) x[j]->f));
		}
		/*if (use.Get_surface_ptr()->edl == TRUE && diffuse_layer_x == FALSE) { */
		if ((use.Get_surface_ptr()->Get_type() == cxxSurface::DDL || use.Get_surface_ptr()->Get_type() == cxxSurface::CCM) && dl_type_x == cxxSurface::NO_DL)
		{
			charge = x[j]->f;
		}
		else
		{
			charge = calc_surface_charge(name.c_str());
		}
		output_msg(sformatf("\t%11.3e  Surface charge, eq\n",
				   (double) charge));
		if (x[j]->type == SURFACE_CB)
		{
			cxxSurfaceCharge * charge_ptr = use.Get_surface_ptr()->Find_charge(x[j]->surface_charge);
			if ((charge_ptr->Get_specific_area() *
				 charge_ptr->Get_grams()) > 0)
			{
#ifdef NO_UTF8_ENCODING
				output_msg(sformatf("\t%11.3e  sigma, C/m2\n",
#else
				output_msg(sformatf("\t%11.3e  sigma, C/m²\n",
#endif
						   (double) (charge * F_C_MOL /
									 (charge_ptr->Get_specific_area() *
									  charge_ptr->Get_grams()))));
			}
			else
			{
#ifdef NO_UTF8_ENCODING
				output_msg(sformatf("\tundefined  sigma, C/m2\n"));
#else
				output_msg(sformatf("\tundefined  sigma, C/m²\n"));
#endif
			}
			if (use.Get_surface_ptr()->Get_type() == cxxSurface::CCM)
			{			
				output_msg(sformatf("\t%11.3e  capacitance, F/m^2\n",
					   (double) (charge_ptr->Get_capacitance0())));
			}
			output_msg(sformatf("\t%11.3e  psi, V\n",
					   (double) (x[j]->master[0]->s->la * 2 * R_KJ_DEG_MOL *
								 tk_x * LOG_10 / F_KJ_V_EQ)));
			output_msg(sformatf("\t%11.3e  -F*psi/RT\n",
					   (double) (x[j]->master[0]->s->la * (-2) * LOG_10)));
			output_msg(sformatf("\t%11.3e  exp(-F*psi/RT)\n",
					   exp(x[j]->master[0]->s->la * (-2) * LOG_10)));
			cxxSurfaceComp * comp_ptr = surface_ptr->Find_comp(x[j]->surface_comp);
			if (comp_ptr->Get_phase_name().size() > 0)
			{
				output_msg(sformatf(
#ifdef NO_UTF8_ENCODING
						   "\t%11.3e  specific area, m2/mol %s\n",
#else
						   "\t%11.3e  specific area, m²/mol %s\n",
#endif
						   (double) charge_ptr->Get_specific_area(),
						   comp_ptr->Get_phase_name().c_str()));
				output_msg(sformatf(
#ifdef NO_UTF8_ENCODING
						   "\t%11.3e  m2 for %11.3e moles of %s\n\n",
#else
						   "\t%11.3e  m² for %11.3e moles of %s\n\n",
#endif
						   (double) (charge_ptr->Get_grams() *
									 charge_ptr->Get_specific_area()),
						   (double) charge_ptr->Get_grams(),
						   comp_ptr->Get_phase_name().c_str()));
			}
			else if (comp_ptr->Get_rate_name().size() > 0)
			{
				output_msg(sformatf(
#ifdef NO_UTF8_ENCODING
						   "\t%11.3e  specific area, m2/mol %s\n",
#else
						   "\t%11.3e  specific area, m²/mol %s\n",
#endif
						   (double) charge_ptr->Get_specific_area(),
						   comp_ptr->Get_rate_name().c_str()));
				output_msg(sformatf(
#ifdef NO_UTF8_ENCODING
						   "\t%11.3e  m2 for %11.3e moles of %s\n\n",
#else
						   "\t%11.3e  m² for %11.3e moles of %s\n\n",
#endif
						   (double) (charge_ptr->Get_grams() *
									 charge_ptr->Get_specific_area()),
						   (double) charge_ptr->Get_grams(),
						   comp_ptr->Get_rate_name().c_str()));
			}
			else
			{
				output_msg(sformatf(
#ifdef NO_UTF8_ENCODING
						   "\t%11.3e  specific area, m2/g\n",
#else
						   "\t%11.3e  specific area, m²/g\n",
#endif
						   (double) charge_ptr->Get_specific_area()));
#ifdef NO_UTF8_ENCODING
				output_msg(sformatf("\t%11.3e  m2 for %11.3e g\n\n",
#else
				output_msg(sformatf("\t%11.3e  m² for %11.3e g\n\n",
#endif
						   (double) (charge_ptr->Get_specific_area() *
									 charge_ptr->Get_grams()),
						   (double) charge_ptr->Get_grams()));
			}
			if (dl_type_x != cxxSurface::NO_DL)
				print_diffuse_layer(charge_ptr);
			output_msg(sformatf("\n"));
/*
 *   Heading for species
 */
			for (int k = j - 1; k < count_unknowns; k++)
			{
				if (x[k]->type != SURFACE)
					continue;
				if (x[j] != x[k]->potential_unknown)
					continue;
				master_ptr = x[k]->master[0];
				output_msg(sformatf("%-14s\n",
						   x[k]->master[0]->elt->name));
				output_msg(sformatf("\t%11.3e  moles",
						   (double) x[k]->moles));
				cxxSurfaceComp * comp_k_ptr = surface_ptr->Find_comp(x[k]->surface_comp);
				if (comp_k_ptr->Get_phase_name().size() > 0)
				{
					output_msg(sformatf("\t[%g mol/(mol %s)]\n",
							   (double) comp_k_ptr->Get_phase_proportion(),
							   comp_k_ptr->Get_phase_name().c_str()));
				}
				else if (comp_k_ptr->Get_rate_name().size() > 0)
				{
					output_msg(sformatf(
							   "\t[%g mol/(mol kinetic reactant %s)]\n",
							   (double) comp_k_ptr->Get_phase_proportion(),
							   comp_k_ptr->Get_rate_name().c_str()));
				}
				else
				{
					output_msg(sformatf("\n"));
				}
				output_msg(sformatf("\t%-15s%12s%12s%12s%12s\n", " ",
						   " ", "Mole", " ", "Log"));
				output_msg(sformatf("\t%-15s%12s%12s%12s%12s\n\n",
						   "Species", "Moles", "Fraction", "Molality",
						   "Molality"));
				for (int i = 0; i < count_species_list; i++)
				{
					if (species_list[i].master_s != master_ptr->s)
						continue;
/*
 *   Print species data
 */
					if (x[k]->moles >= MIN_RELATED_SURFACE)
					{
						molfrac =
							(LDBLE) (species_list[i].s->moles) / x[k]->moles *
							species_list[i].s->equiv;
					}
					else
					{
						molfrac = 0.0;
					}
					output_msg(sformatf(
							   "\t%-15s%12.3e%12.3f%12.3e%12.3f\n",
							   species_list[i].s->name,
							   (double) species_list[i].s->moles,
							   (double) molfrac,
							   (double) (species_list[i].s->moles /
										 mass_water_aq_x),
							   log10(species_list[i].s->moles /
									 mass_water_aq_x)));
				}
				output_msg(sformatf("\n"));
			}
		}
		else
		{
			int k = j;
			master_ptr = x[k]->master[0];
			output_msg(sformatf("%-14s\n", x[k]->master[0]->elt->name));
			output_msg(sformatf("\t%11.3e  moles\n",
					   (double) x[k]->moles));
			output_msg(sformatf("\t%-15s%12s%12s%12s%12s\n", " ", " ",
					   "Mole", " ", "Log"));
			output_msg(sformatf("\t%-15s%12s%12s%12s%12s\n\n",
					   "Species", "Moles", "Fraction", "Molality",
					   "Molality"));
			for (int i = 0; i < count_species_list; i++)
			{
				if (species_list[i].master_s != master_ptr->s)
					continue;
/*
 *   Print species data
 */
				if (x[k]->moles >= MIN_RELATED_SURFACE)
				{
					molfrac =
						(double) (species_list[i].s->moles) / x[k]->moles *
						species_list[i].s->equiv;
				}
				else
				{
					molfrac = 0.0;
				}
				output_msg(sformatf(
						   "\t%-15s%12.3e%12.3f%12.3e%12.3f\n",
						   species_list[i].s->name,
						   (double) species_list[i].s->moles,
						   (double) molfrac,
						   (double) (species_list[i].s->moles /
									 mass_water_aq_x),
						   log10(species_list[i].s->moles / mass_water_aq_x)));
			}
			output_msg(sformatf("\n"));
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
print_surface_cd_music(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints description of cd music surfaces, including charge and potential,
 *   grams and specific area, moles of each species on surface sites,
 *   and description of diffuse layer if applicable.
 */
	cxxSurface *surface_ptr;
	std::string name;
	struct master *master_ptr, *master_ptr0, *master_ptr1, *master_ptr2;
	struct unknown *unknown_ptr0, *unknown_ptr1, *unknown_ptr2;
	LDBLE molfrac, charge0, charge1, charge2, sum;
/*
 *  Print surface speciation
 */
	surface_ptr = use.Get_surface_ptr();
	if (surface_ptr == NULL || pr.surface == FALSE || pr.all == FALSE)
		return (OK);

	if (state >= REACTION)
	{
		print_centered("Surface composition");
	}
/*
 *   Print list of species
 */

	s_h2o->lm = s_h2o->la;
	for (int j = 0; j < count_unknowns; j++)
	{
		if (x[j]->type != SURFACE_CB)
			continue;
		name = x[j]->master[0]->elt->name;
		Utilities::replace("_psi", "", name);
		output_msg(sformatf("%-14s\n", name.c_str()));
		cxxSurfaceCharge * charge_ptr = use.Get_surface_ptr()->Find_charge(x[j]->surface_charge);
/*
 *   Description of surface
 */
		if (dl_type_x != cxxSurface::NO_DL)
		{
			output_msg(sformatf(
					   "\t%11.3e  Surface + diffuse layer charge, eq\n\n",
					   (double) (x[j + 2]->f + (charge_ptr->Get_sigma0() + charge_ptr->Get_sigma1()) * (charge_ptr->Get_specific_area() * charge_ptr->Get_grams()) / F_C_MOL)));
		}
		master_ptr0 =
			surface_get_psi_master(charge_ptr->Get_name().c_str(), SURF_PSI);
		master_ptr1 =
			surface_get_psi_master(charge_ptr->Get_name().c_str(), SURF_PSI1);
		master_ptr2 =
			surface_get_psi_master(charge_ptr->Get_name().c_str(), SURF_PSI2);
		unknown_ptr0 = x[master_ptr0->unknown->number];
		unknown_ptr1 = x[master_ptr1->unknown->number];
		unknown_ptr2 = x[master_ptr2->unknown->number];

		charge0 = unknown_ptr0->f;
		charge1 = unknown_ptr1->f;
		if (dl_type_x != cxxSurface::NO_DL)
		{
			charge2 =
				charge_ptr->Get_sigma2() *
				(charge_ptr->Get_specific_area() *
				 charge_ptr->Get_grams()) / F_C_MOL;
		}
		else
		{
			charge2 = unknown_ptr2->f;
		}
		sum = 0;
		for (int k = 0; k < x[j]->count_comp_unknowns; k++)
		{
			sum +=
				x[j]->comp_unknowns[k]->moles *
				x[j]->comp_unknowns[k]->master[0]->s->z;
		}
		output_msg(sformatf("\t%11.3e  Surface charge, plane 0, eq\n",
				   (double) (charge0 + sum)));
		output_msg(sformatf("\t%11.3e  Surface charge, plane 1, eq\n",
				   (double) charge1));
		output_msg(sformatf("\t%11.3e  Surface charge, plane 2, eq\n",
				   (double) charge2));
		output_msg(sformatf(
				   "\t%11.3e  Sum of surface charge, all planes, eq\n\n",
				   (double) (charge0 + sum + charge1 + charge2)));
		if (x[j]->type == SURFACE_CB)
		{
			if ((charge_ptr->Get_specific_area() *
				 charge_ptr->Get_grams()) > 0)
			{
				output_msg(sformatf(
#ifdef NO_UTF8_ENCODING
						   "\t%11.3e  sigma, plane 0, C/m2\n",
#else
						   "\t%11.3e  sigma, plane 0, C/m²\n",
#endif
						   (double) charge_ptr->Get_sigma0()));
				output_msg(sformatf(
#ifdef NO_UTF8_ENCODING
						   "\t%11.3e  sigma, plane 1, C/m2\n",
#else
						   "\t%11.3e  sigma, plane 1, C/m²\n",
#endif
						   (double) charge_ptr->Get_sigma1()));
				output_msg(sformatf(
#ifdef NO_UTF8_ENCODING
						   "\t%11.3e  sigma, plane 2, C/m2\n",
#else
						   "\t%11.3e  sigma, plane 2, C/m²\n",
#endif
						   (double) charge_ptr->Get_sigma2()));
				output_msg(sformatf(
#ifdef NO_UTF8_ENCODING
						   "\t%11.3e  sigma, diffuse layer, C/m2\n\n",
#else
						   "\t%11.3e  sigma, diffuse layer, C/m²\n\n",
#endif
						   (double) charge_ptr->Get_sigmaddl()));
			}
			else
			{
#ifdef NO_UTF8_ENCODING
				output_msg(sformatf("\tundefined  sigma, C/m2\n"));
#else
				output_msg(sformatf("\tundefined  sigma, C/m²\n"));
#endif
			}
			output_msg(sformatf("\t%11.3e  psi, plane 0, V\n",
					   (double) (-master_ptr0->s->la * LOG_10 * R_KJ_DEG_MOL * tk_x / F_KJ_V_EQ)));
			output_msg(sformatf("\t%11.3e  psi, plane 1, V\n",
					   (double) (-master_ptr1->s->la * LOG_10 * R_KJ_DEG_MOL * tk_x / F_KJ_V_EQ)));
			output_msg(sformatf("\t%11.3e  psi, plane 2, V\n\n",
					   (double) (-master_ptr2->s->la * LOG_10 * R_KJ_DEG_MOL * tk_x / F_KJ_V_EQ)));
			output_msg(sformatf("\t%11.3e  exp(-F*psi/RT), plane 0\n",
					   (double) (exp(master_ptr0->s->la * LOG_10))));
			output_msg(sformatf("\t%11.3e  exp(-F*psi/RT), plane 1\n",
					   (double) (exp(master_ptr1->s->la * LOG_10))));
			output_msg(sformatf(
					   "\t%11.3e  exp(-F*psi/RT), plane 2\n\n",
					   (double) (exp(master_ptr2->s->la * LOG_10))));

			output_msg(sformatf("\t%11.3e  capacitance 0-1, F/m^2\n",
					   (double) (charge_ptr->Get_capacitance0())));
			output_msg(sformatf("\t%11.3e  capacitance 1-2, F/m^2\n",
					   (double) (charge_ptr->Get_capacitance1())));
			cxxSurfaceComp * comp_ptr = surface_ptr->Find_comp(x[j]->surface_comp);
			if (comp_ptr->Get_phase_name().size() > 0)
			{
				output_msg(sformatf(
						   "\t%11.3e  specific area, m^2/mol %s\n",
						   (double) charge_ptr->Get_specific_area(),
						   comp_ptr->Get_phase_name().c_str()));
				output_msg(sformatf(
						   "\t%11.3e  m^2 for %11.3e moles of %s\n\n",
						   (double) (charge_ptr->Get_grams() *
									 charge_ptr->Get_specific_area()),
						   (double) charge_ptr->Get_grams(),
						   comp_ptr->Get_phase_name().c_str()));
			}
			else if (comp_ptr->Get_rate_name().size() > 0)
			{
				output_msg(sformatf(
						   "\t%11.3e  specific area, m^2/mol %s\n",
						   (double) charge_ptr->Get_specific_area(),
						   comp_ptr->Get_rate_name().c_str()));
				output_msg(sformatf(
						   "\t%11.3e  m^2 for %11.3e moles of %s\n\n",
						   (double) (charge_ptr->Get_grams() *
									 charge_ptr->Get_specific_area()),
						   (double) charge_ptr->Get_grams(),
						   comp_ptr->Get_rate_name().c_str()));
			}
			else
			{
				output_msg(sformatf(
						   "\t%11.3e  specific area, m^2/g\n",
						   (double) charge_ptr->Get_specific_area()));
				output_msg(sformatf("\t%11.3e  m^2 for %11.3e g\n\n",
						   (double) (charge_ptr->Get_specific_area() *
									 charge_ptr->Get_grams()),
						   (double) charge_ptr->Get_grams()));
			}
			if (dl_type_x != cxxSurface::NO_DL)
				print_diffuse_layer(charge_ptr);
			output_msg(sformatf("\n"));
/*
 *   Heading for species
 */
			for (int k = j - 1; k < count_unknowns; k++)
			{
				if (x[k]->type != SURFACE)
					continue;
				if (x[j] != x[k]->potential_unknown)
					continue;
				master_ptr = x[k]->master[0];
				output_msg(sformatf("%-14s\n",
						   x[k]->master[0]->elt->name));
				output_msg(sformatf("\t%11.3e  moles",
						   (double) x[k]->moles));
				cxxSurfaceComp * comp_k_ptr = surface_ptr->Find_comp(x[k]->surface_comp);
				if (comp_k_ptr->Get_phase_name().size() > 0)
				{
					output_msg(sformatf("\t[%g mol/(mol %s)]\n",
							   (double) comp_k_ptr->Get_phase_proportion(),
							   comp_k_ptr->Get_phase_name().c_str()));
				}
				else if (comp_k_ptr->Get_rate_name().size() > 0)
				{
					output_msg(sformatf(
							   "\t[%g mol/(mol kinetic reactant %s)]\n",
							   (double) comp_k_ptr->Get_phase_proportion(),
							   comp_k_ptr->Get_rate_name().c_str()));
				}
				else
				{
					output_msg(sformatf("\n"));
				}
				output_msg(sformatf("\t%-20s%12s%12s%12s%12s\n", " ",
						   " ", "Mole", " ", "Log"));
				output_msg(sformatf("\t%-20s%12s%12s%12s%12s\n\n",
						   "Species", "Moles", "Fraction", "Molality",
						   "Molality"));
				for (int i = 0; i < count_species_list; i++)
				{
					if (species_list[i].master_s != master_ptr->s)
						continue;
/*
 *   Print species data
 */
					if (x[k]->moles >= MIN_RELATED_SURFACE)
					{
						molfrac =
							(LDBLE) (species_list[i].s->moles) / x[k]->moles *
							species_list[i].s->equiv;
					}
					else
					{
						molfrac = 0.0;
					}
					output_msg(sformatf(
							   "\t%-20s%12.3e%12.3f%12.3e%12.3f\n",
							   species_list[i].s->name,
							   (double) species_list[i].s->moles,
							   (double) molfrac,
							   (double) (species_list[i].s->moles /
										 mass_water_aq_x),
							   log10(species_list[i].s->moles /
									 mass_water_aq_x)));
				}
				output_msg(sformatf("\n"));
			}
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
print_totals(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print total concentrations of elements, molality and moles.
 */
	int i, pure_water;
	LDBLE EC, dens;

	if (pr.totals == FALSE || pr.all == FALSE)
		return (OK);
	print_centered("Solution composition");
	pure_water = TRUE;
	output_msg(sformatf("\t%-15s%12s%12s\n\n", "Elements", "Molality",
			   "Moles"));
	for (i = 0; i < count_unknowns; i++)
	{
		if (x[i] == alkalinity_unknown)
		{
			output_msg(sformatf("\t%-15s%12.3e%12.3e\n",
					   "Alkalinity",
					   (double) (x[i]->f / mass_water_aq_x),
					   (double) x[i]->f));
			pure_water = FALSE;
		}
		if (x[i] == ph_unknown)
			continue;
		if (x[i] == pe_unknown)
			continue;
		if (x[i] == charge_balance_unknown)
		{
			output_msg(sformatf("\t%-15s%12.3e%12.3e",
					   x[i]->description,
					   (double) (x[i]->sum / mass_water_aq_x),
					   (double) x[i]->sum));
			output_msg(sformatf("  Charge balance\n"));
			pure_water = FALSE;
			continue;
		}
		if (x[i]->type == SOLUTION_PHASE_BOUNDARY)
		{
			output_msg(sformatf("\t%-15s%12.3e%12.3e",
					   x[i]->description,
					   (double) (x[i]->sum / mass_water_aq_x),
					   (double) x[i]->sum));
			output_msg(sformatf("  Equilibrium with %s\n",
					   x[i]->phase->name));
			pure_water = FALSE;
			continue;
		}
		if (x[i]->type == MB)
		{
			output_msg(sformatf("\t%-15s%12.3e%12.3e\n",
					   x[i]->description,
					   (double) (x[i]->sum / mass_water_aq_x),
					   (double) x[i]->sum));
			pure_water = FALSE;
		}
	}

	if (pure_water == TRUE)
	{
		output_msg(sformatf("\t%-15s\n", "Pure water"));
	}
/*
 *   Description of solution
 */
	output_msg(sformatf("\n"));
	print_centered("Description of solution");
/*
 *   pH
 */
	output_msg(sformatf("%45s%7.3f    ", "pH  = ",
			   (double) (-(s_hplus->la))));
	if (ph_unknown == NULL)
	{
		output_msg(sformatf("\n"));
	}
	else if (ph_unknown == charge_balance_unknown)
	{
		output_msg(sformatf("  Charge balance\n"));
	}
	else if (ph_unknown->type == SOLUTION_PHASE_BOUNDARY)
	{
		output_msg(sformatf("  Equilibrium with %s\n",
				   ph_unknown->phase->name));
	}
	else if (ph_unknown->type == ALK)
	{
		output_msg(sformatf("  Adjust alkalinity\n"));
	}
/*
 *    pe
 */
	output_msg(sformatf("%45s%7.3f    ", "pe  = ",
			   (double) (-(s_eminus->la))));
	if (pe_unknown == NULL)
	{
		output_msg(sformatf("\n"));
	}
	else if (pe_unknown == charge_balance_unknown)
	{
		output_msg(sformatf("  Charge balance\n"));
	}
	else if (pe_unknown->type == SOLUTION_PHASE_BOUNDARY)
	{
		output_msg(sformatf("  Equilibrium with %s\n",
				   pe_unknown->phase->name));
	}
	else if (pe_unknown->type == MH)
	{
		output_msg(sformatf("  Adjusted to redox equilibrium\n"));
	}
/*
 *   Others
 */
	EC = calc_SC();
	if (EC > 0)
	{
		//output_msg(sformatf("%36s%i%7s%i\n",
		output_msg(sformatf("%35s%3.0f%7s%i\n",
#ifdef NO_UTF8_ENCODING
				   "Specific Conductance (uS/cm, ", tc_x, "oC)  = ", (int) EC));
#else
				   "Specific Conductance (µS/cm, ", tc_x, "°C)  = ", (int) EC));
#endif
	}
/* VP: Density Start */
	if (print_density)
	{
		dens = calc_dens();
#ifdef NO_UTF8_ENCODING
		output_msg(sformatf("%45s%9.5f", "Density (g/cm3)  = ",
#else
		output_msg(sformatf("%45s%9.5f", "Density (g/cm³)  = ",
#endif
			   (double) dens));
		if (dens > 1.999) output_msg(sformatf("%18s\n", " (Program limit)"));
		else output_msg(sformatf("\n"));
		output_msg(sformatf("%45s%9.5f\n", "     Volume (L)  = ",
			   (double) calc_solution_volume()));
	}
/* VP: Density End */
#ifdef NPP
	if (print_viscosity)
	{
		output_msg(sformatf("%45s%9.5f", "Viscosity (mPa s)  = ",
			   (double) viscos));
		if (tc_x > 200 && !pure_water) 
		{
			output_msg(sformatf("%18s\n", 
#ifdef NO_UTF8_ENCODING
				   " (solute contributions limited to 200 oC)"));
#else
				   " (solute contributions limited to 200 °C)"));
#endif
		}
		else output_msg(sformatf("\n"));
	}
#endif
	output_msg(sformatf("%45s%7.3f\n", "Activity of water  = ",
			   exp(s_h2o->la * LOG_10)));
	output_msg(sformatf("%45s%11.3e\n", "Ionic strength (mol/kgw)  = ",
			   (double) mu_x));
	output_msg(sformatf("%45s%11.3e\n", "Mass of water (kg)  = ",
			   (double) mass_water_aq_x));
	if (alkalinity_unknown == NULL)
	{
		output_msg(sformatf("%45s%11.3e\n",
				   "Total alkalinity (eq/kg)  = ",
				   (double) (total_alkalinity / mass_water_aq_x)));
	}
	if (carbon_unknown == NULL)
	{
		output_msg(sformatf("%45s%11.3e\n",
				   "Total carbon (mol/kg)  = ",
				   (double) (total_carbon / mass_water_aq_x)));
	}
	output_msg(sformatf("%45s%11.3e\n", "Total CO2 (mol/kg)  = ",
			   (double) (total_co2 / mass_water_aq_x)));
#ifdef NO_UTF8_ENCODING
	output_msg(sformatf("%45s%6.2f\n", "Temperature (oC)  = ",
#else
	output_msg(sformatf("%45s%6.2f\n", "Temperature (°C)  = ",
#endif
			   (double) tc_x));

	if (patm_x != 1.0)
	{
		/* only print if different than default */
		output_msg(sformatf("%45s%5.2f\n", "Pressure (atm)  = ",
			(double) patm_x));
	}

	output_msg(sformatf("%45s%11.3e\n", "Electrical balance (eq)  = ",
			   (double) cb_x));
	output_msg(sformatf("%45s%6.2f\n",
			   "Percent error, 100*(Cat-|An|)/(Cat+|An|)  = ",
			   (double) (100 * cb_x / total_ions_x)));
	output_msg(sformatf("%45s%3d\n", "Iterations  = ", iterations));
	if (pitzer_model == TRUE || sit_model == TRUE)
	{
		if (always_full_pitzer == FALSE)
		{
			output_msg(sformatf("%45s%3d\n", "Gamma iterations  = ",
				   gamma_iterations));
		}
		else
		{
			output_msg(sformatf("%45s%3d\n", "Gamma iterations  = ",
				  iterations));
		}
		output_msg(sformatf("%45s%9.5f\n", "Osmotic coefficient  = ",
				    (double) COSMOT));
		if (print_density) output_msg(sformatf("%45s%9.5f\n", "Density of water  = ",
				   (double) DW0));
	}
	output_msg(sformatf("%45s%e\n", "Total H  = ", (double) total_h_x));
	output_msg(sformatf("%45s%e\n", "Total O  = ", (double) total_o_x));
	output_msg(sformatf("\n"));

	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
print_user_print(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print with user defined BASIC print routine
 */
	cxxKinetics *kinetics_ptr;

	char l_command[] = "run";

	if (pr.user_print == FALSE || pr.all == FALSE)
		return (OK);
	if (user_print->commands == NULL)
		return (OK);
	kinetics_ptr = NULL;
	if (use.Get_kinetics_in() == TRUE)
	{
		kinetics_ptr = use.Get_kinetics_ptr();
		if (state == TRANSPORT || state == PHAST || state == ADVECTION)
		{
			use.Set_kinetics_ptr(Utilities::Rxn_find(Rxn_kinetics_map, use.Get_n_kinetics_user()));
		}
		else
		{
			use.Set_kinetics_ptr(Utilities::Rxn_find(Rxn_kinetics_map, -2));
		}
	}
	print_centered("User print");
	if (user_print->new_def == TRUE)
	{
		/*      basic_renumber(user_print->commands, &user_print->linebase, &user_print->varbase, &user_print->loopbase); */
		if (basic_compile
			(user_print->commands, &user_print->linebase,
			 &user_print->varbase, &user_print->loopbase) != 0)
		{
			error_msg("Fatal Basic error in USER_PRINT.", STOP);
		}
		user_print->new_def = FALSE;
	}
	if (basic_run
		(l_command, user_print->linebase, user_print->varbase,
		 user_print->loopbase) != 0)
	{
		error_msg("Fatal Basic error in USER_PRINT.", STOP);
	}
	output_msg(sformatf("\n"));
	if (use.Get_kinetics_in() == TRUE)
	{
		use.Set_kinetics_ptr(kinetics_ptr);
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
print_using(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print entities used in calculation
 */
	cxxMix * mix_ptr;
	cxxSolution *solution_ptr;

	if (pr.use == FALSE || pr.all == FALSE)
		return (OK);
	if (state < REACTION || phast == TRUE)
		return (OK);
/*
 *   Mixture or Solution
 */
	if (use.Get_mix_in() == TRUE)
	{
		if (state == TRANSPORT)
		{
			mix_ptr = Utilities::Rxn_find(Rxn_mix_map, use.Get_n_mix_user());
		}
		else
		{
			mix_ptr = Utilities::Rxn_find(Rxn_mix_map, use.Get_n_mix_user_orig());
		}
		if (mix_ptr == NULL)
		{
			mix_ptr = use.Get_mix_ptr();
		}
		if (mix_ptr != NULL)
		{
			if (state == TRANSPORT)
			{
				output_msg(sformatf("Using mix %d.\t%s\n",
						   use.Get_n_mix_user(), mix_ptr->Get_description().c_str()));
			}
			else
			{
				output_msg(sformatf("Using mix %d.\t%s\n",
						   use.Get_n_mix_user_orig(), mix_ptr->Get_description().c_str()));
			}

		}
	}
	else
	{
		solution_ptr = Utilities::Rxn_find(Rxn_solution_map, use.Get_n_solution_user());
		output_msg(sformatf("Using solution %d.\t%s\n",
				   use.Get_n_solution_user(), solution_ptr->Get_description().c_str()));
	}
/*
 *   Exchange and surface
 */
	if (use.Get_exchange_in())
	{
		cxxExchange *exchange_ptr = Utilities::Rxn_find(Rxn_exchange_map, use.Get_n_exchange_user());
		output_msg(sformatf("Using exchange %d.\t%s\n",
				   use.Get_n_exchange_user(), exchange_ptr->Get_description().c_str()));
	}
	if (use.Get_surface_in())
	{
		cxxSurface *surface_ptr = Utilities::Rxn_find(Rxn_surface_map, use.Get_n_surface_user());
		output_msg(sformatf("Using surface %d.\t%s\n",
				   use.Get_n_surface_user(), surface_ptr->Get_description().c_str()));
	}
	if (use.Get_pp_assemblage_in() == TRUE)
	{
		cxxPPassemblage * pp_assemblage_ptr = Utilities::Rxn_find(Rxn_pp_assemblage_map, use.Get_n_pp_assemblage_user());
		output_msg(sformatf("Using pure phase assemblage %d.\t%s\n",
				   use.Get_n_pp_assemblage_user(), pp_assemblage_ptr->Get_description().c_str()));
	}
	if (use.Get_ss_assemblage_in() == TRUE)
	{
		cxxSSassemblage * ss_assemblage_ptr = Utilities::Rxn_find(Rxn_ss_assemblage_map, use.Get_n_ss_assemblage_user());
		output_msg(sformatf(
				   "Using solid solution assemblage %d.\t%s\n",
				   use.Get_n_ss_assemblage_user(),
				   ss_assemblage_ptr->Get_description().c_str()));
	}
	if (use.Get_gas_phase_in())
	{
		cxxGasPhase * gas_phase_ptr = Utilities::Rxn_find(Rxn_gas_phase_map, use.Get_n_gas_phase_user());
		output_msg(sformatf("Using gas phase %d.\t%s\n",
				   use.Get_n_gas_phase_user(), gas_phase_ptr->Get_description().c_str()));
	}
	if (use.Get_temperature_in())
	{
		cxxTemperature *temperature_ptr = Utilities::Rxn_find(Rxn_temperature_map, use.Get_n_temperature_user());
		output_msg(sformatf("Using temperature %d.\t%s\n",
				   use.Get_n_temperature_user(), temperature_ptr->Get_description().c_str()));
	}
	if (use.Get_pressure_in())
	{
		cxxPressure *pressure_ptr = Utilities::Rxn_find(Rxn_pressure_map, use.Get_n_pressure_user());
		output_msg(sformatf("Using pressure %d.\t%s\n",
				   use.Get_n_pressure_user(), pressure_ptr->Get_description().c_str()));
	}
	if (use.Get_reaction_in())
	{
		if (state != TRANSPORT || transport_step > 0)
		{
			cxxReaction *reaction_ptr = Utilities::Rxn_find(Rxn_reaction_map, use.Get_n_reaction_user());
			output_msg(sformatf("Using reaction %d.\t%s\n",
					   use.Get_n_reaction_user(), reaction_ptr->Get_description().c_str()));
		}
	}
	if (use.Get_kinetics_in())
	{
		cxxKinetics * kinetics_ptr;
		if (state == TRANSPORT || state == PHAST || state == ADVECTION)
		{
			kinetics_ptr = Utilities::Rxn_find(Rxn_kinetics_map, use.Get_n_kinetics_user());
		}
		else
		{
			kinetics_ptr = Utilities::Rxn_find(Rxn_kinetics_map, -2);
		}
		output_msg(sformatf("Using kinetics %d.\t%s\n",
				   use.Get_n_kinetics_user(), kinetics_ptr->Get_description().c_str()));
	}
	output_msg(sformatf("\n"));
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
punch_gas_phase(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints selected gas phase data
 */
	//int i;
	LDBLE p, total_moles, volume;
	LDBLE moles;
	bool PR = false;

	//if (punch.count_gases <= 0)
	if (current_selected_output->Get_gases().size() == 0)
		return (OK);
	p = 0.0;
	total_moles = 0.0;
	volume = 0.0;
	cxxGasPhase * gas_phase_ptr = use.Get_gas_phase_ptr();

	if (gas_unknown != NULL && use.Get_gas_phase_ptr() != NULL)
	{
		if (gas_phase_ptr->Get_v_m() >= 0.01)
		{
			PR = true;
		}
		if (gas_phase_ptr->Get_type() == cxxGasPhase::GP_PRESSURE)
		{
			if (gas_unknown->moles >= 1e-12)
			{
				gas_phase_ptr->Set_total_moles(gas_unknown->moles);
				gas_phase_ptr->Set_volume(
					gas_phase_ptr->Get_total_moles() * R_LITER_ATM * tk_x /
					gas_phase_ptr->Get_total_p());
				if (PR)
				{
					gas_phase_ptr->Set_volume(gas_phase_ptr->Get_v_m() * gas_unknown->moles);
				}
			}
			else
			{
				gas_phase_ptr->Set_volume(0);
			}
		}
		p = gas_phase_ptr->Get_total_p();
		total_moles = gas_phase_ptr->Get_total_moles();
		//volume = total_moles * R_LITER_ATM * tk_x / gas_phase_ptr->Get_total_p();
 		//if (gas_phase_ptr->Get_v_m() > 0.03)
 		//	volume = 0.03 * gas_phase_ptr->Get_total_moles();
		volume = gas_phase_ptr->Get_volume();

	}
	if (!current_selected_output->Get_high_precision())
	{
		fpunchf("pressure", "%12.4e\t", (double) p);
		fpunchf("total mol", "%12.4e\t", (double) total_moles);
		fpunchf("volume", "%12.4e\t", (double) volume);
	}
	else
	{
		fpunchf("pressure", "%20.12e\t", (double) p);
		fpunchf("total mol", "%20.12e\t", (double) total_moles);
		fpunchf("volume", "%20.12e\t", (double) volume);
	}
	for (size_t i = 0; i < current_selected_output->Get_gases().size(); i++)
	{
		moles = 0.0;
		if (gas_phase_ptr != NULL && current_selected_output->Get_gases()[i].second != NULL)
		{
			for (size_t j = 0; j < gas_phase_ptr->Get_gas_comps().size(); j++)
			{
				cxxGasComp *gc_ptr = &(gas_phase_ptr->Get_gas_comps()[j]);
				int k;
				struct phase *phase_ptr = phase_bsearch(gc_ptr->Get_phase_name().c_str() , &k, FALSE);
				if (phase_ptr != current_selected_output->Get_gases()[i].second)
					continue;
				moles = phase_ptr->moles_x;
				if (moles <= MIN_TOTAL)
					moles = 0.0;
				break;
			}
		}
		if (!current_selected_output->Get_high_precision())
		{
			fpunchf(sformatf("g_%s", current_selected_output->Get_gases()[i].first.c_str()), "%12.4e\t", (double) moles);
		}
		else
		{
			fpunchf(sformatf("g_%s", current_selected_output->Get_gases()[i].first.c_str()), "%20.12e\t",
					(double) moles);
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
punch_ss_assemblage(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints solid solution composition if present
 */
	//int j, k;
	int found;
	LDBLE moles;

/*
 *   Print solid solutions
 */
	for (size_t k = 0; k < current_selected_output->Get_s_s().size(); k++)
	{
		found = FALSE;
		if (use.Get_ss_assemblage_ptr() != NULL)
		{
			std::vector<cxxSS *> ss_ptrs = use.Get_ss_assemblage_ptr()->Vectorize();
			for (int j = 0; j < (int) ss_ptrs.size(); j++)
			{
				cxxSS * ss_ptr = ss_ptrs[j];
				for (int i = 0; i < (int) ss_ptr->Get_ss_comps().size(); i++)
				{
					cxxSScomp *comp_ptr = &(ss_ptr->Get_ss_comps()[i]);

					if (strcmp_nocase(current_selected_output->Get_s_s()[k].first.c_str(), comp_ptr->Get_name().c_str()) == 0)
					{
						if (ss_ptr->Get_ss_in())
						{
							moles =	comp_ptr->Get_moles();
						}
						else
						{
							moles = 0;
						}
						if (!current_selected_output->Get_high_precision())
						{
							fpunchf(sformatf("s_%s", current_selected_output->Get_s_s()[k].first.c_str()),
									"%12.4e\t", (double) moles);
						}
						else
						{
							fpunchf(sformatf("s_%s", current_selected_output->Get_s_s()[k].first.c_str()),
									"%20.12e\t", (double) moles);
						}
						found = TRUE;
						break;
					}
				}
				if (found == TRUE)
					break;
			}
		}
		if (found == FALSE)
		{
			if (!current_selected_output->Get_high_precision())
			{
				fpunchf(sformatf("s_%s", current_selected_output->Get_s_s()[k].first.c_str()), "%12.4e\t", (double) 0.0);
			}
			else
			{
				fpunchf(sformatf("s_%s", current_selected_output->Get_s_s()[k].first.c_str()), "%20.12e\t",
						(double) 0.0);
			}
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
punch_totals(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print total concentrations of elements, molality and moles.
 */
	//int j;
	LDBLE molality;

	for (size_t j = 0; j < current_selected_output->Get_totals().size(); j++)
	{
		if (current_selected_output->Get_totals()[j].second == NULL)
		{
			molality = 0.0;
		}
		else if (((struct master *) current_selected_output->Get_totals()[j].second)->primary == TRUE)
		{
			if (strncmp(current_selected_output->Get_totals()[j].first.c_str(), "Alkalinity", 20) == 0)
			{
				molality = total_alkalinity / mass_water_aq_x;
			} else
			{
				molality = ((struct master *) current_selected_output->Get_totals()[j].second)->total_primary / mass_water_aq_x;
			}
		}
		else
		{
			molality = ((struct master *) current_selected_output->Get_totals()[j].second)->total / mass_water_aq_x;
		}
		if (!current_selected_output->Get_high_precision())
		{
			fpunchf(sformatf("%s(mol/kgw)", current_selected_output->Get_totals()[j].first.c_str()),
					"%12.4e\t", (double) molality);
		}
		else
		{
			fpunchf(sformatf("%s(mol/kgw)", current_selected_output->Get_totals()[j].first.c_str()),
					"%20.12e\t", (double) molality);
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
punch_molalities(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print concentrations of species (aq, ex, surf)
 */
	//int j;
	LDBLE molality;

	for (size_t j = 0; j < current_selected_output->Get_molalities().size(); j++)
	{
		molality = 0.0;
		if (current_selected_output->Get_molalities()[j].second != NULL
			&& ((struct species *) current_selected_output->Get_molalities()[j].second)->in == TRUE)
		{
			molality = ((struct species *) current_selected_output->Get_molalities()[j].second)->moles / mass_water_aq_x;
		}
		if (!current_selected_output->Get_high_precision())
		{
			fpunchf(sformatf("m_%s(mol/kgw)", current_selected_output->Get_molalities()[j].first.c_str()),
					"%12.4e\t", (double) molality);
		}
		else
		{
			fpunchf(sformatf("m_%s(mol/kgw)", current_selected_output->Get_molalities()[j].first.c_str()),
					"%20.12e\t", (double) molality);
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
punch_activities(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print concentrations of species (aq, ex, surf)
 */
	//int j;
	LDBLE la;

	for (size_t j = 0; j < current_selected_output->Get_activities().size(); j++)
	{
		la = -999.999;
		if (current_selected_output->Get_activities()[j].second != NULL
			&& ((struct species *) current_selected_output->Get_activities()[j].second)->in == TRUE)
		{
			/*la = punch.activities[j].s->lm + punch.activities[j].s->lg; */
			la = log_activity(current_selected_output->Get_activities()[j].first.c_str());
		}
		if (!current_selected_output->Get_high_precision())
		{
			fpunchf(sformatf("la_%s", current_selected_output->Get_activities()[j].first.c_str()), "%12.4e\t",
					(double) la);
		}
		else
		{
			fpunchf(sformatf("la_%s", current_selected_output->Get_activities()[j].first.c_str()),
					"%20.12e\t", (double) la);
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
punch_pp_assemblage(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints masses of selected pure_phases in pp_assemblage
 */
	//int i, j;
	LDBLE moles, delta_moles;
	for (size_t i = 0; i < current_selected_output->Get_pure_phases().size(); i++)
	{
		delta_moles = 0;
		moles = 0.0;
		if (current_selected_output->Get_pure_phases()[i].second != NULL)
		{
			for (int j = 0; j < count_unknowns; j++)
			{
				if (x == NULL || x[j]->type != PP)
					continue;
				//cxxPPassemblageComp * comp_ptr = pp_assemblage_ptr->Find(x[j]->pp_assemblage_comp_name);
				cxxPPassemblageComp * comp_ptr = (cxxPPassemblageComp * ) x[j]->pp_assemblage_comp_ptr;
/*
 *   Print pure phase assemblage data
 */
				if (current_selected_output->Get_pure_phases()[i].second != x[j]->phase)
					continue;
				if (state != TRANSPORT && state != PHAST)
				{
					moles = x[j]->moles;
					delta_moles =
						x[j]->moles - comp_ptr->Get_moles() -
						comp_ptr->Get_delta();
				}
				else
				{
					moles = x[j]->moles;
					delta_moles =
						x[j]->moles - comp_ptr->Get_initial_moles();
				}
				break;
			}
		}
		if (!current_selected_output->Get_high_precision())
		{
			fpunchf(current_selected_output->Get_pure_phases()[i].first.c_str(), "%12.4e\t", (double) moles);
			fpunchf(sformatf("d_%s", current_selected_output->Get_pure_phases()[i].first.c_str()), "%12.4e\t",
					(double) delta_moles);
		}
		else
		{
			fpunchf(current_selected_output->Get_pure_phases()[i].first.c_str(), "%20.12e\t", (double) moles);
			fpunchf(sformatf("d_%s", current_selected_output->Get_pure_phases()[i].first.c_str()),
					"%20.12e\t", (double) delta_moles);
		}
	}
	return (OK);
}

#define PHAST_NULL(x)	  (phast ? NULL : x)
/* ---------------------------------------------------------------------- */
int Phreeqc::
punch_identifiers(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   prints series of integers to identify simulation number,
 *   state of calculations, reaction or transport step number,
 *   and temp, ph, pe, and mass of water for each line
 *   of selected output.
 */
	const char *sformat;
	const char *dformat;
	const char *gformat;
	int i;
	char token[MAX_LENGTH];

	//if (punch.in == FALSE)
	//	return (OK);
	if (!current_selected_output->Get_high_precision())
	{
		sformat = "%12s\t";
		dformat = "%12d\t";
		gformat = "%12g\t";
	}
	else
	{
		sformat = "%20s\t";
		dformat = "%20d\t";
		gformat = "%20g\t";
	}

/*
 *   simulation or simul_tr
 */

	if (current_selected_output->Get_sim())
	{
		if (state != TRANSPORT && state != PHAST)
		{
			fpunchf(PHAST_NULL("sim"), dformat, simulation);
		}
		else
		{
			fpunchf(PHAST_NULL("sim"), dformat, simul_tr);
		}
	}
	if (current_selected_output->Get_state())
	{
		switch (state)
		{
		case 0:
			strcpy(token, "init");
			break;
		case 1:
			strcpy(token, "i_soln");
			break;
		case 2:
			strcpy(token, "i_exch");
			break;
		case 3:
			strcpy(token, "i_surf");
			break;
		case 4:
			strcpy(token, "i_gas");
			break;
		case 5:
			strcpy(token, "react");
			break;
		case 6:
			strcpy(token, "inverse");
			break;
		case 7:
			strcpy(token, "advect");
			break;
		case 8:
			strcpy(token, "transp");
			break;
		default:
			strcpy(token, "unknown");
			break;
		}
		fpunchf(PHAST_NULL("state"), sformat, token);

	}
/*
 *   solution number or cell number and time
 */
	if (current_selected_output->Get_soln())
	{
		if (state == TRANSPORT || state == PHAST)
		{
			fpunchf(PHAST_NULL("soln"), dformat, cell);
		}
		else if (state == ADVECTION)
		{
			fpunchf(PHAST_NULL("soln"), dformat, use.Get_n_solution_user());
		}
		else if (state < REACTION)
		{
			fpunchf(PHAST_NULL("soln"), dformat, use.Get_solution_ptr()->Get_n_user());
		}
		else
		{
			if (use.Get_mix_in() == TRUE)
			{
				if (state != TRANSPORT)
				{
					fpunchf(PHAST_NULL("soln"), dformat, use.Get_n_mix_user_orig());
				}
				else
				{
					fpunchf(PHAST_NULL("soln"), dformat, use.Get_n_mix_user());
				}

			}
			else
			{
				fpunchf(PHAST_NULL("soln"), dformat, use.Get_n_solution_user());
			}
		}
	}
	if (current_selected_output->Get_dist())
	{
		if (state == ADVECTION)
		{
			fpunchf(PHAST_NULL("dist_x"), gformat, (double)use.Get_n_solution_user());
		}
		else if (state == TRANSPORT)
		{
			fpunchf(PHAST_NULL("dist_x"), gformat,
					(double) cell_data[cell - 1].mid_cell_x);
		}
		else
		{
			fpunchf(PHAST_NULL("dist_x"), gformat, (double)-99);
		}
	}
	if (current_selected_output->Get_time())
	{
		LDBLE reaction_time = kin_time_x;
		if (state == REACTION && incremental_reactions == TRUE
			&& use.Get_kinetics_ptr() != NULL)
		{
			if (!use.Get_kinetics_ptr()->Get_equalIncrements())
			{
				reaction_time = 0.0;
				for (i = 0; i < reaction_step; i++)
				{
					if (i < (int) use.Get_kinetics_ptr()->Get_steps().size())
					{
						reaction_time += use.Get_kinetics_ptr()->Get_steps()[i];
					}
					else
					{
						reaction_time +=
							use.Get_kinetics_ptr()->Get_steps().back();
					}
				}
			}
			else
			{
				if (reaction_step > use.Get_kinetics_ptr()->Get_count())
				{
					reaction_time = use.Get_kinetics_ptr()->Get_steps().front();
				}
				else
				{
					reaction_time =
						reaction_step * use.Get_kinetics_ptr()->Get_steps().front() /
						((LDBLE) (use.Get_kinetics_ptr()->Get_count()));
				}
			}
		}
		if (state == REACTION)
		{
			fpunchf(PHAST_NULL("time"), gformat, (double) reaction_time);
		}
		else if (state == TRANSPORT || state == PHAST)
		{
			fpunchf(PHAST_NULL("time"), gformat,
					(double) (initial_total_time + rate_sim_time));
		}
		else if (state == ADVECTION)
		{
			if (advection_kin_time_defined == TRUE)
			{
				fpunchf(PHAST_NULL("time"), gformat,
						(double) (initial_total_time + rate_sim_time));
			}
			else
			{
				fpunchf(PHAST_NULL("time"), gformat, (double)advection_step);
			}
		}
		else
		{
			fpunchf(PHAST_NULL("time"), gformat, (double)-99);
		}
	}

/*
 *   reaction or transport step
 */
	if (current_selected_output->Get_step())
	{
		if (state == REACTION)
		{
			fpunchf(PHAST_NULL("step"), dformat, reaction_step);
		}
		else if (state == ADVECTION)
		{
			fpunchf(PHAST_NULL("step"), dformat, advection_step);
		}
		else if (state == TRANSPORT)
		{
			fpunchf(PHAST_NULL("step"), dformat, transport_step);
		}
		else
		{
			fpunchf(PHAST_NULL("step"), dformat, -99);
		}
	}
	if (current_selected_output->Get_ph())
	{
		if (!current_selected_output->Get_high_precision())
		{
			fpunchf("pH", "%12g\t", (double) (-s_hplus->la));
		}
		else
		{
			fpunchf("pH", "%20.12e\t", (double) (-s_hplus->la));
		}
	}
	if (current_selected_output->Get_pe())
	{
		if (!current_selected_output->Get_high_precision())
		{
			fpunchf("pe", "%12g\t", (double) (-s_eminus->la));
		}
		else
		{
			fpunchf("pe", "%20.12e\t", (double) (-s_eminus->la));
		}
	}
	if (current_selected_output->Get_rxn())
	{
		if (state >= REACTION && use.Get_reaction_in() == TRUE)
		{
			if (!current_selected_output->Get_high_precision())
			{
				fpunchf("reaction", "%12.4e\t", (double) step_x);
			}
			else
			{
				fpunchf("reaction", "%20.12e\t", (double) step_x);
			}
		}
		else
		{
			if (!current_selected_output->Get_high_precision())
			{
				fpunchf("reaction", "%12d\t", -99);
			}
			else
			{
				fpunchf("reaction", "%20d\t", -99);
			}
		}
	}
	if (current_selected_output->Get_temp())
	{
		if (!current_selected_output->Get_high_precision())
		{
			fpunchf("temp(C)", "%12.3f\t", (double) tc_x);
		}
		else
		{
			fpunchf("temp(C)", "%20.12e\t", (double) tc_x);
		}
	}
	if (current_selected_output->Get_alk())
	{
		if (!current_selected_output->Get_high_precision())
		{
			fpunchf("Alk(eq/kgw)", "%12g\t",
					(double) (total_alkalinity / mass_water_aq_x));
		}
		else
		{
			fpunchf("Alk(eq/kgw)", "%20.12e\t",
					(double) (total_alkalinity / mass_water_aq_x));
		}
	}
	if (current_selected_output->Get_mu())
	{
		if (!current_selected_output->Get_high_precision())
		{
			fpunchf("mu", "%12g\t", (double) mu_x);
		}
		else
		{
			fpunchf("mu", "%20.12e\t", (double) mu_x);
		}
	}
	if (current_selected_output->Get_water())
	{
		if (!current_selected_output->Get_high_precision())
		{
			fpunchf("mass_H2O", "%12g\t", (double) mass_water_aq_x);
		}
		else
		{
			fpunchf("mass_H2O", "%20.12e\t", (double) mass_water_aq_x);
		}
	}
	if (current_selected_output->Get_charge_balance())
	{
		if (!current_selected_output->Get_high_precision())
		{
			fpunchf("charge(eq)", "%12g\t", (double) cb_x);
		}
		else
		{
			fpunchf("charge(eq)", "%20.12e\t", (double) cb_x);
		}
	}
	if (current_selected_output->Get_percent_error())
	{
		if (!current_selected_output->Get_high_precision())
		{
			fpunchf("pct_err", "%12g\t",
					(double) (100 * cb_x / total_ions_x));
		}
		else
		{
			fpunchf("pct_err", "%20.12e\t",
					(double) (100 * cb_x / total_ions_x));
		}
	}
	punch_flush();
	return (OK);
}

#undef PHAST_NULL
/* ---------------------------------------------------------------------- */
int Phreeqc::
punch_saturation_indices(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints saturation indices of selected phases
 */
	//int i;
	LDBLE si, iap;
	struct rxn_token *rxn_ptr;

	for (size_t i = 0; i < current_selected_output->Get_si().size(); i++)
	{
		if (current_selected_output->Get_si()[i].second == NULL || ((struct phase *) current_selected_output->Get_si()[i].second)->in == FALSE)
		{
			si = -999.999;
		}
		else
		{
/*
 *   Print saturation index
 */
			iap = 0.0;
			for (rxn_ptr = ((struct phase *) current_selected_output->Get_si()[i].second)->rxn_x->token + 1;
				 rxn_ptr->s != NULL; rxn_ptr++)
			{
				iap += rxn_ptr->s->la * rxn_ptr->coef;
			}
			si = -((struct phase *) current_selected_output->Get_si()[i].second)->lk + iap;
		}
		if (!current_selected_output->Get_high_precision())
		{
			fpunchf(sformatf("si_%s", current_selected_output->Get_si()[i].first.c_str()), "%12.4f\t", (double) si);
		}
		else
		{
			fpunchf(sformatf("si_%s", current_selected_output->Get_si()[i].first.c_str()), "%20.12e\t", (double) si);
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
punch_kinetics(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   prints kinetic reaction,
 *   should be called only on final kinetic step
 */
	cxxKinetics *kinetics_ptr;
	LDBLE moles, delta_moles;

	kinetics_ptr = NULL;
	if (use.Get_kinetics_in() == TRUE)
	{
		if (state == TRANSPORT || state == PHAST || state == ADVECTION)
		{
			kinetics_ptr = Utilities::Rxn_find(Rxn_kinetics_map, use.Get_n_kinetics_user());
		}
		else
		{
			kinetics_ptr = Utilities::Rxn_find(Rxn_kinetics_map, -2);
		}
	}
	for (size_t i = 0; i < current_selected_output->Get_kinetics().size(); i++)
	{
		moles = 0.0;
		delta_moles = 0.0;
		if (kinetics_ptr != NULL)
		{
			for (size_t j = 0; j < kinetics_ptr->Get_kinetics_comps().size(); j++)
			{
				cxxKineticsComp * kinetics_comp_ptr = &(kinetics_ptr->Get_kinetics_comps()[j]);
				if (strcmp_nocase
					(current_selected_output->Get_kinetics()[i].first.c_str(),
					 kinetics_comp_ptr->Get_rate_name().c_str()) == 0)
				{
					if (state != TRANSPORT && state != PHAST)
					{
						moles = kinetics_comp_ptr->Get_m();
						delta_moles = - kinetics_comp_ptr->Get_moles();
					}
					else
					{
						moles =  kinetics_comp_ptr->Get_m();
						delta_moles =
							 kinetics_comp_ptr->Get_m() -
							 kinetics_comp_ptr->Get_initial_moles();
					}
					break;
				}
			}
		}
		if (!current_selected_output->Get_high_precision())
		{
			fpunchf(sformatf("k_%s", current_selected_output->Get_kinetics()[i].first.c_str()), "%12.4e\t",
					(double) moles);
			fpunchf(sformatf("dk_%s", current_selected_output->Get_kinetics()[i].first.c_str()), "%12.4e\t",
					(double) delta_moles);
		}
		else
		{
			fpunchf(sformatf("k_%s", current_selected_output->Get_kinetics()[i].first.c_str()), "%20.12e\t",
					(double) moles);
			fpunchf(sformatf("dk_%s", current_selected_output->Get_kinetics()[i].first.c_str()), "%20.12e\t",
					(double) delta_moles);
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
punch_user_punch(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Punch with user defined BASIC print routine
 */
	char l_command[] = "run";

	n_user_punch_index = 0;
	//if (punch.user_punch == FALSE)
	//	return (OK);
	if (current_user_punch == NULL || !current_selected_output->Get_user_punch())
		return OK;

	struct rate * user_punch = current_user_punch->Get_rate();

	if (user_punch->commands == NULL)
		return (OK);
	if (user_punch->new_def == TRUE)
	{
		if (basic_compile
			(user_punch->commands, &user_punch->linebase,
			 &user_punch->varbase, &user_punch->loopbase) != 0)
		{
			error_msg("Fatal Basic error in USER_PUNCH.", STOP);
		}
		user_punch->new_def = FALSE;
	}
	if (basic_run
		(l_command, user_punch->linebase, user_punch->varbase,
		 user_punch->loopbase) != 0)
	{
		error_msg("Fatal Basic error in USER_PUNCH.", STOP);
	}
	return (OK);
}

#if defined PHREEQ98
/* ---------------------------------------------------------------------- */
int Phreeqc::
punch_user_graph(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Graph with user defined BASIC print routine
 */
	char command[] = "run";

	colnr = 0;
/*    if (pr.user_graph == FALSE || pr.all == FALSE) return(OK); */
/*    if (punch.user_punch == FALSE) return(OK); */
/*	  if (punch.in == FALSE) return(OK); */
	if (user_graph->commands == NULL)
		return (OK);
	if (((state == INITIAL_SOLUTION) || (state == INITIAL_EXCHANGE)
		 || (state == INITIAL_SURFACE) || (state == INITIAL_GAS_PHASE))
		&& (graph_initial_solutions == FALSE))
		return (OK);
	if (FirstCallToUSER_GRAPH)
		AddSeries = TRUE;
	if (state == REACTION)
	{
		/*if (reaction_step == 1) AddSeries = TRUE;
		   else AddSeries = FALSE; */
		if (reaction_step == 1 && !connect_simulations)
			AddSeries = TRUE;
		if (reaction_step > 1)
			AddSeries = FALSE;
	}
	if (state == ADVECTION)
	{
		if (advection_step == 0 && graph_initial_solutions == FALSE)
			return (OK);
		if (((chart_type == 1) && (advection_step == punch_ad_modulus)) ||
			((chart_type == 0) && (advection_step != prev_advection_step)))
			AddSeries = TRUE;
		else
			AddSeries = FALSE;
	}
	if (state == TRANSPORT)
	{
		if (transport_step == 0 && graph_initial_solutions == FALSE)
			return (OK);
		if (((chart_type == 1) && (transport_step == punch_modulus)) ||
			((chart_type == 0) && (transport_step != prev_transport_step)))
			AddSeries = TRUE;
		else
			AddSeries = FALSE;
	}
	if (user_graph->new_def == TRUE)
	{
		if (basic_compile
			(user_graph->commands, &user_graph->linebase,
			 &user_graph->varbase, &user_graph->loopbase) != 0)
		{
			error_msg("Fatal Basic error in USER_GRAPH.", STOP);
		}
		user_graph->new_def = FALSE;
	}
	if (basic_run
		(command, user_graph->linebase, user_graph->varbase,
		 user_graph->loopbase) != 0)
	{
		error_msg("Fatal Basic error in USER_GRAPH.", STOP);
	}
	if (state == ADVECTION)
		prev_advection_step = advection_step;
	if (state == TRANSPORT)
		prev_transport_step = transport_step;
	/*if (state == REACTION) prev_reaction_step = reaction_step; */

	if (FirstCallToUSER_GRAPH)
	{
		start_chart(0);
	}

	FirstCallToUSER_GRAPH = FALSE;

	return (OK);
}
#endif
#if defined(MULTICHART)
/* ---------------------------------------------------------------------- */
int Phreeqc::
punch_user_graph(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Graph with user defined BASIC print routine
 */
	char command[] = "run";
	ChartObject *chart = chart_handler.Get_current_chart();
	if (chart == NULL) return OK;

	chart->Set_AddSeries(false);

	if (chart->Get_rate_command_list().size() == 0)
		return (OK);

	// Skip initial calculations if initial_solutions == false
	if (((state == INITIAL_SOLUTION) || (state == INITIAL_EXCHANGE)
		 || (state == INITIAL_SURFACE) || (state == INITIAL_GAS_PHASE))
		&& (chart->Get_graph_initial_solutions() == false))
		return (OK);

	if (state == REACTION)
	{
		/*if (reaction_step == 1) AddSeries = TRUE;
		   else AddSeries = FALSE; */
		if (reaction_step == 1)
			chart->Set_AddSeries(true);
		if (reaction_step > 1)
			chart->Set_AddSeries(false);
	}

	if (state == ADVECTION)
	{

		if (advection_step == 0 && chart->Get_graph_initial_solutions() == false)
			return (OK);
		if (
				((chart->Get_chart_type() == 1) && (advection_step == punch_ad_modulus)) ||
				((chart->Get_chart_type() == 0) && (advection_step != chart->Get_prev_advection_step()))
			)
		{
			chart->Set_AddSeries(true);
		}
		else
			chart->Set_AddSeries(false);
	}
	if (state == TRANSPORT)
	{
		if (transport_step == 0 && chart->Get_graph_initial_solutions() == FALSE)
			return (OK);
		if (
				((chart->Get_chart_type() == 1) && (transport_step == punch_modulus)) ||
				((chart->Get_chart_type() == 0) && (transport_step != chart->Get_prev_transport_step()))
			)
		{
			chart->Set_AddSeries(true);
		}
		else
		{
			chart->Set_AddSeries(false);
		}
	}

	// From cmdplot_xy merged into transport and advection above

	// From plotXY
	if (chart->Get_FirstCallToUSER_GRAPH())
		chart->Set_prev_sim_no(simulation);
	else
	{
		if (simulation != chart->Get_prev_sim_no())
		{
			chart->Set_AddSeries(true);
		}
	}
	chart->Set_prev_sim_no(simulation);
	if (chart->Get_AddSeries() && !chart->Get_connect_simulations())
	{
		chart->Add_new_series();
	}

	chart->Set_colnr(chart->Get_ColumnOffset());
	chart->Initialize_graph_pts();
	if (chart->Get_rate_new_def())
	{
		if (basic_compile
			(chart->Get_user_graph()->commands, &chart->Get_user_graph()->linebase,
			 &chart->Get_user_graph()->varbase, &chart->Get_user_graph()->loopbase) != 0)
		{
			error_msg("Fatal Basic error in USER_GRAPH.", STOP);
		}
		chart->Set_rate_new_def(false);
	}

	// basic_run calculates points for all graph and plotxy curves
	// colnr identifies the curve and is incremented as each Y/Y2 curve point is added
	if (basic_run
		(command, chart->Get_user_graph()->linebase,
			 chart->Get_user_graph()->varbase, chart->Get_user_graph()->loopbase) != 0)
	{
		error_msg("Fatal Basic error in USER_GRAPH.", STOP);
	}
	chart->Finalize_graph_pts();

	if (state == ADVECTION)
		chart->Set_prev_advection_step(advection_step);
	if (state == TRANSPORT)
		chart->Set_prev_transport_step(transport_step);

	if (chart->Get_FirstCallToUSER_GRAPH())
	{
		chart->start_chart();
	}
	chart->Set_new_ug(false);

	chart->Set_FirstCallToUSER_GRAPH(false);

	return (OK);
}
#endif // MULTICHART

char * Phreeqc::
sformatf(const char *format, ...)
{
	bool success = false;
	do
	{
		va_list args;
		va_start(args, format);
		int j = vsnprintf(sformatf_buffer, sformatf_buffer_size, format, args);
		success = (j > 0 && j < (int) sformatf_buffer_size);
		va_end(args);
		if (!success)
		{
			sformatf_buffer_size *= 2;
			sformatf_buffer = (char *) PHRQ_realloc(sformatf_buffer, sformatf_buffer_size * sizeof(char));
			if (sformatf_buffer == NULL) malloc_error();
		}
	}
	while (!success);

	return sformatf_buffer;
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
print_alkalinity(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints description of solution, uses array species_list for
 *   order of aqueous species.
 */
	int i, j;
	struct species_list *alk_list;
	int count_alk_list;
	LDBLE min;

	if (pr.alkalinity == FALSE || pr.all == FALSE)
		return (OK);
	print_centered("Distribution of alkalinity");
	alk_list =
		(struct species_list *)
		PHRQ_malloc((size_t) (count_s_x * sizeof(struct species_list)));
	if (alk_list == NULL)
	{
		malloc_error();
		return (OK);
	}
	j = 0;
	for (i = 0; i < count_s_x; i++)
	{
		if (s_x[i]->alk == 0.0)
			continue;
		alk_list[j].master_s = s_hplus;
		alk_list[j].s = s_x[i];
		alk_list[j].coef = s_x[i]->alk;
		j++;
	}
	count_alk_list = j;
	min = fabs(censor * total_alkalinity / mass_water_aq_x);
	if (count_alk_list > 0)
	{
		output_msg(sformatf("\t%26s%11.3e\n\n",
				   "Total alkalinity (eq/kgw)  = ",
				   (double) (total_alkalinity / mass_water_aq_x)));
		output_msg(sformatf("\t%-15s%12s%12s%10s\n\n", "Species",
				   "Alkalinity", "Molality", "Alk/Mol"));
		qsort(&alk_list[0], (size_t) count_alk_list,
			  (size_t) sizeof(struct species_list), species_list_compare_alk);
		for (i = 0; i < count_alk_list; i++)
		{
			if (fabs
				(alk_list[i].s->alk * (alk_list[i].s->moles) /
				 mass_water_aq_x) < min)
				continue;
			output_msg(sformatf("\t%-15s%12.3e%12.3e%10.2f\n",
					   alk_list[i].s->name,
					   (double) (alk_list[i].s->alk *
								 (alk_list[i].s->moles) / mass_water_aq_x),
					   (double) ((alk_list[i].s->moles) / mass_water_aq_x),
					   (double) (alk_list[i].s->alk)));
		}
	}

	output_msg(sformatf("\n"));
	alk_list = (struct species_list *) free_check_null(alk_list);
	return (OK);
}

