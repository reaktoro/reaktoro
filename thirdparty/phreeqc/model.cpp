#include "Utils.h"
#include "Phreeqc.h"
#include "phqalloc.h"
#include "cxxMix.h"
#include "Exchange.h"
#include "GasPhase.h"
#include "PPassemblage.h"
#include "SSassemblage.h"
#include "Solution.h"

/* ---------------------------------------------------------------------- */
int Phreeqc::
model(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   model is called after the equations have been set up by prep
 *   and initial guesses have been made in set.
 *
 *   Here is the outline of the calculation sequence:
 *	  residuals--residuals are calculated, if small we are done
 *	  sum_jacobian--jacobian is calculated
 *	  ineq--inequality solver is called
 *	  reset--estimates of unknowns revised, if changes are small solution
 *		 has been found, usually convergence is found in residuals.
 *	  gammas--new activity coefficients
 *	  molalities--calculate molalities
 *	  mb_sums--calculate mass-balance sums
 *	  mb_gases--decide if gas_phase exists
 *	  mb_ss--decide if solid_solutions exists
 *	  switch_bases--check to see if new basis species is needed
 *		 reprep--rewrite equations with new basis species if needed
 *		 revise_guesses--revise unknowns to get initial mole balance
 *	  check_residuals--check convergence one last time
 *		 sum_species--calculate sums of elements from species concentrations
 *
 *	  An additional pass through may be needed if unstable phases still exist
 *		 in the phase assemblage.
 */
	int l_kode, return_kode;
	int r;
	int count_infeasible, count_basis_change;
	int debug_model_save;
	int mass_water_switch_save;

	set_inert_moles();
/*	debug_model = TRUE; */
/*	debug_prep = TRUE; */
/*	debug_set = TRUE; */
	/* mass_water_switch == TRUE, mass of water is constant */
	if (pitzer_model == TRUE && sit_model == TRUE)
	{
	  input_error++;
	  error_msg("Cannot use PITZER and SIT data blocks in same run (database + input file).", STOP);
	}
	if (pitzer_model == TRUE)
	{

		l_kode = model_pz();
		unset_inert_moles();
		return l_kode;
	}
	if (sit_model == TRUE)
	{

		l_kode = model_sit();
		unset_inert_moles();
		return l_kode;
	}
	mass_water_switch_save = mass_water_switch;
	if (mass_water_switch_save == FALSE && delay_mass_water == TRUE)
	{
		mass_water_switch = TRUE;
	}
	debug_model_save = debug_model;
	pe_step_size_now = pe_step_size;
	step_size_now = step_size;
#ifdef NPP
	if (!use.Get_kinetics_in()) status(0, NULL);
#else
	status(0, NULL);
#endif
	iterations = 0;
	count_basis_change = count_infeasible = 0;
	stop_program = FALSE;
	remove_unstable_phases = FALSE;
	for (;;)
	{
		mb_gases();
		mb_ss();
		l_kode = 1;
		while ((r = residuals()) != CONVERGED
			   || remove_unstable_phases == TRUE)
		{
#if defined(PHREEQCI_GUI)
			PhreeqcIWait(this);
#endif
			iterations++;
			if (iterations > itmax - 1 && debug_model == FALSE
				&& pr.logfile == TRUE)
			{
				set_forward_output_to_log(TRUE);
				debug_model = TRUE;
			}
			if (debug_model == TRUE)
			{
				output_msg(sformatf(
						   "\nIteration %d\tStep_size = %f\n", iterations,
						   (double) step_size_now));
				output_msg(sformatf( "\t\tPe_step_size = %f\n\n",
						   (double) pe_step_size_now));
			}
/*
 *   Iterations exceeded
 */
			if (iterations > itmax)
			{
				error_string = sformatf( "Maximum iterations exceeded, %d\n",
						itmax);
				warning_msg(error_string);
				stop_program = TRUE;
				break;
			}
/*
 *   Calculate jacobian
 */
			if (state >= REACTION && numerical_deriv)
			{
				//jacobian_sums();
				numerical_jacobian();
			}
			else /* hmm */
			{
				jacobian_sums();
				numerical_jacobian();
			}
/*
 *   Full matrix with pure phases
 */
			if (r == OK || remove_unstable_phases == TRUE)
			{
				return_kode = ineq(l_kode);
				if (return_kode != OK)
				{
					if (debug_model == TRUE)
					{
						output_msg(sformatf(
								   "Ineq had infeasible solution, "
								   "kode %d, iteration %d\n", return_kode,
								   iterations));
					}
					log_msg(sformatf("Ineq had infeasible solution, "
							   "kode %d, iteration %d\n", return_kode,
							   iterations));
					count_infeasible++;
				}
				if (return_kode == 2)
				{
					ineq(0);
				}
				reset();
			}
			gammas(mu_x);
			if (molalities(FALSE) == ERROR)
			{
				revise_guesses();
/*				adjust_step_size(); */
			}
			if (use.Get_surface_ptr() != NULL &&
				use.Get_surface_ptr()->Get_dl_type() != cxxSurface::NO_DL &&
				use.Get_surface_ptr()->Get_related_phases())
				initial_surface_water();
			mb_sums();
			mb_gases();
			mb_ss();
/*
 *   Switch bases if necessary
 */

			if (switch_bases() == TRUE)
			{
				count_basis_change++;
				reprep();
				gammas(mu_x);
				molalities(TRUE);
				if (use.Get_surface_ptr() != NULL &&
					use.Get_surface_ptr()->Get_dl_type() != cxxSurface::NO_DL &&
					use.Get_surface_ptr()->Get_related_phases())
					initial_surface_water();
				revise_guesses();
				mb_sums();
				mb_gases();
				mb_ss();
			}
/* debug
						species_list_sort();
						sum_species();
			print_species();
			print_exchange();
			print_surface();
 */
			if (stop_program == TRUE)
			{
				break;
			}
		}
/*
 *   Check for stop_program
 */

		if (stop_program == TRUE)
		{
			break;
		}
		if (check_residuals() == ERROR)
		{
			stop_program = TRUE;
			break;
		}
		if (remove_unstable_phases == FALSE && mass_water_switch_save == FALSE
			&& mass_water_switch == TRUE)
		{
			log_msg(sformatf(
					   "\nChanging water switch to FALSE. Iteration %d.\n",
					   iterations));
			mass_water_switch = FALSE;
			continue;
		}
		if (remove_unstable_phases == FALSE)
			break;
		if (debug_model == TRUE)
		{
			output_msg(sformatf(
					   "\nRemoving unstable phases. Iteration %d.\n",
					   iterations));
		}
		log_msg(sformatf( "\nRemoving unstable phases. Iteration %d.\n",
				   iterations));
	}
	log_msg(sformatf( "\nNumber of infeasible solutions: %d\n",
			   count_infeasible));
	log_msg(sformatf( "Number of basis changes: %d\n\n",
			   count_basis_change));
	log_msg(sformatf( "Number of iterations: %d\n\n", iterations));
	debug_model = debug_model_save;
	set_forward_output_to_log(FALSE);
	unset_inert_moles();
	if (stop_program == TRUE)
	{
		return (ERROR);
	}
	return (OK);
}

#ifdef SKIP
/* ---------------------------------------------------------------------- */
int Phreeqc::
adjust_step_size(void)
/* ---------------------------------------------------------------------- */
{
/*
 *  Step sizes are cut down if overflow occurs in molalities
 */
	pe_step_size_now -= (pe_step_size_now - 1.0) / 2.0;
	step_size_now -= (step_size_now - 1.) / 2.0;
	if (pe_step_size_now < 1.5)
		pe_step_size_now = 1.5;
	if (step_size_now < 1.5)
		step_size_now = 1.5;
	log_msg(sformatf( "\tNew step sizes: %f\t%f\t%d\n",
			   step_size_now, pe_step_size_now, iterations));
	return (OK);
}
#endif
/* ---------------------------------------------------------------------- */
int Phreeqc::
check_residuals(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks for convergence of all equations, prints any nonconvergence
 *   Sets variable remove_unstable_phases if a phase is present,
 *   but undersaturated (i.e. aragonite in calcite-saturated solution).
 */
	int i, return_value;
	LDBLE epsilon;
	epsilon = convergence_tolerance;

	return_value = OK;
	if (stop_program == TRUE)
	{
		warning_msg
			("The program has failed to converge to a numerical solution.\n\nThe following equations were not satisfied:");
		/*error_msg("The program has failed to converge to a numerical solution.\n\nThe following equations were not satisfied:", CONTINUE); */
	}
	for (i = 0; i < count_unknowns; i++)
	{
		if (x[i]->type == MB || x[i]->type == ALK)
		{
			if (fabs(residual[i]) >= epsilon * x[i]->moles && fabs(residual[i]) > sqrt(fabs(x[i]->moles) * MIN_TOTAL)
				&& x[i]->moles > MIN_TOTAL /* || stop_program == TRUE */ )
			{
				error_string = sformatf(
						"%20s has not converged. Total: %e\tCalculated: "
						"%e\tResidual: %e\n", x[i]->description,
						(double) x[i]->moles, (double) x[i]->f,
						(double) residual[i]);
				error_msg(error_string, CONTINUE);
				if (x[i]->type == ALK)
				{
					error_msg("Is non-carbonate alkalinity "
							  "greater than total alkalinity?\n", CONTINUE);
				}
				return_value = ERROR;
			}
		}
		else if (x[i]->type == SOLUTION_PHASE_BOUNDARY)
		{
			if (fabs(residual[i]) >= epsilon /* || stop_program == TRUE */ )
			{
				error_string = sformatf(
						"%20s solution phase boundary has not converged. "
						"\tResidual: %e\n", x[i]->description,
						(double) residual[i]);
				error_msg(error_string, CONTINUE);
			}
		}
		else if (x[i]->type == CB)
		{
			if (fabs(residual[i]) >=
				epsilon * mu_x *
				mass_water_aq_x /* || stop_program == TRUE */ )
			{
				error_string = sformatf(
						"%20s Charge balance has not converged. "
						"\tResidual: %e\n", x[i]->description,
						(double) residual[i]);
				error_msg(error_string, CONTINUE);
			}
		}
		else if (x[i]->type == MU /*&& pitzer_model == FALSE && sit_model == FALSE*/)
		{
			if (fabs(residual[i]) >=
				epsilon * mu_x *
				mass_water_aq_x /* || stop_program == TRUE */ )
			{
				error_string = sformatf(
						"%20s Ionic strength has not converged. "
						"\tResidual: %e\n", x[i]->description,
						(double) residual[i]);
				error_msg(error_string, CONTINUE);
			}
		}
		else if (x[i]->type == AH2O && pitzer_model == FALSE && sit_model == FALSE)
		{
			if (fabs(residual[i]) >= epsilon /* || stop_program == TRUE */ )
			{
				error_string = sformatf(
						"%20s Activity of water has not converged. "
						"\tResidual: %e\n", x[i]->description,
						(double) residual[i]);
				error_msg(error_string, CONTINUE);
			}
		}
		else if ((x[i]->type == MH
				  && (pitzer_model == FALSE || pitzer_pe == TRUE)))
		{
#define COMBINE
			/*#define COMBINE_CHARGE */
#ifdef COMBINE
#ifndef COMBINE_CHARGE
			if (fabs(residual[i]) >
				epsilon * (x[i]->moles + 2 * mass_oxygen_unknown->moles))
#else
			if (fabs(residual[i]) >
				epsilon * (x[i]->moles + 2 * mass_oxygen_unknown->moles +
						   charge_balance_unknown->moles))
#endif
#else
			if (fabs(residual[i]) >=
				epsilon * x[i]->moles /* || stop_program == TRUE */ )
#endif
			{
				error_string = sformatf(
						"%20s Mass of hydrogen has not converged. "
						"\tResidual: %e\n", x[i]->description,
						(double) residual[i]);
				error_msg(error_string, CONTINUE);
			}
		}
		else if (x[i]->type == MH2O)
		{
			if (mass_water_switch == TRUE)
				continue;
			if (fabs(residual[i]) >=
				0.01 * epsilon * x[i]->moles /* || stop_program == TRUE */ )
			{
				error_string = sformatf(
						"%20s Mass of oxygen has not converged. "
						"\tResidual: %e\n", x[i]->description,
						(double) residual[i]);
				error_msg(error_string, CONTINUE);
			}
		}
		else if (x[i]->type == PP)
		{
			//cxxPPassemblageComp * comp_ptr = pp_assemblage_ptr->Find(x[i]->pp_assemblage_comp_name);
			cxxPPassemblageComp * comp_ptr = (cxxPPassemblageComp *) x[i]->pp_assemblage_comp_ptr;
			if (comp_ptr->Get_add_formula().size() == 0)
			{
				if (x[i]->dissolve_only == TRUE)
				{
					if ((residual[i] > epsilon && x[i]->moles > 0.0)
						||
						((residual[i] < -epsilon
						  && (comp_ptr->Get_initial_moles() - x[i]->moles) >
						  0)))
					{
						log_msg(sformatf(
								   "%20s Dissolve_only pure phase has not converged. \tResidual: %e\n",
								   x[i]->description, (double) residual[i]));
					}
				}
				else
				{
					if ((residual[i] >= epsilon * 100
						 && x[i]->moles > 0.0) /* || stop_program == TRUE */ )
					{
						remove_unstable_phases = TRUE;
						log_msg(sformatf(
								   "%20s Pure phase has not converged. \tResidual: %e\n",
								   x[i]->description, (double) residual[i]));
					}
					else if (residual[i] <= -epsilon)
					{
						error_string = sformatf(
								"%20s Pure phase has not converged. "
								"\tResidual: %e\n", x[i]->description,
								(double) residual[i]);
						error_msg(error_string, CONTINUE);
					}
				}
			}
			else
			{
				if ((fabs(residual[i]) >= epsilon
					 && x[i]->moles > 0.0) /* || stop_program == TRUE */ )
				{
					log_msg(sformatf(
							   "%s, Pure phase has not converged. \tResidual: %e\n",
							   x[i]->description, (double) residual[i]));
					error_string = sformatf(
							"%s, Pure phase with add formula has not converged.\n\t SI may be a local minimum."
							"\tResidual: %e\n", x[i]->description,
							(double) residual[i]);
					warning_msg(error_string);
				}
			}
		}
		else if (x[i]->type == EXCH)
		{
			if (				/* stop_program == TRUE || */
				   (x[i]->moles <= MIN_RELATED_SURFACE
					&& fabs(residual[i]) > epsilon)
				   || (x[i]->moles > MIN_RELATED_SURFACE
					   && (fabs(residual[i]) > epsilon * x[i]->moles)))
			{
				error_string = sformatf(
						"%20s Exchanger mass balance has not converged. "
						"\tResidual: %e\n", x[i]->description,
						(double) residual[i]);
				error_msg(error_string, CONTINUE);
			}
		}
		else if (x[i]->type == SURFACE)
		{
			if (fabs(residual[i]) < ineq_tol && fabs(residual[i]) < 1e-2*x[i]->moles) continue;
			if (				/* stop_program == TRUE || */
				   (x[i]->moles <= MIN_RELATED_SURFACE
					&& fabs(residual[i]) > epsilon)
				   || (x[i]->moles > MIN_RELATED_SURFACE
					   && (fabs(residual[i]) > epsilon * x[i]->moles)))
			{
				error_string = sformatf(
						"%20s Surface mass balance has not converged. "
						"\tResidual: %e\n", x[i]->description,
						(double) residual[i]);
				error_msg(error_string, CONTINUE);
			}
		}
		else if (x[i]->type == SURFACE_CB || x[i]->type == SURFACE_CB1
				 || x[i]->type == SURFACE_CB2)
		{
			cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[i]->surface_charge);
			if ((charge_ptr->Get_grams() > MIN_RELATED_SURFACE
				 && fabs(residual[i]) >
				 epsilon) /* || stop_program == TRUE */ )
			{
				error_string = sformatf(
						"%20s Surface charge/potential has not converged. "
						"\tResidual: %e\n", x[i]->description,
						(double) residual[i]);
				error_msg(error_string, CONTINUE);
			}
		}
		else if (x[i]->type == GAS_MOLES)
		{
			if (gas_in == FALSE)
				continue;
			if (residual[i] >= epsilon
				|| residual[i] <= -epsilon /* || stop_program == TRUE */ )
			{
				error_string = sformatf(
						"%20s Total moles in gas phase has not converged. "
						"\tResidual: %e\n", x[i]->description,
						(double) residual[i]);
				error_msg(error_string, CONTINUE);
			}
		}
		else if (x[i]->type == PITZER_GAMMA)
		{
			if (fabs(residual[i]) > epsilon)
			{
				error_string = sformatf(
						"%20s log gamma not converged.\tResidual: %e\n",
						x[i]->description, (double) residual[i]);
			}
		}
		else if (x[i]->type == SS_MOLES)
		{
			if (x[i]->ss_in == FALSE)
				continue;
			if (x[i]->moles <= MIN_TOTAL_SS)
				continue;
			if (residual[i] >= epsilon
				|| residual[i] <= -epsilon /* || stop_program == TRUE */ )
			{
				error_string = sformatf(
						"%20s Total moles in solid solution has not converged. "
						"\tResidual: %e\n", x[i]->description,
						(double) residual[i]);
				error_msg(error_string, CONTINUE);
			}
		}
	}
	if (remove_unstable_phases == TRUE)
	{
		log_msg(sformatf( "%20sRemoving unstable phases, iteration %d.",
				   " ", iterations));
	}
	return (return_value);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
gammas(LDBLE mu)
/* ---------------------------------------------------------------------- */
{
/*
 *   Calculates gammas and [moles * d(ln gamma)/d mu] for all aqueous
 *   species.
 */
	int i, j;
	int ifirst, ilast;
	LDBLE d1, d2, d3, f, a_llnl, b_llnl, bdot_llnl, log_g_co2, dln_g_co2, c2_llnl;

	LDBLE c1, c2, a, b;
	LDBLE muhalf, equiv;
	/* Initialize */
	if (mu <= 0) mu = 1e-10;
	if (pitzer_model == TRUE)
		return gammas_pz();
	if (sit_model == TRUE)
		return gammas_sit();
	a_llnl = b_llnl = bdot_llnl = log_g_co2 = dln_g_co2 = c2_llnl = 0;
/*
 *   compute temperature dependence of a and b for debye-huckel
 */
	// a and b are calc'd in calc_dielectrics(tc_x, patm_x);
	k_temp(tc_x, patm_x);
	a = DH_A;
	b = DH_B;

	/*
	 *   LLNL temperature dependence
	 */
	if (llnl_count_temp > 0)
	{
		ifirst = 0;
		ilast = llnl_count_temp;
		if (tc_x < llnl_temp[0] || tc_x > llnl_temp[llnl_count_temp - 1])
		{
			error_msg
				("Temperature out of range of LLNL_AQUEOUS_MODEL parameters",
				 STOP);
		}
		for (i = 0; i < llnl_count_temp; i++)
		{
			if (tc_x >= llnl_temp[i])
				ifirst = i;
			if (tc_x <= llnl_temp[i])
			{
				ilast = i;
				break;
			}
		}
		if (ilast == ifirst)
		{
			f = 1;
		}
		else
		{
			f = (tc_x - llnl_temp[ifirst]) / (llnl_temp[ilast] -
											  llnl_temp[ifirst]);
		}
		a_llnl = (1 - f) * llnl_adh[ifirst] + f * llnl_adh[ilast];
		b_llnl = (1 - f) * llnl_bdh[ifirst] + f * llnl_bdh[ilast];
		bdot_llnl = (1 - f) * llnl_bdot[ifirst] + f * llnl_bdot[ilast];
		/*
		 * CO2 activity coefficient
		 */
		log_g_co2 =
			(llnl_co2_coefs[0] + llnl_co2_coefs[1] * tk_x +
			 llnl_co2_coefs[2] / tk_x) * mu - (llnl_co2_coefs[3] +
											   llnl_co2_coefs[4] * tk_x) *
			(mu / (mu + 1));
		log_g_co2 /= LOG_10;
		dln_g_co2 =
			(llnl_co2_coefs[0] + llnl_co2_coefs[1] * tk_x +
			 llnl_co2_coefs[2] / tk_x) - (llnl_co2_coefs[3] +
										  llnl_co2_coefs[4] * tk_x) * (1 /
																	   ((mu +
																		 1) *
																		(mu +
																		 1)));
	}

/*
 *   constants for equations
 */
	muhalf = sqrt(mu);
	c1 = (-a) * LOG_10 * (1.0 /
						  (2 * muhalf * (muhalf + 1.0) * (muhalf + 1.0)) -
						  0.3);
	c2 = -a / (2 * muhalf);
	if (llnl_count_temp > 0)
	{
		c2_llnl = -a_llnl / (2 * muhalf);
	}

/*
 *   Calculate activity coefficients
 */
	for (i = 0; i < count_s_x; i++)
	{
		switch (s_x[i]->gflag)
		{
		case 0:				/* uncharged */
			s_x[i]->lg = s_x[i]->dhb * mu;
			s_x[i]->dg = s_x[i]->dhb * LOG_10 * s_x[i]->moles;
			break;
		case 1:				/* Davies */
			s_x[i]->lg = -s_x[i]->z * s_x[i]->z * a *
				(muhalf / (1.0 + muhalf) - 0.3 * mu);
			s_x[i]->dg = c1 * s_x[i]->z * s_x[i]->z * s_x[i]->moles;
			break;
		case 2:				/* Extended D-H, WATEQ D-H */
			s_x[i]->lg = -a * muhalf * s_x[i]->z * s_x[i]->z /
				(1.0 + s_x[i]->dha * b * muhalf) + s_x[i]->dhb * mu;
			s_x[i]->dg = (c2 * s_x[i]->z * s_x[i]->z /
						  ((1.0 + s_x[i]->dha * b * muhalf) * (1.0 +
															   s_x[i]->dha *
															   b * muhalf)) +
						  s_x[i]->dhb) * LOG_10 * s_x[i]->moles;
/*			if (mu_x < 1e-6) s_x[i]->dg = 0.0; */
			break;
		case 3:				/* Always 1.0 */
			s_x[i]->lg = 0.0;
			s_x[i]->dg = 0.0;
			break;
		case 4:				/* Exchange */
/*
 *   Find CEC
 *   z contains valence of cation for exchange species, alk contains cec
 */
/* !!!!! */
			for (j = 1; s_x[i]->rxn_x->token[j].s != NULL; j++)
			{
				if (s_x[i]->rxn_x->token[j].s->type == EX)
				{
					s_x[i]->alk =
						s_x[i]->rxn_x->token[j].s->primary->unknown->moles;
					break;
				}
			}
			if (s_x[i]->exch_gflag == 1 && s_x[i]->alk > 0)
			{
				/* Davies */
				d1 = s_x[i]->lg;
				s_x[i]->lg = -s_x[i]->equiv * s_x[i]->equiv * a *
					(muhalf / (1.0 + muhalf) - 0.3 * mu) +
					log10(fabs(s_x[i]->equiv) / s_x[i]->alk);
				if (s_x[i]->a_f && s_x[i]->primary == NULL)
				{
					d2 = s_x[i]->moles * s_x[i]->equiv / s_x[i]->alk;
					if (d2 > 1) d2 = 1;
					d2 = s_x[i]->lg - s_x[i]->a_f * (1 - d2);
					d3 = 0.89;
					if (iterations < 10) d3 = 0.7; else d3 = 0.89;
					s_x[i]->lg = d3 * d1 + (1 - d3) * d2;
				}
				s_x[i]->dg =
					c1 * s_x[i]->equiv * s_x[i]->equiv * s_x[i]->moles;
			}
			else if (s_x[i]->exch_gflag == 2 && s_x[i]->alk > 0)
			{
				/* Extended D-H, WATEQ D-H */
				d1 = s_x[i]->lg;
				s_x[i]->lg = -a * muhalf * s_x[i]->equiv * s_x[i]->equiv /
					(1.0 + s_x[i]->dha * b * muhalf) + s_x[i]->dhb * mu +
					log10(fabs(s_x[i]->equiv) / s_x[i]->alk);
				if (s_x[i]->a_f && s_x[i]->primary == NULL)
				{
					d2 = s_x[i]->moles * s_x[i]->equiv / s_x[i]->alk;
					if (d2 > 1) d2 = 1;
					d2 = s_x[i]->lg - s_x[i]->a_f * (1 - d2);
					d3 = 0.89;
					if (iterations < 10) d3 = 0.7; else d3 = 0.89;
					s_x[i]->lg = d3 * d1 + (1 - d3) * d2;
				}
				s_x[i]->dg = (c2 * s_x[i]->equiv * s_x[i]->equiv /
							  ((1.0 + s_x[i]->dha * b * muhalf) * (1.0 +
																   s_x[i]->
																   dha * b *
																   muhalf)) +
							  s_x[i]->dhb) * LOG_10 * s_x[i]->moles;
			}
			else if (s_x[i]->exch_gflag == 7 && s_x[i]->alk > 0)
			{
				if (llnl_count_temp > 0)
				{
					s_x[i]->lg =
						-a_llnl * muhalf * s_x[i]->equiv * s_x[i]->equiv /
						(1.0 + s_x[i]->dha * b_llnl * muhalf) +
						bdot_llnl * mu +
						log10(fabs(s_x[i]->equiv) / s_x[i]->alk);
					s_x[i]->dg =
						(c2_llnl * s_x[i]->equiv * s_x[i]->equiv /
						 ((1.0 + s_x[i]->dha * b_llnl * muhalf) * (1.0 +
																   s_x[i]->
																   dha *
																   b_llnl *
																   muhalf)) +
						 bdot_llnl) * LOG_10 * s_x[i]->moles;
				}
				else
				{
					error_msg("LLNL_AQUEOUS_MODEL_PARAMETERS not defined.",
							  STOP);
				}
			}
			else
			{
/*
 *   Master species is a dummy variable with meaningless activity and mass
 */
				if (s_x[i]->primary != NULL)
				{
					s_x[i]->lg = 0.0;
					s_x[i]->dg = 0.0;
				}
				else
				{
					if (s_x[i]->alk <= 0)
					{
						s_x[i]->lg = 0.0;
					}
					else
					{
						s_x[i]->lg = log10(fabs(s_x[i]->equiv) / s_x[i]->alk);
					}
					s_x[i]->dg = 0.0;
				}
			}
			break;
		case 5:				/* Always 1.0 */
			s_x[i]->lg = 0.0;
			s_x[i]->dg = 0.0;
			break;
		case 6:				/* Surface */
/*
 *   Find moles of sites.
 *   s_x[i]->equiv is stoichiometric coefficient of sites in species
 */
			for (j = 1; s_x[i]->rxn_x->token[j].s != NULL; j++)
			{
				if (s_x[i]->rxn_x->token[j].s->type == SURF)
				{
					s_x[i]->alk =
						s_x[i]->rxn_x->token[j].s->primary->unknown->moles;
					break;
				}
			}
			if (use.Get_surface_ptr()->Get_type() == cxxSurface::CD_MUSIC)
			{
				/*  mole fraction */
				equiv = 1.0;
			}
			else
			{
				equiv = s_x[i]->equiv;
			}
			if (s_x[i]->alk > 0)
			{
				s_x[i]->lg = log10(equiv / s_x[i]->alk);
				s_x[i]->dg = 0.0;
			}
			else
			{
				s_x[i]->lg = 0.0;
				s_x[i]->dg = 0.0;
			}
			break;
		case 7:				/* LLNL */
			if (llnl_count_temp > 0)
			{
				if (s_x[i]->z == 0)
				{
					s_x[i]->lg = 0.0;
					s_x[i]->dg = 0.0;
				}
				else
				{
					s_x[i]->lg = -a_llnl * muhalf * s_x[i]->z * s_x[i]->z /
						(1.0 + s_x[i]->dha * b_llnl * muhalf) +
						bdot_llnl * mu;
					s_x[i]->dg =
						(c2_llnl * s_x[i]->z * s_x[i]->z /
						 ((1.0 + s_x[i]->dha * b_llnl * muhalf) * (1.0 +
																   s_x[i]->
																   dha *
																   b_llnl *
																   muhalf)) +
						 bdot_llnl) * LOG_10 * s_x[i]->moles;
					break;
				}
			}
			else
			{
				error_msg("LLNL_AQUEOUS_MODEL_PARAMETERS not defined.", STOP);
			}
			break;
		case 8:				/* LLNL CO2 */
			if (llnl_count_temp > 0)
			{
				s_x[i]->lg = log_g_co2;
				s_x[i]->dg = dln_g_co2 * s_x[i]->moles;
			}
			else
			{
				error_msg("LLNL_AQUEOUS_MODEL_PARAMETERS not defined.", STOP);
			}
			break;
		case 9:				/* activity water */
			s_x[i]->lg = log10(exp(s_h2o->la * LOG_10) * gfw_water);
			s_x[i]->dg = 0.0;
			break;
		}
/*
		if (mu_unknown != NULL) {
			if (fabs(residual[mu_unknown->number]) > 0.1 &&
				fabs(residual[mu_unknown->number])/mu_x > 0.5) {
				s_x[i]->dg = 0.0;
			}
		}
 */
	}
	return (OK);
}
/* ------------------------------------------------------------------------------- */
int Phreeqc::
ineq(int in_kode)
/* ------------------------------------------------------------------------------- */
{
/*
 *	Sets up equations and inequalities for Cl1.
 *	Scales columns if necessary.
 *	Eliminates equations that are not necessary because
 *		gas_phase, s_s, or phase equation is not needed
 *	Mallocs space
 *	Calls Cl1
 *	Rescales results if necessary
 */
	int i, j;
	int return_code;
	int l_count_rows;
	int l_count_optimize, count_equal;
	int k, l, m, n;
	int l_klmd, l_nklmd, l_n2d;
	int l_iter;
	LDBLE l_error;
	LDBLE max;
	int l_kode;
	LDBLE min;
#ifdef SLNQ
	LDBLE *slnq_array;
	LDBLE *slnq_delta1;
#endif
/*   Debug
	if (debug_model == TRUE) {
		output_msg(sformatf( "\narray:\n\n");
		array_print(array, count_unknowns, count_unknowns + 1, count_unknowns + 1));
	}
 */
/*
 *   Special case for removing unstable phases
 */
	cxxPPassemblage * pp_assemblage_ptr = use.Get_pp_assemblage_ptr();
	cxxPPassemblageComp * comp_ptr;

	if (remove_unstable_phases == TRUE)
	{
		if (debug_model == TRUE)
		{
			output_msg(sformatf(
					   "\nSolution vector for removing unstable phases:\n"));
		}
		for (i = 0; i < count_unknowns; i++)
		{

			if (x[i]->type == PP)
			{
				//std::map<std::string, cxxPPassemblageComp>::iterator it;
				//it =  pp_assemblage_ptr->Get_pp_assemblage_comps().find(x[i]->pp_assemblage_comp_name);
				//assert(it != pp_assemblage_ptr->Get_pp_assemblage_comps().end());
				cxxPPassemblageComp * comp_ptr = (cxxPPassemblageComp *) x[i]->pp_assemblage_comp_ptr;
				if (residual[i] > 0e-8 && x[i]->moles > 0 &&
					//it->second.Get_add_formula().size() == 0
					comp_ptr->Get_add_formula().size() == 0
					&& x[i]->dissolve_only == FALSE)
				{
					/*
					*   Set mass transfer to all of phase
					*/
					delta[i] = x[i]->moles;
				}
				else
				{
					delta[i] = 0.0;
				}
				if (debug_model == TRUE)
				{
					output_msg(sformatf( "%6d  %-12.12s %10.2e\n", i,
						x[i]->description, (double) delta[i]));
				}
			}
		}
		remove_unstable_phases = FALSE;
		return (OK);
	}
/*
 *   Pitzer model does not have activity of water //or mu
 */
	if (pitzer_model == TRUE || sit_model == TRUE)
	{
		for (i = 0; i < count_unknowns; i++)
		{
			if ((x[i]->type == AH2O && full_pitzer == FALSE) ||
				(x[i]->type == MH && pitzer_model == TRUE && pitzer_pe == FALSE) ||
				/*x[i]->type == MU ||*/
				(x[i]->type == PITZER_GAMMA && full_pitzer == FALSE))
			{
				for (j = 0; j < count_unknowns; j++)
				{
					array[j * (count_unknowns + 1) + i] = 0.0;
				}
				for (j = 0; j < count_unknowns + 1; j++)
				{
					array[i * (count_unknowns + 1) + j] = 0.0;
				}
			}
		}
	}

/*
 * Initialize space if necessary
 */
	ineq_init(3 * count_unknowns, 3 * count_unknowns);
/*
 *   Normalize column
 */
	space((void **) ((void *) &normal), count_unknowns, &normal_max,
		  sizeof(LDBLE));

	for (i = 0; i < count_unknowns; i++)
		normal[i] = 1.0;


	for (i = 0; i < count_unknowns; i++)
	{
		max = 0.0;

		if (x[i]->type == MB || x[i]->type == ALK || x[i]->type == EXCH
			|| x[i]->type == SURFACE || x[i]->type == SURFACE_CB
			|| x[i]->type == SURFACE_CB1 || x[i]->type == SURFACE_CB2)
		{
			if (x[i]->moles <= MIN_RELATED_SURFACE
			   && (x[i]->type == EXCH || x[i]->type == SURFACE))
				continue;
			for (j = 0; j < count_unknowns; j++)
			{
				if (x[i]->type == SURFACE && x[j]->type == SURFACE_CB)
					continue;
				if (x[i]->type == SURFACE_CB1 && x[j]->type == SURFACE_CB2)
					continue;

				if (fabs(array[j * (count_unknowns + 1) + i]) > max)
				{
					max = fabs(array[j * (count_unknowns + 1) + i]);
					if (max > min_value)
						break;
				}
			}
			if (diagonal_scale == TRUE)
			{
				if (fabs(array[i * (count_unknowns + 1) + i]) < min_value)
				{
					max = fabs(array[i * (count_unknowns + 1) + i]);
				}
			}

			if (max == 0)
			{
				array[i * (count_unknowns + 1) + i] = 1e-5 * x[i]->moles;
				max = fabs(1e-5 * x[i]->moles);
			}
		}

		if (x[i]->type == MH && (pitzer_model == FALSE || pitzer_pe == TRUE))
		{
			/* make absolute value of diagonal at least 1e-12 */

			min = 1e-12;
			min = MIN_TOTAL;
			array[x[i]->number * (count_unknowns + 1) + x[i]->number] += min;
			if (fabs
				(array[x[i]->number * (count_unknowns + 1) + x[i]->number]) <
				min)
				array[x[i]->number * (count_unknowns + 1) + x[i]->number] =
					min;
			max = 0.0;

			for (j = 0; j < count_unknowns; j++)
			{
				if (x[j]->type != MB &&
					x[j]->type != SURFACE &&
					x[j]->type != SURFACE_CB &&
					x[j]->type != SURFACE_CB1 &&
					x[j]->type != SURFACE_CB2 &&
					x[j]->type != EXCH && x[j]->type != MH
					&& x[j]->type != MH2O)
					continue;
				if (fabs(array[j * (count_unknowns + 1) + i]) > max)
				{
					max = fabs(array[j * (count_unknowns + 1) + i]);
					if (max > min_value)
						break;
				}
			}
		}

		if (max > 0.0 && max < min_value)
		{
			if (debug_model == TRUE)
			{
				output_msg(sformatf(
						   "Scaling column for %s, max= %e\n",
						   x[i]->description, (double) max));
			}
			for (j = 0; j < count_unknowns; j++)
			{
				array[j * (count_unknowns + 1) + i] *= min_value / max;
			}
			normal[i] = min_value / max;
		}
	}

/*
 *   Allocate arrays for inequality solver
 */
	max_row_count = 2 * count_unknowns + 2;
	max_column_count = count_unknowns + 2;
	space((void **) ((void *) &ineq_array), max_row_count * max_column_count,
		  &ineq_array_max, sizeof(LDBLE));

	space((void **) ((void *) &back_eq), max_row_count, &back_eq_max,
		  sizeof(int));

	space((void **) ((void *) &zero), max_row_count, &zero_max,
		  sizeof(LDBLE));
	zero_double(zero, max_row_count);

	space((void **) ((void *) &res), max_row_count, &res_max, sizeof(LDBLE));
	zero_double(res, max_row_count);

	space((void **) ((void *) &delta1), max_column_count, &delta1_max,
		  sizeof(LDBLE));
	zero_double(delta1, max_column_count);

/*
 *   Copy equations to optimize into ineq_array
 */
	l_count_rows = 0;
	for (i = 0; i < count_unknowns; i++)
	{
		if (iterations < aqueous_only)
			continue;
/*
 *   Pure phases
 */
		if (x[i]->type == PP)
		{
			//std::map<std::string, cxxPPassemblageComp>::iterator it;
			//it =  pp_assemblage_ptr->Get_pp_assemblage_comps().find(x[i]->pp_assemblage_comp_name);
			
			/* not in model, ignore */
			if (x[i]->phase->in == FALSE)
				continue;		
			// delay removing phase
			if (x[i]->moles > 0.0 || x[i]->f <= 0.0 || iterations == 0 || equi_delay == 0)
			{
				x[i]->iteration = iterations;
			}
			cxxPPassemblageComp * comp_ptr = (cxxPPassemblageComp *) x[i]->pp_assemblage_comp_ptr;
			//if (it->second.Get_force_equality())
			if (comp_ptr->Get_force_equality())
				continue;
			/*   Undersaturated and no mass, ignore */
			if (x[i]->f > 0e-8 && x[i]->moles <= 0
				&& iterations >= x[i]->iteration + equi_delay
				&& comp_ptr->Get_add_formula().size() == 0)
			{
				continue;
			}
			else if (x[i]->f < 0e-8 && x[i]->dissolve_only == TRUE
					 //&& (x[i]->moles - it->second.Get_initial_moles() >= 0))
					 && (x[i]->moles - comp_ptr->Get_initial_moles() >= 0))
			{
				continue;
			}
			else
			{
				/*   Copy in saturation index equation (has mass or supersaturated) */
				memcpy((void *) &(ineq_array[l_count_rows * max_column_count]),
					   (void *) &(array[i * (count_unknowns + 1)]),
					   (size_t) (count_unknowns + 1) * sizeof(LDBLE));
				back_eq[l_count_rows] = i;
				//if (it->second.Get_add_formula().size() == 0
				if (comp_ptr->Get_add_formula().size() == 0
					&& x[i]->dissolve_only == FALSE)
				{
					res[l_count_rows] = 1.0;
				}
/*
 *   If infeasible solution on first attempt, remove constraints on IAP
 */
				if (pp_scale != 1)
				{
					for (j = 0; j < count_unknowns + 1; j++)
					{
						ineq_array[l_count_rows * max_column_count + j] *=
							pp_scale;
					}
				}

				if (in_kode != 1)
				{
					res[l_count_rows] = 0.0;
				}
				l_count_rows++;
			}
		}
		else if (x[i]->type == ALK || x[i]->type == SOLUTION_PHASE_BOUNDARY)
		{
/*
 *   Alkalinity and solution phase boundary
 */
			memcpy((void *) &(ineq_array[l_count_rows * max_column_count]),
				   (void *) &(array[i * (count_unknowns + 1)]),
				   (size_t) (count_unknowns + 1) * sizeof(LDBLE));
			back_eq[l_count_rows] = i;
			l_count_rows++;
/*
 *   Gas phase
 */
		}
		else if (x[i]->type == GAS_MOLES && gas_in == TRUE)
		{
			memcpy((void *) &(ineq_array[l_count_rows * max_column_count]),
				   (void *) &(array[i * (count_unknowns + 1)]),
				   (size_t) (count_unknowns + 1) * sizeof(LDBLE));
			back_eq[l_count_rows] = i;

			res[l_count_rows] = 1.0;
			if (in_kode != 1)
			{
				res[l_count_rows] = 0.0;
			}
			l_count_rows++;
/*
 *   Solid solution
 */
		}
		else if (x[i]->type == SS_MOLES && x[i]->ss_in == TRUE)
		{
			memcpy((void *) &(ineq_array[l_count_rows * max_column_count]),
				   (void *) &(array[i * (count_unknowns + 1)]),
				   (size_t) (count_unknowns + 1) * sizeof(LDBLE));
			back_eq[l_count_rows] = i;
			res[l_count_rows] = 1.0;
			if (in_kode != 1)
			{
				res[l_count_rows] = 0.0;
			}
			l_count_rows++;
		}
	}
	l_count_optimize = l_count_rows;
/*
 *   Copy equality equations into ineq_array
 */
	for (i = 0; i < count_unknowns; i++)
	{
		comp_ptr = NULL;
		cxxPPassemblageComp *comp_ptr1 = NULL;
		if (x[i]->type == PP)
		{
			//comp_ptr = pp_assemblage_ptr->Find(x[i]->pp_assemblage_comp_name);
			comp_ptr = (cxxPPassemblageComp *) x[i]->pp_assemblage_comp_ptr;
		}
		if (x[i]->type == SURFACE && x[i]->phase_unknown != NULL)
		{
			comp_ptr = pp_assemblage_ptr->Find(x[i]->phase_unknown->phase->name);
		}
		if ((x[i]->type == SURFACE_CB || x[i]->type == SURFACE_CB1
				 || x[i]->type == SURFACE_CB2)
				&& x[i - 1]->phase_unknown != NULL)
		{
			comp_ptr1 = pp_assemblage_ptr->Find(x[i-1]->phase_unknown->phase->name);
		}
		if (x[i]->type != SOLUTION_PHASE_BOUNDARY &&
			x[i]->type != ALK &&
			x[i]->type != GAS_MOLES && x[i]->type != SS_MOLES &&
			x[i]->type != PITZER_GAMMA
			/* && x[i]->type != PP */
			)
		{
			if (x[i]->type == PP && !comp_ptr->Get_force_equality())
				continue;
			if (x[i]->type == MH && pitzer_model == TRUE && pitzer_pe == FALSE)
				continue;
			if (mass_water_switch == TRUE && x[i] == mass_oxygen_unknown)
				continue;
/*
 *   Mass balance, CB, MU, AH2O, MH, MH2O, others
 */
			if (x[i]->type == EXCH && x[i]->moles <= MIN_RELATED_SURFACE)
				continue;
			if (x[i]->type == SURFACE &&
				x[i]->phase_unknown == NULL
				&& x[i]->moles <= MIN_RELATED_SURFACE)
				continue;
			if ((x[i]->type == SURFACE_CB || x[i]->type == SURFACE_CB1
				 || x[i]->type == SURFACE_CB2) )
			{
				cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[i]->surface_charge);
				if (charge_ptr->Get_grams() <= MIN_RELATED_SURFACE)
				{
					continue;
				}
			}
			if (x[i]->type == SURFACE && x[i]->phase_unknown != NULL &&
				x[i]->phase_unknown->moles <= MIN_RELATED_SURFACE &&
				comp_ptr->Get_add_formula().size() == 0)
				continue;
			if ((x[i]->type == SURFACE_CB || x[i]->type == SURFACE_CB1
				|| x[i]->type == SURFACE_CB2)
				&& x[i - 1]->phase_unknown != NULL)
			{
				cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[i]->surface_charge);
				if (charge_ptr->Get_grams() <= MIN_RELATED_SURFACE &&
					comp_ptr1->Get_add_formula().size() == 0)
				{
					continue;
				}
			}
			memcpy((void *) &(ineq_array[l_count_rows * max_column_count]),
				   (void *) &(array[i * (count_unknowns + 1)]),
				   (size_t) (count_unknowns + 1) * sizeof(LDBLE));
			back_eq[l_count_rows] = i;
			if (mass_water_switch == TRUE && x[i] == mass_hydrogen_unknown)
			{
				k = mass_oxygen_unknown->number;
				for (j = 0; j < count_unknowns; j++)
				{
					ineq_array[l_count_rows * max_column_count + j] -=
						2 * array[k * (count_unknowns + 1) + j];
				}
			}
			l_count_rows++;
		}
		else if (x[i]->type == PITZER_GAMMA && full_pitzer == TRUE)
		{
			memcpy((void *) &(ineq_array[l_count_rows * max_column_count]),
				   (void *) &(array[i * (count_unknowns + 1)]),
				   (size_t) (count_unknowns + 1) * sizeof(LDBLE));
			back_eq[l_count_rows] = i;
			l_count_rows++;
		}
	}
	count_equal = l_count_rows - l_count_optimize;
/*
 *   Copy inequality constraints into ineq
 */
	if (pure_phase_unknown != NULL)
	{
		for (i = 0; i < count_unknowns; i++)
		{
			if (x[i]->type == PP)
			{
				//comp_ptr = pp_assemblage_ptr->Find(x[i]->pp_assemblage_comp_name);
				comp_ptr = (cxxPPassemblageComp *) x[i]->pp_assemblage_comp_ptr;
				/* not in model, ignore */
				if (x[i]->phase->in == FALSE)
					continue;
				/*   No moles and undersaturated, ignore */
				if (x[i]->moles <= 0.0 && x[i]->f > 0e-8 &&
					comp_ptr->Get_add_formula().size() == 0)
				{
					continue;
					/*   No moles of pure phase present, must precipitate */
				}
				else if (x[i]->moles <= 0.0)
				{
					delta1[i] = -1.0;
				}
				else if (x[i]->f < 0e-8 && x[i]->dissolve_only == TRUE
						 && (x[i]->moles - comp_ptr->Get_initial_moles() >=
							 0))
				{
					continue;
				}
				else
				{

					/*   Pure phase is present, force Mass transfer to be <= amount of mineral remaining */
					memcpy((void *)
						   &(ineq_array[l_count_rows * max_column_count]),
						   (void *) &(zero[0]),
						   (size_t) (count_unknowns + 1) * sizeof(LDBLE));
					ineq_array[l_count_rows * max_column_count + i] = 1.0;
					ineq_array[l_count_rows * max_column_count +
							   count_unknowns] = x[i]->moles;
					back_eq[l_count_rows] = i;
					l_count_rows++;
				}
				/*   Pure phase is present and dissolve_only, force ppt to be <= amount of dissolved so far */
				if (x[i]->dissolve_only == TRUE)
				{
					memcpy((void *)
						   &(ineq_array[l_count_rows * max_column_count]),
						   (void *) &(zero[0]),
						   (size_t) (count_unknowns + 1) * sizeof(LDBLE));
					ineq_array[l_count_rows * max_column_count + i] = -1.0;
					ineq_array[l_count_rows * max_column_count +
							   count_unknowns] =
						comp_ptr->Get_initial_moles() - x[i]->moles;
					back_eq[l_count_rows] = i;
					l_count_rows++;
				}
			}
		}
	}
/*
 *   Add inequality for mass of oxygen greater than zero
 */
	if (pitzer_model || sit_model)
	{
		for (i = 0; i < count_unknowns; i++)
		{
			if (x[i]->type == MH2O)
			{
				memcpy((void *) &(ineq_array[l_count_rows * max_column_count]),
					   (void *) &(array[i * (count_unknowns + 1)]),
					   (size_t) (count_unknowns + 1) * sizeof(LDBLE));
				back_eq[l_count_rows] = i;
				for (j = 0; j < count_unknowns; j++)
				{
					if (x[j]->type < PP)
					{
						ineq_array[l_count_rows * max_column_count + j] = 0.0;
					}
					else
					{
						/*ineq_array[l_count_rows*max_column_count + j] = -ineq_array[l_count_rows*max_column_count + j]; */
					}
				}
				ineq_array[l_count_rows * max_column_count + count_unknowns] =
					0.5 * x[i]->moles;
				l_count_rows++;
			}
		}
	}


/*
 *   Hydrogen mass balance is good
 */
/*
 *   No moles and undersaturated, mass transfer must be zero
 */
	if (pure_phase_unknown != NULL)
	{
		for (i = 0; i < count_unknowns; i++)
		{
			if (x[i]->type == PP)
			{
				//comp_ptr = pp_assemblage_ptr->Find(x[i]->pp_assemblage_comp_name);			
				comp_ptr = (cxxPPassemblageComp *) x[i]->pp_assemblage_comp_ptr;
				if ((x[i]->moles <= 0.0 && x[i]->f > 0e-8 &&
					 comp_ptr->Get_add_formula().size() == 0)
					|| x[i]->phase->in == FALSE)
				{
					for (j = 0; j < l_count_rows; j++)
					{
						ineq_array[j * max_column_count + i] = 0.0;
					}
				}
				if (x[i]->dissolve_only == TRUE)
				{
					if (x[i]->f < 0e-8
						&& (x[i]->moles - comp_ptr->Get_initial_moles() >=
							0))
					{
						for (j = 0; j < l_count_rows; j++)
						{
							ineq_array[j * max_column_count + i] = 0.0;
						}
					}
				}
			}
		}
	}
/*
 *   No moles of exchanger
 */
	if (use.Get_exchange_ptr() != NULL
		&& (use.Get_exchange_ptr()->Get_related_phases() ||
		use.Get_exchange_ptr()->Get_related_rate()))
	{
		for (i = 0; i < count_unknowns; i++)
		{
			if (x[i]->type == EXCH && x[i]->moles <= 0)
			{
				for (j = 0; j < l_count_rows; j++)
				{
					ineq_array[j * max_column_count + i] = 0.0;
				}
			}
		}
	}
/*
 *   No moles of surface
 */
	if (use.Get_surface_ptr() != NULL
		&& (use.Get_surface_ptr()->Get_related_phases()
			|| use.Get_surface_ptr()->Get_related_rate()))
	{
		for (i = 0; i < count_unknowns; i++)
		{
			comp_ptr = NULL;
			cxxPPassemblageComp * comp_ptr1 = NULL;
			if (x[i]->type == SURFACE && x[i]->phase_unknown != NULL)
			{
				comp_ptr = pp_assemblage_ptr->Find(x[i]->phase_unknown->phase->name);
			}
			if ((x[i]->type == SURFACE_CB || x[i]->type == SURFACE_CB1
				  || x[i]->type == SURFACE_CB2)
				 && (x[i - 1]->phase_unknown != NULL))
			{
				comp_ptr1 = pp_assemblage_ptr->Find(x[i-1]->phase_unknown->phase->name);
			}
			LDBLE grams = 0;
			if ((x[i]->type == SURFACE_CB || x[i]->type == SURFACE_CB1
				  || x[i]->type == SURFACE_CB2)
				 && x[i - 1]->phase_unknown != NULL)
			{
				cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[i]->surface_charge);
				grams = charge_ptr->Get_grams();
			}
			if ((x[i]->type == SURFACE && x[i]->phase_unknown != NULL &&
				 x[i]->phase_unknown->moles <= MIN_RELATED_SURFACE &&
				 comp_ptr->Get_add_formula().size() == 0) ||
				((x[i]->type == SURFACE_CB || x[i]->type == SURFACE_CB1
				  || x[i]->type == SURFACE_CB2)
				 && x[i - 1]->phase_unknown != NULL &&
				 grams <= MIN_RELATED_SURFACE &&
				 comp_ptr1->Get_add_formula().size() == 0))
			{
				for (j = 0; j < l_count_rows; j++)
				{
					ineq_array[j * max_column_count + i] = 0.0;
				}
			}
		}
	}
/*
 *   No moles of surface
 */
	if (use.Get_surface_ptr() != NULL)
	{
		LDBLE grams = 0;;
		for (i = 0; i < count_unknowns; i++)
		{
			if ((x[i]->type == SURFACE_CB || x[i]->type == SURFACE_CB1
				  || x[i]->type == SURFACE_CB2) && x[i - 1]->phase_unknown == NULL)
			{
				cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[i]->surface_charge);
				grams = charge_ptr->Get_grams();
			}
			if ((x[i]->type == SURFACE &&
				 x[i]->phase_unknown == NULL &&
				 x[i]->moles <= MIN_RELATED_SURFACE) ||
				((x[i]->type == SURFACE_CB || x[i]->type == SURFACE_CB1
				  || x[i]->type == SURFACE_CB2)
				 && x[i - 1]->phase_unknown == NULL
				 && grams <= MIN_RELATED_SURFACE))
			{
				for (j = 0; j < l_count_rows; j++)
				{
					ineq_array[j * max_column_count + i] = 0.0;
				}
			}
		}
	}
/*
 *   Moles of gas must be >= zero
 */
	if (gas_in == TRUE)
	{
		for (i = gas_unknown->number; i < count_unknowns; i++)
		{
			if (x[i]->type == GAS_MOLES)
			{
				memcpy((void *) &(ineq_array[l_count_rows * max_column_count]),
					   (void *) &(zero[0]),
					   (size_t) (count_unknowns + 1) * sizeof(LDBLE));
				ineq_array[l_count_rows * max_column_count + i] = -1.0;
				ineq_array[l_count_rows * max_column_count + count_unknowns] =
					x[i]->moles;
				back_eq[l_count_rows] = i;
				l_count_rows++;
			}
			else
			{
				break;
			}
		}
	}
	else if (use.Get_gas_phase_ptr() != NULL && gas_in == FALSE)
	{
/*
 *   Moles of gas small and sum p < ptotal
 */
		i = gas_unknown->number;
		for (j = 0; j < l_count_rows; j++)
		{
			ineq_array[j * max_column_count + i] = 0.0;
		}
	}
/*
 *   Phase must be "in" and moles of solid solution must be >= zero
 */

	if (ss_unknown != NULL)
	{
		for (i = ss_unknown->number; i < count_unknowns; i++)
		{
			if (x[i]->type != SS_MOLES)
				break;
			if (x[i]->phase->in == TRUE && x[i]->ss_in == TRUE)
			{
				memcpy((void *) &(ineq_array[l_count_rows * max_column_count]),
					   (void *) &(zero[0]),
					   (size_t) (count_unknowns + 1) * sizeof(LDBLE));
				ineq_array[l_count_rows * max_column_count + i] = 1.0;
				ineq_array[l_count_rows * max_column_count + count_unknowns] =
					0.99 * x[i]->moles - MIN_TOTAL_SS;
				back_eq[l_count_rows] = i;
				l_count_rows++;
			}
			else
			{
				for (j = 0; j < l_count_rows; j++)
				{
					ineq_array[j * max_column_count + i] = 0.0;
				}
			}
		}
	}
/*
 *   Add inequality if total moles of element is less than zero
 */
	if (negative_concentrations)
	{
		for (i = 0; i < count_unknowns; i++)
		{
			if (x[i]->type == MB && x[i]->moles < 0.0)
			{
				memcpy((void *) &(ineq_array[l_count_rows * max_column_count]),
					   (void *) &(array[i * (count_unknowns + 1)]),
					   (size_t) (count_unknowns + 1) * sizeof(LDBLE));
				back_eq[l_count_rows] = i;
				for (j = 0; j < count_unknowns; j++)
				{
					if (x[j]->type < PP)
					{
						ineq_array[l_count_rows * max_column_count + j] = 0.0;
					}
				}
				l_count_rows++;
			}
		}
	}
/*
 *   Zero column for mass of water
 */
	if (mass_oxygen_unknown != NULL && mass_water_switch == TRUE)
	{
		k = mass_oxygen_unknown->number;
		for (j = 0; j < l_count_rows + 1; j++)
		{
			ineq_array[j * max_column_count + k] = 0;
		}
	}
/*
 *   Scale column for pure phases
 */
	for (i = 0; i < count_unknowns; i++)
	{
		if ((x[i]->type == PP || x[i]->type == SS_MOLES)
			&& pp_column_scale != 1.0)
		{
			for (j = 0; j < l_count_rows; j++)
			{
				ineq_array[j * max_column_count + i] *= pp_column_scale;
			}
			normal[i] = pp_column_scale;
		}

	}
	if (debug_model == TRUE)
	{
		output_msg(sformatf( "\nA and B arrays:\n\n"));
		array_print(ineq_array, l_count_rows, count_unknowns + 1,
					max_column_count);
	}
/*
 *   Calculate dimensions
 */
	k = l_count_optimize;			/* rows in A */
	l = count_equal;			/* rows in C */
	m = l_count_rows - l - k;		/* rows in E */
	if (m < 0)
		m = 0;

	if (debug_model == TRUE)
	{
		output_msg(sformatf( "k, l, m\t%d\t%d\t%d\n", k, l, m));
	}
#define SHRINK_ARRAY
#ifdef SHRINK_ARRAY
	if (sit_model && full_pitzer == FALSE)
	{
		n = count_unknowns - (int) s_list.size();
		for (int i = 0; i < l_count_rows; i++)
		{
			//for (int j = 0; j < n; j++)
			//{
			//	ineq_array[i*(n+2) + j] = ineq_array[i*(count_unknowns+2) +j];
			//}
			if (i > 0)
			{
				memcpy((void *) &ineq_array[i*(n+2)], (void *) &ineq_array[i*(count_unknowns+2)], (size_t) (n) * sizeof(LDBLE));
			}
			ineq_array[i*(n+2) + n] = ineq_array[i*(count_unknowns+2) + count_unknowns];
		}
	}
	else
	{
		n = count_unknowns;			/* columns in A, C, E */
	}
#else
	n = count_unknowns;			/* columns in A, C, E */
#endif
	l_klmd = max_row_count - 2;
	l_nklmd = n + l_klmd;
	l_n2d = n + 2;
/*
 *   Retain constraints on mineral mass transfers, even if infeasible on
 *   first attempt.
 */
	l_kode = 1;

	if (in_kode == 2)
	{
		l_kode = 1;
	}
	l_iter = 2*(n + l_count_rows);
/*
 *   Allocate space for arrays
 */
	space((void **) ((void *) &cu), 2 * l_nklmd, &cu_max, sizeof(LDBLE));

	space((void **) ((void *) &iu), 2 * l_nklmd, &iu_max, sizeof(int));

	space((void **) ((void *) &is), l_klmd, &is_max, sizeof(int));

#ifdef SLNQ
	slnq_array =
		(LDBLE *) PHRQ_malloc((size_t) count_unknowns *
							  (count_unknowns + 1) * sizeof(LDBLE));
	if (slnq_array == NULL)
		malloc_error();
	for (i = 0; i < k + l; i++)
	{
		for (j = 0; j <= count_unknowns; j++)
		{
			slnq_array[i * (count_unknowns + 1) + j] =
				ineq_array[i * max_column_count + j];
		}
	}
	slnq_delta1 =
		(LDBLE *) PHRQ_malloc((size_t) max_column_count * sizeof(LDBLE));
	if (slnq_delta1 == NULL)
		malloc_error();
	memcpy((void *) &(slnq_delta1[0]), (void *) &(zero[0]),
		   (size_t) max_column_count * sizeof(LDBLE));
#endif
/*
 *   Call CL1
 */
	cl1(k, l, m, n,
		l_nklmd, l_n2d, ineq_array,
		&l_kode, ineq_tol, &l_iter, delta1, res, &l_error, cu, iu, is, FALSE);
/*   Set return_kode */
	if (l_kode == 1)
	{
		return_code = ERROR;
	}
	else if (l_kode == 2)
	{
		return_code = 2;
	}
	else if (l_kode == 3)
	{
		return_code = ERROR;
		return_code = 2;
		warning_msg("Too many iterations in Cl1. Should not have done this.");
	}
	else
	{
		return_code = OK;
	}
#ifdef SLNQ
/*	if (l_kode > 0 && ((k + l) == count_unknowns)) { */
	if (l_kode > 0 && ((k + l) <= count_unknowns))
	{
		if (add_trivial_eqns(k + l, count_unknowns, slnq_array) == TRUE)
		{
			if (debug_model == TRUE)
				output_msg(sformatf( "Calling SLNQ, iteration %d\n",
						   iterations));
			log_msg(sformatf( "Calling SLNQ, iteration %d\n",
					   iterations));
			if (slnq
				(count_unknowns, slnq_array, slnq_delta1, count_unknowns + 1,
				 debug_model) == OK)
			{
				memcpy((void *) &(delta1[0]), (void *) &(slnq_delta1[0]),
					   (size_t) count_unknowns * sizeof(LDBLE));
				if (debug_model == TRUE)
					output_msg(sformatf( "Using SLNQ results.\n"));
				log_msg(sformatf( "Using SLNQ results.\n"));
				return_code = OK;
			}
			else
			{
				if (debug_model == TRUE)
					output_msg(sformatf(
							   "Could not use SLNQ results.\n"));
				log_msg(sformatf( "Could not use SLNQ results.\n"));
			}
		}
		else
		{
			log_msg(sformatf(
					   "Could not call SLNQ, row %d, unknowns %d, iteration %d\n",
					   k + l, count_unknowns, iterations));
		}
	}
	else if (l_kode > 0)
	{
		log_msg(sformatf( "Could not call SLNQ, row %d, unknowns %d\n",
				   k + l, count_unknowns));
	}
#endif
/*   Copy delta1 into delta and scale */
#ifdef SHRINK_ARRAY	
	memcpy((void *) &(delta[0]), (void *) &(zero[0]),
		   (size_t) count_unknowns * sizeof(LDBLE));
#endif
	memcpy((void *) &(delta[0]), (void *) &(delta1[0]),
		   (size_t) n * sizeof(LDBLE));
	for (i = 0; i < n; i++)
		delta[i] *= normal[i];
/*
 *   Rescale columns of array
 */
	for (i = 0; i < count_unknowns; i++)
	{
		if (normal[i] != 1.0)
		{
			for (j = 0; j < count_unknowns; j++)
			{
				array[j * (count_unknowns + 1) + i] /= normal[i];
			}
		}
	}

	/* adjust unknowns where cl1 failed */
	if (l_kode != 0)
	{
		for (i = 0; i < l_count_rows; i++)
		{
			/*
			output_msg(sformatf( "%6d  %-12.12s %10.2e\n", i,
					   x[back_eq[i]]->description, (double) res[i);
			*/
			j = back_eq[i];
			if (x[j]->type == MB && delta[j] == 0.0 && fabs(res[i]) > ineq_tol)
			{
				delta[j] = res[i]/fabs(res[i]) * 1;
			}
		}
	}

/*
 *   Debug, write results of ineq
 */

	if (debug_model == TRUE)
	{
		output_msg(sformatf( "kode: %d\titer: %d\terror: %e\n", l_kode,
				   l_iter, (double) l_error));
		output_msg(sformatf( "\nsolution vector:\n"));
		for (i = 0; i < count_unknowns; i++)
		{
			output_msg(sformatf( "%6d  %-12.12s %10.2e", i,
					   x[i]->description, (double) delta[i]));
			if (x[i]->type == PP)
			{
				output_msg(sformatf( "   -SI %10.2e   Moles %10.2e",
						   (double) x[i]->f, (double) x[i]->moles));
				if (x[i]->f < 0e-8 || x[i]->moles > 0.0)
				{
					output_msg(sformatf( " **"));
				}
			}
			output_msg(sformatf( "\n"));
		}

		output_msg(sformatf( "\nresidual vector:\n"));
		for (i = 0; i < l_count_rows; i++)
		{
			output_msg(sformatf( "%6d  %-12.12s %10.2e\n", i,
					   x[back_eq[i]]->description, (double) res[i]));
		}
	}
#ifdef SLNQ
	slnq_array = (LDBL *) free_check_null(slnq_array);
	slnq_delta1 = (LDBL *) free_check_null(slnq_delta1);
#endif
	return (return_code);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
jacobian_sums(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Fills in jacobian array, uses arrays sum_jacob0, sum_jacob1, and
 *   sum_jacob2.
 */
	int i, j, k;
	LDBLE sinh_constant;
/*
 *   Clear array, note residuals are in array[i, count_unknowns+1]
 */
	for (i = 0; i < count_unknowns; i++)
	{
		array[i] = 0.0;
	}
	for (i = 1; i < count_unknowns; i++)
	{
		memcpy((void *) &(array[i * (count_unknowns + 1)]),
			   (void *) &(array[0]), (size_t) count_unknowns * sizeof(LDBLE));
	}
/*
 *   Add constant terms
 */
	for (k = 0; k < count_sum_jacob0; k++)
	{
		*sum_jacob0[k].target += sum_jacob0[k].coef;
	}
/*
 *   Add terms with coefficients of 1.0
 */
	for (k = 0; k < count_sum_jacob1; k++)
	{
		*sum_jacob1[k].target += *sum_jacob1[k].source;
	}
/*
 *   Add terms with coefficients != 1.0
 */
	for (k = 0; k < count_sum_jacob2; k++)
	{
		*sum_jacob2[k].target += *sum_jacob2[k].source * sum_jacob2[k].coef;
	}
/*
 *   Make final adustments to jacobian array
 */
/*
 *   Ionic strength
 */
	if (mu_unknown != NULL)
	{
		for (i = 0; i < count_unknowns; i++)
		{
			// using straight mu equation
			array[mu_unknown->number * (count_unknowns + 1) + i] *= 0.5;
		}
		array[mu_unknown->number * (count_unknowns + 1) +
			  mu_unknown->number] -= mass_water_aq_x;
	}
/*
 *   Mass of oxygen
 */
	if (mass_oxygen_unknown != NULL && mu_unknown != NULL)
	{
		array[mu_unknown->number * (count_unknowns + 1) +
			  mass_oxygen_unknown->number] -= mu_x * mass_water_aq_x;
	}
/*
 *   Activity of water
 */
	if (ah2o_unknown != NULL)
	{
		if (dampen_ah2o)
		{
			// factors for tanh attenuation to make ah2o positive
			LDBLE y_sum = ah2o_unknown->f;
			LDBLE x_h2o = mass_water_aq_x;
			LDBLE a = AH2O_FACTOR;
			LDBLE lim = 95.0 - 100.0*a*y_sum/x_h2o;
			LDBLE factor = -a*(x_h2o*(-0.5*tanh(lim) - 0.5) + (50.0*a*y_sum - 47.5 * x_h2o) / (cosh(lim)*cosh(lim))) /
				(x_h2o*x_h2o);

			for (i = 0; i < count_unknowns; i++)
			{
				array[ah2o_unknown->number * (count_unknowns + 1) + i] *= factor;
			}
			// activity of water term
			array[ah2o_unknown->number * (count_unknowns + 1) +
				  ah2o_unknown->number] -= exp(s_h2o->la * LOG_10);

			// mass of water term
			if (mass_oxygen_unknown != NULL)
			{
				array[ah2o_unknown->number * (count_unknowns + 1) + mass_oxygen_unknown->number] -=
					  a*y_sum*(x_h2o*(0.5*tanh(lim) + 0.5) + (47.5*x_h2o - 50.0*a*y_sum)/(cosh(lim)*cosh(lim))) / 
					  (x_h2o*x_h2o*x_h2o);
			}
		}
		else
		{
			for (i = 0; i < count_unknowns; i++)
			{
				array[ah2o_unknown->number * (count_unknowns + 1) + i] *= -AH2O_FACTOR;
			}
			array[ah2o_unknown->number * (count_unknowns + 1) + ah2o_unknown->number] -=
				mass_water_aq_x * exp(s_h2o->la * LOG_10);
			if (mass_oxygen_unknown != NULL)
			{
				array[ah2o_unknown->number * (count_unknowns + 1) + mass_oxygen_unknown->number] -=
					(exp(s_h2o->la * LOG_10) - 1) * mass_water_aq_x;
			}
		}
	}

/*
 *   Surface charge balance
 */
	if (surface_unknown != NULL && dl_type_x == cxxSurface::NO_DL)
	{
		if (use.Get_surface_ptr()->Get_type() != cxxSurface::CCM)
		{
			sinh_constant =
				//sqrt(8 * EPSILON * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000) * tk_x *
				//	 1000);
				sqrt(8 * eps_r * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000) * tk_x *
				1000);
			for (i = 0; i < count_unknowns; i++)
			{
				cxxSurfaceCharge *charge_ptr = NULL;
				if (x[i]->type == SURFACE_CB)
				{
					charge_ptr = use.Get_surface_ptr()->Find_charge(x[i]->surface_charge);
				}
				if (x[i]->type == SURFACE_CB && charge_ptr->Get_grams() > 0)
				{
					for (j = 0; j < count_unknowns; j++)
					{
						array[x[i]->number * (count_unknowns + 1) + j] *=
							F_C_MOL / (charge_ptr->Get_specific_area() *
							charge_ptr->Get_grams());
					}
					array[x[i]->number * (count_unknowns + 1) + x[i]->number] -=
						sinh_constant * sqrt(mu_x) *
						cosh(x[i]->master[0]->s->la * LOG_10);
					if (mu_unknown != NULL)
					{
						array[x[i]->number * (count_unknowns + 1) +
							mu_unknown->number] -=
							0.5 * sinh_constant / sqrt(mu_x) *
							sinh(x[i]->master[0]->s->la * LOG_10);
					}
				}
			}
		}
		else
		{
			for (i = 0; i < count_unknowns; i++)
			{
				cxxSurfaceCharge *charge_ptr = NULL;
				if (x[i]->type == SURFACE_CB)
				{
					charge_ptr = use.Get_surface_ptr()->Find_charge(x[i]->surface_charge);
				}
				if (x[i]->type == SURFACE_CB && charge_ptr->Get_grams() > 0)
				{
					for (j = 0; j < count_unknowns; j++)
					{
						array[x[i]->number * (count_unknowns + 1) + j] *=
							F_C_MOL / (charge_ptr->Get_specific_area() *
							charge_ptr->Get_grams());
					}
					array[x[i]->number * (count_unknowns + 1) + x[i]->number] -=
						charge_ptr->Get_capacitance0() * 2 * R_KJ_DEG_MOL *
								 tk_x * LOG_10 / F_KJ_V_EQ;
				}
			}
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
mb_sums(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Calculates sums of species for calculation of mass balances, charge
 *   balance. Also calculates saturation indices for solution_phase_boundaries
 *   and pure_phases. After this routine total calcium calculated from all
 *   calcium species in solution is stored in x[i]->f.  Also calculates
 *   x[i]->sum for some types of unknowns. Uses arrays sum_mb1 and
 *   sum_mb1, which are generated in prep and reprep.
 */
	int k;
/*
 *   Clear functions in unknowns
 */
	for (k = 0; k < count_unknowns; k++)
	{
		x[k]->f = 0.0;
		x[k]->sum = 0.0;
	}
/*
 *   Add terms with coefficients of 1.0
 */
	for (k = 0; k < count_sum_mb1; k++)
	{
		*sum_mb1[k].target += *sum_mb1[k].source;
/*		{ k += 1; k -= 1;} */
	}
/*
 *   Add terms with coefficients != 1.0
 */
	for (k = 0; k < count_sum_mb2; k++)
	{
		*sum_mb2[k].target += *sum_mb2[k].source * sum_mb2[k].coef;
/*		{ k += 1; k -= 1;} */
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
mb_gases(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Determines whether gas_phase equation is needed
 */
	gas_in = FALSE;
	if (gas_unknown == NULL || use.Get_gas_phase_ptr() == NULL)
		return (OK);
	cxxGasPhase *gas_phase_ptr = use.Get_gas_phase_ptr();
	if (gas_phase_ptr->Get_type() == cxxGasPhase::GP_PRESSURE)
	{
		if (gas_unknown->f > gas_phase_ptr->Get_total_p() + 1e-7 ||
			gas_unknown->moles > MIN_TOTAL)
		{
			gas_in = TRUE;
			//patm_x = gas_phase_ptr->Get_total_p();

		}
	}
	else
	{
		if (numerical_fixed_volume && (gas_phase_ptr->Get_pr_in() || force_numerical_fixed_volume))
		{
			gas_in = TRUE;
		}
		else
		{
			gas_in = FALSE;
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
mb_ss(void)
/* ---------------------------------------------------------------------- */
{
	LDBLE lp, log10_iap, total_moles;
	LDBLE iapc, iapb, l_kc, l_kb, lc, lb, xcaq, xbaq, xb, xc;
	LDBLE sigmapi_aq, sigmapi_solid;
	LDBLE total_p;
	struct rxn_token *rxn_ptr;
/*
 *   Determines whether solid solution equation is needed
 */
	if (ss_unknown == NULL || use.Get_ss_assemblage_ptr() == NULL)
		return (OK);
	std::vector<cxxSS *> ss_ptrs = use.Get_ss_assemblage_ptr()->Vectorize();
	for (size_t i = 0; i < ss_ptrs.size(); i++)
	{
		cxxSS *ss_ptr = ss_ptrs[i];
		total_moles = 0;
		for (size_t j = 0; j < ss_ptr->Get_ss_comps().size(); j++)
		{
			cxxSScomp *comp_ptr = &(ss_ptr->Get_ss_comps()[j]);
			total_moles += comp_ptr->Get_moles();
		}
		if (total_moles > 1e-13)
		{
			ss_ptr->Set_ss_in(true);
		}
		else if (ss_ptr->Get_a0() != 0.0 || ss_ptr->Get_a1() != 0.0)
		{
			int l;
			struct phase *phase0_ptr = phase_bsearch(ss_ptr->Get_ss_comps()[0].Get_name().c_str(), &l, FALSE);
			struct phase *phase1_ptr = phase_bsearch(ss_ptr->Get_ss_comps()[1].Get_name().c_str(), &l, FALSE);
			/*
			 *  Calculate IAPc and IAPb
			 */
			if (phase0_ptr->rxn_x != NULL)
			{
				log10_iap = 0;
				for (rxn_ptr = phase0_ptr->rxn_x->token + 1;
					 rxn_ptr->s != NULL; rxn_ptr++)
				{
					log10_iap += rxn_ptr->s->la * rxn_ptr->coef;
				}
				iapc = exp(log10_iap * LOG_10);
			}
			else
			{
				iapc = 1e-99;
			}
			if (phase1_ptr->rxn_x != NULL)
			{
				log10_iap = 0;
				for (rxn_ptr = phase1_ptr->rxn_x->token + 1;
					 rxn_ptr->s != NULL; rxn_ptr++)
				{
					log10_iap += rxn_ptr->s->la * rxn_ptr->coef;
				}
				iapb = exp(log10_iap * LOG_10);
			}
			else
			{
				iapb = 1e-99;
			}
			/*
			 *  Calculate sigma pi, aq
			 */
			sigmapi_aq = iapc + iapb;
			/*
			 *  Calculate xc,aq and xb, aq
			 */
			xcaq = iapc / (iapb + iapc);
			xbaq = iapb / (iapb + iapc);
			/*
			 *  Get Kc and Kb
			 */
			l_kc = exp(phase0_ptr->lk * LOG_10);
			l_kb = exp(phase1_ptr->lk * LOG_10);
			/*
			 *  Solve for xb
			 */
			xb = ss_root(ss_ptr->Get_a0(), ss_ptr->Get_a1(), l_kc, l_kb, xcaq, xbaq);
			/*
			 *  Calculate lambdac and lambdab
			 */
			xc = 1 - xb;
			lc = exp((ss_ptr->Get_a0() - ss_ptr->Get_a1() * (-4 * xb + 3)) * xb * xb);
			lb = exp((ss_ptr->Get_a0() + ss_ptr->Get_a1() * (4 * xb - 1)) * xc * xc);
			/*
			 *  Calculate sigma pi, solid
			 */
			sigmapi_solid = xb * lb * l_kb + xc * lc * l_kc;
			/*
			 * If Sigma pi, solid < sigma pi, aq, then use eqns
			 */
			if (sigmapi_solid < sigmapi_aq)
			{
				ss_ptr->Set_ss_in(true);
			}
			else
			{
				ss_ptr->Set_ss_in(false);
			}
		}
		else
		{
			/*
			 *  Calculate total mole fraction from solution activities
			 */
			total_p = 0;
			for (size_t j = 0; j < ss_ptr->Get_ss_comps().size(); j++)
			{
				cxxSScomp *comp_ptr = &(ss_ptr->Get_ss_comps()[j]);
				int l;
				struct phase *phase_ptr = phase_bsearch(comp_ptr->Get_name().c_str(), &l, FALSE);
				if (phase_ptr->in == TRUE)
				{
					lp = -phase_ptr->lk;
					for (rxn_ptr = phase_ptr->rxn_x->token + 1;
						 rxn_ptr->s != NULL; rxn_ptr++)
					{
						lp += rxn_ptr->s->la * rxn_ptr->coef;
					}
					total_p += exp(lp * LOG_10);
				}
			}
			if (total_p > 1.0)
			{
				ss_ptr->Set_ss_in(true);
			}
			else
			{
				ss_ptr->Set_ss_in(false);
			}
		}
	}
	for (int i = ss_unknown->number; i < count_unknowns; i++)
	{
		if (x[i]->type != SS_MOLES)
			break;
		//cxxSS *ss_ptr = use.Get_ss_assemblage_ptr()->Find(x[i]->ss_name);
		cxxSS *ss_ptr = (cxxSS *) x[i]->ss_ptr;
		x[i]->ss_in = ss_ptr->Get_ss_in() ? TRUE : FALSE;
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
molalities(int allow_overflow)
/* ---------------------------------------------------------------------- */
{
/*
 *   Calculates la for master species
 *   Calculates lm and moles from lk, lg, and la's of master species
 *   Adjusts lm of h2 and o2.
 */
	int i, j;
	LDBLE total_g;
	struct rxn_token *rxn_ptr;
/*
 *   la for master species
 */
	for (i = 0; i < count_master; i++)
	{
		if (master[i]->in == REWRITE)
		{
			master[i]->s->la = master[i]->s->lm + master[i]->s->lg;
		}
	}
	if (dl_type_x != cxxSurface::NO_DL)
	{
		s_h2o->tot_g_moles = s_h2o->moles;
		s_h2o->tot_dh2o_moles = 0.0;
	}
	for (i = 0; i < count_s_x; i++)
	{
		if (s_x[i]->type > HPLUS && s_x[i]->type != EX
			&& s_x[i]->type != SURF)
			continue;
/*
 *   lm and moles for all aqueous species
 */
		s_x[i]->lm = s_x[i]->lk - s_x[i]->lg;
		for (rxn_ptr = s_x[i]->rxn_x->token + 1; rxn_ptr->s != NULL;
			 rxn_ptr++)
		{
			s_x[i]->lm += rxn_ptr->s->la * rxn_ptr->coef;
			/*
			if (isnan(rxn_ptr->s->la))
			{
				fprintf(stderr,"molalities la %s %e\n", rxn_ptr->s->name, rxn_ptr->s->la);
			}
			*/
		}
		if (s_x[i]->type == EX)
		{
			s_x[i]->moles = Utilities::safe_exp(s_x[i]->lm * LOG_10);

		}
		else if (s_x[i]->type == SURF)
		{
			s_x[i]->moles = Utilities::safe_exp(s_x[i]->lm * LOG_10);

		}
		else
		{
			s_x[i]->moles = under(s_x[i]->lm) * mass_water_aq_x;
			if (s_x[i]->moles / mass_water_aq_x > 100)
			{
				log_msg(sformatf( "Overflow: %s\t%e\t%e\t%d\n",
						   s_x[i]->name,
						   (double) (s_x[i]->moles / mass_water_aq_x),
						   (double) s_x[i]->lm, iterations));

				if (iterations >= 0 && allow_overflow == FALSE)
				{
					return (ERROR);
				}
			}

		}
	}
/*
 *   other terms for diffuse layer model
 */
	if (use.Get_surface_ptr() != NULL 
		&& use.Get_surface_ptr()->Get_type() == cxxSurface::CD_MUSIC
		&& dl_type_x != cxxSurface::NO_DL)
	{
		calc_all_donnan();
	}

	struct species *s_ptr = NULL;
	for (i = 0; i < count_s_x; i++)
	{
		s_ptr = s_x[i];
		if (s_ptr->type > HPLUS && s_ptr->type != EX && s_ptr->type != SURF)
			continue;
		if (use.Get_surface_ptr() != NULL && dl_type_x != cxxSurface::NO_DL	&& s_ptr->type <= HPLUS)
		{
			total_g = 0.0;
			s_ptr->tot_dh2o_moles = 0.0;
			for (j = 0; j < (int) use.Get_surface_ptr()->Get_surface_charges().size(); j++)
			{
				cxxSurfaceCharge & charge_ref = use.Get_surface_ptr()->Get_surface_charges()[j];
				cxxSpeciesDL & dl_ref = s_diff_layer[s_ptr->number][charge_ref.Get_name()];
				cxxSurfDL & surf_dl_ref = charge_ref.Get_g_map()[s_ptr->z];

				//s_diff_layer[is][charge_ref.Get_name()] = dl_ref
				//charge_ref.Get_g_map()[s_ptr->z] = surf_dl
/*
 *   partially corrected formulation assumes mass of water in diffuse layer
 *   is insignificant. Excess is calculated on the basis of moles_water_aq_x
 *   instead of moles_water_bulk.
 */
				/* revised eq. 61 */
				dl_ref.Set_g_moles(s_ptr->moles * s_ptr->erm_ddl *
					(surf_dl_ref.Get_g() + charge_ref.Get_mass_water() / mass_water_aq_x));
				if (s_ptr->moles > 1e-30)
				{
					dl_ref.Set_dg_g_moles(s_ptr->dg * dl_ref.Get_g_moles() / s_ptr->moles);
				}

				/*
				 *  first term of 63 is summed for all surfaces in
				 *  s_ptr->tot_g_moles. This sum is then used in
				 *  the jacobian for species i
				 */
				total_g += surf_dl_ref.Get_g() + charge_ref.Get_mass_water() / mass_water_aq_x;
				/* revised eq. 63, second term */
				/* g.dg is dg/dx(-2y**2) or dg/d(ln y) */
				dl_ref.Set_dx_moles(s_ptr->moles * s_ptr->erm_ddl * surf_dl_ref.Get_dg());
				/* revised eq. 63, third term */
				dl_ref.Set_dh2o_moles(-s_ptr->moles * s_ptr->erm_ddl *
					charge_ref.Get_mass_water() / mass_water_aq_x);
				s_ptr->tot_dh2o_moles += dl_ref.Get_dh2o_moles();

				/* surface related to phase */
				dl_ref.Set_drelated_moles(s_ptr->moles * s_ptr->erm_ddl * charge_ref.Get_specific_area() *
					use.Get_surface_ptr()->Get_thickness() / mass_water_aq_x);
			}
			s_ptr->tot_g_moles = s_ptr->moles * (1 + total_g * s_ptr->erm_ddl);

			/* note that dg is for cb, act water, mu eqns */
			/* dg_total_g for mole balance eqns */
			/* dg_g_moles for surface cb */

			if (s_ptr->moles > 1e-30)
			{
				s_ptr->dg_total_g =	s_ptr->dg * s_ptr->tot_g_moles / s_ptr->moles;
			}
			else
			{
				s_ptr->dg_total_g = 0.0;
			}
			if (debug_diffuse_layer == TRUE)
			{
				output_msg(sformatf( "%s\t%e\t%e\n", s_ptr->name,
						   (double) s_ptr->moles,
						   (double) s_ptr->tot_g_moles));
				output_msg(sformatf( "\tg\n"));
				for (j = 0; j < (int) use.Get_surface_ptr()->Get_surface_charges().size(); j++)
				{
					cxxSurfaceCharge &charge_ref = use.Get_surface_ptr()->Get_surface_charges()[j];
					output_msg(sformatf( "\t%e",
							   (double) charge_ref.Get_g_map()[s_ptr->z].Get_g()));
				}
				output_msg(sformatf( "\n\tg_moles\n"));
				for (j = 0; j < (int) use.Get_surface_ptr()->Get_surface_charges().size(); j++)
				{
					cxxSurfaceCharge &charge_ref = use.Get_surface_ptr()->Get_surface_charges()[j];
					int is = s_ptr->number;
					output_msg(sformatf( "\t%e",
						(double) s_diff_layer[is][charge_ref.Get_name()].Get_g_moles()));
				}
				output_msg(sformatf( "\n\tdg\n"));
				for (j = 0; j < (int) use.Get_surface_ptr()->Get_surface_charges().size(); j++)
				{
					cxxSurfaceCharge &charge_ref = use.Get_surface_ptr()->Get_surface_charges()[j];
					output_msg(sformatf( "\t%e",
							   (double) charge_ref.Get_g_map()[s_ptr->z].Get_dg()));
				}
				output_msg(sformatf( "\n\tdx_moles\n"));
				for (j = 0; j < (int) use.Get_surface_ptr()->Get_surface_charges().size(); j++)
				{
					int is = s_ptr->number;
					cxxSurfaceCharge &charge_ref = use.Get_surface_ptr()->Get_surface_charges()[j];
					output_msg(sformatf( "\t%e",
						(double) s_diff_layer[is][charge_ref.Get_name()].Get_dx_moles()));
				}
				output_msg(sformatf( "\n\tdh2o_moles\t%e\n",
						   (double) s_ptr->tot_dh2o_moles));
				for (j = 0; j < (int) use.Get_surface_ptr()->Get_surface_charges().size(); j++)
				{
					cxxSurfaceCharge &charge_ref = use.Get_surface_ptr()->Get_surface_charges()[j];
					int is = s_ptr->number;
					output_msg(sformatf( "\t%e",
						s_diff_layer[is][charge_ref.Get_name()].Get_dh2o_moles()));
				}
				output_msg(sformatf( "\n"));
			}
		}
	}
	calc_gas_pressures();
	calc_ss_fractions();

	return (OK);
}
#ifdef SKIP
/* ---------------------------------------------------------------------- */
int Phreeqc::
molalities(int allow_overflow)
/* ---------------------------------------------------------------------- */
{
/*
 *   Calculates la for master species
 *   Calculates lm and moles from lk, lg, and la's of master species
 *   Adjusts lm of h2 and o2.
 */
	int i, j;
	LDBLE total_g;
	struct rxn_token *rxn_ptr;
/*
 *   la for master species
 */
	for (i = 0; i < count_master; i++)
	{
		if (master[i]->in == REWRITE)
		{
			master[i]->s->la = master[i]->s->lm + master[i]->s->lg;
		}
	}
	if (dl_type_x != cxxSurface::NO_DL)
	{
		s_h2o->tot_g_moles = s_h2o->moles;
		s_h2o->tot_dh2o_moles = 0.0;
	}
	for (i = 0; i < count_s_x; i++)
	{
		if (s_x[i]->type > HPLUS && s_x[i]->type != EX
			&& s_x[i]->type != SURF)
			continue;
/*
 *   lm and moles for all aqueous species
 */
		s_x[i]->lm = s_x[i]->lk - s_x[i]->lg;
		for (rxn_ptr = s_x[i]->rxn_x->token + 1; rxn_ptr->s != NULL;
			 rxn_ptr++)
		{
			s_x[i]->lm += rxn_ptr->s->la * rxn_ptr->coef;
			/*
			if (isnan(rxn_ptr->s->la))
			{
				fprintf(stderr,"molalities la %s %e\n", rxn_ptr->s->name, rxn_ptr->s->la);
			}
			*/
		}
		if (s_x[i]->type == EX)
		{
			s_x[i]->moles = Utilities::safe_exp(s_x[i]->lm * LOG_10);

		}
		else if (s_x[i]->type == SURF)
		{
			s_x[i]->moles = Utilities::safe_exp(s_x[i]->lm * LOG_10);

		}
		else
		{
			s_x[i]->moles = under(s_x[i]->lm) * mass_water_aq_x;
			if (s_x[i]->moles / mass_water_aq_x > 100)
			{
				log_msg(sformatf( "Overflow: %s\t%e\t%e\t%d\n",
						   s_x[i]->name,
						   (double) (s_x[i]->moles / mass_water_aq_x),
						   (double) s_x[i]->lm, iterations));

				if (iterations >= 0 && allow_overflow == FALSE)
				{
					return (ERROR);
				}
			}

		}
	}
/*
 *   other terms for diffuse layer model
 */
	if (use.Get_surface_ptr() != NULL && use.Get_surface_ptr()->Get_type() == cxxSurface::CD_MUSIC
		&& dl_type_x != cxxSurface::NO_DL)
		calc_all_donnan();

	for (i = 0; i < count_s_x; i++)
	{
		if (s_x[i]->type > HPLUS && s_x[i]->type != EX
			&& s_x[i]->type != SURF)
			continue;
		if (use.Get_surface_ptr() != NULL && dl_type_x != cxxSurface::NO_DL
			&& s_x[i]->type <= HPLUS)
		{
			total_g = 0.0;
			s_x[i]->tot_dh2o_moles = 0.0;
			for (j = 0; j < (int) use.Get_surface_ptr()->Get_surface_charges().size(); j++)
			{
				int is = s_x[i]->number;
				cxxSurfaceCharge & charge_ref = use.Get_surface_ptr()->Get_surface_charges()[j];
/*
 *   partially corrected formulation assumes mass of water in diffuse layer
 *   is insignificant. Excess is calculated on the basis of moles_water_aq_x
 *   instead of moles_water_bulk.
 */
				/* revised eq. 61 */
				s_diff_layer[is][charge_ref.Get_name()].Set_g_moles(s_x[i]->moles * s_x[i]->erm_ddl *
					(charge_ref.Get_g_map()[s_x[i]->z].Get_g() +
					charge_ref.Get_mass_water() /
					mass_water_aq_x));
				if (s_x[i]->moles > 1e-30)
				{
					s_diff_layer[is][charge_ref.Get_name()].Set_dg_g_moles(s_x[i]->dg * s_diff_layer[is][charge_ref.Get_name()].Get_g_moles() /
						s_x[i]->moles);
				}

				/*
				 *  first term of 63 is summed for all surfaces in
				 *  s_x[i]->tot_g_moles. This sum is then used in
				 *  the jacobian for species i
				 */
				total_g +=
					charge_ref.Get_g_map()[s_x[i]->z].Get_g() +
					charge_ref.Get_mass_water() /
					mass_water_aq_x;
				/* revised eq. 63, second term */
				/* g.dg is dg/dx(-2y**2) or dg/d(ln y) */
				s_diff_layer[is][charge_ref.Get_name()].Set_dx_moles(
					s_x[i]->moles * s_x[i]->erm_ddl *
					charge_ref.Get_g_map()[s_x[i]->z].Get_dg());
				/* revised eq. 63, third term */
				s_diff_layer[is][charge_ref.Get_name()].Set_dh2o_moles(
					-s_x[i]->moles * s_x[i]->erm_ddl *
					charge_ref.Get_mass_water() /
					mass_water_aq_x);
				s_x[i]->tot_dh2o_moles += s_diff_layer[is][charge_ref.Get_name()].Get_dh2o_moles();

				/* surface related to phase */
				s_diff_layer[is][charge_ref.Get_name()].Set_drelated_moles(
					s_x[i]->moles * s_x[i]->erm_ddl *
					charge_ref.Get_specific_area() *
					use.Get_surface_ptr()->Get_thickness() / mass_water_aq_x);
			}
			s_x[i]->tot_g_moles =
				s_x[i]->moles * (1 + total_g /* s_x[i]->erm_ddl */ );

			/* note that dg is for cb, act water, mu eqns */
			/* dg_total_g for mole balance eqns */
			/* dg_g_moles for surface cb */

			if (s_x[i]->moles > 1e-30)
			{
				s_x[i]->dg_total_g =
					s_x[i]->dg * s_x[i]->tot_g_moles / s_x[i]->moles;
			}
			else
			{
				s_x[i]->dg_total_g = 0.0;
			}
			if (debug_diffuse_layer == TRUE)
			{
				output_msg(sformatf( "%s\t%e\t%e\n", s_x[i]->name,
						   (double) s_x[i]->moles,
						   (double) s_x[i]->tot_g_moles));
				output_msg(sformatf( "\tg\n"));
				for (j = 0; j < (int) use.Get_surface_ptr()->Get_surface_charges().size(); j++)
				{
					cxxSurfaceCharge &charge_ref = use.Get_surface_ptr()->Get_surface_charges()[j];
					output_msg(sformatf( "\t%e",
							   (double) charge_ref.Get_g_map()[s_x[i]->z].Get_g()));
				}
				output_msg(sformatf( "\n\tg_moles\n"));
				for (j = 0; j < (int) use.Get_surface_ptr()->Get_surface_charges().size(); j++)
				{
					cxxSurfaceCharge &charge_ref = use.Get_surface_ptr()->Get_surface_charges()[j];
					int is = s_x[i]->number;
					output_msg(sformatf( "\t%e",
						(double) s_diff_layer[is][charge_ref.Get_name()].Get_g_moles()));
				}
				output_msg(sformatf( "\n\tdg\n"));
				for (j = 0; j < (int) use.Get_surface_ptr()->Get_surface_charges().size(); j++)
				{
					cxxSurfaceCharge &charge_ref = use.Get_surface_ptr()->Get_surface_charges()[j];
					output_msg(sformatf( "\t%e",
							   (double) charge_ref.Get_g_map()[s_x[i]->z].Get_dg()));
				}
				output_msg(sformatf( "\n\tdx_moles\n"));
				for (j = 0; j < (int) use.Get_surface_ptr()->Get_surface_charges().size(); j++)
				{
					int is = s_x[i]->number;
					cxxSurfaceCharge &charge_ref = use.Get_surface_ptr()->Get_surface_charges()[j];
					output_msg(sformatf( "\t%e",
						(double) s_diff_layer[is][charge_ref.Get_name()].Get_dx_moles()));
				}
				output_msg(sformatf( "\n\tdh2o_moles\t%e\n",
						   (double) s_x[i]->tot_dh2o_moles));
				for (j = 0; j < (int) use.Get_surface_ptr()->Get_surface_charges().size(); j++)
				{
					cxxSurfaceCharge &charge_ref = use.Get_surface_ptr()->Get_surface_charges()[j];
					int is = s_x[i]->number;
					output_msg(sformatf( "\t%e",
						s_diff_layer[is][charge_ref.Get_name()].Get_dh2o_moles()));
				}
				output_msg(sformatf( "\n"));
			}
		}
	}
	calc_gas_pressures();
	calc_ss_fractions();

	return (OK);
}
#endif
/* ---------------------------------------------------------------------- */
int Phreeqc::
calc_gas_pressures(void)
/* ---------------------------------------------------------------------- */
{
	int n_g = 0;
	LDBLE lp, V_m = 0;
	struct rxn_token *rxn_ptr;
	std::vector<struct phase *> phase_ptrs;
	bool PR = false, pr_done = false;
/*
 *   moles and partial pressures for gases
 */
	if (use.Get_gas_phase_ptr() == NULL)
		return (OK);
	cxxGasPhase *gas_phase_ptr = use.Get_gas_phase_ptr();
	if (gas_phase_ptr->Get_type() == cxxGasPhase::GP_VOLUME && (gas_phase_ptr->Get_pr_in() || force_numerical_fixed_volume) && numerical_fixed_volume)
	{
		if (iterations > 2)
			return calc_fixed_volume_gas_pressures();
		else
			return OK;
	}
	if (iterations > 2 && gas_phase_ptr->Get_type() == cxxGasPhase::GP_VOLUME)
	{
		gas_phase_ptr->Set_total_moles(0);
	}

	for (size_t i = 0; i < gas_phase_ptr->Get_gas_comps().size(); i++)
	{
		const cxxGasComp * gas_comp_ptr = &(gas_phase_ptr->Get_gas_comps()[i]);
		int j;
		struct phase *phase_ptr = phase_bsearch(gas_comp_ptr->Get_phase_name().c_str(), &j, FALSE);
		if (phase_ptr->in == TRUE)
		{
			phase_ptrs.push_back(phase_ptr);
			if (!PR && phase_ptr->t_c > 0 && phase_ptr->p_c > 0)
				PR = true;
			n_g++;
		}
		if (iterations > 2 && gas_phase_ptr->Get_type() == cxxGasPhase::GP_VOLUME)
		{
			gas_phase_ptr->Set_total_moles(gas_phase_ptr->Get_total_moles() + phase_ptr->moles_x);
		}
	}
	if (gas_phase_ptr->Get_type() == cxxGasPhase::GP_PRESSURE)
	{
		if (PR /*&& gas_unknown->gas_phase->total_p > 1 */ && iterations > 0)
		{
			calc_PR(phase_ptrs, gas_phase_ptr->Get_total_p(), tk_x, 0);
		}
	} else
	{
		if (PR)
		{
			if (gas_phase_ptr->Get_total_moles() > 0)
			{
				V_m = gas_phase_ptr->Get_volume() / gas_phase_ptr->Get_total_moles();
				if (V_m < 0.016)
				{
						V_m = 0.016;
				} else if (V_m > 1e4)
				{
					V_m = 1e4;
				}
				if (V_m < 0.02)
					V_m = (8. * gas_phase_ptr->Get_v_m() + V_m) / 9;
				else if (V_m < 0.03)
					V_m = (6. * gas_phase_ptr->Get_v_m() + V_m) / 7;
				else if (V_m < 0.05)
					V_m = (4. * gas_phase_ptr->Get_v_m() + V_m) / 5;
				else if (V_m < 0.07)
					V_m = (2. * gas_phase_ptr->Get_v_m() + V_m) / 3;
				else
					V_m = (1. * gas_phase_ptr->Get_v_m() + V_m) / 2;
				if (iterations > 99 && numerical_fixed_volume == false)
				{
					//V_m *= 1; /* debug */
					numerical_fixed_volume = true;
					//switch_numerical = true;
					warning_msg
						("Numerical method failed, switching to numerical derivatives.");
					prep();
					//switch_numerical = false;
				}
			} else
				V_m = 1.0;
			calc_PR(phase_ptrs, 0, tk_x, V_m);
			pr_done = true;
		} else
		{
			gas_phase_ptr->Set_total_p(0);
		}
	}

	gas_phase_ptr->Set_total_moles(0);
	std::vector<cxxGasComp> gas_comps;
	for (size_t i = 0; i < gas_phase_ptr->Get_gas_comps().size(); i++)
	{
		const cxxGasComp *gas_comp = &(gas_phase_ptr->Get_gas_comps()[i]);
		int j;
		struct phase *phase_ptr = phase_bsearch(gas_comp->Get_phase_name().c_str(), &j, FALSE);
		if (phase_ptr->in == TRUE)
		{
			lp = -phase_ptr->lk;
			for (rxn_ptr = phase_ptr->rxn_x->token + 1; rxn_ptr->s != NULL;
				 rxn_ptr++)
			{
				lp += rxn_ptr->s->la * rxn_ptr->coef;
			}
			phase_ptr->p_soln_x = exp(LOG_10 * (lp - phase_ptr->pr_si_f));
			if (!strcmp(phase_ptr->name, "H2O(g)") && phase_ptr->p_soln_x > 90)
					phase_ptr->p_soln_x = 90;

			if (gas_phase_ptr->Get_type() == cxxGasPhase::GP_PRESSURE)
			{
				phase_ptr->moles_x = phase_ptr->p_soln_x *
					gas_unknown->moles / gas_phase_ptr->Get_total_p();
				phase_ptr->fraction_x =
					phase_ptr->moles_x / gas_unknown->moles;
			}
			else
			{
				if (pr_done)
				{
					lp = phase_ptr->p_soln_x / gas_phase_ptr->Get_total_p() *
						gas_phase_ptr->Get_volume() / V_m;
					if (lp > 0)
						phase_ptr->moles_x = lp;
					if (iterations > 50)
					{
						lp *= 1.0; /* debug */
					}
				} else
				{
					phase_ptr->moles_x = phase_ptr->p_soln_x *
						gas_phase_ptr->Get_volume() / (R_LITER_ATM * tk_x);
					gas_phase_ptr->Set_total_p(gas_phase_ptr->Get_total_p() + phase_ptr->p_soln_x);
				}
				gas_phase_ptr->Set_total_moles(gas_phase_ptr->Get_total_moles() +
					phase_ptr->moles_x);
			}
		}
		else
		{
			phase_ptr->moles_x = 0;
			phase_ptr->fraction_x = 0;
		}
	}

	if (gas_phase_ptr->Get_type() == cxxGasPhase::GP_VOLUME && !PR)
	{
		/*
		 * Fixed-volume gas phase reacting with a solution
		 * Change pressure used in logK to pressure of gas phase
		 */
		if (gas_phase_ptr->Get_total_p() > 1500)
		{
			gas_phase_ptr->Set_total_moles(0);
			for (size_t i = 0; i < gas_phase_ptr->Get_gas_comps().size(); i++)
			{
				const cxxGasComp *gas_comp = &(gas_phase_ptr->Get_gas_comps()[i]);
				int j;
				struct phase *phase_ptr = phase_bsearch(gas_comp->Get_phase_name().c_str(), &j, FALSE);
				if (phase_ptr->in == TRUE)
				{
					phase_ptr->moles_x *= 1500.0 / gas_phase_ptr->Get_total_p();
					gas_phase_ptr->Set_total_moles(gas_phase_ptr->Get_total_moles() +
						phase_ptr->moles_x);
				}
			}
			gas_phase_ptr->Set_total_p(1500.0);
		}
	}

	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
calc_ss_fractions(void)
/* ---------------------------------------------------------------------- */
{
	LDBLE moles, n_tot;
/*
 *   moles and lambdas for solid solutions
 */
	if (ss_unknown == NULL)
		return (OK);
/*
 *  Calculate mole fractions and log lambda and derivative factors
 */
	if (use.Get_ss_assemblage_ptr() == NULL)
		return (OK);
	std::vector<cxxSS *> ss_ptrs = use.Get_ss_assemblage_ptr()->Vectorize();
	for (size_t i = 0; i < ss_ptrs.size(); i++)
	{
		cxxSS *ss_ptr = ss_ptrs[i];
		n_tot = 0;
		for (size_t k = 0; k < ss_ptr->Get_ss_comps().size(); k++)
		{
			cxxSScomp *comp_ptr = &(ss_ptr->Get_ss_comps()[k]);
			moles = comp_ptr->Get_moles();
			if (moles < 0)
			{
				moles = MIN_TOTAL_SS;
				comp_ptr->Set_initial_moles(moles);
			}
			n_tot += moles;
		}
		ss_ptr->Set_total_moles(n_tot);
		for (size_t k = 0; k < ss_ptr->Get_ss_comps().size(); k++)
		{
			cxxSScomp *comp_ptr = &(ss_ptr->Get_ss_comps()[k]);
			int l;
			struct phase *phase_ptr = phase_bsearch(comp_ptr->Get_name().c_str(), &l, FALSE);
			moles = comp_ptr->Get_moles();
			if (moles < 0)
			{
				moles = MIN_TOTAL_SS;
			}
			comp_ptr->Set_fraction_x(moles / n_tot);
			comp_ptr->Set_log10_fraction_x(log10(moles / n_tot));

			/* all mb and jacobian items must be in x or phase to be static between models */
			phase_ptr->log10_fraction_x = comp_ptr->Get_log10_fraction_x();
		}
		if (ss_ptr->Get_a0() != 0.0 || ss_ptr->Get_a1() != 0)
		{
			ss_binary(ss_ptr);
		}
		else
		{
			ss_ideal(ss_ptr);
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
ss_binary(cxxSS *ss_ptr)
/* ---------------------------------------------------------------------- */
{
	LDBLE nb, nc, n_tot, xb, xc, dnb, dnc, l_a0, l_a1;
	LDBLE xb2, xb3, xb4, xc2, xc3;
	LDBLE xb1, xc1;
/*
 * component 0 is major component
 * component 1 is minor component
 * xb is the mole fraction of second component (formerly trace)
 * xc is the mole fraction of first component (formerly major)
*/
/*
 *  Calculate mole fractions and log lambda and derivative factors
 */
	n_tot = ss_ptr->Get_total_moles();
	cxxSScomp *comp0_ptr = &(ss_ptr->Get_ss_comps()[0]);
	cxxSScomp *comp1_ptr = &(ss_ptr->Get_ss_comps()[1]);
	int l;
	struct phase *phase0_ptr = phase_bsearch(comp0_ptr->Get_name().c_str(), &l, FALSE);
	struct phase *phase1_ptr = phase_bsearch(comp1_ptr->Get_name().c_str(), &l, FALSE);

	nc = comp0_ptr->Get_moles();
	xc = nc / n_tot;
	nb = comp1_ptr->Get_moles();
	xb = nb / n_tot;
/*
 *   In miscibility gap
 */
	l_a0 = ss_ptr->Get_a0();
	l_a1 = ss_ptr->Get_a1();
	if (ss_ptr->Get_miscibility() && xb > ss_ptr->Get_xb1()
		&& xb < ss_ptr->Get_xb2())
	{
		xb1 = ss_ptr->Get_xb1();
		xc1 = 1.0 - xb1;
		comp0_ptr->Set_fraction_x(xc1);
		comp0_ptr->Set_log10_fraction_x(log10(xc1));
		phase0_ptr->log10_fraction_x =
			comp0_ptr->Get_log10_fraction_x();

		comp1_ptr->Set_fraction_x(xb1);
		comp1_ptr->Set_log10_fraction_x(log10(xb1));
		phase1_ptr->log10_fraction_x =
			comp1_ptr->Get_log10_fraction_x();

		comp0_ptr->Set_log10_lambda(
			xb1 * xb1 * (l_a0 - l_a1 * (3 - 4 * xb1)) / LOG_10);
		phase0_ptr->log10_lambda =
			comp0_ptr->Get_log10_lambda();

		comp1_ptr->Set_log10_lambda(
			xc1 * xc1 * (l_a0 + l_a1 * (4 * xb1 - 1)) / LOG_10);
		phase1_ptr->log10_lambda =
			comp1_ptr->Get_log10_lambda();

		comp0_ptr->Set_dnb(0);
		comp0_ptr->Set_dnc(0);
		comp1_ptr->Set_dnb(0);
		comp1_ptr->Set_dnc(0);
		phase0_ptr->dnb = 0;
		phase0_ptr->dnc = 0;
		phase1_ptr->dnb = 0;
		phase1_ptr->dnc = 0;
	}
	else
	{
/*
 *   Not in miscibility gap
 */
		comp0_ptr->Set_fraction_x(xc);
		comp0_ptr->Set_log10_fraction_x(log10(xc));
		phase0_ptr->log10_fraction_x =
			comp0_ptr->Get_log10_fraction_x();

		comp1_ptr->Set_fraction_x(xb);
		comp1_ptr->Set_log10_fraction_x(log10(xb));
		phase1_ptr->log10_fraction_x =
			comp1_ptr->Get_log10_fraction_x();

		comp0_ptr->Set_log10_lambda(
			xb * xb * (l_a0 - l_a1 * (3 - 4 * xb)) / LOG_10);
		phase0_ptr->log10_lambda =
			comp0_ptr->Get_log10_lambda();

		comp1_ptr->Set_log10_lambda(
			xc * xc * (l_a0 + l_a1 * (4 * xb - 1)) / LOG_10);
		phase1_ptr->log10_lambda =
			comp1_ptr->Get_log10_lambda();

		xc2 = xc * xc;
		xc3 = xc2 * xc;
		xb2 = xb * xb;
		xb3 = xb2 * xb;
		xb4 = xb3 * xb;
		xb4 = xb4;
		xc3 = xc3;

		/* used derivation that did not substitute x2 = 1-x1 */

		/* first component, df1/dn1 */
		dnc = 2 * l_a0 * xb2 + 12 * l_a1 * xc * xb2 + 6 * l_a1 * xb2;
		phase0_ptr->dnc = -xb / nc + dnc / n_tot;


		/* first component, df1/dn2 */
		dnb =
			1 - 2 * l_a0 * xb + 2 * l_a0 * xb2 + 8 * l_a1 * xc * xb -
			12 * l_a1 * xc * xb2 - 2 * l_a1 * xb + 2 * l_a1 * xb2;
		phase0_ptr->dnb = dnb / n_tot;

		/* second component, df2/dn1 */
		dnc =
			1 - 2 * l_a0 * xc + 2 * l_a0 * xc2 - 8 * l_a1 * xb * xc +
			12 * l_a1 * xb * xc2 + 2 * l_a1 * xc - 2 * l_a1 * xc2;
		phase1_ptr->dnc = dnc / n_tot;

		/* second component, df2/dn2 */
		dnb = 2 * l_a0 * xc2 + 12 * l_a1 * xb * xc2 - 6 * l_a1 * xc2;
		phase1_ptr->dnb = -xc / nb + dnb / n_tot;

	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
ss_ideal(cxxSS *ss_ptr)
/* ---------------------------------------------------------------------- */
{
	LDBLE n_tot, n_tot1;

/*
 * component 0 is major component
 * component 1 is minor component
 * xb is the mole fraction of second component (formerly trace)
 * xc is the mole fraction of first component (formerly major)
*/
/*
 *  Calculate mole fractions and log lambda and derivative factors
 */
	n_tot = ss_ptr->Get_total_moles();

/*
 *   Ideal solid solution
 */
	ss_ptr->Set_dn(1.0 / n_tot);
	for (size_t k = 0; k < ss_ptr->Get_ss_comps().size(); k++)
	{
		cxxSScomp *compk_ptr = &(ss_ptr->Get_ss_comps()[k]);
		int l;
		struct phase *phasek_ptr = phase_bsearch(compk_ptr->Get_name().c_str(), &l, FALSE);
		n_tot1 = 0;
		for (size_t j = 0; j < ss_ptr->Get_ss_comps().size(); j++)
		{
			cxxSScomp *compj_ptr = &(ss_ptr->Get_ss_comps()[j]);
			if (j != k)
			{
				n_tot1 += compj_ptr->Get_moles();
			}
		}
		compk_ptr->Set_log10_lambda(0);
		compk_ptr->Set_log10_lambda(0);

		compk_ptr->Set_dnb(-(n_tot1) / (compk_ptr->Get_moles() * n_tot));
		phasek_ptr->dnb = compk_ptr->Get_dnb();

		compk_ptr->Set_dn(ss_ptr->Get_dn());
		phasek_ptr->dn = ss_ptr->Get_dn();
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
reset(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks deltas (changes to unknowns) to make sure they are reasonable
 *   Scales deltas if necessary
 *   Updates unknowns with deltas
 */

	int i;
	int converge;
	LDBLE up, down;
	LDBLE d;
	LDBLE factor, f0;
	LDBLE sum_deltas;
	LDBLE step_up;
	LDBLE mu_calc;
	LDBLE old_moles;
	int warning;
/*
 *   Calculate interphase mass transfers
 */
	cxxPPassemblageComp * comp_ptr;

	step_up = log(step_size_now);
	factor = 1.;

	if ((pure_phase_unknown != NULL || ss_unknown != NULL)
		&& calculating_deriv == FALSE)
	{
/*
 *   Don`t take out more mineral than is present
 */
		for (i = 0; i < count_unknowns; i++)
		{
			if (x[i]->type == PP || x[i]->type == SS_MOLES)
			{

				if (delta[i] < -1e8)
				{
					delta[i] = -10.;
				}
				else if (delta[i] > 1e8)
				{
					delta[i] = 10;
				}
				if (x[i]->dissolve_only == TRUE)
				{
					assert (x[i]->type == PP);
					//comp_ptr = pp_assemblage_ptr->Find(x[i]->pp_assemblage_comp_name);
					comp_ptr = (cxxPPassemblageComp *) x[i]->pp_assemblage_comp_ptr;
					assert(comp_ptr);
					if ((delta[i] < 0.0)
						&& (-delta[i] >
							(comp_ptr->Get_initial_moles() - x[i]->moles)))
					{
						if ((comp_ptr->Get_initial_moles() - x[i]->moles) !=
							0.0)
						{
							f0 = fabs(delta[i] /
									  (comp_ptr->Get_initial_moles() -
									   x[i]->moles));
							if (f0 > factor)
							{
								if (debug_model == TRUE)
								{
									output_msg(sformatf(
											   "%-10.10s, Precipitating too much dissolve_only mineral.\tDelta %e\tCurrent %e\tInitial %e\n",
											   x[i]->description,
											   (double) delta[i],
											   (double) x[i]->moles,
											   (double) comp_ptr->Get_initial_moles()));
								}
								factor = f0;
							}
						}
						else
						{
							if (debug_model == TRUE)
							{
								output_msg(sformatf(
										   "%-10.10s, Precipitating dissolve_only mineral.\tDelta %e\n",
										   x[i]->description,
										   (double) delta[i]));
							}
							delta[i] = 0;
						}
					}
				}
				if ( /* delta[i] > 0.0 && */ x[i]->moles > 0.0
					&& delta[i] > x[i]->moles)
				{
					f0 = delta[i] / x[i]->moles;
					if (f0 > factor)
					{
						if (debug_model == TRUE)
						{
							output_msg(sformatf(
									   "%-10.10s, Removing more than total mineral.\t%f\n",
									   x[i]->description, (double) f0));
						}
						factor = f0;
					}
				}
				else if (delta[i] > 0.0 && x[i]->moles <= 0.0)
				{
					if (debug_model == TRUE)
					{
						output_msg(sformatf(
								   "%-10.10s\tDelta: %e\tMass: %e   "
								   "Dissolving mineral with 0.0 mass.\n ",
								   x[i]->description, (double) delta[i],
								   (double) x[i]->moles));
					}
					delta[i] = 0.0;
				}
				else if (x[i]->ss_comp_name != NULL && delta[i] < -x[i]->phase->delta_max)
				// Uses delta_max computed in step
				// delta_max is the maximum amount of the mineral that could form based
				// on the limiting element in the system
				{
					f0 = -delta[i] / x[i]->phase->delta_max;
					if (f0 > factor)
					{
						if (debug_model == TRUE)
						{
							output_msg(sformatf(
									   "%-10.10s, Precipitating too much mineral.\t%f\n",
									   x[i]->description, (double) f0));
						}
						factor = f0;
					}
				}
			}
		}
	}
/*
 *   Calculate change in element concentrations due to pure phases and gases
 */
	warning = 0;
	for (i = 0; i < count_unknowns; i++)
	{

		/*if (isnan(delta[i]))*/
		if (!PHR_ISFINITE((double) delta[i]))
		{
			warning ++;
			delta[i] = 0;
		}

	}
	if (warning > 0)
	{
		error_string = sformatf( "%d delta equal NaN\n", warning);
		warning_msg(error_string);
	}
	if (pure_phase_unknown != NULL || gas_unknown != NULL
		|| ss_unknown != NULL)
	{
		for (i = 0; i < count_unknowns; i++)
		{
			x[i]->delta = 0.0;
		}

		for (i = 0; i < count_sum_delta; i++)
		{
			*sum_delta[i].target += *sum_delta[i].source * sum_delta[i].coef;
		}

/*
 *   Apply factor from minerals to deltas
 */


		for (i = 0; i < count_unknowns; i++)
		{
			x[i]->delta /= factor;
			if (x[i]->type == PP || x[i]->type == SS_MOLES)
				delta[i] /= factor;
		}

	}


/*
 *   Calc factor for mass balance equations for aqueous unknowns
 */
	factor = 1.0;
	sum_deltas = 0.0;
	for (i = 0; i < count_unknowns; i++)

	{
		/* fixes underflow problem on Windows */
		if (delta[i] > 0)
		{
			sum_deltas += delta[i];
		}
		else
		{
			sum_deltas -= delta[i];

		}
		/*sum_deltas += fabs(delta[i]); */
		if (calculating_deriv == FALSE)
		{
			up = step_up;
			down = up;
			if (x[i]->type <= SOLUTION_PHASE_BOUNDARY)
			{
				up = step_up;
				down = 1.3 * up;
			}
			else if (x[i]->type == MU)
			{
				up = 100 * mu_x;
				down = mu_x;
			}
			else if (x[i]->type == AH2O)
			{
				down = up;
				if (pitzer_model || sit_model)
				  {
				    up = 0.05;
				    down = -0.03;
				  }
			}
			else if (x[i]->type == MH)
			{
				up = log(pe_step_size_now);
				down = 1.3 * up;
			}
			else if (x[i]->type == MH2O)
			{
				/* ln gH2O + delta; ln(gH2O*delta); */
				/*
				   up = log(10.);
				   down = log(4.);
				 */
				up = log(1.3);
				down = log(1.2);

			}
			else if (x[i]->type == PP)
			{
				continue;
			}
			else if (x[i]->type == GAS_MOLES)
			{
				up = 1000. * x[i]->moles;
				if (up <= 0.0)
					up = 1e-1;
				if (up >= 1.0)
					up = 1.;
				down = x[i]->moles;
			}
			else if (x[i]->type == SS_MOLES)
			{
				continue;
			}
			else if (x[i]->type == EXCH)
			{
				up = step_up;
				down = 1.3 * up;
			}
			else if (x[i]->type == SURFACE)
			{
				up = step_up;
				down = 1.3 * up;
			}
			else if (x[i]->type == PITZER_GAMMA)
			{
				up = step_up;
				if (up > 1)
				{
					up = 0.7;
				}
				down = 1.3 * up;
			}
			else if (x[i]->type == SURFACE_CB || x[i]->type == SURFACE_CB1
					 || x[i]->type == SURFACE_CB2)
			{
				up = step_up;
				down = 1.3 * up;
				/*
				   up = 1.3;
				   down = 1.2;
				 */
			}

			if (delta[i] > 0.0)
			{
				f0 = delta[i] / up;
				if (f0 > factor)
				{
					if (debug_model == TRUE)
					{
						output_msg(sformatf( "%-10.10s\t%f\n",
								   x[i]->description, (double) f0));
					}
					factor = f0;
				}
			}
			else
			{
				f0 = delta[i] / (-down);
				if (f0 > factor)
				{
					if (debug_model == TRUE)
					{
						output_msg(sformatf( "%-10.10s\t%f\n",
								   x[i]->description, (double) f0));
					}
					factor = f0;
				}
			}
		}
	}

	/*converge=TRUE; */

	if (debug_model == TRUE)
	{
		output_msg(sformatf( "\nSum of deltas: %12.6f\n",
				   (double) sum_deltas));
		output_msg(sformatf( "Factor: %12.4e\n", (double) factor));
	}
	factor = 1.0 / factor;

	for (i = 0; i < count_unknowns; i++)
	{
		if (x[i]->type != PP && x[i]->type != SS_MOLES)
			delta[i] *= factor;
	}

	warning = 0;
	for (i = 0; i < count_unknowns; i++)
	{

		/*if (isnan(delta[i]))*/
		if (!PHR_ISFINITE((double) delta[i]))
		{
			warning ++;
			delta[i] = 0;
		}

	}
	if (warning > 0)
	{
		error_string = sformatf( "%d delta equal NaN after scaling\n", warning);
		warning_msg(error_string);
	}
/*
 *   Solution mass balances: MB, ALK, CB, SOLUTION_PHASE_BOUNDARY
 */
	for (i = 0; i < count_unknowns; i++)
	{
		if (x[i]->type == MB || x[i]->type == ALK || x[i]->type == EXCH
			|| x[i]->type == SURFACE)
		{
			/*if ( fabs(delta[i]) >= epsilon ) converge = FALSE; */
			d = delta[i] / LOG_10;
			/* surface */
			if (x[i]->type == SURFACE)
			{
				cxxSurfaceComp *comp_ptr = use.Get_surface_ptr()->Find_comp(x[i]->surface_comp);
				old_moles = x[i]->moles;
				if (x[i]->phase_unknown != NULL)
				{
					x[i]->moles = comp_ptr->Get_phase_proportion() *
						(x[i]->phase_unknown->moles -
						 delta[x[i]->phase_unknown->number]);
					if (x[i]->phase_unknown->moles -
						delta[x[i]->phase_unknown->number] <=
						MIN_RELATED_SURFACE)
					{
						x[i]->moles = 0.0;
						if (fabs(x[i]->f) > MIN_RELATED_SURFACE)
						{
							x[i]->master[0]->s->la -= 5.;
						}
					}
					if (old_moles <= 0 && x[i]->moles > 0)
					{
						x[i]->master[0]->s->la = log10(x[i]->moles) - 5.;
					}
				}
				else if (comp_ptr->Get_phase_name().size() > 0)
				{
					/* probably initial surface calculation */
					if (x[i]->moles <= MIN_RELATED_SURFACE)
					{
						x[i]->moles = 0.0;
						if (fabs(x[i]->f) > MIN_RELATED_SURFACE)
						{
							x[i]->master[0]->s->la -= 5.;
						}
					}
				}
			}
			/* exch */
			if (x[i]->type == EXCH && x[i]->moles <= MIN_RELATED_SURFACE)
			{
				x[i]->moles = 0.0;
				if (fabs(x[i]->f) > MIN_RELATED_SURFACE)
				{
					x[i]->master[0]->s->la -= 5.;
				}
			}
			if (debug_model == TRUE)
			{
				output_msg(sformatf(
						   "%-10.10s %-9s%10.5f   %-9s%10.5f   %-6s%10.2e   "
						   "%-8s%10.2e\n", x[i]->description, "old la",
						   (double) x[i]->master[0]->s->la, "new la",
						   (double) x[i]->master[0]->s->la + (double) d,
						   "delta", (double) delta[i], "delta/c", (double) d));
			}
			x[i]->master[0]->s->la += d;
			//if (x[i]->master[0]->s->la < (double) (DBL_MIN_10_EXP + 10))
			//	x[i]->master[0]->s->la = (double) (DBL_MIN_10_EXP + 10);

/*
 * Surface charge balance
 */

		}
		else if (x[i]->type == SURFACE_CB || x[i]->type == SURFACE_CB1
				 || x[i]->type == SURFACE_CB2)
		{
			d = delta[i] / LOG_10;
			cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[i]->surface_charge);
			if (x[i]->phase_unknown != NULL)
			{
				charge_ptr->Set_grams(	
					(x[i]->phase_unknown->moles -
					 delta[x[i]->phase_unknown->number]));
				if (charge_ptr->Get_grams() <= MIN_RELATED_SURFACE)
				{
					charge_ptr->Set_grams(0.0);
				}
			}
			if (charge_ptr->Get_grams() <= MIN_RELATED_SURFACE)
			{
				charge_ptr->Set_grams(0.0);
			}
			x[i]->related_moles = charge_ptr->Get_grams();
			if (debug_model == TRUE)
			{
				output_msg(sformatf(
						   "%-10.10s %-9s%10.5f   %-9s%10.5f   %-6s%10.2e\n",
						   x[i]->description, "old f*psi",
						   (double) x[i]->master[0]->s->la, "new f*psi",
						   (double) x[i]->master[0]->s->la + (double) d,
						   "delta", (double) d));
			}

			x[i]->master[0]->s->la += d;

			/* recalculate g's for component */
			if (dl_type_x != cxxSurface::NO_DL
				&& (use.Get_surface_ptr()->Get_type() == cxxSurface::DDL || use.Get_surface_ptr()->Get_type() == cxxSurface::CCM
					|| (use.Get_surface_ptr()->Get_type() == cxxSurface::CD_MUSIC
						&& x[i]->type == SURFACE_CB2)))
			{
				if (debug_diffuse_layer == TRUE)
				{
					output_msg(sformatf(
							   "\ncharge, old g, new g, dg*delta,"
							   " dg, delta\n"));
				}
				std::map<LDBLE, cxxSurfDL>::iterator jit;
				for (jit = charge_ptr->Get_g_map().begin(); jit != charge_ptr->Get_g_map().end(); jit++)
				{
					if (debug_diffuse_layer == TRUE)
					{
						output_msg(sformatf(
								   "\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\n",
								   
								   (double) jit->second.Get_g(),
								   (double) jit->second.Get_g() +
								   (double) (jit->second.Get_dg() *
											 delta[i]),
								   (double) (jit->second.Get_dg() *
											 delta[i]),
								   (double) jit->second.Get_dg(),
								   (double) delta[i]));
					}
					if (use.Get_surface_ptr()->Get_dl_type() != cxxSurface::DONNAN_DL)
					{
						jit->second.Set_g(jit->second.Get_g() +
							jit->second.Get_dg() * delta[i]);
					}
				}
				if (use.Get_surface_ptr()->Get_dl_type() == cxxSurface::DONNAN_DL)
				{
					calc_all_donnan();
				}
			}

/*   Solution phase boundary */
		}
		else if (x[i]->type == SOLUTION_PHASE_BOUNDARY)
		{
			/*if (fabs(delta[i]) > epsilon) converge=FALSE; */
			d = delta[i] / LOG_10;
			if (debug_model == TRUE)
			{
				output_msg(sformatf(
						   "%-10.10s %-9s%10.5f   %-9s%10.5f   %-6s%10.2e   %-8s%10.2e\n",
						   x[i]->description, "old la",
						   (double) x[i]->master[0]->s->la, "new la",
						   (double) (x[i]->master[0]->s->la + d), "delta",
						   (double) delta[i], "delta/c", (double) d));
			}
			x[i]->master[0]->s->la += d;
/*   Charge balance */
		}
		else if (x[i]->type == CB)
		{
			/*if (fabs(delta[i]) > epsilon * mu_x * mass_water_aq_x ) converge=FALSE; */
			d = delta[i] / LOG_10;
			if (debug_model == TRUE)
			{
				output_msg(sformatf(
						   "%-10.10s %-9s%10.5f   %-9s%10.5f   %-6s%10.2e   %-8s%10.2e\n",
						   x[i]->description, "old la",
						   (double) x[i]->master[0]->s->la, "new la",
						   (double) (x[i]->master[0]->s->la + d), "delta",
						   (double) delta[i], "delta/c", (double) d));
			}
			x[i]->master[0]->s->la += d;
/*   Ionic strength */
		}
		else if (x[i]->type == MU)
		{
			// using straight ionic strength equation
			{
				mu_calc = 0.5 * mu_unknown->f / mass_water_aq_x;
			}
			if (debug_model == TRUE)
			{
				output_msg(sformatf( "Calculated mu: %e\n",
						   (double) mu_calc));
				output_msg(sformatf(
						   "%-10.10s %-9s%10.5f   %-9s%10.5f   %-6s%10.2e\n",
						   x[i]->description, "old mu", (double) mu_x,
						   "new mu", (double) (mu_x + delta[i]), "delta",
						   (double) delta[i]));
			}
			d = mu_x + delta[i];
			if (d < 1e-11)
			{
				delta[i] = sqrt(mu_calc * mu_x) - mu_x;
				mu_x = sqrt(mu_calc * mu_x);
			}
			else
			{
				mu_x += delta[i];
			}
			//if (mu_x <= 1e-8)
			//{
			//	mu_x = 1e-8;
			//}
/*   Activity of water */
		}
		else if (x[i]->type == AH2O)
		{
			/*if (pitzer_model == TRUE && full_pitzer == FALSE) continue; */
			/*if (fabs(delta[i]) > epsilon) converge=FALSE; */
			d = delta[i] / LOG_10;
			if (debug_model == TRUE)
			{
				output_msg(sformatf(
						   "%-10.10s %-9s%10.5f   %-9s%10.5f   %-6s%10.2e   %-8s%10.2e\n",
						   x[i]->description, "old la",
						   (double) x[i]->master[0]->s->la, "new la",
						   (double) (x[i]->master[0]->s->la + d), "delta",
						   (double) delta[i], "delta/c", (double) d));
			}
			s_h2o->la += d;
			ah2o_x = exp(s_h2o->la * LOG_10);
/*   pe */
		}
		else if (x[i]->type == MH)
		{
			/*if (fabs(delta[i]) > epsilon) converge=FALSE; */
			d = delta[i] / LOG_10;
			if (debug_model == TRUE)
			{
				output_msg(sformatf(
						   "%-10.10s %-9s%10.5f   %-9s%10.5f   %-6s%10.2e   %-8s%10.2e\n",
						   x[i]->description, "old pe",
						   (double) x[i]->master[0]->s->la, "new pe",
						   (double) (x[i]->master[0]->s->la + d), "delta",
						   (double) delta[i], "delta/c", (double) d));
			}
			s_eminus->la += d;
/*   Mass of water */
		}
		else if (x[i]->type == MH2O)
		{
			if (mass_water_switch == TRUE)
				continue;
			/*if (fabs(delta[i]) > epsilon * mass_water_aq_x) converge=FALSE; */
			/* ln(gh2o) + delta, log(gh2o) + d, gh2o * 10**d */
			d = exp(delta[i]);
			if (debug_model == TRUE)
			{
				output_msg(sformatf(
						   "%-10.10s %-9s%10.2e   %-9s%10.2e   %-6s%10.2e   %-8s%10.2e\n",
						   x[i]->description, "old MH2O",
						   (double) mass_water_aq_x, "new MH2O",
						   (double) (mass_water_aq_x * d), "delta",
						   (double) delta[i], "10**d/c", (double) d));
			}
			mass_water_aq_x *= d;

			mass_water_bulk_x = mass_water_aq_x + mass_water_surfaces_x;
			if (debug_model == TRUE && dl_type_x != cxxSurface::NO_DL)
			{
				output_msg(sformatf(
						   "mass_water bulk: %e\taq: %e\tsurfaces: %e\n",
						   (double) mass_water_bulk_x,
						   (double) mass_water_aq_x,
						   (double) mass_water_surfaces_x));
			}
			x[i]->master[0]->s->moles = mass_water_aq_x / gfw_water;
/*appt */
			if (use.Get_surface_ptr() != NULL)
			{
				if (use.Get_surface_ptr()->Get_debye_lengths() > 0)
					x[i]->master[0]->s->moles = mass_water_bulk_x / gfw_water;
			}

			if (mass_water_aq_x < 1e-10)
			{
				error_string = sformatf(
						"Mass of water is less than 1e-10 kilogram.\n"
						"The aqueous phase may not be stable relative to given masses of minerals.");
				warning_msg(error_string);
				stop_program = TRUE;
				return (TRUE);
			}
/*   Pure phases */
		}
		else if (x[i]->type == PP)
		{
			//comp_ptr = pp_assemblage_ptr->Find(x[i]->pp_assemblage_comp_name);
			comp_ptr = (cxxPPassemblageComp *) x[i]->pp_assemblage_comp_ptr;
			/*if (fabs(delta[i]) > epsilon) converge=FALSE; */
			if (debug_model == TRUE)
			{
				output_msg(sformatf(
						   "%-10.10s %-9s%10.2e   %-9s%10.2e   %-6s%10.2e\n",
						   x[i]->description, "old mass",
						   (double) x[i]->moles, "new mass",
						   (double) (x[i]->moles - delta[i]), "delta",
						   (double) delta[i]));
			}
			if (equal(x[i]->moles, delta[i], ineq_tol))
			{
				x[i]->moles = 0.0;
			}
			else
			{
				x[i]->moles -= delta[i];
			}

			if (x[i]->dissolve_only == TRUE)
			{
				if (equal
					(x[i]->moles, comp_ptr->Get_initial_moles(), ineq_tol))
					x[i]->moles = comp_ptr->Get_initial_moles();
			}
			/*if (fabs(x[i]->moles) < MIN_RELATED_SURFACE) x[i]->moles = 0.0; */
		}
		else if (x[i]->type == GAS_MOLES)
		{
			cxxGasPhase *gas_phase_ptr = use.Get_gas_phase_ptr();
			/*if (fabs(delta[i]) > epsilon) converge=FALSE; */
			/*if (gas_in == TRUE && fabs(residual[i]) > epsilon) converge=FALSE; */
			if (debug_model == TRUE)
			{
				output_msg(sformatf(
						   "%-10.10s %-9s%10.2e   %-9s%10.2e   %-6s%10.2e\n",
						   x[i]->description, "old mol",
						   (double) x[i]->moles, "new mol",
						   (double) (x[i]->moles + delta[i]), "delta",
						   (double) delta[i]));
			}
			x[i]->moles += delta[i];
			if (x[i]->moles < MIN_TOTAL)
				x[i]->moles = MIN_TOTAL;
			if (x[i] == gas_unknown && gas_phase_ptr->Get_type() == cxxGasPhase::GP_VOLUME && 
				 !calculating_deriv)
			{					
				if (debug_model == TRUE)
				{
					output_msg(sformatf(
								"%-10.10s %-9s%10.2e   %-9s%10.2e   %-6s%10.2e\n",
								"Pressure", "old P",
								(double) last_patm_x, "new P",
								(double) gas_phase_ptr->Get_total_p(), "iter P",
								(double) patm_x));
				}
				patm_x = gas_phase_ptr->Get_total_p();
				//if (patm_x < 1e-10 && patm_x < p_sat)
				//{
				//	patm_x = ( 1 * patm_x + p_sat) / 2.0;
				//}
				if (patm_x > 1500)
				  patm_x = 1500;
			}
			last_patm_x = patm_x;
		}
		else if (x[i]->type == SS_MOLES)
		{

			/*if (fabs(delta[i]) > epsilon) converge=FALSE; */
			/*if (x[i]->ss_in == TRUE && fabs(residual[i]) > epsilon) converge=FALSE; */
			if (debug_model == TRUE)
			{
				output_msg(sformatf(
						   "%-10.10s %-9s%10.2e   %-9s%10.2e   %-6s%10.2e\n",
						   x[i]->description, "old mol",
						   (double) x[i]->moles, "new mol",
						   (double) (x[i]->moles - delta[i]), "delta",
						   (double) delta[i]));
			}
			x[i]->moles -= delta[i];
			if (x[i]->moles < MIN_TOTAL_SS && calculating_deriv == FALSE)
				x[i]->moles = MIN_TOTAL_SS;
			cxxSScomp *comp_ptr = (cxxSScomp *) x[i]->ss_comp_ptr;
			comp_ptr->Set_moles(x[i]->moles);
/*   Pitzer gamma */
		}
		else if (x[i]->type == PITZER_GAMMA)
		{
			if (full_pitzer == FALSE)
				continue;
			d = delta[i];
			if (debug_model == TRUE)
			{
				output_msg(sformatf(
						   "%-10.10s %-9s%10.5f   %-9s%10.5f   %-6s%10.2e   %-8s%10.2e\n",
						   x[i]->description, "old lg", (double) x[i]->s->lg,
						   "new lg", (double) (x[i]->s->lg + d), "delta",
						   (double) delta[i], "delta", (double) d));
			}
			x[i]->s->lg += d;
		}
	}
/*
 *   Reset total molalities in mass balance equations
 */
	if (pure_phase_unknown != NULL || gas_unknown != NULL
		|| ss_unknown != NULL)
	{
		for (i = 0; i < count_unknowns; i++)
		{
			if (x[i]->type == MB || x[i]->type == MH ||
				x[i]->type == MH2O ||
				x[i]->type == CB || x[i]->type == EXCH
				|| x[i]->type == SURFACE)
			{
				/*if (fabs(x[i]->delta) > epsilon*x[i]->moles) converge = FALSE; */
				if (x[i]->type == SURFACE)
					x[i]->delta = 0.0;
				if (debug_model == TRUE)
				{
					output_msg(sformatf(
							   "%-10.10s %-9s%10.2e   %-9s%10.2e   %-6s%10.2e\n",
							   x[i]->description, "old mole",
							   (double) x[i]->moles, "new mole",
							   (double) (x[i]->moles + x[i]->delta), "delta",
							   (double) x[i]->delta));
				}
				x[i]->moles += x[i]->delta;
			}
		}
	}
	converge = FALSE;
	return (converge);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
residuals(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Calculates residuals for all equations
 */
	int i, j;
	int converge;

	LDBLE l_toler;
	LDBLE sum_residual;
	LDBLE sinh_constant;
	LDBLE sum, sum1;
	struct master *master_ptr, *master_ptr1, *master_ptr2;
	LDBLE sigmaddl, negfpsirt;
	int print_fail;
	std::vector<LDBLE> cd_psi;
	print_fail = FALSE;
	sum_residual = 0.0;
	sigmaddl = 0;
	sum = 0;
/*
 *   Calculate residuals
 */
	converge = TRUE;
	l_toler = convergence_tolerance;

	for (i = 0; i < count_unknowns; i++)
	{
		if (x[i]->type == MB)
		{
			residual[i] = x[i]->moles - x[i]->f;
			if ((fabs(residual[i]) > l_toler * x[i]->moles && 
				fabs(residual[i]) > sqrt(fabs(x[i]->moles) * MIN_TOTAL) && 
				x[i]->moles > MIN_TOTAL) 
				|| 
				x[i]->moles < 0)
			{
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual %d: %s %d %e\n", iterations,
							   x[i]->description, i, residual[i]));
				converge = FALSE;
			}
		}
		else if (x[i]->type == ALK)
		{
			residual[i] = x[i]->moles - x[i]->f;
			if (fabs(residual[i]) > l_toler * x[i]->moles)
			{
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual %d: %s %d %e\n", iterations,
							   x[i]->description, i, residual[i]));
				converge = FALSE;
			}
		}
		else if (x[i]->type == SOLUTION_PHASE_BOUNDARY)
		{
			residual[i] = x[i]->f * LOG_10;
			if (fabs(residual[i]) > l_toler)
			{
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual %d: %s %d %e\n", iterations,
							   x[i]->description, i, residual[i]));
				converge = FALSE;
			}
		}
		else if (x[i]->type == CB)
		{
			residual[i] = -x[i]->f;
			if (ph_unknown == charge_balance_unknown)
			{
				residual[i] += x[i]->moles;
			}
			if (fabs(residual[i]) >= l_toler * mu_x * mass_water_aq_x)
			{
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual %d: %s %d %e\n", iterations,
							   x[i]->description, i, residual[i]));
				converge = FALSE;
			}
		}
		else if (x[i]->type == MU /*&& pitzer_model == FALSE && sit_model == FALSE*/)
		{
			// Using straight ionic strength equation
			{
				residual[i] = mass_water_aq_x * mu_x - 0.5 * x[i]->f;
			}
			if (fabs(residual[i]) > l_toler * mu_x * mass_water_aq_x)
			{
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual %d: %s %d %e\n", iterations,
							   x[i]->description, i, residual[i]));
				converge = FALSE;
			}
		}
		else if (x[i]->type == AH2O)
		{
			if (dampen_ah2o)
			{
				// a = 0.017; Y = sum(m(i)); X = Mw (mass of water)
				// Aw = (1 - a y/x) 0.5 (tanh(100 (0.95 - ay/x)) + 1) + 0.05 (0.5 (1 - tanh(100(0.95 - ay/x))))
				//residual[i] =  exp(s_h2o->la * LOG_10) - ((1.0 - 0.017 * x[i]->f/mass_water_aq_x) *
				//	0.5*(tanh(100.0*(0.95 - 0.017*x[i]->f/mass_water_aq_x)) + 1) +
				//	0.05*0.5*(1.0 - tanh(100.0*(0.95 - 0.017*x[i]->f/mass_water_aq_x)))) ;
				residual[i] =  exp(s_h2o->la * LOG_10) - ((1.0 - AH2O_FACTOR * x[i]->f/mass_water_aq_x) *
					0.5*(tanh(100.0*(0.95 - AH2O_FACTOR*x[i]->f/mass_water_aq_x)) + 1) +
					0.05*0.5*(1.0 - tanh(100.0*(0.95 - AH2O_FACTOR*x[i]->f/mass_water_aq_x)))) ;
			}
			else
			{
				//residual[i] = mass_water_aq_x * exp(s_h2o->la * LOG_10) - mass_water_aq_x +
				//	0.017 * x[i]->f;
				residual[i] = mass_water_aq_x * exp(s_h2o->la * LOG_10) - mass_water_aq_x +
					AH2O_FACTOR * x[i]->f;
			}

			if (pitzer_model || sit_model)
			{
				residual[i] = pow((LDBLE) 10.0, s_h2o->la) - AW;
				if (full_pitzer == FALSE)
				{
					residual[i] = 0.0;
				}
			}
			if (fabs(residual[i]) > l_toler)
			{
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual %d: %s %d %e\n", iterations,
							   x[i]->description, i, residual[i]));
				converge = FALSE;
			}
		}
		else if (x[i]->type == MH
				 && (pitzer_model == FALSE || pitzer_pe == TRUE))
		{
#ifdef COMBINE
			residual[i] = x[i]->moles - x[i]->f;
#else
			residual[i] = (x[i]->moles - 2 * s_h2o->moles) - x[i]->f;
			x[i]->f += 2 * s_h2o->moles;
#endif
			if (mass_water_switch == TRUE)
			{
				residual[i] -=
					2 * (mass_oxygen_unknown->moles - mass_oxygen_unknown->f);
			}
#ifdef COMBINE
#ifndef COMBINE_CHARGE
			if (fabs(residual[i]) >
				l_toler * (x[i]->moles + 2 * mass_oxygen_unknown->moles))
			{
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual %d: %s %d %e\n", iterations,
							   x[i]->description, i, residual[i]));
				converge = FALSE;
			}
#else
			if (fabs(residual[i]) >
				l_toler * (x[i]->moles + 2 * mass_oxygen_unknown->moles +
						 charge_balance_unknown->moles))
			{
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual %d: %s %d %e\n", iterations,
							   x[i]->description, i, residual[i]));
				converge = FALSE;
			}
#endif
#else
			if (fabs(residual[i]) > l_toler * x[i]->moles)
			{
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual %d: %s %d %e\n", iterations,
							   x[i]->description, i, residual[i]));
				converge = FALSE;
			}
#endif
		}
		else if (x[i]->type == MH2O)
		{
			if (mass_water_switch == TRUE)
				continue;
#ifdef COMBINE
			residual[i] = x[i]->moles - x[i]->f;
#else
			residual[i] = (x[i]->moles - s_h2o->moles) - x[i]->f;
			x[i]->f += s_h2o->moles;
#endif
			if (fabs(residual[i]) > 0.01 * l_toler * x[i]->moles)
			{
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual %d: %s %d %e\n", iterations,
							   x[i]->description, i, residual[i]));
				converge = FALSE;
			}
		}
		else if (x[i]->type == PP)
		{
			//cxxPPassemblageComp * comp_ptr = pp_assemblage_ptr->Find(x[i]->pp_assemblage_comp_name);
			cxxPPassemblageComp * comp_ptr = (cxxPPassemblageComp *) x[i]->pp_assemblage_comp_ptr;
			residual[i] = x[i]->f * LOG_10;
			if (comp_ptr->Get_add_formula().size() == 0)
			{
				if (x[i]->dissolve_only == TRUE)
				{
					if ((residual[i] > l_toler && x[i]->moles > 0.0)
						|| (residual[i] < -l_toler
							&& (comp_ptr->Get_initial_moles() -
								x[i]->moles) > 0))
					{
						if (print_fail)
							output_msg(sformatf(
									   "Failed Residual %d: %s %d %e\n",
									   iterations, x[i]->description, i,
									   residual[i]));
						converge = FALSE;
					}
				}
				else
				{
					if (residual[i] < -l_toler || iterations < 1)
					{
						if (print_fail)
							output_msg(sformatf(
									   "Failed Residual %d: %s %d %e\n",
									   iterations, x[i]->description, i,
									   residual[i]));
						converge = FALSE;
					}
				}
			}
			else
			{
				/* if (x[i]->moles > 0.0 && fabs(residual[i]) > l_toler) converge = FALSE; */
				if (residual[i] < -l_toler || iterations < 1)
				{
					if (print_fail)
						output_msg(sformatf(
								   "Failed Residual %d: %s %d %e\n",
								   iterations, x[i]->description, i,
								   residual[i]));
					converge = FALSE;
				}
			}
		}
		else if (x[i]->type == GAS_MOLES)
		{
			cxxGasPhase *gas_phase_ptr = use.Get_gas_phase_ptr();
			if (gas_phase_ptr->Get_type() == cxxGasPhase::GP_VOLUME && 
				(gas_phase_ptr->Get_pr_in() || force_numerical_fixed_volume) && numerical_fixed_volume)
			{
				residual[i] = x[i]->moles - x[i]->phase->moles_x;
			}
			else
			{
				residual[i] = gas_phase_ptr->Get_total_p() - x[i]->f;
			}

			if (fabs(residual[i]) > l_toler && gas_in == TRUE)
			{
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual %d: %s %d %e\n", iterations,
							   x[i]->description, i, residual[i]));
				converge = FALSE;
			}
			if (gas_phase_ptr->Get_type() == cxxGasPhase::GP_VOLUME &&
				(fabs(last_patm_x - patm_x) > 0.001 || fabs(last_patm_x - gas_phase_ptr->Get_total_p()) > 0.001)
				&& !calculating_deriv)
			{
				if (print_fail)
					output_msg(sformatf("Failed pressure test %d: %e %e %e\n", iterations, last_patm_x, patm_x, 
						gas_phase_ptr->Get_total_p()));
				converge = FALSE;
			}
		}
		else if (x[i]->type == SS_MOLES)
		{
			residual[i] = x[i]->f * LOG_10;
			if (fabs(residual[i]) > l_toler && x[i]->ss_in == TRUE)
			{
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual %d: %s %d %e\n", iterations,
							   x[i]->description, i, residual[i]));
				converge = FALSE;
			}
		}
		else if (x[i]->type == EXCH)
		{
			residual[i] = x[i]->moles - x[i]->f;
			if (x[i]->moles <= MIN_RELATED_SURFACE)
			{
				if (fabs(residual[i]) > l_toler)
				{
					if (print_fail)
						output_msg(sformatf(
								   "Failed Residual %d: %s %d %e\n",
								   iterations, x[i]->description, i,
								   residual[i]));
					converge = FALSE;
				}
			}
			else if (fabs(residual[i]) > l_toler * x[i]->moles)
			{
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual %d: %s %d %e\n", iterations,
							   x[i]->description, i, residual[i]));
				converge = FALSE;
			}
		}
		else if (x[i]->type == SURFACE)
		{
			residual[i] = x[i]->moles - x[i]->f;
			if (x[i]->moles <= MIN_RELATED_SURFACE)
			{
				if (fabs(residual[i]) > l_toler)
				{
					if (print_fail)
						output_msg(sformatf(
								   "Failed Residual %d: %s %d %e\n",
								   iterations, x[i]->description, i,
								   residual[i]));
					converge = FALSE;
				}
			}
			else if (fabs(residual[i]) < ineq_tol && fabs(residual[i]) < 1e-2*x[i]->moles)
			{
			}
			else if (fabs(residual[i]) > l_toler * x[i]->moles)
			{
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual %d: %s %d %e\n", iterations,
							   x[i]->description, i, residual[i]));
				converge = FALSE;
			}
		}
		else if (x[i]->type == PITZER_GAMMA)
		{
			if (full_pitzer == FALSE)
				continue;
			residual[i] = x[i]->s->lg - x[i]->s->lg_pitzer;
			if (fabs(residual[i]) > l_toler)
			{
				/*
				   fprintf(stderr,"Residuals %d: %s %d %e\n", iterations, x[i]->description, i, residual[i]);
				 */
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual %d: %s %d %e\n", iterations,
							   x[i]->description, i, residual[i]));
				converge = FALSE;
			}
		}
		else if (x[i]->type == SURFACE_CB && use.Get_surface_ptr()->Get_type() == cxxSurface::DDL)
		{
			/*sinh_constant = 0.1174; */
			sinh_constant =
				//sqrt(8 * EPSILON * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000) *
				//	 tk_x * 1000);
				sqrt(8 * eps_r * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000) *
					 tk_x * 1000);
/*			if (x[i]->surface_charge->grams <= MIN_RELATED_SURFACE) { */
			cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[i]->surface_charge);
			if (charge_ptr->Get_grams() == 0)
			{
				residual[i] = 0.0;
			}
			else if (dl_type_x != cxxSurface::NO_DL)
			{
				residual[i] = -x[i]->f;
			}
			else
			{
/*
 *   sinh_constant is (8 e e0 R T 1000)**1/2
 *				 = sqrt(8*EPSILON*EPSILON_ZERO*(R_KJ_DEG_MOL*1000)*t_x*1000)
 *				 ~ 0.1174 at 25C
 */
				residual[i] =
					sinh_constant * sqrt(mu_x) *
					sinh(x[i]->master[0]->s->la * LOG_10) -
					x[i]->f * F_C_MOL / (charge_ptr->Get_specific_area() *
										 charge_ptr->Get_grams());
			}
			if (debug_model == TRUE)
			{
				output_msg(sformatf( "Charge/Potential\n"));
				if (charge_ptr->Get_grams() > 0)
				{
					output_msg(sformatf(
							   "\tSum of surface charge %e eq\n",
							   (double) (x[i]->f
										 /* F_C_MOL / (x[i]->surface_charge->specific_area * x[i]->surface_charge->grams) */
							   )));
				}
				else
				{
					output_msg(sformatf( "\tResidual %e\n",
							   (double) x[i]->f));
				}
				output_msg(sformatf( "\t				grams %g\n",
						   (double) charge_ptr->Get_grams()));
				output_msg(sformatf( "\tCharge from potential %e eq\n",
						   (double) (charge_ptr->Get_specific_area() *
									 charge_ptr->Get_grams() / F_C_MOL *
									 sinh_constant * sqrt(mu_x) *
									 sinh(x[i]->master[0]->s->la * LOG_10))));
				output_msg(sformatf( "\t			 FPsi/2RT %e\n",
						   (double) (x[i]->master[0]->s->la * LOG_10)));
				output_msg(sformatf( "\t	   Sinh(FPsi/2RT) %e\n",
						   sinh(x[i]->master[0]->s->la * LOG_10)));
				output_msg(sformatf( "\t	   Cosh(FPsi/2RT) %e\n",
						   cosh(x[i]->master[0]->s->la * LOG_10)));
				output_msg(sformatf( "\t		   Sqrt(mu_x) %e\n",
						   sqrt(mu_x)));
			}
			if (charge_ptr->Get_grams() > MIN_RELATED_SURFACE
				&& fabs(residual[i]) > l_toler)
			{
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual %d: %s %d %e\n", iterations,
							   x[i]->description, i, residual[i]));
				converge = FALSE;
			}
		}
		else if (x[i]->type == SURFACE_CB
				 && use.Get_surface_ptr()->Get_type() == cxxSurface::CD_MUSIC)
		{
			cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[i]->surface_charge);
			if (charge_ptr->Get_grams() == 0)
			{
				residual[i] = 0.0;
				cd_psi.clear();
				cd_psi.push_back(0.0);
				cd_psi.push_back(0.0);
				cd_psi.push_back(0.0);
			}
			else
			{
				/* sum is in moles of charge */
				/*psi = pow(10, x[i]->surface_charge->psi_master->s->la); */ /* = exp(-Fpsi/RT) */
				master_ptr =
					surface_get_psi_master(charge_ptr->Get_name().c_str(),
										   SURF_PSI);
				master_ptr1 =
					surface_get_psi_master(charge_ptr->Get_name().c_str(),
										   SURF_PSI1);
				master_ptr2 =
					surface_get_psi_master(charge_ptr->Get_name().c_str(),
										   SURF_PSI2);
				
				cd_psi.clear();
				
				cd_psi.push_back(-(master_ptr->s->la * LOG_10) * R_KJ_DEG_MOL * tk_x /
					F_KJ_V_EQ);
				cd_psi.push_back(-(master_ptr1->s->la * LOG_10) * R_KJ_DEG_MOL * tk_x /
					F_KJ_V_EQ);
				cd_psi.push_back(-(master_ptr2->s->la * LOG_10) * R_KJ_DEG_MOL * tk_x /
					F_KJ_V_EQ);
				sum = 0;
				for (j = 0; j < x[i]->count_comp_unknowns; j++)
				{
					sum +=
						x[i]->comp_unknowns[j]->moles *
						x[i]->comp_unknowns[j]->master[0]->s->z;
				}
				charge_ptr->Set_sigma0(
					(x[i]->f +
					 sum) * F_C_MOL / (charge_ptr->Get_specific_area() *
									   charge_ptr->Get_grams()));
				/* f is in moles */
				/* eqns A-3 */
				residual[i] =
					charge_ptr->Get_sigma0() -
					charge_ptr->Get_capacitance0() *
					(cd_psi[0] - cd_psi[1]);
			}
			if (charge_ptr->Get_grams() > MIN_RELATED_SURFACE
				&& fabs(residual[i]) > l_toler)
			{
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual A %d: %s %d %e\n",
							   iterations, x[i]->description, i, residual[i]));
				converge = FALSE;
			}
		}		
		else if (x[i]->type == SURFACE_CB && use.Get_surface_ptr()->Get_type() == cxxSurface::CCM)
		{
			cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[i]->surface_charge);
			if (charge_ptr->Get_grams() == 0)
			{
				residual[i] = 0.0;
			}
			else if (dl_type_x != cxxSurface::NO_DL)
			{
				residual[i] = -x[i]->f;
			}
			else
			{
				residual[i] =
					charge_ptr->Get_capacitance0() * x[i]->master[0]->s->la * 2 * R_KJ_DEG_MOL *
								 tk_x * LOG_10 / F_KJ_V_EQ -
					x[i]->f * F_C_MOL / (charge_ptr->Get_specific_area() *
										 charge_ptr->Get_grams());
			}
			if (debug_model == TRUE)
			{
				output_msg(sformatf( "Charge/Potential\n"));
				if (charge_ptr->Get_grams() > 0)
				{
					output_msg(sformatf(
							   "\tSum of surface charge %e eq\n",
							   (double) (x[i]->f
										 /* F_C_MOL / (x[i]->surface_charge->specific_area * x[i]->surface_charge->grams) */
							   )));
				}
				else
				{
					output_msg(sformatf( "\tResidual %e\n",
							   (double) x[i]->f));
				}
				output_msg(sformatf( "\t				grams %g\n",
						   (double) charge_ptr->Get_grams()));
				output_msg(sformatf( "\tCharge from potential %e eq\n",
						   (double) (charge_ptr->Get_capacitance0() * x[i]->master[0]->s->la * 2 * R_KJ_DEG_MOL *
								 tk_x * LOG_10 / F_KJ_V_EQ)));
				output_msg(sformatf( "\t			      Psi %e\n",
						   (double) (x[i]->master[0]->s->la * 2 * R_KJ_DEG_MOL *
								 tk_x * LOG_10 / F_KJ_V_EQ)));
			}
			if (charge_ptr->Get_grams() > MIN_RELATED_SURFACE
				&& fabs(residual[i]) > l_toler)
			{
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual %d: %s %d %e\n", iterations,
							   x[i]->description, i, residual[i]));
				converge = FALSE;
			}
		}
		else if (x[i]->type == SURFACE_CB1)
		{
			cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[i]->surface_charge);
			if (charge_ptr->Get_grams() == 0)
			{
				residual[i] = 0.0;
			}
			else
			{
				/* eqns A-4 */
				/*psi = pow(10, x[i]->surface_charge->psi_master1->s->la); */ /* = exp(-Fpsi/RT) */
				charge_ptr->Set_sigma1(
					x[i]->f * F_C_MOL / (charge_ptr->Get_specific_area() *
										 charge_ptr->Get_grams()));
				residual[i] =
					(charge_ptr->Get_sigma0() +
					 charge_ptr->Get_sigma1()) -
					charge_ptr->Get_capacitance1() *
					(cd_psi[1] - cd_psi[2]);
			}
			if (charge_ptr->Get_grams() > MIN_RELATED_SURFACE
				&& fabs(residual[i]) > l_toler)
			{
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual B %d: %s %d %e\n",
							   iterations, x[i]->description, i, residual[i]));
				converge = FALSE;
			}
		}
		else if (x[i]->type == SURFACE_CB2)
		{
			cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[i]->surface_charge);
			if (charge_ptr->Get_grams() == 0)
			{
				residual[i] = 0.0;
			}
			else if (dl_type_x != cxxSurface::NO_DL)
			{
				sum = 0;
				sum1 = 0;
				for (j = 0; j < count_s_x; j++)
				{
					if (s_x[j]->type == SURF)
					{
						sum += under(s_x[j]->lm) * s_x[j]->dz[2];
					}
					if (s_x[j]->type < H2O)
					{
						int is = s_x[j]->number;
						sum1 += s_x[j]->z * s_x[j]->z * s_diff_layer[is][charge_ptr->Get_name()].Get_g_moles();
					}
				}
				charge_ptr->Set_sigma2(
					sum * F_C_MOL / (charge_ptr->Get_specific_area() *
									 charge_ptr->Get_grams()));
				charge_ptr->Set_sigmaddl(
					(x[i]->f -
					 sum) * F_C_MOL / (charge_ptr->Get_specific_area() *
									   charge_ptr->Get_grams()));

				residual[i] =
					x[i]->f + (charge_ptr->Get_sigma0() +
							   charge_ptr->Get_sigma1()) *
					(charge_ptr->Get_specific_area() *
					 charge_ptr->Get_grams()) / F_C_MOL;
				/* residual[i] = sum + (x[i]->surface_charge->sigma0 + x[i]->surface_charge->sigma1) * (x[i]->surface_charge->specific_area * x[i]->surface_charge->grams) / F_C_MOL */
			}
			else
			{
				/* eqns A-6 and A-7 */
				sinh_constant =
					//sqrt(8 * EPSILON * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000) *
					//	 tk_x * 1000);
					sqrt(8 * eps_r * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000) *
						 tk_x * 1000);
				/*
				 *   sinh_constant is (8 e e0 R T 1000)**1/2
				 *       = sqrt(8*EPSILON*EPSILON_ZERO*(R_KJ_DEG_MOL*1000)*t_x*1000)
				 *       ~ 0.1174 at 25C
				 */
				master_ptr2 =
					surface_get_psi_master(charge_ptr->Get_name().c_str(),
										   SURF_PSI2);
				negfpsirt = master_ptr2->s->la * LOG_10;
				sum = 0;
				sum1 = 0;
				for (j = 0; j < count_s_x; j++)
				{
					if (s_x[j]->type < H2O)
					{
						sum +=
							under(s_x[j]->lm) *
							(exp(s_x[j]->z * negfpsirt) - 1);
						sum1 += under(s_x[j]->lm) * s_x[j]->z;
					}
				}

				/* add fictitious monovalent ion that balances charge */
				//sum += fabs(sum1) * (exp(-sum1 / fabs(sum1) * negfpsirt) - 1);
				if (sum1 >= 0) 
				{
					sum += fabs(sum1) * (exp(-negfpsirt) - 1);
				}
				else
				{
					sum += fabs(sum1) * (exp(negfpsirt) - 1);
				}

				if (sum < 0)
				{
					sum = -sum;
					{
						if (print_fail)
							output_msg(sformatf(
									   "Failed Residual C %d: %s %d %e %e\n",
									   iterations, x[i]->description, i, sum,
									   l_toler));
						converge = FALSE;
					}
					/*output_msg(sformatf( "Negative sum, iteration %d\n", iterations)); */
				}
				charge_ptr->Set_sigma2(
					x[i]->f * F_C_MOL / (charge_ptr->Get_specific_area() *
										 charge_ptr->Get_grams()));
				if ((negfpsirt) < 0)
				{
					sigmaddl = -0.5 * sinh_constant * sqrt(sum);
				}
				else
				{
					sigmaddl = 0.5 * sinh_constant * sqrt(sum);
				}
				charge_ptr->Set_sigmaddl(sigmaddl);
				residual[i] =
					(charge_ptr->Get_sigma0() +
					 charge_ptr->Get_sigma1() +
					 charge_ptr->Get_sigma2()) + sigmaddl;
			}
			if (debug_model == TRUE)
			{
				master_ptr =
					surface_get_psi_master(charge_ptr->Get_name().c_str(),
										   SURF_PSI);
				master_ptr1 =
					surface_get_psi_master(charge_ptr->Get_name().c_str(),
										   SURF_PSI1);
				master_ptr2 =
					surface_get_psi_master(charge_ptr->Get_name().c_str(),
										   SURF_PSI2);
				output_msg(sformatf( "CD_music Charge/Potential 2\n"));
				output_msg(sformatf( "\tgrams	      %g\n",
						   (double) charge_ptr->Get_grams()));
				output_msg(sformatf( "\tCapacitances       %g\t%g\n",
						   (double) charge_ptr->Get_capacitance0(),
						   charge_ptr->Get_capacitance1()));
				output_msg(sformatf( "\t-F/(RT)	    %g\n",
						   (double) -F_KJ_V_EQ / (R_KJ_DEG_MOL * tk_x)));
				output_msg(sformatf( "\tResidual 0	 %14e\n",
						   (double) residual[master_ptr->unknown->number]));
				output_msg(sformatf( "\tResidual 1	 %14e\n",
						   (double) residual[master_ptr1->unknown->number]));
				output_msg(sformatf( "\tResidual 2	 %14e\n",
						   (double) residual[master_ptr2->unknown->number]));
				output_msg(sformatf( "\texp(-FPsi0/RT)     %14e",
						   (double) pow((LDBLE) 10., master_ptr->s->la)));
				output_msg(sformatf( "\tPsi0	       %14e\n",
						   (double) cd_psi[0]));
				output_msg(sformatf( "\texp(-FPsi1/RT)     %14e",
						   (double) pow((LDBLE) 10., master_ptr1->s->la)));
				output_msg(sformatf( "\tPsi1	       %14e\n",
						   (double) cd_psi[1]));
				output_msg(sformatf( "\texp(-FPsi2/RT)     %14e",
						   (double) pow((LDBLE) 10., master_ptr2->s->la)));
				output_msg(sformatf( "\tPsi2	       %14e\n",
						   (double) cd_psi[2]));
				output_msg(sformatf( "\tf 0		%14e",
						   (double) master_ptr->unknown->f));
				output_msg(sformatf( "\tsigma 0	    %14e\n",
						   (double) charge_ptr->Get_sigma0()));
				output_msg(sformatf( "\tf 1		%14e",
						   (double) master_ptr1->unknown->f));
				output_msg(sformatf( "\tsigma 1	    %14e\n",
						   (double) charge_ptr->Get_sigma1()));
				output_msg(sformatf( "\tf 2		%14e",
						   (double) master_ptr2->unknown->f));
				output_msg(sformatf( "\tsigma 2	    %14e\n",
						   (double) charge_ptr->Get_sigma2()));
				output_msg(sformatf( "\tsigma ddl	  %14e\n",
						   (double) sigmaddl));
				output_msg(sformatf( "\texp sum	    %14e\n",
						   (double) sum));

			}
			if (charge_ptr->Get_grams() > MIN_RELATED_SURFACE
				&& fabs(residual[i]) > l_toler)
			{
				if (print_fail)
					output_msg(sformatf(
							   "Failed Residual D %d: %s %d %e\n",
							   iterations, x[i]->description, i, residual[i]));
				converge = FALSE;
			}
		}
/*
 *   Store residuals in array
 */
		array[(i + 1) * (count_unknowns + 1) - 1] = residual[i];
		sum_residual += fabs(residual[i]);
	}
/*
 *   Return
 */
	if ((pitzer_model == TRUE || sit_model == TRUE) && iterations < 1)
		return (OK);
	if (converge == TRUE)
	{
		return (CONVERGED);
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
set(int initial)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sets initial guesses for unknowns if initial == TRUE
 *   Revises guesses whether initial is true or not
 */
	int i;
	cxxSolution *solution_ptr;
/*
 *   Set initial log concentrations to zero
 */
	if (pitzer_model == TRUE)
		return (set_pz(initial));
	if (sit_model == TRUE)
		return (set_sit(initial));
	iterations = -1;
	solution_ptr = use.Get_solution_ptr();
	for (i = 0; i < count_s_x; i++)
	{
		s_x[i]->lm = LOG_ZERO_MOLALITY;
		s_x[i]->lg = 0.0;
	}
/*
 *   Set master species activities
 */

	tc_x = solution_ptr->Get_tc();
	tk_x = tc_x + 273.15;

	patm_x = solution_ptr->Get_patm();  // done in calc_rho_0(tc, pa)

/*
 *   H+, e-, H2O
 */
	mass_water_aq_x = solution_ptr->Get_mass_water();
	mu_x = solution_ptr->Get_mu();
	s_h2o->moles = mass_water_aq_x / gfw_water;
	s_h2o->la = log10(solution_ptr->Get_ah2o());
	s_hplus->la = -solution_ptr->Get_ph();
	s_hplus->lm = s_hplus->la;
	s_hplus->moles = exp(s_hplus->lm * LOG_10) * mass_water_aq_x;
	s_eminus->la = -solution_ptr->Get_pe();
	if (initial == TRUE)
		initial_guesses();
	if (dl_type_x != cxxSurface::NO_DL)
		initial_surface_water();
	revise_guesses();
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
initial_guesses(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Make initial guesses for activities of master species and
 *   ionic strength
 */
	int i;
	cxxSolution *solution_ptr;

	solution_ptr = use.Get_solution_ptr();
	mu_x =
		s_hplus->moles +
		exp((solution_ptr->Get_ph() - 14.) * LOG_10) * mass_water_aq_x;
	mu_x /= mass_water_aq_x;
	s_h2o->la = 0.0;
	for (i = 0; i < count_unknowns; i++)
	{
		if (x[i] == ph_unknown || x[i] == pe_unknown)
			continue;
		if (x[i]->type < CB)
		{
			mu_x +=
				x[i]->moles / mass_water_aq_x * 0.5 * x[i]->master[0]->s->z *
				x[i]->master[0]->s->z;
			x[i]->master[0]->s->la = log10(x[i]->moles / mass_water_aq_x);
		}
		else if (x[i]->type == CB)
		{
			x[i]->master[0]->s->la =
				log10(0.001 * x[i]->moles / mass_water_aq_x);
		}
		else if (x[i]->type == SOLUTION_PHASE_BOUNDARY)
		{
			x[i]->master[0]->s->la =
				log10(0.001 * x[i]->moles / mass_water_aq_x);
		}
		else if (x[i]->type == EXCH)
		{
			if (x[i]->moles <= 0)
			{
				x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
			}
			else
			{
				x[i]->master[0]->s->la = log10(x[i]->moles);
			}
		}
		else if (x[i]->type == SURFACE)
		{
			if (x[i]->moles <= 0)
			{
				x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
			}
			else
			{
				x[i]->master[0]->s->la = log10(0.1 * x[i]->moles);
			}
		}
		else if (x[i]->type == SURFACE_CB)
		{
			x[i]->master[0]->s->la = 0.0;
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
revise_guesses(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Revise la's of master species
 */
	int i;
	int l_iter, max_iter, repeat, fail;
	LDBLE weight, f;
	LDBLE d;

	max_iter = 10;
	gammas(mu_x);
	l_iter = 0;
	repeat = TRUE;
	fail = FALSE;;
	while (repeat == TRUE)
	{
		l_iter++;
		if (debug_set == TRUE)
		{
			output_msg(sformatf( "\nBeginning set iteration %d.\n",
					   l_iter));
		}
		if (l_iter == max_iter + 1)
		{
			log_msg(sformatf(
					   "Did not converge in set, iteration %d.\n",
					   iterations));
			fail = TRUE;
		}
		if (l_iter > 2 * max_iter)
		{
			log_msg(sformatf(
					   "Did not converge with relaxed criteria in set.\n"));
			return (OK);
		}
		molalities(TRUE);
		mb_sums();
		if (state < REACTION)
		{
			sum_species();
		}
		else
		{
			for (i = 0; i < count_unknowns; i++)
			{
				x[i]->sum = x[i]->f;
			}
		}
/* debug
		if (debug_set == TRUE) {
			pr.species = TRUE;
			pr.all = TRUE;
			print_species();
		}
 */
		repeat = FALSE;
		for (i = 0; i < count_unknowns; i++)
		{
			if (x[i] == ph_unknown || x[i] == pe_unknown)
				continue;
			if (x[i]->type == MB ||
/*				x[i]->type == ALK || */
				x[i]->type == CB ||
				x[i]->type == SOLUTION_PHASE_BOUNDARY ||
				x[i]->type == EXCH || x[i]->type == SURFACE)
			{

				if (debug_set == TRUE)
				{
					output_msg(sformatf(
							   "\n\t%5s  at beginning of set %d: %e\t%e\t%e\n",
							   x[i]->description, l_iter, (double) x[i]->sum,
							   (double) x[i]->moles,
							   (double) x[i]->master[0]->s->la));
				}
				if (fabs(x[i]->moles) < 1e-30)
					x[i]->moles = 0;

				f = fabs(x[i]->sum);
				/*if (isnan(f) || !_finite(f))*/
				if (!PHR_ISFINITE((double) f))
				{
					f = 0;
				}
				if (f == 0 && x[i]->moles == 0)
				{
					x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
					continue;
				}
				else if (f == 0)
				{
					repeat = TRUE;
					x[i]->master[0]->s->la += 5;
					if (x[i]->master[0]->s->la < -999.)
						x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
				}
				else if (fail == TRUE && f < 1.5 * fabs(x[i]->moles))
				{
					continue;
				}
				else if (f > 1.5 * fabs(x[i]->moles)
						 || f < 1e-5 * fabs(x[i]->moles))
				{
					weight = (f < 1e-5 * fabs(x[i]->moles)) ? 0.3 : 1.0;
					if (x[i]->moles <= 0)
					{
						x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
					}
					else
					{
						repeat = TRUE;
						d = 0;
#ifdef SKIP
						d = weight * log10(fabs(x[i]->moles / x[i]->sum));
						double d1 = d;
#else
						// avoid underflows and overflows
						if (x[i]->moles > 1e101 || x[i]->moles < 1e-101 ||
							x[i]->sum > 1e101 || x[i]->sum < 1e-101)
						{
							LDBLE d1 = log10(x[i]->moles);
							LDBLE d2 = log10(x[i]->sum);
							LDBLE d3 = d1 - d2;
							if (d3 > DBL_MAX_10_EXP/2)
							{
								d = pow(10.0, DBL_MAX_10_EXP/2.);
							}
							else if (d3 < DBL_MIN_10_EXP/2.)
							{
								d = pow(10.0, DBL_MIN_10_EXP/2.);
							}
						}
						else
						{
							d = fabs(x[i]->moles / x[i]->sum);
						}
						LDBLE d1;
						if (d > 0)
						{
							d1 = weight * log10(d);
							if (PHR_ISFINITE((double) d1) /*&& d1 < 5.0*/)
							{
								x[i]->master[0]->s->la += d1;
							}
							else
							{
								warning_msg("Adjustment to la in revise_guesses was NaN\n");
								x[i]->master[0]->s->la += 5.0;
							}
						}
						else
						{
							x[i]->master[0]->s->la += 5.0;
						}
#endif
					}
					if (debug_set == TRUE)
					{
						output_msg(sformatf(
								   "\t%5s not converged in set %d: %e\t%e\t%e\n",
								   x[i]->description, l_iter,
								   (double) x[i]->sum, (double) x[i]->moles,
								   (double) x[i]->master[0]->s->la));
					}
				}
			}
			else if (x[i]->type == ALK)
			{
				f = total_co2;
				if (fail == TRUE && f < 1.5 * fabs(x[i]->moles))
				{
					continue;
				}
				if (f > 1.5 * fabs(x[i]->moles)
					|| f < 1e-5 * fabs(x[i]->moles))
				{
					repeat = TRUE;
					weight = (f < 1e-5 * fabs(x[i]->moles)) ? 0.3 : 1.0;
					x[i]->master[0]->s->la += weight *
						log10(fabs(x[i]->moles / x[i]->sum));
					if (debug_set == TRUE)
					{
						output_msg(sformatf(
								   "%s not converged in set. %e\t%e\t%e\n",
								   x[i]->description, (double) x[i]->sum,
								   (double) x[i]->moles,
								   (double) x[i]->master[0]->s->la));
					}
				}
			}
		}
	}
	log_msg(sformatf( "Iterations in revise_guesses: %d\n", l_iter));
	mu_x = mu_unknown->f * 0.5 / mass_water_aq_x;
	if (mu_x <= 1e-8)
	{
		mu_x = 1e-8;
	}
	gammas(mu_x);
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
sum_species(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Calculates total alk, total carbon, total co2, electrical balance,
 *   total hydrogen, and total oxygen.
 *
 *   Sorts species for summing and printing based on valence state and
 *   concentrations.
 *
 *   Sums total valence states and stores in master[i]->total.
 */
	int i, j;
	struct master *master_ptr;
/*
 *   Set global variables
 */
	ph_x = -s_hplus->la;
	solution_pe_x = -s_eminus->la;
	ah2o_x = exp(s_h2o->la * LOG_10);

	density_x = 1.0;
	if (s_o2 != NULL)
		s_o2->moles = under(s_o2->lm) * mass_water_aq_x;
	if (s_h2 != NULL)
		s_h2->moles = under(s_h2->lm) * mass_water_aq_x;

/*
 *   Calculate sums
 */
	total_alkalinity = 0.0;
	total_carbon = 0.0;
	total_co2 = 0.0;
	cb_x = 0.0;
	total_ions_x = 0.0;
	total_o_x = 0.0;
	total_h_x = 0.0;
	for (i = 0; i < count_s_x; i++)
	{
		if (s_x[i]->type == EX)
			continue;
		if (s_x[i]->type == SURF)
			continue;
		cb_x += s_x[i]->z * s_x[i]->moles;
		total_ions_x += fabs(s_x[i]->z * s_x[i]->moles);
		total_alkalinity += s_x[i]->alk * s_x[i]->moles;
		total_carbon += s_x[i]->carbon * s_x[i]->moles;
		total_co2 += s_x[i]->co2 * s_x[i]->moles;

		total_h_x += s_x[i]->h * s_x[i]->moles;
		total_o_x += s_x[i]->o * s_x[i]->moles;

		if (use.Get_surface_ptr() != NULL)
		{
			if (use.Get_surface_ptr()->Get_debye_lengths() > 0 && state >= REACTION
				&& s_x[i]->type == H2O)
			{
				total_h_x -= 2 * mass_water_surfaces_x / gfw_water;
				total_o_x -= mass_water_surfaces_x / gfw_water;
			}
		}
	}
/*
 *   Sum valence states, put in master->total
 */
	for (i = 0; i < count_master; i++)
	{
		master[i]->total = 0.0;
		master[i]->total_primary = 0.0;
	}
	for (i = 0; i < count_species_list; i++)
	{
		if (species_list[i].master_s->secondary != NULL)
		{
			master_ptr = species_list[i].master_s->secondary;
		}
		else
		{
			master_ptr = species_list[i].master_s->primary;
		}
		master_ptr->total += species_list[i].s->moles * species_list[i].coef;
	}
/*
 *   Calculate mass-balance sums
 */
	for (i = 0; i < count_unknowns; i++)
	{
		if (x[i]->type == MB ||
			x[i]->type == SOLUTION_PHASE_BOUNDARY ||
			x[i]->type == EXCH ||
			x[i]->type == SURFACE ||
			(x[i]->type == CB && x[i] != ph_unknown && x[i] != pe_unknown))
		{
			x[i]->sum = 0.0;
			for (j = 0; x[i]->master[j] != NULL; j++)
			{
				x[i]->sum += x[i]->master[j]->total;
			}
		}
		else if (x[i]->type == ALK)
		{
			x[i]->sum = total_co2;
		}
	}
/*
 *   Calculate total element concentrations
 */
	for (i = 0; i < count_master; i++)
	{
		master[i]->elt->primary->total_primary += master[i]->total;
	}
	/*
	 *  Calculate isotope ratios
	 */
	calculate_values();
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
surface_model(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Use extra iterative loop to converge on g_factors
 */
	int i, debug_diffuse_layer_save, debug_model_save;
	cxxSolution *solution_ptr;
	LDBLE prev_aq_x;
/*
 *   Allocate space for g factors for diffuse layer in surface complexation
 */
	debug_diffuse_layer_save = debug_diffuse_layer;
	debug_model_save = debug_model;
	if (last_model.force_prep == TRUE)
	{
		same_model = FALSE;
	}
	else
	{
		same_model = check_same_model();
	}
	if (dl_type_x != cxxSurface::NO_DL && same_model == FALSE)
	{
		s_diff_layer.clear();
		for (i = 0; i < count_s; i++)
		{
			std::map < std::string, cxxSpeciesDL > dl;
			s_diff_layer.push_back(dl);
			for (size_t j = 0; j < use.Get_surface_ptr()->Get_surface_charges().size(); j++)
			{
				cxxSpeciesDL species_dl;
				std::string name = use.Get_surface_ptr()->Get_surface_charges()[j].Get_name();
				s_diff_layer.back()[name] = species_dl;
			}
		}
	}
	if (state >= REACTION && dl_type_x != cxxSurface::NO_DL)
	{
		if (use.Get_mix_ptr() != NULL)
		{
			mass_water_bulk_x = 0.0;
			std::map<int, LDBLE>::const_iterator cit;
			for (cit = use.Get_mix_ptr()->Get_mixComps().begin(); cit != use.Get_mix_ptr()->Get_mixComps().end(); cit++)
			{
				solution_ptr = Utilities::Rxn_find(Rxn_solution_map, cit->first);
				mass_water_bulk_x += solution_ptr->Get_mass_water() * cit->second;
			}
		}
		else
		{
			mass_water_bulk_x = use.Get_solution_ptr()->Get_mass_water();
		}
		for (i = 0; i < (int) use.Get_surface_ptr()->Get_surface_charges().size(); i++)
		{
			cxxSurfaceCharge *charge_ptr = &(use.Get_surface_ptr()->Get_surface_charges()[i]);
			mass_water_bulk_x += charge_ptr->Get_mass_water();
			if (use.Get_surface_ptr()->Get_debye_lengths() > 0)
				charge_ptr->Set_mass_water(0.0);
		}
	}
	prep();
	k_temp(tc_x, patm_x);
	if (use.Get_surface_ptr()->Get_dl_type() == cxxSurface::DONNAN_DL)
	{
		initial_surface_water();
		calc_init_donnan();
	}
	else
	{
		calc_init_g();
	}
	if (state >= REACTION && !use.Get_surface_ptr()->Get_new_def())
	{
		set(FALSE);
	}
	else
	{
		set(TRUE);
	}
	if (model() == ERROR)
		return (ERROR);
	g_iterations = 0;

	if (use.Get_surface_ptr()->Get_dl_type() == cxxSurface::DONNAN_DL)
	{
		do
		{
			g_iterations++;
			prev_aq_x = mass_water_aq_x;
			k_temp(tc_x, patm_x);
			gammas(mu_x);
			molalities(TRUE);
			mb_sums();
			if (model() == ERROR)
				return (ERROR);
			if (!use.Get_surface_ptr()->Get_related_phases()
				&& !use.Get_surface_ptr()->Get_related_rate())
				initial_surface_water();
			if (debug_model == TRUE)
			{
				output_msg(sformatf(
						   "Surface_model (Donnan approximation): %d g_iterations, %d model iterations\n",
						   g_iterations, iterations));
			}
		}
		while ((calc_all_donnan() == FALSE
				|| fabs(1 - prev_aq_x / mass_water_aq_x) > 1e-6)
			   && g_iterations < itmax);
	}
	else
	{
		do
		{
			g_iterations++;
			if (g_iterations > itmax - 10)
			{
				debug_model = TRUE;
				debug_diffuse_layer = TRUE;
			}

			gammas(mu_x);
			molalities(TRUE);
			mb_sums();
			if (model() == ERROR)
				return (ERROR);
			if (!use.Get_surface_ptr()->Get_related_phases()
				&& !use.Get_surface_ptr()->Get_related_rate())
				initial_surface_water();
			if (debug_model == TRUE)
			{
				output_msg(sformatf(
						   "Surface_model (full integration): %d g_iterations, %d iterations\n",
						   g_iterations, iterations));
			}
		}
		while (calc_all_g() == FALSE && g_iterations < itmax);
	}
	if (g_iterations >= itmax)
	{
		pr.all = TRUE;
		pr.surface = TRUE;
		pr.species = TRUE;
		pr.use = TRUE;
		print_all();
		error_msg("Did not converge on g (diffuse layer excess)", STOP);
	}
	debug_diffuse_layer = debug_diffuse_layer_save;
	debug_model = debug_model_save;
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
free_model_allocs(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   free space allocated in model
 */
	int i;
	if (x != NULL)
	{
		for (i = 0; i < max_unknowns; i++)
		{
			unknown_free(x[i]);
		}
	}
	x = (struct unknown **) free_check_null(x);
	max_unknowns = 0;
	array = (LDBLE *) free_check_null(array);
	delta = (LDBLE *) free_check_null(delta);
	residual = (LDBLE *) free_check_null(residual);
	s_x = (struct species **) free_check_null(s_x);
	count_s_x = 0;
	sum_mb1 = (struct list1 *) free_check_null(sum_mb1);
	count_sum_mb1 = 0;
	sum_mb2 = (struct list2 *) free_check_null(sum_mb2);
	count_sum_mb2 = 0;
	sum_jacob0 = (struct list0 *) free_check_null(sum_jacob0);
	count_sum_jacob0 = 0;
	sum_jacob1 = (struct list1 *) free_check_null(sum_jacob1);
	count_sum_jacob1 = 0;
	sum_jacob2 = (struct list2 *) free_check_null(sum_jacob2);
	count_sum_jacob2 = 0;
	sum_delta = (struct list2 *) free_check_null(sum_delta);
	count_sum_delta = 0;
	return (OK);
}

#ifdef SLNQ
/* ---------------------------------------------------------------------- */
int Phreeqc::
add_trivial_eqns(int rows, int cols, LDBLE * matrix)
/* ---------------------------------------------------------------------- */
{
	int r, i, j;
	r = rows;
	if (rows == cols)
		return (OK);
	if (rows > cols)
		return (ERROR);
	for (i = 0; i < cols; i++)
	{
		for (j = 0; j < rows; j++)
		{
			if (matrix[j * (cols + 1) + i] != 0.0)
				break;
		}
		if (j < rows)
			continue;
		for (j = 0; j < cols + 1; j++)
			matrix[r * (cols + 1) + j] = 0.0;
		matrix[r * (cols + 1) + i] = 1.0;
		r++;
	}
	if (r == cols)
		return (OK);
	return (ERROR);
}
#endif
#define ZERO_TOL 1.0e-30
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
ss_root(LDBLE l_a0, LDBLE l_a1, LDBLE l_kc, LDBLE l_kb, LDBLE xcaq, LDBLE xbaq)
/* ---------------------------------------------------------------------- */
{
	int i;
	LDBLE x0, y0, x1, y1, xb, miny;

/*
 *  Bracket answer
 */
	x0 = 0.0;
	x1 = 0.0;
	y0 = ss_f(x0, l_a0, l_a1, l_kc, l_kb, xcaq, xbaq);
	miny = fabs(y0);
	for (i = 1; i <= 10; i++)
	{
		x1 = (LDBLE) i / 10;
		y1 = ss_f(x1, l_a0, l_a1, l_kc, l_kb, xcaq, xbaq);
		if (fabs(y1) < miny)
		{
			miny = fabs(y1);
		}
		if (y0 * y1 < 0)
		{
			break;
		}
		else
		{
			x0 = x1;
			y0 = y1;
		}
	}
/*
 *  Interval halve
 */
	if (i > 10)
	{
		xb = 0.0;
	}
	else
	{
		xb = ss_halve(l_a0, l_a1, x0, x1, l_kc, l_kb, xcaq, xbaq);
	}
	return (xb);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
ss_halve(LDBLE l_a0, LDBLE l_a1, LDBLE x0, LDBLE x1, LDBLE l_kc, LDBLE l_kb,
		  LDBLE xcaq, LDBLE xbaq)
/* ---------------------------------------------------------------------- */
{
	int i;
	LDBLE l_x, y0, dx, y;

	y0 = ss_f(x0, l_a0, l_a1, l_kc, l_kb, xcaq, xbaq);
	dx = (x1 - x0);
/*
 *  Loop for interval halving
 */
	for (i = 0; i < 100; i++)
	{
		dx *= 0.5;
		l_x = x0 + dx;
		y = ss_f(l_x, l_a0, l_a1, l_kc, l_kb, xcaq, xbaq);
		if (dx < 1e-8 || y == 0)
		{
			break;
		}
		if (y0 * y >= 0)
		{
			x0 = l_x;
			y0 = y;
		}
	}
	return (x0 + dx);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
ss_f(LDBLE xb, LDBLE l_a0, LDBLE l_a1, LDBLE l_kc, LDBLE l_kb, LDBLE xcaq,
	  LDBLE xbaq)
/* ---------------------------------------------------------------------- */
{
/*
 *  Need root of this function to determine xb
 */
	LDBLE lb, lc, f, xc, r;
	xc = 1 - xb;
	if (xb == 0)
		xb = 1e-20;
	if (xc == 0)
		xc = 1e-20;
	lc = exp((l_a0 - l_a1 * (-4 * xb + 3)) * xb * xb);
	lb = exp((l_a0 + l_a1 * (4 * xb - 1)) * xc * xc);
	r = lc * l_kc / (lb * l_kb);
	f = xcaq * (xb / r + xc) + xbaq * (xb + r * xc) - 1;
	return (f);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
numerical_jacobian(void)
/* ---------------------------------------------------------------------- */
{
	LDBLE *base;
	LDBLE d, d1, d2;
	int i, j;
	cxxGasPhase *gas_phase_ptr = use.Get_gas_phase_ptr();
	if (!
		(numerical_deriv || 
		(use.Get_surface_ptr() != NULL && use.Get_surface_ptr()->Get_type() == cxxSurface::CD_MUSIC) ||
		(gas_phase_ptr != NULL && gas_phase_ptr->Get_type() == cxxGasPhase::GP_VOLUME && 
		(gas_phase_ptr->Get_pr_in() || force_numerical_fixed_volume) && numerical_fixed_volume) 
		))
		return(OK);

	calculating_deriv = TRUE;
	gammas(mu_x);
	molalities(TRUE);
	mb_sums();
	residuals();
/*
 *   Clear array, note residuals are in array[i, count_unknowns+1]
 */

	for (i = 0; i < count_unknowns; i++)
	{
		array[i] = 0.0;
	}
	for (i = 1; i < count_unknowns; i++)
	{
		memcpy((void *) &(array[i * (count_unknowns + 1)]),
			   (void *) &(array[0]), (size_t) count_unknowns * sizeof(LDBLE));
	}

	base = (LDBLE *) PHRQ_malloc((size_t) count_unknowns * sizeof(LDBLE));
	if (base == NULL)
		malloc_error();
	for (i = 0; i < count_unknowns; i++)
	{
		base[i] = residual[i];
	}
	d = 0.0001;
	d1 = d * log(10.0);
	d2 = 0;
	for (i = 0; i < count_unknowns; i++)
	{
		switch (x[i]->type)
		{
		case MB:
		case ALK:
		case CB:
		case SOLUTION_PHASE_BOUNDARY:
		case EXCH:
		case SURFACE:
		case SURFACE_CB:
		case SURFACE_CB1:
		case SURFACE_CB2:
			x[i]->master[0]->s->la += d;
			d2 = d1;
			break;
		case MH:
			s_eminus->la += d;
			d2 = d1;
			break;
		case AH2O:
			x[i]->master[0]->s->la += d;
			d2 = d1;
			break;
		case PITZER_GAMMA:
			x[i]->s->lg += d;
			d2 = d;
			break;
		case MH2O:
			mass_water_aq_x *= (1.0 + d);
			x[i]->master[0]->s->moles = mass_water_aq_x / gfw_water;
			d2 = log(1.0 + d);
			break;
		case MU:
			d2 = d * mu_x;
			mu_x += d2;
			gammas(mu_x);
			break;
		case PP:
			for (j = 0; j < count_unknowns; j++)
			{
				delta[j] = 0.0;
			}
			d2 = -1e-8;
			delta[i] = d2;
			reset();
			d2 = delta[i];
			break;
		case SS_MOLES:
			if (x[i]->ss_in == FALSE)
				continue;
			for (j = 0; j < count_unknowns; j++)
			{
				delta[j] = 0.0;
			}
			/*d2 = -1e-8; */
			d2 = -d * x[i]->moles;
			d2 = -.1 * x[i]->moles;
			/*
			   if (d2 > -1e-10) d2 = -1e-10;
			   calculating_deriv = FALSE;
			 */
			delta[i] = d2;
			/*fprintf (stderr, "delta before reset %e\n", delta[i]); */
			reset();
			d2 = delta[i];
			/*fprintf (stderr, "delta after reset %e\n", delta[i]); */
			break;
		case GAS_MOLES:
			if (gas_in == FALSE)
				continue;

			d2 = d * x[i]->moles;
			if (d2 < 1e-14)
				d2 = 1e-14;
			x[i]->moles += d2;
			break;
		}
		molalities(TRUE);
		mb_sums();
		/*
		   mb_ss();
		   mb_gases();
		 */
		residuals();
		//output_msg(sformatf( "%d\n", i));
		for (j = 0; j < count_unknowns; j++)
		{
			// avoid overflow
			if (residual[j] > 1.0e101)
			{
			  LDBLE t = (LDBLE) pow((LDBLE) 10.0, (LDBLE) (DBL_MAX_10_EXP - 50.0));
				if (residual[j]  > t)
				{
					array[j * (count_unknowns + 1) + i] = -pow(10.0, DBL_MAX_10_EXP - 50.0);
				}
				else
				{
					array[j * (count_unknowns + 1) + i] = -(residual[j] - base[j]) / d2;
				}
			}
			else if (residual[j] < -1.0e101)
			{
				LDBLE t = pow((LDBLE) 10.0, (LDBLE) (DBL_MIN_10_EXP + 50.0));
				if (residual[j]  < -t)
				{
					array[j * (count_unknowns + 1) + i] = pow(10.0, DBL_MIN_10_EXP + 50.0);
				}
				else
				{
					array[j * (count_unknowns + 1) + i] = -(residual[j] - base[j]) / d2;
				}
			}
			else
			{
				array[j * (count_unknowns + 1) + i] = -(residual[j] - base[j]) / d2;
				if (!PHR_ISFINITE(array[j * (count_unknowns + 1) + i]))
				{
					//fprintf(stderr, "oops, got NaN: %e, %e, %e, %e\n", residual[j], base[j], d2, array[j * (count_unknowns + 1) + i]);
				}
			}

			//output_msg(sformatf( "\t%d %e %e %e %e\n", j, array[j*(count_unknowns + 1) + i] , residual[j], base[j], d2));
		}
		switch (x[i]->type)
		{
		case MB:
		case ALK:
		case CB:
		case SOLUTION_PHASE_BOUNDARY:
		case EXCH:
		case SURFACE:
		case SURFACE_CB:
		case SURFACE_CB1:
		case SURFACE_CB2:
		case AH2O:
			x[i]->master[0]->s->la -= d;
			break;
		case MH:
			s_eminus->la -= d;
			if (array[i * (count_unknowns + 1) + i] == 0)
			{
				/*output_msg(sformatf( "Zero diagonal for MH\n")); */
				array[i * (count_unknowns + 1) + i] =
				  under(s_h2->lm) * 2;
			}
			break;
		case PITZER_GAMMA:
			x[i]->s->lg -= d;
			break;
		case MH2O:
			mass_water_aq_x /= (1 + d);
			x[i]->master[0]->s->moles = mass_water_aq_x / gfw_water;
			break;
		case MU:
			mu_x -= d2;
			gammas(mu_x);
			break;
		case PP:
			delta[i] = -d2;
			reset();
			break;
		case SS_MOLES:
			delta[i] = -d2;
			reset();
			break;
		case GAS_MOLES:
			x[i]->moles -= d2;
			break;
		}
	}
	molalities(TRUE);
	mb_sums();
	mb_gases();
	mb_ss();
	residuals();
	free_check_null(base);
	calculating_deriv = FALSE;
	return OK;
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
ineq_init(int l_max_row_count, int l_max_column_count)
/* ---------------------------------------------------------------------- */
{
	if (normal == NULL)
	{
		normal =
			(LDBLE *) PHRQ_malloc((size_t) count_unknowns * sizeof(LDBLE));
		normal_max = count_unknowns;
		if (normal == NULL)
			malloc_error();
	}
	if (ineq_array == NULL)
	{
		ineq_array =
			(LDBLE *) PHRQ_malloc((size_t) l_max_row_count * l_max_column_count *
								  sizeof(LDBLE));
		if (ineq_array == NULL)
			malloc_error();
		ineq_array_max = l_max_row_count * l_max_column_count;
	}
	if (back_eq == NULL)
	{
		back_eq = (int *) PHRQ_malloc((size_t) l_max_row_count * sizeof(int));
		if (back_eq == NULL)
			malloc_error();
		back_eq_max = l_max_row_count;
	}
	if (zero == NULL)
	{
		zero = (LDBLE *) PHRQ_malloc((size_t) l_max_row_count * sizeof(LDBLE));
		if (zero == NULL)
			malloc_error();
		zero_max = l_max_row_count;
	}
	if (res == NULL)
	{
		res = (LDBLE *) PHRQ_malloc((size_t) l_max_row_count * sizeof(LDBLE));
		if (res == NULL)
			malloc_error();
		res_max = l_max_row_count;
	}
	if (delta1 == NULL)
	{
		delta1 =
			(LDBLE *) PHRQ_malloc((size_t) l_max_column_count * sizeof(LDBLE));
		if (delta1 == NULL)
			malloc_error();
		delta1_max = l_max_column_count;
	}
	if (cu == NULL)
	{
		cu = (LDBLE *) PHRQ_malloc((size_t) 3 * l_max_row_count *
								   sizeof(LDBLE));
		if (cu == NULL)
			malloc_error();
		cu_max = 3 * l_max_row_count;
	}
	if (iu == NULL)
	{
		iu = (int *) PHRQ_malloc((size_t) 3 * l_max_row_count * sizeof(int));
		if (iu == NULL)
			malloc_error();
		iu_max = 3 * l_max_row_count;
	}
	if (is == NULL)
	{
		is = (int *) PHRQ_malloc((size_t) 3 * l_max_row_count * sizeof(int));
		if (is == NULL)
			malloc_error();
		is_max = 3 * l_max_row_count;
	}
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
set_inert_moles(void)
/* ---------------------------------------------------------------------- */
{
	int j;
	if (use.Get_pp_assemblage_ptr() == NULL) return;
	
	for (j = 0; j < count_unknowns; j++)
	{
		if (x[j]->type != PP) continue;
		//cxxPPassemblageComp * comp_ptr = pp_assemblage_ptr->Find(x[j]->pp_assemblage_comp_name); 
		cxxPPassemblageComp * comp_ptr = (cxxPPassemblageComp *) x[j]->pp_assemblage_comp_ptr;
		if (comp_ptr->Get_precipitate_only())
		{
			x[j]->inert_moles = x[j]->moles;
			x[j]->moles = 0;
		}
	}
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
unset_inert_moles()
/* ---------------------------------------------------------------------- */
{
	int j;
	if (use.Get_pp_assemblage_ptr() == NULL) return;
	for (j = 0; j < count_unknowns; j++)
	{
		if (x[j]->type != PP) continue;
		//cxxPPassemblageComp * comp_ptr = pp_assemblage_ptr->Find(x[j]->pp_assemblage_comp_name); 		
		cxxPPassemblageComp * comp_ptr = (cxxPPassemblageComp *) x[j]->pp_assemblage_comp_ptr;
		if (comp_ptr->Get_precipitate_only())
		{
			x[j]->moles += x[j]->inert_moles;
			x[j]->inert_moles = 0;
		}
	}
}
