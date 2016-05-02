#include "Phreeqc.h"
#include "GasPhase.h"
/* ---------------------------------------------------------------------- */
int Phreeqc::
setup_fixed_volume_gas(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Fill in data for gas phase unknown (sum of partial pressures)
 *   in unknown structure
 */
	if (use.Get_gas_phase_ptr() == NULL)
		return (OK);
	cxxGasPhase *gas_phase_ptr = use.Get_gas_phase_ptr();
/*
 *   One for each gas component
 */
	gas_unknowns.clear();
	gas_unknown = NULL;
	gas_phase_ptr->Set_total_moles(0);
	for (size_t i = 0; i < gas_phase_ptr->Get_gas_comps().size(); i++)
	{
		const cxxGasComp *comp_ptr = &(gas_phase_ptr->Get_gas_comps()[i]);
		int j;
		struct phase *phase_ptr = phase_bsearch(comp_ptr->Get_phase_name().c_str(), &j, FALSE);
		x[count_unknowns]->type = GAS_MOLES;
		x[count_unknowns]->description = phase_ptr->name;
		x[count_unknowns]->phase = phase_ptr;
		x[count_unknowns]->moles = comp_ptr->Get_moles();
		if (x[count_unknowns]->moles <= 0)
		{
			x[count_unknowns]->moles = MIN_TOTAL;
		}
		x[count_unknowns]->ln_moles = log(x[count_unknowns]->moles);
		gas_unknowns.push_back(x[count_unknowns]);
		gas_phase_ptr->Set_total_moles(gas_phase_ptr->Get_total_moles() + x[count_unknowns]->moles);
		x[count_unknowns]->phase->moles_x = x[count_unknowns]->moles;
		count_unknowns++;
	}
	if (gas_unknowns.size() > 0)
	{
		gas_unknown = gas_unknowns[0];
	}

	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
build_fixed_volume_gas(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Put coefficients into lists to sum iaps to test for equilibrium
 *   Put coefficients into lists to build jacobian for 
 *      sum of partial pressures equation and
 *      mass balance equations for elements contained in gases
 */
	int row, col;
	struct master *master_ptr;
	struct rxn_token *rxn_ptr;
	struct unknown *unknown_ptr;
	LDBLE coef, coef_elt;

	if (gas_unknown == NULL)
		return (OK);
	cxxGasPhase *gas_phase_ptr = use.Get_gas_phase_ptr();
	for (size_t i = 0; i < gas_phase_ptr->Get_gas_comps().size(); i++)
	{
		const cxxGasComp *comp_ptr = &(gas_phase_ptr->Get_gas_comps()[i]);
		int j;
		struct phase *phase_ptr = phase_bsearch(comp_ptr->Get_phase_name().c_str(), &j, FALSE);
/*
 *   Determine elements in gas component
 */
		count_elts = 0;
		paren_count = 0;
		if (phase_ptr->rxn_x == NULL)
			continue;
		add_elt_list(phase_ptr->next_elt, 1.0);
#define COMBINE
#ifdef COMBINE
		change_hydrogen_in_elt_list(0);
#endif
/*
 *   Build mass balance sums for each element in gas
 */
		if (debug_prep == TRUE)
		{
			output_msg(sformatf( "\n\tMass balance summations %s.\n\n",
					   phase_ptr->name));
		}

		/* All elements in gas */
		for (j = 0; j < count_elts; j++)
		{
			unknown_ptr = NULL;
			if (strcmp(elt_list[j].elt->name, "H") == 0)
			{
				unknown_ptr = mass_hydrogen_unknown;
			}
			else if (strcmp(elt_list[j].elt->name, "O") == 0)
			{
				unknown_ptr = mass_oxygen_unknown;
			}
			else
			{
				if (elt_list[j].elt->primary->in == TRUE)
				{
					unknown_ptr = elt_list[j].elt->primary->unknown;
				}
				else if (elt_list[j].elt->primary->s->secondary != NULL)
				{
					unknown_ptr =
						elt_list[j].elt->primary->s->secondary->unknown;
				}
			}
			if (unknown_ptr != NULL)
			{
				coef = elt_list[j].coef;
				store_mb(&(gas_unknowns[i]->moles), &(unknown_ptr->f), coef);
				if (debug_prep == TRUE)
				{
					output_msg(sformatf( "\t\t%-24s%10.3f\n",
							   unknown_ptr->description, (double) coef));
				}
			}
		}
		if (gas_phase_ptr->Get_type() == cxxGasPhase::GP_PRESSURE)
		{
			/* Total pressure of gases */
			store_mb(&(phase_ptr->p_soln_x), &(gas_unknown->f),
					 1.0);
		}
/*
 *   Build jacobian sums for mass balance equations
 */
		if (debug_prep == TRUE)
		{
			output_msg(sformatf( "\n\tJacobian summations %s.\n\n",
					   phase_ptr->name));
		}
		for (j = 0; j < count_elts; j++)
		{
			unknown_ptr = NULL;
			if (strcmp(elt_list[j].elt->name, "H") == 0)
			{
				unknown_ptr = mass_hydrogen_unknown;
			}
			else if (strcmp(elt_list[j].elt->name, "O") == 0)
			{
				unknown_ptr = mass_oxygen_unknown;
			}
			else
			{
				if (elt_list[j].elt->primary->in == TRUE)
				{
					unknown_ptr = elt_list[j].elt->primary->unknown;
				}
				else if (elt_list[j].elt->primary->s->secondary != NULL)
				{
					unknown_ptr =
						elt_list[j].elt->primary->s->secondary->unknown;
				}
			}
			if (unknown_ptr == NULL)
			{
				continue;
			}
			if (debug_prep == TRUE)
			{
				output_msg(sformatf( "\n\t%s.\n",
						   unknown_ptr->description));
			}
			row = unknown_ptr->number * (count_unknowns + 1);
			coef_elt = elt_list[j].coef;
			for (rxn_ptr = phase_ptr->rxn_x->token + 1;
				 rxn_ptr->s != NULL; rxn_ptr++)
			{

				if (rxn_ptr->s->secondary != NULL
					&& rxn_ptr->s->secondary->in == TRUE)
				{
					master_ptr = rxn_ptr->s->secondary;
				}
				else if (rxn_ptr->s->primary != NULL && rxn_ptr->s->primary->in == TRUE)
				{
					master_ptr = rxn_ptr->s->primary;
				}
				else
				{
					master_ptr = master_bsearch_primary(rxn_ptr->s->name);
					master_ptr->s->la = -999.0;
				}
				if (debug_prep == TRUE)
				{
					output_msg(sformatf( "\t\t%s\n",
							   master_ptr->s->name));
				}
				if (master_ptr->unknown == NULL)
				{
					continue;
				}
				if (master_ptr->in == FALSE)
				{
					error_string = sformatf(
							"Element, %s, in phase, %s, is not in model.",
							master_ptr->elt->name, phase_ptr->name);
					error_msg(error_string, CONTINUE);
					input_error++;
				}
				col = master_ptr->unknown->number;
				coef = coef_elt * rxn_ptr->coef;
				store_jacob(&(gas_unknowns[i]->moles),
							&(array[row + col]), coef);
				if (debug_prep == TRUE)
				{
					output_msg(sformatf( "\t\t%-24s%10.3f\t%d\t%d\n",
							   master_ptr->s->name, (double) coef,
							   row / (count_unknowns + 1), col));
				}
			}
			if (gas_phase_ptr->Get_type() == cxxGasPhase::GP_PRESSURE)
			{
				/* derivative wrt total moles of gas */
				store_jacob(&(phase_ptr->fraction_x),
							&(array[row + gas_unknown->number]), coef_elt);
				if (debug_prep == TRUE)
				{
					output_msg(sformatf( "\t\t%-24s%10.3f\t%d\t%d\n",
							   "gas moles", (double) elt_list[j].coef,
							   row / (count_unknowns + 1),
							   gas_unknown->number));
				}
			}
		}
/*
 *   Build jacobian sums for sum of partial pressures equation
 */
		if (gas_phase_ptr->Get_type() != cxxGasPhase::GP_PRESSURE)
			continue;
		if (debug_prep == TRUE)
		{
			output_msg(sformatf( "\n\tPartial pressure eqn %s.\n\n",
					   phase_ptr->name));
		}
		unknown_ptr = gas_unknown;
		row = unknown_ptr->number * (count_unknowns + 1);
		for (rxn_ptr = phase_ptr->rxn_x->token + 1; rxn_ptr->s != NULL; rxn_ptr++)
		{
			if (rxn_ptr->s != s_eminus && rxn_ptr->s->in == FALSE)
			{
				error_string = sformatf(
					"Element in species, %s, in phase, %s, is not in model.",
					rxn_ptr->s->name, phase_ptr->name);
				warning_msg(error_string);
			}
			else
			{
				if (rxn_ptr->s->secondary != NULL
					&& rxn_ptr->s->secondary->in == TRUE)
				{
					master_ptr = rxn_ptr->s->secondary;
				}
				else if (rxn_ptr->s->primary != NULL && rxn_ptr->s->primary->in == TRUE)
				{
					master_ptr = rxn_ptr->s->primary;
				}
				else
				{
					master_ptr = master_bsearch_primary(rxn_ptr->s->name);
					if (master_ptr && master_ptr->s)
					{
						master_ptr->s->la = -999.0;
					}
				}

				if (master_ptr == NULL)
				{
					error_string = sformatf(
						"Master species for %s, in phase, %s, is not in model.",
						rxn_ptr->s->name, phase_ptr->name);
					error_msg(error_string, CONTINUE);
					input_error++;
				}
				else
				{
					if (debug_prep == TRUE)
					{
						output_msg(sformatf( "\t\t%s\n", master_ptr->s->name));
					}
					if (master_ptr->unknown == NULL)
					{
						assert(false);
						continue;
					}
					if (master_ptr->in == FALSE)
					{
						error_string = sformatf(
							"Element, %s, in phase, %s, is not in model.",
							master_ptr->elt->name, phase_ptr->name);
						warning_msg(error_string);
					}
					col = master_ptr->unknown->number;
					coef = rxn_ptr->coef;
					store_jacob(&(phase_ptr->p_soln_x), &(array[row + col]), coef);
					if (debug_prep == TRUE)
					{
						output_msg(sformatf( "\t\t%-24s%10.3f\t%d\t%d\n",
							master_ptr->s->name, (double) coef,
							row / (count_unknowns + 1), col));
					}
				}
			}
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
calc_PR(void)
/* ---------------------------------------------------------------------- */
/*  Calculate fugacity and fugacity coefficient for gas pressures if critical T and P
    are defined.
  1) Solve molar volume V_m or total pressure P from Peng-Robinson's EOS:
  P = R * T / (V_m - b) - a * aa / (V_m^2 + 2 * b * V_m - b^2)
     a = 0.457235 * (R * T_c)^2 / P_c
     b = 0.077796 * R * T_c / P_c
     aa = (1 + kk * (1 - T_r^0.5))^2
     kk = 0.37464 + 1.54226 * omega - 0.26992 * omega^2
     T_r = T / T_c
  multicomponent gas phase:
     use: b_sum = Sum(x_i * b), x_i is mole-fraction
          a_aa_sum = Sum_i( Sum_j(x_i * x_j * (a_i * aa_i * a_j * aa_j)^0.5) )
  2) Find the fugacity coefficient phi for gas i:
  log(phi_i) = B_ratio * (z - 1) - log(z - B) + A / (2.8284 * B) * (B_ratio - 2 / a_aa_sum * a_aa_sum2) *\
           log((z + 2.4142 * B) / (z - 0.4142 * B))
     B_ratio = b_i / b_sum
     A = a_aa_sum * P / R_TK^2
     B = b_sum * P / R_TK
     a_aa_sum2 = Sum_j(x_j * (a_aa_i * a_aa_j)^0.5
  3) correct the solubility of gas i with:
  pr_si_f = log10(phi_i) -  Delta_V_i * (P - 1) / (2.303 * R * TK);
*/
{
	LDBLE T_c, P_c;
	LDBLE A, B, B_r, /*b2,*/ kk, oo, a_aa, T_r;
	LDBLE m_sum, /*b_sum, a_aa_sum,*/ a_aa_sum2;
	LDBLE phi;
	LDBLE /*R_TK,*/ R = R_LITER_ATM; /* L atm / (K mol) */
	LDBLE r3[4], r3_12, rp, rp3, rq, rz, ri, ri1, one_3 = 0.33333333333333333;
	LDBLE disct, vinit, v1, ddp, dp_dv, dp_dv2;
	int it;
	struct phase *phase_ptr;
	LDBLE V_m = 0, P = 0;

	LDBLE TK = tk_x;
	bool halved;
	R_TK = R * TK;
	m_sum = b_sum = a_aa_sum = 0.0;
	size_t i;

	if (gas_unknowns.size() == 0)
	{
		error_msg("No gas unknowns.", STOP);
	}
	cxxGasPhase * gas_phase_ptr = use.Get_gas_phase_ptr();
	
	for (i = 0; i < gas_unknowns.size(); i++)
	{
		m_sum += gas_unknowns[i]->moles;
		phase_ptr = gas_unknowns[i]->phase;
		if (phase_ptr->t_c == 0.0 || phase_ptr->p_c == 0.0)
			error_msg("Cannot calculate a mixture of ideal and Peng_Robinson gases,\n       please define Tc and Pc for the active gases in PHASES.", STOP);
			//continue;
		if (!phase_ptr->pr_a)
		{
			T_c = phase_ptr->t_c;
			P_c = phase_ptr->p_c;
			phase_ptr->pr_a = 0.457235 * R * R * T_c * T_c / P_c;
			phase_ptr->pr_b = 0.077796 * R * T_c / P_c;
			T_r = TK / T_c;
			oo = phase_ptr->omega;
			kk = 0.37464 + oo * (1.54226 - 0.26992 * oo);
			phase_ptr->pr_alpha = pow(1 + kk * (1 - sqrt(T_r)), 2);
			phase_ptr->pr_tk = TK;
			phase_ptr->pr_in = true;
		}
		if (phase_ptr->pr_tk != TK)
		{
			T_r = TK / phase_ptr->t_c;
			oo = phase_ptr->omega;
			kk = 0.37464 + oo * (1.54226 - 0.26992 * oo);
			phase_ptr->pr_alpha = pow(1 + kk * (1 - sqrt(T_r)), 2);
			phase_ptr->pr_tk = TK;
			phase_ptr->pr_in = true;
		}
	}
	if (m_sum == 0)
			return (OK);
	gas_phase_ptr->Set_v_m(gas_phase_ptr->Get_volume() / m_sum);
	for (i = 0; i < gas_unknowns.size(); i++)
	{
		phase_ptr = gas_unknowns[i]->phase;
		phase_ptr->fraction_x = gas_unknowns[i]->moles / m_sum;							// phase_ptr->fraction_x updated
	}

	for (i = 0; i < gas_unknowns.size(); i++)
	{
		a_aa_sum2 = 0.0;
		phase_ptr = gas_unknowns[i]->phase;
		//if (phase_ptr->t_c == 0.0 || phase_ptr->p_c == 0.0)
		//	continue;
		b_sum += phase_ptr->fraction_x * phase_ptr->pr_b;
		size_t i1;
		struct phase *phase_ptr1;
		for (i1 = 0; i1 <  gas_unknowns.size(); i1++)
		{
			phase_ptr1 = gas_unknowns[i1]->phase;
			//if (phase_ptr1->t_c == 0.0 || phase_ptr1->p_c == 0.0)
			//	continue;
			if (phase_ptr1->fraction_x == 0)
				continue;
			a_aa = sqrt(phase_ptr->pr_a * phase_ptr->pr_alpha *
				        phase_ptr1->pr_a * phase_ptr1->pr_alpha);
			if (!strcmp(phase_ptr->name, "H2O(g)"))
			{
				if (!strcmp(phase_ptr1->name, "CO2(g)"))
					a_aa *= 0.81; // Soreide and Whitson, 1992, FPE 77, 217
				else if (!strcmp(phase_ptr1->name, "H2S(g)"))
					a_aa *= 0.81;
				else if (!strcmp(phase_ptr1->name, "CH4(g)"))
					a_aa *= 0.51;
				else if (!strcmp(phase_ptr1->name, "N2(g)"))
					a_aa *= 0.51;
			}
			if (!strcmp(phase_ptr1->name, "H2O(g)"))
			{
				if (!strcmp(phase_ptr->name, "CO2(g)"))
					a_aa *= 0.81;
				else if (!strcmp(phase_ptr->name, "H2S(g)"))
					a_aa *= 0.81;
				else if (!strcmp(phase_ptr->name, "CH4(g)"))
					a_aa *= 0.51;
				else if (!strcmp(phase_ptr->name, "N2(g)"))
					a_aa *= 0.51;
			}
			a_aa_sum += phase_ptr->fraction_x * phase_ptr1->fraction_x * a_aa;
			a_aa_sum2 += phase_ptr1->fraction_x * a_aa;
		}
		phase_ptr->pr_aa_sum2 = a_aa_sum2;
	}
	b2 = b_sum * b_sum;

	if (gas_phase_ptr->Get_type() == cxxGasPhase::GP_VOLUME)
	{
		V_m = gas_phase_ptr->Get_volume() / m_sum;
		P = 0.0;
		while (P <= 0)
		{
			P = R_TK / (V_m - b_sum) - a_aa_sum / (V_m * (V_m + 2 * b_sum) - b2);
			if (P <= 0.0)
			{
				V_m *= 2.0;
			//a_aa_sum /= 2.0;
			}
		}
		if (iterations > 0 && P < 150 && V_m < 1.01)
		{
			// check for 3-roots...
			r3[1] = b_sum - R_TK / P;
			r3[2] = -3.0 * b2 + (a_aa_sum - R_TK * 2.0 * b_sum) / P;
			r3[3] = b2 * b_sum + (R_TK * b2 - b_sum * a_aa_sum) / P;
			// the discriminant of the cubic eqn...
			disct = 18. * r3[1] * r3[2] * r3[3] -
				4. * pow(r3[1], 3) * r3[3] + 
				r3[1] * r3[1] * r3[2] * r3[2] -
				4. * pow(r3[2], 3) - 
				27. * r3[3] * r3[3];
			if (disct > 0)
			{
				// 3-roots, find the largest P...
				it = 0;
				halved = false;
				ddp = 1e-9;
				v1 = vinit = 0.729;
				dp_dv = f_Vm(v1, this);
				while (fabs(dp_dv) > 1e-11 && it < 40)
				{
					it +=1;
					dp_dv2 = f_Vm(v1 - ddp, this);
					v1 -= (dp_dv * ddp / (dp_dv - dp_dv2));
					if (!halved && (v1 > vinit || v1 < 0.03))
					{
						if (vinit > 0.329)
							vinit -= 0.1;
						else
							vinit -=0.05;
						if (vinit < 0.03)
						{
							vinit = halve(f_Vm, 0.03, 1.0, 1e-3);
							if (f_Vm(vinit - 2e-3, this) < 0)
								vinit = halve(f_Vm, vinit + 2e-3, 1.0, 1e-3);
							halved = true;
						}
						v1 = vinit;
					}
					dp_dv = f_Vm(v1, this);
					if (fabs(dp_dv) < 1e-11)
					{
						if (f_Vm(v1 - 1e-4, this) < 0)
						{
							v1 = halve(f_Vm, v1 + 1e-4, 1.0, 1e-3);
							dp_dv = f_Vm(v1, this);
						}
					}
				}
				if (it == 40)
				{
// accept a (possible) whobble in the curve...
//					error_msg("No convergence when calculating P in Peng-Robinson.", STOP);
				}
				if (V_m < v1 && it < 40)
					P = R_TK / (v1 - b_sum) - a_aa_sum / (v1 * (v1 + 2 * b_sum) - b2);
			}
		}
		if (P <= 0) // iterations = -1
			P = 1.;
		gas_phase_ptr->Set_total_p(P);												// phase_ptr->total_p updated
		gas_phase_ptr->Set_v_m(V_m);
	}
	else
	{
		assert(false);
		P = gas_phase_ptr->Get_total_p();
		r3[1] = b_sum - R_TK / P;
		r3_12 = r3[1] * r3[1];
		r3[2] = -3.0 * b2 + (a_aa_sum - R_TK * 2.0 * b_sum) / P;
		r3[3] = b2 * b_sum + (R_TK * b2 - b_sum * a_aa_sum) / P;
		// solve t^3 + rp*t + rq = 0.
		// molar volume V_m = t - r3[1] / 3... 
		rp = r3[2] - r3_12 / 3;
		rp3 = rp * rp * rp;
		rq = (2.0 * r3_12 * r3[1] - 9.0 * r3[1] * r3[2]) / 27 + r3[3];
		rz = rq * rq / 4 + rp3 / 27;
		if (rz >= 0) // Cardono's method...
		{
			ri = sqrt(rz);
			if (ri + rq / 2 <= 0)
			{
				V_m = pow(ri - rq / 2, one_3) + pow(- ri - rq / 2, one_3) - r3[1] / 3;
			} else
			{
				ri = - pow(ri + rq / 2, one_3);
				V_m = ri - rp / (3.0 * ri) - r3[1] / 3;
			}
		} else // use complex plane...
		{
			ri = sqrt(- rp3 / 27); // rp < 0
			ri1 = acos(- rq / 2 / ri);
			V_m = 2.0 * pow(ri, one_3) * cos(ri1 / 3) - r3[1] / 3;
		}
		gas_phase_ptr->Set_v_m(V_m);												 // phase_ptr->fraction_x updated
	}
 // calculate the fugacity coefficients...
	for (i = 0; i < gas_unknowns.size(); i++)
	{
		phase_ptr = gas_unknowns[i]->phase;
		if (phase_ptr->fraction_x == 0.0)
		{
			phase_ptr->pr_p = 0;
			phase_ptr->pr_phi = 1;
			phase_ptr->pr_si_f = 0.0;
			//phase_ptr->logk[vm0] = 0.0; // ***
			continue;
		}	
		phase_ptr->pr_p = phase_ptr->fraction_x * P;

		//if (phase_ptr->t_c == 0.0 || phase_ptr->p_c == 0.0)
		//{
		//	phase_ptr->pr_phi = 1;
		//	continue;
		//}
		rz = P * V_m / R_TK;
		A = a_aa_sum * P / (R_TK * R_TK);
		B = b_sum * P / R_TK;
		B_r = phase_ptr->pr_b / b_sum;
		if (rz > B)
		{
			phi = B_r * (rz - 1) - log(rz - B) + A / (2.828427 * B) * (B_r - 2.0 * phase_ptr->pr_aa_sum2 / a_aa_sum) *
				  log((rz + 2.41421356 * B) / (rz - 0.41421356 * B));
			if (phi > 4.44)
				phi = 4.44;
		}
		else
			phi = -3.0; // fugacity coefficient > 0.05
		phase_ptr->pr_phi = exp(phi);
		phase_ptr->pr_si_f = phi / LOG_10;										// pr_si_f updated
		// ****
		//phase_ptr->logk[vm0] = V_m * 1e3 * phase_ptr->fraction_x;
		phase_ptr->pr_in = true;
	}
	return (V_m);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
calc_fixed_volume_gas_pressures(void)
/* ---------------------------------------------------------------------- */
{
	int n_g = 0;
	LDBLE lp;
	struct rxn_token *rxn_ptr;
	struct phase *phase_ptr;
	bool PR = false, pr_done = false;
	size_t i;
/*
 *   moles and partial pressures for gases
 */
	if (use.Get_gas_phase_ptr() == NULL)
		return (OK);
	cxxGasPhase *gas_phase_ptr = use.Get_gas_phase_ptr();
	gas_phase_ptr->Set_total_moles(0);

	for (i = 0; i < gas_unknowns.size(); i++)
	{
		phase_ptr = gas_unknowns[i]->phase;
		if (phase_ptr->in == TRUE)
		{
			if (!PR && phase_ptr->t_c > 0 && phase_ptr->p_c > 0)
				PR = true;
			n_g++;
		}
		gas_phase_ptr->Set_total_moles(gas_phase_ptr->Get_total_moles() + gas_unknowns[i]->moles);
	}
	if (PR && gas_phase_ptr->Get_total_moles() > 0)
	{
		calc_PR();
		pr_done = true;
		gas_phase_ptr->Set_total_moles(0);
	} else
	{
		gas_phase_ptr->Set_total_p(0);
		gas_phase_ptr->Set_total_moles(0);
	}
	for (i = 0; i < gas_unknowns.size(); i++)
	{
		phase_ptr = gas_unknowns[i]->phase;
		if (phase_ptr->in == TRUE)
		{
			lp = -phase_ptr->lk;
			//lp = -k_calc(phase_ptr->rxn_x->logk, tk_x, use.Get_gas_phase_ptr()->total_p * PASCAL_PER_ATM);
			for (rxn_ptr = phase_ptr->rxn_x->token + 1; rxn_ptr->s != NULL;
				 rxn_ptr++)
			{
				lp += rxn_ptr->s->la * rxn_ptr->coef;
			}
			phase_ptr->p_soln_x = exp(LOG_10 * (lp - phase_ptr->pr_si_f));

			if (pr_done)
			{
				lp = phase_ptr->p_soln_x / gas_phase_ptr->Get_total_p() *
					gas_phase_ptr->Get_volume() / gas_phase_ptr->Get_v_m();
				phase_ptr->moles_x = lp;
			}
			else
			{
				phase_ptr->moles_x = phase_ptr->p_soln_x * gas_phase_ptr->Get_volume() / (R_LITER_ATM * tk_x);
				gas_phase_ptr->Set_total_p(gas_phase_ptr->Get_total_p() + phase_ptr->p_soln_x);
			}
			gas_phase_ptr->Set_total_moles(gas_phase_ptr->Get_total_moles() + phase_ptr->moles_x);
		}
		else
		{
			phase_ptr->moles_x = 0;
			phase_ptr->fraction_x = 0;
		}
	}

	return (OK);
}
