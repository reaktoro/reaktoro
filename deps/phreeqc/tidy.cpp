#include "Utils.h"
#include "Phreeqc.h"
#include "phqalloc.h"
#include "Exchange.h"
#include "GasPhase.h"
#include "PPassemblage.h"
#include "SSassemblage.h"
#include "cxxKinetics.h"
#include "Solution.h"

#define ZERO_TOL 1.0e-30

/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_model(void)
/* ---------------------------------------------------------------------- */
{
	int n_user, last;
	int new_named_logk;
	/*
	 * Determine if any new elements, species, phases have been read
	 */
	state = INITIALIZE;
	new_model = FALSE;
	new_pp_assemblage = FALSE;
	new_surface = FALSE;
	new_exchange = FALSE;
	new_reaction = FALSE;
	new_temperature = FALSE;
	new_mix = FALSE;
	new_solution = FALSE;
	new_gas_phase = FALSE;
	new_inverse = FALSE;
	new_punch = FALSE;
	new_surface = FALSE;
	new_ss_assemblage = FALSE;
	new_kinetics = FALSE;
	new_pitzer = FALSE;
	new_named_logk = FALSE;

	if (keycount[Keywords::KEY_SOLUTION_SPECIES] > 0				||	/*"species" */
		keycount[Keywords::KEY_SOLUTION_MASTER_SPECIES] > 0			||	/*"master" */
		keycount[Keywords::KEY_PHASES] > 0							||	/*"phases" */
		keycount[Keywords::KEY_EXCHANGE_SPECIES] > 0				||	/*"exchange_species" */
		keycount[Keywords::KEY_EXCHANGE_MASTER_SPECIES] > 0			||	/*"master_exchange_species" */
		keycount[Keywords::KEY_SURFACE_SPECIES] > 0					||	/*"surface_species" */
		keycount[Keywords::KEY_SURFACE_MASTER_SPECIES] > 0			||	/*"master_surface_species" */
		keycount[Keywords::KEY_RATES] > 0							||	/*"rates" */
		keycount[Keywords::KEY_LLNL_AQUEOUS_MODEL_PARAMETERS] > 0	||	/*"llnl_aqueous_model_parameters" */
		(keycount[Keywords::KEY_DATABASE] > 0 && simulation == 0)	||	/*"database" */
		keycount[Keywords::KEY_NAMED_EXPRESSIONS] > 0				||	/*"named_analytical_expressions" */
		keycount[Keywords::KEY_ISOTOPES] > 0						||	/*"isotopes" */
		keycount[Keywords::KEY_CALCULATE_VALUES] > 0				||	/*"calculate_values" */
		keycount[Keywords::KEY_ISOTOPE_RATIOS] > 0					||	/*"isotopes_ratios", */
		keycount[Keywords::KEY_ISOTOPE_ALPHAS] > 0					||	/*"isotopes_alphas" */
		keycount[Keywords::KEY_PITZER] > 0							||	/*"pitzer" */
		keycount[Keywords::KEY_SIT] > 0								/*"sit" */
		)
	{							
		new_model = TRUE;
	}
	if (keycount[Keywords::KEY_EQUILIBRIUM_PHASES] > 0		|| 
		keycount[Keywords::KEY_EQUILIBRIUM_PHASES_RAW] > 0	||
		keycount[Keywords::KEY_EQUILIBRIUM_PHASES_MODIFY])
	{
		new_pp_assemblage = TRUE;					/*"pure_phases" */
	}
	if (keycount[Keywords::KEY_SURFACE] > 0					||
		keycount[Keywords::KEY_SURFACE_RAW] > 0				||
		keycount[Keywords::KEY_SURFACE_MODIFY])
	{
		new_surface = TRUE;							/*"surface" */
	}
	if (keycount[Keywords::KEY_EXCHANGE] > 0				||
		keycount[Keywords::KEY_EXCHANGE_RAW] > 0				||
		keycount[Keywords::KEY_EXCHANGE_MODIFY])
	{
		new_exchange = TRUE;						/*"exchange" */
	}
	if (keycount[Keywords::KEY_REACTION] > 0				/*||
		keycount[Keywords::KEY_REACTION_RAW] > 0			||
		keycount[Keywords::KEY_REACTION_MODIFY]*/)
	{
		new_reaction = TRUE;						/*"reaction" */
	}
	if (keycount[Keywords::KEY_REACTION_TEMPERATURE] > 0		/*||
		keycount[Keywords::KEY_REACTION_TEMPERATURE_RAW] > 0	||
		keycount[Keywords::KEY_REACTION_TEMPERATURE_MODIFY]*/)
	{
		new_temperature = TRUE;						/*"reacton_temperature" */
	}
	if (keycount[Keywords::KEY_MIX] > 0						||
		keycount[Keywords::KEY_MIX_RAW] > 0)	
	{
		new_mix = TRUE;								/*"mix" */
	}
	if (keycount[Keywords::KEY_SOLUTION] > 0 ||			
		keycount[Keywords::KEY_SOLUTION_SPREAD] > 0			||
		keycount[Keywords::KEY_SOLUTION_RAW] > 0			||
		keycount[Keywords::KEY_SOLUTION_MODIFY])
	{												/*"solution" */
		new_solution = TRUE;
	}
	if (keycount[Keywords::KEY_GAS_PHASE]  > 0				||
		keycount[Keywords::KEY_GAS_PHASE_RAW] > 0			||
		keycount[Keywords::KEY_GAS_PHASE_MODIFY])
	{
		new_gas_phase = TRUE;						/*"gas_phase" */
	}
	if (keycount[Keywords::KEY_SOLID_SOLUTIONS] > 0			||
		keycount[Keywords::KEY_SOLID_SOLUTIONS_RAW] > 0		||
		keycount[Keywords::KEY_SOLID_SOLUTIONS_MODIFY])
	{
		new_ss_assemblage = TRUE;					/*"solid_solutions" */
	}
	if (keycount[Keywords::KEY_KINETICS] > 0				/*||
		keycount[Keywords::KEY_KINETICS_RAW] > 0			||
		keycount[Keywords::KEY_KINETICS_MODIFY]*/)
	{
		new_kinetics = TRUE;						/*"kinetics" */
	}
	if (keycount[Keywords::KEY_INVERSE_MODELING] > 0)
	{
		new_inverse = TRUE;							/*"inverse_modeling" */
	}
	if (keycount[Keywords::KEY_SELECTED_OUTPUT] > 0			||		/*"selected_output" */
		keycount[Keywords::KEY_USER_PUNCH] > 0)						/*"user_punch" */
	{
		new_punch = TRUE;
	}

	if (keycount[Keywords::KEY_COPY] > 0)
	{
		new_copy = TRUE;							/*"copy" */
	}
	if (keycount[Keywords::KEY_PITZER] > 0)
	{
		new_pitzer = TRUE;							/*"pitzer" */
	}
	if (keycount[Keywords::KEY_NAMED_EXPRESSIONS] > 0)
	{
		new_named_logk = TRUE;						/*"named_log_k" */
	}

/*
 *   Sort arrays
 */

/* species */
	if (new_model == TRUE)
	{
		qsort(s, (size_t) count_s, (size_t) sizeof(struct species *), s_compare);

/* master species */
		qsort(master, (unsigned) count_master, sizeof(struct master *), master_compare);

/* elements */
		qsort(elements, (size_t) count_elements, (size_t) sizeof(struct element *), element_compare);
/* phases */
		qsort(phases, (size_t) count_phases, (size_t) sizeof(struct phase *), phase_compare);

	}

/* named_log_k */
	if (new_named_logk)
	{
		tidy_logk();
	}
/*
 *   Check pointers, write reactions for species
 */
	if (new_model)
	{
		sum_species_map.clear();

		tidy_species();

		tidy_phases();

		tidy_master_isotope();
/*
 *   calculate gfw of water, kg/mole
 */
		compute_gfw("H2O", &gfw_water);
		gfw_water *= 0.001;
	}
/*
 *   tidy surface data
 */
	if (new_model || new_surface)
	{
		tidy_surface();
	}
/*
 *   tidy inverse data
 */
	if (new_inverse)
	{
		tidy_inverse();
	}
/*
 *   tidy gas phase data
 */
	if (new_gas_phase)
	{
		tidy_gas_phase();
	}
/*
 *   tidy pp_assemblage data
 */
	if (new_model || new_pp_assemblage)
	{
		tidy_pp_assemblage();
	}
/*
 *   tidy ss_assemblage data
 */
	if (new_model || new_ss_assemblage)
	{
		tidy_ss_assemblage();
	}
/*
 *   tidy exchange data, after pp_assemblages
 */
	if (new_exchange)
	{
		tidy_exchange();
		tidy_min_exchange();
		tidy_kin_exchange();
	}
/*
 *   tidy surface data
 */
	if (new_surface)
	{
		tidy_min_surface();
		tidy_kin_surface();
	}
/*
 *   tidy solution isotope data
 */
	if (new_solution)
	{
		tidy_isotopes();
	}
	if (new_model)
	{
		tidy_isotope_ratios();
		tidy_isotope_alphas();
	}
/*
 *   Duplicate kinetics
 */
	if (new_kinetics)
	{
		std::map<int, cxxKinetics>::iterator it;
		for (it = Rxn_kinetics_map.begin(); it != Rxn_kinetics_map.end(); it++)
		{
			n_user = it->second.Get_n_user();
			last = it->second.Get_n_user_end();
			it->second.Set_n_user_end(n_user);
			Utilities::Rxn_copies(Rxn_kinetics_map, n_user, last);
		}
	}

/*
 *   Tidy pitzer information
 */
	if (pitzer_model && new_model)
	{
		pitzer_tidy();
	}
/*
 *   Tidy SIT information
 */
	if (sit_model && new_model)
	{
		sit_tidy();
	}
/*
 *   Tidy punch information
 */
	if (get_input_errors() == 0 && (new_punch || new_model))
	{
		tidy_punch();
	}
/*
 *   Tidy solution information
 */
	if (new_solution)
	{
		tidy_solutions();
	}

	/*      if (new_model || new_exchange || new_pp_assemblage || new_surface || new_gas_phase || new_kinetics) reset_last_model(); */
	if (new_model)
	{
		reset_last_model();
	}
/*
 *   make sure essential species are defined
 */
	//if (new_model)
	{
		if (s_h2o == NULL)
		{
			input_error++;
			//error_msg("H2O not defined.", STOP);
			error_msg("H2O not defined.", CONTINUE);
		}
		else
		{
			if (s_h2o->primary == NULL)
			{
				input_error++;
				error_msg("H2O, primary master species for O, not defined.",
					CONTINUE);
			}
			if (s_h2o->secondary == NULL)
			{
				input_error++;
				error_msg("H2O, secondary master species for O(-2), not defined.",
					CONTINUE);
			}
			if (s_h2o->type != H2O)
			{
				input_error++;
				error_msg("H2O can only be defined in SOLUTION_SPECIES.",
					CONTINUE);
			}
		}
		if (s_hplus == NULL && s_h3oplus == NULL)
		{
			input_error++;
			error_msg("Neither H+ nor H3O+ are defined in solution_species.",
				STOP);
		}
		else if (s_hplus == NULL && s_h3oplus != NULL)
		{
			s_hplus = s_h3oplus;
			s_h3oplus = NULL;
		}
		else if (s_hplus != NULL && s_h3oplus == NULL)
		{
		}
		else if (s_hplus != NULL && s_h3oplus != NULL)
		{
			input_error++;
			error_msg("Cannot define both H+ and H3O+ in solution_species.",
				STOP);
		}
		if (s_hplus->primary == NULL)
		{
			input_error++;
			error_msg("H3O+, primary master species for H, not defined.",
				CONTINUE);
		}
		if (s_hplus->secondary == NULL)
		{
			input_error++;
			error_msg("H3O+, secondary master species for H(1), not defined.",
				CONTINUE);
		}
		if (s_eminus == NULL)
		{
			input_error++;
			error_msg("e- not defined in solution_species.", CONTINUE);
		}
		if (s_eminus->primary == NULL)
		{
			input_error++;
			error_msg("e-, primary master species for E-, not defined.",
				CONTINUE);
		}
		if (pitzer_model == FALSE || pitzer_pe == TRUE)
		{
			if (s_h2 == NULL)
			{
				input_error++;
				error_msg("H2(aq) not defined in solution_species.", CONTINUE);
			}
			if (s_o2 == NULL)
			{
				input_error++;
				error_msg("O2(aq) not defined in solution_species.", CONTINUE);
			}
		}
		element_h_one = element_store("H(1)");
		if (element_h_one == NULL)
		{
			input_error++;
			error_msg("H(1) not defined in solution_master_species.", CONTINUE);
		}
	}
/*
 *   Error check, program termination
 */
	if (get_input_errors() > 0 || parse_error > 0)
	{
		error_msg("Calculations terminating due to input errors.", STOP);
	}

	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
check_species_input(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Check species data for completeness
 */
	int i;
	int return_value;

	return_value = OK;
	for (i = 0; i < count_s; i++)
	{
		if (s[i]->next_elt == NULL)
		{
			input_error++;
			return_value = ERROR;
			error_string = sformatf(
					"Elements in species have not been tabulated, %s.",
					s[i]->name);
			error_msg(error_string, CONTINUE);
		}
		if (s[i]->rxn == NULL)
		{
			input_error++;
			return_value = ERROR;
			error_string = sformatf(
					"Reaction for species has not been defined, %s.",
					s[i]->name);
			error_msg(error_string, CONTINUE);
		}
		else
		{
			select_log_k_expression(s[i]->logk, s[i]->rxn->logk);
			add_other_logk(s[i]->rxn->logk, s[i]->count_add_logk,
						   s[i]->add_logk);
		}
	}
	return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
select_log_k_expression(LDBLE * source_k, LDBLE * target_k)
/* ---------------------------------------------------------------------- */
{
	int j, analytic;

	analytic = FALSE;
	for (j = T_A1; j <= T_A6; j++)
	{
		if (source_k[j] != 0.0)
		{
			analytic = TRUE;
			break;
		}
	}
	if (analytic == TRUE)
	{
		target_k[logK_T0] = 0.0;
		target_k[delta_h] = 0.0;
		for (j = T_A1; j <= T_A6; j++)
		{
			target_k[j] = source_k[j];
		}
	}
	else
	{
		target_k[logK_T0] = source_k[logK_T0];
		target_k[delta_h] = source_k[delta_h];
		for (j = T_A1; j <= T_A6; j++)
		{
			target_k[j] = 0.0;
		}
	}
	for (j = delta_v; j < MAX_LOG_K_INDICES; j++)
	{
		target_k[j] = source_k[j];
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_logk(void)
/* ---------------------------------------------------------------------- */
/*
 *  Picks log k expression
 */
{
	int i;
	for (i = 0; i < count_logk; i++)
	{
		select_log_k_expression(logk[i]->log_k_original, logk[i]->log_k);
		logk[i]->done = FALSE;
	}
	for (i = 0; i < count_logk; i++)
	{
		if (logk[i]->done == FALSE)
		{
			add_logks(logk[i], 0);
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
add_other_logk(LDBLE * source_k, int count_add_logk,
			   struct name_coef *add_logk)
/* ---------------------------------------------------------------------- */
{
	int i, j, analytic;
	struct logk *logk_ptr;
	char token[MAX_LENGTH];
	LDBLE coef;
	ENTRY item, *found_item;

	if (count_add_logk == 0)
		return (OK);
	for (i = 0; i < count_add_logk; i++)
	{
		coef = add_logk[i].coef;
		strcpy(token, add_logk[i].name);
		str_tolower(token);
		item.key = token;
		item.data = NULL;
		found_item = hsearch_multi(logk_hash_table, item, FIND);
		if (found_item == NULL)
		{
			input_error++;
			error_string = sformatf(
					"Could not find named temperature expression, %s\n",
					token);
			error_msg(error_string, CONTINUE);
			return (ERROR);
		}
		logk_ptr = (struct logk *) found_item->data;
		analytic = FALSE;
		for (j = T_A1; j <= T_A6; j++)
		{
			if (logk_ptr->log_k[j] != 0.0)
			{
				analytic = TRUE;
				break;
			}
		}
		if (analytic == TRUE)
		{
			for (j = T_A1; j <= T_A6; j++)
			{
				source_k[j] += logk_ptr->log_k[j] * coef;
			}
		}
		else
		{
			source_k[logK_T0] += logk_ptr->log_k[logK_T0] * coef;
			source_k[delta_h] += logk_ptr->log_k[delta_h] * coef;
		}
		for (j = delta_v; j < MAX_LOG_K_INDICES; j++)
		{
			source_k[j] += logk_ptr->log_k[j] * coef;
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
add_logks(struct logk *logk_ptr, int repeats)
/* ---------------------------------------------------------------------- */
{
	int i, j;
	struct logk *next_logk_ptr;
	char token[MAX_LENGTH];
	LDBLE coef;
	ENTRY item, *found_item;
	/*
	 *  Adds in other named_expressions to get complete log K
	 *  Evaluates others recursively if necessary
	 */
	if (repeats > 15)
	{
		input_error++;
		error_string = sformatf( "Circular definition of named_logk? %s\n",
				logk_ptr->name);
		error_msg(error_string, CONTINUE);
		return (ERROR);
	}
	for (i = 0; i < logk_ptr->count_add_logk; i++)
	{
		coef = logk_ptr->add_logk[i].coef;
		strcpy(token, logk_ptr->add_logk[i].name);
		str_tolower(token);
		item.key = token;
		item.data = NULL;
		found_item = hsearch_multi(logk_hash_table, item, FIND);
		if (found_item == NULL)
		{
			input_error++;
			error_string = sformatf(
					"Could not find named temperature expression, %s\n",
					token);
			error_msg(error_string, CONTINUE);
			return (ERROR);
		}
		next_logk_ptr = (struct logk *) found_item->data;
		if (next_logk_ptr->done == FALSE)
		{
			/*output_msg(sformatf( "Done == FALSE\n", token)); */
			if (add_logks(next_logk_ptr, repeats + 1) == ERROR)
			{
				return (ERROR);
			}
		}
		for (j = 0; j < MAX_LOG_K_INDICES; j++)
		{
			logk_ptr->log_k[j] += next_logk_ptr->log_k[j] * coef;
		}
	}
	logk_ptr->done = TRUE;
	return (OK);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
coef_in_master(struct master * master_ptr)
/* ---------------------------------------------------------------------- */
{
	int l;
	LDBLE coef;
	char *ptr;
	char elt_name[MAX_LENGTH];
	struct elt_list *next_elt;

	coef = 0.0;
	char * temp_name = string_duplicate(master_ptr->elt->name);
	ptr = temp_name;
	get_elt(&ptr, elt_name, &l);
	free_check_null(temp_name);
	for (next_elt = master_ptr->s->next_elt; next_elt->elt != NULL;
		 next_elt++)
	{
		if (strcmp(elt_name, next_elt->elt->name) == 0)
		{
			coef = next_elt->coef;
			break;
		}
	}
	return (coef);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
rewrite_eqn_to_secondary(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Write equation for species in terms of secondary species
 *   Result is in trxn.
 */
	LDBLE coef;
	int repeat, i, add_count;
	struct rxn_token_temp *token_ptr;
/*
 *
 */
	add_count = 0;
	repeat = TRUE;
/*
 *   Reduce reaction equation to primary and secondary species
 */
	while (repeat == TRUE)
	{
		repeat = FALSE;
		/*   Check for too many iterations */
		if (++add_count > MAX_ADD_EQUATIONS)
		{
			parse_error++;
			error_string = sformatf(
					"Could not reduce equation to secondary master species, %s.",
					trxn.token[0].name);
			error_msg(error_string, CONTINUE);
			break;
		}

		for (i = 1; i < count_trxn; i++)
		{
			token_ptr = &(trxn.token[i]);
			if (token_ptr->s == NULL)
			{
				error_string = sformatf(
						"NULL species pointer for species, %s.",
						token_ptr->name);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			if (token_ptr->s->secondary == NULL
				&& token_ptr->s->primary == NULL)
			{
				coef = token_ptr->coef;
				trxn_add(token_ptr->s->rxn, coef, TRUE);
				repeat = TRUE;
				break;
			}
		}
	}
	trxn_combine();
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
replace_solids_gases(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Write equation for species in terms of secondary species
 *   Result is in trxn.
 */
	LDBLE coef;
	int n;
	int repeat, i, add_count;
	struct rxn_token_temp *token_ptr;
	struct phase *phase_ptr;
	int replaced;
	char token[MAX_LENGTH];
/*
 *
 */
	add_count = 0;
	repeat = TRUE;
	replaced = FALSE;
/*
 *   Reduce reaction equation to primary and secondary species
 */
	while (repeat == TRUE)
	{
		repeat = FALSE;
		/*   Check for too many iterations */
		if (++add_count > MAX_ADD_EQUATIONS)
		{
			parse_error++;
			error_string = sformatf(
					"Could not remove all solids and gases from equation, %s.",
					trxn.token[0].name);
			error_msg(error_string, CONTINUE);
			break;
		}

		for (i = 1; i < count_trxn; i++)
		{
			token_ptr = &(trxn.token[i]);
			if (token_ptr->s == NULL)
			{
				phase_ptr = phase_bsearch(token_ptr->name, &n, FALSE);
				/* try phase name without (g) or  (s) */
				if (phase_ptr == NULL)
				{
					strcpy(token, token_ptr->name);
					replace("(g)", "", token);
					replace("(s)", "", token);
					replace("(G)", "", token);
					replace("(S)", "", token);
					phase_ptr = phase_bsearch(token, &n, FALSE);
				}
				if (phase_ptr == NULL)
				{
					input_error++;
					error_string = sformatf( "Phase not found, %s.",
							token_ptr->name);
					error_msg(error_string, CONTINUE);
					break;
				}
				coef = token_ptr->coef;
				/* add reaction for solid/gas */
				/* debug
				   output_msg(sformatf( "Reaction to add.\n"));
				   rxn_print(phase_ptr->rxn);
				 */
				trxn_add_phase(phase_ptr->rxn, coef, FALSE);

				/* remove solid/gas from trxn list */
				trxn.token[i].name = phase_ptr->rxn->token[0].name;
				trxn.token[i].s = phase_ptr->rxn->token[0].s;
				trxn.token[i].coef = -coef * phase_ptr->rxn->token[0].coef;
				repeat = TRUE;
				replaced = TRUE;
				/* debug
				   output_msg(sformatf( "Before combined.\n"));
				   trxn_print();
				 */
				/* combine */
				trxn_combine();
				/* debug
				   output_msg(sformatf( "Combined.\n"));
				   trxn_print();
				 */
				break;
			}
		}
	}
	trxn_combine();
	return (replaced);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
rewrite_eqn_to_primary(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Write equation for secondary master species in terms of primary master species
 *   Store result in reaction structure for master species
 *   rewrite if necessary.
 *
 */
	int repeat, j, add_count;

/*
 *   Check secondary master species
 */
	repeat = TRUE;
	add_count = 0;
/*
 *   Check if reaction contains only primary master species
 */
	while (repeat == TRUE)
	{
		repeat = FALSE;
/*
 *   Check for too many iterations
 */
		if (++add_count > MAX_ADD_EQUATIONS)
		{
			parse_error++;
			error_string = sformatf(
					"Could not reduce equation to primary master species, %s.",
					trxn.token[0].s->name);
			error_msg(error_string, CONTINUE);
			break;
		}
/*
 *   Go through species in reaction for secondary master species, look for non-primary
 *   species as reactants, rewrite
 */
		for (j = 1; j < count_trxn; j++)
		{
			if (trxn.token[j].s->primary == NULL)
			{
				trxn_add(trxn.token[j].s->rxn, trxn.token[j].coef, TRUE);
				repeat = TRUE;
				break;
			}
		}
	}
	trxn_combine();
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_gas_phase(void)
/* ---------------------------------------------------------------------- */
{
	int n_user, last;
	LDBLE P, V_m;
	bool PR;
/*
 *   Find all gases for each gas_phase in phase list
 */
	for (std::set<int>::const_iterator nit = Rxn_new_gas_phase.begin(); nit != Rxn_new_gas_phase.end(); nit++)
	{
		std::map<int, cxxGasPhase>::iterator it = Rxn_gas_phase_map.find(*nit);
		if (it == Rxn_gas_phase_map.end())
		{
			assert(false);
		}
		cxxGasPhase *gas_phase_ptr = &(it->second);
		PR = false;
		P = 0.0;
		for (size_t j = 0; j < gas_phase_ptr->Get_gas_comps().size(); j++)
		{
			int k;
			struct phase *phase_ptr = phase_bsearch(gas_phase_ptr->Get_gas_comps()[j].Get_phase_name().c_str(), &k, FALSE);
			if (phase_ptr == NULL)
			{
				input_error++;
				error_string = sformatf(
					"Gas not found in PHASES database, %s.",
					gas_phase_ptr->Get_gas_comps()[j].Get_phase_name().c_str());
				error_msg(error_string, CONTINUE);
				continue;
			}
			else
			{
				if (phase_ptr->t_c > 0 && phase_ptr->p_c > 0)
					PR = true;
			}
		}
		gas_phase_ptr->Set_pr_in(PR);

		if (gas_phase_ptr->Get_new_def())
		{
			for (size_t j = 0; j < gas_phase_ptr->Get_gas_comps().size(); j++)
			{
				/*
				*   Fixed pressure
				*/
				if (gas_phase_ptr->Get_type() == cxxGasPhase::GP_PRESSURE)
				{
					if (gas_phase_ptr->Get_solution_equilibria())
					{
						input_error++;
						error_string = sformatf(
							"Gas phase %d: cannot use '-equilibrium' option with fixed pressure gas phase.",
							gas_phase_ptr->Get_n_user());
						error_msg(error_string, CONTINUE);
					}
					/* calculate moles */
					if (gas_phase_ptr->Get_gas_comps()[j].Get_p_read() != NAN)
					{
						P += gas_phase_ptr->Get_gas_comps()[j].Get_p_read();
						if (!PR)
							gas_phase_ptr->Get_gas_comps()[j].Set_moles(
							gas_phase_ptr->Get_gas_comps()[j].Get_p_read() * gas_phase_ptr->Get_volume() /
							R_LITER_ATM / gas_phase_ptr->Get_temperature());
					}
					else
					{
						input_error++;
						error_string = sformatf(
							"Gas phase %d: partial pressure of gas component %s not defined.",
							gas_phase_ptr->Get_n_user(), gas_phase_ptr->Get_gas_comps()[j].Get_phase_name().c_str());
						error_msg(error_string, CONTINUE);
					}
				}
				else
				{
					/*
					*   Fixed volume
					*/
					if (!gas_phase_ptr->Get_solution_equilibria())
					{
						if (gas_phase_ptr->Get_gas_comps()[j].Get_p_read() != NAN)
						{
							P += gas_phase_ptr->Get_gas_comps()[j].Get_p_read();
							if (!PR)
								gas_phase_ptr->Get_gas_comps()[j].Set_moles (
								gas_phase_ptr->Get_gas_comps()[j].Get_p_read() *
								gas_phase_ptr->Get_volume() / R_LITER_ATM /
								gas_phase_ptr->Get_temperature());
						}
						else
						{
							input_error++;
							error_string = sformatf(
								"Gas phase %d: moles of gas component %s not defined.",
								gas_phase_ptr->Get_n_user(),
								gas_phase_ptr->Get_gas_comps()[j].Get_phase_name().c_str());
							error_msg(error_string, CONTINUE);
						}
					}
				}
			}

			if (PR && P > 0)
			{
				std::vector<struct phase *> phase_ptrs;
				size_t j_PR;
				std::vector<cxxGasComp> &gc = gas_phase_ptr->Get_gas_comps();
				for (j_PR = 0; j_PR < gas_phase_ptr->Get_gas_comps().size(); j_PR++)
				{
					int k;
					struct phase *phase_ptr = phase_bsearch(gas_phase_ptr->Get_gas_comps()[j_PR].Get_phase_name().c_str(), &k, FALSE);
					if (gc[j_PR].Get_p_read() == 0)
						continue;
					if (phase_ptr)
					{
						phase_ptr->moles_x = gc[j_PR].Get_p_read() / P;
						phase_ptrs.push_back(phase_ptr);
					}
				}
				V_m = calc_PR(phase_ptrs, P, gas_phase_ptr->Get_temperature(), 0);
				gas_phase_ptr->Set_v_m(V_m);
				if (gas_phase_ptr->Get_type() == cxxGasPhase::GP_VOLUME)
				{
					gas_phase_ptr->Set_total_p(P);
				}
				for (j_PR = 0; j_PR < gas_phase_ptr->Get_gas_comps().size(); j_PR++)
				{
					int k;
					struct phase *phase_ptr = phase_bsearch(gc[j_PR].Get_phase_name().c_str(), &k, FALSE);
					if (gc[j_PR].Get_p_read() == 0)
					{
						gc[j_PR].Set_moles(0.0);
					} else
					{
						if (phase_ptr)
						{
							gc[j_PR].Set_moles(phase_ptr->moles_x *	gas_phase_ptr->Get_volume() / V_m);
							gas_phase_ptr->Set_total_moles(gas_phase_ptr->Get_total_moles() + gc[j_PR].Get_moles());
						}
					}
				}
			}
			/* 
			*   Duplicate gas phase, only if not solution equilibria
			*/
			if (!gas_phase_ptr->Get_solution_equilibria())
			{
				gas_phase_ptr->Set_new_def(false);
				n_user = gas_phase_ptr->Get_n_user();
				last = gas_phase_ptr->Get_n_user_end();
				gas_phase_ptr->Set_n_user_end(n_user);
				for (int j1 = n_user + 1; j1 <= last; j1++)
				{
					Utilities::Rxn_copy(Rxn_gas_phase_map, n_user, j1);
				}
			}
			else
			{
				gas_phase_ptr->Set_new_def(true);
			}
		}
	}
	return (OK);
}
#ifdef SKIP
/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_gas_phase(void)
/* ---------------------------------------------------------------------- */
{
	int n_user, last;
	LDBLE P, V_m;
	bool PR;
/*
 *   Find all gases for each gas_phase in phase list
 */
	for (std::set<int>::const_iterator nit = Rxn_new_gas_phase.begin(); nit != Rxn_new_gas_phase.end(); nit++)
	{
		std::map<int, cxxGasPhase>::iterator it = Rxn_gas_phase_map.find(*nit);
		if (it == Rxn_gas_phase_map.end())
		{
			assert(false);
		}
		cxxGasPhase *gas_phase_ptr = &(it->second);
		PR = false;
		P = 0.0;
		std::vector<cxxGasComp> gc = gas_phase_ptr->Get_gas_comps();
		for (size_t j = 0; j < gc.size(); j++)
		{
			int k;
			struct phase *phase_ptr = phase_bsearch(gc[j].Get_phase_name().c_str(), &k, FALSE);
			if (phase_ptr == NULL)
			{
				input_error++;
				error_string = sformatf(
						"Gas not found in PHASES database, %s.",
						gc[j].Get_phase_name().c_str());
				error_msg(error_string, CONTINUE);
				continue;
			}
			else
			{
				if (phase_ptr->t_c > 0 && phase_ptr->p_c > 0)
					PR = true;
			}
			gas_phase_ptr->Set_pr_in(PR);
			if (gas_phase_ptr->Get_new_def())
			{
				if (j == gc.size() - 1)
					gas_phase_ptr->Set_new_def(false);
				/*
				*   Fixed pressure
				*/
				if (gas_phase_ptr->Get_type() == cxxGasPhase::GP_PRESSURE)
				{
					if (gas_phase_ptr->Get_solution_equilibria())
					{
						input_error++;
						error_string = sformatf(
							"Gas phase %d: cannot use '-equilibrium' option with fixed pressure gas phase.",
							gas_phase_ptr->Get_n_user());
						error_msg(error_string, CONTINUE);
					}
					/* calculate moles */
					if (gc[j].Get_p_read() != NAN)
					{
						P += gc[j].Get_p_read();
						if (!PR)
							gc[j].Set_moles(
								gc[j].Get_p_read() * gas_phase_ptr->Get_volume() /
								R_LITER_ATM / gas_phase_ptr->Get_temperature());
					}
					else
					{
						input_error++;
						error_string = sformatf(
							"Gas phase %d: partial pressure of gas component %s not defined.",
							gas_phase_ptr->Get_n_user(), gc[j].Get_phase_name().c_str());
						error_msg(error_string, CONTINUE);
					}
				}
				else
				{
					/*
					*   Fixed volume
					*/
					if (!gas_phase_ptr->Get_solution_equilibria())
					{
						if (gc[j].Get_p_read() != NAN)
						{
							P += gc[j].Get_p_read();
							if (!PR)
								gc[j].Set_moles (
									gc[j].Get_p_read() *
									gas_phase_ptr->Get_volume() / R_LITER_ATM /
									gas_phase_ptr->Get_temperature());
						}
						else
						{
							input_error++;
							error_string = sformatf(
								"Gas phase %d: moles of gas component %s not defined.",
								gas_phase_ptr->Get_n_user(),
								gc[j].Get_phase_name().c_str());
							error_msg(error_string, CONTINUE);
						}
					}
				}


				gas_phase_ptr->Set_gas_comps(gc);

				if (PR && P > 0 && j == gc.size() - 1)
				{
					std::vector<struct phase *> phase_ptrs;
					size_t j_PR;
					for (j_PR = 0; j_PR < gas_phase_ptr->Get_gas_comps().size(); j_PR++)
					{
						int k;
						struct phase *phase_ptr = phase_bsearch(gas_phase_ptr->Get_gas_comps()[j_PR].Get_phase_name().c_str(), &k, FALSE);
						if (gc[j_PR].Get_p_read() == 0)
							continue;
						phase_ptr->moles_x = gc[j_PR].Get_p_read() / P;
						phase_ptrs.push_back(phase_ptr);
					}
					V_m = calc_PR(phase_ptrs, P, gas_phase_ptr->Get_temperature(), 0);
					gas_phase_ptr->Set_v_m(V_m);
					if (gas_phase_ptr->Get_type() == cxxGasPhase::GP_VOLUME)
					{
						gas_phase_ptr->Set_total_p(P);
					}
					std::vector<cxxGasComp> gc = gas_phase_ptr->Get_gas_comps();
					for (j_PR = 0; j_PR < gas_phase_ptr->Get_gas_comps().size(); j_PR++)
					{
						int k;
						struct phase *phase_ptr = phase_bsearch(gc[j_PR].Get_phase_name().c_str(), &k, FALSE);
						if (gc[j_PR].Get_p_read() == 0)
						{
							gc[j_PR].Set_moles(0.0);
						} else
						{
							gc[j_PR].Set_moles(phase_ptr->moles_x *
								gas_phase_ptr->Get_volume() / V_m);
							gas_phase_ptr->Set_total_moles(gas_phase_ptr->Get_total_moles() + gc[j_PR].Get_moles());
						}
					}
					gas_phase_ptr->Set_gas_comps(gc);
				}
				/* 
				*   Duplicate gas phase, only if not solution equilibria
				*/
				if (!gas_phase_ptr->Get_solution_equilibria())
				{
					n_user = gas_phase_ptr->Get_n_user();
					last = gas_phase_ptr->Get_n_user_end();
					gas_phase_ptr->Set_n_user_end(n_user);
					for (int j1 = n_user + 1; j1 <= last; j1++)
					{
						Utilities::Rxn_copy(Rxn_gas_phase_map, n_user, j1);
					}
				}
				else
				{
					gas_phase_ptr->Set_new_def(true);
				}
			}
		}
	}
	return (OK);
}
#endif
/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_inverse(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   After all of data are read, fill in data for an inverse structure,
 *   including master species pointers, phase pointers, and uncertainties
 *   and a list of all elements from phases or -balance input.
 */
	int i, j, k, l;
	int count_in;
	LDBLE value;
	struct inv_elts *inv_elts;
	struct master *master_ptr;
	struct master *master_alk_ptr;
	struct elt_list *elt_list_ptr;
	master_alk_ptr = master_bsearch("Alkalinity");
	for (i = 0; i < count_inverse; i++)
	{
		if (inverse[i].new_def != TRUE)
			continue;
/*
 *   Set default uncertainties for all solutions, if necessary
 */
		if (inverse[i].count_uncertainties < inverse[i].count_solns)
		{
			inverse[i].uncertainties =
				(LDBLE *) PHRQ_realloc(inverse[i].uncertainties,
									   (size_t) inverse[i].count_solns *
									   sizeof(LDBLE));
			if (inverse[i].uncertainties == NULL)
				malloc_error();
			for (j = inverse[i].count_uncertainties;
				 j < inverse[i].count_solns; j++)
			{
				inverse[i].uncertainties[j] =
					inverse[i].uncertainties[inverse[i].count_uncertainties -
											 1];
			}
		}
/*
 *   Set default ph uncertainties for all solutions, if necessary
 */
		if (inverse[i].count_ph_uncertainties < inverse[i].count_solns)
		{
			inverse[i].ph_uncertainties =
				(LDBLE *) PHRQ_realloc(inverse[i].ph_uncertainties,
									   (size_t) inverse[i].count_solns *
									   sizeof(LDBLE));
			if (inverse[i].ph_uncertainties == NULL)
				malloc_error();
			for (j = inverse[i].count_ph_uncertainties;
				 j < inverse[i].count_solns; j++)
			{
				inverse[i].ph_uncertainties[j] =
					inverse[i].ph_uncertainties[inverse[i].
												count_ph_uncertainties - 1];
			}
		}
/*
 *   Set default force for all solutions
 */
		if (inverse[i].count_force_solns < inverse[i].count_solns)
		{
			inverse[i].force_solns =
				(int *) PHRQ_realloc(inverse[i].force_solns,
									 (size_t) inverse[i].count_solns *
									 sizeof(int));
			if (inverse[i].force_solns == NULL)
				malloc_error();
			for (j = inverse[i].count_force_solns; j < inverse[i].count_solns;
				 j++)
			{
				inverse[i].force_solns[j] = FALSE;
			}
		}
/*
 *   Find master species for element, set uncertainties
 */
		for (j = 0; j < inverse[i].count_elts; j++)
		{
			inverse[i].elts[j].master =
				master_bsearch_primary(inverse[i].elts[j].name);
			if (inverse[i].elts[j].master == NULL)
			{
				input_error++;
				error_string = sformatf( "No master species for element, %s.",
						inverse[i].elts[j].name);
				error_msg(error_string, CONTINUE);
				continue;
			}
			inverse[i].elts[j].uncertainties =
				(LDBLE *) PHRQ_realloc(inverse[i].elts[j].uncertainties,
									   (size_t) inverse[i].count_solns *
									   sizeof(LDBLE));
			if (inverse[i].elts[j].uncertainties == NULL)
				malloc_error();
			if (inverse[i].elts[j].count_uncertainties == 0)
			{
/* use default uncertainties for element */
				for (k = 0; k < inverse[i].count_solns; k++)
				{
					inverse[i].elts[j].uncertainties[k] =
						inverse[i].uncertainties[k];
				}
			}
			else if (inverse[i].elts[j].count_uncertainties <
					 inverse[i].count_solns)
			{
/* use input uncertainties, fill in any missing at end */
				value =
					inverse[i].elts[j].uncertainties[inverse[i].elts[j].
													 count_uncertainties - 1];
				for (k = inverse[i].elts[j].count_uncertainties;
					 k < inverse[i].count_solns; k++)
				{
					inverse[i].elts[j].uncertainties[k] = value;
				}
			}
		}
/*
 *   Find phase
 */
		count_elts = 0;
		paren_count = 0;
		for (j = 0; j < inverse[i].count_phases; j++)
		{
			inverse[i].phases[j].phase =
				phase_bsearch(inverse[i].phases[j].name, &k, FALSE);
			if (inverse[i].phases[j].phase == NULL)
			{
				input_error++;
				error_string = sformatf( "Could not find phase, %s.",
						inverse[i].phases[j].name);
				error_msg(error_string, CONTINUE);
				continue;
			}
/*
 *   Find isotope elements
 */
			if (inverse[i].phases[j].count_isotopes > 0)
			{
				for (k = 0; k < inverse[i].phases[j].count_isotopes; k++)
				{
					inverse[i].phases[j].isotopes[k].primary = NULL;
					inverse[i].phases[j].isotopes[k].master = NULL;
					master_ptr =
						master_bsearch(inverse[i].phases[j].isotopes[k].
									   elt_name);
					if (master_ptr == NULL)
					{
						input_error++;
						error_string = sformatf(
								"Element not found for isotope calculation: %s.",
								inverse[i].phases[j].isotopes[k].elt_name);
						error_msg(error_string, CONTINUE);
						continue;
					}
					if (master_ptr->primary != TRUE)
					{
						input_error++;
						error_string = sformatf(
								"Isotope ratio may only be used"
								" for total element in phase.\n"
								"Secondary species not allowed: %s.",
								master_ptr->elt->name);
						error_msg(error_string, CONTINUE);
						continue;
					}
					inverse[i].phases[j].isotopes[k].primary = master_ptr;
					inverse[i].phases[j].isotopes[k].master = master_ptr;
					/* find coefficient for element */
					for (elt_list_ptr = inverse[i].phases[j].phase->next_elt;
						 elt_list_ptr->elt != NULL; elt_list_ptr++)
					{
						if (elt_list_ptr->elt == master_ptr->elt)
						{
							inverse[i].phases[j].isotopes[k].coef =
								elt_list_ptr->coef;
							break;
						}
					}
					if (elt_list_ptr == NULL)
					{
						input_error++;
						error_string = sformatf(
								"Element, %s,for which isotope ratio was defined is not found in phase, %s",
								master_ptr->elt->name,
								inverse[i].phases[j].phase->name);
						error_msg(error_string, CONTINUE);
						continue;
					}
				}
				qsort(inverse[i].phases[j].isotopes,
					  (size_t) inverse[i].phases[j].count_isotopes,
					  (size_t) sizeof(struct isotope), isotope_compare);
			}
			add_elt_list(inverse[i].phases[j].phase->next_elt, 1.0);

		}
		if (get_input_errors() > 0)
			return (ERROR);
/*
 *   Sort elements in reaction and combine
 */
		if (count_elts > 0)
		{
			qsort(elt_list, (size_t) count_elts,
				  (size_t) sizeof(struct elt_list), elt_list_compare);
			elt_list_combine();
		}
/*
 *   Mark master species list
 */
		for (j = 0; j < count_master; j++)
			master[j]->in = FALSE;
		for (j = 0; j < count_elts; j++)
		{
			elt_list[j].elt->master->in = TRUE;
		}
		/* Include all input elements */
		for (j = 0; j < inverse[i].count_elts; j++)
		{
			inverse[i].elts[j].master->in = TRUE;
		}
		s_eminus->primary->in = TRUE;	/* Include electrons */
		master_alk_ptr->in = TRUE;	/* Include alkalinity */
/*
 *   Unmark primary and mark secondary master species for redox elements
 */
		count_in = 0;
		inverse[i].count_redox_rxns = 0;
		for (j = 0; j < count_master; j++)
		{
			/*   skip all secondary master species in this loop */
			if (master[j]->primary == FALSE || master[j]->in == FALSE)
				continue;
			count_in++;
			if (j + 1 == count_master)
				continue;
			/*   if next master species is secondary, mark all 
			   secondary master species until a primary is found */
			if (master[j + 1]->primary == FALSE)
			{
				master[j]->in = FALSE;
				count_in--;
				for (k = j + 1; k < count_master; k++)
				{
					if (master[k]->primary == FALSE)
					{
						count_in++;
						master[k]->in = TRUE;
						if (master[k]->s->primary == NULL)
						{
							inverse[i].count_redox_rxns++;
						}
					}
					else
					{
						break;
					}
				}
			}
		}
/*
 *   Save list of master species in inv_elts structure
 */
		inv_elts =
			(struct inv_elts *) PHRQ_malloc((size_t) (count_in) *
											sizeof(struct inv_elts));
		if (inv_elts == NULL)
			malloc_error();
		count_in = 0;
		for (j = 0; j < count_master; j++)
		{
			/* skip H(1) and O(-2) */
			if (master[j]->s == s_hplus || master[j]->s == s_h2o)
				continue;
			if (master[j]->in == TRUE)
			{
				/* set master */
				inv_elts[count_in].master = master[j];
				/* alloc uncertainties and set default */
				inv_elts[count_in].uncertainties =
					(LDBLE *) PHRQ_malloc((size_t) inverse[i].count_solns *
										  sizeof(LDBLE));
				if (inv_elts[count_in].uncertainties == NULL)
					malloc_error();
				for (k = 0; k < inverse[i].count_solns; k++)
				{
					inv_elts[count_in].uncertainties[k] =
						inverse[i].uncertainties[k];
				}
				count_in++;
			}
		}
		if (s_co3->secondary->in == TRUE)
		{
			inverse[i].carbon = TRUE;
		}
		else
		{
			inverse[i].carbon = FALSE;
		}
/*
 *   copy in input uncertainties 
 */
		/* copy primary redox to all secondary redox */
		for (j = 0; j < inverse[i].count_elts; j++)
		{
			master_ptr = master_bsearch(inverse[i].elts[j].name);
			if (master_ptr == NULL)
			{
				input_error++;
				error_string = sformatf( "Element not found, %s.",
						inverse[i].elts[j].name);
				error_msg(error_string, CONTINUE);
				continue;
			}
			if (master_ptr->primary == FALSE
				|| master_ptr->s->secondary == NULL)
				continue;
			for (k = 0; k < count_in; k++)
			{
				if (master_ptr == inv_elts[k].master->elt->primary)
				{
					for (l = 0; l < inverse[i].count_solns; l++)
					{
						inv_elts[k].uncertainties[l] =
							inverse[i].elts[j].uncertainties[l];
					}
				}
			}
			inverse[i].elts[j].uncertainties =
				(LDBLE *) free_check_null(inverse[i].elts[j].uncertainties);
		}
		/* copy masters that are not primary redox */
		for (j = 0; j < inverse[i].count_elts; j++)
		{
			master_ptr = master_bsearch(inverse[i].elts[j].name);
			if (master_ptr == NULL)
			{
				input_error++;
				error_string = sformatf( "Element not found, %s.",
						inverse[i].elts[j].name);
				error_msg(error_string, CONTINUE);
				continue;
			}
			if (master_ptr->primary == TRUE
				&& master_ptr->s->secondary != NULL)
				continue;
			for (k = 0; k < count_in; k++)
			{
				if (master_ptr == inv_elts[k].master)
				{
					for (l = 0; l < inverse[i].count_solns; l++)
					{
						inv_elts[k].uncertainties[l] =
							inverse[i].elts[j].uncertainties[l];
					}
					break;
				}
			}
			inverse[i].elts[j].uncertainties =
				(LDBLE *) free_check_null(inverse[i].elts[j].uncertainties);
		}
/*
 *   replace elts in inverse struct
 */
		inverse[i].elts =
			(struct inv_elts *) free_check_null(inverse[i].elts);
		inverse[i].elts = inv_elts;
		inverse[i].count_elts = count_in;
		for (j = 0; j < inverse[i].count_elts; j++)
		{
/* debug
			output_msg(sformatf( "\t%d\t%s", j, inverse[i].elts[j].master->elt->name));
			for (k = 0; k < inverse[i].count_solns; k++) {
				output_msg(sformatf( "\t%f", inverse[i].elts[j].uncertainties[k]));
			}
			output_msg(sformatf("\n"));
 */
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_phases(void)
/* ---------------------------------------------------------------------- */
{
	int i;
	int replaced;
	/*
	 *  Fix log Ks first, so they can possibly be added to other phase equations
	 */
	for (i = 0; i < count_phases; i++)
	{
		select_log_k_expression(phases[i]->logk, phases[i]->rxn->logk);
		add_other_logk(phases[i]->rxn->logk, phases[i]->count_add_logk,
					   phases[i]->add_logk);
		phases[i]->rxn->token[0].name = phases[i]->name;
		phases[i]->rxn->token[0].s = NULL;
	}
	/*
	 *   Rewrite all phases to secondary species
	 */
	for (i = 0; i < count_phases; i++)
	{
		/*
		 *   Rewrite equation
		 */
		count_trxn = 0;
		trxn_add_phase(phases[i]->rxn, 1.0, FALSE);
		trxn.token[0].name = phases[i]->name;
		/* debug 
		   output_msg(sformatf( "%s PHASE.\n", phases[i]->name));
		   trxn_print();
		 */
		replaced = replace_solids_gases();
		phases[i]->replaced = replaced;
		/*  save rxn */
		/*
		rxn_free(phases[i]->rxn);
		phases[i]->rxn = rxn_alloc(count_trxn + 1);
		trxn_copy(phases[i]->rxn);
		*/
		/*  save rxn_s */
		trxn_reverse_k();
		rewrite_eqn_to_secondary();
		trxn_reverse_k();
		rxn_free(phases[i]->rxn_s);
		phases[i]->rxn_s = rxn_alloc(count_trxn + 1);
		trxn_copy(phases[i]->rxn_s);
		/*
		 *   Check equation
		 */
		if (phases[i]->check_equation == TRUE)
		{
			if (replaced == FALSE)
			{
				phase_rxn_to_trxn(phases[i], phases[i]->rxn);
			}
			else
			{
				phase_rxn_to_trxn(phases[i], phases[i]->rxn_s);
			}
			if (check_eqn(FALSE) == ERROR)
			{
				input_error++;
				error_string = sformatf(
						"Equation for phase %s does not balance.",
						phases[i]->name);
				error_msg(error_string, CONTINUE);
			}
		}
	}

	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_pp_assemblage(void)
/* ---------------------------------------------------------------------- */
{
	LDBLE coef;
	char *ptr;
/*
 *   Find pointers for pure phases
 */
	//std::map<int, cxxPPassemblage>::iterator it;
	//it = Rxn_pp_assemblage_map.begin();
	//for ( ; it != Rxn_pp_assemblage_map.end(); it++)
	//{
	//for (size_t nn = 0; nn < Rxn_new_pp_assemblage.size(); nn++)
	//{
		//std::map<int, cxxPPassemblage>::iterator kit = Rxn_pp_assemblage_map.find(Rxn_new_pp_assemblage[nn]);
	for (std::set<int>::const_iterator nit = Rxn_new_pp_assemblage.begin(); nit != Rxn_new_pp_assemblage.end(); nit++)
	{
		std::map<int, cxxPPassemblage>::iterator kit = Rxn_pp_assemblage_map.find(*nit);
		if (kit == Rxn_pp_assemblage_map.end())
		{
			assert(false);
		}
		//if (!kit->second.Get_new_def()) continue;
		cxxPPassemblage *pp_assemblage_ptr = &(kit->second);
		count_elts = 0;
		paren_count = 0;
		coef = 1.0;
		pp_assemblage_ptr->Set_new_def(false);

		std::map<std::string, cxxPPassemblageComp>::iterator it;
		it =  pp_assemblage_ptr->Get_pp_assemblage_comps().begin();
		for ( ; it != pp_assemblage_ptr->Get_pp_assemblage_comps().end(); it++)
		{
			int k;
			struct phase *phase_ptr = phase_bsearch(it->first.c_str(), &k, FALSE);
			if (phase_ptr == NULL)
			{
				input_error++;
				error_string = sformatf( "Phase not found in database, %s.",
						it->first.c_str());
				error_msg(error_string, CONTINUE);
				continue;
			}
			else
			{
				add_elt_list(phase_ptr->next_elt, coef);
			}
			if (it->second.Get_add_formula().size() > 0)
			{
				int first = count_elts;
				phase_ptr =	phase_bsearch(it->second.Get_add_formula().c_str(), &k, FALSE);
				if (phase_ptr != NULL)
				{
					it->second.Set_add_formula(phase_ptr->formula);
				}
				{
					char * temp_add = string_duplicate(it->second.Get_add_formula().c_str());
					ptr = temp_add;
					get_elts_in_species(&ptr, coef);
					free_check_null(temp_add);
				}
				/* check that all elements are in the database */
				for (int l = first; l < count_elts; l++)
				{
					if (elt_list[l].elt->master == NULL)
					{
						input_error++;
						error_string = sformatf(
								"Element \"%s\" in alternative phase for \"%s\" in EQUILIBRIUM_PHASES not found in database.",
								elt_list[l].elt->name,
								it->first.c_str());
						error_msg(error_string, CONTINUE);
					}
				}
			}
		}

/*
 *   Store list with all elements in phases and add formulae
 */
		cxxNameDouble nd = elt_list_NameDouble();
		pp_assemblage_ptr->Set_eltList(nd);
/*
 *   Duplicate pure phases if necessary
 */

		int n_user = pp_assemblage_ptr->Get_n_user();
		int n_user_end = pp_assemblage_ptr->Get_n_user_end();
		pp_assemblage_ptr->Set_n_user_end(n_user);
		Utilities::Rxn_copies(Rxn_pp_assemblage_map, n_user, n_user_end);

	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_ss_assemblage(void)
/* ---------------------------------------------------------------------- */
{
	struct phase *phase_ptr;
	LDBLE nb, nc, n_tot, xb, xc, dnb, dnc, l_a0, l_a1;
	LDBLE xb2, xb3, xb4, xc2, xc3;
	LDBLE moles;
/*
 *   Find pointers for pure phases
 */
	//std::map<int, cxxSSassemblage>::iterator it;
	//for (it = Rxn_ss_assemblage_map.begin(); it != Rxn_ss_assemblage_map.end(); it++)
	//{
	//for (size_t nn = 0; nn < Rxn_new_ss_assemblage.size(); nn++)
	//{
		//std::map<int, cxxSSassemblage>::iterator it = Rxn_ss_assemblage_map.find(Rxn_new_ss_assemblage[nn]);
	for (std::set<int>::const_iterator nit = Rxn_new_ss_assemblage.begin(); nit != Rxn_new_ss_assemblage.end(); nit++)
	{
		std::map<int, cxxSSassemblage>::iterator it = Rxn_ss_assemblage_map.find(*nit);
		if (it == Rxn_ss_assemblage_map.end())
		{
			assert(false);
		}
		//if (!it->second.Get_new_def()) continue;
		count_elts = 0;
		paren_count = 0;
		cxxSSassemblage *ss_assemblage_ptr = &(it->second);
		std::vector<cxxSS *> ss_ptrs = ss_assemblage_ptr->Vectorize();
		for (size_t j = 0; j < ss_ptrs.size(); j++)
		{
			cxxSS *ss_ptr = ss_ptrs[j];
			for (size_t k = 0; k < ss_ptr->Get_ss_comps().size(); k++)
			{
				cxxSScomp * comp_ptr = &(ss_ptr->Get_ss_comps()[k]);
				int k1;
				phase_ptr =	phase_bsearch(comp_ptr->Get_name().c_str(), &k1, FALSE);
				if (phase_ptr == NULL)
				{
					input_error++;
					error_string = sformatf(
							"Phase not found in database, %s, assemblage %d.",
							comp_ptr->Get_name().c_str(),
							ss_assemblage_ptr->Get_n_user());
					error_msg(error_string, CONTINUE);
					continue;
				}
				else
				{
					phase_ptr->moles_x = 0;
					phase_ptr->fraction_x = 0;
				}
				if (comp_ptr->Get_moles() == NAN)
				{
					input_error++;
					error_string = sformatf(
							"Moles of solid solution component not defined, %s, assemblage %d.",
							comp_ptr->Get_name().c_str(),
							ss_assemblage_ptr->Get_n_user());
					error_msg(error_string, CONTINUE);
					continue;
				}
			}

			if (ss_assemblage_ptr->Get_new_def())
			{
				/*
				 *  Calculate a0 and a1 first
				 */
				ss_calc_a0_a1(ss_ptr);
				
				n_tot = 0;
				for (size_t k = 0; k < ss_ptr->Get_ss_comps().size(); k++)
				{
					cxxSScomp * comp_ptr = &(ss_ptr->Get_ss_comps()[k]);
					moles = comp_ptr->Get_moles();
					if (moles <= 0.0)
					{
						moles = MIN_TOTAL_SS;
						comp_ptr->Set_initial_moles(moles);
					}
					n_tot += moles;
				}

				for (size_t k = 0; k < ss_ptr->Get_ss_comps().size(); k++)
				{
					cxxSScomp * comp_ptr = &(ss_ptr->Get_ss_comps()[k]);
					moles = comp_ptr->Get_moles();
					if (moles <= 0.0)
					{
						moles = MIN_TOTAL_SS;
					}
					comp_ptr->Set_fraction_x(moles / n_tot);
					comp_ptr->Set_log10_fraction_x(log10(moles / n_tot));
				}
				l_a0 = ss_ptr->Get_a0();
				l_a1 = ss_ptr->Get_a1();

/*
 *   Binary solid solution
 */
				if (l_a0 != 0.0 || l_a1 != 0)
				{
					ss_ptr->Set_dn(1.0 / n_tot);
					nc = ss_ptr->Get_ss_comps()[0].Get_moles();
					if (nc == 0)
						nc = MIN_TOTAL_SS;
					nb = ss_ptr->Get_ss_comps()[1].Get_moles();
					if (nb == 0)
						nb = MIN_TOTAL_SS;
					xc = nc / n_tot;
					xb = nb / n_tot;

					/* lambdas */
					ss_ptr->Get_ss_comps()[0].Set_log10_lambda(xb * xb * (l_a0 - l_a1 * (3 - 4 * xb)) / LOG_10);
					ss_ptr->Get_ss_comps()[1].Set_log10_lambda(xc * xc * (l_a0 + l_a1 * (4 * xb - 1)) / LOG_10);

					/* derivatives wrt nc and nb */
					xc2 = xc * xc;
					xc3 = xc2 * xc;
					xb2 = xb * xb;
					xb3 = xb2 * xb;
					xb4 = xb3 * xb;

					/* component 1 */
					dnb =
						-2 * l_a0 * xb * xc2 - 8 * l_a1 * xb2 * xc2 +
						6 * l_a1 * xb * xc2 - 4 * l_a1 * xc * xb4 -
						8 * l_a1 * xb3 * xc2 - 4 * l_a1 * xb2 * xc3 -
						2 * l_a0 * xc * xb2 - 8 * l_a1 * xc * xb3 +
						6 * l_a1 * xc * xb2 + 1;
					ss_ptr->Get_ss_comps()[0].Set_dnb(dnb / n_tot);
					dnc =
						2 * l_a0 * xb3 + 2 * l_a0 * xc * xb2 + 8 * l_a1 * xb4 +
						8 * l_a1 * xc * xb3 - 2 * l_a1 * xb3 - 6 * l_a1 * xc * xb2;
					ss_ptr->Get_ss_comps()[0].Set_dnc(-xb / nc + dnc / n_tot);
					ss_ptr->Get_ss_comps()[0].Set_dn(1.0 / n_tot);

					/* component 2 */
					dnb =
						2 * l_a0 * xb * xc2 + 2 * l_a0 * xc3 +
						8 * l_a1 * xb2 * xc2 + 8 * l_a1 * xb * xc3 -
						2 * l_a1 * xb * xc2 - 6 * l_a1 * xc3;
					ss_ptr->Get_ss_comps()[1].Set_dnb(-xc / nb + dnb / n_tot);
					dnc =
						-2 * l_a0 * xc * xb2 - 8 * l_a1 * xc * xb3 +
						2 * l_a1 * xc * xb2 - 2 * l_a0 * xb * xc2 -
						8 * l_a1 * xb2 * xc2 + 6 * l_a1 * xb * xc2 + 1;
					ss_ptr->Get_ss_comps()[1].Set_dnc(dnc / n_tot);
					ss_prep(ss_ptr->Get_tk(), ss_ptr, TRUE);
					ss_ptr->Get_ss_comps()[1].Set_dn(1.0 / n_tot);
/*
 *   Ideal solid solution
 */
				}
				else
				{
					ss_ptr->Set_dn(1.0 / n_tot);
					for (size_t k = 0; k < ss_ptr->Get_ss_comps().size(); k++)
					{
						cxxSScomp * comp_ptr = &(ss_ptr->Get_ss_comps()[k]);
						comp_ptr->Set_log10_lambda(0);
						moles = comp_ptr->Get_moles();
						if (moles <= 0.0)
							moles = MIN_TOTAL_SS;
						comp_ptr->Set_dnb((n_tot - moles) / (moles * n_tot));
						comp_ptr->Set_dn(1.0 / n_tot);
					}
				}
			}
		}
		ss_assemblage_ptr->Set_new_def(false);
/*
 *   Duplicate ss_assemblage if necessary
 */
		int n_user = ss_assemblage_ptr->Get_n_user();
		int n_user_end = ss_assemblage_ptr->Get_n_user_end();
		Utilities::Rxn_copies(Rxn_ss_assemblage_map, n_user, n_user_end);
		ss_assemblage_ptr->Set_n_user_end(n_user);
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_punch(void)
/* ---------------------------------------------------------------------- */
{
	//int i, j, l;
	int punch_save;
	//char token[MAX_LENGTH];
/*
 *   tidy punch information
 */
	std::map < int, SelectedOutput >::iterator so_it = SelectedOutput_map.begin(); 
	for ( ; so_it != SelectedOutput_map.end(); so_it++)
	{
		current_selected_output = &(so_it->second);
		if (current_selected_output == NULL)
			continue;


		/* totals */

		for (size_t i = 0; i < current_selected_output->Get_totals().size(); i++)
		{
			std::pair< std::string, void *> &pair_ptr = current_selected_output->Get_totals()[i];
			pair_ptr.second = master_bsearch(pair_ptr.first.c_str());
		}

		/* molalities */
		for (size_t i = 0; i < current_selected_output->Get_molalities().size(); i++)
		{
			std::pair< std::string, void *> &pair_ptr = current_selected_output->Get_molalities()[i];
			pair_ptr.second = s_search(pair_ptr.first.c_str());
		}

		/* log activities */

		//for (i = 0; i < punch.count_activities; i++)
		for (size_t i = 0; i < current_selected_output->Get_activities().size(); i++)
		{
			std::pair< std::string, void *> &pair_ptr = current_selected_output->Get_activities()[i];
			pair_ptr.second = s_search(pair_ptr.first.c_str());
		}

		/* equilibrium phases */

		//for (i = 0; i < punch.count_pure_phases; i++)
		for (size_t i = 0; i < current_selected_output->Get_pure_phases().size(); i++)
		{
			int j;
			std::pair< std::string, void *> &pair_ptr = current_selected_output->Get_pure_phases()[i];
			pair_ptr.second = phase_bsearch(pair_ptr.first.c_str(), &j, FALSE);
		}

		/* saturation indices */

		//for (i = 0; i < punch.count_si; i++)
		for (size_t i = 0; i < current_selected_output->Get_si().size(); i++)
		{
			int j;
			std::pair< std::string, void *> &pair_ptr = current_selected_output->Get_si()[i];
			pair_ptr.second = phase_bsearch(pair_ptr.first.c_str(), &j, FALSE);
		}

		/* gases */

		//for (i = 0; i < punch.count_gases; i++)
		for (size_t i = 0; i < current_selected_output->Get_gases().size(); i++)
		{
			int j;
			std::pair< std::string, void *> &pair_ptr = current_selected_output->Get_gases()[i];
			pair_ptr.second = phase_bsearch(pair_ptr.first.c_str(), &j, FALSE);
		}
	}
	/*
	*  Always write new headings when SELECTED_OUTPUT is read
	*/
	so_it = SelectedOutput_map.begin(); 
	for ( ; so_it != SelectedOutput_map.end(); so_it++)
	{
		current_selected_output = &(so_it->second);
		if (current_selected_output == NULL || 
			!current_selected_output->Get_new_def())
			continue;
		phrq_io->Set_punch_ostream(current_selected_output->Get_punch_ostream());


		int l;
		if (current_selected_output->Get_high_precision() == false)
		{
			l = 12;
		}
		else
		{
			l = 20;
		}
		// UserPunch
		std::map < int, UserPunch >::iterator up_it = UserPunch_map.find(current_selected_output->Get_n_user());
		current_user_punch = up_it == UserPunch_map.end() ? NULL : &(up_it->second);

		punch_save = pr.punch;
		pr.punch = TRUE;
		phrq_io->Set_punch_on(true);

		/* constant stuff, sim, pH, etc. */

		if (current_selected_output->Get_sim() == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "sim"));
		}
		if (current_selected_output->Get_state() == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "state"));
		}
		if (current_selected_output->Get_soln() == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "soln"));
		}
		if (current_selected_output->Get_dist() == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "dist_x"));
		}
		if (current_selected_output->Get_time() == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "time"));
		}
		if (current_selected_output->Get_step() == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "step"));
		}
		if (current_selected_output->Get_ph() == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "pH"));
		}
		if (current_selected_output->Get_pe() == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "pe"));
		}
		if (current_selected_output->Get_rxn() == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "reaction"));
		}
		if (current_selected_output->Get_temp() == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "temp"));
		}
		if (current_selected_output->Get_alk() == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "Alk"));
		}
		if (current_selected_output->Get_mu() == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "mu"));
		}
		if (current_selected_output->Get_water() == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "mass_H2O"));
		}
		if (current_selected_output->Get_charge_balance() == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "charge"));
		}
		if (current_selected_output->Get_percent_error() == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "pct_err"));
		}
		/* totals */

		//for (i = 0; i < punch.count_totals; i++)
		for (size_t i = 0; i < current_selected_output->Get_totals().size(); i++)
		{
			std::pair< std::string, void *> &pair_ref = current_selected_output->Get_totals()[i];
			fpunchf_heading(sformatf("%*s\t", l, pair_ref.first.c_str()));
			if (pair_ref.second == NULL)
			{
				error_string = sformatf( "Did not find master species,"
					" %s.", pair_ref.first.c_str());
				warning_msg(error_string);
			}
			//fpunchf_heading(sformatf("%*s\t", l, punch.totals[i].name));
			//if (punch.totals[i].master == NULL)
			//{
			//	error_string = sformatf( "Did not find master species,"
			//		" %s.", punch.totals[i].name);
			//	warning_msg(error_string);
			//}
		}

		/* molalities */

		//for (i = 0; i < punch.count_molalities; i++)
		for (size_t i = 0; i < current_selected_output->Get_molalities().size(); i++)
		{
			std::pair< std::string, void *> &pair_ref = current_selected_output->Get_molalities()[i];
			std::string name = "m_";
			name.append(pair_ref.first);
			fpunchf_heading(sformatf("%*s\t", l, name.c_str()));
			if (pair_ref.second == NULL)
			{
				error_string = sformatf( "Did not find species,"
					" %s.", pair_ref.first.c_str());
				warning_msg(error_string);
			}
			//strcpy(token, "m_");
			//strcat(token, punch.molalities[i].name);
			//fpunchf_heading(sformatf("%*s\t", l, token));
			//if (punch.molalities[i].s == NULL)
			//{
			//	error_string = sformatf( "Did not find species,"
			//		" %s.", punch.molalities[i].name);
			//	warning_msg(error_string);
			//}
		}

		/* log activities */

		//for (i = 0; i < punch.count_activities; i++)
		for (size_t i = 0; i < current_selected_output->Get_activities().size(); i++)
		{
			std::pair< std::string, void *> &pair_ref = current_selected_output->Get_activities()[i];
			std::string name = "la_";
			name.append(pair_ref.first);
			fpunchf_heading(sformatf("%*s\t", l, name.c_str()));
			if (pair_ref.second == NULL)
			{
				error_string = sformatf( "Did not find species,"
					" %s.", pair_ref.first.c_str());
				warning_msg(error_string);
			}
			//strcpy(token, "la_");
			//strcat(token, punch.activities[i].name);
			//fpunchf_heading(sformatf("%*s\t", l, token));
			//if (punch.activities[i].s == NULL)
			//{
			//	error_string = sformatf( "Did not find species, "
			//		"%s.", punch.activities[i].name);
			//	warning_msg(error_string);
			//}
		}

		/* equilibrium phases */

		//for (i = 0; i < punch.count_pure_phases; i++)
		for (size_t i = 0; i < current_selected_output->Get_pure_phases().size(); i++)
		{
			std::pair< std::string, void *> &pair_ref = current_selected_output->Get_pure_phases()[i];
			fpunchf_heading(sformatf("%*s\t", l, pair_ref.first.c_str()));
			std::string name = "d_";
			name.append(pair_ref.first);
			fpunchf_heading(sformatf("%*s\t", l, name.c_str()));
			if (pair_ref.second == NULL)
			{
				error_string = sformatf( "Did not find phase,"
					" %s.", pair_ref.first.c_str());
				warning_msg(error_string);
			}
			//strcpy(token, "d_");
			//strcat(token, punch.pure_phases[i].name);
			//fpunchf_heading(sformatf("%*s\t", l, punch.pure_phases[i].name));
			//fpunchf_heading(sformatf("%*s\t", l, token));
			//if (punch.pure_phases[i].phase == NULL)
			//{
			//	error_string = sformatf( "Did not find phase, "
			//		"%s.", punch.pure_phases[i].name);
			//	warning_msg(error_string);
			//}
		}

		/* saturation indices */

		//for (i = 0; i < punch.count_si; i++)
		for (size_t i = 0; i < current_selected_output->Get_si().size(); i++)
		{
			std::pair< std::string, void *> &pair_ref = current_selected_output->Get_si()[i];
			std::string name = "si_";
			name.append(pair_ref.first);
			fpunchf_heading(sformatf("%*s\t", l, name.c_str()));
			if (pair_ref.second == NULL)
			{
				error_string = sformatf( "Did not find phase,"
					" %s.", pair_ref.first.c_str());
				warning_msg(error_string);
			}
			//strcpy(token, "si_");
			//strcat(token, punch.si[i].name);
			//fpunchf_heading(sformatf("%*s\t", l, token));
			//if (punch.si[i].phase == NULL)
			//{
			//	error_string = sformatf( "Did not find phase, "
			//		"%s.", punch.si[i].name);
			//	warning_msg(error_string);
			//}
		}

		/* gases */

		//if (punch.count_gases > 0)
		if (current_selected_output->Get_gases().size() > 0)
		{
			fpunchf_heading(sformatf("%*s\t", l, "pressure"));
			fpunchf_heading(sformatf("%*s\t", l, "total mol"));
			fpunchf_heading(sformatf("%*s\t", l, "volume"));
		}
		//for (i = 0; i < punch.count_gases; i++)
		for (size_t i = 0; i < current_selected_output->Get_gases().size(); i++)
		{
			std::pair< std::string, void *> &pair_ref = current_selected_output->Get_gases()[i];
			std::string name = "g_";
			name.append(pair_ref.first);
			fpunchf_heading(sformatf("%*s\t", l, name.c_str()));
			if (pair_ref.second == NULL)
			{
				error_string = sformatf( "Did not find phase,"
					" %s.", pair_ref.first.c_str());
				warning_msg(error_string);
			}
			//strcpy(token, "g_");
			//strcat(token, punch.gases[i].name);
			//fpunchf_heading(sformatf("%*s\t", l, token));
			//if (punch.gases[i].phase == NULL)
			//{
			//	error_string = sformatf( "Did not find phase, "
			//		"%s.", punch.gases[i].name);
			//	warning_msg(error_string);
			//}
		}

		/* kinetics */

		//for (i = 0; i < punch.count_kinetics; i++)
		for (size_t i = 0; i < current_selected_output->Get_kinetics().size(); i++)
		{
			std::pair< std::string, void *> &pair_ref = current_selected_output->Get_kinetics()[i];
			std::string name = "k_";
			name.append(pair_ref.first);
			fpunchf_heading(sformatf("%*s\t", l, name.c_str()));
			name = "dk_";
			name.append(pair_ref.first);
			fpunchf_heading(sformatf("%*s\t", l, name.c_str()));
			//strcpy(token, "k_");
			//strcat(token, punch.kinetics[i].name);
			//fpunchf_heading(sformatf("%*s\t", l, token));
			//strcpy(token, "dk_");
			//strcat(token, punch.kinetics[i].name);
			//fpunchf_heading(sformatf("%*s\t", l, token));
		}

		/* solid solutions */

		//for (i = 0; i < punch.count_s_s; i++)
		for (size_t i = 0; i < current_selected_output->Get_s_s().size(); i++)
		{
			std::pair< std::string, void *> &pair_ref = current_selected_output->Get_s_s()[i];
			std::string name = "s_";
			name.append(pair_ref.first);
			fpunchf_heading(sformatf("%*s\t", l, name.c_str()));
			//strcpy(token, "s_");
			//strcat(token, punch.s_s[i].name);
			//fpunchf_heading(sformatf("%*s\t", l, token));
		}

		/* isotopes */

		//for (i = 0; i < punch.count_isotopes; i++)
		for (size_t i = 0; i < current_selected_output->Get_isotopes().size(); i++)
		{
			std::pair< std::string, void *> &pair_ref = current_selected_output->Get_isotopes()[i];
			if (isotope_ratio_search(pair_ref.first.c_str()) == NULL)
			{
				error_string = sformatf(
					"Did not find isotope_ratio definition for "
					"%s in -isotopes of SELECTED_OUTPUT.\n%s must be defined in ISOTOPE_RATIO data block.",
					pair_ref.first.c_str(), pair_ref.first.c_str());
				warning_msg(error_string);
			}
			std::string name = "I_";
			name.append(pair_ref.first);
			fpunchf_heading(sformatf("%*s\t", l, name.c_str()));
			//if (isotope_ratio_search(punch.isotopes[i].name) == NULL)
			//{
			//	error_string = sformatf(
			//		"Did not find isotope_ratio definition for "
			//		"%s in -isotopes of SELECTED_OUTPUT.\n%s must be defined in ISOTOPE_RATIO data block.",
			//		punch.isotopes[i].name, punch.isotopes[i].name);
			//	warning_msg(error_string);
			//}
			//strcpy(token, "I_");
			//strcat(token, punch.isotopes[i].name);
			//fpunchf_heading(sformatf("%*s\t", l, token));
		}

		/* calculate_values */

		//for (i = 0; i < punch.count_calculate_values; i++)
		for (size_t i = 0; i < current_selected_output->Get_calculate_values().size(); i++)
		{
			std::pair< std::string, void *> &pair_ref = current_selected_output->Get_calculate_values()[i];
			if (calculate_value_search(pair_ref.first.c_str()) == NULL)
			{
				error_string = sformatf(
					"Did not find calculate_values definition for "
					"%s in -calculate_values of SELECTED_OUTPUT.\n%s must be defined in CALCULATE_VALUES data block.",
					pair_ref.first.c_str(),
					pair_ref.first.c_str());
				warning_msg(error_string);
			}
			std::string name = "V_";
			name.append(pair_ref.first);
			fpunchf_heading(sformatf("%*s\t", l, name.c_str()));
			//if (calculate_value_search(punch.calculate_values[i].name) == NULL)
			//{
			//	error_string = sformatf(
			//		"Did not find calculate_values definition for "
			//		"%s in -calculate_values of SELECTED_OUTPUT.\n%s must be defined in CALCULATE_VALUES data block.",
			//		punch.calculate_values[i].name,
			//		punch.calculate_values[i].name);
			//	warning_msg(error_string);
			//}
			//strcpy(token, "V_");
			//strcat(token, punch.calculate_values[i].name);
			//fpunchf_heading(sformatf("%*s\t", l, token));
		}

		/* user_punch */
		if (current_user_punch != NULL && current_selected_output->Get_user_punch())
		{
			for (size_t i = 0; i < current_user_punch->Get_headings().size(); i++)
			{
				fpunchf_heading(sformatf("%*s\t", l, current_user_punch->Get_headings()[i].c_str()));
			}
		}
		fpunchf_heading("\n");
		//if (punch.user_punch == TRUE)
		//{
		//	for (i = 0; i < user_punch_count_headings; i++)
		//	{
		//		fpunchf_heading(sformatf("%*s\t", l, user_punch_headings[i]));
		//	}
		//}
		//fpunchf_heading("\n");

		current_selected_output->Set_new_def(false);
		pr.punch = punch_save;
		phrq_io->Set_punch_on(pr.punch == TRUE);

		punch_flush();
	}

	current_selected_output = NULL;
	current_user_punch = NULL;
	phrq_io->Set_punch_ostream(NULL);
	return (OK);
}
#ifdef SKIP
/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_punch(void)
/* ---------------------------------------------------------------------- */
{
	int i, j, l;
	int punch_save;
	char token[MAX_LENGTH];
/*
 *   tidy punch information
 */
	if (punch.high_precision == FALSE)
	{
		l = 12;
	}
	else
	{
		l = 20;
	}
	if (punch.in == TRUE)
	{
		/* totals */

		for (i = 0; i < punch.count_totals; i++)
		{
			punch.totals[i].master = master_bsearch(punch.totals[i].name);
		}

		/* molalities */

		for (i = 0; i < punch.count_molalities; i++)
		{
			punch.molalities[i].s = s_search(punch.molalities[i].name);
		}

		/* log activities */

		for (i = 0; i < punch.count_activities; i++)
		{
			punch.activities[i].s = s_search(punch.activities[i].name);
		}

		/* equilibrium phases */

		for (i = 0; i < punch.count_pure_phases; i++)
		{
			punch.pure_phases[i].phase =
				phase_bsearch(punch.pure_phases[i].name, &j, FALSE);
		}

		/* saturation indices */

		for (i = 0; i < punch.count_si; i++)
		{
			punch.si[i].phase = phase_bsearch(punch.si[i].name, &j, FALSE);
		}

		/* gases */

		for (i = 0; i < punch.count_gases; i++)
		{
			punch.gases[i].phase =
				phase_bsearch(punch.gases[i].name, &j, FALSE);
		}
	}
	/*
	 *  Always write new headings when SELECTED_OUTPUT is read
	 */
	if (punch.new_def == TRUE && punch.in == TRUE)
	{
		punch_save = pr.punch;
		pr.punch = TRUE;
		phrq_io->Set_punch_on(true);

		/* constant stuff, sim, pH, etc. */

		if (punch.sim == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "sim"));
		}
		if (punch.state == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "state"));
		}
		if (punch.soln == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "soln"));
		}
		if (punch.dist == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "dist_x"));
		}
		if (punch.time == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "time"));
		}
		if (punch.step == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "step"));
		}
		if (punch.ph == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "pH"));
		}
		if (punch.pe == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "pe"));
		}
		if (punch.rxn == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "reaction"));
		}
		if (punch.temp == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "temp"));
		}
		if (punch.alk == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "Alk"));
		}
		if (punch.mu == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "mu"));
		}
		if (punch.water == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "mass_H2O"));
		}
		if (punch.charge_balance == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "charge"));
		}
		if (punch.percent_error == TRUE)
		{
			fpunchf_heading(sformatf("%*s\t", l, "pct_err"));
		}
		/* totals */

		for (i = 0; i < punch.count_totals; i++)
		{
			fpunchf_heading(sformatf("%*s\t", l, punch.totals[i].name));
			if (punch.totals[i].master == NULL)
			{
				error_string = sformatf( "Did not find master species,"
						" %s.", punch.totals[i].name);
				warning_msg(error_string);
			}
		}

		/* molalities */

		for (i = 0; i < punch.count_molalities; i++)
		{
			strcpy(token, "m_");
			strcat(token, punch.molalities[i].name);
			fpunchf_heading(sformatf("%*s\t", l, token));
			if (punch.molalities[i].s == NULL)
			{
				error_string = sformatf( "Did not find species,"
						" %s.", punch.molalities[i].name);
				warning_msg(error_string);
			}
		}

		/* log activities */

		for (i = 0; i < punch.count_activities; i++)
		{
			strcpy(token, "la_");
			strcat(token, punch.activities[i].name);
			fpunchf_heading(sformatf("%*s\t", l, token));
			if (punch.activities[i].s == NULL)
			{
				error_string = sformatf( "Did not find species, "
						"%s.", punch.activities[i].name);
				warning_msg(error_string);
			}
		}

		/* equilibrium phases */

		for (i = 0; i < punch.count_pure_phases; i++)
		{
			strcpy(token, "d_");
			strcat(token, punch.pure_phases[i].name);
			fpunchf_heading(sformatf("%*s\t", l, punch.pure_phases[i].name));
			fpunchf_heading(sformatf("%*s\t", l, token));
			if (punch.pure_phases[i].phase == NULL)
			{
				error_string = sformatf( "Did not find phase, "
						"%s.", punch.pure_phases[i].name);
				warning_msg(error_string);
			}
		}

		/* saturation indices */

		for (i = 0; i < punch.count_si; i++)
		{
			strcpy(token, "si_");
			strcat(token, punch.si[i].name);
			fpunchf_heading(sformatf("%*s\t", l, token));
			if (punch.si[i].phase == NULL)
			{
				error_string = sformatf( "Did not find phase, "
						"%s.", punch.si[i].name);
				warning_msg(error_string);
			}
		}

		/* gases */

		if (punch.count_gases > 0)
		{
			fpunchf_heading(sformatf("%*s\t", l, "pressure"));
			fpunchf_heading(sformatf("%*s\t", l, "total mol"));
			fpunchf_heading(sformatf("%*s\t", l, "volume"));
		}
		for (i = 0; i < punch.count_gases; i++)
		{
			strcpy(token, "g_");
			strcat(token, punch.gases[i].name);
			fpunchf_heading(sformatf("%*s\t", l, token));
			if (punch.gases[i].phase == NULL)
			{
				error_string = sformatf( "Did not find phase, "
						"%s.", punch.gases[i].name);
				warning_msg(error_string);
			}
		}

		/* kinetics */

		for (i = 0; i < punch.count_kinetics; i++)
		{
			strcpy(token, "k_");
			strcat(token, punch.kinetics[i].name);
			fpunchf_heading(sformatf("%*s\t", l, token));
			strcpy(token, "dk_");
			strcat(token, punch.kinetics[i].name);
			fpunchf_heading(sformatf("%*s\t", l, token));
		}

		/* solid solutions */

		for (i = 0; i < punch.count_s_s; i++)
		{
			strcpy(token, "s_");
			strcat(token, punch.s_s[i].name);
			fpunchf_heading(sformatf("%*s\t", l, token));
		}

		/* isotopes */

		for (i = 0; i < punch.count_isotopes; i++)
		{
			if (isotope_ratio_search(punch.isotopes[i].name) == NULL)
			{
				error_string = sformatf(
						"Did not find isotope_ratio definition for "
						"%s in -isotopes of SELECTED_OUTPUT.\n%s must be defined in ISOTOPE_RATIO data block.",
						punch.isotopes[i].name, punch.isotopes[i].name);
				warning_msg(error_string);
			}
			strcpy(token, "I_");
			strcat(token, punch.isotopes[i].name);
			fpunchf_heading(sformatf("%*s\t", l, token));
		}

		/* calculate_values */

		for (i = 0; i < punch.count_calculate_values; i++)
		{
			if (calculate_value_search(punch.calculate_values[i].name) == NULL)
			{
				error_string = sformatf(
						"Did not find calculate_values definition for "
						"%s in -calculate_values of SELECTED_OUTPUT.\n%s must be defined in CALCULATE_VALUES data block.",
						punch.calculate_values[i].name,
						punch.calculate_values[i].name);
				warning_msg(error_string);
			}
			strcpy(token, "V_");
			strcat(token, punch.calculate_values[i].name);
			fpunchf_heading(sformatf("%*s\t", l, token));
		}

		/* user_punch */
		if (punch.user_punch == TRUE)
		{
			for (i = 0; i < user_punch_count_headings; i++)
			{
				fpunchf_heading(sformatf("%*s\t", l, user_punch_headings[i]));
			}
		}
		fpunchf_heading("\n");

		punch.new_def = FALSE;
		pr.punch = punch_save;
		phrq_io->Set_punch_on(pr.punch == TRUE);
	}
	punch_flush();
	return (OK);
}
#endif
/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_species(void)
/* ---------------------------------------------------------------------- */
{
	int i, j;
	struct master *master_ptr;
	char c, *ptr;
/*
 *   Make sure species pointers are ok
 */
	if (check_species_input() == ERROR)
	{
		error_msg("Calculations terminating due to input errors.", STOP);
	}
/*
 *   Set secondary and primary pointers in species structures
 */
	for (i = 0; i < count_s; i++)
	{
		s[i]->number = i;
		s[i]->primary = NULL;
		s[i]->secondary = NULL;
		if (s[i]->check_equation == TRUE)
		{
			species_rxn_to_trxn(s[i]);
			if (check_eqn(TRUE) == ERROR)
			{
				input_error++;
				error_string = sformatf(
						"Equation for species %s does not balance.",
						s[i]->name);
				error_msg(error_string, CONTINUE);
			}
		}
	}
	for (i = 0; i < count_master; i++)
	{
		char * temp_name = string_duplicate(master[i]->elt->name);
		ptr = temp_name;
		if (ptr[0] != '[')
		{
			while ((c = (int) *(++ptr)) != '\0')
			{
				if (isupper((int) c))
				{
					input_error++;
					error_string = sformatf(
							"Element or valence name in SOLUTION_MASTER_SPECIES should include only one element, %s.",
							master[i]->elt->name);
					error_msg(error_string, CONTINUE);
					break;
				}
			}
		}
		free_check_null(temp_name);
		/* store sequence number in master structure */
		master[i]->number = i;
		if (strcmp(master[i]->elt->name, "Alkalinity") != 0)
		{
			if (master[i]->primary == TRUE)
			{
				master[i]->s->primary = master[i];
			}
			else
			{
				master[i]->s->secondary = master[i];
			}
		}
		if (strcmp(master[i]->elt->name, "C") == 0)
		{
			s_co3 = master[i]->s;
		}
		if (master[i]->gfw_formula != NULL)
		{
			if (compute_gfw(master[i]->gfw_formula, &master[i]->gfw) == ERROR)
			{
				input_error++;
				error_string = sformatf(
						"Calculating gfw for master species, %s, formula %s.",
						master[i]->elt->name, master[i]->gfw_formula);
				error_msg(error_string, CONTINUE);
			}
		}
	}
/*
 *   Write equations for all master species in terms of primary
 *   master species, set coefficient of element in master species
 */
	for (i = 0; i < count_master; i++)
	{
		count_trxn = 0;
		if (master[i]->s->primary != NULL)
		{
			trxn_add(master[i]->s->rxn, 1.0, FALSE);
			trxn_add(master[i]->s->rxn, -1.0, TRUE);
		}
		else
		{
			trxn_add(master[i]->s->rxn, 1.0, FALSE);
			rewrite_eqn_to_primary();
		}
		rxn_free(master[i]->rxn_primary);
		master[i]->rxn_primary = rxn_alloc(count_trxn + 1);
		trxn_copy(master[i]->rxn_primary);
		master[i]->coef = coef_in_master(master[i]);
	}
/*
 *   Rewrite all species to secondary species
 */
	for (i = 0; i < count_s; i++)
	{
		count_trxn = 0;
		if (s[i]->primary != NULL || s[i]->secondary != NULL)
		{
			trxn_add(s[i]->rxn, 1.0, FALSE);
			trxn_add(s[i]->rxn, -1.0, TRUE);
		}
		else
		{
			trxn_add(s[i]->rxn, 1.0, FALSE);
			rewrite_eqn_to_secondary();
		}
		rxn_free(s[i]->rxn_s);
		s[i]->rxn_s = rxn_alloc(count_trxn + 1);
		trxn_copy(s[i]->rxn_s);
		/* calculate alkalinity */
		s[i]->alk = calc_alk(s[i]->rxn_s);
		/* set co2 coefficient */
		s[i]->co2 = 0.0;
		for (j = 1; j < count_trxn; j++)
		{
			if (trxn.token[j].s == s_co3)
			{
				s[i]->co2 = trxn.token[j].coef;
				break;
			}
		}
	}
/*
 *   Set pointer in element to master species
 */
	for (i = 0; i < count_elements; i++)
	{
		elements[i]->master = master_bsearch(elements[i]->name);
		if (elements[i]->master == NULL)
		{
			input_error++;
			error_string = sformatf( "No master species for element %s.",
					elements[i]->name);
			error_msg(error_string, CONTINUE);
		}
		elements[i]->primary = master_bsearch_primary(elements[i]->name);
		if (elements[i]->primary == NULL)
		{
			input_error++;
			error_string = sformatf( "No master species for element %s.",
					elements[i]->name);
			error_msg(error_string, CONTINUE);
		}
	}
/*
 *   Make sure all primary master species for redox elements
 *   are also secondary master species
 */
	for (i = 0; i < count_master; i++)
	{
		if (master[i]->primary == FALSE)
		{
			master_ptr = master[i]->s->secondary->elt->primary;
			if (master_ptr == NULL)
			{
				input_error++;
				error_string = sformatf(
						"Every primary master species for a redox element\n"
						"\tmust also be a secondary master species.\n"
						"\tError in definitions related to %s .\n",
						master[i]->s->name);
				error_msg(error_string, CONTINUE);

			}
			else if (master_ptr->s->secondary == NULL)
			{
				input_error++;
				error_string = sformatf(
						"Every primary master species for a redox element\n"
						"\tmust also be a secondary master species.\n"
						"\t%s is the primary master species for element %s.\n"
						"\tAnother entry in SOLUTION_MASTER_SPECIES is needed.\n"
						"\tDefine species %s as a secondary master species for a valence state.\n"
						"\tFor example: \n" "\t%s(0)\t%s alk gfw",
						master_ptr->s->name, master_ptr->elt->name,
						master_ptr->s->name, master_ptr->elt->name,
						master_ptr->s->name);
				error_msg(error_string, CONTINUE);
			}
		}
	}
/*
 *   Calculate H and O if alternate mass balance is given
 */
	for (i = 0; i < count_s; i++)
	{
		if (s[i]->next_secondary != NULL)
		{
			s[i]->h = 0.0;
			s[i]->o = 0.0;
			for (j = 0; s[i]->next_secondary[j].elt != NULL; j++)
			{
				if (s[i]->next_secondary[j].elt->primary == NULL)
					continue;
				if (s[i]->next_secondary[j].elt->primary->s == s_hplus || s[i]->next_secondary[j].elt->primary->s == s_h3oplus)
				{
					s[i]->h += s[i]->next_secondary[j].coef;
				}
				else if (s[i]->next_secondary[j].elt->primary->s == s_h2o)
				{
					s[i]->o += s[i]->next_secondary[j].coef;
				}
				else if (s[i]->mole_balance != NULL)
				{
					master_ptr = s[i]->next_secondary[j].elt->master;
					if (master_ptr != NULL)
					{
						if (master_ptr->primary == TRUE)
						{
							if (master_ptr->s->secondary != NULL)
							{
								master_ptr = master_ptr->s->secondary;
							}
						}
					}
					else
					{
						input_error++;
						error_string = sformatf(
							"Element in -mole_balance %s not defined for species %s.\n", s[i]->mole_balance, s[i]->name);
						error_msg(error_string, CONTINUE);
						continue;
					}
					if (master_ptr->coef != 1)
					{
						s[i]->next_secondary[j].coef /= master_ptr->coef;
					}
				}
			}
			if (s[i]->type == EX)
			{
				for (j = 0; s[i]->next_secondary[j].elt != NULL; j++)
				{
					if (s[i]->next_secondary[j].elt->primary->s->type == EX)
					{
						s[i]->equiv = s[i]->next_secondary[j].coef;
						break;
					}
				}
			}
		}
		if (s[i]->type == EX)
		{
/*
 *   Find valence of cation from coefficients of reaction components
 *   Changed to be coefficient of exchanger
 */
			LDBLE exchange_coef = 0.0;
			for (j = 1; s[i]->rxn_s->token[j].s != NULL; j++)
			{
				if (s[i]->rxn_s->token[j].s->type == EX)
				{
					exchange_coef = s[i]->rxn_s->token[j].coef;
					break;
				}
			}
			if (exchange_coef == 0.0)
			{
				input_error++;
				error_string = sformatf(
					"No exchange species found in equation for %s.\n", s[i]->name);
				error_msg(error_string, CONTINUE);
				continue;
			}
			s[i]->equiv = exchange_coef;
		}
		if (s[i]->type == SURF)
		{
			LDBLE surface_coef = 0.0;
			/*
			 *   Find coefficient of surface in rxn, store in equiv
			 */
			for (j = 1; s[i]->rxn_s->token[j].s != NULL; j++)
			{
				if (s[i]->rxn_s->token[j].s->type == SURF)
				{
					surface_coef = s[i]->rxn_s->token[j].coef;
					break;
				}
			}
			if (surface_coef == 0.0)
			{
				input_error++;
				error_string = sformatf(
					"No surface species found in equation for %s.\n", s[i]->name);
				error_msg(error_string, CONTINUE);
				continue;
			}
			s[i]->equiv = surface_coef;
		}
	}
	
	for (i = 0; i < count_master; i++)
	{
		if (master[i]->gfw <= 0.0)
		{
			if (master[i]->type >= EMINUS) continue;
			if ((strcmp(master[i]->elt->name, "E") != 0) && 
			    (strcmp(master[i]->elt->name, "e") != 0) &&
			    (strcmp(master[i]->elt->name, "H(1)") != 0) &&
			    (strcmp(master[i]->elt->name, "O(-2)") != 0)
			    )
			{
				input_error++;
				error_string = sformatf(
					"Gram formula wt in SOLUTION_MASTER_SPECIES should not be <= 0.0, %s.\n", master[i]->elt->name);
				error_msg(error_string, CONTINUE);
			}
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_surface(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   After all of data are read, fill in master species for surface comps
 *   Sort surface
 */
	char *ptr1;
	cxxSurface *surface_ptr;
	//std::map<int, cxxSurface>::iterator kit;
	//for (kit = Rxn_surface_map.begin(); kit != Rxn_surface_map.end(); kit++)
	//{
	//for (size_t nn = 0; nn < Rxn_new_surface.size(); nn++)
	//{
		//std::map<int, cxxSurface>::iterator kit = Rxn_surface_map.find(Rxn_new_surface[nn]);
	for (std::set<int>::const_iterator nit = Rxn_new_surface.begin(); nit != Rxn_new_surface.end(); nit++)
	{
		std::map<int, cxxSurface>::iterator kit = Rxn_surface_map.find(*nit);
		if (kit == Rxn_surface_map.end())
		{
			assert(false);
		}
		//if (!kit->second.Get_new_def()) continue;
		surface_ptr = &(kit->second);
		// ccm incompatible with Donnan or diffuse_layer
		if (surface_ptr->Get_type() == cxxSurface::CCM)
		{
			if (surface_ptr->Get_dl_type() == cxxSurface::BORKOVEK_DL || surface_ptr->Get_dl_type() == cxxSurface::DONNAN_DL)
			{
					input_error++;
					error_string = "Cannot use -diffuse_layer or -donnan calculation with Constant Capacity Model.";
					error_msg(error_string, CONTINUE);
					continue;
			}
		}
		for (size_t i = 0; i < surface_ptr->Get_surface_comps().size(); i++)
		{
			cxxSurfaceComp *comp_ptr = &(surface_ptr->Get_surface_comps()[i]);
/*
 *   Find master species for each surface
 */
			cxxNameDouble::iterator jit = comp_ptr->Get_totals().begin();
			for ( ; jit != comp_ptr->Get_totals().end(); jit++ )
			{
				struct element *elt_ptr = element_store(jit->first.c_str());
				struct master *master_ptr = elt_ptr->master;
				if (master_ptr == NULL)
				{
					input_error++;
					error_string = sformatf(
							"Master species not in database for %s, "
							"skipping element.",
							elt_ptr->name);
					error_msg(error_string, CONTINUE);
					continue;
				}
				if (master_ptr->type != SURF)
					continue;
				comp_ptr->Set_master_element(elt_ptr->name);
/*
 *   Set flags
 */
				cxxSurfaceCharge *charge_ptr = surface_ptr->Find_charge(comp_ptr->Get_charge_name());
				/*
				 * Calculate moles of sites
				 */
				if (surface_ptr->Get_new_def()
					&& surface_ptr->Get_sites_units() == cxxSurface::SITES_DENSITY
					&& comp_ptr->Get_phase_name().size() == 0)
				{
					if (charge_ptr == NULL)
					{
						input_error++;
						error_string = sformatf(
								"Surface type is incompatible with site units for %s.",
								comp_ptr->Get_formula().c_str());
						error_msg(error_string, CONTINUE);
						continue;
					}
					comp_ptr->Set_moles(
						comp_ptr->Get_moles() * 1.0e18 *
						charge_ptr->Get_specific_area() *
						charge_ptr->Get_grams() / AVOGADRO);
					/*
					 *  Calculate totals
					 */
					count_elts = 0;
					paren_count = 0;
					{
						char * temp_formula = string_duplicate(comp_ptr->Get_formula().c_str());
						ptr1 = temp_formula;
						get_elts_in_species(&ptr1, comp_ptr->Get_moles());
						free_check_null(temp_formula);
					}
					{
						cxxNameDouble nd = elt_list_NameDouble();
						comp_ptr->Set_totals(nd);
					}
				}
				if (surface_ptr->Get_type() == cxxSurface::CD_MUSIC)
				{
					charge_ptr->Set_charge_balance(charge_ptr->Get_charge_balance() +
						comp_ptr->Get_moles() *
						comp_ptr->Get_formula_z());
				}
				break;
			}
#ifdef SKIP_MUSIC
			/*
			 *  If charge of formula is non zero
			 */
			if (surface_ptr->type == CD_MUSIC)
			{
				surface_ptr->comps[i].cb =
					surface_ptr->comps[i].formula_z *
					surface_ptr->comps[i].moles;
			}
#endif
		}
		/*
		 * Check that all surface comps have a corresponding master
		 */
		for (size_t i = 0; i < surface_ptr->Get_surface_comps().size(); i++)
		{
			if (surface_ptr->Get_surface_comps()[i].Get_master_element().size() == 0)
			{
				input_error++;
				error_string = sformatf(
						"No surface master species for surface component %s, ",
						surface_ptr->Get_surface_comps()[i].Get_formula().c_str());
				error_msg(error_string, CONTINUE);
			}
		}
/*
 *   Sort components
 */
		std::map<std::string, cxxSurfaceComp> comp_map;
		for (size_t i = 0; i < surface_ptr->Get_surface_comps().size(); i++)
		{
			cxxSurfaceComp *comp_ptr = &(surface_ptr->Get_surface_comps()[i]);
			comp_map[comp_ptr->Get_formula()] = *comp_ptr;
		}
		std::map<std::string, cxxSurfaceComp>::iterator it = comp_map.begin();
		surface_ptr->Get_surface_comps().clear();
		for ( ; it != comp_map.end(); it++)
		{
			surface_ptr->Get_surface_comps().push_back(it->second);
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_solutions(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Define n_user for any solutions read by solution_spread that
 *   don`t have n_user defined
 */
	struct master *master_ptr;

	/*
	 *  Calculate solution numbers
	 */
	if (unnumbered_solutions.size() > 0)
	{
		int last = 0;			
		std::map<int, cxxSolution>::iterator jit;
		for (jit = Rxn_solution_map.begin(); jit != Rxn_solution_map.end(); jit++)
		{
			if (jit->second.Get_n_user() > last)
				last = jit->second.Get_n_user();
			if (jit->second.Get_n_user_end() > last)
				last = jit->second.Get_n_user_end();
		}
		if (save.solution == TRUE)
		{
			if (save.n_solution_user > last)
				last = save.n_solution_user;
			if (save.n_solution_user_end > last)
				last = save.n_solution_user_end;
		}
		// put unnumbered solutions in map
		for (size_t i = 0; i < unnumbered_solutions.size(); i++)
		{
			if (use.Get_n_solution_user() < 0)
			{
				use.Set_n_solution_user(last + 1);
			}
			unnumbered_solutions[i].Set_n_user_both(++last);
			Rxn_solution_map[last] = unnumbered_solutions[i];
			Rxn_new_solution.insert(last);
		}
		unnumbered_solutions.clear();
	}
	/*
	 * Check that elements are in database
	 */
	//std::map<int, cxxSolution>::iterator it;
	//for (it = Rxn_solution_map.begin(); it != Rxn_solution_map.end(); it++)
	//for(size_t n = 0; n < Rxn_new_solution.size(); n++)
	//{
		//std::map<int, cxxSolution>::iterator it = Rxn_solution_map.find(Rxn_new_solution[n]);
	for (std::set<int>::const_iterator nit = Rxn_new_solution.begin(); nit != Rxn_new_solution.end(); nit++)
	{
		std::map<int, cxxSolution>::iterator it = Rxn_solution_map.find(*nit);
		if (it == Rxn_solution_map.end())
		{
			assert(false);
			continue;
		}
		cxxSolution &solution_ref = it->second;
		//if (solution_ref.Get_new_def())
		{
			cxxISolution *initial_data_ptr = solution_ref.Get_initial_data();
			if (initial_data_ptr != NULL)
			{
				std::map<std::string, cxxISolutionComp>::iterator iit = initial_data_ptr->Get_comps().begin();
				for ( ; iit != initial_data_ptr->Get_comps().end(); iit++)
				{
					cxxISolutionComp &comp_ref = iit->second;
					if (strcmp(comp_ref.Get_description().c_str(), "H(1)") == 0 ||
						strcmp(comp_ref.Get_description().c_str(), "E") == 0)
					{
						comp_ref.Set_moles(0.0);
						continue;
					}
					std::string token;
					std::string description = comp_ref.Get_description();
					std::string::iterator b = description.begin();
					std::string::iterator e = description.end();
					CParser::copy_token(token, b, e);

					master_ptr = master_bsearch(token.c_str());
					if (master_ptr == NULL)
					{
						error_string = sformatf(
							"Could not find element in database, %s.\n\tConcentration is set to zero.",
							comp_ref.Get_description().c_str());
						warning_msg(error_string);
						comp_ref.Set_input_conc(0.0);
						continue;
					}
				}
			}
		}
	}
	
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
species_rxn_to_trxn(struct species *s_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copy reaction from reaction structure to 
 *   temp reaction structure.
 */
	int i;

	for (i = 0; s_ptr->rxn->token[i].s != NULL; i++)
	{
		trxn.token[i].name = s_ptr->rxn->token[i].s->name;
		trxn.token[i].z = s_ptr->rxn->token[i].s->z;
		trxn.token[i].s = s_ptr->rxn->token[i].s;
		trxn.token[i].unknown = NULL;
		trxn.token[i].coef = s_ptr->rxn->token[i].coef;
		count_trxn = i + 1;
		if (count_trxn + 1 >= max_trxn)
		{
			space((void **) ((void *) &(trxn.token)), count_trxn + 1,
				  &max_trxn, sizeof(struct rxn_token_temp));
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
phase_rxn_to_trxn(struct phase *phase_ptr, struct reaction *rxn_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copy reaction from reaction structure to 
 *   temp reaction structure.
 */
	int i, l;
	char *ptr;
	char token[MAX_LENGTH];
	LDBLE l_z;

	trxn.token[0].name = phase_ptr->formula;
	/* charge */
	char * temp_formula = string_duplicate(phase_ptr->formula);
	ptr = temp_formula;
	get_token(&ptr, token, &l_z, &l);
	free_check_null(temp_formula);
	trxn.token[0].z = l_z;
	trxn.token[0].s = NULL;
	trxn.token[0].unknown = NULL;
	/*trxn.token[0].coef = -1.0; */
	/* check for leading coefficient of 1.0 for phase did not work */
	trxn.token[0].coef = phase_ptr->rxn->token[0].coef;
	for (i = 1; rxn_ptr->token[i].s != NULL; i++)
	{
		trxn.token[i].name = rxn_ptr->token[i].s->name;
		trxn.token[i].z = rxn_ptr->token[i].s->z;
		trxn.token[i].s = NULL;
		trxn.token[i].unknown = NULL;
		trxn.token[i].coef = rxn_ptr->token[i].coef;
		count_trxn = i + 1;
		if (count_trxn + 1 >= max_trxn)
		{
			space((void **) ((void *) &(trxn.token)), count_trxn + 1,
				  &max_trxn, sizeof(struct rxn_token_temp));
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_isotopes(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Isotope ratios for each element or element valence state
 */
	LDBLE isotope_number;
	struct master *master_ptr, *primary_ptr;

	int primary_number = 0;
	primary_ptr = NULL;
	std::map<int, cxxSolution>::iterator it;
	for (it = Rxn_solution_map.begin(); it != Rxn_solution_map.end(); it++)
	{
		std::map<std::string, cxxSolutionIsotope> new_isotopes;
		cxxSolution &solution_ref = it->second;
		if (!solution_ref.Get_new_def())
			continue;
		if (solution_ref.Get_isotopes().size() == 0)
			continue;
	
		std::map<std::string, cxxSolutionIsotope> primary_isotopes;
/*
 *   Make list of primary master species for isotopes
 */
		std::map < std::string, cxxSolutionIsotope >::iterator kit = solution_ref.Get_isotopes().begin();
		for ( ; kit != solution_ref.Get_isotopes().end(); kit++)
		{
			cxxSolutionIsotope &isotope_ref = kit->second;
			master_ptr = master_bsearch_primary(isotope_ref.Get_elt_name().c_str());
			isotope_number = isotope_ref.Get_isotope_number();
			if (master_ptr == NULL)
			{
				input_error++;
				error_string = sformatf(
						"In isotope calculation: element not defined: %s.",
						isotope_ref.Get_elt_name().c_str());
				error_msg(error_string, CONTINUE);
				continue;
			}
			std::ostringstream iso_name_str;
			iso_name_str << (int) isotope_number << master_ptr->elt->name;

			std::map < std::string, cxxSolutionIsotope >::iterator jit;
			jit = primary_isotopes.find(iso_name_str.str().c_str());
			if (jit == primary_isotopes.end())
			{
				cxxSolutionIsotope temp_isotope;
				temp_isotope.Set_isotope_name(iso_name_str.str().c_str());
				temp_isotope.Set_elt_name(master_ptr->elt->name);
				temp_isotope.Set_isotope_number(isotope_number);
				primary_isotopes[iso_name_str.str().c_str()] = temp_isotope;
			}
		}
		if (get_input_errors() > 0)
			return (ERROR);
/*
 *   Go through all redox states of the list of primary species and isotope number
 */
		for (kit = primary_isotopes.begin(); kit != primary_isotopes.end(); kit++)
		{
			/* find index number of master species, set flag to FALSE */
			master_ptr = master_bsearch(kit->second.Get_elt_name().c_str());
			isotope_number = kit->second.Get_isotope_number();
			for (int k = 0; k < count_master; k++)
			{
				master[k]->isotope = FALSE;
			}
			primary_number = master_ptr->number;
			primary_ptr = master_ptr;

			/* go through isotopes of solution and fill in master species */
			std::map < std::string, cxxSolutionIsotope >::iterator lit = solution_ref.Get_isotopes().begin();
			for ( ; lit != solution_ref.Get_isotopes().end(); lit++)
			{
				master_ptr = master_bsearch(lit->second.Get_elt_name().c_str());
				if (master_ptr == NULL)
				{
					input_error++;
					error_string = sformatf(
							"In isotope calculation: element not defined: %s.",
							lit->second.Get_elt_name().c_str());
					error_msg(error_string, CONTINUE);
					continue;
				}

				/* only fill for pertinent isotope */
				if (master_ptr->elt->primary != primary_ptr)
					continue;
				if (lit->second.Get_isotope_number() !=	isotope_number)
					continue;

				/* for primary, fill in ratio for all secondary species */
				if (master_ptr->primary == TRUE	&& master_ptr->s->secondary != NULL)
				{
					for (int k = primary_number + 1; k < count_master; k++)
					{
						if (master[k]->elt->primary != primary_ptr)
							break;
						master[k]->isotope_ratio = lit->second.Get_ratio();
						master[k]->isotope_ratio_uncertainty = lit->second.Get_ratio_uncertainty();
						if (master[k]->isotope == TRUE)
						{
							error_string = sformatf(
									"In isotope calculation: redefinition of isotope ratio for %s.",
									lit->second.Get_elt_name().c_str());
							error_msg(error_string, CONTINUE);
						}
						master[k]->isotope = TRUE;
					}
				}
				/* for secondary and non redox, set ratio */
				else
				{
					master_ptr->isotope_ratio =	lit->second.Get_ratio();
					master_ptr->isotope_ratio_uncertainty =	lit->second.Get_ratio_uncertainty();
					if (master_ptr->isotope == TRUE)
					{
						error_string = sformatf(
								"In isotope calculation: redefinition of isotope ratio for %s.",
								lit->second.Get_elt_name().c_str());
						error_msg(error_string, CONTINUE);
					}
					master_ptr->isotope = TRUE;
				}
			}
/*
 *   Write new isotope structure
 */
			for (int k = 0; k < count_master; k++)
			{
				/* skip primary master species of redox elements */
				if (master[k]->primary == TRUE && master[k]->s->secondary != NULL)
					continue;
				if (master[k]->elt->primary == primary_ptr && master[k]->isotope == FALSE)
				{
					input_error++;
					error_string = sformatf(
							"Isotopic ratio not defined for element or valence state %g%s, using 0.",
							(double) isotope_number, master[k]->elt->name);
					warning_msg(error_string);
					master[k]->isotope = TRUE;
					master[k]->isotope_ratio = 0.0;
					master[k]->isotope_ratio_uncertainty = 0.001;
				}
				if (master[k]->isotope == FALSE)
					continue;
				cxxSolutionIsotope temp_iso;
				temp_iso.Set_isotope_number(isotope_number);
				temp_iso.Set_elt_name(master[k]->elt->name);
				temp_iso.Set_total(0);
				temp_iso.Set_ratio(master[k]->isotope_ratio);
				temp_iso.Set_ratio_uncertainty(master[k]->isotope_ratio_uncertainty);
				if (master[k]->isotope_ratio_uncertainty != NAN)
				{
					temp_iso.Set_ratio_uncertainty_defined(true);
				}
				std::string token = sformatf("%d%s", (int) isotope_number,
						master[k]->elt->name);
				temp_iso.Set_isotope_name(token.c_str());
				new_isotopes[token] = temp_iso;
			}
		}
		solution_ref.Set_isotopes(new_isotopes);
	}

	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_kin_exchange(void)
/* ---------------------------------------------------------------------- */
/*
 *  If exchanger is related to mineral, exchanger amount is 
 *  set in proportion
 */
{
	cxxKinetics *kinetics_ptr;
	char *ptr;
	LDBLE conc;

	//std::map<int, cxxExchange>::iterator it = Rxn_exchange_map.begin();
	//for ( ; it != Rxn_exchange_map.end(); it++)
	//for (size_t nn = 0; nn < Rxn_new_exchange.size(); nn++)
	//{
		//std::map<int, cxxExchange>::iterator it = Rxn_exchange_map.find(Rxn_new_exchange[nn]);
	for (std::set<int>::const_iterator nit = Rxn_new_exchange.begin(); nit != Rxn_new_exchange.end(); nit++)
	{
		std::map<int, cxxExchange>::iterator it = Rxn_exchange_map.find(*nit);
		if (it == Rxn_exchange_map.end())
		{
			assert(false);
		}
		cxxExchange * exchange_ptr = &(it->second);
		//if (!exchange_ptr->Get_new_def())
		//	continue;
		//if (exchange_ptr->Get_n_user() < 0)
		//	continue;
		// check elements
		for (size_t j = 0; j < exchange_ptr->Get_exchange_comps().size(); j++)
		{
			cxxExchComp & comp_ref = exchange_ptr->Get_exchange_comps()[j];
			if (comp_ref.Get_rate_name().size() == 0)
				continue;
			/* First find exchange master species */
			cxxNameDouble nd = comp_ref.Get_totals();
			cxxNameDouble::iterator kit = nd.begin();
			bool found_exchange = false;
			for (; kit != nd.end(); kit++)
			{
				/* Find master species */
				struct element *elt_ptr = element_store(kit->first.c_str());
				if (elt_ptr == NULL || elt_ptr->master == NULL)
				{
					input_error++;
					error_string = sformatf( "Master species not in database "
							"for %s, skipping element.",
							kit->first.c_str());
					error_msg(error_string, CONTINUE);
					continue;
				}
				if (elt_ptr->master->type == EX)
					found_exchange = true;;
			}
			if (!found_exchange)
			{
				input_error++;
				error_string = sformatf(
						"Exchange formula does not contain an exchange master species, %s",
						comp_ref.Get_formula().c_str());
				error_msg(error_string, CONTINUE);
				continue;
			}

			/* Now find associated kinetic reaction ...  */
			if ((kinetics_ptr = Utilities::Rxn_find(Rxn_kinetics_map, exchange_ptr->Get_n_user())) == NULL)
			{
				input_error++;
				error_string = sformatf(
						"Kinetics %d must be defined to use exchange related to kinetic reaction, %s",
						exchange_ptr->Get_n_user(), comp_ref.Get_formula().c_str());
				error_msg(error_string, CONTINUE);
				continue;
			}
			size_t k;
			for (k = 0; k < kinetics_ptr->Get_kinetics_comps().size(); k++)
			{
				if (strcmp_nocase
					(comp_ref.Get_rate_name().c_str(),
					 kinetics_ptr->Get_kinetics_comps()[k].Get_rate_name().c_str()) == 0)
				{
					break;
				}
			}
			if (k == kinetics_ptr->Get_kinetics_comps().size())
			{
				input_error++;
				error_string = sformatf(
						"Kinetic reaction, %s, related to exchanger, %s, not found in KINETICS %d",
						comp_ref.Get_rate_name().c_str(), comp_ref.Get_formula().c_str(), exchange_ptr->Get_n_user());
				error_msg(error_string, CONTINUE);
				continue;
			}

			/* use database name for phase */
			comp_ref.Set_rate_name(kinetics_ptr->Get_kinetics_comps()[k].Get_rate_name().c_str());

			/* make exchanger concentration proportional to mineral ... */
			conc = kinetics_ptr->Get_kinetics_comps()[k].Get_m() * comp_ref.Get_phase_proportion();

			count_elts = 0;
			paren_count = 0;
			{
				char * temp_formula = string_duplicate(comp_ref.Get_formula().c_str());
				ptr = temp_formula;
				get_elts_in_species(&ptr, conc);
				free_check_null(temp_formula);
			}
			comp_ref.Set_totals(elt_list_NameDouble());
/*
 *   No check on availability of exchange elements 
 */
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_min_exchange(void)
/* ---------------------------------------------------------------------- */
/*
 *  If exchanger is related to mineral, exchanger amount is 
 *  set in proportion
 */
{
	int n, jj;
	char *ptr;
	LDBLE conc;

	//std::map<int, cxxExchange>::iterator it = Rxn_exchange_map.begin();
	//for ( ; it != Rxn_exchange_map.end(); it++)
	//{
	//for (size_t nn = 0; nn < Rxn_new_exchange.size(); nn++)
	//{
		//std::map<int, cxxExchange>::iterator it = Rxn_exchange_map.find(Rxn_new_exchange[nn]);
	for (std::set<int>::const_iterator nit = Rxn_new_exchange.begin(); nit != Rxn_new_exchange.end(); nit++)
	{
		std::map<int, cxxExchange>::iterator it = Rxn_exchange_map.find(*nit);
		if (it == Rxn_exchange_map.end())
		{
			assert(false);
		}
		cxxExchange * exchange_ptr = &(it->second);
		//if (!exchange_ptr->Get_new_def())
		//		continue;
		//if (exchange_ptr->Get_n_user() < 0)
		//	continue;
		n = exchange_ptr->Get_n_user();
		// check elements
		for (size_t j = 0; j < exchange_ptr->Get_exchange_comps().size(); j++)
		{
			cxxExchComp & comp_ref = exchange_ptr->Get_exchange_comps()[j];
			if (comp_ref.Get_phase_name().size() == 0)
				continue;
			/* First find exchange master species */
			cxxNameDouble nd = comp_ref.Get_totals();
			cxxNameDouble::iterator kit = nd.begin();
			bool found_exchange = false;
			for (; kit != nd.end(); kit++)
			{
				/* Find master species */
				struct element *elt_ptr = element_store(kit->first.c_str());
				if (elt_ptr == NULL || elt_ptr->master == NULL)
				{
					input_error++;
					error_string = sformatf( "Master species not in database "
							"for %s, skipping element.",
							kit->first.c_str());
					error_msg(error_string, CONTINUE);
					continue;
				}
				if (elt_ptr->master->type == EX)
				{
					found_exchange = true;;
				}
			}
			if (!found_exchange)
			{
				input_error++;
				error_string = sformatf(
						"Exchange formula does not contain an exchange master species, %s",
						comp_ref.Get_formula().c_str());
				error_msg(error_string, CONTINUE);
				continue;
			}
			cxxPPassemblage *pp_assemblage_ptr = Utilities::Rxn_find(Rxn_pp_assemblage_map, n);
			/* Now find the mineral on which exchanger depends...  */
			if (pp_assemblage_ptr == NULL)
			{
				input_error++;
				error_string = sformatf(
						"Equilibrium_phases %d must be defined to use exchange related to mineral phase, %s",
						n, comp_ref.Get_formula().c_str());
				error_msg(error_string, CONTINUE);
				continue;
			}
			std::map<std::string, cxxPPassemblageComp>::iterator jit;
			jit =  pp_assemblage_ptr->Get_pp_assemblage_comps().begin();
			for ( ; jit != pp_assemblage_ptr->Get_pp_assemblage_comps().end(); jit++)
			{
				if (strcmp_nocase(comp_ref.Get_phase_name().c_str(), jit->first.c_str()) == 0)
				{
					break;
				}
			}
			if (jit == pp_assemblage_ptr->Get_pp_assemblage_comps().end() )
			{
				input_error++;
				error_string = sformatf(
						"Mineral, %s, related to exchanger, %s, not found in Equilibrium_Phases %d",
						comp_ref.Get_phase_name().c_str(), comp_ref.Get_formula().c_str(), n);
				error_msg(error_string, CONTINUE);
				continue;
			}
			/* use database name for phase */
			comp_ref.Set_phase_name(jit->first.c_str());
			/* make exchanger concentration proportional to mineral ... */
			conc = jit->second.Get_moles() * comp_ref.Get_phase_proportion();
			count_elts = 0;
			paren_count = 0;
			{
				char * temp_formula = string_duplicate(comp_ref.Get_formula().c_str());
				ptr = temp_formula;
				get_elts_in_species(&ptr, conc);
				free_check_null(temp_formula);
			}
			comp_ref.Set_totals(elt_list_NameDouble());
/*
 *   make sure exchange elements are in phase
 */
			count_elts = 0;
			paren_count = 0;
			{
				char * temp_formula = string_duplicate(comp_ref.Get_formula().c_str());
				ptr = temp_formula;
				get_elts_in_species(&ptr, -comp_ref.Get_phase_proportion());
				free_check_null(temp_formula);
			}
			int l;
			struct phase *phase_ptr = phase_bsearch(jit->first.c_str(), &l, FALSE);
			{
				char * temp_formula = string_duplicate(phase_ptr->formula);
				ptr = temp_formula;
				get_elts_in_species(&ptr, 1.0);
				free_check_null(temp_formula);
			}
			qsort(elt_list, (size_t) count_elts,
				  (size_t) sizeof(struct elt_list), elt_list_compare);
			elt_list_combine();
			for (jj = 0; jj < count_elts; jj++)
			{
				if (elt_list[jj].elt->primary->s->type != EX
					&& elt_list[jj].coef < 0)
				{
					input_error++;
					error_string = sformatf(
							"Stoichiometry of exchanger, %s * %g mol sites/mol phase,\n\tmust be a subset of the related phase %s, %s.",
							comp_ref.Get_formula().c_str(),
							(double) comp_ref.Get_phase_proportion(),
							phase_ptr->name,
							phase_ptr->formula);
					error_msg(error_string, CONTINUE);
					break;
				}
			}
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_min_surface(void)
/* ---------------------------------------------------------------------- */
/*
 *  If surface is related to mineral, surface amount is 
 *  set in proportion
 */
{
	//std::map<int, cxxSurface>::iterator kit;
	//for (kit = Rxn_surface_map.begin(); kit != Rxn_surface_map.end(); kit++)
	//{
	//for (size_t nn = 0; nn < Rxn_new_surface.size(); nn++)
	//{
	//	std::map<int, cxxSurface>::iterator kit = Rxn_surface_map.find(Rxn_new_surface[nn]);
	for (std::set<int>::const_iterator nit = Rxn_new_surface.begin(); nit != Rxn_new_surface.end(); nit++)
	{
		std::map<int, cxxSurface>::iterator kit = Rxn_surface_map.find(*nit);
		if (kit == Rxn_surface_map.end())
		{
			assert(false);
		}
		cxxSurface *surface_ptr = &(kit->second);
		if (!surface_ptr->Get_new_def()) continue;
		if (!surface_ptr->Get_new_def())
			continue;
		if (surface_ptr->Get_n_user() < 0)
			continue;
		for (size_t j = 0; j < surface_ptr->Get_surface_comps().size(); j++)
		{
			cxxSurfaceComp *surface_comp_ptr = &(surface_ptr->Get_surface_comps()[j]);
			cxxSurfaceCharge *surface_charge_ptr = surface_ptr->Find_charge(surface_comp_ptr->Get_charge_name());
			if (surface_comp_ptr->Get_phase_name().size() == 0)
				continue;
			int n = surface_ptr->Get_n_user();

			/* First find surface master species */

			cxxNameDouble::iterator it;
			for (it = surface_comp_ptr->Get_totals().begin(); it != surface_comp_ptr->Get_totals().end(); it++)
			{
				/* Find master species */
				struct element *elt_ptr = element_store(it->first.c_str());
				struct master *master_ptr = elt_ptr->master;
				if (master_ptr == NULL)
				{
					input_error++;
					error_string = sformatf( "Master species not in database "
							"for %s, skipping element.",
							elt_ptr->name);
					error_msg(error_string, CONTINUE);
					continue;
				}
				if (master_ptr->type != SURF)
					continue;
				surface_comp_ptr->Set_master_element(elt_ptr->name);
				break;
			}
			if (surface_comp_ptr->Get_master_element().size() == 0)
			{
				input_error++;
				error_string = sformatf(
						"Surface formula does not contain a surface master species, %s",
						surface_comp_ptr->Get_formula().c_str());
				error_msg(error_string, CONTINUE);
				continue;
			}

			/* Now find the mineral on which surface depends...  */
			cxxPPassemblage * pp_assemblage_ptr = Utilities::Rxn_find(Rxn_pp_assemblage_map, n);
			if (pp_assemblage_ptr == NULL)
			{
				input_error++;
				error_string = sformatf(
						"Equilibrium_phases %d must be defined to use surface related to mineral phase, %s",
						n, surface_comp_ptr->Get_formula().c_str());
				error_msg(error_string, CONTINUE);
				continue;
			}
			std::map<std::string, cxxPPassemblageComp>::iterator jit;
			jit =  pp_assemblage_ptr->Get_pp_assemblage_comps().begin();
			for ( ; jit != pp_assemblage_ptr->Get_pp_assemblage_comps().end(); jit++)
			{
				if (strcmp_nocase(surface_comp_ptr->Get_phase_name().c_str(),
					 jit->first.c_str()) == 0)
				{
					break;
				}
			}
			if (jit == pp_assemblage_ptr->Get_pp_assemblage_comps().end())
			{
				input_error++;
				error_string = sformatf(
						"Mineral, %s, related to surface, %s, not found in Equilibrium_Phases %d",
						surface_comp_ptr->Get_phase_name().c_str(), surface_comp_ptr->Get_formula().c_str(), n);
				error_msg(error_string, CONTINUE);
				continue;
			}
			int l;
			struct phase *phase_ptr = phase_bsearch(jit->first.c_str(), &l, FALSE);
			if (phase_ptr == NULL)
			{
				input_error++;
				error_string = sformatf(
						"Mineral, %s, related to surface, %s, not found in database.",
						jit->first.c_str(), surface_comp_ptr->Get_formula().c_str());
				error_msg(error_string, CONTINUE);
				continue;
			}
			/* use database name for phase */
			surface_comp_ptr->Set_phase_name(phase_ptr->name);
			/* make surface concentration proportional to mineral ... */
			LDBLE conc = jit->second.Get_moles() * surface_comp_ptr->Get_phase_proportion();
#ifdef SKIP_MUSIC
			comp_ptr->cb = conc * comp_ptr->formula_z;
#endif
/*			if (conc < MIN_RELATED_SURFACE) conc = 0.0; */
			{
				char * temp_formula = string_duplicate(surface_comp_ptr->Get_formula().c_str());
				char *ptr = temp_formula;
				count_elts = 0;
				paren_count = 0;
				get_elts_in_species(&ptr, conc);
				free_check_null(temp_formula);
			}
			{
				if (surface_ptr->Get_new_def())
				{
					cxxNameDouble nd = elt_list_NameDouble();
					surface_comp_ptr->Set_totals(nd);
				}
				else
				{
					surface_comp_ptr->Get_totals()[surface_comp_ptr->Get_master_element()] = conc;
				}
			}

			/* area */
			if (surface_charge_ptr)
			{
				surface_charge_ptr->Set_grams(jit->second.Get_moles());
			}
/*
 *   make sure surface elements are in phase
 *   logically necessary for mass balance and to avoid negative concentrations when dissolving phase
 */
			count_elts = 0;
			paren_count = 0;
			{
				char * temp_formula = string_duplicate(phase_ptr->formula);
				char * ptr = temp_formula;
				get_elts_in_species(&ptr, 1.0);
				free_check_null(temp_formula);
			}
			// Revise logic for surface related to mineral
			for (size_t jj = 0; jj < surface_ptr->Get_surface_comps().size(); jj++)
			{
				cxxSurfaceComp *comp_jj_ptr = &(surface_ptr->Get_surface_comps()[jj]);
				// Use formula for all types of surfaces
				{
					char * temp_formula = string_duplicate(comp_jj_ptr->Get_formula().c_str());
					char *ptr = temp_formula;
					get_elts_in_species(&ptr,
										-comp_jj_ptr->Get_phase_proportion());

					if (surface_ptr->Get_type() != cxxSurface::CD_MUSIC)
					{

						// Warn if not master species and charge balanced
						struct element *elt_ptr = element_store(comp_jj_ptr->Get_master_element().c_str());
						if (elt_ptr->master == NULL)
						{
							input_error++;
							error_string = sformatf("Unknown element definition in SURFACE \n\t for surface related to equilibrium_phase: SURFACE %d.", 
								surface_ptr->Get_n_user());
							error_msg(error_string);
							continue;
						}
						if (elt_ptr->master->s == NULL || elt_ptr->master->s->name == NULL)
						{
							input_error++;
							error_string = sformatf("Unknown master species definition in SURFACE \n\t for surface related to equilibrium_phase: SURFACE %d.", 
								surface_ptr->Get_n_user());
							error_msg(error_string);
							continue;
						}
						if (strcmp(elt_ptr->master->s->name, temp_formula) != 0)
						{
							error_string = sformatf("Suggest using master species formula in SURFACE \n\t for surface related to equilibrium_phase: %s.", 
								elt_ptr->master->s->name);
							warning_msg(error_string);
						}
						if (elt_ptr->master->s->z != 0.0)
						{
							error_string = sformatf(
								"Suggest master species of surface, %s, be uncharged for surface related to equilibrium_phase.",
								elt_ptr->master->s->name);
							warning_msg(error_string);
						}	
					}
					free_check_null(temp_formula);
				}
			}
			qsort(elt_list, (size_t) count_elts,
				  (size_t) sizeof(struct elt_list), elt_list_compare);
			elt_list_combine();
			/* Makes no sense: sorbed species need not be in mineral structure... */
			/* But elements that can desorb into solution must be in mineral */
			/* If you precipitate Ca-Mont, and make SurfMg (assuming this is the 
			   formula in SURFACE), where does the Mg come from? 
			   Further, if you precipitate Ca-Mont, make SurfCa, desorb
			   all the Ca, then dissolve the "Ca-Mont", you must remove SurfCa, or you
			   will end up with Ca in solution. H and O are excluded */
			for (int jj = 0; jj < count_elts; jj++)
			{
				if (elt_list[jj].elt->primary->s->type != SURF
					&& elt_list[jj].coef < 0
					//&& elt_list[jj].elt->primary->s != s_hplus
					//&& elt_list[jj].elt->primary->s != s_h2o
					)
				{
					struct element *elt_ptr = element_store(surface_comp_ptr->Get_master_element().c_str());
					error_string = sformatf(
							"Element %s in sum of surface sites,\n"
							"\t including %s * %g mol sites/mol phase,\n"
							"\t exceeds stoichiometry in the related phase %s, %s.",
							elt_list[jj].elt->name,
							elt_ptr->master->s->name,
							(double) surface_comp_ptr->Get_phase_proportion(),
							phase_ptr->name,
							phase_ptr->formula);
					warning_msg(error_string);
					warning_msg("The mismatch in stoichiometry may cause mass-balance errors or unwanted redox reactions.");
					break;
				}
			}
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_kin_surface(void)
/* ---------------------------------------------------------------------- */
/*
 *  If surface is related to mineral, surface amount is 
 *  set in proportion
 */
{
	cxxKinetics *kinetics_ptr;
	struct phase *phase_ptr;
	struct elt_list *elt_list_kinetics;
	int count_elts_kinetics;

	//std::map<int, cxxSurface>::iterator it;
	//for (it = Rxn_surface_map.begin(); it != Rxn_surface_map.end(); it++)
	//{
	//for (size_t nn = 0; nn < Rxn_new_surface.size(); nn++)
	//{
	//	std::map<int, cxxSurface>::iterator it = Rxn_surface_map.find(Rxn_new_surface[nn]);
	for (std::set<int>::const_iterator nit = Rxn_new_surface.begin(); nit != Rxn_new_surface.end(); nit++)
	{
		std::map<int, cxxSurface>::iterator it = Rxn_surface_map.find(*nit);
		if (it == Rxn_surface_map.end())
		{
			assert(false);
		}
		cxxSurface *surface_ptr = &(it->second);
		if (!surface_ptr->Get_new_def())
			continue;
		if (surface_ptr->Get_n_user() < 0)
			continue;
		int n = surface_ptr->Get_n_user();
		for (size_t j = 0; j < surface_ptr->Get_surface_comps().size(); j++)
		{
			cxxSurfaceComp *comp_ptr = &(surface_ptr->Get_surface_comps()[j]);
			if (comp_ptr->Get_rate_name().size() == 0)
				continue;
			comp_ptr->Set_master_element("");

			/* First find surface master species */
			int k;
			cxxNameDouble::iterator kit;
			for (kit = comp_ptr->Get_totals().begin(); kit != comp_ptr->Get_totals().end(); kit++)
			{
				/* Find master species */
				struct element *elt_ptr = element_store(kit->first.c_str());
				struct master *master_ptr = elt_ptr->master;
				if (master_ptr == NULL)
				{
					input_error++;
					error_string = sformatf( "Master species not in database "
							"for %s, skipping element.",
							elt_ptr->name);
					error_msg(error_string, CONTINUE);
					continue;
				}
				if (master_ptr->type != SURF)
					continue;
				comp_ptr->Set_master_element(elt_ptr->name);
				break;
			}
			if (comp_ptr->Get_master_element().size() == 0)
			{
				input_error++;
				error_string = sformatf(
						"Surface formula does not contain a surface master species, %s",
						comp_ptr->Get_formula().c_str());
				error_msg(error_string, CONTINUE);
				continue;
			}

			/* Now find the kinetic reaction on which surface depends...  */
			if ((kinetics_ptr = Utilities::Rxn_find(Rxn_kinetics_map, n)) == NULL)
			{
				input_error++;
				error_string = sformatf(
						"Kinetics %d must be defined to use surface related to kinetic reaction, %s",
						n, comp_ptr->Get_formula().c_str());
				error_msg(error_string, CONTINUE);
				continue;
			}
			for (k = 0; k < (int) kinetics_ptr->Get_kinetics_comps().size(); k++)
			{
				cxxKineticsComp *kin_comp_ptr = &(kinetics_ptr->Get_kinetics_comps()[k]);
				if (strcmp_nocase
					(comp_ptr->Get_rate_name().c_str(),
					 kin_comp_ptr->Get_rate_name().c_str()) == 0)
				{
					break;
				}
			}
			if (k == (int) kinetics_ptr->Get_kinetics_comps().size())
			{
				input_error++;
				error_string = sformatf(
						"Kinetic reaction, %s, related to surface, %s, not found in Kinetics %d",
						comp_ptr->Get_rate_name().c_str(), comp_ptr->Get_formula().c_str(), n);
				error_msg(error_string, CONTINUE);
				continue;
			}
			cxxKineticsComp *kin_comp_ptr = &(kinetics_ptr->Get_kinetics_comps()[k]);
			/* use database name for phase */
			comp_ptr->Set_rate_name(kin_comp_ptr->Get_rate_name().c_str());

			/* make surface concentration proportional to mineral ... */
			LDBLE conc = kin_comp_ptr->Get_m() * comp_ptr->Get_phase_proportion();

/*			if (conc < MIN_RELATED_SURFACE) conc = 0.0; */
			{
				char * temp_formula = string_duplicate(comp_ptr->Get_formula().c_str());
				char *ptr = temp_formula;
				count_elts = 0;
				paren_count = 0;
				get_elts_in_species(&ptr, conc);
				free_check_null(temp_formula);
			}
			{
				if (surface_ptr->Get_new_def())
				{
					cxxNameDouble nd = elt_list_NameDouble();
					comp_ptr->Set_totals(nd);
				}
				else
				{
					comp_ptr->Get_totals()[comp_ptr->Get_master_element()] = conc;
				}
			}

			/* area */

			cxxSurfaceCharge *charge_ptr = surface_ptr->Find_charge(comp_ptr->Get_charge_name());
			charge_ptr->Set_grams(kin_comp_ptr->Get_m());
		}
/*
 *   check on elements
 */
		/* Go through each kinetic reaction, add all related surface compositions
		 * check for negative values
		 */
		if (!surface_ptr->Get_related_rate())
			continue;
		kinetics_ptr = Utilities::Rxn_find(Rxn_kinetics_map, n);
		if (kinetics_ptr == NULL)
		{
				input_error++;
				error_string = sformatf(
						"Error in SURFACE related to KINETICS. ");
				error_msg(error_string, CONTINUE);
				continue;
		}
		for (size_t k = 0; k < kinetics_ptr->Get_kinetics_comps().size(); k++)
		{
			cxxKineticsComp *kin_comp_ptr = &(kinetics_ptr->Get_kinetics_comps()[k]);
			count_elts = 0;
			paren_count = 0;

			/* added in kinetics formula */
			cxxNameDouble::iterator jit = kin_comp_ptr->Get_namecoef().begin();
			for (; jit != kin_comp_ptr->Get_namecoef().end(); jit++)
			{
				std::string name = jit->first;
				LDBLE coef = jit->second;
				phase_ptr = NULL;
				int jj;
				phase_ptr = phase_bsearch(name.c_str(), &jj, FALSE);
				if (phase_ptr != NULL)
				{
					add_elt_list(phase_ptr->next_elt, 1.0);
				}
				else
				{
					char * temp_name = string_duplicate(name.c_str());
					char * ptr = temp_name;
					get_elts_in_species(&ptr, coef);
					free_check_null(temp_name);
				}
			}
			/* save kinetics formula */
			if (count_elts > 0)
			{
				qsort(elt_list, (size_t) count_elts,
					  (size_t) sizeof(struct elt_list), elt_list_compare);
				elt_list_combine();
			}
			elt_list_kinetics = elt_list_save();
			count_elts_kinetics = count_elts;

			/* get surface formulas */
			count_elts = 0;
			paren_count = 0;
			cxxSurfaceComp *comp_ptr_save = NULL;
			for (size_t j = 0; j < surface_ptr->Get_surface_comps().size(); j++)
			{
				cxxSurfaceComp *comp_ptr = &(surface_ptr->Get_surface_comps()[j]);
				if (comp_ptr->Get_rate_name().size() == 0)
					continue;
				comp_ptr_save = comp_ptr;
				if (strcmp_nocase
					(comp_ptr->Get_rate_name().c_str(),
					 kin_comp_ptr->Get_rate_name().c_str()) == 0)
				{
					char * temp_formula = string_duplicate( comp_ptr->Get_formula().c_str());
					char *ptr = temp_formula;
					get_elts_in_species(&ptr, -1 * comp_ptr->Get_phase_proportion());
					free_check_null(temp_formula);
				}
			}
			if (count_elts > 0)
			{
				qsort(elt_list, (size_t) count_elts,
					  (size_t) sizeof(struct elt_list), elt_list_compare);
				elt_list_combine();
			}
			for (int j = 0; j < count_elts; j++)
			{
				if (elt_list[j].elt == NULL)
				{
					input_error++;
					error_string = sformatf(
						"Cannot identify elements in kinetics component %s.",
						comp_ptr_save->Get_formula().c_str());
					error_msg(error_string, CONTINUE);
					continue;
				}
				if (elt_list[j].elt->primary == NULL )
				{
					input_error++;
					error_string = sformatf(
						"Cannot identify primary element in kinetics component %s.",
						comp_ptr_save->Get_formula().c_str());
					error_msg(error_string, CONTINUE);
					continue;
				}
				if (elt_list[j].elt->primary->s == NULL)
				{
					input_error++;
					error_string = sformatf(
						"Cannot identify primary species for an element in kinetics component %s.",
						comp_ptr_save->Get_formula().c_str());
					error_msg(error_string, CONTINUE);
					continue;
				}

				if (elt_list[j].elt->primary->s->type <= H2O)
				{
					int l;
					for (l = 0; l < count_elts_kinetics; l++)
					{
						if (elt_list[j].elt == elt_list_kinetics[l].elt)
						{
							break;
						}
					}
					if (l == count_elts_kinetics)
					{
						input_error++;
						error_string = sformatf(
								"Stoichiometry of surface, %s * %g mol sites/mol reactant,\n\tmust be a subset of the formula defined for the related reactant %s.\n\tElement %s is not present in reactant formula.",
								comp_ptr_save->Get_formula().c_str(),
								(double) comp_ptr_save->Get_phase_proportion(),
								comp_ptr_save->Get_rate_name().c_str(), elt_list[j].elt->name);
						error_msg(error_string, CONTINUE);
					}
					else if (fabs(elt_list[j].coef) >
							 fabs(elt_list_kinetics[l].coef))
					{
						input_error++;
						error_string = sformatf(
								"Stoichiometry of surface, %s * %g mol sites/mol reactant,\n\tmust be a subset of the formula defined for the related reactant %s.\n\tCoefficient of element %s in surface exceeds amount present in reactant formula.",
								comp_ptr_save->Get_formula().c_str(),
								(double) comp_ptr_save->Get_phase_proportion(),
								comp_ptr_save->Get_rate_name().c_str(), elt_list[j].elt->name);
						error_msg(error_string, CONTINUE);
					}
				}
			}
			elt_list_kinetics =
				(struct elt_list *) free_check_null(elt_list_kinetics);
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
ss_prep(LDBLE t, cxxSS *ss_ptr, int print)
/* ---------------------------------------------------------------------- */
{
	int i, j, k, converged, divisions;
	LDBLE r, rt, ag0, ag1, crit_pt;
	LDBLE xc, tc;
	LDBLE l_x, x0, x1, xsm1, xsm2, xb1, xb2;
	LDBLE xc1, xc2;
	LDBLE facb1, faca1, spim1, xblm1, acrae, acrael, xliapt, xliapm;
	LDBLE xaly, xaly1, xaly2;
	LDBLE faca, facb, spialy, facal, facbl;
	LDBLE tol;

	if (pr.ss_assemblage == FALSE)
		print = FALSE;
	tol = 1e-6;
	r = R_KJ_DEG_MOL;
	rt = r * t;
	a0 = ss_ptr->Get_ag0() / rt;
	a1 = ss_ptr->Get_ag1() / rt;
	ss_ptr->Set_a0(a0);
	ss_ptr->Set_a1(a1);
	ag0 = a0 * rt;
	ag1 = a1 * rt;
	cxxSScomp *comp0_ptr = &(ss_ptr->Get_ss_comps()[0]);
	cxxSScomp *comp1_ptr = &(ss_ptr->Get_ss_comps()[1]);
	struct phase *phase0_ptr = phase_bsearch(comp0_ptr->Get_name().c_str(), &k, FALSE);
	struct phase *phase1_ptr = phase_bsearch(comp1_ptr->Get_name().c_str(), &k, FALSE);
	kc = exp(k_calc(phase0_ptr->rxn->logk, t, REF_PRES_PASCAL) * LOG_10);
	kb = exp(k_calc(phase1_ptr->rxn->logk, t, REF_PRES_PASCAL) * LOG_10);
	crit_pt = fabs(a0) + fabs(a1);
/*
 *   Default, no miscibility or spinodal gaps
 */
	ss_ptr->Set_miscibility(false);
	ss_ptr->Set_spinodal(false);
	xsm1 = 0.5;
	xsm2 = 0.5;
	xb1 = 0.5;
	xb2 = 0.5;
	xc1 = 0;
	xc2 = 0;

	if (crit_pt >= tol)
	{
/*
 *   Miscibility gap information
 */
		if (fabs(a1) < tol)
		{
			xc = 0.5;
			tc = ag0 / (2 * r);
		}
		else
		{
			xc = 0.5 + (pow((ag0 * ag0 + 27 * ag1 * ag1), (LDBLE) 0.5) -
						ag0) / (18 * ag1);
			tc = (12 * ag1 * xc - 6 * ag1 + 2 * ag0) * (xc - xc * xc) / r;
		}
		if (print == TRUE)
		{
			error_string = sformatf( "Description of Solid Solution %s",
					ss_ptr->Get_name().c_str());
			dup_print(error_string, TRUE);
		}
		if (print == TRUE)
		{
			output_msg(sformatf(
					   "\t                              Temperature: %g kelvin\n",
					   (double) t));
			output_msg(sformatf(
					   "\t                       A0 (dimensionless): %g\n",
					   (double) a0));
			output_msg(sformatf(
					   "\t                       A1 (dimensionless): %g\n",
					   (double) a1));
			output_msg(sformatf(
					   "\t                              A0 (kJ/mol): %g\n",
					   (double) ag0));
			output_msg(sformatf(
					   "\t                              A1 (kJ/mol): %g\n\n",
					   (double) ag1));
		}
		if (xc < 0 || xc > 1)
		{
			if (print == TRUE)
				output_msg(sformatf(
						   "No miscibility gap above 0 degrees kelvin.\n"));
		}
		else
		{
			if (print == TRUE)
			{
				output_msg(sformatf(
						   "\t    Critical mole-fraction of component 2: %g\n",
						   (double) xc));
				output_msg(sformatf(
						   "\t                     Critical temperature: %g kelvin\n",
						   (double) tc));
				output_msg(sformatf(
						   "\n(The critical temperature calculation assumes that the Guggenheim model\ndefined at %g kelvin is valid at the critical temperature.)\n\n\n",
						   (double) t));
			}
		}
/*
 *   Calculate miscibility and spinodal gaps
 */
		if (tc >= t)
		{

			/* search for sign changes */
			x0 = 0;
			x1 = 1;
			if (scan(f_spinodal, &x0, &x1) == TRUE)
			{

				/* find first spinodal pt */
				xsm1 = halve(f_spinodal, x0, x1, tol);
				xsm1 = halve(f_spinodal, x0, x1, tol);
				ss_ptr->Set_spinodal(true);

				/* find second spinodal pt */
				x0 = x1;
				x1 = 1;
				if (scan(f_spinodal, &x0, &x1) == TRUE)
				{
					xsm2 = halve(f_spinodal, x0, x1, tol);
				}
				else
				{
					error_msg("Failed to find second spinodal point.", STOP);
				}
			}
		}
	}
/*
 *   Now find Miscibility gap
 */
	if (ss_ptr->Get_spinodal())
	{
		if (print == TRUE)
			output_msg(sformatf(
					   "\t Spinodal-gap mole fractions, component 2: %g\t%g\n",
					   (double) xsm1, (double) xsm2));
		converged = FALSE;
		if (converged == FALSE)
		{
			for (i = 1; i < 3; i++)
			{
				divisions = (int) pow(10., i);
				for (j = 0; j < divisions; j++)
				{
					for (k = divisions; k > 0; k--)
					{
						xc1 = (LDBLE) j / divisions + 0.001;
						xc2 = (LDBLE) k / divisions;
						converged = solve_misc(&xc1, &xc2, tol);
						if (converged == TRUE)
							break;
					}
					if (converged == TRUE)
						break;
				}
				if (converged == TRUE)
					break;
			}
		}
		if (converged == FALSE)
		{
			error_msg("Failed to find miscibility gap.", STOP);
		}
		ss_ptr->Set_miscibility(true);
		if (xc1 < xc2)
		{
			xb1 = 1 - xc2;
			xb2 = 1 - xc1;
			xc1 = 1 - xb1;
			xc2 = 1 - xb2;
		}
		else
		{
			xb1 = 1 - xc1;
			xb2 = 1 - xc2;
		}
		facb1 = kb * xb1 * exp(xc1 * xc1 * (a0 + a1 * (4 * xb1 - 1)));
		faca1 = kc * xc1 * exp(xb1 * xb1 * (a0 - a1 * (3 - 4 * xb1)));
		spim1 = log10(faca1 + facb1);
		xblm1 = 1. / (1. + faca1 / facb1);
		acrae = facb1 / faca1;
		acrael = log10(acrae);
		xliapt = log10(facb1);
		xliapm = log10(faca1);

		if (print == TRUE)
		{
			output_msg(sformatf(
					   "\t   Miscibility-gap fractions, component 2: %g\t%g\n",
					   (double) xb1, (double) xb2));
			output_msg(sformatf(
					   "\n\t\t\tEutectic Point Calculations\n\n"));
			output_msg(sformatf(
					   "\t     Aqueous activity ratio (comp2/comp1): %g\n",
					   (double) acrae));
			output_msg(sformatf(
					   "\t Log aqueous activity ratio (comp2/comp1): %g\n",
					   (double) acrael));
			output_msg(sformatf(
					   "\t Aqueous activity fraction of component 2: %g\n",
					   (double) xblm1));
			output_msg(sformatf(
					   "\t                    Log IAP (component 2): %g\n",
					   (double) xliapt));
			output_msg(sformatf(
					   "\t                    Log IAP (component 1): %g\n",
					   (double) xliapm));
			output_msg(sformatf(
					   "\t                               Log Sum Pi: %g\n",
					   (double) spim1));
		}
		ss_ptr->Set_tk(t);
		ss_ptr->Set_xb1(xb1);
		ss_ptr->Set_xb2(xb2);
	}
/*
 *   Alyotropic point calculation
 */
	xaly = -1.0;
	l_x = a0 * a0 + 3 * a1 * a1 + 6 * a1 * log(kb / kc);
	if (l_x > 0)
	{
		if (fabs(l_x - a0 * a0) >= tol)
		{
			xaly1 = (-(a0 - 3 * a1) + pow(l_x, (LDBLE) 0.5)) / (6 * a1);
			xaly2 = (-(a0 - 3 * a1) - pow(l_x, (LDBLE) 0.5)) / (6 * a1);
			if (xaly1 >= 0 && xaly1 <= 1)
			{
				xaly = xaly1;
			}
			if (xaly2 >= 0 && xaly2 <= 1)
			{
				xaly = xaly2;
			}
		}
		else
		{
			xaly = 0.5 + log(kb / kc) / (2 * a0);
		}
		if (xaly > 0 && xaly < 1)
		{
			faca =
				kc * (1 -
					  xaly) * exp(xaly * xaly * (a0 - a1 * (3 - 4 * xaly)));
			facb =
				kb * xaly * exp((1 - xaly) * (1 - xaly) *
								(a0 + a1 * (4 * xaly - 1.0)));
			spialy = log10(faca + facb);
			facal = log10(faca);
			facbl = log10(facb);
			if (xaly > xb1 && xaly < xb2)
			{
				if (print == TRUE)
					output_msg(sformatf(
							   "\nLocal minimum in the solidus curve coresponding to a maximum\nin the minimum stoichiometric saturation curve.\n\n"));
			}
			else
			{
				if (print == TRUE)
					output_msg(sformatf(
							   "\n\t\t\tAlyotropic Point\n\n"));
			}
			if (print == TRUE)
			{
				output_msg(sformatf(
						   "\t       Solid mole fraction of component 2: %g\n",
						   (double) xaly));
				output_msg(sformatf(
						   "\t                    Log IAP (component 2): %g\n",
						   (double) facbl));
				output_msg(sformatf(
						   "\t                    Log IAP (component 1): %g\n",
						   (double) facal));
				output_msg(sformatf(
						   "\t                               Log Sum Pi: %g\n",
						   (double) spialy));
			}
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
halve(LDBLE f(LDBLE x, void *), LDBLE x0, LDBLE x1, LDBLE tol)
/* ---------------------------------------------------------------------- */
{
	int i;
	LDBLE l_x, y, y0, dx;

	y0 = f(x0, this);
	dx = (x1 - x0);
/*
 *  Loop for interval halving
 */
	for (i = 0; i < 100; i++)
	{
		dx *= 0.5;
		l_x = x0 + dx;
		y = f(l_x, this);
		if (dx < tol || y == 0)
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
 int Phreeqc::
scan(LDBLE f(LDBLE x, void *), LDBLE * xx0, LDBLE * xx1)
/* ---------------------------------------------------------------------- */
{
	int i, j;
	LDBLE fx0, fx1, divisions;
	LDBLE x0, x1, diff;

	x0 = *xx0;
	x1 = *xx1;
	diff = x1 - x0;
	for (j = 0; j < 3; j++)
	{
		fx0 = f(x0, this);
		divisions = (int) pow((LDBLE) 10, (LDBLE) j);
		for (i = 1; i < divisions; i++)
		{
			x1 = *xx0 + diff * (LDBLE) i / divisions;
			fx1 = f(x1, this);
			if (fx0 * fx1 <= 0)
			{
				*xx0 = x0;
				*xx1 = x1;
				return (TRUE);
			}
			x0 = x1;
			fx0 = fx1;
		}
	}
	return (FALSE);
}


/* ---------------------------------------------------------------------- */
 LDBLE Phreeqc::
f_spinodal(LDBLE x, void * cookie)
/* ---------------------------------------------------------------------- */
{
	LDBLE fx;
	Phreeqc * pThis;
	pThis = (Phreeqc *) cookie;
	fx = -12 * pThis->a1 * x * x * x + (18 * pThis->a1 - 
		2 * pThis->a0) * x * x + (2 * pThis->a0 -
		6 * pThis->a1) * x - 1.0;
	return (fx);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
slnq(int n, LDBLE * a, LDBLE * l_delta, int ncols, int print)
/* ---------------------------------------------------------------------- */
{
	int i, j, k, m;
/* debug
 */
	int row;

	LDBLE b;
/* Debug
*/
	if (print == TRUE)
	{
		output_msg(sformatf( "\nArray in slnq: \n\n"));
		for (i = 0; i < ncols - 1; i++)
		{
			row = i * (n + 1);
			for (j = 0; j < ncols; j++)
			{
				output_msg(sformatf( "%10.2e", (double) a[row + j]));
			}
			output_msg(sformatf( "\n"));
		}
		output_msg(sformatf( "\n"));
	}

	if (n == 0)
		return (OK);
/*   Trivial case */
	if (n == 1)
	{
		if (fabs(a[0]) < ZERO_TOL)
			goto slnq_error;
		l_delta[0] = a[1] / a[0];
		return (OK);
	}

/*   Reduction loop */
	for (i = 0; i < n - 1; i++)
	{
		b = fabs(a[i * ncols + i]);
		m = i;

/*   Find maximum value in column */
		for (j = i + 1; j < n; j++)
		{
			if (fabs(a[j * ncols + i]) > b)
			{
				b = fabs(a[j * ncols + i]);
				m = j;
			}
		}

/*   Check for singularity */
		if (b < ZERO_TOL)
			goto slnq_error;

/*   Exchange rows if necessary */
		if (m != i)
		{
			for (j = i; j <= n; j++)
			{
				b = a[i * ncols + j];
				a[i * ncols + j] = a[m * ncols + j];
				a[m * ncols + j] = b;
			}
		}

/*   Make a[i][i]=1.0 */
		for (j = n; j >= i; j--)
		{
			a[i * ncols + j] /= a[i * ncols + i];
		}

/*   Reduction step */
		for (j = i + 1; j < n; j++)
		{
			if (a[j * ncols + i] == 0.0)
				continue;
			b = -a[j * ncols + i];
			for (k = i + 1; k <= n; k++)
			{
				a[j * ncols + k] += b * a[i * ncols + k];
			}
		}
	}

/*   Calculation of l_delta[n] */
	if (fabs(a[(n - 1) * ncols + n - 1]) > ZERO_TOL)
	{
		l_delta[n - 1] = a[(n - 1) * ncols + n] / a[(n - 1) * ncols + n - 1];
	}
	else
	{
		output_msg(sformatf( "Error: Divide by zero in slnq.\n"));
		l_delta[n] = 0.0;
		goto slnq_error;
	}

/*   Back substitution for other delta values */
	for (i = n - 2; i >= 0; i--)
	{
		l_delta[i] = a[i * ncols + n];
		for (j = i + 1; j < n; j++)
		{
			l_delta[i] -= a[i * ncols + j] * l_delta[j];
		}
	}
	if (print == TRUE)
	{
		output_msg(sformatf( "\nResults from slnq: \n\n"));
		for (i = 0; i < n; i++)
		{
			output_msg(sformatf( "%10.2e", (double) l_delta[i]));
		}
		output_msg(sformatf( "\n"));
	}
	return (OK);

  slnq_error:{
		error_string = sformatf(
				"Error: Singular matrix in subroutine slnq. \n");
		warning_msg(error_string);
	}
	return (ERROR);
}

/* ---------------------------------------------------------------------- */
 int Phreeqc::
solve_misc(LDBLE * xxc1, LDBLE * xxc2, LDBLE tol)
/* ---------------------------------------------------------------------- */
{
	int i, repeat, converged, max_iter;
	LDBLE x1, x2, xb1, xb2;
	LDBLE xc1, xc1_2, xc1_3, xc2, xc2_2, xc2_3;
	LDBLE lc1, lc2, lb1, lb2;
	LDBLE a[6], d[2];
	LDBLE t;

	d[0] = d[1] = 0;
	xc1 = *xxc1;
	xc2 = *xxc2;
	x1 = 0;
	x2 = 0;
	converged = TRUE;
	max_iter = 25;
	for (i = 0; i < max_iter; i++)
	{
		/*
		   output_msg(sformatf( "%d  xc1: %g\txc2 %g\n", i, xc1, xc2));
		 */
		xb1 = 1 - xc1;
		xb2 = 1 - xc2;
		xc1_2 = xc1 * xc1;
		xc1_3 = xc1_2 * xc1;
		xc2_2 = xc2 * xc2;
		xc2_3 = xc2_2 * xc2;

		lc1 = exp(xb1 * xb1 * (a0 - a1 * (3 - 4 * xb1)));
		lb1 = exp(xc1 * xc1 * (a0 + a1 * (4 * xb1 - 1)));
		lc2 = exp(xb2 * xb2 * (a0 - a1 * (3 - 4 * xb2)));
		lb2 = exp(xc2 * xc2 * (a0 + a1 * (4 * xb2 - 1)));

		/* -fb */
		a[2] = -(xb1 * lb1 - xb2 * lb2);

		/* -fc */
		a[5] = -(xc1 * lc1 - xc2 * lc2);

		if (fabs(a[2]) < tol && fabs(a[5]) < tol)
			break;

		/* dfb/dxc1 */
		t = exp(a0 * xc1_2 - 4 * a1 * xc1_3 + 3 * a1 * xc1_2);
		a[0] =
			(2 * a0 * xc1 + 6 * a1 * xc1 - 2 * a0 * xc1_2 + 12 * a1 * xc1_3 -
			 18 * a1 * xc1_2 - 1) * t;

		/* dfb/dxc2 */
		t = exp(a0 * xc2_2 - 4 * a1 * xc2_3 + 3 * a1 * xc2_2);
		a[1] =
			(2 * a0 * xc2_2 - 12 * a1 * xc2_3 - 2 * a0 * xc2 +
			 18 * a1 * xc2_2 - 6 * a1 * xc2 + 1) * t;


		/* dfc/dxc1 */
		t = exp(a0 * xc1_2 - 2 * a0 * xc1 + a0 - 4 * a1 * xc1_3 +
				9 * a1 * xc1_2 - 6 * a1 * xc1 + a1);
		a[3] =
			(2 * a0 * xc1_2 - 2 * a0 * xc1 - 12 * a1 * xc1_3 +
			 18 * a1 * xc1_2 - 6 * a1 * xc1 + 1) * t;

		/* dfc/dxc2 */
		t = exp(a0 * xc2_2 - 2 * a0 * xc2 + a0 - 4 * a1 * xc2_3 +
				9 * a1 * xc2_2 - 6 * a1 * xc2 + a1);
		a[4] =
			(-2 * a0 * xc2_2 + 2 * a0 * xc2 + 12 * a1 * xc2_3 -
			 18 * a1 * xc2_2 + 6 * a1 * xc2 - 1) * t;


		/* solve for dxc1 and dxc2 */
		slnq(2, a, d, 3, FALSE);

		repeat = TRUE;
		while (repeat == TRUE)
		{
			x1 = xc1 + d[0];
			x2 = xc2 + d[1];
			if (x1 > 1 || x1 < 0 || x2 > 1 || x2 < 0)
			{
				d[0] *= 0.5;
				d[1] *= 0.5;
			}
			else
			{
				repeat = FALSE;
			}
		};
		xc1 = x1;
		xc2 = x2;

		if (fabs(xc1 - xc2) < .01)
		{
			converged = FALSE;
			break;
		}
	}
	if (i == max_iter)
		converged = FALSE;
	*xxc1 = xc1;
	*xxc2 = xc2;
	return (converged);
}
/* ---------------------------------------------------------------------- */
 int Phreeqc::
ss_calc_a0_a1(cxxSS *ss_ptr)
/* ---------------------------------------------------------------------- */
{
	int i, done;
	LDBLE r, rt;
	std::vector<LDBLE> p;
	LDBLE q1, q2, xbq1, xbq2, xb1, xb2, xc1, xc2;
	LDBLE r1, r2, pa1, pb1, pa2, pb2, xsm1, xsm2;
	LDBLE pn9, pn10, c5, c6, pl9, pl10, pj9, pj10;
	LDBLE xc, tc;
	LDBLE spialy, azero, phi1, phi2, test;
	LDBLE dq1, dq2, denom, ratio, dr1, dr2, x21, x22, x61, x62;
	LDBLE l_a0, l_a1, ag0, ag1;
	LDBLE wg2, wg1, alpha2, alpha3;
	LDBLE l_kc, l_kb;
	LDBLE xaly, xcaly, alpha0, alpha1, fx, fx1;
	LDBLE tol;

	tol = 1e-6;
	rt = ss_ptr->Get_tk() * R_KJ_DEG_MOL;
	if (ss_ptr->Get_ss_comps().size() < 2)
	{
		input_error++;
		error_string = sformatf(
				"Two components not defined for solid solution ",
				ss_ptr->Get_name().c_str());
		error_msg(error_string, CONTINUE);
		return (ERROR);
	}
	cxxSScomp *comp0_ptr = &(ss_ptr->Get_ss_comps()[0]);
	cxxSScomp *comp1_ptr = &(ss_ptr->Get_ss_comps()[1]);
	int k;
	struct phase *phase0_ptr = phase_bsearch(comp0_ptr->Get_name().c_str(), &k, FALSE);
	struct phase *phase1_ptr = phase_bsearch(comp1_ptr->Get_name().c_str(), &k, FALSE);
	if (phase0_ptr == NULL || phase1_ptr == NULL)
	{
		input_error++;
		error_string = sformatf(
				"Two components were not defined for %s solid solution",
				ss_ptr->Get_name().c_str());
		error_msg(error_string, CONTINUE);
		return (ERROR);
	}
	l_kc = exp(k_calc(phase0_ptr->rxn->logk, ss_ptr->Get_tk(), REF_PRES_PASCAL) *
			 LOG_10);
	l_kb = exp(k_calc(phase1_ptr->rxn->logk, ss_ptr->Get_tk(), REF_PRES_PASCAL) *
			 LOG_10);

	p = ss_ptr->Get_p();

	l_a0 = 0;
	l_a1 = 0;
	ag0 = 0;
	ag1 = 0;
	dq2 = 0;
	switch (ss_ptr->Get_input_case())
	{
		/*
		 *  dimensionless a0 and a1
		 */
	case cxxSS::SS_PARM_A0_A1:
		l_a0 = p[0];
		l_a1 = p[1];
		ag0 = l_a0 * rt;
		ag1 = l_a1 * rt;
		break;
		/*
		 *  two activity coefficients
		 *  q1, q2, xbq1, xbq2
		 */
	case cxxSS::SS_PARM_GAMMAS:
		q1 = p[0];
		q2 = p[1];
		xbq1 = p[2];
		xbq2 = p[3];
		done = FALSE;
		if (fabs(1 - xbq1) > 0 && q1 > 0)
		{
			dq1 = log(q1) / ((1 - xbq1) * (1 - xbq1));
			if (xbq2 <= 0 || xbq2 > 1)
			{
				l_a0 = dq1;
				l_a1 = 0;
				done = TRUE;
			}
		}
		if (done == FALSE)
		{
			if (fabs(xbq2) < 0 || q2 <= 0)
			{
				input_error++;
				error_string = sformatf(
						"No solution possible for A0 and A1 calculation from two activity coefficients, %s.\n",
						ss_ptr->Get_name().c_str());
				error_msg(error_string, CONTINUE);
				done = TRUE;
			}
		}
		if (done == FALSE)
		{
			dq2 = log(q2) / (xbq2 * xbq2);
			if (xbq1 < 0. || xbq2 > 1.)
			{
				l_a0 = dq2;
				l_a1 = 0;
				done = TRUE;
			}
		}
		if (done == FALSE)
		{
			denom = 4 * (xbq1 - xbq2) + 2;
			if (fabs(denom) >= tol)
			{
				if (fabs(1 - xbq1) > 0 && q1 > 0)
				{
					dq1 = log(q1) / ((1 - xbq1) * (1 - xbq1));
					l_a0 = (dq1 * (3 - 4 * xbq2) +
						  dq2 * (4 * xbq1 - 1)) / denom;
					l_a1 = (dq1 - dq2) / denom;
					done = TRUE;
				}
			}
		}
		if (done == FALSE)
		{
			input_error++;
			error_string = sformatf(
					"No solution possible for A0 and A1 calculation from two activity coefficients, %s.\n",
					ss_ptr->Get_name().c_str());
			error_msg(error_string, CONTINUE);
		}
		/* io = 1 */
		ag0 = l_a0 * rt;
		ag1 = l_a1 * rt;
		break;
		/*
		 *  two distribution coefficients
		 *  q1, q2, xbq1, xbq2
		 */
	case cxxSS::SS_PARM_DIST_COEF:
		q1 = p[0];
		q2 = p[1];
		xbq1 = p[2];
		xbq2 = p[3];
		ratio = l_kc / l_kb;
		dr1 = log(q1 / ratio);
		x21 = 2 * xbq1 - 1;
		if (fabs(xbq1 - xbq2) < tol || xbq2 < 0)
		{
			l_a0 = dr1 / x21;
			l_a1 = 0;
		}
		else
		{
			dr2 = log(q2 / ratio);
			x22 = 2 * xbq2 - 1;
			if (xbq1 < 0.)
			{
				l_a0 = dr2 / x22;
				l_a1 = 0;
			}
			else
			{
				x61 = 6 * xbq1 * xbq1 - 6 * xbq1 + 1;
				x62 = 6 * xbq2 * xbq2 - 6 * xbq2 + 1;
				if (fabs(x22 * x61 - x21 * x62) < tol)
				{
					input_error++;
					error_string = sformatf(
							"No solution possible for A0 and A1 calculation from two distribution coefficients, %s.\n",
							ss_ptr->Get_name().c_str());
					error_msg(error_string, CONTINUE);
				}
				l_a0 = (x61 * dr2 - x62 * dr1) / (x22 * x61 - x21 * x62);
				l_a1 = (x21 * dr2 - x22 * dr1) / (x21 * x62 - x22 * x61);
			}
		}

		/* io = 1 */
		ag0 = l_a0 * rt;
		ag1 = l_a1 * rt;
		break;
		/*
		 *  from miscibility gap fractions
		 *  q1, q2
		 */
	case cxxSS::SS_PARM_MISCIBILITY:
		q1 = p[0];
		q2 = p[1];
		xb1 = q1;
		xb2 = q2;
		xc1 = 1 - xb1;
		xc2 = 1 - xb2;
		r1 = log(xb1 / xb2);
		r2 = log(xc1 / xc2);
		pa1 = xc2 * xc2 - xc1 * xc1;
		pb1 =
			3 * (xc2 * xc2 - xc1 * xc1) - 4 * (xc2 * xc2 * xc2 -
											   xc1 * xc1 * xc1);
		pa2 = xb2 * xb2 - xb1 * xb1;
		pb2 =
			-(3 * (xb2 * xb2 - xb1 * xb1) -
			  4 * (xb2 * xb2 * xb2 - xb1 * xb1 * xb1));
		l_a0 = (r1 - pb1 / pb2 * r2) / (pa1 - pa2 * pb1 / pb2);
		l_a1 = (r1 - pa1 / pa2 * r2) / (pb1 - pb2 * pa1 / pa2);

		/* io = 1 */
		ag0 = l_a0 * rt;
		ag1 = l_a1 * rt;
		break;
		/*
		 *  from spinodal gap fractions
		 *  q1, q2
		 */
	case cxxSS::SS_PARM_SPINODAL:
		q1 = p[0];
		q2 = p[1];
		xsm1 = q1;
		xsm2 = q2;
		pn9 = 1 / xsm1;
		pn10 = 1 / xsm2;
		c5 = 1 - xsm1;
		c6 = 1 - xsm2;
		pl9 = 6 * c5 - 12 * c5 * c5;
		pl10 = 6 * c6 - 12 * c6 * c6;
		pj9 = 2 * c5;
		pj10 = 2 * c6;
		l_a0 = (pn9 - pl9 / pl10 * pn10) / (pj9 - pl9 / pl10 * pj10);
		l_a1 = (pn9 - pj9 / pj10 * pn10) / (pl9 - pj9 / pj10 * pl10);

		/* io = 1 */
		ag0 = l_a0 * rt;
		ag1 = l_a1 * rt;
		break;
		/*
		 *  from critical point
		 *  q1, q2
		 */
	case cxxSS::SS_PARM_CRITICAL:
		xc = p[0];
		tc = p[1];
		r = R_KJ_DEG_MOL;
		ag1 = r * tc * (2 * xc - 1) / (12 * xc * xc * (1 - xc) * (1 - xc));
		ag0 = (r * tc / (xc * (1 - xc)) - (12 * xc - 6) * ag1) / 2;

		/* io = 0 */
		l_a0 = ag0 / rt;
		l_a1 = ag1 / rt;
		break;
		/*
		 *  from alyotropic point
		 *  q1, q2
		 */
	case cxxSS::SS_PARM_ALYOTROPIC:
		q1 = p[0];
		q2 = p[1];
		xaly = q1;
		r = log(l_kb / l_kc);
		alpha0 = 2 * xaly - 1;
		alpha1 = 6 * xaly * (xaly - 1) + 1;
		spialy = pow((LDBLE) 10., q2);
		l_a0 = -999.;
		l_a1 = -999.;
		if (fabs(alpha0) < tol)
		{
			input_error++;
			error_string = sformatf(
					"No solution possible for A0 and A1 calculation from alyotropic point, %s.\n",
					ss_ptr->Get_name().c_str());
			error_msg(error_string, CONTINUE);
		}
		else
		{
			azero = 1;
			if (fabs(alpha0) > tol)
				azero = r / alpha0;
			xcaly = 1 - xaly;
/*
 *  Solve for a0 by Newton's method
 */
			for (i = 0; i < 50; i++)
			{
				phi1 =
					xcaly * xcaly * (azero +
									 (r - azero * alpha0) * (4 * xaly -
															 1) / alpha1);
				phi2 =
					xaly * xaly * (azero +
								   (3 - 4 * xaly) * (azero * alpha0 -
													 r) / alpha1);
				phi1 = xaly * l_kb * exp(phi1);
				phi2 = xcaly * l_kc * exp(phi2);
				fx = phi1 + phi2 - spialy;
				fx1 =
					xcaly * xcaly * (1 -
									 alpha0 * (4 * xaly -
											   1) / alpha1) * phi1 +
					xaly * xaly * (1 +
								   alpha0 * (3 - 4 * xaly) / alpha1) * phi2;
				if (fabs(fx1) < 1e-10)
				{
					input_error++;
					error_string = sformatf(
							"Could not find A0 and A1 calculation from alyotropic point, %s.\n",
							ss_ptr->Get_name().c_str());
					error_msg(error_string, CONTINUE);
					break;
				}
				l_a0 = azero - fx / fx1;
				test = fabs(l_a0 - azero) + fabs(fx);
				azero = l_a0;
				if (test < tol)
					break;
			}
			if (i == 50)
			{
				input_error++;
				error_string = sformatf(
						"Too many iterations, could not find A0 and A1 calculation from alyotropic point, %s.\n",
						ss_ptr->Get_name().c_str());
				error_msg(error_string, CONTINUE);
			}
			else
			{
				l_a1 = (r - l_a0 * alpha0) / alpha1;

				/* io = 0 */
				ag0 = l_a0 * rt;
				ag1 = l_a1 * rt;
			}

		}
		break;
		/*
		 *  dimensional (kJ/mol) Guggenheim parameters
		 *  ag0, ag1
		 */
	case cxxSS::SS_PARM_DIM_GUGG:
		ag0 = p[0];
		ag1 = p[1];
		l_a0 = ag0 / rt;
		l_a1 = ag1 / rt;
		break;
		/*
		 *  Waldbaum-Thompson
		 *  wg2, wg1
		 */
	case cxxSS::SS_PARM_WALDBAUM:
		wg2 = p[0];
		wg1 = p[1];
		ag0 = (wg2 + wg1) / 2;
		ag1 = (wg2 - wg1) / 2;
		l_a0 = ag0 / rt;
		l_a1 = ag1 / rt;
		break;
		/*
		 *  Margules
		 *  alpha2, alpha3
		 */
	case cxxSS::SS_PARM_MARGULES:
		alpha2 = p[0];
		alpha3 = p[1];
		l_a0 = alpha2 + 3 * alpha3 / 4;
		l_a1 = alpha3 / 4;
		ag0 = l_a0 * rt;
		ag1 = l_a1 * rt;
		break;
	case cxxSS::SS_PARM_NONE:
		break;
	}
	ss_ptr->Set_ag0(ag0);
	ss_ptr->Set_ag1(ag1);
	ss_ptr->Set_a0(l_a0);
	ss_ptr->Set_a1(l_a1);
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_master_isotope(void)
/* ---------------------------------------------------------------------- */
{
	int i;
	struct master *master_ptr;

	for (i = 0; i < count_master_isotope; i++)
	{
		/*
		 * Mark master species list as minor isotope
		 */
		if (master_isotope[i]->minor_isotope == TRUE)
		{
			master_ptr = master_bsearch(master_isotope[i]->name);
			if (master_ptr == NULL)
			{
				input_error++;
				error_string = sformatf(
						"Did not find master species for isotope, %s",
						master_isotope[i]->name);
				error_msg(error_string, CONTINUE);
				master_isotope[i]->master = NULL;
				continue;
			}
			else
			{
				master_isotope[i]->master = master_ptr;
			}
			master_ptr->minor_isotope = TRUE;
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_isotope_ratios(void)
/* ---------------------------------------------------------------------- */
{
	int i;
	struct master *master_ptr;
	struct master_isotope *master_isotope_ptr;
	struct calculate_value *calculate_value_ptr;

	for (i = 0; i < count_isotope_ratio; i++)
	{
		/*
		 * Mark master species list as minor isotope
		 */
		master_isotope_ptr =
			master_isotope_search(isotope_ratio[i]->isotope_name);
		if (master_isotope_ptr == NULL)
		{
			input_error++;
			error_string = sformatf(
					"For ISOTOPE_RATIO %s, did not find ISOTOPE definition for this isotope, %s",
					isotope_ratio[i]->name, isotope_ratio[i]->isotope_name);
			error_msg(error_string, CONTINUE);
		}
		master_ptr = master_bsearch(isotope_ratio[i]->isotope_name);
		if (master_ptr == NULL)
		{
			input_error++;
			error_string = sformatf(
					"For ISOTOPE_RATIO %s, did not find SOLUTION_MASTER_SPECIES for isotope, %s",
					isotope_ratio[i]->name, isotope_ratio[i]->isotope_name);
			error_msg(error_string, CONTINUE);
		}
		calculate_value_ptr = calculate_value_search(isotope_ratio[i]->name);
		if (calculate_value_ptr == NULL)
		{
			input_error++;
			error_string = sformatf(
					"For ISOTOPE_RATIOS %s, did not find corresponding CALCULATE_VALUE definition",
					isotope_ratio[i]->name);
			error_msg(error_string, CONTINUE);
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_isotope_alphas(void)
/* ---------------------------------------------------------------------- */
{
	int i;
	struct calculate_value *calculate_value_ptr;
	struct logk *logk_ptr;

	for (i = 0; i < count_isotope_alpha; i++)
	{
		/*
		 * Mark master species list as minor isotope
		 */
		calculate_value_ptr = calculate_value_search(isotope_alpha[i]->name);
		if (calculate_value_ptr == NULL)
		{
			input_error++;
			error_string = sformatf(
					"For ISOTOPE_ALPHAS %s, did not find corresponding CALCULATE_VALUE definition",
					isotope_alpha[i]->name);
			error_msg(error_string, CONTINUE);
		}
		if (isotope_alpha[i]->named_logk != NULL)
		{
			logk_ptr = logk_search(isotope_alpha[i]->named_logk);
			if (logk_ptr == NULL)
			{
				input_error++;
				error_string = sformatf(
						"For ISOTOPE_ALPHAS %s, did not find corresponding NAMED_EXPRESSION definition %s.",
						isotope_alpha[i]->name, isotope_alpha[i]->named_logk);
				error_msg(error_string, CONTINUE);
			}
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
reset_last_model(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Initialize model
 */
	last_model.force_prep = TRUE;
	last_model.count_exchange = 0;
	last_model.exchange =
		(struct master **) free_check_null(last_model.exchange);
	last_model.count_gas_phase = 0;
	last_model.gas_phase =
		(struct phase **) free_check_null(last_model.gas_phase);
	last_model.count_ss_assemblage = 0;
	last_model.ss_assemblage =
		(const char **) free_check_null(last_model.ss_assemblage);
	last_model.count_pp_assemblage = 0;
	last_model.pp_assemblage =
		(struct phase **) free_check_null(last_model.pp_assemblage);
	last_model.add_formula =
		(const char **) free_check_null(last_model.add_formula);
	last_model.si = (LDBLE *) free_check_null(last_model.si);
	last_model.dl_type = cxxSurface::NO_DL;
	last_model.count_surface_comp = 0;
	last_model.surface_comp =
		(const char **) free_check_null(last_model.surface_comp);
	last_model.count_surface_charge = 0;
	last_model.surface_charge =
		(const char **) free_check_null(last_model.surface_charge);
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
tidy_exchange(void)
/* ---------------------------------------------------------------------- */
/*
 *  If exchanger is related to mineral, exchanger amount is 
 *  set in proportion
 */
{
	//std::map<int, cxxExchange>::iterator it = Rxn_exchange_map.begin();
	//for ( ; it != Rxn_exchange_map.end(); it++)
	//for (size_t nn = 0; nn < Rxn_new_exchange.size(); nn++)
	//{
	//	std::map<int, cxxExchange>::iterator it = Rxn_exchange_map.find(Rxn_new_exchange[nn]);
	for (std::set<int>::const_iterator nit = Rxn_new_exchange.begin(); nit != Rxn_new_exchange.end(); nit++)
	{
		std::map<int, cxxExchange>::iterator it = Rxn_exchange_map.find(*nit);
		if (it == Rxn_exchange_map.end())
		{
			assert(false);
		}
		//std::map<int, cxxExchange>::iterator it = Rxn_exchange_map.begin();
		cxxExchange * exchange_ptr = &(it->second);
		//if (!exchange_ptr->Get_new_def())
		//	continue;
		//if (exchange_ptr->Get_n_user() < 0)
		//	continue;
		// check elements
		for (size_t j = 0; j < exchange_ptr->Get_exchange_comps().size(); j++)
		{
			cxxExchComp & comp_ref = exchange_ptr->Get_exchange_comps()[j];
			if (comp_ref.Get_phase_name().size() > 0)
				continue;
			if (comp_ref.Get_rate_name().size() > 0)
				continue;

			/* Check elements */
			cxxNameDouble nd = comp_ref.Get_totals();
			cxxNameDouble::iterator kit = nd.begin();
			for (; kit != nd.end(); kit++)
			{
				/* Find master species */
				struct element *elt_ptr = element_store(kit->first.c_str());
				if (elt_ptr == NULL || elt_ptr->master == NULL)
				{
					input_error++;
					error_string = sformatf( "Master species not in database "
							"for %s, skipping element.",
							kit->first.c_str());
					error_msg(error_string, CONTINUE);
					break;
				}
			}
		}
	}
	return (OK);
}

