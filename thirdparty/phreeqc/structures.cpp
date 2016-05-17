#include "Utils.h"
#include "Phreeqc.h"
#include <iostream>

#include "phqalloc.h"
#include "Temperature.h"
#include "cxxMix.h"
#include "Exchange.h"
#include "GasPhase.h"
#include "Reaction.h"
#include "PPassemblage.h"
#include "Use.h"
#include "SSassemblage.h"
#include "cxxKinetics.h"
#include "Surface.h"
#include "Solution.h"


/* ---------------------------------------------------------------------- */
int Phreeqc::
clean_up(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Free all allocated memory, except strings
 */
	int i, j;
#if defined MULTICHART
	chart_handler.End_timer();
	output_flush();
#if 0
	// Wait for charts to end
	while (0 != this->chart_handler.Get_active_charts())
	{
		System::Threading::Thread::Sleep(60);
	}
#endif
#endif

	description_x = (char *) free_check_null(description_x);
	isotopes_x.clear();
	moles_per_kilogram_string =
		(char *) free_check_null(moles_per_kilogram_string);
	pe_string = (char *) free_check_null(pe_string);
/* model */
	last_model.exchange =
		(struct master **) free_check_null(last_model.exchange);
	last_model.gas_phase =
		(struct phase **) free_check_null(last_model.gas_phase);
	last_model.pp_assemblage =
		(struct phase **) free_check_null(last_model.pp_assemblage);
	last_model.ss_assemblage =
		(const char **) free_check_null(last_model.ss_assemblage);
	last_model.add_formula =
		(const char **) free_check_null(last_model.add_formula);
	last_model.si = (LDBLE *) free_check_null(last_model.si);
	last_model.surface_comp =
		(const char **) free_check_null(last_model.surface_comp);
	last_model.surface_charge =
		(const char **) free_check_null(last_model.surface_charge);

	/* model */
	free_model_allocs();

/* species */

	for (j = 0; j < count_s; j++)
	{
		s_free(s[j]);
		s[j] = (struct species *) free_check_null(s[j]);
	}
	s = (struct species **) free_check_null(s);

/* master species */

	for (j = 0; j < count_master; j++)
	{
		master_free(master[j]);
	}
	master = (struct master **) free_check_null(master);

/* elements */

	for (j = 0; j < count_elements; j++)
	{
		elements[j] = (struct element *) free_check_null(elements[j]);
	}
	elements = (struct element **) free_check_null(elements);

/* solutions */
	Rxn_solution_map.clear();

/* surfaces */
	Rxn_surface_map.clear();

/* exchange */
	Rxn_exchange_map.clear();

/* pp assemblages */
	Rxn_pp_assemblage_map.clear();

/* s_s assemblages */
	Rxn_ss_assemblage_map.clear();

/* irreversible reactions */
	Rxn_reaction_map.clear();

/* temperature */
	Rxn_temperature_map.clear();

/* pressure */
	Rxn_pressure_map.clear();

/* unknowns */

	for (j = 0; j < max_unknowns; j++)
	{
		unknown_free(x[j]);
	}
	x = (struct unknown **) free_check_null(x);

/* mixtures */
	Rxn_mix_map.clear();

/* phases */

	for (j = 0; j < count_phases; j++)
	{
		phase_free(phases[j]);
		phases[j] = (struct phase *) free_check_null(phases[j]);
	}
	phases = (struct phase **) free_check_null(phases);

/* inverse */
	for (j = 0; j < count_inverse; j++)
	{
		inverse_free(&(inverse[j]));
	}
	inverse = (struct inverse *) free_check_null(inverse);

/* gases */
	Rxn_gas_phase_map.clear();

/* kinetics */
	Rxn_kinetics_map.clear();
	x0_moles = (LDBLE *) free_check_null(x0_moles);
	m_temp = (LDBLE *) free_check_null(m_temp);
	m_original = (LDBLE *) free_check_null(m_original);
	rk_moles = (LDBLE *) free_check_null(rk_moles);

/* rates */
	for (j = 0; j < count_rates; j++)
	{
		rate_free(&rates[j]);
	}
	rates = (struct rate *) free_check_null(rates);

/* logk hash table */
	for (j = 0; j < count_logk; j++)
	{
		free_check_null(logk[j]->add_logk);
		logk[j] = (struct logk *) free_check_null(logk[j]);
	}
	logk = (struct logk **) free_check_null(logk);

/* save_values */
	for (j = 0; j < count_save_values; j++)
	{
		save_values[j].subscripts =
			(int *) free_check_null(save_values[j].subscripts);
	}
	save_values = (struct save_values *) free_check_null(save_values);

/* model */

/* global solution */

	pe_x.clear();

/* species_list */

	species_list = (struct species_list *) free_check_null(species_list);

/* transport data */

	stag_data = (struct stag_data *) free_check_null(stag_data);
	cell_data = (struct cell_data *) free_check_null(cell_data);

/* punch */
#ifdef SKIP
	punch.totals = (struct name_master *) free_check_null(punch.totals);
	punch.molalities =
		(struct name_species *) free_check_null(punch.molalities);
	punch.activities =
		(struct name_species *) free_check_null(punch.activities);
	punch.pure_phases =
		(struct name_phase *) free_check_null(punch.pure_phases);
	punch.si = (struct name_phase *) free_check_null(punch.si);
	punch.gases = (struct name_phase *) free_check_null(punch.gases);
	punch.s_s = (struct name_phase *) free_check_null(punch.s_s);
	punch.kinetics = (struct name_phase *) free_check_null(punch.kinetics);
#endif
	advection_punch = (int *) free_check_null(advection_punch);
	advection_print = (int *) free_check_null(advection_print);
#ifdef SKIP
	punch.isotopes = (struct name_master *) free_check_null(punch.isotopes);
	punch.calculate_values =
		(struct name_master *) free_check_null(punch.calculate_values);
#endif
	SelectedOutput_map.clear();
	UserPunch_map.clear();

/*  user_print and user_punch */
	rate_free(user_print);
	user_print = (struct rate *) free_check_null(user_print);
#ifdef SKIP
	rate_free(user_punch);
	user_print = (struct rate *) free_check_null(user_print);

	user_punch = (struct rate *) free_check_null(user_punch);
	user_punch_headings = (const char **) free_check_null(user_punch_headings);
#endif

	/*
	   Free llnl aqueous model parameters
	 */
	llnl_temp = (LDBLE *) free_check_null(llnl_temp);
	llnl_adh = (LDBLE *) free_check_null(llnl_adh);
	llnl_bdh = (LDBLE *) free_check_null(llnl_bdh);
	llnl_bdot = (LDBLE *) free_check_null(llnl_bdot);
	llnl_co2_coefs = (LDBLE *) free_check_null(llnl_co2_coefs);
	/*
	 * Copier space
	 */
	copier_free(&copy_solution);
	copier_free(&copy_pp_assemblage);
	copier_free(&copy_exchange);
	copier_free(&copy_surface);
	copier_free(&copy_ss_assemblage);
	copier_free(&copy_gas_phase);
	copier_free(&copy_kinetics);
	copier_free(&copy_mix);
	copier_free(&copy_reaction);
	copier_free(&copy_temperature);
	copier_free(&copy_pressure);

#if defined PHREEQ98 
	rate_free(user_graph);
	user_graph = (struct rate *) free_check_null(user_graph);
	user_graph_headings = (char **) free_check_null(user_graph_headings);
#endif
	/* master_isotope */
	for (i = 0; i < count_master_isotope; i++)
	{
		master_isotope[i] =
			(struct master_isotope *) free_check_null(master_isotope[i]);
	}
	master_isotope =
		(struct master_isotope **) free_check_null(master_isotope);
	hdestroy_multi(master_isotope_hash_table);
	master_isotope_hash_table = NULL;

	/* calculate_value */
	for (i = 0; i < count_calculate_value; i++)
	{
		calculate_value_free(calculate_value[i]);
		calculate_value[i] =
			(struct calculate_value *) free_check_null(calculate_value[i]);
	}
	calculate_value =
		(struct calculate_value **) free_check_null(calculate_value);
	hdestroy_multi(calculate_value_hash_table);
	calculate_value_hash_table = NULL;

	/* isotope_ratio */
	for (i = 0; i < count_isotope_ratio; i++)
	{
		isotope_ratio[i] =
			(struct isotope_ratio *) free_check_null(isotope_ratio[i]);
	}
	isotope_ratio = (struct isotope_ratio **) free_check_null(isotope_ratio);
	hdestroy_multi(isotope_ratio_hash_table);
	isotope_ratio_hash_table = NULL;

	/* isotope_alpha */
	for (i = 0; i < count_isotope_alpha; i++)
	{
		isotope_alpha[i] =
			(struct isotope_alpha *) free_check_null(isotope_alpha[i]);
	}
	isotope_alpha = (struct isotope_alpha **) free_check_null(isotope_alpha);
	hdestroy_multi(isotope_alpha_hash_table);
	isotope_alpha_hash_table = NULL;

	free_tally_table();

	/* CVODE memory */
	free_cvode();

	pitzer_clean_up();

	sit_clean_up();


/* hash tables */

	hdestroy_multi(elements_hash_table);
	hdestroy_multi(species_hash_table);
	hdestroy_multi(logk_hash_table);
	hdestroy_multi(phases_hash_table);

	elements_hash_table = NULL;
	species_hash_table = NULL;
	logk_hash_table = NULL;
	phases_hash_table = NULL;

/* strings */
#ifdef HASH
	strings_hash_clear();
#else
	strings_map_clear();
#endif

/* delete basic interpreter */
	basic_free();
	change_surf = (struct Change_Surf *) free_check_null(change_surf);

/* miscellaneous work space */

	elt_list = (struct elt_list *) free_check_null(elt_list);
	trxn.token = (struct rxn_token_temp *) free_check_null(trxn.token);
	mb_unknowns = (struct unknown_list *) free_check_null(mb_unknowns);
	line = (char *) free_check_null(line);
	line_save = (char *) free_check_null(line_save);

	zeros = (LDBLE *) free_check_null(zeros);
	scratch = (LDBLE *) free_check_null(scratch);
	x_arg = (LDBLE *) free_check_null(x_arg);
	res_arg = (LDBLE *) free_check_null(res_arg);

	normal = (LDBLE *) free_check_null(normal);
	ineq_array = (LDBLE *) free_check_null(ineq_array);
	back_eq = (int *) free_check_null(back_eq);
	zero = (LDBLE *) free_check_null(zero);
	res = (LDBLE *) free_check_null(res);
	delta1 = (LDBLE *) free_check_null(delta1);
	cu = (LDBLE *) free_check_null(cu);
	iu = (int *) free_check_null(iu);
	is = (int *) free_check_null(is);

/*  x_arg = res_arg = scratch = NULL; */
	x_arg_max = res_arg_max = scratch_max = 0;


/* free user database name if defined */
	user_database = (char *) free_check_null(user_database);
	//selected_output_file_name =
	//	(char *) free_check_null(selected_output_file_name);
	dump_file_name = (char *) free_check_null(dump_file_name);
#ifdef PHREEQCI_GUI
	free_spread();
#endif
	title_x = (char *) free_check_null(title_x);

	count_elements = 0;
	count_master = 0;
	count_phases = 0;
	count_s = 0;
	count_logk = 0;
	count_master_isotope = 0;
	count_rates = 0;
	count_inverse = 0;
	count_save_values = 0;

	llnl_count_temp = 0;
	llnl_count_adh = 0;
	llnl_count_bdh = 0;
	llnl_count_bdot = 0;
	llnl_count_co2_coefs = 0;

	count_calculate_value = 0;
	count_isotope_ratio = 0;
	count_isotope_alpha = 0;

	default_data_base = (char *) free_check_null(default_data_base);
	sformatf_buffer = (char *) free_check_null(sformatf_buffer);
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
reinitialize(void)
/* ---------------------------------------------------------------------- */
{
/* solutions */
	Rxn_solution_map.clear();

/* surfaces */
	Rxn_surface_map.clear();

/* exchange */
	Rxn_exchange_map.clear();

/* pp assemblages */
	Rxn_pp_assemblage_map.clear();

/* s_s assemblages */
	Rxn_ss_assemblage_map.clear();

/* gases */
	Rxn_gas_phase_map.clear();

/* kinetics */
	Rxn_kinetics_map.clear();

/* irreversible reactions */
	Rxn_reaction_map.clear();

	// Temperature
	Rxn_temperature_map.clear();

	// Pressure
	Rxn_pressure_map.clear();
	return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "element"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
int Phreeqc::
element_compare(const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
	const struct element *element_ptr1, *element_ptr2;
	element_ptr1 = *(const struct element **) ptr1;
	element_ptr2 = *(const struct element **) ptr2;
/*      return(strcmp_nocase(element_ptr1->name, element_ptr2->name)); */
	return (strcmp(element_ptr1->name, element_ptr2->name));

}

/* ---------------------------------------------------------------------- */
struct element * Phreeqc::
element_store(const char *element)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "element" in the hash table for elements.
 *
 *   If found, pointer to the appropriate element structure is returned.
 *
 *   If the string is not found, a new entry is made at the end of
 *   the elements array (position count_elements) and count_elements is
 *   incremented. A new entry is made in the hash table. Pointer to
 *   the new structure is returned.
 *
 *   Arguments:
 *      element    input, character string to be located or stored.
 *
 *   Returns:
 *      The address of an elt structure that contains the element data.
 */
	int n;
	struct element *elts_ptr;
	ENTRY item, *found_item;
	char token[MAX_LENGTH];
/*
 *   Search list
 */
	strcpy(token, element);

	item.key = token;
	item.data = NULL;
	found_item = hsearch_multi(elements_hash_table, item, FIND);
	if (found_item != NULL)
	{
		elts_ptr = (struct element *) (found_item->data);
		return (elts_ptr);
	}
/*
 *   Save new elt structure and return pointer to it
 */
	/* make sure there is space in elements */
	elements[count_elements] =
		(struct element *) PHRQ_malloc((size_t) sizeof(struct element));
	if (elements[count_elements] == NULL)
		malloc_error();
	/* set name pointer in elements structure */
	elements[count_elements]->name = string_hsave(token);
	/* set return value */
	elements[count_elements]->master = NULL;
	elements[count_elements]->primary = NULL;
	elements[count_elements]->gfw = 0.0;
	n = count_elements++;
	if (count_elements >= max_elements)
	{
		space((void **) ((void *) &elements), count_elements, &max_elements,
			  sizeof(struct element *));
	}
/*
 *   Update hash table
 */
	item.key = elements[n]->name;
	item.data = (void *) elements[n];
	found_item = hsearch_multi(elements_hash_table, item, ENTER);
	if (found_item == NULL)
	{
		error_string = sformatf( "Hash table error in element_store.");
		error_msg(error_string, CONTINUE);
	}
	return (elements[n]);
}

/* **********************************************************************
 *
 *   Routines related to structure "elt_list"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
int Phreeqc::
elt_list_combine(void)
/* ---------------------------------------------------------------------- */
/*
 *      Function goes through the list of elements pointed to by elt_list
 *      and combines the coefficients of elements that are the same.
 *      Assumes elt_list has been sorted by element name.
 */
{
	int i, j;

	if (count_elts < 1)
	{
		output_msg("elt_list_combine: How did this happen?\n");
		return (ERROR);
	}
	if (count_elts == 1)
	{
		return (OK);
	}
	j = 0;
	for (i = 1; i < count_elts; i++)
	{
		if (elt_list[i].elt == elt_list[j].elt)
		{
			elt_list[j].coef += elt_list[i].coef;
		}
		else
		{
			j++;
			if (i != j)
			{
				elt_list[j].elt = elt_list[i].elt;
				elt_list[j].coef = elt_list[i].coef;
			}
		}
	}
	count_elts = j + 1;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
elt_list_compare(const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
	const struct elt_list *a, *b;

	a = (const struct elt_list *) ptr1;
	b = (const struct elt_list *) ptr2;
	return (strncmp(a->elt->name, b->elt->name, MAX_LENGTH));
}

/* ---------------------------------------------------------------------- */
struct elt_list * Phreeqc::
elt_list_dup(struct elt_list *elt_list_ptr_old)
/* ---------------------------------------------------------------------- */
{
/*
 *  Duplicates the elt_list structure pointed to by elt_list_ptr_old.
 */
	int i, count_totals;
	struct elt_list *elt_list_ptr_new;
/*
 *   Count totals data and copy
 */
	if (elt_list_ptr_old == NULL)
		return (NULL);
	for (i = 0; elt_list_ptr_old[i].elt != NULL; i++);
	count_totals = i;
/*
 *   Malloc space and store element data
 */
	elt_list_ptr_new =
		(struct elt_list *) PHRQ_malloc((size_t) (count_totals + 1) *
										sizeof(struct elt_list));
	if (elt_list_ptr_new == NULL)
		malloc_error();
	memcpy(elt_list_ptr_new, elt_list_ptr_old,
		   (size_t) (count_totals + 1) * sizeof(struct elt_list));
	return (elt_list_ptr_new);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
elt_list_print(struct elt_list *elt_list_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *  Duplicates the elt_list structure pointed to by elt_list_ptr_old.
 */
	int i;
/*
 *   Debug print for element list
 */
	if (elt_list_ptr == NULL)
		return (ERROR);
	output_msg(sformatf( "Elt_list\n"));
	for (i = 0; elt_list_ptr[i].elt != NULL; i++)
	{
		output_msg(sformatf( "\t%s\t%e\n", elt_list_ptr[i].elt->name,
				   (double) elt_list_ptr[i].coef));
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
cxxNameDouble Phreeqc::
elt_list_NameDouble(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Takes data from work space elt_list, makes NameDouble
 */
	cxxNameDouble nd;
	for(int i = 0; i < count_elts; i++)
	{
		nd.add(elt_list[i].elt->name, elt_list[i].coef);
	}
	return (nd);
}
/* ---------------------------------------------------------------------- */
struct elt_list * Phreeqc::
elt_list_save(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Takes data from work space elt_list, allocates a new elt_list structure,
 *   copies data from work space to new structure, and returns pointer to
 *   new structure.
 */
	int j;
	struct elt_list *elt_list_ptr;
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
 *   Malloc space and store element data
 */
	elt_list_ptr =
		(struct elt_list *) PHRQ_malloc((size_t) (count_elts + 1) *
										sizeof(struct elt_list));
	if (elt_list_ptr == NULL)
	{
		malloc_error();
	}
	else
	{
		for (j = 0; j < count_elts; j++)
		{
			elt_list_ptr[j].elt = elt_list[j].elt;
			elt_list_ptr[j].coef = elt_list[j].coef;
		}
		elt_list_ptr[count_elts].elt = NULL;
	}
	return (elt_list_ptr);
}
/* ---------------------------------------------------------------------- */
struct elt_list * Phreeqc::
NameDouble2elt_list(const cxxNameDouble &nd)
/* ---------------------------------------------------------------------- */
{
/*
 *   Takes NameDouble allocates space and fills new elt_list struct
 */
	struct elt_list *elt_list_ptr = (struct elt_list *) PHRQ_malloc((nd.size() + 1) * sizeof(struct elt_list));
	if (elt_list_ptr == NULL)
	{
		malloc_error();
	}
	else
	{
		cxxNameDouble::const_iterator it = nd.begin();
		int i = 0;
		for( ; it != nd.end(); it++)
		{
			elt_list_ptr[i].elt = element_store(it->first.c_str());
			elt_list_ptr[i].coef = it->second;
			i++;
		}
		elt_list_ptr[i].elt = NULL;
		elt_list_ptr[i].coef = 0;
	}
	return (elt_list_ptr);
}
/* **********************************************************************
 *
 *   Routines related to structure "inverse"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct inverse * Phreeqc::
inverse_alloc(void)
/* ---------------------------------------------------------------------- */
/*
 *   Allocates space for a new inverse structure at position count_inverse.
 *   Initializes structure.
 *      arguments
 *      input:  none
 *      output: pointer to an inverse structure
 *      return: OK
 */
{
	struct inverse *inverse_ptr = NULL;

	count_inverse++;
	inverse =
		(struct inverse *) PHRQ_realloc(inverse,
										(size_t) count_inverse *
										sizeof(struct inverse));
	if (inverse == NULL)
	{
		malloc_error();
		return inverse_ptr;
	}
	inverse_ptr = &(inverse[count_inverse - 1]);
/*
 *   Initialize variables
 */
	inverse_ptr->description = NULL;
	inverse_ptr->count_uncertainties = 0;
	inverse_ptr->count_solns = 0;
	inverse_ptr->count_elts = 0;
	inverse_ptr->count_isotopes = 0;
	inverse_ptr->count_i_u = 0;
	inverse_ptr->count_phases = 0;
	inverse_ptr->count_force_solns = 0;
/*
 *   allocate space for pointers in structure to NULL
 */

	inverse_ptr->uncertainties =
		(LDBLE *) PHRQ_malloc((size_t) sizeof(LDBLE));
	if (inverse_ptr->uncertainties == NULL)
	{
		malloc_error();
		return inverse_ptr;
	}

	inverse_ptr->ph_uncertainties =
		(LDBLE *) PHRQ_malloc((size_t) sizeof(LDBLE));
	if (inverse_ptr->ph_uncertainties == NULL)
	{
		malloc_error();
		return inverse_ptr;
	}

	inverse_ptr->force_solns = (int *) PHRQ_malloc((size_t) sizeof(int));
	if (inverse_ptr->force_solns == NULL)
	{
		malloc_error();
		return inverse_ptr;
	}

	inverse_ptr->dalk_dph = NULL;
	inverse_ptr->dalk_dc = NULL;

	inverse_ptr->solns = NULL;

	inverse_ptr->elts =
		(struct inv_elts *) PHRQ_malloc((size_t) sizeof(struct inv_elts));
	if (inverse_ptr->elts == NULL)
	{
		malloc_error();
		return inverse_ptr;
	}
	inverse_ptr->elts[0].name = NULL;
	inverse_ptr->elts[0].uncertainties = NULL;

	inverse_ptr->isotopes =
		(struct inv_isotope *) PHRQ_malloc((size_t)
										   sizeof(struct inv_isotope));
	if (inverse_ptr->isotopes == NULL)
	{
		malloc_error();
		return inverse_ptr;
	}
	inverse_ptr->isotopes[0].isotope_name = NULL;
	inverse_ptr->isotopes[0].isotope_number = 0;
	inverse_ptr->isotopes[0].elt_name = NULL;

	inverse_ptr->i_u =
		(struct inv_isotope *) PHRQ_malloc((size_t)
										   sizeof(struct inv_isotope));
	if (inverse_ptr->i_u == NULL)
	{
		malloc_error();
		return inverse_ptr;
	}
	inverse_ptr->i_u[0].isotope_name = NULL;
	inverse_ptr->i_u[0].isotope_number = 0;
	inverse_ptr->i_u[0].elt_name = NULL;

	inverse_ptr->phases =
		(struct inv_phases *) PHRQ_malloc((size_t) sizeof(struct inv_phases));
	if (inverse_ptr->phases == NULL)
	{
		malloc_error();
		return inverse_ptr;
	}

	return (inverse_ptr);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
inverse_compare(const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compare inverse values for n_user
 */
	const struct inverse *nptr1;
	const struct inverse *nptr2;

	nptr1 = (const struct inverse *) ptr1;
	nptr2 = (const struct inverse *) ptr2;
	if (nptr1->n_user > nptr2->n_user)
		return (1);
	if (nptr1->n_user < nptr2->n_user)
		return (-1);
	return (0);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
inverse_delete(int i)
/* ---------------------------------------------------------------------- */
{
/*
 *   Deletes inverse i from list (i is not user number),
 *   Frees memory allocated to inverse struct
 *   Input: i, number of inverse struct to delete
 *   Return: OK
 */
	int j;

	inverse_free(&(inverse[i]));
	for (j = i; j < (count_inverse - 1); j++)
	{
		memcpy((void *) &(inverse[j]), (void *) &(inverse[j + 1]),
			   (size_t) sizeof(struct inverse));
	}
	count_inverse--;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
inverse_free(struct inverse *inverse_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Free all memory for an inverse structure.
 */
	int i;

	inverse_ptr->description =
		(char *) free_check_null(inverse_ptr->description);
/*   Free solns */
	inverse_ptr->solns = (int *) free_check_null(inverse_ptr->solns);

/*   Free uncertainties */
	inverse_ptr->uncertainties =
		(LDBLE *) free_check_null(inverse_ptr->uncertainties);
	inverse_ptr->ph_uncertainties =
		(LDBLE *) free_check_null(inverse_ptr->ph_uncertainties);

/*   Free force_solns */
	inverse_ptr->force_solns =
		(int *) free_check_null(inverse_ptr->force_solns);

/*   Free elts */
	for (i = 0; i < inverse_ptr->count_elts; i++)
	{
		inverse_ptr->elts[i].uncertainties =
			(LDBLE *) free_check_null(inverse_ptr->elts[i].uncertainties);
	};
	inverse_ptr->elts =
		(struct inv_elts *) free_check_null(inverse_ptr->elts);

/*   Free isotopes */
	for (i = 0; i < inverse_ptr->count_isotopes; i++)
	{
		inverse_ptr->isotopes[i].uncertainties =
			(LDBLE *) free_check_null(inverse_ptr->isotopes[i].uncertainties);
	};
	inverse_ptr->isotopes =
		(struct inv_isotope *) free_check_null(inverse_ptr->isotopes);

	for (i = 0; i < inverse_ptr->count_i_u; i++)
	{
		inverse_ptr->i_u[i].uncertainties =
			(LDBLE *) free_check_null(inverse_ptr->i_u[i].uncertainties);
	};
	inverse_ptr->i_u =
		(struct inv_isotope *) free_check_null(inverse_ptr->i_u);

/*   Free phases */
	for (i = 0; i < inverse_ptr->count_phases; i++)
	{
		inverse_ptr->phases[i].isotopes =
			(struct isotope *) free_check_null(inverse_ptr->phases[i].
											   isotopes);
	}
	inverse_ptr->phases =
		(struct inv_phases *) free_check_null(inverse_ptr->phases);

/*   Free carbon derivatives */
	inverse_ptr->dalk_dph = (LDBLE *) free_check_null(inverse_ptr->dalk_dph);
	inverse_ptr->dalk_dc = (LDBLE *) free_check_null(inverse_ptr->dalk_dc);

	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
inverse_isotope_compare(const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
	int i;
	const struct inv_isotope *iso_ptr1, *iso_ptr2;

	iso_ptr1 = (const struct inv_isotope *) ptr1;
	iso_ptr2 = (const struct inv_isotope *) ptr2;
	i = strcmp_nocase(iso_ptr1->elt_name, iso_ptr2->elt_name);
	if (i != 0)
		return (i);
	if (iso_ptr1->isotope_number < iso_ptr2->isotope_number)
	{
		return (-1);
	}
	else if (iso_ptr1->isotope_number > iso_ptr2->isotope_number)
	{
		return (1);
	}
	return (0);
}

/* ---------------------------------------------------------------------- */
struct inverse * Phreeqc::
inverse_search(int n_user, int *n)
/* ---------------------------------------------------------------------- */
{
/*   Linear search of the structure array "inverse" for user number n_user.
 *
 *   Arguments:
 *      n_user  input, user number
 *      n       output, position in inverse
 *
 *   Returns:
 *      if found, the address of the inverse element
 *      if not found, NULL
 *
 */
	int i;
	for (i = 0; i < count_inverse; i++)
	{
		if (inverse[i].n_user == n_user)
		{
			*n = i;
			return (&(inverse[i]));
		}
	}
/*
 *   An inverse structure with n_user was not found
 */
	return (NULL);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
inverse_sort(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sort array of inverse structures
 */
	if (count_inverse > 0)
	{
		qsort(inverse, (size_t) count_inverse,
			  (size_t) sizeof(struct inverse), inverse_compare);
	}
	return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "master"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct master * Phreeqc::
master_alloc(void)
/* ---------------------------------------------------------------------- */
/*
 *   Allocates space to a master structure and initializes the space.
 *      arguments: void
 *      return: pointer to a master structure
 */
{
	struct master *ptr;
	ptr = (struct master *) PHRQ_malloc(sizeof(struct master));
	if (ptr == NULL)
		malloc_error();
/*
 *   set pointers in structure to NULL
 */
	ptr->in = FALSE;
	ptr->number = -1;
	ptr->last_model = -1;
	ptr->type = 0;
	ptr->primary = FALSE;
	ptr->coef = 0.0;
	ptr->total = 0.0;
	ptr->isotope_ratio = 0;
	ptr->isotope_ratio_uncertainty = 0;
	ptr->isotope = 0;
	ptr->total_primary = 0;
	ptr->elt = NULL;
	ptr->alk = 0.0;
	ptr->gfw = 0.0;
	ptr->gfw_formula = NULL;
	ptr->unknown = NULL;
	ptr->s = NULL;
	ptr->rxn_primary = NULL;
	ptr->rxn_secondary = NULL;
	ptr->pe_rxn = NULL;
	ptr->minor_isotope = FALSE;
	return (ptr);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
master_delete(char *ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Delete master species:  Free memory of master species structure, free
 *   the structure, and remove from array master.
 *
 *   Input
 *	ptr  character string with name of element or valence state
 *   Returns
 *	TRUE if master species was deleted.
 *	FALSE if master species was not found.
 */
	int j, n;

	if (master_search(ptr, &n) == NULL)
		return (FALSE);
	master_free(master[n]);
	for (j = n; j < (count_master - 1); j++)
	{
		master[j] = master[j + 1];
	}
	count_master--;
	return (TRUE);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
master_free(struct master *master_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Free memory pointed to by master species pointer, master_ptr.
 *   Frees master_ptr itself.
 */
	if (master_ptr == NULL)
		return (ERROR);
	rxn_free(master_ptr->rxn_primary);
	rxn_free(master_ptr->rxn_secondary);
	master_ptr = (struct master *) free_check_null(master_ptr);
	return (OK);
}

/* ---------------------------------------------------------------------- */
struct master * Phreeqc::
master_bsearch(const char *ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Uses binary search. Assumes master is in sort order.
 *   Find master species for string (*ptr) containing name of element or valence state.
 *
 *   Input: ptr    pointer to string containing element name
 *
 *   Return: pointer to master structure containing name ptr or NULL.
 */
	void *void_ptr;
	if (count_master == 0)
	{
		return (NULL);
	}
	void_ptr = bsearch((const char *) ptr,
					   (char *) master,
					   (unsigned) count_master,
					   sizeof(struct master *), master_compare_string);
	if (void_ptr == NULL)
	{
		char * dup = string_duplicate(ptr);
		replace("(+","(", dup);
		void_ptr = bsearch((const char *) dup,
			(char *) master,
			(unsigned) count_master,
			sizeof(struct master *), master_compare_string);
		dup = (char *) free_check_null(dup);
	}
	if (void_ptr == NULL)
	{
		return (NULL);
	}
	else
	{
		return (*(struct master **) void_ptr);
	}
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
master_compare_string(const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
	const char *string_ptr;
	const struct master *master_ptr;

	string_ptr = (const char *) ptr1;
	master_ptr = *(const struct master **) ptr2;
	return (strcmp_nocase(string_ptr, master_ptr->elt->name));
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
master_compare(const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
	const struct master *master_ptr1, *master_ptr2;
	master_ptr1 = *(const struct master **) ptr1;
	master_ptr2 = *(const struct master **) ptr2;
	return (strcmp_nocase(master_ptr1->elt->name, master_ptr2->elt->name));
}

/* ---------------------------------------------------------------------- */
struct master * Phreeqc::
master_bsearch_primary(const char *ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Find primary master species for first element in the string, ptr.
 *   Uses binary search. Assumes master is in sort order.
 */
	int l;
	char *ptr1;
	char elt[MAX_LENGTH];
	struct master *master_ptr_primary;
/*
 *   Find element name
 */
	char * temp_name = string_duplicate(ptr);
	ptr1 = temp_name;
	get_elt(&ptr1, elt, &l);
	free_check_null(temp_name);
/*
 *   Search master species list
 */
	master_ptr_primary = master_bsearch(elt);
	if (master_ptr_primary == NULL)
	{
		input_error++;
		error_string = sformatf(
				"Could not find primary master species for %s.", ptr);
		error_msg(error_string, CONTINUE);
	}
	return (master_ptr_primary);
}
/* ---------------------------------------------------------------------- */
struct master * Phreeqc::
master_bsearch_secondary(char *ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Find secondary master species that corresponds to the primary master species.
 *   i.e. S(6) for S.
 */
	int l;
	char *ptr1;
	char elt[MAX_LENGTH];
	struct master *master_ptr_primary, *master_ptr=NULL, *master_ptr_secondary=NULL;
	int j;
/*
 *   Find element name
 */
	ptr1 = ptr;
	get_elt(&ptr1, elt, &l);
/*
 *   Search master species list
 */
	master_ptr_primary = master_bsearch(elt);
	if (master_ptr_primary == NULL)
	{
		input_error++;
		error_string = sformatf(
				"Could not find primary master species for %s.", ptr);
		error_msg(error_string, CONTINUE);
	}
/*
 *  If last in list or not redox
*/
	if (master_ptr_primary)
	{
		if ((master_ptr_primary->number >= count_master - 1) || 
			(master[master_ptr_primary->number + 1]->elt->primary != master_ptr_primary))
		{
			return(master_ptr_primary);
		}
		/*
		*  Find secondary master with same species as primary
		*/
		master_ptr = NULL;
		for (j = master_ptr_primary->number + 1; j < count_master; j++)
		{
			if (master[j]->s == master_ptr_primary->s)
			{
				master_ptr = master[j];
			}
		}
	}
/*
 *
 */
	if (master_ptr != NULL && master_ptr->elt != NULL && (master_ptr->elt->primary == master_ptr_primary))
	{
		master_ptr_secondary = master_ptr;
	}
	else
	{		
		input_error++;
		error_string = sformatf(
				"Could not find secondary master species for %s.", ptr);
		error_msg(error_string, STOP);
	}


	return (master_ptr_secondary);
}
/* ---------------------------------------------------------------------- */
struct master * Phreeqc::
master_search(char *ptr, int *n)
/* ---------------------------------------------------------------------- */
{
/*
 *   Linear search of master to find master species in string, ptr.
 *   Returns pointer if found. n contains position in array master.
 *   Returns NULL if not found.
 */
	int i;
	struct master *master_ptr;
/*
 *   Search master species list
 */
	*n = -999;
	for (i = 0; i < count_master; i++)
	{
		if (strcmp(ptr, master[i]->elt->name) == 0)
		{
			*n = i;
			master_ptr = master[i];
			return (master_ptr);
		}
	}
	return (NULL);
}
/* **********************************************************************
 *
 *   Routines related to structure "phases"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct phase * Phreeqc::
phase_alloc(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Allocates space to a phase structure and initializes
 *      arguments: void
 *      return: pointer to new phase structure
 */
	struct phase *phase_ptr;
/*
 *   Allocate space
 */
	phase_ptr = (struct phase *) PHRQ_malloc(sizeof(struct phase));
	if (phase_ptr == NULL)
		malloc_error();
/*
 *   Initialize space
 */
	phase_init(phase_ptr);
	return (phase_ptr);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
phase_compare(const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compares names of phases for sort
 */
	const struct phase *phase_ptr1, *phase_ptr2;
	phase_ptr1 = *(const struct phase **) ptr1;
	phase_ptr2 = *(const struct phase **) ptr2;
	return (strcmp_nocase(phase_ptr1->name, phase_ptr2->name));
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
phase_compare_string(const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
	const char *char_ptr;
	const struct phase *phase_ptr;
	char_ptr = (const char *) ptr1;
	phase_ptr = *(const struct phase **) ptr2;
	return (strcmp_nocase(char_ptr, phase_ptr->name));
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
phase_delete(int i)
/* ---------------------------------------------------------------------- */
{
/*
 *   Deletes phase i from list, phases
 *   Frees memory allocated to phase[i] and renumbers phases to remove number i.
 *   Input: i, number of phase
 *   Return: OK
 */
	int j;

	phase_free(phases[i]);
	phases[i] = (struct phase *) free_check_null(phases[i]);
	for (j = i; j < (count_phases - 1); j++)
	{
		phases[j] = phases[j + 1];
	}
	count_phases--;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
phase_free(struct phase *phase_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Frees memory allocated within phase[i], does not free phase structure
 *   Input: i, number of phase
 *   Return: OK
 */
	if (phase_ptr == NULL)
		return (ERROR);
	phase_ptr->next_elt =
		(struct elt_list *) free_check_null(phase_ptr->next_elt);
	phase_ptr->next_sys_total =
		(struct elt_list *) free_check_null(phase_ptr->next_sys_total);
	rxn_free(phase_ptr->rxn);
	rxn_free(phase_ptr->rxn_s);
	rxn_free(phase_ptr->rxn_x);
	phase_ptr->add_logk =
		(struct name_coef *) free_check_null(phase_ptr->add_logk);
	return (OK);
}

/* ---------------------------------------------------------------------- */
struct phase * Phreeqc::
phase_bsearch(const char *ptr, int *j, int print)
/* ---------------------------------------------------------------------- */
{
/*   Binary search the structure array "phases" for a name that is equal to
 *   ptr. Assumes array phases is in sort order.
 *
 *   Arguments:
 *      name  input, a character string to be located in phases.
 *      j	    index number in array phases.
 *
 *   Returns:
 *      if found, pointer to phase structure.
 *      if not found, NULL
 *
 */
	void *void_ptr;

	void_ptr = NULL;
	if (count_phases > 0)
	{
		void_ptr = (void *)
			bsearch((char *) ptr,
					(char *) phases,
					(size_t) count_phases,
					(size_t) sizeof(struct phase *), phase_compare_string);
	}
	if (void_ptr == NULL && print == TRUE)
	{
		error_string = sformatf( "Could not find phase in list, %s.", ptr);
		error_msg(error_string, CONTINUE);
	}

	if (void_ptr == NULL)
	{
		*j = -1;
		return (NULL);
	}

	*j = (int) ((struct phase **) void_ptr - phases);
	return (*(struct phase **) void_ptr);
}

/* ---------------------------------------------------------------------- */
 int Phreeqc::
phase_init(struct phase *phase_ptr)
/* ---------------------------------------------------------------------- */
/*
 *   set pointers in phase structure to NULL
 */
{
	int i;

	phase_ptr->name = NULL;
	phase_ptr->formula = NULL;
	phase_ptr->in = FALSE;
	phase_ptr->lk = 0.0;
	for (i = 0; i < MAX_LOG_K_INDICES; i++)
		phase_ptr->logk[i] = 0.0;
	phase_ptr->original_units = kjoules;
	phase_ptr->count_add_logk = 0;
	phase_ptr->add_logk = NULL;
	phase_ptr->moles_x = 0;
	phase_ptr->delta_max = 0;
	phase_ptr->p_soln_x = 0;
	phase_ptr->fraction_x = 0;
	phase_ptr->log10_lambda = 0;
	phase_ptr->log10_fraction_x = 0;
	phase_ptr->dn = 0;
	phase_ptr->dnb = 0;
	phase_ptr->dnc = 0;
	phase_ptr->gn = 0;
	phase_ptr->gntot = 0;
	phase_ptr->t_c = 0.0;
	phase_ptr->p_c = 0.0;
	phase_ptr->omega = 0.0;
	phase_ptr->pr_a = 0.0;
	phase_ptr->pr_b = 0.0;
	phase_ptr->pr_alpha = 0.0;
	phase_ptr->pr_tk = 0;
	phase_ptr->pr_p = 0;
	phase_ptr->pr_phi = 1.0;
	phase_ptr->pr_aa_sum2 = 0;
	for (i = 0; i < 9; i++)
		phase_ptr->delta_v[i] = 0.0;
	phase_ptr->pr_si_f = 0;
	phase_ptr->pr_in = false;
	phase_ptr->type = SOLID;
	phase_ptr->next_elt = NULL;
	phase_ptr->next_sys_total = NULL;
	phase_ptr->check_equation = TRUE;
	phase_ptr->rxn = NULL;
	phase_ptr->rxn_s = NULL;
	phase_ptr->rxn_x = NULL;
	phase_ptr->replaced = 0;
	phase_ptr->in_system = 1;
	phase_ptr->original_deltav_units = cm3_per_mol;
	return (OK);
}

/* ---------------------------------------------------------------------- */
struct phase * Phreeqc::
phase_store(const char *name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "name" in the hash table for phases.
 *
 *   If found, pointer to the appropriate phase structure is returned.
 *
 *   If the string is not found, a new entry is made at the end of
 *   the phases array (position count_phases) and count_phases is
 *   incremented. A new entry is made in the hash table. Pointer to
 *   the new structure is returned.
 *
 *   Arguments:
 *      name    input, character string to be located or stored.
 *
 *   Returns:
 *      The address of a phase structure that contains the phase data.
 *      If phase existed, it is reinitialized. The structure returned
 *      contains only the name of the phase.
 */
	int n;
	struct phase *phase_ptr;
	ENTRY item, *found_item;
	char token[MAX_LENGTH];
	const char *ptr;
/*
 *   Search list
 */

	strcpy(token, name);
	str_tolower(token);
	ptr = string_hsave(token);

	item.key = ptr;
	item.data = NULL;
	found_item = hsearch_multi(phases_hash_table, item, FIND);
	if (found_item != NULL)
	{
		phase_ptr = (struct phase *) (found_item->data);
		phase_free(phase_ptr);
		phase_init(phase_ptr);
		phase_ptr->name = string_hsave(name);
		return (phase_ptr);
	}
/*
 *   Make new phase structure and return pointer to it
 */
	/* make sure there is space in phases */
	n = count_phases++;
	if (count_phases >= max_phases)
	{
		space((void **) ((void *) &phases), count_phases, &max_phases,
			  sizeof(struct phase *));
	}
	phases[n] = phase_alloc();
	/* set name in phase structure */
	phases[n]->name = string_hsave(name);
/*
 *   Update hash table
 */
	item.key = ptr;
	item.data = (void *) phases[n];
	found_item = hsearch_multi(phases_hash_table, item, ENTER);
	if (found_item == NULL)
	{
		error_string = sformatf( "Hash table error in phase_store.");
		error_msg(error_string, CONTINUE);
	}

	return (phases[n]);
}
/* **********************************************************************
 *
 *   Routines related to structure "rates"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct rate * Phreeqc::
rate_bsearch(char *ptr, int *j)
/* ---------------------------------------------------------------------- */
{
/*   Binary search the structure array "rates" for a name that is equal to
 *   ptr. Assumes array rates is in sort order.
 *
 *   Arguments:
 *      name  input, a character string to be located in rates.
 *      j	    index number in array rates.
 *
 *   Returns:
 *      if found, pointer to rate structure.
 *      if not found, NULL
 *
 */
	void *void_ptr;

	if (count_rates == 0)
	{
		*j = -1;
		return (NULL);
	}
	void_ptr = (void *)
		bsearch((char *) ptr,
				(char *) rates,
				(size_t) count_rates,
				(size_t) sizeof(struct rate *), rate_compare_string);

	if (void_ptr == NULL)
	{
		*j = -1;
		return (NULL);
	}

	*j = (int) ((struct rate *) void_ptr - rates);
	return ((struct rate *) void_ptr);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
rate_compare(const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compares names of rates for sort
 */
	const struct rate *rate_ptr1, *rate_ptr2;
	rate_ptr1 = *(const struct rate **) ptr1;
	rate_ptr2 = *(const struct rate **) ptr2;
	return (strcmp_nocase(rate_ptr1->name, rate_ptr2->name));
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
rate_compare_string(const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
	const char *char_ptr;
	const struct rate *rate_ptr;
	char_ptr = (const char *) ptr1;
	rate_ptr = *(const struct rate **) ptr2;
	return (strcmp_nocase(char_ptr, rate_ptr->name));
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
rate_free(struct rate *rate_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Frees memory allocated within rate[i], does not free rate structure
 *   Input: i, number of rate
 *   Return: OK
 */
	

	if (rate_ptr == NULL)
		return (ERROR);
	rate_ptr->commands = (char *) free_check_null(rate_ptr->commands);
	if (rate_ptr->linebase != NULL)
	{
		char cmd[] = "new; quit";
		basic_run(cmd, rate_ptr->linebase, rate_ptr->varbase, rate_ptr->loopbase);
		rate_ptr->linebase = NULL;
		rate_ptr->varbase = NULL;
		rate_ptr->loopbase = NULL;
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
struct rate * Phreeqc::
rate_search(const char *name_in, int *n)
/* ---------------------------------------------------------------------- */
{
/*   Linear search of the structure array "rates" for name.
 *
 *   Arguments:
 *     name     input, name of rate
 *      n       output, position in rates
 *
 *   Returns:
 *      if found, the address of the pp_assemblage element
 *      if not found, NULL
 */
	std::map<const char *, int>::iterator it;

	const char * name;
	name = string_hsave(name_in);

	it = rates_map.find(name);
	if (it != rates_map.end())
	{
		*n = it->second;
		if (*n >= 0)
		{
			return &(rates[it->second]);
		}
		return NULL;
	}

	int i;
	*n = -1;
	for (i = 0; i < count_rates; i++)
	{
		if (strcmp_nocase(rates[i].name, name) == 0)
		{
			*n = i;
			rates_map[name] = i;
			return (&(rates[i]));
		}
	}
/*
 *   rate name not found
 */
	rates_map[name] = *n;
	return (NULL);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
rate_sort(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sort array of rate structures
 */
	if (count_rates > 0)
	{
		qsort(rates, (size_t) count_rates, (size_t) sizeof(struct rate),
			  rate_compare);
	}
	return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "reaction", balanced chemical reactions
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct reaction * Phreeqc::
rxn_alloc(int ntokens)
/* ---------------------------------------------------------------------- */
{
	int i;
/*
 *   Allocates space to a rxn structure
 *      input: ntokens, number of tokens in reaction
 *      return: pointer to a species structure
 */
	struct reaction *rxn_ptr;
/*
 *   Malloc reaction structure
 */
	rxn_ptr = (struct reaction *) PHRQ_malloc(sizeof(struct reaction));
	if (rxn_ptr == NULL)
		malloc_error();
/*
 *   zero log k data
 */
	for (i = 0; i < MAX_LOG_K_INDICES; i++)
	{
		rxn_ptr->logk[i] = 0.0;
	}
/*
 *   zero dz data
 */
	for (i = 0; i < 3; i++)
	{
		rxn_ptr->dz[i] = 0.0;
	}
/*
 *   Malloc rxn_token structure
 */
	rxn_ptr->token =
		(struct rxn_token *) PHRQ_malloc((size_t) ntokens *
										 sizeof(struct rxn_token));
	for (i = 0; i < ntokens; i++)
	{
		rxn_ptr->token[i].s = NULL;
		rxn_ptr->token[i].name = NULL;
		rxn_ptr->token[i].coef = 0.0;
	}

	if (rxn_ptr->token == NULL)
		malloc_error();
	return (rxn_ptr);
}

/* ---------------------------------------------------------------------- */
struct reaction * Phreeqc::
rxn_dup(struct reaction *rxn_ptr_old)
/* ---------------------------------------------------------------------- */
{
/*
 *   mallocs space for a reaction and copies the reaction
 *   input: rxn_ptr_old, pointer to a reaction structure to copy
 *
 *   Return: rxn_ptr_new,  pointer to duplicated structure to copy
 */
	int i;
	struct reaction *rxn_ptr_new;

	if (rxn_ptr_old == NULL)
		return (NULL);
	for (i = 0; rxn_ptr_old->token[i].s != NULL; i++);

	rxn_ptr_new = rxn_alloc(i + 1);
/*
 *   Copy logk data
 */
	memcpy(rxn_ptr_new->logk, rxn_ptr_old->logk, (size_t) MAX_LOG_K_INDICES * sizeof(LDBLE));
/*
 *   Copy dz data
 */
	memcpy(rxn_ptr_new->dz, rxn_ptr_old->dz, (size_t) (3 * sizeof(LDBLE)));
/*
 *   Copy tokens
 */
	memcpy(rxn_ptr_new->token, rxn_ptr_old->token,
		   (size_t) (i + 1) * sizeof(struct rxn_token));

	return (rxn_ptr_new);
}
/* ---------------------------------------------------------------------- */
struct reaction * Phreeqc::
cxxChemRxn2rxn(cxxChemRxn &cr)
/* ---------------------------------------------------------------------- */
{
/*
 *   mallocs space for a reaction and copies the cxxChemRxn to a struct reaction
 *
 *   Return: rxn_ptr_new,  pointer to new structure 
 */
	for (int i = 0; i < (int) cr.Get_tokens().size(); i++)
	{
		if (cr.Get_tokens()[i].s != NULL)
		{
			cr.Get_tokens()[i].s = s_store(cr.Get_tokens()[i].s->name, cr.Get_tokens()[i].s->z, FALSE);
		}
		if (cr.Get_tokens()[i].name != NULL)
		{
			cr.Get_tokens()[i].name = string_hsave(cr.Get_tokens()[i].name);
		}
		else
		{
			if (cr.Get_tokens()[i].s != NULL)
			{
				cr.Get_tokens()[i].name = string_hsave(cr.Get_tokens()[i].s->name);
			}
			else
			{
				cr.Get_tokens()[i].name=NULL;
			}
		}
	}

	count_trxn = 0;
	trxn_add(cr, 1.0, 1);

	struct reaction *rxn_ptr_new = rxn_alloc(count_trxn + 1);
	trxn_copy(rxn_ptr_new);

	// cleanup pointers for copy operator name, and s may point into another instance

	for (int i = 0; rxn_ptr_new->token[i].s != NULL; i++)
	{
		rxn_ptr_new->token[i].name = string_hsave(rxn_ptr_new->token[i].name);
		LDBLE  z = rxn_ptr_new->token[i].s->z;
		rxn_ptr_new->token[i].s = s_store(rxn_ptr_new->token[i].name, z, false);
	}
	return (rxn_ptr_new);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
rxn_find_coef(struct reaction * r_ptr, const char *str)
/* ---------------------------------------------------------------------- */
{
/*
 *   Finds coefficient of token in reaction.
 *   input: r_ptr, pointer to a reaction structure
 *	  str, string to find as reaction token
 *
 *   Return: 0.0, if token not found
 *	   coefficient of token, if found.
 */
	struct rxn_token *r_token;
	LDBLE coef;

	r_token = r_ptr->token + 1;
	coef = 0.0;
	while (r_token->s != NULL)
	{
		if (strcmp(r_token->s->name, str) == 0)
		{
			coef = r_token->coef;
			break;
		}
		r_token++;
	}
	return (coef);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
rxn_free(struct reaction *rxn_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Frees space allocated for a reaction structure
 *      input: rxn_ptr, pointer to reaction structure
 *      return: ERROR, if pointer is NULL
 *	      OK, otherwise.
 */
	if (rxn_ptr == NULL)
		return (ERROR);
	rxn_ptr->token = (struct rxn_token *) free_check_null(rxn_ptr->token);
	rxn_ptr = (struct reaction *) free_check_null(rxn_ptr);
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
rxn_print(struct reaction *rxn_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Frees space allocated for a reaction structure
 *      input: rxn_ptr, pointer to reaction structure
 *      return: ERROR, if pointer is NULL
 *	      OK, otherwise.
 */
	struct rxn_token *next_token;
	int i;
	if (rxn_ptr == NULL)
		return (ERROR);
	next_token = rxn_ptr->token;
	output_msg(sformatf( "log k data:\n"));
	for (i = 0; i < MAX_LOG_K_INDICES; i++)
	{
		output_msg(sformatf( "\t%f\n", (double) rxn_ptr->logk[i]));
	}
	output_msg(sformatf( "Reaction definition\n"));
	while (next_token->s != NULL || next_token->name != NULL)
	{
		output_msg(sformatf( "\tcoef %f ", next_token->coef));
		if (next_token->s != NULL)
		{
			output_msg(sformatf( "\tspecies token: %s ",
					   next_token->s->name));
		}
		if (next_token->name != NULL)
		{
			output_msg(sformatf( "\tname token: %s", next_token->name));
		}
		output_msg(sformatf( "\n"));
		next_token++;
	}
	output_msg(sformatf( "dz data\n"));
	for (i = 0; i < 3; i++)
	  {
	    output_msg(sformatf( "\t%d %e\n", i, (double) rxn_ptr->dz[i]));
	    
	  }
	return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "species"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct species * Phreeqc::
s_alloc(void)
/* ---------------------------------------------------------------------- */
/*
 *   Allocates space to a species structure, initializes
 *      arguments: void
 *      return: pointer to a species structure
 */
{
	struct species *s_ptr;
	s_ptr = (struct species *) PHRQ_malloc(sizeof(struct species));
	if (s_ptr == NULL)
		malloc_error();
/*
 *   set pointers in structure to NULL, variables to zero
 */
	s_init(s_ptr);

	return (s_ptr);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
s_compare(const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
	const struct species *s_ptr1, *s_ptr2;
	s_ptr1 = *(const struct species **) ptr1;
	s_ptr2 = *(const struct species **) ptr2;
	return (strcmp(s_ptr1->name, s_ptr2->name));

}

/* ---------------------------------------------------------------------- */
int Phreeqc::
s_delete(int i)
/* ---------------------------------------------------------------------- */
{
/*
 *   Delete species i: free memory and renumber array of pointers, s.
 */
	int j;

	s_free(s[i]);
	s[i] = (struct species *) free_check_null(s[i]);
	for (j = i; j < (count_s - 1); j++)
	{
		s[j] = s[j + 1];
	}
	count_s--;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
s_free(struct species *s_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Free space allocated for species structure, s_ptr. Does not free s_ptr.
 */
	if (s_ptr == NULL)
		return (ERROR);
	s_ptr->next_elt = (struct elt_list *) free_check_null(s_ptr->next_elt);
	s_ptr->next_secondary =
		(struct elt_list *) free_check_null(s_ptr->next_secondary);
	s_ptr->next_sys_total =
		(struct elt_list *) free_check_null(s_ptr->next_sys_total);
	s_ptr->add_logk = (struct name_coef *) free_check_null(s_ptr->add_logk);
	rxn_free(s_ptr->rxn);
	rxn_free(s_ptr->rxn_s);
	rxn_free(s_ptr->rxn_x);
	return (OK);
}

/* ---------------------------------------------------------------------- */
 int Phreeqc::
s_init(struct species *s_ptr)
/* ---------------------------------------------------------------------- */
/*
 *      return: pointer to a species structure
 */
{
	int i;
/*
 *   set pointers in structure to NULL
 */
	s_ptr->name = NULL;
	s_ptr->mole_balance = NULL;
	s_ptr->in = FALSE;
	s_ptr->number = 0;
	s_ptr->primary = NULL;
	s_ptr->secondary = NULL;
	s_ptr->gfw = 0.0;
	s_ptr->z = 0.0;
	s_ptr->dw = 0.0;
	s_ptr->erm_ddl = 1.0;
	s_ptr->equiv = 0;
	s_ptr->alk = 0.0;
	s_ptr->carbon = 0.0;
	s_ptr->co2 = 0.0;
	s_ptr->h = 0.0;
	s_ptr->o = 0.0;
	s_ptr->dha = 0.0;
	s_ptr->dhb = 0.0;
	s_ptr->a_f = 0.0;
	s_ptr->lk = 0.0;
	for (i = 0; i < MAX_LOG_K_INDICES; i++)
	{
		s_ptr->logk[i] = 0.0;
	}
/* VP: Density Start */
	for (i = 0; i < 6; i++)
	{
		s_ptr->millero[i] = 0.0;
	}
/* VP: Density End */
	s_ptr->original_units = kjoules;
	s_ptr->count_add_logk = 0;
	s_ptr->add_logk = NULL;
	s_ptr->lg = 0.0;
	s_ptr->lg_pitzer = 0.0;
	s_ptr->lm = 0.0;
	s_ptr->la = 0.0;
	s_ptr->dg = 0.0;
	s_ptr->dg_total_g = 0;
	s_ptr->moles = 0.0;
	s_ptr->type = 0;
	s_ptr->gflag = 0;
	s_ptr->exch_gflag = 0;
	s_ptr->next_elt = NULL;
	s_ptr->next_secondary = NULL;
	s_ptr->next_sys_total = NULL;
	s_ptr->check_equation = TRUE;
	s_ptr->original_deltav_units = cm3_per_mol;

	s_ptr->rxn = NULL;
	s_ptr->rxn_s = NULL;
	s_ptr->rxn_x = NULL;
	s_ptr->tot_g_moles = 0;
	s_ptr->tot_dh2o_moles = 0;
	for (i = 0; i < 5; i++)
	{
		s_ptr->cd_music[i] = 0.0;
	}
	for (i = 0; i < 3; i++)
	{
		s_ptr->dz[i] = 0.0;
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
struct species * Phreeqc::
s_search(const char *name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "name" in the hash table for species.
 *
 *   Arguments:
 *      name  input, a character string to be located in species.
 *      i    is obsolete.
 *
 *   Returns:
 *   If found, pointer to the appropriate species structure is returned.
 *       else, NULL pointer is returned.
 */
	struct species *s_ptr;
	ENTRY item, *found_item;
	char safe_name[MAX_LENGTH];

	strcpy(safe_name, name);
	item.key = safe_name;
	item.data = NULL;
	found_item = hsearch_multi(species_hash_table, item, FIND);
	if (found_item != NULL)
	{
		s_ptr = (struct species *) (found_item->data);
		return (s_ptr);
	}
	return (NULL);
}

/* ---------------------------------------------------------------------- */
struct species * Phreeqc::
s_store(const char *name, LDBLE l_z, int replace_if_found)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "name" in the hash table for species.
 *
 *   Pointer to a species structure is always returned.
 *
 *   If the string is not found, a new entry is made at the end of
 *      the elements array (position count_elements) and count_elements is
 *      incremented. A new entry is made in the hash table. Pointer to
 *      the new structure is returned.
 *   If "name" is found and replace is true, pointers in old species structure
 *      are freed and replaced with additional input.
 *   If "name" is found and replace is false, the old species structure is not
 *      modified and a pointer to it is returned.
 *
 *   Arguments:
 *      name    input, character string to be found in "species".
 *      l_z      input, charge on "name"
 *      replace_if_found input, TRUE means reinitialize species if found
 *		     FALSE means just return pointer if found.
 *
 *   Returns:
 *      pointer to species structure "s" where "name" can be found.
 */
	int n;
	struct species *s_ptr;
	ENTRY item, *found_item;
/*
 *   Search list
 */
	item.key = name;
	item.data = NULL;
	found_item = hsearch_multi(species_hash_table, item, FIND);

	if (found_item != NULL && replace_if_found == FALSE)
	{
		s_ptr = (struct species *) (found_item->data);
		return (s_ptr);
	}
	else if (found_item != NULL && replace_if_found == TRUE)
	{
		s_ptr = (struct species *) (found_item->data);
		s_free(s_ptr);
		s_init(s_ptr);
	}
	else
	{
		n = count_s++;
		/* make sure there is space in s */
		if (count_s >= max_s)
		{
			space((void **) ((void *) &s), count_s, &max_s,
				  sizeof(struct species *));
		}
		/* Make new species structure */
		s[n] = s_alloc();
		s_ptr = s[n];
	}
	/* set name and z in pointer in species structure */
	s_ptr->name = string_hsave(name);
	s_ptr->z = l_z;
/*
 *   Update hash table
 */
	item.key = s_ptr->name;
	item.data = (void *) s_ptr;
	found_item = hsearch_multi(species_hash_table, item, ENTER);
	if (found_item == NULL)
	{
		error_string = sformatf( "Hash table error in species_store.");
		error_msg(error_string, CONTINUE);
	}

	return (s_ptr);
}
/* **********************************************************************
 *
 *   Routines related to structure "save_values"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct save_values * Phreeqc::
save_values_bsearch(struct save_values *k, int *n)
/* ---------------------------------------------------------------------- */
{
/*
 *   Binary search save_values to find if one exists with given coefficients
 *   Save_Values is assumed to be in sort order by count_subscripts and
 *   values of subscripts
 */
	void *void_ptr;
	if (count_save_values == 0)
	{
		*n = -999;
		return (NULL);
	}
	void_ptr = (void *)
		bsearch((char *) k,
				(char *) save_values,
				(size_t) count_save_values,
				(size_t) sizeof(struct save_values), save_values_compare);
	if (void_ptr == NULL)
	{
		*n = -999;
		return (NULL);
	}
	*n = (int) ((struct save_values *) void_ptr - save_values);
	return ((struct save_values *) void_ptr);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
save_values_compare(const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
	int i;
	const struct save_values *save_values_ptr1, *save_values_ptr2;
	save_values_ptr1 = (const struct save_values *) ptr1;
	save_values_ptr2 = (const struct save_values *) ptr2;
	if (save_values_ptr1->count_subscripts <
		save_values_ptr2->count_subscripts)
	{
		return (-1);
	}
	else if (save_values_ptr1->count_subscripts >
			 save_values_ptr2->count_subscripts)
	{
		return (1);
	}
	else
	{
		for (i = 0; i < save_values_ptr1->count_subscripts; i++)
		{
			if (save_values_ptr1->subscripts[i] <
				save_values_ptr2->subscripts[i])
			{
				return (-1);
			}
			else if (save_values_ptr1->subscripts[i] >
					 save_values_ptr2->subscripts[i])
			{
				return (1);
			}
		}
	}
	return (0);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
save_values_sort(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sort array of save_values structures
 */
	if (count_save_values > 0)
	{
		qsort(save_values, (size_t) count_save_values,
			  (size_t) sizeof(struct save_values), save_values_compare);
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
save_values_store(struct save_values *s_v)
/* ---------------------------------------------------------------------- */
{
/*
 *   Look for subscripts
 */
	int n, i;
	struct save_values *s_v_ptr;

	s_v_ptr = save_values_bsearch(s_v, &n);
	if (s_v_ptr != NULL)
	{
		s_v_ptr->value = s_v->value;
	}
	else
	{
		save_values =
			(struct save_values *) PHRQ_realloc(save_values,
												(size_t) (count_save_values +
														  1) *
												sizeof(struct save_values));
		if (save_values == NULL)
			malloc_error();
		save_values[count_save_values].value = s_v->value;
		save_values[count_save_values].count_subscripts =
			s_v->count_subscripts;
		i = s_v->count_subscripts;
		if (i == 0)
			i = 1;
		save_values[count_save_values].subscripts =
			(int *) PHRQ_malloc((size_t) i * sizeof(int));
		if (save_values[count_save_values].subscripts == NULL)
			malloc_error();
		save_values[count_save_values].subscripts =
			(int *) memcpy(save_values[count_save_values].subscripts,
						   s_v->subscripts, (size_t) i * sizeof(int));
		count_save_values++;
		save_values_sort();
	}

	if (count_save_values > 0)
	{
		qsort(save_values, (size_t) count_save_values,
			  (size_t) sizeof(struct save_values), save_values_compare);
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
isotope_compare(const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
	int i;
	const struct isotope *iso_ptr1, *iso_ptr2;

	iso_ptr1 = (const struct isotope *) ptr1;
	iso_ptr2 = (const struct isotope *) ptr2;
	i = strcmp_nocase(iso_ptr1->elt_name, iso_ptr2->elt_name);
	if (i != 0)
		return (i);
	if (iso_ptr1->isotope_number < iso_ptr2->isotope_number)
	{
		return (-1);
	}
	else if (iso_ptr1->isotope_number > iso_ptr2->isotope_number)
	{
		return (1);
	}
	return (0);
}
/* **********************************************************************
 *
 *   Routines related to structure "species_list"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
 int Phreeqc::
species_list_compare(const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
	int j;
	const char *name1, *name2;
	const struct species_list *nptr1, *nptr2;

	nptr1 = (const struct species_list *) ptr1;
	nptr2 = (const struct species_list *) ptr2;

/*
 *   Put H+ first
 */
	if (nptr1->master_s != nptr2->master_s)
	{
		/*
		if (nptr1->master_s == s_hplus)
			return (-1);
		if (nptr2->master_s == s_hplus)
			return (1);
		*/
		if ((strcmp(nptr1->master_s->name,"H+") == 0) || (strcmp(nptr1->master_s->name,"H3O+") == 0))
			return (-1);
		if ((strcmp(nptr2->master_s->name,"H+") == 0) || (strcmp(nptr2->master_s->name,"H3O+") == 0))
			return (1);
	}
/*
 *   Other element valence states
 */
	if (nptr1->master_s->secondary != NULL)
	{
		name1 = nptr1->master_s->secondary->elt->name;
	}
	else
	{
		name1 = nptr1->master_s->primary->elt->name;
	}
	if (nptr2->master_s->secondary != NULL)
	{
		name2 = nptr2->master_s->secondary->elt->name;
	}
	else
	{
		name2 = nptr2->master_s->primary->elt->name;
	}
/*
 *   Compare name of primary or secondary master species; log molality
 */

	j = strcmp(name1, name2);

/*
 *   Different master species
 */
	if (j != 0)
		return (j);

/*
 *   Else, descending order by log molality
 */
	if (nptr1->s->lm > nptr2->s->lm)
	{
		return (-1);
	}
	else if (nptr1->s->lm < nptr2->s->lm)
	{
		return (1);
	}
	else
	{
		return (0);
	}
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
species_list_compare_alk(const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
	const struct species_list *nptr1, *nptr2;
	LDBLE alk1, alk2;

	nptr1 = (const struct species_list *) ptr1;
	nptr2 = (const struct species_list *) ptr2;
/*
 *   Else, descending order by log molality
 */
	alk1 = fabs(under(nptr1->s->lm) * nptr1->s->alk);
	alk2 = fabs(under(nptr2->s->lm) * nptr2->s->alk);

	if (alk1 > alk2)
	{
		return (-1);
	}
	else if (alk1 < alk2)
	{
		return (1);
	}
	else
	{
		return (0);
	}
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
species_list_compare_master(const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
	const char *name1, *name2;
	const struct species_list *nptr1, *nptr2;

	nptr1 = (const struct species_list *) ptr1;
	nptr2 = (const struct species_list *) ptr2;

/*
 *   Put H+ first
 */
	if (nptr1->master_s != nptr2->master_s)
	{
		/*
		if (nptr1->master_s == s_hplus)
			return (-1);
		if (nptr2->master_s == s_hplus)
			return (1);
		*/
		if ((strcmp(nptr1->master_s->name,"H+") == 0) || (strcmp(nptr1->master_s->name,"H3O+") == 0))
			return (-1);
		if ((strcmp(nptr2->master_s->name,"H+") == 0) || (strcmp(nptr2->master_s->name,"H3O+") == 0))
			return (1);
	}
/*
 *   Other element valence states
 */
	if (nptr1->master_s->secondary != NULL)
	{
		name1 = nptr1->master_s->secondary->elt->name;
	}
	else
	{
		name1 = nptr1->master_s->primary->elt->name;
	}
	if (nptr2->master_s->secondary != NULL)
	{
		name2 = nptr2->master_s->secondary->elt->name;
	}
	else
	{
		name2 = nptr2->master_s->primary->elt->name;
	}
/*
 *   Compare name of primary or secondary master species; log molality
 */

	return (strcmp(name1, name2));
}


/* ---------------------------------------------------------------------- */
int Phreeqc::
species_list_sort(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sort list using rules in species_list_compare
 */
	if (count_species_list > 0)
	{
		qsort(&species_list[0], (size_t) count_species_list,
			  (size_t) sizeof(struct species_list), species_list_compare);
	}
	return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "surface"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct Change_Surf * Phreeqc::
change_surf_alloc(int count)
/* ---------------------------------------------------------------------- */
{
	if (count == 1)
		return (change_surf);
	change_surf =
		(struct Change_Surf *) PHRQ_realloc(change_surf,
											(size_t) count *
											sizeof(struct Change_Surf));
	if (change_surf == NULL)
		malloc_error();
	change_surf[count - 1].cell_no = -99;
	change_surf[count - 1].next = FALSE;
	change_surf[count - 2].next = TRUE;

	return (change_surf);
}
/* ---------------------------------------------------------------------- */
struct master * Phreeqc::
surface_get_psi_master(const char *name, int plane)
/* ---------------------------------------------------------------------- */
{
	struct master *master_ptr;
	std::string token;

	if (name == NULL)
		return (NULL);
	token = name;
	token.append("_psi");
	switch (plane)
	{
	case SURF_PSI:
		break;
	case SURF_PSI1:
		token.append("b");
		break;
	case SURF_PSI2:
		token.append("d");
		break;
	default:
		error_msg("Unknown plane for surface_get_psi_master", STOP);
	}
	master_ptr = master_bsearch(token.c_str());
	return (master_ptr);
}
/* **********************************************************************
 *
 *   Routines related to structure "trxn"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
int Phreeqc::
rxn_token_temp_compare(const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
	const struct rxn_token_temp *rxn_token_temp_ptr1, *rxn_token_temp_ptr2;
	rxn_token_temp_ptr1 = (const struct rxn_token_temp *) ptr1;
	rxn_token_temp_ptr2 = (const struct rxn_token_temp *) ptr2;
	return (strcmp(rxn_token_temp_ptr1->name, rxn_token_temp_ptr2->name));
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
trxn_add(cxxChemRxn &r_ptr, LDBLE coef, int combine)
/* ---------------------------------------------------------------------- */
{
/*
 *   Adds reactions together.
 *
 *   Global variable count_trxn determines which position in trxn is used.
 *      If count_trxn=0, then the equation effectively is copied into trxn.
 *      If count_trxn>0, then new equation is added to existing equation.
 *
 *   Arguments:
 *      *r_ptr	 points to rxn structure to add.
 *
 *       coef	  added equation is multiplied by coef.
 *       combine       if TRUE, reaction is reaction is sorted and
 *		     like terms combined.
 */
/*
 *   Accumulate log k for reaction
 */
	if (count_trxn == 0)
	{
		for (int i = 0; i < MAX_LOG_K_INDICES; i++)
		{
			trxn.logk[i] = r_ptr.Get_logk()[i];
		}
		for (int i = 0; i < 3; i++)
		{
			trxn.dz[i] = r_ptr.Get_dz()[i];
		}
	}
	else
	{
		for (int i = 0; i < MAX_LOG_K_INDICES; i++)
		{
			trxn.logk[i] += coef * (r_ptr.Get_logk()[i]);
		}
		for (int i = 0; i < 3; i++)
		{
			trxn.dz[i] += coef * r_ptr.Get_dz()[i];
		}
	}
/*
 *   Copy  equation into work space
 */
	for (size_t j = 0; j < r_ptr.Get_tokens().size(); j++)
	{
		if (count_trxn + 1 >= max_trxn)
		{
			space((void **) ((void *) &(trxn.token)), count_trxn + 1,
				  &max_trxn, sizeof(struct rxn_token_temp));
		}
		trxn.token[count_trxn].name = r_ptr.Get_tokens()[j].name;
		trxn.token[count_trxn].s = r_ptr.Get_tokens()[j].s;
		trxn.token[count_trxn].coef = coef * r_ptr.Get_tokens()[j].coef;
		count_trxn++;
	}
	if (combine == TRUE)
		trxn_combine();
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
trxn_add(struct reaction *r_ptr, LDBLE coef, int combine)
/* ---------------------------------------------------------------------- */
{
/*
 *   Adds reactions together.
 *
 *   Global variable count_trxn determines which position in trxn is used.
 *      If count_trxn=0, then the equation effectively is copied into trxn.
 *      If count_trxn>0, then new equation is added to existing equation.
 *
 *   Arguments:
 *      *r_ptr	 points to rxn structure to add.
 *
 *       coef	  added equation is multiplied by coef.
 *       combine       if TRUE, reaction is reaction is sorted and
 *		     like terms combined.
 */
	int i;
	struct rxn_token *next_token;
/*
 *   Accumulate log k for reaction
 */
	if (count_trxn == 0)
	{
		memcpy((void *) trxn.logk, (void *) r_ptr->logk,
			(size_t) MAX_LOG_K_INDICES * sizeof(LDBLE));
		for (i = 0; i < 3; i++)
		{
			trxn.dz[i] = r_ptr->dz[i];
		}
	}
	else
	{
		for (i = 0; i < MAX_LOG_K_INDICES; i++)
		{
			trxn.logk[i] += coef * (r_ptr->logk[i]);
		}
		for (i = 0; i < 3; i++)
		{
			trxn.dz[i] += coef * r_ptr->dz[i];
		}
	}
/*
 *   Copy  equation into work space
 */
	next_token = r_ptr->token;
	while (next_token->s != NULL)
	{
		if (count_trxn + 1 >= max_trxn)
		{
			space((void **) ((void *) &(trxn.token)), count_trxn + 1,
				  &max_trxn, sizeof(struct rxn_token_temp));
		}
		trxn.token[count_trxn].name = next_token->s->name;
		trxn.token[count_trxn].s = next_token->s;
		trxn.token[count_trxn].coef = coef * next_token->coef;
		count_trxn++;
		next_token++;
	}
	if (combine == TRUE)
		trxn_combine();
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
trxn_add_phase(struct reaction *r_ptr, LDBLE coef, int combine)
/* ---------------------------------------------------------------------- */
{
/*
 *   Adds reactions together.
 *
 *   Global variable count_trxn determines which position in trxn is used.
 *      If count_trxn=0, then the equation effectively is copied into trxn.
 *      If count_trxn>0, then new equation is added to existing equation.
 *
 *   Arguments:
 *      *r_ptr	 points to rxn structure to add.
 *
 *       coef	  added equation is multiplied by coef.
 *       combine       if TRUE, reaction is reaction is sorted and
 *		     like terms combined.
 */
	int i;
	struct rxn_token *next_token;
/*
 *   Accumulate log k for reaction
 */
	if (count_trxn == 0)
	{
		memcpy((void *) trxn.logk, (void *) r_ptr->logk,
			(size_t) MAX_LOG_K_INDICES * sizeof(LDBLE));
	}
	else
	{
		for (i = 0; i < MAX_LOG_K_INDICES; i++)
		{
			trxn.logk[i] += coef * (r_ptr->logk[i]);
		}
	}
/*
 *   Copy  equation into work space
 */
	next_token = r_ptr->token;
	while (next_token->s != NULL || next_token->name != NULL)
	{
		if (count_trxn + 1 >= max_trxn)
		{
			space((void **) ((void *) &(trxn.token)), count_trxn + 1,
				  &max_trxn, sizeof(struct rxn_token_temp));
		}
		if (next_token->s != NULL)
		{
			trxn.token[count_trxn].name = next_token->s->name;
			trxn.token[count_trxn].s = next_token->s;
		}
		else
		{
			trxn.token[count_trxn].name = next_token->name;
			trxn.token[count_trxn].s = NULL;
		}
		trxn.token[count_trxn].coef = coef * next_token->coef;
		count_trxn++;
		next_token++;
	}
	if (combine == TRUE)
		trxn_combine();
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
trxn_combine(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Combines coefficients of tokens that are equal in temporary
 *   reaction structure, trxn.
 */
	int j, k;
/*
 *   Sort trxn species
 */
	trxn_sort();
/*
 *   Combine trxn tokens
 */
	j = 1;
	for (k = 2; k < count_trxn; k++)
	{
		if (trxn.token[k].s != NULL)
		{
			if ((j > 0) && (trxn.token[k].s == trxn.token[j].s))
			{
				trxn.token[j].coef += trxn.token[k].coef;
				if (equal(trxn.token[j].coef, 0.0, 1e-5))
					j--;
			}
			else
			{
				j++;
				if (k != j)
				{
					trxn.token[j].name = trxn.token[k].name;
					trxn.token[j].s = trxn.token[k].s;
					trxn.token[j].coef = trxn.token[k].coef;
				}
			}
		}
		else
		{
			if ((j > 0) && (trxn.token[k].s == trxn.token[j].s)
				&& (trxn.token[k].name == trxn.token[j].name))
			{
				trxn.token[j].coef += trxn.token[k].coef;
				if (equal(trxn.token[j].coef, 0.0, 1e-5))
					j--;
			}
			else
			{
				j++;
				if (k != j)
				{
					trxn.token[j].name = trxn.token[k].name;
					trxn.token[j].s = trxn.token[k].s;
					trxn.token[j].coef = trxn.token[k].coef;
				}
			}
		}
	}
	count_trxn = j + 1;			/* number excluding final NULL */
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
trxn_copy(struct reaction *rxn_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies trxn to a reaction structure.
 *
 *   Input: rxn_ptr, pointer to reaction structure to copy trxn to.
 *
 */
	int i;
/*
 *   Copy logk data
 */
	for (i = 0; i < MAX_LOG_K_INDICES; i++)
	{
		rxn_ptr->logk[i] = trxn.logk[i];
	}
/*
 *   Copy dz data
 */
	for (i = 0; i < 3; i++)
	{
		rxn_ptr->dz[i] = trxn.dz[i];
	}
/*
 *   Copy tokens
 */
	for (i = 0; i < count_trxn; i++)
	{
		rxn_ptr->token[i].s = trxn.token[i].s;
		rxn_ptr->token[i].name = trxn.token[i].name;
		rxn_ptr->token[i].coef = trxn.token[i].coef;
	}
	rxn_ptr->token[count_trxn].s = NULL;

	return (OK);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
trxn_find_coef(const char *str, int start)
/* ---------------------------------------------------------------------- */
{
/*
 *   Finds coefficient of specified token in trxn.
 *   Input: str, token name in reaction.
 *
 *   Return: 0.0, if token not found.
 *	   coefficient of token, if token found.
 */
	int i;
	LDBLE coef;

	coef = 0.0;
	for (i = start; i < count_trxn; i++)
	{
		if (strcmp(trxn.token[i].s->name, str) == 0)
		{
			coef = trxn.token[i].coef;
			break;
		}
	}
	return (coef);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
trxn_multiply(LDBLE coef)
/* ---------------------------------------------------------------------- */
{
/*
 *   Multiplies temporary reaction, trxn,  by a constant
 *
 *   Arguments:
 *       input: coef	  multiplier.
 */
	int i;
/*
 *   Multiply log k for reaction
 */
	for (i = 0; i < MAX_LOG_K_INDICES; i++)
	{
		trxn.logk[i] *= coef;
	}
/*
 *   Multiply dz for reaction
 */
	for (i = 0; i < 3; i++)
	{
		trxn.dz[i] *= coef;
	}
/*
 *   Multiply coefficients of reaction
 */
	for (i = 0; i < count_trxn; i++)
	{
		trxn.token[i].coef *= coef;
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
trxn_print(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints trxn
 */
	int i;
/*
 *   Print log k for reaction
 */

	output_msg(sformatf( "\tlog k data:\n"));
	for (i = 0; i < MAX_LOG_K_INDICES; i++)
	{
		output_msg(sformatf( "\t\t%f\n", (double) trxn.logk[i]));
	}

/*
 *   Print dz for reaction
 */
	output_msg(sformatf( "\tdz data:\n"));
	for (i = 0; i < 3; i++)
	{
		output_msg(sformatf( "\t\t%f\n", (double) trxn.dz[i]));
	}
/*
 *   Print stoichiometry
 */
	output_msg(sformatf( "\tReaction stoichiometry\n"));
	for (i = 0; i < count_trxn; i++)
	{
		output_msg(sformatf( "\t\t%-20s\t%10.2f\n", trxn.token[i].name,
				   (double) trxn.token[i].coef));
	}
	output_msg(sformatf( "\n"));
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
trxn_reverse_k(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Changes K from dissociation to association and back
 */
	int i;
/*
 *   Accumulate log k for reaction
 */
   for (i = 0; i < MAX_LOG_K_INDICES; i++)
	{
		trxn.logk[i] = -trxn.logk[i];
	}
	for (i = 0; i < 3; i++)
	{
		trxn.dz[i] = -trxn.dz[i];
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
trxn_sort(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compare names in tokens in trxn array for sorting
 */
	if (count_trxn - 1 > 0)
	{
		qsort(&trxn.token[1],
			  (size_t) count_trxn - 1,
			  (size_t) sizeof(struct rxn_token_temp), rxn_token_temp_compare);
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
trxn_swap(const char *token)
/* ---------------------------------------------------------------------- */
{
/*
 *   Moves specified token to initial position in reaction.
 *   Input: token, token name to move to initial position.
 *
 *   Return: ERROR, if token not found.
 *	   OK, if token moved to initial position.
 */
	int i, j;
	LDBLE coef;
/*
 *   Locate token
 */
	for (j = 0; j < count_trxn; j++)
	{
		if (strcmp(trxn.token[j].s->name, token) == 0)
			break;
	}
	if (j >= count_trxn)
	{
		input_error++;
		error_string = sformatf( "Could not find token in equation, %s.", token);
		error_msg(error_string, CONTINUE);
		for (i = 0; i < count_trxn; i++)
		{
			output_msg(sformatf( "%f\t%s\t",
					   (double) trxn.token[i].coef, trxn.token[i].name));
		}
		output_msg(sformatf( "\n"));
		return (ERROR);
	}
/*
 *   Swap token to first position
 */
	trxn.token[count_trxn].name = trxn.token[0].name;
	trxn.token[count_trxn].s = trxn.token[0].s;
	trxn.token[count_trxn].coef = trxn.token[0].coef;

	trxn.token[0].name = trxn.token[j].name;
	trxn.token[0].s = trxn.token[j].s;
	trxn.token[0].coef = trxn.token[j].coef;

	trxn.token[j].name = trxn.token[count_trxn].name;
	trxn.token[j].s = trxn.token[count_trxn].s;
	trxn.token[j].coef = trxn.token[count_trxn].coef;
/*
 *   Make coefficient of token -1.0
 */
	coef = -1.0 / trxn.token[0].coef;
	trxn_multiply(coef);
	return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "unknown"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct unknown * Phreeqc::
unknown_alloc(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Allocates space to an "unknown" structure
 *      arguments: void
 *      return: pointer to an "unknown" structure
 */
	struct unknown *unknown_ptr;
/*
 *   Allocate space
 */
	unknown_ptr = (struct unknown *) PHRQ_malloc(sizeof(struct unknown));
	if (unknown_ptr == NULL)
		malloc_error();
/*
 *   set pointers in structure to NULL
 */
	unknown_ptr->type = 0;
	unknown_ptr->moles = 0.0;
	unknown_ptr->ln_moles = 0.0;
	unknown_ptr->f = 0.0;
	unknown_ptr->sum = 0.0;
	unknown_ptr->delta = 0.0;
	unknown_ptr->la = 0.0;
	unknown_ptr->number = 0;
	unknown_ptr->description = NULL;
	unknown_ptr->master = NULL;
	unknown_ptr->phase = NULL;
	unknown_ptr->si = 0.0;
	unknown_ptr->s = NULL;
	unknown_ptr->exch_comp = NULL;
	unknown_ptr->pp_assemblage_comp_ptr = NULL;
	unknown_ptr->ss_name = NULL;
	unknown_ptr->ss_ptr = NULL;
	unknown_ptr->ss_comp_name = NULL;
	unknown_ptr->ss_comp_ptr = NULL;
	unknown_ptr->ss_comp_number = 0;
	unknown_ptr->ss_in = FALSE;
	unknown_ptr->surface_comp = NULL;
	unknown_ptr->related_moles = 0.0;
	unknown_ptr->potential_unknown = NULL;
	unknown_ptr->potential_unknown1 = NULL;
	unknown_ptr->potential_unknown2 = NULL;
	unknown_ptr->count_comp_unknowns = 0;
	unknown_ptr->comp_unknowns = NULL;
	unknown_ptr->phase_unknown = NULL;
	unknown_ptr->surface_charge = NULL;
	unknown_ptr->mass_water = 0.0;
	unknown_ptr->dissolve_only = FALSE;

	return (unknown_ptr);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
unknown_delete(int i)
/* ---------------------------------------------------------------------- */
{
/*
 *   Delete unknow from list x
 */
	int j;

	unknown_free(x[i]);
	for (j = i; j < (count_unknowns); j++)
	{
		x[j] = x[j + 1];
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
unknown_free(struct unknown *unknown_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Frees space allocated to an unknown structure, frees unknown_ptr.
 */
	if (unknown_ptr == NULL)
		return (ERROR);
	unknown_ptr->master =
		(struct master **) free_check_null(unknown_ptr->master);
	if (unknown_ptr->type == SURFACE_CB)
	{
		/*
		   surface_charge_free(unknown_ptr->surface_charge);
		   unknown_ptr->surface_charge = (struct surface_charge *) free_check_null(unknown_ptr->surface_charge);
		 */
	}
	unknown_ptr->comp_unknowns =
		(struct unknown **) free_check_null(unknown_ptr->comp_unknowns);
	unknown_ptr = (struct unknown *) free_check_null(unknown_ptr);
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
system_duplicate(int i, int save_old)
/* ---------------------------------------------------------------------- */
{
	Utilities::Rxn_copy(Rxn_solution_map, i, save_old);

	Utilities::Rxn_copy(Rxn_pp_assemblage_map, i, save_old);

	Utilities::Rxn_copy(Rxn_exchange_map, i, save_old);

	Utilities::Rxn_copy(Rxn_surface_map, i, save_old);

	Utilities::Rxn_copy(Rxn_gas_phase_map, i, save_old);

	Utilities::Rxn_copy(Rxn_kinetics_map, i, save_old);

	Utilities::Rxn_copy(Rxn_ss_assemblage_map, i, save_old);

	return (OK);
}

/* ---------------------------------------------------------------------- */
struct logk * Phreeqc::
logk_store(char *name, int replace_if_found)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "name" in the hash table for logk.
 *
 *   Pointer to a logk structure is always returned.
 *
 *   If the string is not found, a new entry is made in the hash table. Pointer to
 *      the new structure is returned.
 *   If "name" is found and replace is true, pointers in old logk structure
 *      are freed and replaced with additional input.
 *   If "name" is found and replace is false, the old logk structure is not
 *      modified and a pointer to it is returned.
 *
 *   Arguments:
 *      name    input, character string to be found in "logk".
 *      replace_if_found input, TRUE means reinitialize logk structure if found
 *		     FALSE means just return pointer if found.
 *
 *   Returns:
 *      pointer to logk structure "logk" where "name" can be found.
 */
	int n;
	struct logk *logk_ptr;
	ENTRY item, *found_item;
/*
 *   Search list
 */
	str_tolower(name);
	item.key = name;
	item.data = NULL;
	found_item = hsearch_multi(logk_hash_table, item, FIND);

	if (found_item != NULL && replace_if_found == FALSE)
	{
		logk_ptr = (struct logk *) (found_item->data);
		return (logk_ptr);
	}
	else if (found_item != NULL && replace_if_found == TRUE)
	{
		logk_ptr = (struct logk *) (found_item->data);
		logk_init(logk_ptr);
	}
	else
	{
		n = count_logk++;
		/* make sure there is space in s */
		if (count_logk >= max_logk)
		{
			space((void **) ((void *) &logk), count_logk, &max_logk,
				  sizeof(struct logk *));
		}
		/* Make new logk structure */
		logk[n] = logk_alloc();
		logk_ptr = logk[n];
	}
	/* set name and z in pointer in logk structure */
	logk_ptr->name = string_hsave(name);
/*
 *   Update hash table
 */
	item.key = logk_ptr->name;
	item.data = (void *) logk_ptr;
	found_item = hsearch_multi(logk_hash_table, item, ENTER);
	if (found_item == NULL)
	{
		error_string = sformatf( "Hash table error in logk_store.");
		error_msg(error_string, CONTINUE);
	}

	return (logk_ptr);
}

/* ---------------------------------------------------------------------- */
struct logk * Phreeqc::
logk_alloc(void)
/* ---------------------------------------------------------------------- */
/*
 *   Allocates space to a logk structure, initializes
 *      arguments: void
 *      return: pointer to a logk structure
 */
{
	struct logk *logk_ptr;
	logk_ptr = (struct logk *) PHRQ_malloc(sizeof(struct logk));
	if (logk_ptr == NULL)
		malloc_error();
/*
 *   set pointers in structure to NULL, variables to zero
 */
	logk_init(logk_ptr);

	return (logk_ptr);
}

/* ---------------------------------------------------------------------- */
 int Phreeqc::
logk_init(struct logk *logk_ptr)
/* ---------------------------------------------------------------------- */
/*
 *      return: pointer to a logk structure
 */
{
	int i;
/*
 *   set pointers in structure to NULL
 */
	logk_ptr->name = NULL;
/*
 *   set varibles = 0
 */
	logk_ptr->lk = 0.0;
	for (i = 0; i < MAX_LOG_K_INDICES; i++)
	{
		logk_ptr->log_k[i] = 0.0;
		logk_ptr->log_k_original[i] = 0.0;
	}
	logk_ptr->count_add_logk = 0;
	logk_ptr->add_logk = NULL;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
logk_copy2orig(struct logk *logk_ptr)
/* ---------------------------------------------------------------------- */
/*
 *   Copies log k data to logk_original
 */
{
	int i;
	for (i = 0; i < MAX_LOG_K_INDICES; i++)
	{
		logk_ptr->log_k_original[i] = logk_ptr->log_k[i];
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
struct logk * Phreeqc::
logk_search(const char *name_in)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "name" in the hash table for logk.
 *
 *   Arguments:
 *      name    input, character string to be found in "logk".
 *
 *   Returns:
 *      pointer to logk structure "logk" where "name" can be found.
 *      or NULL if not found.
 */
	struct logk *logk_ptr;
	ENTRY item, *found_item;
/*
 *   Search list
 */
	char * name = string_duplicate(name_in);
	str_tolower(name);
	item.key = name;
	item.data = NULL;
	found_item = hsearch_multi(logk_hash_table, item, FIND);
	free_check_null(name);
	if (found_item != NULL)
	{
		logk_ptr = (struct logk *) (found_item->data);
		return (logk_ptr);
	}
	return (NULL);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
entity_exists(const char *name, int n_user)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads solution,		   0 Solution
 *	 reaction,		   1 Reaction
 *	 exchange,		   2 Exchange
 *	 surface,		    3 Surface
 *	 gas_phase,		  4 Gas_phase
 *	 equilibrium_phases,	 5 Pure_phase
 *	 solid_solution,	     6 Ss_phase
 *	 kinetics,		   7 Kinetics
 *	 mix,			8 Mix
 *	 reaction_temperature	9 Temperature
 *	 unknown		     10 UnKnown
 */
	int return_value;
	char token[MAX_LENGTH];
	enum entity_type type;
/*
 *   Read keyword
 */
	strncpy(token, name, MAX_LENGTH-1);
	token[MAX_LENGTH-1] = '\0';
	type = get_entity_enum(token);
	return_value = TRUE;
	switch (type)
	{
	case UnKnown:
		warning_msg
			("EXISTS expecting keyword solution, mix, kinetics, reaction, reaction_temperature, equilibrium_phases, exchange, surface, gas_phase, or solid_solutions.");
		return_value = 2;
		break;
	case Solution:				/* Solution */
		if (Utilities::Rxn_find(Rxn_solution_map, n_user) == NULL)
		{
			return_value = FALSE;
		}
		break;
	case Pure_phase:			/* Pure phases */
		if (Utilities::Rxn_find(Rxn_pp_assemblage_map, n_user) == NULL)
		{
			return_value = FALSE;
		}
		break;
	case Reaction:				/* Reaction */
		if (Utilities::Rxn_find(Rxn_reaction_map, n_user) == NULL)
		{
			return_value = FALSE;
		}
		break;
	case Mix:					/* Mix */
		if (Utilities::Rxn_find(Rxn_mix_map, n_user) == NULL)
		{
			return_value = FALSE;
		}
		break;
	case Exchange:				/* Ex */
		if (Utilities::Rxn_find(Rxn_exchange_map, n_user) == NULL)
		{
			return_value = FALSE;
		}
		break;
	case Surface:				/* Surface */
		if (Utilities::Rxn_find(Rxn_surface_map, n_user) == NULL)
		{
			return_value = FALSE;
		}
		break;
	case Temperature:
		if (Utilities::Rxn_find(Rxn_temperature_map, n_user) == NULL)
		{
			return_value = FALSE;
		}
	case Pressure:
		if (Utilities::Rxn_find(Rxn_pressure_map, n_user) == NULL)
		{
			return_value = FALSE;
		}
	case Gas_phase:			/* Gas */
		if (Utilities::Rxn_find(Rxn_gas_phase_map, n_user) == NULL)
		{
			return_value = FALSE;
		}
		break;
	case Kinetics:				/* Kinetics */
		if (Utilities::Rxn_find(Rxn_kinetics_map, n_user) == NULL)
		{
			return_value = FALSE;
		}
		break;
	case Ss_phase:				/* solid_solutions */
		if (Utilities::Rxn_find(Rxn_ss_assemblage_map, n_user) == NULL)
		{
			return_value = FALSE;
		}
		break;
	}
	return (return_value);
}

/* ---------------------------------------------------------------------- */
enum entity_type Phreeqc::
get_entity_enum(char *name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads solution,		   0 Solution
 *	 reaction,		   1 Reaction
 *	 exchange,		   2 Exchange
 *	 surface,		    3 Surface
 *	 gas_phase,		  4 Gas_phase
 *	 equilibrium_phases,	 5 Pure_phase
 *	 solid_solution,	     6 Ss_phase
 *	 kinetics,		   7 Kinetics
 *	 mix,			8 Mix
 *	 reaction_temperature	9 Temperature
 *   reaction_pressure
 *	 unknown		     10 UnKnown
 *
 */
	int i;
	char *ptr;
	char token[MAX_LENGTH];
/*
 *   Read keyword
 */
	ptr = name;
	copy_token(token, &ptr, &i);
	check_key(token);

	switch (next_keyword)
	{
	case Keywords::KEY_SOLUTION:					/* Solution */
		return (Solution);
		break;
	case Keywords::KEY_EQUILIBRIUM_PHASES:		/* Pure phases */
		return (Pure_phase);
		break;
	case Keywords::KEY_REACTION:					/* Reaction */
		return (Reaction);
		break;
	case Keywords::KEY_MIX:						/* Mix */
		return (Mix);
		break;
	case Keywords::KEY_EXCHANGE:					/* Ex */
		return (Exchange);
		break;
	case Keywords::KEY_SURFACE:					/* Surface */
		return (Surface);
		break;
	case Keywords::KEY_REACTION_TEMPERATURE:		/* Temperature */
		return (Temperature);
		break;
	case Keywords::KEY_REACTION_PRESSURE:		/* Pressure */
		return (Pressure);
		break;
	case Keywords::KEY_GAS_PHASE:					/* Gas */
		return (Gas_phase);
		break;
	case Keywords::KEY_KINETICS:					/* Kinetics */
		return (Kinetics);
		break;
	case Keywords::KEY_SOLID_SOLUTIONS:			/* solid_solutions */
		return (Ss_phase);
		break;
	default:
		warning_msg
			("EXISTS expecting keyword solution, mix, kinetics, reaction, reaction_temperature, equilibrium_phases, exchange, surface, gas_phase, or solid_solutions.");
		break;
	}
	return (UnKnown);
}

/*
 * copier routines
 */
/* ---------------------------------------------------------------------- */
int Phreeqc::
copier_add(struct copier *copier_ptr, int n_user, int start, int end)
/* ---------------------------------------------------------------------- */
/*
 *   add new set of copy instructions
 */
{

	if (copier_ptr->count >= copier_ptr->max)
	{
		copier_ptr->max = copier_ptr->count * 2;
		copier_ptr->n_user =
			(int *) PHRQ_realloc(copier_ptr->n_user,
								 (size_t) (copier_ptr->max * sizeof(int)));
		if (copier_ptr->n_user == NULL)
		{
			malloc_error();
			return (OK);
		}
		copier_ptr->start =
			(int *) PHRQ_realloc(copier_ptr->start,
								 (size_t) (copier_ptr->max * sizeof(int)));
		if (copier_ptr->start == NULL)
		{
			malloc_error();
			return (OK);
		}
		copier_ptr->end =
			(int *) PHRQ_realloc(copier_ptr->end,
								 (size_t) (copier_ptr->max * sizeof(int)));
		if (copier_ptr->end == NULL)
		{
			malloc_error();
			return (OK);
		}
	}
	copier_ptr->n_user[copier_ptr->count] = n_user;
	copier_ptr->start[copier_ptr->count] = start;
	copier_ptr->end[copier_ptr->count] = end;
	copier_ptr->count++;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
copier_free(struct copier *copier_ptr)
/* ---------------------------------------------------------------------- */
/*
 *   initialize copier structure
 */
{

	copier_ptr->n_user = (int *) free_check_null(copier_ptr->n_user);
	copier_ptr->start = (int *) free_check_null(copier_ptr->start);
	copier_ptr->end = (int *) free_check_null(copier_ptr->end);
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
copier_init(struct copier *copier_ptr)
/* ---------------------------------------------------------------------- */
/*
 *   initialize copier structure
 */
{

	copier_ptr->count = 0;
	copier_ptr->max = 10;
	copier_ptr->n_user =
		(int *) PHRQ_malloc((size_t) (copier_ptr->max * sizeof(int)));
	copier_ptr->start =
		(int *) PHRQ_malloc((size_t) (copier_ptr->max * sizeof(int)));
	copier_ptr->end =
		(int *) PHRQ_malloc((size_t) (copier_ptr->max * sizeof(int)));
	return (OK);
}
#include "StorageBin.h"

void Phreeqc::
Use2cxxStorageBin(cxxStorageBin & sb)
{
	//Add everything from use structure to storagebin sb

	sb.Get_system().Set_io(sb.Get_io());
	if (use.Get_mix_in())
	{
		cxxMix *entity = Utilities::Rxn_find(Rxn_mix_map, use.Get_n_mix_user());
		if (entity != NULL)
		{
			sb.Set_Mix(use.Get_n_mix_user(), entity);
		}
	}
	else if (use.Get_solution_in())
	{
		cxxSolution *entity = Utilities::Rxn_find(Rxn_solution_map, use.Get_n_solution_user());
		if (entity != NULL)
		{
			sb.Set_Solution(use.Get_n_solution_user(), entity);
		}
	}
	if (use.Get_pp_assemblage_in())
	{
		cxxPPassemblage *entity_ptr = Utilities::Rxn_find(Rxn_pp_assemblage_map, use.Get_n_pp_assemblage_user());
		if (entity_ptr != NULL)
		{
			sb.Set_PPassemblage(use.Get_n_pp_assemblage_user(), entity_ptr);
		}
	}
	if (use.Get_exchange_in())
	{
		cxxExchange *entity_ptr = Utilities::Rxn_find(Rxn_exchange_map, use.Get_n_exchange_user());
		if (entity_ptr != NULL)
		{
			sb.Set_Exchange(use.Get_n_exchange_user(), entity_ptr);
		}
	}
	if (use.Get_surface_in())
	{
		cxxSurface *entity_ptr = Utilities::Rxn_find(Rxn_surface_map, use.Get_n_surface_user());
		if (entity_ptr != NULL)
		{
			sb.Set_Surface(use.Get_n_surface_user(), entity_ptr);
		}
	}
	if (use.Get_gas_phase_in())
	{
		cxxGasPhase *entity_ptr = Utilities::Rxn_find(Rxn_gas_phase_map, use.Get_n_gas_phase_user());
		if (entity_ptr != NULL)
		{
			sb.Set_GasPhase(use.Get_n_gas_phase_user(), entity_ptr);
		}
	}
	if (use.Get_ss_assemblage_in())
	{
		cxxSSassemblage *entity_ptr = Utilities::Rxn_find(Rxn_ss_assemblage_map, use.Get_n_ss_assemblage_user());
		if (entity_ptr != NULL)
		{
			sb.Set_SSassemblage(use.Get_n_ss_assemblage_user(), entity_ptr);
		}
	}
	if (use.Get_kinetics_in())
	{
		cxxKinetics *entity_ptr = Utilities::Rxn_find(Rxn_kinetics_map, use.Get_n_kinetics_user());
		if (entity_ptr != NULL)
		{
			sb.Set_Kinetics(use.Get_n_kinetics_user(), entity_ptr);
		}
	}
	if (use.Get_reaction_in())
	{
		cxxReaction *entity = Utilities::Rxn_find(Rxn_reaction_map, use.Get_n_reaction_user());
		if (entity != NULL)
		{
			sb.Set_Reaction(use.Get_n_reaction_user(), entity);
		}
	}
	if (use.Get_temperature_in())
	{
		cxxTemperature *entity = Utilities::Rxn_find(Rxn_temperature_map, use.Get_n_temperature_user());
		if (entity != NULL)
		{
			sb.Set_Temperature(use.Get_n_temperature_user(), entity);
		}
	}
	if (use.Get_pressure_in())
	{
		cxxPressure *entity = Utilities::Rxn_find(Rxn_pressure_map, use.Get_n_pressure_user());
		if (entity != NULL)
		{
			sb.Set_Pressure(use.Get_n_pressure_user(), entity);
		}
	}
}

void Phreeqc::
phreeqc2cxxStorageBin(cxxStorageBin & sb)
	//
	// Fills StorageBin sb with all reactants from phreeqc instance.
	// equivalent to old import_phreeqc.
	//
{
	// Solutions
	{
		std::map<int, cxxSolution>::iterator it;
		for (it = Rxn_solution_map.begin(); it != Rxn_solution_map.end(); it++)
		{
			sb.Set_Solution(it->second.Get_n_user(), &(it->second));	
		}
	}
	// Exchangers
	{
		std::map<int, cxxExchange>::iterator it;
		for (it = Rxn_exchange_map.begin(); it != Rxn_exchange_map.end(); it++)
		{
			sb.Set_Exchange(it->second.Get_n_user(), &(it->second));	
		}
	}
	// GasPhases
	{
		std::map<int, cxxGasPhase>::iterator it;
		for (it = Rxn_gas_phase_map.begin(); it != Rxn_gas_phase_map.end(); it++)
		{
			sb.Set_GasPhase(it->second.Get_n_user(), &(it->second));	
		}
	}

	// Kinetics
	{
		std::map<int, cxxKinetics>::iterator it;
		for (it = Rxn_kinetics_map.begin(); it != Rxn_kinetics_map.end(); it++)
		{
			sb.Set_Kinetics(it->second.Get_n_user(), &(it->second));	
		}
	}
	// PPassemblages
	{
		std::map<int, cxxPPassemblage>::iterator it;
		for (it = Rxn_pp_assemblage_map.begin(); it != Rxn_pp_assemblage_map.end(); it++)
		{
			sb.Set_PPassemblage(it->second.Get_n_user(), &(it->second));	
		}
	}
	// SSassemblages
	{
		std::map<int, cxxSSassemblage>::iterator it;
		for (it = Rxn_ss_assemblage_map.begin(); it != Rxn_ss_assemblage_map.end(); it++)
		{
			sb.Set_SSassemblage(it->second.Get_n_user(), &(it->second));	
		}
	}
	// Surfaces
	{
		std::map<int, cxxSurface>::iterator it;
		for (it = Rxn_surface_map.begin(); it != Rxn_surface_map.end(); it++)
		{
			sb.Set_Surface(it->second.Get_n_user(), &(it->second));	
		}
	}
	// Mixes
	{
		std::map<int, cxxMix>::iterator it;
		for (it = Rxn_mix_map.begin(); it != Rxn_mix_map.end(); it++)
		{
			sb.Set_Mix(it->second.Get_n_user(), &(it->second));	
		}
	}

	// Reactions
	{
		std::map<int, cxxReaction>::iterator it;
		for (it = Rxn_reaction_map.begin(); it != Rxn_reaction_map.end(); it++)
		{
			sb.Set_Reaction(it->second.Get_n_user(), &(it->second));	
		}
	}

	// Temperatures
	{
		std::map<int, cxxTemperature>::iterator it;
		for (it = Rxn_temperature_map.begin(); it != Rxn_temperature_map.end(); it++)
		{
			sb.Set_Temperature(it->second.Get_n_user(), &(it->second));	
		}
	}

	// Pressures
	{
		std::map<int, cxxPressure>::iterator it;
		for (it = Rxn_pressure_map.begin(); it != Rxn_pressure_map.end(); it++)
		{
			sb.Set_Pressure(it->second.Get_n_user(), &(it->second));	
		}
	}
}

void Phreeqc::
phreeqc2cxxStorageBin(cxxStorageBin & sb, int n)
		//
		// copy phreeqc reactants numbered n to StorageBin sb
		//
{
	// Solutions
	{
		cxxSolution *entity_ptr = Utilities::Rxn_find(Rxn_solution_map, n);
		if (entity_ptr != NULL)
		{
			sb.Set_Solution(n, entity_ptr);
		}
	}
	// Exchangers
	{
		cxxExchange *entity_ptr = Utilities::Rxn_find(Rxn_exchange_map, n);
		if (entity_ptr != NULL)
		{
			sb.Set_Exchange(n, entity_ptr);
		}
	}

	// GasPhases
	{
		cxxGasPhase *entity_ptr = Utilities::Rxn_find(Rxn_gas_phase_map, n);
		if (entity_ptr != NULL)
		{
			sb.Set_GasPhase(n, entity_ptr);
		}
	}

	// Kinetics
	{
		cxxKinetics *entity_ptr = Utilities::Rxn_find(Rxn_kinetics_map, n);
		if (entity_ptr != NULL)
		{
			sb.Set_Kinetics(n, entity_ptr);
		}
	}
	// PPassemblages
	{
		cxxPPassemblage *entity_ptr = Utilities::Rxn_find(Rxn_pp_assemblage_map, n);
		if (entity_ptr != NULL)
		{
			sb.Set_PPassemblage(n, entity_ptr);
		}
	}
	// SSassemblages
	{
		cxxSSassemblage *entity_ptr = Utilities::Rxn_find(Rxn_ss_assemblage_map, n);
		if (entity_ptr != NULL)
		{
			sb.Set_SSassemblage(n, entity_ptr);
		}
	}
	// Surfaces
	{
		cxxSurface *entity_ptr = Utilities::Rxn_find(Rxn_surface_map, n);
		if (entity_ptr != NULL)
		{
			sb.Set_Surface(n, entity_ptr);
		}
	}
}
void Phreeqc::
cxxStorageBin2phreeqc(cxxStorageBin & sb, int n)
//
// copy all reactants from storage bin number n to phreeqc
// replaces any existing reactants in phreeqc
//
{
	// Solutions
	{
		std::map < int, cxxSolution >::const_iterator it = sb.Get_Solutions().find(n);
		if (it != sb.Get_Solutions().end())
		{
			Rxn_solution_map[n] = it->second;
		}
	}
	// Exchangers
	{
		std::map < int, cxxExchange >::const_iterator it = sb.Get_Exchangers().find(n);
		if (it != sb.Get_Exchangers().end())
		{
			Rxn_exchange_map[n] = it->second;
		}
	}

	// GasPhases
	{
		std::map < int, cxxGasPhase >::const_iterator it = sb.Get_GasPhases().find(n);
		if (it != sb.Get_GasPhases().end())
		{
			Rxn_gas_phase_map[n] = it->second;
		}
	}

	// Kinetics
	{
		std::map < int, cxxKinetics >::const_iterator it = sb.Get_Kinetics().find(n);
		if (it != sb.Get_Kinetics().end())
		{
			Rxn_kinetics_map[n] = it->second;
		}
	}
	// PPassemblages
	{
		std::map < int, cxxPPassemblage >::const_iterator it = sb.Get_PPassemblages().find(n);
		if (it != sb.Get_PPassemblages().end())
		{
			Rxn_pp_assemblage_map[n] = it->second;
		}
	}
	// SSassemblages
	{
		std::map < int, cxxSSassemblage >::const_iterator it = sb.Get_SSassemblages().find(n);
		if (it != sb.Get_SSassemblages().end())
		{
			Rxn_ss_assemblage_map[n] = it->second;
		}
	}
	// Surfaces
	{
		std::map < int, cxxSurface >::const_iterator it = sb.Get_Surfaces().find(n);
		if (it != sb.Get_Surfaces().end())
		{
			Rxn_surface_map[n] = it->second;
		}
	}
	// Mixes
	{
		std::map < int, cxxMix >::const_iterator it = sb.Get_Mixes().find(n);
		if (it != sb.Get_Mixes().end())
		{
			Rxn_mix_map[n] = it->second;
		}
	}

	// Reactions
	{
		std::map < int, cxxReaction >::const_iterator it = sb.Get_Reactions().find(n);
		if (it != sb.Get_Reactions().end())
		{
			Rxn_reaction_map[n] = it->second;
		}
	}
	// Temperatures
	{
		std::map < int, cxxTemperature >::const_iterator it = sb.Get_Temperatures().find(n);
		if (it != sb.Get_Temperatures().end())
		{
			Rxn_temperature_map[n] = it->second;
		}
	}
	// Pressures
	{
		std::map < int, cxxPressure >::const_iterator it = sb.Get_Pressures().find(n);
		if (it != sb.Get_Pressures().end())
		{
			Rxn_pressure_map[n] = it->second;
		}
	}
}
void Phreeqc::
cxxStorageBin2phreeqc(cxxStorageBin & sb)
//
// copy data from storage bin to phreeqc
// replaces any existing reactants in phreeqc
//
{
	// Solutions
	{
		std::map < int, cxxSolution >::const_iterator it = sb.Get_Solutions().begin();
		for ( ; it != sb.Get_Solutions().end(); it++)
		{
			Rxn_solution_map[it->first] = it->second;
		}
	}
	// Exchangers
	{
		std::map < int, cxxExchange >::const_iterator it = sb.Get_Exchangers().begin();
		for ( ; it != sb.Get_Exchangers().end(); it++)
		{
			Rxn_exchange_map[it->first] = it->second;
		}
	}

	// GasPhases
	{
		std::map < int, cxxGasPhase >::const_iterator it = sb.Get_GasPhases().begin();
		for ( ; it != sb.Get_GasPhases().end(); it++)
		{
			Rxn_gas_phase_map[it->first] = it->second;
		}
	}

	// Kinetics
	{
		std::map < int, cxxKinetics >::const_iterator it = sb.Get_Kinetics().begin();
		for ( ; it != sb.Get_Kinetics().end(); it++)
		{
			Rxn_kinetics_map[it->first] = it->second;
		}
	}
	// PPassemblages
	{
		std::map < int, cxxPPassemblage >::const_iterator it = sb.Get_PPassemblages().begin();
		for ( ; it != sb.Get_PPassemblages().end(); it++)
		{
			Rxn_pp_assemblage_map[it->first] = it->second;
		}
	}
	// SSassemblages
	{
		std::map < int, cxxSSassemblage >::const_iterator it = sb.Get_SSassemblages().begin();
		for ( ; it != sb.Get_SSassemblages().end(); it++)
		{
			Rxn_ss_assemblage_map[it->first] = it->second;
		}
	}
	// Surfaces
	{
		std::map < int, cxxSurface >::const_iterator it = sb.Get_Surfaces().begin();
		for ( ; it != sb.Get_Surfaces().end(); it++)
		{
			Rxn_surface_map[it->first] = it->second;
		}
	}
	// Mixes
	{
		std::map < int, cxxMix >::const_iterator it = sb.Get_Mixes().begin();
		for ( ; it != sb.Get_Mixes().end(); it++)
		{
			Rxn_mix_map[it->first] = it->second;
		}
	}

	// Reactions
	{
		std::map < int, cxxReaction >::const_iterator it = sb.Get_Reactions().begin();
		for ( ; it != sb.Get_Reactions().end(); it++)
		{
			Rxn_reaction_map[it->first] = it->second;
		}
	}
	// Temperatures
	{
		std::map < int, cxxTemperature >::const_iterator it = sb.Get_Temperatures().begin();
		for ( ; it != sb.Get_Temperatures().end(); it++)
		{
			Rxn_temperature_map[it->first] = it->second;
		}
	}
	// Pressures
	{
		std::map < int, cxxPressure >::const_iterator it = sb.Get_Pressures().begin();
		for ( ; it != sb.Get_Pressures().end(); it++)
		{
			Rxn_pressure_map[it->first] = it->second;
		}
	}
}


