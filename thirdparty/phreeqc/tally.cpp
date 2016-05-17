#include "Utils.h"
#include "Phreeqc.h"
#include "phqalloc.h"
#include "Temperature.h"
#include "Exchange.h"
#include "GasPhase.h"
#include "Reaction.h"
#include "PPassemblage.h"
#include "SSassemblage.h"
#include "cxxKinetics.h"
#include "Solution.h"
/*   
     Calling sequence 

Initialization:
---------------

build_tally_table();                   finds all elements (rows),
                                       possible return values (columns),
				       allocates space
get_tally_table_rows_columns(int *rows, int *columns)
				       returns number of rows and columns in table
get_tally_table_row_heading(int row, char *string)
				       row is C row number
				       returns row descripter for row
get_tally_table_column_heading(int column, int *type, char *string)
				       column is C column number
				       returns column heading for column

Each call to phreeqc:
---------------------

zero_tally_table();                    initialize table to 0s
set_reaction_moles(n_user, moles)      n_user is reservoir number
                                       moles is number of moles of reaction to add

int set_reaction_temperature(n_user, tc)

fill_tally_table(int *n_user, int index_conservative, int n_buffer)
                                       n_user is reservoir number
				       index_conservative is solution number
				           where conservative mixing is stored
                                       slot is 0 for initial

   run phreeqc here.

fill_tally_table(int *n_user, int index_conservative, int n_buffer)
                                       n_user is reservoir number
				       index_conservative is solution number
				           where conservative mixing is stored
                                       slot is 1 for final
store_tally_table(LDBLE *array, int row_dim, int col_dim, LDBLE fill_factor) 
                                       row_dim is Fortran dimension
                                       col_dim is Fortran dimension
				       array is space from Fortran
				       stores conservative mixing (column 0)
				       stores reaction (column 1)
				       difference between slot 1 and slot 0 for
				       all other entities (columns 2-n)

Finalization:
-------------
int free_tally_table(void);       Frees space

*/
/* ---------------------------------------------------------------------- */
int Phreeqc::
get_all_components(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Counts components in any defined solution, gas_phase, exchanger,
 *   surface, or pure_phase_assemblage
 *
 *   Returns n_comp, which is total, including H, O, elements, and Charge
 *           names contains character strings with names of components
 */
	int i, j;
/*
 *   Accumulate all aqueous components
 */
	add_all_components_tally();

	// add secondary master species
	for (i = 0; i < count_master; i++)
	{
		if (master[i]->total > 0.0 && master[i]->s->type == AQ && master[i]->primary == TRUE)
		{
			for (int j = i + 1; j < count_master; j++)
			{
				if (master[j]->elt->primary == master[i])
				{
					master[j]->total = 1.0;
				}
				else
				{
					break;
				}
			}
		}
	}


/*
 *   Count components + Alkalinity + total_h + total_o
 */
	tally_count_component = 3;
	for (i = 0; i < count_master; i++)
	{
		if (master[i]->total > 0.0 && master[i]->s->type == AQ)
		{
			tally_count_component++;
		}
	}
/*
 *   Put information in buffer.
 *   Buffer contains an entry for every primary master
 *   species that can be used in the transport problem.
 */
	t_buffer =
		(struct tally_buffer *) PHRQ_malloc((size_t) tally_count_component *
											sizeof(struct tally_buffer));

	// store alkalinity
	j = 0;
	t_buffer[j].name = string_hsave("Alkalinity");
	t_buffer[j].master = master_bsearch("Alkalinity");
	t_buffer[j].gfw = t_buffer[j].master->elt->gfw;
	j++;		

	// store total_h
	t_buffer[j].name = string_hsave("Total_H");
	t_buffer[j].master = NULL;
	compute_gfw("H", &(t_buffer[j].gfw));
	j++;

	// store total_o
	t_buffer[j].name = string_hsave("Total_O");
	t_buffer[j].master = NULL;
	compute_gfw("O", &(t_buffer[j].gfw));
	j++;

	for (i = 0; i < count_master; i++)
	{
		if (master[i]->total > 0.0 && master[i]->s->type == AQ)
		{
			t_buffer[j].name = master[i]->elt->name;
			t_buffer[j].master = master[i];
			t_buffer[j].gfw = master[i]->elt->gfw;
			j++;
		}
	}
	/*
	 *  Return value
	 */
	/**n_comp = count_component;*/
	count_tally_table_rows = tally_count_component;
	return (OK);
}
#ifdef SKIP
/* ---------------------------------------------------------------------- */
int Phreeqc::
get_all_components(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Counts components in any defined solution, gas_phase, exchanger,
 *   surface, or pure_phase_assemblage
 *
 *   Returns n_comp, which is total, including H, O, elements, and Charge
 *           names contains character strings with names of components
 */
	int i, j;
/*
 *   Accumulate all aqueous components
 */
	add_all_components_tally();
/*
 *   Count components, 2 for hydrogen, oxygen,  + others,
 */
	tally_count_component = 0;
	for (i = 0; i < count_master; i++)
	{
		if (master[i]->total > 0.0 && master[i]->s->type == AQ)
		{
			tally_count_component++;
		}
	}
/*
 *   Put information in buffer.
 *   Buffer contains an entry for every primary master
 *   species that can be used in the transport problem.
 *   Each entry in buffer is sent to HST for transort.
 */
	t_buffer =
		(struct tally_buffer *) PHRQ_malloc((size_t) tally_count_component *
											sizeof(struct tally_buffer));
	j = 0;
	for (i = 0; i < count_master; i++)
	{
		if (master[i]->total > 0.0 && master[i]->s->type == AQ)
		{
			t_buffer[j].name = master[i]->elt->name;
			t_buffer[j].master = master[i];
			t_buffer[j].gfw = master[i]->elt->gfw;
			j++;
		}
	}
	/*
	 *  Return value
	 */
	/**n_comp = count_component;*/
	count_tally_table_rows = tally_count_component;
	return (OK);
}
#endif
/* ---------------------------------------------------------------------- */
int Phreeqc::
store_tally_table(LDBLE * l_array, int row_dim_in, int col_dim, LDBLE fill_factor)
/* ---------------------------------------------------------------------- */
{
	int i, j;
	int row_dim = row_dim_in + 1;
	if (tally_table == NULL)
	{
		input_error++;
		error_msg("Tally table not defined, get_tally_table_rows_columns",
				  CONTINUE);
		return (ERROR);
	}
	if (count_tally_table_rows > row_dim)
	{
		input_error++;
		error_msg
			("Too many tally table rows for Fortran storage, store_tally_table",
			 CONTINUE);
		return (ERROR);
	}
	if (count_tally_table_columns > col_dim)
	{
		input_error++;
		error_msg
			("Too many tally table columns for Fortran storage, store_tally_table",
			 CONTINUE);
		return (ERROR);
	}
	/*
	 * store conservative mixing solution
	 */
	for (j = 0; j < count_tally_table_rows; j++)
	{
		l_array[j] = tally_table[0].total[1][j].moles;
	}
	/*
	 * store reaction solution
	 */
	for (j = 0; j < count_tally_table_rows; j++)
	{
		l_array[row_dim + j] = tally_table[1].total[1][j].moles;
	}
	/*
	 *   Calculate deltas
	 */

	diff_tally_table();

	/*
	 * store deltas for all other columns
	 */
	for (i = 2; i < count_tally_table_columns; i++)
	{
		for (j = 0; j < count_tally_table_rows; j++)
		{
			l_array[i * row_dim + j] =
				tally_table[i].total[2][j].moles / fill_factor;
		}
	}

	/*
	 * Add row for total moles of reactant
	 */
	for (i = 0; i < count_tally_table_columns; i++)
	{
		l_array[i * row_dim + count_tally_table_rows] =
				tally_table[i].moles / fill_factor;
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
get_tally_table_rows_columns(int *rows, int *columns)
/* ---------------------------------------------------------------------- */
{
	*rows = 0;
	*columns = 0;
	if (tally_table == NULL)
	{
		input_error++;
		error_msg("tally table not defined, get_tally_table_rows_columns",
				  CONTINUE);
		return (ERROR);
	}
	*rows = count_tally_table_rows;
	*columns = count_tally_table_columns;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
get_tally_table_row_heading(int row, char *string)
/* ---------------------------------------------------------------------- */
{
	/*
	 *  row is C row number
	 */
	strcpy(string, "");
	if (tally_table == NULL)
	{
		input_error++;
		error_msg("Tally table not defined, get_tally_table row_heading",
				  CONTINUE);
		return (ERROR);
	}
	if (row >= count_tally_table_rows)
	{
		input_error++;
		error_msg
			("Row exceeds tally table size, get_tally_table row_heading",
			 CONTINUE);
		return (ERROR);
	}
	strcpy(string, t_buffer[row].name);
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
get_tally_table_column_heading(int column, int *type, char *string)
/* ---------------------------------------------------------------------- */
{
	/*
	 *  column is C column number
	 */
	*type = -1;
	strcpy(string, "");
	if (tally_table == NULL)
	{
		input_error++;
		error_msg("tally table not defined, get_tally_table_column_heading",
				  CONTINUE);
		return (ERROR);
	}
	if (column >= count_tally_table_columns)
	{
		input_error++;
		error_msg
			("column exceeds tally table size, get_tally_table_column_heading",
			 CONTINUE);
		return (ERROR);
	}
	strcpy(string, tally_table[column].name);
	*type = tally_table[column].type;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
free_tally_table(void)
/* ---------------------------------------------------------------------- */
{
	int i, k;
	if (tally_table == NULL)
		return (OK);
	for (i = 0; i < count_tally_table_columns; i++)
	{
		if (tally_table[i].formula != NULL)
			tally_table[i].formula =
				(struct elt_list *) free_check_null(tally_table[i].formula);
		for (k = 0; k < 3; k++)
		{
			tally_table[i].total[k] =
				(struct tally_buffer *) free_check_null(tally_table[i].
														total[k]);
		}
	}
	tally_table = (struct tally *) free_check_null(tally_table);
	t_buffer = (struct tally_buffer *) free_check_null(t_buffer);
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
zero_tally_table(void)
/* ---------------------------------------------------------------------- */
{
	int i, j, k;
	for (i = 0; i < count_tally_table_columns; i++)
	{
		for (j = 0; j < count_tally_table_rows; j++)
		{
			for (k = 0; k < 3; k++)
			{
				tally_table[i].total[k][j].moles = 0;
			}
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
diff_tally_table(void)
/* ---------------------------------------------------------------------- */
{
	int i, j;
	/*
	   output_msg("Difference\n\n");
	 */
	for (i = 0; i < count_tally_table_columns; i++)
	{
		for (j = 0; j < count_tally_table_rows; j++)
		{
			tally_table[i].total[2][j].moles =
				tally_table[i].total[1][j].moles -
				tally_table[i].total[0][j].moles;
		}

		/*
		   output_msg(sformatf( "Column %d\t%s\tType: %d\n", i, tally_table[i].name, tally_table[i].type));
		   for (j = 0; j < count_tally_table_rows; j++) {
		   output_msg(sformatf( "\t%d\t%s\t%e\n", j, tally_table[i].total[2][j].name, (double) tally_table[i].total[2][j].moles));
		   }
		 */
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
print_tally_table(void)
/* ---------------------------------------------------------------------- */
{
	int i, j;
	output_msg(sformatf( "Tally_table\n\n"));
	for (i = 0; i < count_tally_table_columns; i++)
	{
		output_msg(sformatf( "%s\tType: %d\n", tally_table[i].name,
				   tally_table[i].type));
		output_msg(sformatf( "\n"));
		output_msg(sformatf( "\t%15s\t%15s\t%15s\n", "Initial",
				   "Final", "Difference"));
		for (j = 0; j < count_tally_table_rows; j++)
		{
			output_msg(sformatf( "%5s\t%15g\t%15g\t%15g\n",
					   t_buffer[j].name, tally_table[i].total[0][j].moles,
					   tally_table[i].total[1][j].moles,
					   tally_table[i].total[2][j].moles));
		}
		output_msg(sformatf( "\n"));
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
fill_tally_table(int *n_user, int index_conservative, int n_buffer)
/* ---------------------------------------------------------------------- */
{
/*
 *   Routine accumulates elements from all solutions, phases, gas phases,
 *   exchangers, and surfaces. 
 */
	int found;
	LDBLE moles;
	//char *ptr;
	/*
	 *  Cycle through tally table columns
	 */
	for (int i = 0; i < count_tally_table_columns; i++)
	{
		switch (tally_table[i].type)
		{
		case Solution:
/*
 *   fill solution
 */
			if (n_user[Solution] < 0 || n_buffer == 0)
				break;
			{
				cxxSolution *solution_ptr = NULL;;
				if (i == 0)
				{
					solution_ptr = Utilities::Rxn_find(Rxn_solution_map, index_conservative);
				}
				else if (i == 1)
				{
					solution_ptr = Utilities::Rxn_find(Rxn_solution_map, n_user[Solution]);
				}
				else
				{
					error_msg("Solution is not in first two columns of tally_table", STOP);
				}
				if (solution_ptr == NULL)
					break;
				/*
				 *   Add secondary master species
				 */


				xsolution_zero();

				// adds primary master species
				add_solution(solution_ptr, 1.0, 1.0);

				// adds secondary master species
				cxxNameDouble::iterator jit = solution_ptr->Get_totals().begin();
				for ( ; jit != solution_ptr->Get_totals().end(); jit++)
				{
					struct master *master_ptr = master_bsearch(jit->first.c_str());
					master_ptr->total = jit->second;
				}

				// Fill table
				master_to_tally_table(tally_table[i].total[n_buffer]);
				
				// Add alkalinity
				tally_table[i].total[n_buffer][0].moles = solution_ptr->Get_total_alkalinity();

				// Add total_h
				tally_table[i].total[n_buffer][1].moles = solution_ptr->Get_total_h();

				// Add total_o
				tally_table[i].total[n_buffer][2].moles = solution_ptr->Get_total_o();
			}
			break;
#ifdef SKIP
		case Solution:
/*
 *   fill solution
 */
			if (n_user[Solution] < 0 || n_buffer == 0)
				break;
			{
				cxxSolution *solution_ptr;
				if (i == 0)
				{
					solution_ptr = Utilities::Rxn_find(Rxn_solution_map, index_conservative);
				}
				else if (i == 1)
				{
					solution_ptr = Utilities::Rxn_find(Rxn_solution_map, n_user[Solution]);
				}
				else
				{
					solution_ptr = NULL;
					error_msg
						("Solution is not in first two columns of tally_table",
						STOP);
				}
				if (solution_ptr == NULL)
					break;
				xsolution_zero();
				add_solution(solution_ptr, 1.0, 1.0);
				count_elts = 0;
				paren_count = 0;
				for (int j = 0; j < count_master; j++)
				{
					if (master[j]->total > 0.0)
					{
						char * temp_name = string_duplicate(master[j]->elt->primary->elt->name);
						ptr = temp_name;
						get_elts_in_species(&ptr, master[j]->total);
						free_check_null(temp_name);
					}
				}
				qsort(elt_list, (size_t) count_elts,
					(size_t) sizeof(struct elt_list), elt_list_compare);
				elt_list_combine();
				elt_list_to_tally_table(tally_table[i].total[n_buffer]);
			}
			break;
#endif
		case Reaction:
			/*
			 *   fill reaction
			 */
			if (n_user[Reaction] < 0)
				break;
			{
				cxxReaction *reaction_ptr = Utilities::Rxn_find(Rxn_reaction_map, n_user[Reaction]);
				if (reaction_ptr == NULL)
					break;
				count_elts = 0;
				paren_count = 0;
				if (n_buffer == 1)
				{
					moles = reaction_ptr->Get_steps()[0];
				}
				else
				{
					moles = 0.0;
				}
				reaction_calc(reaction_ptr);
				add_elt_list(reaction_ptr->Get_elementList(), moles);
				elt_list_to_tally_table(tally_table[i].total[n_buffer]);
			}
			break;
		case Pure_phase:
			/*
			 *   fill an equilibrium phase
			 */
			if (n_user[Pure_phase] < 0)
				break;
			{
				cxxPPassemblage * pp_assemblage_ptr = Utilities::Rxn_find(Rxn_pp_assemblage_map, n_user[Pure_phase]);
				if (pp_assemblage_ptr == NULL)
					break;
				std::map<std::string, cxxPPassemblageComp>::iterator it;
				it =  pp_assemblage_ptr->Get_pp_assemblage_comps().begin();
				for ( ; it != pp_assemblage_ptr->Get_pp_assemblage_comps().end(); it++)
				{
					if (string_hsave(it->second.Get_name().c_str()) ==
						tally_table[i].name)
						break;
					if (strcmp_nocase(it->second.Get_name().c_str(),
						tally_table[i].name) == 0)
						break;
				}
				if (it == pp_assemblage_ptr->Get_pp_assemblage_comps().end())
					break;
				count_elts = 0;
				paren_count = 0;
				moles = it->second.Get_moles();
				tally_table[i].moles = moles;
				add_elt_list(tally_table[i].formula, moles);
				elt_list_to_tally_table(tally_table[i].total[n_buffer]);
			}
			break;
		case Exchange:
			{
				/*
				*   fill exchange
				*/
				if (n_user[Exchange] < 0)
					break;
				cxxExchange * exchange_ptr = Utilities::Rxn_find(Rxn_exchange_map, n_user[Exchange]);
				if (exchange_ptr == NULL)
					break;
				count_elts = 0;
				paren_count = 0;
				for (size_t j = 0; j < exchange_ptr->Get_exchange_comps().size(); j++)
				{
					add_elt_list(exchange_ptr->Get_exchange_comps()[j].Get_totals(), 1.0);
				}
				qsort(elt_list, (size_t) count_elts,
					(size_t) sizeof(struct elt_list), elt_list_compare);
				elt_list_combine();
				elt_list_to_tally_table(tally_table[i].total[n_buffer]);
			}
			break;
		case Surface:
			/*
			 *   fill surface
			 */
			if (n_user[Surface] < 0)
				break;
			{
				cxxSurface * surface_ptr = Utilities::Rxn_find(Rxn_surface_map, n_user[Surface]);
				if (surface_ptr == NULL)
					break;
				count_elts = 0;
				paren_count = 0;
				for (size_t j = 0; j < surface_ptr->Get_surface_comps().size(); j++)
				{
					add_elt_list(surface_ptr->Get_surface_comps()[j].Get_totals(), 1.0);
				}
				qsort(elt_list, (size_t) count_elts,
					(size_t) sizeof(struct elt_list), elt_list_compare);
				elt_list_combine();
				elt_list_to_tally_table(tally_table[i].total[n_buffer]);
			}
			break;
		case Ss_phase:
			if (n_user[Ss_phase] < 0)
				break;
			{
				/*
				*   fill an solid solution phase
			    */

				cxxSSassemblage * ss_assemblage_ptr = Utilities::Rxn_find(Rxn_ss_assemblage_map, n_user[Ss_phase]);
				if (ss_assemblage_ptr == NULL)
					break;
				found = FALSE;
				moles = 0.0;
				std::vector<cxxSS *> ss_ptrs = ss_assemblage_ptr->Vectorize();
				for (size_t j = 0; j < ss_ptrs.size(); j++)
				{
					cxxSS * ss_ptr = ss_ptrs[j];
					cxxSScomp * comp_ptr = NULL;
					size_t k;
					for (k = 0; k < ss_ptr->Get_ss_comps().size(); k++)
					{
						comp_ptr = &(ss_ptr->Get_ss_comps()[k]);
						int l;
						struct phase *phase_ptr = phase_bsearch(comp_ptr->Get_name().c_str(), &l, FALSE);
						if (phase_ptr->name == tally_table[i].name)
							break;
						if (strcmp_nocase(phase_ptr->name, tally_table[i].name) == 0)
							break;
					}
					if (k < ss_ptr->Get_ss_comps().size() && comp_ptr)
					{
						moles = comp_ptr->Get_moles();
						found = TRUE;
						break;
					}
				}
				if (found == FALSE)
					break;
				count_elts = 0;
				paren_count = 0;
				tally_table[i].moles = moles;
				add_elt_list(tally_table[i].formula, moles);
				elt_list_to_tally_table(tally_table[i].total[n_buffer]);
			}
			break;
		case Gas_phase:
			/*
			 *   fill in gas phase
			 */
			if (n_user[Gas_phase] < 0)
				break;
			{
				cxxGasPhase * gas_phase_ptr = Utilities::Rxn_find(Rxn_gas_phase_map, n_user[Gas_phase]);
				if (gas_phase_ptr == NULL)
					break;
				count_elts = 0;
				paren_count = 0;
				const std::vector<cxxGasComp> *gc = &(gas_phase_ptr->Get_gas_comps());
				for (size_t l = 0; l < gc->size(); l++)
				{
					int k;
					struct phase *phase_ptr = phase_bsearch((*gc)[l].Get_phase_name().c_str(), &k, FALSE);

					add_elt_list(phase_ptr->next_elt, (*gc)[l].Get_moles());
				}
				qsort(elt_list, (size_t) count_elts,
					(size_t) sizeof(struct elt_list), elt_list_compare);
				elt_list_combine();
				elt_list_to_tally_table(tally_table[i].total[n_buffer]);
				break;
			}
		case Kinetics:
			{
				/*
				*   fill in kinetics
				*/
				if (n_user[Kinetics] < 0)
					break;
				cxxKinetics *kinetics_ptr = Utilities::Rxn_find(Rxn_kinetics_map, n_user[Kinetics]);
				if (kinetics_ptr == NULL)
					break;
				cxxKineticsComp * kinetics_comp_ptr = NULL;
				size_t j;
				for (j = 0; j < kinetics_ptr->Get_kinetics_comps().size(); j++)
				{
					kinetics_comp_ptr = &(kinetics_ptr->Get_kinetics_comps()[j]);
					if (string_hsave(kinetics_comp_ptr->Get_rate_name().c_str()) == tally_table[i].name)
						break;
					if (strcmp_nocase
						(kinetics_comp_ptr->Get_rate_name().c_str(), tally_table[i].name) == 0)
						break;
				}
				if (j >= kinetics_ptr->Get_kinetics_comps().size())
					break;
				moles = 0.0;
				if (kinetics_comp_ptr)
				{
					moles = kinetics_comp_ptr->Get_m();
				}
				tally_table[i].moles = moles;
				count_elts = 0;
				paren_count = 0;
				add_elt_list(tally_table[i].formula, moles);
				elt_list_to_tally_table(tally_table[i].total[n_buffer]);
			}
			break;
		case Mix:
			break;
		case Temperature:
		case Pressure:
			break;
		case UnKnown:
			break;
		}
	}

	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
elt_list_to_tally_table(struct tally_buffer *buffer_ptr)
/* ---------------------------------------------------------------------- */
{
	int i, j;
	for (i = 0; i < count_tally_table_rows; i++)
	{
		buffer_ptr[i].moles = 0.0;
	}
	/*
	 * copy element list amounts to buffer in tally table
	 * for column number
	 */

	for (j = 0; j < count_elts; j++)
	{
		if (elt_list[j].elt->primary->s == s_h2o)
			continue;
		if (elt_list[j].elt->primary->s == s_hplus)
			continue;
		if (elt_list[j].elt->primary->s == s_h3oplus)
			continue;
		if (elt_list[j].elt->primary->type != AQ)
			continue;
		for (i = 0; i < count_tally_table_rows; i++)
		{
			if (buffer_ptr[i].master != NULL)
			{
				if (elt_list[j].elt->primary ==
					buffer_ptr[i].master->elt->primary)
				{
					buffer_ptr[i].moles = elt_list[j].coef;
					break;
				}
			}
		}
		if (i >= count_tally_table_rows)
		{
			error_msg("Should not be here in elt_list_to_tally_table", STOP);
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
master_to_tally_table(struct tally_buffer *buffer_ptr)
/* ---------------------------------------------------------------------- */
{
	int i, j;
	for (i = 0; i < count_tally_table_rows; i++)
	{
		buffer_ptr[i].moles = 0.0;
	}
	/*
	 * copy element list amounts to buffer in tally table
	 * for column number
	 */
	for (j  = 0; j < count_master; j++)
	{
		if (master[j]->total <= 0)
			continue;
		if (master[j]->elt->primary->s == s_h2o)
			continue;
		if (master[j]->elt->primary->s == s_hplus)
			continue;
		if (master[j]->elt->primary->s == s_h3oplus)
			continue;
		if (master[j]->elt->primary->type != AQ)
			continue;
		for (i = 0; i < count_tally_table_rows; i++)
		{
			if (master[j] ==  buffer_ptr[i].master)
			{
				buffer_ptr[i].moles = master[j]->total;
				break;
			}
		}
		if (i >= count_tally_table_rows)
		{
			error_msg("Should not be here in master_to_tally_table", STOP);
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
build_tally_table(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Routine accumulates elements from all solutions, phases, gas phases,
 *   exchangers, and surfaces. Counts number of aqueous components
 *   to transport. Stores in global variable tally_count_component.
 *   Also calculates a number greater than all user numbers and
 *   stores in global variable first_user_number.
 */
	int j, k, l, n, p, save_print_use;
	int count_tt_pure_phase, count_tt_ss_phase, count_tt_kinetics;
	struct phase *phase_ptr;
	char token[MAX_LENGTH];
	char *ptr;
/*
 *  make list of all elements in all entitites
 *  defines the number of rows in the table
 */
	get_all_components();

	save_print_use = pr.use;
	pr.use = FALSE;
/*
 *  find nuber of columns
 */
	count_tally_table_columns = 0;
/*
 *   add one for conservative mixing
 */
	n = count_tally_table_columns;
	extend_tally_table();
	tally_table[n].name = string_hsave("Solution_conservative");
	tally_table[n].type = Solution;
/*
 *   add one for mixing plus reaction
 */
	n = count_tally_table_columns;
	extend_tally_table();
	tally_table[n].name = string_hsave("Solution_reaction");
	tally_table[n].type = Solution;
/*
 *   add one for reactions
 */
	if (Rxn_reaction_map.size() > 0)
	{
		n = count_tally_table_columns;
		extend_tally_table();
		tally_table[n].name = string_hsave("Reaction");
		tally_table[n].type = Reaction;
	}
/*
 *   add one for exchangers
 */
	if (Rxn_exchange_map.size() > 0)
	{
		n = count_tally_table_columns;
		extend_tally_table();
		tally_table[n].name = string_hsave("Exchange");
		tally_table[n].type = Exchange;
	}
/*
 *   add one for surface
 */
	if (Rxn_surface_map.size() > 0)
	{
		n = count_tally_table_columns;
		extend_tally_table();
		tally_table[n].name = string_hsave("Surface");
		tally_table[n].type = Surface;
	}
/*
 *   add one for gases
 */
	if (Rxn_gas_phase_map.size() > 0)
	{
		n = count_tally_table_columns;
		extend_tally_table();
		tally_table[n].name = string_hsave("Gas_phase");
		tally_table[n].type = Gas_phase;
	}
/*
 *   Count pure phases
 */
	count_tt_pure_phase = 0;
	if (Rxn_pp_assemblage_map.size() > 0)
	{
		/* 
		 * Go through all pure phases in pure phase assemblages
		 */
		std::map<int, cxxPPassemblage>::iterator it;
		for (it = Rxn_pp_assemblage_map.begin(); it != Rxn_pp_assemblage_map.end(); it++)
		{
			cxxPPassemblage * pp_assemblage_ptr = &(it->second);
			std::map<std::string, cxxPPassemblageComp>::iterator jit;
			jit =  pp_assemblage_ptr->Get_pp_assemblage_comps().begin();
			for ( ; jit != pp_assemblage_ptr->Get_pp_assemblage_comps().end(); jit++)
			{
				cxxPPassemblageComp * comp_ptr = &(jit->second);
				int l;
				struct phase * phase_ptr = phase_bsearch(jit->first.c_str(), &l, FALSE);
				/* 
				 * check if already in tally_table
				 */
				for (k = 1; k < count_tally_table_columns; k++)
				{
					if (tally_table[k].type == Pure_phase &&
						tally_table[k].name == phase_ptr->name &&
						tally_table[k].add_formula ==
						string_hsave(comp_ptr->Get_add_formula().c_str()))
						break;
				}
				if (k < count_tally_table_columns)
					continue;
				/*
				 * Add to table
				 */
				count_tt_pure_phase++;
				n = count_tally_table_columns;
				extend_tally_table();
				tally_table[n].name = phase_ptr->name;
				tally_table[n].type = Pure_phase;
				tally_table[n].add_formula = string_hsave(comp_ptr->Get_add_formula().c_str());
				count_elts = 0;
				paren_count = 0;
				if (comp_ptr->Get_add_formula().size() > 0)
				{
					strcpy(token, comp_ptr->Get_add_formula().c_str());
					ptr = &(token[0]);
					get_elts_in_species(&ptr, 1.0);
				}
				else
				{
					strcpy(token, phase_ptr->formula);
					add_elt_list(phase_ptr->next_elt, 1.0);
				}
				qsort(elt_list, (size_t) count_elts,
					  (size_t) sizeof(struct elt_list), elt_list_compare);
				elt_list_combine();
				tally_table[n].formula = elt_list_save();
			}
		}
	}
/*
 *   Add solid-solution pure phases
 */
	count_tt_ss_phase = 0;
	if (Rxn_ss_assemblage_map.size() > 0)
	{
		/* 
		 * Go through all components of all solid solutions in solid-solution assemblages
		 */
		std::map<int, cxxSSassemblage>::iterator it;
		for (it = Rxn_ss_assemblage_map.begin(); it != Rxn_ss_assemblage_map.end(); it++)
		{
			cxxSSassemblage *ss_assemblage_ptr = &(it->second);
			std::vector<cxxSS *> ss_ptrs = ss_assemblage_ptr->Vectorize();
			for (j = 0; j < (int) ss_ptrs.size(); j++)
			{
				cxxSS * ss_ptr = ss_ptrs[j];
				for (k = 0; k < (int) ss_ptr->Get_ss_comps().size(); k++)
				{
					cxxSScomp *comp_ptr = &(ss_ptr->Get_ss_comps()[k]);
					int l;
					struct phase *phase_ptr = phase_bsearch(comp_ptr->Get_name().c_str(), &l, FALSE);
					/* 
					 * check if already in tally_table
					 */
					for (l = 1; l < count_tally_table_columns; l++)
					{
						if (tally_table[l].type == Ss_phase &&
							tally_table[l].name == phase_ptr->name)
							break;
					}
					if (l < count_tally_table_columns)
						continue;
					/*
					 * Add to table
					 */
					count_tt_ss_phase++;
					n = count_tally_table_columns;
					extend_tally_table();
					tally_table[n].name = phase_ptr->name;
					tally_table[n].type = Ss_phase;
					count_elts = 0;
					paren_count = 0;
					strcpy(token, phase_ptr->formula);
					add_elt_list(phase_ptr->next_elt, 1.0);
					qsort(elt_list, (size_t) count_elts,
						  (size_t) sizeof(struct elt_list), elt_list_compare);
					elt_list_combine();
					tally_table[n].formula = elt_list_save();
				}
			}
		}
	}
/*
 *   Add kinetic reactants
 */
	count_tt_kinetics = 0;
	if (Rxn_kinetics_map.size() > 0)
	{
		std::map<int, cxxKinetics>::iterator it;
		for (it = Rxn_kinetics_map.begin(); it != Rxn_kinetics_map.end(); it++)
		{
			cxxKinetics *kinetics_ptr = &(it->second);
			for (j = 0; j < (int) kinetics_ptr->Get_kinetics_comps().size(); j++)
			{
				cxxKineticsComp *kinetics_comp_ptr = &(kinetics_ptr->Get_kinetics_comps()[j]);
				/* 
				 * check if already in tally_table
				 */
				for (l = 1; l < count_tally_table_columns; l++)
				{
					if (tally_table[l].type == Kinetics &&
						tally_table[l].name == string_hsave(kinetics_comp_ptr->Get_rate_name().c_str()))
						break;
				}
				if (l < count_tally_table_columns)
					continue;
				/*
				 * Add to table
				 */
				count_tt_kinetics++;
				n = count_tally_table_columns;
				extend_tally_table();
				tally_table[n].name = string_hsave(kinetics_comp_ptr->Get_rate_name().c_str());
				tally_table[n].type = Kinetics;
				/*
				 * get formula for kinetic component
				 */
				count_elts = 0;
				paren_count = 0;
				phase_ptr = NULL;
				if (kinetics_comp_ptr->Get_namecoef().size() == 1)
				{
					strcpy(token, kinetics_comp_ptr->Get_namecoef().begin()->first.c_str());
					phase_ptr = phase_bsearch(token, &p, FALSE);
				}
				if (phase_ptr != NULL)
				{
					add_elt_list(phase_ptr->next_elt, 1.0);
				}
				else
				{
					cxxNameDouble::iterator it = kinetics_comp_ptr->Get_namecoef().begin();
					for ( ; it != kinetics_comp_ptr->Get_namecoef().end(); it++)
					{
						std::string name = it->first;
						LDBLE coef = it->second;
						char * temp_name = string_duplicate(name.c_str());
						ptr = temp_name;
						get_elts_in_species(&ptr, 1.0 * coef);
						free_check_null(temp_name);
					}
				}
				qsort(elt_list, (size_t) count_elts,
					  (size_t) sizeof(struct elt_list), elt_list_compare);
				elt_list_combine();
				tally_table[n].formula = elt_list_save();
			}
		}
	}

#ifdef SKIP
	/*
	 *  Debug print for table definition
	 */
	output_msg(sformatf( "List of rows for tally table\n"));
	for (i = 0; i < count_tally_table_rows; i++)
	{
		output_msg(sformatf( "\t%-s\n", buffer[i].name));
	}
	output_msg(sformatf( "\nList of columns for tally table\n"));
	for (i = 0; i < count_tally_table_columns; i++)
	{
		output_msg(sformatf( "\t%-20s\tType: %d\n",
				   tally_table[i].name, tally_table[i].type));
		if (tally_table[i].formula != NULL)
		{
			for (j = 0; tally_table[i].formula[j].elt != NULL; j++)
			{
				output_msg(sformatf( "\t\t%-10s\t%f\n",
						   tally_table[i].formula[j].elt->name,
						   (double) tally_table[i].formula[j].coef));
			}
		}
	}
#endif
	pr.use = save_print_use;
	return (OK);
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
add_all_components_tally(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Routine accumulates elements from all solutions, phases, gas phases,
 *   exchangers, and surfaces. Counts number of aqueous components
 *   to transport. Stores in global variable tally_count_component.
 *   Also calculates a number greater than all user numbers and
 *   stores in global variable first_user_number.
 */
	int save_print_use;

	save_print_use = pr.use;
	pr.use = FALSE;
/*
 *   Delete solutions less than -1
 */

/*
 *   add all solutions
 */
	xsolution_zero();
	{
		std::map<int, cxxSolution>::iterator it = Rxn_solution_map.begin();
		for (; it != Rxn_solution_map.end(); it++)
		{
			add_solution(&it->second, 1.0 / it->second.Get_mass_water(), 1.0);
		}
	}
/*
 *   add all irrev reactions
 */
	{
		std::map<int, cxxReaction>::iterator it = Rxn_reaction_map.begin();
		for (; it != Rxn_reaction_map.end(); it++)
		{
			add_reaction(&it->second, 1, 1.0);
		}
	}
/*
 *   Add pure phases
 */
	{
		std::map<int, cxxPPassemblage>::iterator it = Rxn_pp_assemblage_map.begin();
		for (; it != Rxn_pp_assemblage_map.end(); it++)
		{
			add_pp_assemblage(&(it->second));
		}
	}
/*
 *   Exchangers
 */
	{
		std::map<int, cxxExchange>::iterator it = Rxn_exchange_map.begin();
		for (; it != Rxn_exchange_map.end(); it++)
		{
			add_exchange(&it->second);
		}
	}
/*
 *   Surfaces
 */
	{
		std::map<int, cxxSurface>::iterator it = Rxn_surface_map.begin();
		for (; it != Rxn_surface_map.end(); it++)
		{
			add_surface(&it->second);
		}
	}
/*
 *   Gases
 */
	{
		std::map<int, cxxGasPhase>::iterator it = Rxn_gas_phase_map.begin();
		for ( ; it != Rxn_gas_phase_map.end(); it++)
		{
			add_gas_phase(&it->second);
		}
	}
/*
 *   Add solid-solution pure phases
 */
	{
		std::map<int, cxxSSassemblage>::iterator it;
		for (it = Rxn_ss_assemblage_map.begin(); it != Rxn_ss_assemblage_map.end(); it++)
		{
			add_ss_assemblage(&(it->second));
		}
	}
/*
 *   Add elements in kinetic reactions
 */
	{
		std::map<int, cxxKinetics>::iterator it = Rxn_kinetics_map.begin();
		for ( ; it != Rxn_kinetics_map.end(); it++)
		{
			calc_dummy_kinetic_reaction_tally(&(it->second));
			add_kinetics(&(it->second));
		}
	}
/*
 *   reset pr.use
 */
	pr.use = save_print_use;
	return;
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
calc_dummy_kinetic_reaction_tally(cxxKinetics *kinetics_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *    Go through kinetic components and add positive amount of each reactant
 */
	LDBLE coef;
	char *ptr;
	struct phase *phase_ptr;
/*
 *   Go through list and generate list of elements and
 *   coefficient of elements in reaction
 */
	kinetics_ptr->Get_totals().clear();
	count_elts = 0;
	paren_count = 0;
	for (size_t i = 0; i < kinetics_ptr->Get_kinetics_comps().size(); i++)
	{
		cxxKineticsComp * kinetics_comp_ptr = &(kinetics_ptr->Get_kinetics_comps()[i]);
		coef = 1.0;
/*
 *   Reactant is a pure phase, copy formula into token
 */
		phase_ptr = NULL;
		if (kinetics_comp_ptr->Get_namecoef().size() == 1)
		{
			std::string name = kinetics_comp_ptr->Get_namecoef().begin()->first;
			int j;
			phase_ptr = phase_bsearch(name.c_str(), &j, FALSE);
		}
		if (phase_ptr != NULL)
		{
			add_elt_list(phase_ptr->next_elt, coef);
		}
		else
		{
			cxxNameDouble::iterator it = kinetics_comp_ptr->Get_namecoef().begin();
			for ( ; it != kinetics_comp_ptr->Get_namecoef().end(); it++)
			{
				std::string name = it->first;
				char * temp_name = string_duplicate(name.c_str());
				ptr = temp_name;
				get_elts_in_species(&ptr, coef);
				free_check_null(temp_name);
			}
		}
	}
	kinetics_ptr->Set_totals(elt_list_NameDouble());

	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
extend_tally_table(void)
/* ---------------------------------------------------------------------- */
{
	int i, j;
	/* 
	 * adds another column to tally_table
	 * increments number of columns
	 */
	tally_table =
		(struct tally *) PHRQ_realloc((void *) tally_table,
									  (size_t) (count_tally_table_columns +
												1) * sizeof(struct tally));
	if (tally_table == NULL)
		malloc_error();
	for (i = 0; i < 3; i++)
	{
		tally_table[count_tally_table_columns].total[i] =
			(struct tally_buffer *)
			PHRQ_malloc((size_t) (count_tally_table_rows) *
						sizeof(struct tally_buffer));
		if (tally_table[count_tally_table_columns].total[i] == NULL)
			malloc_error();
		for (j = 0; j < count_tally_table_rows; j++)
		{
			tally_table[count_tally_table_columns].total[i][j].name =
				t_buffer[j].name;
			tally_table[count_tally_table_columns].total[i][j].master =
				t_buffer[j].master;
		}
	}
	tally_table[count_tally_table_columns].name = NULL;
	tally_table[count_tally_table_columns].type = UnKnown;
	tally_table[count_tally_table_columns].add_formula = NULL;
	tally_table[count_tally_table_columns].moles = 0.0;
	tally_table[count_tally_table_columns].formula = NULL;
	count_tally_table_columns++;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
set_reaction_moles(int n_user, LDBLE moles)
/* ---------------------------------------------------------------------- */
{
	cxxReaction *reaction_ptr = Utilities::Rxn_find(Rxn_reaction_map, n_user);
	if (reaction_ptr == NULL)
		return (ERROR);
	std::vector<LDBLE> v;
	v.push_back(moles);
	reaction_ptr->Set_steps(v);
	reaction_ptr->Set_countSteps(1);
	reaction_ptr->Set_equalIncrements(true);
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
set_reaction_temperature(int n_user, LDBLE tc)
/* ---------------------------------------------------------------------- */
{
	cxxTemperature * temperature_ptr = Utilities::Rxn_find(Rxn_temperature_map, n_user);
	if (temperature_ptr == NULL)
		return (ERROR);
	temperature_ptr->Get_temps().clear();
	temperature_ptr->Get_temps().push_back(tc);
	temperature_ptr->Set_equalIncrements(false);
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
set_kinetics_time(int n_user, LDBLE step)
/* ---------------------------------------------------------------------- */
{
	cxxKinetics *kinetics_ptr = Utilities::Rxn_find(Rxn_kinetics_map, n_user);

	if (kinetics_ptr == NULL)
		return (ERROR);
	kinetics_ptr->Get_steps().clear();
	kinetics_ptr->Get_steps().push_back(step);
	kinetics_ptr->Set_equalIncrements(false);
	return (OK);
}

