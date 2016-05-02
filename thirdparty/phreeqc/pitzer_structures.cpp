#include "Phreeqc.h"
#include "phqalloc.h"

#include <list>
#include <string>

/* **********************************************************************
 *
 *   Routines related to structure "pitz_param"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct pitz_param * Phreeqc::
pitz_param_alloc(void)
/* ---------------------------------------------------------------------- */
{
	struct pitz_param *pitz_param_ptr;
	pitz_param_ptr =
		(struct pitz_param *) PHRQ_malloc(sizeof(struct pitz_param));
	if (pitz_param_ptr == NULL)
		malloc_error();
	return (pitz_param_ptr);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
pitz_param_init(struct pitz_param *pitz_param_ptr)
/* ---------------------------------------------------------------------- */
{
	int i;
/*
 *   Frees all data associated with pitz_param structure.
 */

	if (pitz_param_ptr == NULL)
		return (ERROR);
	pitz_param_ptr->species[0] = NULL;
	pitz_param_ptr->species[1] = NULL;
	pitz_param_ptr->species[2] = NULL;
	pitz_param_ptr->ispec[0] = -1;
	pitz_param_ptr->ispec[1] = -1;
	pitz_param_ptr->ispec[2] = -1;
	pitz_param_ptr->type = TYPE_Other;
	pitz_param_ptr->p = 0.0;
	pitz_param_ptr->U.b0 = 0.0;
	for (i = 0; i < 6; i++)
	{
		pitz_param_ptr->a[i] = 0.0;
	}
	pitz_param_ptr->alpha = 0.0;
	pitz_param_ptr->thetas = NULL;
	pitz_param_ptr->os_coef = 0.;
	for (i = 0; i < 3; i++)
	{
		pitz_param_ptr->ln_coef[i] = 0.0;
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
struct pitz_param * Phreeqc::
pitz_param_read(char *string, int n)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read pitzer parameter info from string
 *   n is number of species (character values)
 *          
 */
	int l, i, j, k;
	char *ptr;
	char token[2 * MAX_LENGTH];
	struct pitz_param pzp, *pzp_ptr;

	if (n != 2 && n != 3)
		return (NULL);
	if (string == NULL)
		return (NULL);

	pitz_param_init(&pzp);
	ptr = string;
	if (copy_token(token, &ptr, &l) == EMPTY)
		return (NULL);
	ptr = string;
	for (i = 0; i < n; i++)
	{
		int j = copy_token(token, &ptr, &l);
		if (j == EMPTY)
			return (NULL);
		if (j != UPPER && token[0] != '(')
		{
			input_error++;
			std::ostringstream err;
			err << "Wrong number of species for a Pitzer parameter.\n"  << line;
			error_msg(err.str().c_str(), CONTINUE);
		}
		pzp.species[i] = string_hsave(token);
	}
	k = 0;
	for (i = 0; i < 6; i++)
	{
		if (copy_token(token, &ptr, &l) == EMPTY)
			break;
		j = sscanf(token, SCANFORMAT, &pzp.a[i]);
		if (j <= 0)
			break;
		k++;
	}
	if (k <= 0)
		return (NULL);
	pzp_ptr = pitz_param_duplicate(&pzp);
	return (pzp_ptr);
}

/* ---------------------------------------------------------------------- */
struct pitz_param * Phreeqc::
pitz_param_duplicate(struct pitz_param *old_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Allocates space and makes duplicate copy of pitz_param structure
 */
	struct pitz_param *new_ptr;

	new_ptr = pitz_param_alloc();
	pitz_param_init(new_ptr);
/*
 *   Copy data
 */
	pitz_param_copy(old_ptr, new_ptr);
	return (new_ptr);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
pitz_param_copy(struct pitz_param *old_ptr, struct pitz_param *new_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies pitz_param data from old_ptr to new location, new_ptr.
 *   Space for the new_ptr structure must already be malloced.
 */
/*
 *   Store data for structure pitz_param
 */
	memcpy(new_ptr, old_ptr, sizeof(struct pitz_param));
	return (OK);
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
pitz_param_store(struct pitz_param *pzp_ptr, bool force_copy)
/* ---------------------------------------------------------------------- */
{
/*
 *  Does linear search of pitz_params for same type and species
 *  Returns -1 if not found, index number in pitz_params if found
 */
	int i;
	if (pzp_ptr == NULL)
		return;
	if (pzp_ptr->type == TYPE_Other)
		return;

	std::set< std::string > header;
	for (i = 0; i < 3; i++)
	{
		if (pzp_ptr->species[i] != NULL) header.insert(pzp_ptr->species[i]);
	}

	std::ostringstream key_str;
	key_str << pzp_ptr->type << " ";
	std::set< std::string >::iterator it = header.begin();
	for(; it != header.end(); ++it)
	{
		key_str << *it << " ";
	}
	std::string key = key_str.str().c_str();
	std::map< std::string, size_t>::iterator jit = pitz_param_map.find(key);
	if (jit != pitz_param_map.end())
	{
		if (pzp_ptr->species[2] != NULL)
		{
			error_string = sformatf( "Redefinition of parameter, %s %s %s\n", 
			pzp_ptr->species[0], pzp_ptr->species[1], pzp_ptr->species[2]);
		}
		else
		{
			error_string = sformatf( "Redefinition of parameter, %s %s\n", 
			pzp_ptr->species[0], pzp_ptr->species[1]);
		}
	    warning_msg(error_string);
		pitz_params[(*jit).second] = (struct pitz_param *) free_check_null(pitz_params[(*jit).second]);
		pitz_params[(*jit).second] = pzp_ptr;
	}
	else
	{
		if (count_pitz_param >= max_pitz_param)
		{
			space((void **) ((void *) &pitz_params),
				count_pitz_param, &max_pitz_param,
				sizeof(struct pitz_param *));
		}
		if (force_copy)
		{
			pitz_params[count_pitz_param] = pitz_param_duplicate(pzp_ptr);
			// clean up pointers
			// species 
			for (i = 0; i < 3; i++)
			{
				if (pzp_ptr->species[i] != NULL) 
				{
					pitz_params[count_pitz_param]->species[i] = string_hsave(pzp_ptr->species[i]);
				}
			}
			// thetas
			pitz_params[count_pitz_param]->thetas = NULL;
		}
		else
		{
			pitz_params[count_pitz_param] = pzp_ptr;
		}
		pitz_param_map[key] = count_pitz_param;
		count_pitz_param++;
	}
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
sit_param_store(struct pitz_param *pzp_ptr, bool force_copy)
/* ---------------------------------------------------------------------- */
{
/*
 *  Does linear search of pitz_params for same type and species
 *  Returns -1 if not found, index number in pitz_params if found
 */
	int i;
	if (pzp_ptr == NULL)
		return;
	if (pzp_ptr->type == TYPE_Other)
		return;

	std::set< std::string > header;
	for (i = 0; i < 3; i++)
	{
		if (pzp_ptr->species[i] != NULL) header.insert(pzp_ptr->species[i]);
	}

	std::ostringstream key_str;
	key_str << pzp_ptr->type << " ";
	std::set< std::string >::iterator it = header.begin();
	for(; it != header.end(); ++it)
	{
		key_str << *it << " ";
	}
	std::string key = key_str.str().c_str();

	std::map< std::string, size_t>::iterator jit = sit_param_map.find(key);
	if (jit != sit_param_map.end())
	{
		if (pzp_ptr->species[2] != NULL)
		{
			error_string = sformatf( "Redefinition of parameter, %s %s %s\n", 
			pzp_ptr->species[0], pzp_ptr->species[1], pzp_ptr->species[2]);
		}
		else
		{
			error_string = sformatf( "Redefinition of parameter, %s %s\n", 
			pzp_ptr->species[0], pzp_ptr->species[1]);
		}
	    warning_msg(error_string);
		sit_params[(*jit).second] = (struct pitz_param *) free_check_null(sit_params[(*jit).second]);
		sit_params[(*jit).second] = pzp_ptr;
	}
	else
	{
		if (count_sit_param >= max_sit_param)
		{
			space((void **) ((void *) &sit_params),
				count_sit_param, &max_sit_param,
				sizeof(struct pitz_param *));
		}
		if (force_copy)
		{
			sit_params[count_sit_param] = pitz_param_duplicate(pzp_ptr);
			// clean up pointers
			// species 
			for (i = 0; i < 3; i++)
			{
				if (pzp_ptr->species[i] != NULL) 
				{
					sit_params[count_sit_param]->species[i] = string_hsave(pzp_ptr->species[i]);
				}
			}
			// thetas
			sit_params[count_sit_param]->thetas = NULL;
		}
		else
		{
			sit_params[count_sit_param] = pzp_ptr;
		}
		sit_param_map[key] = count_sit_param;
		count_sit_param++;
	}
}

/* **********************************************************************
 *
 *   Routines related to structure "theta_parm"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct theta_param * Phreeqc::
theta_param_alloc(void)
/* ---------------------------------------------------------------------- */
{
	struct theta_param *theta_param_ptr;
	theta_param_ptr =
		(struct theta_param *) PHRQ_malloc(sizeof(struct theta_param));
	if (theta_param_ptr == NULL)
		malloc_error();
	return (theta_param_ptr);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
theta_param_init(struct theta_param *theta_param_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Frees all data associated with theta_param structure.
 */

	if (theta_param_ptr == NULL)
		return (ERROR);
	theta_param_ptr->zj = 0;
	theta_param_ptr->zk = 0;
	theta_param_ptr->etheta = 0;
	theta_param_ptr->ethetap = 0;
	return (OK);
}

/* ---------------------------------------------------------------------- */
struct theta_param * Phreeqc::
theta_param_search(LDBLE zj, LDBLE zk)
/* ---------------------------------------------------------------------- */
{
/*
 *  Does linear search of theta_params for same charge
 *  Returns NULL if not found, index number in theta_params if found
 */
	int i;
	for (i = 0; i < count_theta_param; i++)
	{
		if ((theta_params[i]->zj == zj && theta_params[i]->zk == zk) ||
			(theta_params[i]->zj == zk && theta_params[i]->zk == zj))
		{
			return theta_params[i];
		}
	}
	return NULL;
}
