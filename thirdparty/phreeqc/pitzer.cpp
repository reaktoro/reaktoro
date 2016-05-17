#include "Phreeqc.h"
#include "phqalloc.h"
#include "Exchange.h"
#include "Solution.h"
#define PITZER_LISTS
#define PITZER

/* ---------------------------------------------------------------------- */
int Phreeqc::
pitzer_init(void)
/* ---------------------------------------------------------------------- */
{
	int i;
/*
 *      Initialization for pitzer
 */
	pitzer_model = FALSE;
	max_pitz_param = 100;
	count_pitz_param = 0;
	use_etheta = TRUE;
	space((void **) ((void *) &pitz_params), INIT, &max_pitz_param,
		  sizeof(struct pitz_param *));

	max_theta_param = 100;
	count_theta_param = 0;
	space((void **) ((void *) &theta_params), INIT, &max_theta_param,
		  sizeof(struct theta_param *));

	ICON = TRUE;
	OTEMP = -100.;
	OPRESS = -100.;
	for (i = 0; i < 23; i++)
	{
		BK[i] = 0.0;
		DK[i] = 0.0;
	}
	pitzer_pe = FALSE;
	VP = 0;
	DW0 = 0;
	return OK;
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
pitzer_tidy(void)
/* ---------------------------------------------------------------------- */
{
	/*
	 *      Make lists of species for cations, anions, neutral
	 */
	const char *string1, *string2;
	int i, j, order;
	int i0, i1, i2;
	int count_pos, count_neg, count_neut, count[3], jj;
	LDBLE z0, z1;
	struct pitz_param *pzp_ptr;
	struct theta_param *theta_param_ptr;
	/*
	* Ensure new parameters are calculated
	*/
	OTEMP = -100.;
	OPRESS = -100.;
	/*
	 *  allocate pointers to species structures
	 */
	if (spec != NULL)
		spec = (struct species **) free_check_null(spec);
	spec =
		(struct species **)
		PHRQ_malloc((size_t) (3 * count_s * sizeof(struct species *)));
	if (spec == NULL)
		malloc_error();
	for (i = 0; i < 3 * count_s; i++)
		spec[i] = NULL;
	cations = spec;
	neutrals = &(spec[count_s]);
	anions = &(spec[2 * count_s]);
	MAXCATIONS = count_s;
	FIRSTANION = 2 * count_s;
	MAXNEUTRAL = count_s;
	count_cations = 0;
	count_anions = 0;
	count_neutrals = 0;
	if (itmax < 200)
		itmax = 200;
	/*
	 *  allocate other arrays for Pitzer
	 */
	if (IPRSNT != NULL)
		IPRSNT = (int *) free_check_null(IPRSNT);
	IPRSNT = (int *) PHRQ_malloc((size_t) (3 * count_s * sizeof(int)));
	if (IPRSNT == NULL)
		malloc_error();
	if (M != NULL)
		M = (LDBLE *) free_check_null(M);
	M = (LDBLE *) PHRQ_malloc((size_t) (3 * count_s * sizeof(LDBLE)));
	if (M == NULL)
		malloc_error();
	if (LGAMMA != NULL)
		LGAMMA = (LDBLE *) free_check_null(LGAMMA);
	LGAMMA = (LDBLE *) PHRQ_malloc((size_t) (3 * count_s * sizeof(LDBLE)));
	if (LGAMMA == NULL)
		malloc_error();


	for (i = 0; i < count_s; i++)
	{
		if (s[i] == s_eminus)
			continue;
		if (s[i] == s_h2o)
			continue;
		if (s[i]->type == EX || s[i]->type == SURF) 
			continue;
		if (s[i]->z < -.001)
		{
			anions[count_anions++] = s[i];
		}
		else if (s[i]->z > .001)
		{
			cations[count_cations++] = s[i];
		}
		else
		{
			neutrals[count_neutrals++] = s[i];
		}
	}
	/*
	 *  Add etheta to parameter list in case theta not defined for 
	 *  cation-cation or anion-anion pair
	 *  Remove old TYPE_ETHETA definitions
	 */
	j = 0;
	for (i = 0; i < count_pitz_param; i++)
	{
		if (pitz_params[i]->type == TYPE_ETHETA)
		{
			pitz_params[i] =
				(struct pitz_param *) free_check_null(pitz_params[i]);
		}
		else
		{
			pitz_params[j++] = pitz_params[i];
		}
	}
	count_pitz_param = j;
	for (i = 0; i < count_cations - 1; i++)
	{
		for (j = i + 1; j < count_cations; j++)
		{
			sprintf(line, "%s %s 1", spec[i]->name, spec[j]->name);
			pzp_ptr = pitz_param_read(line, 2);
			pzp_ptr->type = TYPE_ETHETA;
			if (count_pitz_param >= max_pitz_param)
			{
				space((void **) ((void *) &pitz_params), count_pitz_param,
					  &max_pitz_param, sizeof(struct pitz_param *));
			}
			pitz_params[count_pitz_param++] = pzp_ptr;

		}
	}
	for (i = 2 * count_s; i < 2 * count_s + count_anions - 1; i++)
	{
		for (j = i + 1; j < 2 * count_s + count_anions; j++)
		{
			sprintf(line, "%s %s 1", spec[i]->name, spec[j]->name);
			pzp_ptr = pitz_param_read(line, 2);
			pzp_ptr->type = TYPE_ETHETA;
			if (count_pitz_param >= max_pitz_param)
			{
				space((void **) ((void *) &pitz_params), count_pitz_param,
					  &max_pitz_param, sizeof(struct pitz_param *));
			}
			pitz_params[count_pitz_param] = pzp_ptr;
			count_pitz_param++;
		}
	}
	/*
	 *  put species numbers in pitz_params
	 */
	for (i = 0; i < count_pitz_param; i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (pitz_params[i]->species[j] == NULL)
				continue;
			pitz_params[i]->ispec[j] = ISPEC(pitz_params[i]->species[j]);
			if ((j < 2 && pitz_params[i]->ispec[j] == -1) ||
				(j == 2
				 && (pitz_params[i]->type == TYPE_PSI
					 || pitz_params[i]->type == TYPE_ZETA)
				 && pitz_params[i]->ispec[j] == -1))
			{
				input_error++;
				error_string = sformatf(
						"Species for Pitzer parameter not defined in SOLUTION_SPECIES, %s",
						pitz_params[i]->species[j]);
				error_msg(error_string, CONTINUE);
				return (ERROR);
			}
		}
	}
	/*
	 * MacInnes data
	 */
	string1 = string_hsave("K+");
	string2 = string_hsave("Cl-");
	IC = ISPEC(string2);
	for (i = 0; i < count_pitz_param; i++)
	{
		if ((pitz_params[i]->species[0] == string1 &&
			pitz_params[i]->species[1] == string2) ||
			(pitz_params[i]->species[0] == string2 &&
			pitz_params[i]->species[1] == string1) )
		{
			switch (pitz_params[i]->type)
			{
			case TYPE_B0:
				mcb0 = pitz_params[i];
				break;
			case TYPE_B1:
				mcb1 = pitz_params[i];
				break;
			case TYPE_C0:
				mcc0 = pitz_params[i];
				break;
			case TYPE_B2:
			case TYPE_THETA:
			case TYPE_LAMDA:
			case TYPE_ZETA:
			case TYPE_PSI:
			case TYPE_ETHETA:
			case TYPE_ALPHAS:
			case TYPE_MU:
			case TYPE_ETA:
			case TYPE_Other:
			default:
				break;
			}
		}
	}
	if (mcb0 == NULL && mcb1 == NULL && mcc0 == NULL && ICON == TRUE)
	{
		error_string = sformatf(
				"No KCl interaction parameters, turning off MacInnis scaling.");
		warning_msg(error_string);
		ICON = FALSE;
	}
	/*
	 * Set alpha values
	 */
	for (i = 0; i < count_pitz_param; i++)
	{
		z0 = fabs(spec[pitz_params[i]->ispec[0]]->z);
		z1 = fabs(spec[pitz_params[i]->ispec[1]]->z);
		if (equal(z0, 1.0, 1e-8) || equal(z1, 1.0, 1e-8))
		{
			order = 1;
		}
		else if (equal(z0, 2.0, 1e-8) && equal(z1, 2.0, 1e-8))
		{
			order = 2;
		}
		else
		{
			order = 3;
		}
		if (pitz_params[i]->type == TYPE_B1)
		{
			switch (order)
			{
			case 1:
			case 3:
				pitz_params[i]->alpha = 2.0;
				break;
			case 2:
				pitz_params[i]->alpha = 1.4;
				break;
			}
		}
		else if (pitz_params[i]->type == TYPE_B2)
		{
			switch (order)
			{
			case 1:
				pitz_params[i]->alpha = 12.0;
				break;
			case 2:
				pitz_params[i]->alpha = 12.0;
				break;
			case 3:
				pitz_params[i]->alpha = 50.0;
				break;
			}
		}
	}
	/*
	 * Add specific alphas
	 */
	for (i = 0; i < count_pitz_param; i++)
	{
		if (pitz_params[i]->type == TYPE_ALPHAS)
		{
			for (j = 0; j < count_pitz_param; j++)
			{
				if (pitz_params[j]->type != TYPE_B1)
					continue;
				if (pitz_params[i]->ispec[0] != pitz_params[j]->ispec[0])
					continue;
				if (pitz_params[i]->ispec[1] != pitz_params[j]->ispec[1])
					continue;
				pitz_params[j]->alpha = pitz_params[i]->a[0];
				break;
			}
			for (j = 0; j < count_pitz_param; j++)
			{
				if (pitz_params[j]->type != TYPE_B2)
					continue;
				if (pitz_params[i]->ispec[0] != pitz_params[j]->ispec[0])
					continue;
				if (pitz_params[i]->ispec[1] != pitz_params[j]->ispec[1])
					continue;
				pitz_params[j]->alpha = pitz_params[i]->a[1];
				break;
			}
		}
	}

	/*
	 *   Add thetas pointer to etheta pitzer parameters
	 */

	if (count_theta_param > 0)
	{
		for (i = 0; i < count_theta_param; i++)
		{
			theta_params[i] =
				(struct theta_param *) free_check_null(theta_params[i]);
		}
	}
	count_theta_param = 0;
	for (i = 0; i < count_pitz_param; i++)
	{
		if (pitz_params[i]->type == TYPE_ETHETA)
		{
			z0 = spec[pitz_params[i]->ispec[0]]->z;
			z1 = spec[pitz_params[i]->ispec[1]]->z;
			theta_param_ptr = theta_param_search(z0, z1);
			if (theta_param_ptr == NULL)
			{
				if (count_theta_param >= max_theta_param)
				{
					space((void **) ((void *) &theta_params),
						  count_theta_param, &max_theta_param,
						  sizeof(struct theta_param *));
				}
				theta_params[count_theta_param] = theta_param_alloc();
				theta_param_init(theta_params[count_theta_param]);
				theta_params[count_theta_param]->zj = z0;
				theta_params[count_theta_param]->zk = z1;
				theta_param_ptr = theta_params[count_theta_param];
				count_theta_param++;
			}
			pitz_params[i]->thetas = theta_param_ptr;
		}
	}
	/*
	 *  Tidy TYPE_MU
	 */

	/* Coef for Osmotic coefficient for TYPE_MU */

	for (i = 0; i < count_pitz_param; i++)
	{
		if (pitz_params[i]->type == TYPE_MU)
		{
			i0 = pitz_params[i]->ispec[0];
			i1 = pitz_params[i]->ispec[1];
			i2 = pitz_params[i]->ispec[2];
			count_pos = count_neg = count_neut = 0;
			for (j = 0; j <= 2; j++)
			{
				if (spec[pitz_params[i]->ispec[j]]->z > 0)
				{
					count_pos++;
				}
				if (spec[pitz_params[i]->ispec[j]]->z == 0)
				{
					count_neut++;
				}
				if (spec[pitz_params[i]->ispec[j]]->z < 0)
				{
					count_neg++;
				}
			}
			/* All neutral */
			if (count_neut == 3)
			{
				if (i0 == i1 && i1 == i2)
				{
					/* type n, n, n */
					pitz_params[i]->os_coef = 1;
					continue;
				}
				else if (i0 == i1 || i1 == i2 || i0 == i2)
				{
					/* type n, n, n' */
					pitz_params[i]->os_coef = 3;
					continue;
				}
				else
				{
					/* type n, n', n'' */
					pitz_params[i]->os_coef = 6;
					continue;
				}
			}
			/* Two neutral, one anion or cation */
			if (i0 == i1 || i1 == i2 || i0 == i2)
			{
				/* type n, n, a|c */
				pitz_params[i]->os_coef = 3;
				continue;
			}
			else
			{
				/* type n, n', a|c */
				pitz_params[i]->os_coef = 6;
				continue;
			}
		}
	}

	/* Coef for gammas for TYPE_MU */

	for (i = 0; i < count_pitz_param; i++)
	{
		if (pitz_params[i]->type == TYPE_MU)
		{
			for (j = 0; j <= 2; j++)
			{
				count[j] = 0;
				for (jj = 0; jj <= 2; jj++)
				{
					if (pitz_params[i]->ispec[j] == pitz_params[i]->ispec[jj])
					{
						count[j]++;
					}
				}
			}
			for (j = 0; j <= 2; j++)
			{
				/* cation or anion */
				if (spec[pitz_params[i]->ispec[j]]->z < 0
					|| spec[pitz_params[i]->ispec[j]]->z > 0)
				{
					if (count[0] > 1 || count[1] > 1)
					{
						pitz_params[i]->ln_coef[j] = 3;
					}
					else
					{
						pitz_params[i]->ln_coef[j] = 6;
					}
					continue;
				}
				/* Neutral */
				if (count[j] == 3)
				{
					pitz_params[i]->ln_coef[j] = 1;
				}
				else if (count[j] == 2)
				{
					pitz_params[i]->ln_coef[j] = 3;
				}
				else if (count[j] == 1)
				{
					if (count[0] > 1 || count[1] > 1)
					{
						pitz_params[i]->ln_coef[j] = 3;
					}
					else
					{
						pitz_params[i]->ln_coef[j] = 6;
					}
				}
			}
		}
	}
	/*  Debug TYPE_MU coefficients */
	/*
	   for (i = 0; i < count_pitz_param; i++)
	   {
	   if (pitz_params[i]->type == TYPE_MU)
	   {
	   fprintf(stderr, "%s\t%s\t%s\n", pitz_params[i]->species[0], pitz_params[i]->species[1], pitz_params[i]->species[2]);
	   fprintf(stderr, "%f\t%f\t%f\n", pitz_params[i]->ln_coef[0], pitz_params[i]->ln_coef[1], pitz_params[i]->ln_coef[2]);
	   fprintf(stderr, "%f\n\n", pitz_params[i]->os_coef);
	   }
	   }
	 */
	/*
	 *  Tidy TYPE_LAMDA
	 */

	/* Coef for Osmotic coefficient for TYPE_LAMDA */

	for (i = 0; i < count_pitz_param; i++)
	{
		if (pitz_params[i]->type == TYPE_LAMDA)
		{
			i0 = pitz_params[i]->ispec[0];
			i1 = pitz_params[i]->ispec[1];
			/* All neutral */
			if (i0 == i1)
			{
				/* type n, n */
				pitz_params[i]->os_coef = 0.5;
				pitz_params[i]->ln_coef[0] = 1;
				pitz_params[i]->ln_coef[1] = 1;
			}
			else
			{
				/* type nn', na, nc */
				pitz_params[i]->os_coef = 1;
				pitz_params[i]->ln_coef[0] = 2;
				pitz_params[i]->ln_coef[1] = 2;
			}
		}
	}
	/*  Debug TYPE_LAMDA coefficients */
	/*
	   for (i = 0; i < count_pitz_param; i++)
	   {
	   if (pitz_params[i]->type == TYPE_LAMDA)
	   {
	   fprintf(stderr, "%s\t%s\n", pitz_params[i]->species[0], pitz_params[i]->species[1]);
	   fprintf(stderr, "%f\t%f\n", pitz_params[i]->ln_coef[0], pitz_params[i]->ln_coef[1]);
	   fprintf(stderr, "%f\n\n", pitz_params[i]->os_coef);
	   }
	   }
	 */
	/* remake map */
	{
		pitz_param_map.clear();
		for (int j = 0; j < count_pitz_param; j++)
		{	
			std::set< std::string > header;
			for (int i = 0; i < 3; i++)
			{
				if (pitz_params[j]->species[i] != NULL) header.insert(pitz_params[j]->species[i]);
			}
			std::ostringstream key_str;
			key_str << pitz_params[j]->type << " ";
			std::set< std::string >::iterator it = header.begin();
			for(; it != header.end(); ++it)
			{
				key_str << *it << " ";
			}
			std::string key = key_str.str().c_str();
			pitz_param_map[key] = j;
		}
		assert ((int) pitz_param_map.size() == count_pitz_param);
	}
	return OK;
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
ISPEC(const char *name)
/* ---------------------------------------------------------------------- */
/*
 *      Find species number in spec for character string species name
 */
{
	int i;
	for (i = 0; i < 3 * count_s; i++)
	{
		if (spec[i] == NULL)
			continue;
		if (name == spec[i]->name)
		{
			return (i);
		}
	}
	return (-1);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_pitzer(void)
/* ---------------------------------------------------------------------- */
{
	/*
	 *      Reads advection information
	 *
	 *      Arguments:
	 *         none
	 *
	 *      Returns:
	 *         KEYWORD if keyword encountered, input_error may be incremented if
	 *                    a keyword is encountered in an unexpected position
	 *         EOF     if eof encountered while reading mass balance concentrations
	 *         ERROR   if error occurred reading data
	 *
	 */
	/*
	 *   Read advection parameters: 
	 *        number of cells;
	 *        number of shifts;
	 */
	int n;
	struct pitz_param *pzp_ptr;
	pitz_param_type pzp_type;

	int return_value, opt, opt_save;
	char *next_char;
	const char *opt_list[] = {
		"b0",					/* 0 */
		"b1",					/* 1 */
		"b2",					/* 2 */
		"c0",					/* 3 */
		"theta",				/* 4 */
		"lamda",				/* 5 */
		"zeta",					/* 6 */
		"psi",					/* 7 */
		"macinnes",				/* 8 */
		"macinnis",				/* 9 */
		"mac",					/* 10 */
		"redox",				/* 11 */
		"pe",					/* 12 */
		"alphas",				/* 13 */
		"mu",					/* 14 */
		"eta",					/* 15 */
		"etheta",				/* 16 */
		"use_etheta",			/* 17 */
		"lambda"                /* 18 */
	};
	int count_opt_list = 19;
	/*
	 *   Read lines
	 */
	opt_save = OPTION_ERROR;
	return_value = UNKNOWN;
	n = -1;
	pzp_type = TYPE_Other;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT)
		{
			opt = opt_save;
		}
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_DEFAULT:
			pzp_ptr = pitz_param_read(line, n);
			if (pzp_ptr != NULL)
			{
				pzp_ptr->type = pzp_type;
				pitz_param_store(pzp_ptr, false);
			}
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in PITZER keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* b0 */
			pzp_type = TYPE_B0;
			n = 2;
			opt_save = OPTION_DEFAULT;
			break;
		case 1:				/* b1 */
			pzp_type = TYPE_B1;
			n = 2;
			opt_save = OPTION_DEFAULT;
			break;
		case 2:				/* b2 */
			pzp_type = TYPE_B2;
			n = 2;
			opt_save = OPTION_DEFAULT;
			break;
		case 3:				/* c0 */
			pzp_type = TYPE_C0;
			n = 2;
			opt_save = OPTION_DEFAULT;
			break;
		case 4:				/* theta */
			pzp_type = TYPE_THETA;
			n = 2;
			opt_save = OPTION_DEFAULT;
			break;
		case 5:				/* lamda */
		case 18:            /* lambda */
			pzp_type = TYPE_LAMDA;
			n = 2;
			opt_save = OPTION_DEFAULT;
			break;
		case 6:				/* zeta */
			pzp_type = TYPE_ZETA;
			n = 3;
			opt_save = OPTION_DEFAULT;
			break;
		case 7:				/* psi */
			pzp_type = TYPE_PSI;
			n = 3;
			opt_save = OPTION_DEFAULT;
			break;
		case 13:				/* alphas */
			pzp_type = TYPE_ALPHAS;
			n = 2;
			opt_save = OPTION_DEFAULT;
			break;
		case 8:				/* macinnes */
		case 9:				/* macinnis */
		case 10:				/* mac */
			opt_save = OPTION_ERROR;
			ICON = get_true_false(next_char, TRUE);
			break;
		case 11:				/* redox */
		case 12:				/* pe */
			opt_save = OPTION_ERROR;
			pitzer_pe = get_true_false(next_char, TRUE);
			break;
		case 14:				/* mu */
			pzp_type = TYPE_MU;
			n = 3;
			opt_save = OPTION_DEFAULT;
			break;
		case 15:				/* eta */
			pzp_type = TYPE_ETA;
			n = 3;
			opt_save = OPTION_DEFAULT;
			break;
		case 16:				/* etheta */
		case 17:				/* use_etheta */
			opt_save = OPTION_ERROR;
			use_etheta = get_true_false(next_char, TRUE);
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
	pitzer_model = TRUE;
	return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
PTEMP(LDBLE TK)
/* ---------------------------------------------------------------------- */
{
/*
C
C     SUBROUTINE TO CALUCLATE TEMPERATURE DEPENDENCE OF PITZER PARAMETER
C
*/
	LDBLE TR = 298.15;

	if (fabs(TK - OTEMP) < 0.001 && fabs(patm_x - OPRESS) < 0.1)
		return OK;
	DW0 = rho_0 = calc_rho_0(TK - 273.15, patm_x);
	VP = patm_x;
#if !defined(PITZER_LISTS)
	int i;
	for (i = 0; i < count_pitz_param; i++)
	{
		calc_pitz_param(pitz_params[i], TK, TR);
	}
#else
	for (size_t j = 0; j < param_list.size(); j++)
	{
		int i = param_list[j];
		calc_pitz_param(pitz_params[i], TK, TR);
	}
	if (mcb0) 
	{
		calc_pitz_param(mcb0, TK, TR);
	}
	if (mcb1)
	{
		calc_pitz_param(mcb1, TK, TR);
	}
	if (mcc0)
	{
		calc_pitz_param(mcc0, TK, TR);
	}
#endif
	calc_dielectrics(TK - 273.15, patm_x);
	OTEMP = TK;
	OPRESS = patm_x;
	return OK;
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
calc_pitz_param(struct pitz_param *pz_ptr, LDBLE TK, LDBLE TR)
/* ---------------------------------------------------------------------- */
{
	LDBLE param;
	/*
	 */
	if (fabs(TK - TR) < 0.001)
	{
		param = pz_ptr->a[0];
	}
	else
	{
		param = (pz_ptr->a[0] +
				 pz_ptr->a[1] * (1.e0 / TK - 1.e0 / TR) +
				 pz_ptr->a[2] * log(TK / TR) +
				 pz_ptr->a[3] * (TK - TR) + 
				 pz_ptr->a[4] * (TK * TK - TR * TR)) +
				 pz_ptr->a[5] * (1.e0 / (TK * TK) - 1.e0 / (TR * TR));
	}
	pz_ptr->p = param;
	switch (pz_ptr->type)
	{
	case TYPE_B0:
		pz_ptr->U.b0 = param;
		break;
	case TYPE_B1:
		pz_ptr->U.b1 = param;
		break;
	case TYPE_B2:
		pz_ptr->U.b2 = param;
		break;
	case TYPE_C0:
		pz_ptr->U.c0 = param;
		break;
	case TYPE_THETA:
		pz_ptr->U.theta = param;
		break;
	case TYPE_LAMDA:
		pz_ptr->U.lamda = param;
		break;
	case TYPE_ZETA:
		pz_ptr->U.zeta = param;
		break;
	case TYPE_ETHETA:
		break;
	case TYPE_PSI:
		pz_ptr->U.psi = param;
		break;
	case TYPE_ALPHAS:
		break;
	case TYPE_MU:
		pz_ptr->U.mu = param;
		break;
	case TYPE_ETA:
		pz_ptr->U.eta = param;
		break;
	case TYPE_Other:
	default:
		error_msg("Should not be TYPE_Other in function calc_pitz_param",
				  STOP);
		break;
	}
	return OK;
}
#if !defined(PITZER_LISTS)
/* ---------------------------------------------------------------------- */
int Phreeqc::
pitzer(void)
/* ---------------------------------------------------------------------- */
{
	int i, i0, i1, i2;
	LDBLE param, l_alpha, z0, z1;
	LDBLE etheta, ethetap;
	/*
	   LDBLE CONV, XI, XX, OSUM, BIGZ, DI, F, XXX, GAMCLM, 
	   CSUM, PHIMAC, OSMOT, BMXP, ETHEAP, CMX, BMX, PHI,
	   BMXPHI, PHIPHI, AW, A, B;
	 */
	LDBLE CONV, XX, OSUM, BIGZ, DI, F, F1, F2, F_var, XXX, GAMCLM, CSUM, PHIMAC, OSMOT,
		B, B1, B2;
	LDBLE I, TK;
	/*
	   C
	   C     INITIALIZE
	   C
	 */
	CONV = 1.0 / log(10.0);
	XX = 0.0;
	OSUM = 0.0;
	/*n
	   I = *I_X;
	   TK = *TK_X;
	 */
	I = mu_x;
	TK = tk_x;
	/*      DH_AB(TK, &A, &B); */
	/*
	   C
	   C     TRANSFER DATA FROM TO M
	   C
	 */
	for (i = 0; i < 3 * count_s; i++)
	{
		IPRSNT[i] = FALSE;
		M[i] = 0.0;
		if (spec[i] != NULL && spec[i]->in == TRUE)
		{
			if (spec[i]->type == EX ||
				spec[i]->type == SURF || spec[i]->type == SURF_PSI)
				continue;
			M[i] = under(spec[i]->lm);
			if (M[i] > MIN_TOTAL)
				IPRSNT[i] = TRUE;
		}
	}
	if (ICON == TRUE)
	{
		IPRSNT[IC] = TRUE;
	}
	/*
	   ICON = 0;
	   M[1] = 1.40070736;
	   M[4] = 2.52131086E-05;
	   M[140] = 4.59985435E-09;
	 */

	/*
	   C
	   C     COMPUTE PITZER COEFFICIENTS' TEMPERATURE DEPENDENCE
	   C
	 */
	PTEMP(TK);
	for (i = 0; i < 2 * count_s + count_anions; i++)
	{
		LGAMMA[i] = 0.0;
		if (IPRSNT[i] == TRUE)
		{
			XX = XX + M[i] * fabs(spec[i]->z);
			OSUM = OSUM + M[i];
		}
	}
	/*
	   C
	   C     EQUATION (8)
	   C
	 */
	BIGZ = XX;
	DI = sqrt(I);
	/*
	   C
	   C     CALCULATE F & GAMCLM
	   C
	 */
	B = 1.2;
	F = F1 = F2 = -A0 * (DI / (1.0 + B * DI) + 2.0 * log(1.0 + B * DI) / B);
	if (patm_x > 1.0)
	{
		LDBLE pap;
		pap = (7e-5 + 1.93e-9 * pow(TK - 250.0, 2.0)) * patm_x;
		B1 = B - (pap > 0.2 ? 0.2 : pap);
		pap = (9.65e-10 * pow(TK - 263.0, 2.773)) * pow(patm_x, 0.623);
		//pap = (-5.22e-4 + 7.19e-8 * pow(TK - 263.0, 2.0)) * pow(patm_x, 0.623);
		B2 = B - (pap > 0.2 ? 0.2 : pap);
		if (B1 != 0)
			F1 = -A0 * (DI / (1.0 + B1 * DI) + 2.0 * log(1.0 + B1 * DI) / B1);
		if (B2 != 0)
			F2 = -A0 * (DI / (1.0 + B2 * DI) + 2.0 * log(1.0 + B2 * DI) / B2);
	}
	XXX = 2.0 * DI;
	XXX =
		(1.0 - (1.0 + XXX - XXX * XXX * 0.5) * exp(-XXX)) / (XXX * XXX);
	/*GAMCLM=F+I*2.0e0*(BCX(1,IK,IC)+BCX(2,IK,IC)*XXX)+1.5e0*BCX(4,IK,IC)*I*I; */
	/*GAMCLM=F+I*2.0e0*(mcb0->U.b0 + mcb1->U.b1*XXX) + 1.5e0*mcc0->U.c0*I*I; */
	/*GAMCLM = F + I * 2.0e0 * (mcb0->p + mcb1->p * XXX) + 1.5e0 * mcc0->p * I * I; */
	GAMCLM = F1;
	if (mcb0 != NULL)
		GAMCLM += I * 2.0 * mcb0->p;
	if (mcb1 != NULL)
		GAMCLM += I * 2.0 * mcb1->p * XXX;
	if (mcc0 != NULL)
		GAMCLM += 1.5 * mcc0->p * I * I;
	CSUM = 0.0;
	OSMOT = -(A0) * pow(I, (LDBLE) 1.5) / (1.0 + B * DI);
	/*
	 *  Calculate ethetas
	 */
	for (i = 0; i < count_theta_param; i++)
	{
		z0 = theta_params[i]->zj;
		z1 = theta_params[i]->zk;
		ETHETAS(z0, z1, I, &etheta, &ethetap);
		theta_params[i]->etheta = etheta;
		theta_params[i]->ethetap = ethetap;
	}
	/*
	 *  Sums for F, LGAMMA, and OSMOT
	 */
	for (i = 0; i < count_pitz_param; i++)
	{
		i0 = pitz_params[i]->ispec[0];
		i1 = pitz_params[i]->ispec[1];
		if (IPRSNT[i0] == FALSE || IPRSNT[i1] == FALSE)
			continue;
		z0 = spec[i0]->z;
		z1 = spec[i1]->z;
		param = pitz_params[i]->p;
		l_alpha = pitz_params[i]->alpha;
		F_var = 0;
		switch (pitz_params[i]->type)
		{
		case TYPE_B0:
			LGAMMA[i0] += M[i1] * 2.0 * param;
			LGAMMA[i1] += M[i0] * 2.0 * param;
			OSMOT += M[i0] * M[i1] * param;
			break;
		case TYPE_B1:
			if (param != 0.0)
			{
				F_var = M[i0] * M[i1] * param * GP(l_alpha * DI) / I;
				LGAMMA[i0] += M[i1] * 2.0 * param * G(l_alpha * DI);
				LGAMMA[i1] += M[i0] * 2.0 * param * G(l_alpha * DI);
				OSMOT += M[i0] * M[i1] * param * exp(-l_alpha * DI);
			}
			break;
		case TYPE_B2:
			if (param != 0.0)
			{
				F_var = M[i0] * M[i1] * param * GP(l_alpha * DI) / I;
				LGAMMA[i0] += M[i1] * 2.0 * param * G(l_alpha * DI);
				LGAMMA[i1] += M[i0] * 2.0 * param * G(l_alpha * DI);
				OSMOT += M[i0] * M[i1] * param * exp(-l_alpha * DI);
			}
			break;
		case TYPE_C0:
			CSUM +=
				M[i0] * M[i1] * pitz_params[i]->p / (2.0 *
													 sqrt(fabs(z0 * z1)));
			LGAMMA[i0] += M[i1] * BIGZ * param / (2.0 * sqrt(fabs(z0 * z1)));
			LGAMMA[i1] += M[i0] * BIGZ * param / (2.0 * sqrt(fabs(z0 * z1)));
			OSMOT +=
				M[i0] * M[i1] * BIGZ * param / (2.0 * sqrt(fabs(z0 * z1)));
			break;
		case TYPE_THETA:
			LGAMMA[i0] += 2.0 * M[i1] * (param /*+ ETHETA(z0, z1, I) */ );
			LGAMMA[i1] += 2.0 * M[i0] * (param /*+ ETHETA(z0, z1, I) */ );
			OSMOT += M[i0] * M[i1] * param;
			break;
		case TYPE_ETHETA:
			/*
			   ETHETAS(z0, z1, I, &etheta, &ethetap);
			 */
			if (use_etheta == TRUE)
			{
				etheta = pitz_params[i]->thetas->etheta;
				ethetap = pitz_params[i]->thetas->ethetap;
				F_var = M[i0] * M[i1] * ethetap;
				LGAMMA[i0] += 2.0 * M[i1] * etheta;
				LGAMMA[i1] += 2.0 * M[i0] * etheta;
				OSMOT += M[i0] * M[i1] * (etheta + I * ethetap);
				/*
				   F += M[i0]*M[i1]*ETHETAP(z0, z1, I);
				   LGAMMA[i0] += 2.0*M[i1]*(ETHETA(z0, z1, I) ); 
				   LGAMMA[i1] += 2.0*M[i0]*(ETHETA(z0, z1, I) ); 
				   OSMOT += M[i0]*M[i1]*(ETHETA(z0, z1, I) + I*ETHETAP(z0, z1, I) ); 
				 */
			}
			break;
		case TYPE_PSI:
			i2 = pitz_params[i]->ispec[2];
			if (IPRSNT[i2] == FALSE)
				continue;
			LGAMMA[i0] += M[i1] * M[i2] * param;
			LGAMMA[i1] += M[i0] * M[i2] * param;
			LGAMMA[i2] += M[i0] * M[i1] * param;
			OSMOT += M[i0] * M[i1] * M[i2] * param;
			break;
		case TYPE_LAMDA:
			LGAMMA[i0] += M[i1] * param * pitz_params[i]->ln_coef[0];
			LGAMMA[i1] += M[i0] * param * pitz_params[i]->ln_coef[1];
			OSMOT += M[i0] * M[i1] * param * pitz_params[i]->os_coef;
			break;
		case TYPE_ZETA:
			i2 = pitz_params[i]->ispec[2];
			if (IPRSNT[i2] == FALSE)
				continue;
			LGAMMA[i0] += M[i1] * M[i2] * param;
			LGAMMA[i1] += M[i0] * M[i2] * param;
			LGAMMA[i2] += M[i0] * M[i1] * param;
			OSMOT += M[i0] * M[i1] * M[i2] * param;
			break;
		case TYPE_MU:
			i2 = pitz_params[i]->ispec[2];
			if (IPRSNT[i2] == FALSE)
				continue;

			LGAMMA[i0] += M[i1] * M[i2] * param * pitz_params[i]->ln_coef[0];
			LGAMMA[i1] += M[i0] * M[i2] * param * pitz_params[i]->ln_coef[1];
			LGAMMA[i2] += M[i0] * M[i1] * param * pitz_params[i]->ln_coef[2];
			OSMOT += M[i0] * M[i1] * M[i2] * param * pitz_params[i]->os_coef;
			break;
		case TYPE_ETA:
			i2 = pitz_params[i]->ispec[2];
			if (IPRSNT[i2] == FALSE)
				continue;
			LGAMMA[i0] += M[i1] * M[i2] * param;
			LGAMMA[i1] += M[i0] * M[i2] * param;
			LGAMMA[i2] += M[i0] * M[i1] * param;
			OSMOT += M[i0] * M[i1] * M[i2] * param;
			break;
		case TYPE_ALPHAS:
			break;
		case TYPE_Other:
		default:
			error_msg("TYPE_Other in pitz_param list.", STOP);
			break;
		}
	F += F_var;
	F1 += F_var;
	F2 += F_var;
	}

	/*
	 *  Add F and CSUM terms to LGAMMA
	 */

	for (i = 0; i < count_cations; i++)
	{
		if (!IPRSNT[i])
			continue;
		z0 = fabs(spec[i]->z);
		F_var = (z0 == 1 ? F1 : (z0 == 2.0 ? F2 : F));
		LGAMMA[i] += z0 * z0 * F_var + z0 * CSUM;
	}
	for (i = 2 * count_s; i < 2 * count_s + count_anions; i++)
	{
		if (!IPRSNT[i])
			continue;
		z0 = fabs(spec[i]->z);
		F_var = (z0 == 1 ? F1 : (z0 == 2.0 ? F2 : F));
		LGAMMA[i] += z0 * z0 * F_var + z0 * CSUM;
	}
	/*
	   C
	   C     CONVERT TO MACINNES CONVENTION
	   C
	 */
	if (ICON == TRUE)
	{
		PHIMAC = LGAMMA[IC] - GAMCLM;
		/*
		   C
		   C     CORRECTED ERROR IN PHIMAC, NOVEMBER, 1989
		   C
		 */
		for (i = 0; i < 2 * count_s + count_anions; i++)
		{
			if (IPRSNT[i] == TRUE)
			{
				LGAMMA[i] = LGAMMA[i] + spec[i]->z * PHIMAC;
			}
		}
	}

	COSMOT = 1.0 + 2.0 * OSMOT / OSUM;
	/*
	   C
	   C     CALCULATE THE ACTIVITY OF WATER
	   C
	 */
	AW = exp(-OSUM * COSMOT / 55.50837);
	/*
	if (AW > 1.0)
		AW = 1.0;
	*/
	/*s_h2o->la=log10(AW); */
	mu_x = I;
	for (i = 0; i < 2 * count_s + count_anions; i++)
	{
		if (IPRSNT[i] == FALSE)
			continue;
		/*spec[i]->lg=LGAMMA[i]*CONV; */
		spec[i]->lg_pitzer = LGAMMA[i] * CONV;
		/*
		   output_msg(sformatf( "%d %s:\t%e\t%e\t%e\t%e \n", i, spec[i]->name, M[i], spec[i]->la, spec[i]->lg_pitzer, spec[i]->lg));
		 */
	}
	/*
	   output_msg(sformatf( "OSUM: %e\n", OSUM));
	   output_msg(sformatf( "OSMOT: %e\n", OSMOT));
	   output_msg(sformatf( "COSMOT: %e\n", COSMOT));
	   output_msg(sformatf( "F: %e\n", F));
	   output_msg(sformatf( "AW: %e\n", AW));
	 */
	/*
	 *I_X = I;
	 *COSMOT_X = COSMOT;
	 */
	return (OK);
}
#else
/* ---------------------------------------------------------------------- */
int Phreeqc::
pitzer(void)
/* ---------------------------------------------------------------------- */
{
	int i, i0, i1, i2;
	LDBLE param, l_alpha, z0, z1;
	LDBLE etheta, ethetap;
	/*
	   LDBLE CONV, XI, XX, OSUM, BIGZ, DI, F, XXX, GAMCLM, 
	   CSUM, PHIMAC, OSMOT, BMXP, ETHEAP, CMX, BMX, PHI,
	   BMXPHI, PHIPHI, AW, A, B;
	 */
	LDBLE CONV, XX, OSUM, BIGZ, DI, F, F1, F2, F_var, XXX, GAMCLM, CSUM, PHIMAC, OSMOT,
		B, B1, B2;
	LDBLE I, TK;
	/*
	   C
	   C     INITIALIZE
	   C
	 */
	CONV = 1.0 / log(10.0);
	XX = 0.0;
	OSUM = 0.0;
	I = mu_x;
	TK = tk_x;
	/*      DH_AB(TK, &A, &B); */
	/*
	   C
	   C     TRANSFER DATA FROM TO M
	   C
	 */
 	for (size_t j = 0; j < s_list.size(); j++)
 	{
 		i = s_list[j];
		IPRSNT[i] = FALSE;
		M[i] = 0.0;
		if (spec[i] != NULL && spec[i]->in == TRUE)
		{
			if (spec[i]->type == EX ||
				spec[i]->type == SURF || spec[i]->type == SURF_PSI)
				continue;
			M[i] = under(spec[i]->lm);
			if (M[i] > MIN_TOTAL)
				IPRSNT[i] = TRUE;
		}
	}	
	if (ICON == TRUE)
	{
		IPRSNT[IC] = TRUE;
	}
	/*
	   C
	   C     COMPUTE PITZER COEFFICIENTS' TEMPERATURE DEPENDENCE
	   C
	 */
	PTEMP(TK);
	for (size_t j = 0; j < s_list.size(); j++)
	{
		int i = s_list[j];
		LGAMMA[i] = 0.0;
		XX = XX + M[i] * fabs(spec[i]->z);
		OSUM = OSUM + M[i];
	}
	/*
	   C
	   C     EQUATION (8)
	   C
	 */
	BIGZ = XX;
	DI = sqrt(I);
	/*
	   C
	   C     CALCULATE F & GAMCLM
	   C
	 */
	B = 1.2;
	F = F1 = F2 = -A0 * (DI / (1.0 + B * DI) + 2.0 * log(1.0 + B * DI) / B);
	if (patm_x > 1.0)
	{
		LDBLE pap = 0.0;
		pap = (7e-5 + 1.93e-9 * pow(TK - 250.0, 2.0)) * patm_x;
		B1 = B - (pap > 0.2 ? 0.2 : pap);
		if (TK > 263.0)
		{
			pap = (9.65e-10 * pow(TK - 263.0, 2.773)) * pow(patm_x, 0.623);
			//pap = (-5.22e-4 + 7.19e-8 * pow(TK - 263.0, 2.0)) * pow(patm_x, 0.623);
		}
		B2 = B - (pap > 0.2 ? 0.2 : pap);
		if (B1 != 0)
			F1 = -A0 * (DI / (1.0 + B1 * DI) + 2.0 * log(1.0 + B1 * DI) / B1);
		if (B2 != 0)
			F2 = -A0 * (DI / (1.0 + B2 * DI) + 2.0 * log(1.0 + B2 * DI) / B2);
	}
	XXX = 2.0 * DI;
	XXX = (1.0 - (1.0 + XXX - XXX * XXX * 0.5) * exp(-XXX)) / (XXX * XXX);
	GAMCLM = F1;
	if (mcb0 != NULL)
		GAMCLM += I * 2.0 * mcb0->p;
	if (mcb1 != NULL)
		GAMCLM += I * 2.0 * mcb1->p * XXX;
	if (mcc0 != NULL)
		GAMCLM += 1.5 * mcc0->p * I * I;
	CSUM = 0.0;
	OSMOT = -(A0) * pow(I, (LDBLE) 1.5) / (1.0 + B * DI);
	/*
	 *  Calculate ethetas
	 */
	if (use_etheta == TRUE)
	{
		for (i = 0; i < count_theta_param; i++)
		{
			z0 = theta_params[i]->zj;
			z1 = theta_params[i]->zk;
			ETHETAS(z0, z1, I, &etheta, &ethetap);
			theta_params[i]->etheta = etheta;
			theta_params[i]->ethetap = ethetap;
		}
	}
	/*
	 *  Sums for F, LGAMMA, and OSMOT
	 */
 	for (size_t j = 0; j < param_list.size(); j++)
 	{
 		int i = param_list[j];
		i0 = pitz_params[i]->ispec[0];
		i1 = pitz_params[i]->ispec[1];
		z0 = spec[i0]->z;
		z1 = spec[i1]->z;
		param = pitz_params[i]->p;
		l_alpha = pitz_params[i]->alpha;
		F_var = 0;
		switch (pitz_params[i]->type)
		{
		case TYPE_B0:
			LGAMMA[i0] += M[i1] * 2.0 * param;
			LGAMMA[i1] += M[i0] * 2.0 * param;
			OSMOT += M[i0] * M[i1] * param;
			break;
		case TYPE_B1:
			if (param != 0.0)
			{
				F_var = M[i0] * M[i1] * param * GP(l_alpha * DI) / I;
				LGAMMA[i0] += M[i1] * 2.0 * param * G(l_alpha * DI);
				LGAMMA[i1] += M[i0] * 2.0 * param * G(l_alpha * DI);
				OSMOT += M[i0] * M[i1] * param * exp(-l_alpha * DI);
			}
			break;
		case TYPE_B2:
			if (param != 0.0)
			{
				F_var = M[i0] * M[i1] * param * GP(l_alpha * DI) / I;
				LGAMMA[i0] += M[i1] * 2.0 * param * G(l_alpha * DI);
				LGAMMA[i1] += M[i0] * 2.0 * param * G(l_alpha * DI);
				OSMOT += M[i0] * M[i1] * param * exp(-l_alpha * DI);
			}
			break;
		case TYPE_C0:
			CSUM +=
				M[i0] * M[i1] * pitz_params[i]->p / (2.0 *
													 sqrt(fabs(z0 * z1)));
			LGAMMA[i0] += M[i1] * BIGZ * param / (2.0 * sqrt(fabs(z0 * z1)));
			LGAMMA[i1] += M[i0] * BIGZ * param / (2.0 * sqrt(fabs(z0 * z1)));
			OSMOT +=
				M[i0] * M[i1] * BIGZ * param / (2.0 * sqrt(fabs(z0 * z1)));
			break;
		case TYPE_THETA:
			LGAMMA[i0] += 2.0 * M[i1] * (param /*+ ETHETA(z0, z1, I) */ );
			LGAMMA[i1] += 2.0 * M[i0] * (param /*+ ETHETA(z0, z1, I) */ );
			OSMOT += M[i0] * M[i1] * param;
			break;
		case TYPE_ETHETA:
			/*
			   ETHETAS(z0, z1, I, &etheta, &ethetap);
			 */
			if (use_etheta == TRUE)
			{
				etheta = pitz_params[i]->thetas->etheta;
				ethetap = pitz_params[i]->thetas->ethetap;
				F_var = M[i0] * M[i1] * ethetap;
				LGAMMA[i0] += 2.0 * M[i1] * etheta;
				LGAMMA[i1] += 2.0 * M[i0] * etheta;
				OSMOT += M[i0] * M[i1] * (etheta + I * ethetap);
			}
			break;
		case TYPE_PSI:
			i2 = pitz_params[i]->ispec[2];
			if (IPRSNT[i2] == FALSE)
				continue;
			LGAMMA[i0] += M[i1] * M[i2] * param;
			LGAMMA[i1] += M[i0] * M[i2] * param;
			LGAMMA[i2] += M[i0] * M[i1] * param;
			OSMOT += M[i0] * M[i1] * M[i2] * param;
			break;
		case TYPE_LAMDA:
			LGAMMA[i0] += M[i1] * param * pitz_params[i]->ln_coef[0];
			LGAMMA[i1] += M[i0] * param * pitz_params[i]->ln_coef[1];
			OSMOT += M[i0] * M[i1] * param * pitz_params[i]->os_coef;
			break;
		case TYPE_ZETA:
			i2 = pitz_params[i]->ispec[2];
			if (IPRSNT[i2] == FALSE)
				continue;
			LGAMMA[i0] += M[i1] * M[i2] * param;
			LGAMMA[i1] += M[i0] * M[i2] * param;
			LGAMMA[i2] += M[i0] * M[i1] * param;
			OSMOT += M[i0] * M[i1] * M[i2] * param;
			break;
		case TYPE_MU:
			i2 = pitz_params[i]->ispec[2];
			if (IPRSNT[i2] == FALSE)
				continue;

			LGAMMA[i0] += M[i1] * M[i2] * param * pitz_params[i]->ln_coef[0];
			LGAMMA[i1] += M[i0] * M[i2] * param * pitz_params[i]->ln_coef[1];
			LGAMMA[i2] += M[i0] * M[i1] * param * pitz_params[i]->ln_coef[2];
			OSMOT += M[i0] * M[i1] * M[i2] * param * pitz_params[i]->os_coef;
			break;
		case TYPE_ETA:
			i2 = pitz_params[i]->ispec[2];
			if (IPRSNT[i2] == FALSE)
				continue;
			LGAMMA[i0] += M[i1] * M[i2] * param;
			LGAMMA[i1] += M[i0] * M[i2] * param;
			LGAMMA[i2] += M[i0] * M[i1] * param;
			OSMOT += M[i0] * M[i1] * M[i2] * param;
			break;
		case TYPE_ALPHAS:
			break;
		case TYPE_Other:
		default:
			error_msg("TYPE_Other in pitz_param list.", STOP);
			break;
		}
		F += F_var;
		F1 += F_var;
		F2 += F_var;
	}

	/*
	 *  Add F and CSUM terms to LGAMMA
	 */
	for (size_t j = 0; j < ion_list.size(); j++)
	{
		int i = ion_list[j];
		z0 = fabs(spec[i]->z);
		F_var = (z0 == 1 ? F1 : (z0 == 2.0 ? F2 : F));
		LGAMMA[i] += z0 * z0 * F_var + z0 * CSUM;
	}
	/*
	   C
	   C     CONVERT TO MACINNES CONVENTION
	   C
	 */
	if (ICON == TRUE)
	{
		PHIMAC = LGAMMA[IC] - GAMCLM;
		/*
		   C
		   C     CORRECTED ERROR IN PHIMAC, NOVEMBER, 1989
		   C
		 */
		for (size_t j = 0; j < s_list.size(); j++)
		{
			int i = s_list[j];
			LGAMMA[i] = LGAMMA[i] + spec[i]->z * PHIMAC;
		}
	}

	COSMOT = 1.0 + 2.0 * OSMOT / OSUM;
	/*
	   C
	   C     CALCULATE THE ACTIVITY OF WATER
	   C
	 */
	AW = exp(-OSUM * COSMOT / 55.50837);
	/*
	if (AW > 1.0)
		AW = 1.0;
	*/
	/*s_h2o->la=log10(AW); */
	mu_x = I;
	for (size_t j = 0; j < s_list.size(); j++)
	{
		int i = s_list[j];
		spec[i]->lg_pitzer = LGAMMA[i] * CONV;
	}
	/*
	 *I_X = I;
	 *COSMOT_X = COSMOT;
	 */
	return (OK);
}
#endif
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
G(LDBLE L_Y)
/* ---------------------------------------------------------------------- */
{
	LDBLE d=0.0;
	if (L_Y != 0.0)
	{
		d = 2.0e0 * (1.0e0 - (1.0e0 + L_Y) * exp(-L_Y)) / (L_Y * L_Y);
	}

	return (d);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
GP(LDBLE L_Y)
/* ---------------------------------------------------------------------- */
{
	LDBLE d=0.0;
	if (L_Y != 0.0)
	{
		d = -2.0e0 * (1.0e0 - (1.0e0 + L_Y + L_Y * L_Y / 2.0e0) * exp(-L_Y)) /
			(L_Y * L_Y);
	}
	return d;
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
ETHETAS(LDBLE ZJ, LDBLE ZK, LDBLE I, LDBLE * etheta, LDBLE * ethetap)
/* ---------------------------------------------------------------------- */
{
	/* Revised ETHETAS code thanks to Wouter Falkena and the MoReS team, June, 2015 */
   *etheta = 0.0;
   *ethetap = 0.0;

   if (ZJ == ZK)
      return (OK);

   const LDBLE XCON = 6.0e0 * A0 * sqrt(I);
   const LDBLE ZZ = ZJ * ZK;
/*
C
C     NEXT 3 ARE EQUATION (A1)
C
*/
   const LDBLE XJK = XCON * ZZ;
   const LDBLE XJJ = XCON * ZJ * ZJ;
   const LDBLE XKK = XCON * ZK * ZK;

/*
C
C     EQUATION (A3)
C
*/
   LDBLE JAY_XJK;
   LDBLE JPRIME_XJK;
   ETHETA_PARAMS( XJK, JAY_XJK, JPRIME_XJK );

   LDBLE JAY_XJJ;
   LDBLE JPRIME_XJJ;
   ETHETA_PARAMS( XJJ, JAY_XJJ, JPRIME_XJJ );

   LDBLE JAY_XKK;
   LDBLE JPRIME_XKK;
   ETHETA_PARAMS( XKK, JAY_XKK, JPRIME_XKK );

   *etheta =
      ZZ * (JAY_XJK - JAY_XJJ / 2.0e0 - JAY_XKK / 2.0e0) / (4.0e0 * I);
   *ethetap =
      ZZ * (JPRIME_XJK - JPRIME_XJJ / 2.0e0 -
            JPRIME_XKK / 2.0e0) / (8.0e0 * I * I) - *etheta / I;

   return (OK);
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
ETHETA_PARAMS(LDBLE X, LDBLE& JAY, LDBLE& JPRIME )
/* ---------------------------------------------------------------------- */
/*
C
C     NUMERICAL APPROXIMATION TO THE INTEGRALS IN THE EXPRESSIONS FOR J0
C     AND J1.  CHEBYSHEV APPROXIMATION IS USED.  THE CONSTANTS 'AK' ARE
C     DEFINED IN BLOCK COMMON.
C
*/
/*
C
C     AK IS USED TO CALCULATE HIGHER ORDER ELECTROSTATIC TERMS IN
C     SUBROUTINE PITZER
C
*/
{
   static const LDBLE AKX[42] = {
      1.925154014814667e0, -.060076477753119e0, -.029779077456514e0,
      -.007299499690937e0, 0.000388260636404e0, 0.000636874599598e0,
      0.000036583601823e0, -.000045036975204e0, -.000004537895710e0,
      0.000002937706971e0, 0.000000396566462e0, -.000000202099617e0,
      -.000000025267769e0, 0.000000013522610e0, 0.000000001229405e0,
      -.000000000821969e0, -.000000000050847e0, 0.000000000046333e0,
      0.000000000001943e0, -.000000000002563e0, -.000000000010991e0,
      0.628023320520852e0, 0.462762985338493e0, 0.150044637187895e0,
      -.028796057604906e0, -.036552745910311e0, -.001668087945272e0,
      0.006519840398744e0, 0.001130378079086e0, -.000887171310131e0,
      -.000242107641309e0, 0.000087294451594e0, 0.000034682122751e0,
      -.000004583768938e0, -.000003548684306e0, -.000000250453880e0,
      0.000000216991779e0, 0.000000080779570e0, 0.000000004558555e0,
      -.000000006944757e0, -.000000002849257e0, 0.000000000237816e0
   };
/*
      LDBLE PRECISION AK, BK, DK
      COMMON / MX8 / AK(0:20,2),BK(0:22),DK(0:22)
*/
   const LDBLE *AK;
   LDBLE L_Z = 0.0;
   LDBLE L_DZ = 0.0;

   if ( X <= 1.0e0 )
   {
      const LDBLE powX0_2 = pow( X, 0.2 );
      L_Z  = 4.0e0 * powX0_2 - 2.0e0;
      L_DZ = 0.8e0 * powX0_2 / 2.0e0;
      AK = &AKX[0];
   }
   else
   {
      const LDBLE powXmin0_1 = pow( X, -0.1 );
      L_Z  = ( 40.0e0 * powXmin0_1 - 22.0e0 ) / 9.0e0;
      L_DZ = -4.0e0 * powXmin0_1 / 18.0e0;
      AK = &AKX[21];
   }

   BK[20] = AK[20];
   BK[19] = L_Z * AK[20] + AK[19];
   DK[19] = AK[20];
   for ( int i = 18; i >= 0; i-- )
   {
      BK[i] = L_Z * BK[i + 1] - BK[i + 2] + AK[i];
      DK[i] = BK[i + 1] + L_Z * DK[i + 1] - DK[i + 2];
   }

   JAY = X / 4.0e0 - 1.0e0 + 0.5e0 * (BK[0] - BK[2]);
   JPRIME = X * .25e0 + L_DZ * (DK[0] - DK[2]);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
pitzer_clean_up(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Free all allocated memory, except strings
 */
	int i;
	for (i = 0; i < count_pitz_param; i++)
	{
		pitz_params[i] =
			(struct pitz_param *) free_check_null(pitz_params[i]);
	}
	count_pitz_param = 0;
	pitz_param_map.clear();
	pitz_params = (struct pitz_param **) free_check_null(pitz_params);
	for (i = 0; i < count_theta_param; i++)
	{
		theta_params[i] =
			(struct theta_param *) free_check_null(theta_params[i]);
	}
	count_theta_param = 0;
	theta_params = (struct theta_param **) free_check_null(theta_params);
	LGAMMA = (LDBLE *) free_check_null(LGAMMA);
	IPRSNT = (int *) free_check_null(IPRSNT);
	spec = (struct species **) free_check_null(spec);
	M = (LDBLE *) free_check_null(M);

	return OK;
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
set_pz(int initial)
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
	iterations = -1;
	solution_ptr = use.Get_solution_ptr();
	for (i = 0; i < count_s_x; i++)
	{
		s_x[i]->lm = LOG_ZERO_MOLALITY;
		s_x[i]->lg_pitzer = 0.0;
	}
	if (initial == TRUE || set_and_run_attempt > 0)
	{
		for (i = 0; i < count_s_x; i++)
		{
			s_x[i]->lg = 0.0;
		}
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
	AW = pow((LDBLE) 10.0, s_h2o->la);
	s_hplus->la = -solution_ptr->Get_ph();
	s_hplus->lm = s_hplus->la;
	s_hplus->moles = exp(s_hplus->lm * LOG_10) * mass_water_aq_x;
	s_eminus->la = -solution_ptr->Get_pe();
	if (initial == TRUE)
		pitzer_initial_guesses();
	if (dl_type_x != cxxSurface::NO_DL)
		initial_surface_water();
	pitzer_revise_guesses();
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
pitzer_initial_guesses(void)
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
pitzer_revise_guesses(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Revise molalities species
 */
	int i;
	int l_iter, max_iter, repeat, fail;
	LDBLE weight, f;

	max_iter = 10;
	/* gammas(mu_x); */
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
		/*pitzer(); */
		/*s_h2o->la = 0.0; */
		/*molalities(TRUE); */
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
		/*n
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
/*			    x[i]->type == ALK || */
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
				if (f == 0 && x[i]->moles == 0)
				{
					x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
					continue;
				}
				else if (f == 0)
				{
					repeat = TRUE;
					x[i]->master[0]->s->la += 5;
/*!!!!*/ if (x[i]->master[0]->s->la < -999.)
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
						x[i]->master[0]->s->la +=
							weight * log10(fabs(x[i]->moles / x[i]->sum));
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
	log_msg(sformatf( "Iterations in pitzer_revise_guesses: %d\n", l_iter));
	/*mu_x = mu_unknown->f * 0.5 / mass_water_aq_x; */
	if (mu_x <= 1e-8)
	{
		mu_x = 1e-8;
	}
	/*gammas(mu_x); */
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
jacobian_pz(void)
/* ---------------------------------------------------------------------- */
{
	LDBLE *base;
	LDBLE d, d1, d2;
	int i, j;

Restart:
	int pz_max_unknowns = max_unknowns;
	//k_temp(tc_x, patm_x);
	if (full_pitzer == TRUE)
	{
		molalities(TRUE);
		pitzer();
		residuals();
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
		case AH2O:
			x[i]->master[0]->s->la += d;
			d2 = d1;
			break;
		case PITZER_GAMMA:
			if (!full_pitzer) 
				continue;
			x[i]->s->lg += d;
			d2 = d;
			break;
		case MH2O:
			mass_water_aq_x *= (1.0 + d);
			x[i]->master[0]->s->moles = mass_water_aq_x / gfw_water;
			d2 = log(1.0 + d);
			break;
		case MH:
			if (pitzer_pe == TRUE)
			{
				s_eminus->la += d;
				d2 = d1;
				break;
			}
			else
			{
				continue;
			}
		case GAS_MOLES:
			if (gas_in == FALSE)
				continue;
			d2 = d * x[i]->moles;
			if (d2 < 1e-14)
				d2 = 1e-14;
			x[i]->moles += d2;
			break;
		case MU:
			//continue;
			d2 = d * mu_x;
			mu_x += d2;
			//k_temp(tc_x, patm_x);
			gammas(mu_x);
			break;
		case PP:
		case SS_MOLES:
			continue;
			break;
		}
		molalities(TRUE);
		if (max_unknowns > pz_max_unknowns) 
		{
			base = (LDBLE *) free_check_null(base);
			gammas_pz();
			jacobian_sums();
			goto Restart;
		}
		if (full_pitzer == TRUE)
			pitzer();
		mb_sums();
		residuals();
		for (j = 0; j < count_unknowns; j++)
		{
			array[j * (count_unknowns + 1) + i] =
				-(residual[j] - base[j]) / d2;
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
				array[i * (count_unknowns + 1) + i] =
					exp(s_h2->lm * LOG_10) * 2;
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
			//k_temp(tc_x, patm_x);
			gammas(mu_x);
			break;
		case GAS_MOLES:
			if (gas_in == FALSE)
				continue;
			x[i]->moles -= d2;
			break;

		}
	}
	molalities(TRUE);
	if (full_pitzer == TRUE)
		pitzer();
	mb_sums();
	residuals();
	free_check_null(base);
	return OK;
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
model_pz(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   model is called after the equations have been set up by prep
 *   and initial guesses have been made in set.
 * 
 *   Here is the outline of the calculation sequence:
 *      residuals--residuals are calculated, if small we are done
 *      sum_jacobian--jacobian is calculated 
 *      ineq--inequality solver is called
 *      reset--estimates of unknowns revised, if changes are small solution
 *         has been found, usually convergence is found in residuals.
 *      gammas--new activity coefficients
 *      molalities--calculate molalities
 *      mb_sums--calculate mass-balance sums
 *      mb_gases--decide if gas_phase exists
 *      mb_ss--decide if solid_solutions exists
 *      switch_bases--check to see if new basis species is needed
 *         reprep--rewrite equations with new basis species if needed
 *         pitzer_revise_guesses--revise unknowns to get initial mole balance
 *      check_residuals--check convergence one last time
 *         sum_species--calculate sums of elements from species concentrations
 *
 *      An additional pass through may be needed if unstable phases still exist
 *         in the phase assemblage. 
 */
	int l_kode, return_kode;
	int r;
	int count_infeasible, count_basis_change;
	int debug_model_save;
	int mass_water_switch_save;

/*	debug_model = TRUE; */
/*	debug_prep = TRUE; */
/*	debug_set = TRUE; */
	/* mass_water_switch == TRUE, mass of water is constant */
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
	gamma_iterations = 0;
	count_basis_change = count_infeasible = 0;
	stop_program = FALSE;
	remove_unstable_phases = FALSE;
	if (always_full_pitzer == TRUE)
	{
		full_pitzer = TRUE;
	}
	else
	{
		full_pitzer = FALSE;
	}
#if defined(PITZER_LISTS)
	//pitzer_make_lists();
#endif
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
			gammas_pz();
			jacobian_sums();
			jacobian_pz();
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
					log_msg(sformatf( "Ineq had infeasible solution, "
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
			gammas_pz();
			if (full_pitzer == TRUE)
				pitzer();
			if (always_full_pitzer == TRUE)
			{
				full_pitzer = TRUE;
			}
			else
			{
				full_pitzer = FALSE;
			}
			molalities(TRUE);
			if (use.Get_surface_ptr() != NULL &&
				use.Get_surface_ptr()->Get_dl_type() != cxxSurface::NO_DL &&
				use.Get_surface_ptr()->Get_related_phases())
				initial_surface_water();
			mb_sums();
			mb_gases();
			mb_ss();
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
		/* remove_unstable_phases is set in check_residuals */
		if (remove_unstable_phases == FALSE && mass_water_switch_save == FALSE
			&& mass_water_switch == TRUE)
		{
			log_msg(sformatf(
					   "\nChanging water switch to FALSE. Iteration %d.\n",
					   iterations));
			mass_water_switch = FALSE;
			continue;
		}
		gamma_iterations++;
		if (gamma_iterations > itmax)
		{
			error_string = sformatf( "Maximum gamma iterations exceeded, %d\n",
					itmax);
			warning_msg(error_string);
			stop_program = TRUE;
			break;
		}
		if (check_gammas_pz() != TRUE)
		{
			full_pitzer = TRUE;
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
	log_msg(sformatf( "Number of iterations: %d\n", iterations));
	log_msg(sformatf( "Number of gamma iterations: %d\n\n", gamma_iterations));
	debug_model = debug_model_save;
	set_forward_output_to_log(FALSE);
	if (stop_program == TRUE)
	{
		return (ERROR);
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
check_gammas_pz(void)
/* ---------------------------------------------------------------------- */
{
	LDBLE /*old_aw,*/ old_mu, tol;
	int converge, i;

	old_mu = mu_x;
	pitzer();
	molalities(TRUE);
	mb_sums();
	converge = TRUE;
	tol = convergence_tolerance * 10.;
	for (i = 0; i < count_unknowns; i++)
	{
		if (x[i]->type != PITZER_GAMMA)
			continue;
		if (fabs(x[i]->s->lg - x[i]->s->lg_pitzer) > tol)
		{
			converge = FALSE;
		}
	}
	if (fabs(old_mu - mu_x) > tol)
		converge = FALSE;

	if ((pow((LDBLE) 10.0, s_h2o->la) - AW) > tol)
		converge = FALSE;
	return converge;
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
gammas_pz()
/* ---------------------------------------------------------------------- */
{
/*
 *   Need exchange gammas for pitzer
 */
	int i, j;
	LDBLE coef;
	/* Initialize */
	k_temp(tc_x, patm_x);
/*
 *   Calculate activity coefficients
 */
	for (i = 0; i < count_s_x; i++)
	{
		switch (s_x[i]->gflag)
		{
		case 0:				/* uncharged */
		case 1:				/* Davies */
		case 2:				/* Extended D-H, WATEQ D-H */
		case 3:				/* Always 1.0 */
			break;
		case 4:				/* Exchange */
			/* Now calculated in next loop */
			break;
		case 5:				/* Always 1.0 */
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
			if (s_x[i]->alk > 0)
			{
				s_x[i]->lg = log10(s_x[i]->equiv / s_x[i]->alk);
				s_x[i]->dg = 0.0;
			}
			else
			{
				s_x[i]->lg = 0.0;
				s_x[i]->dg = 0.0;
			}
			break;
		case 7:				/* LLNL */
			break;
		case 8:				/* LLNL CO2 */
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
	/*
	 *  calculate exchange gammas 
	 */

	if (use.Get_exchange_ptr() != NULL)
	{
		for (i = 0; i < count_s_x; i++)
		{
			switch (s_x[i]->gflag)
			{
			case 0:			/* uncharged */
			case 1:			/* Davies */
			case 2:			/* Extended D-H, WATEQ D-H */
			case 3:			/* Always 1.0 */
			case 5:			/* Always 1.0 */
			case 6:			/* Surface */
			case 7:			/* LLNL */
			case 8:			/* LLNL CO2 */
			case 9:			/* activity water */
				break;
			case 4:			/* Exchange */

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
							s_x[i]->rxn_x->token[j].s->primary->unknown->
							moles;
						break;
					}
				}
				/*
				 *   Master species is a dummy variable with meaningless activity and mass
				 */
				s_x[i]->lg = 0.0;
				s_x[i]->dg = 0.0;
				if (s_x[i]->primary != NULL)
				{
					break;
				}
				/*
				 *   All other species
				 */

				/* modific 29 july 2005... */
				if (s_x[i]->equiv != 0 && s_x[i]->alk > 0)
				{
					s_x[i]->lg = log10(fabs(s_x[i]->equiv) / s_x[i]->alk);
				}
				if (use.Get_exchange_ptr()->Get_pitzer_exchange_gammas())
				{
					/* Assume equal gamma's of solute and exchangeable species...  */
					for (j = 1; s_x[i]->rxn_x->token[j].s != NULL; j++)
					{
						if (s_x[i]->rxn_x->token[j].s->type == EX)
							continue;
						coef = s_x[i]->rxn_x->token[j].coef;
						s_x[i]->lg += coef * s_x[i]->rxn_x->token[j].s->lg;
						s_x[i]->dg += coef * s_x[i]->rxn_x->token[j].s->dg;
					}
				}
			}
		}
	}
/* ...end modific 29 july 2005 */

	return (OK);
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
pitzer_make_lists(void)
/* ---------------------------------------------------------------------- */
{
	double log_min = log10(MIN_TOTAL);
	s_list.clear();
	cation_list.clear();
	neutral_list.clear();
	anion_list.clear();
	ion_list.clear();
	param_list.clear();
	OTEMP = -100.0;
	for (int j = 0; j < 3; j++)
	{
		int min, max;
		switch (j)
		{
		case 0:
			min = 0;
			max = count_cations;
			break;
		case 1:
			min = count_s;
			max = count_s + count_neutrals;
			break;
		case 2:
			min = 2*count_s;
			max = 2*count_s + count_anions;
			break;
		}
		for (int i = min; i < max; i++)
		{
			IPRSNT[i] = FALSE;
			M[i] = 0.0;
		if ((spec[i] != NULL && spec[i]->in == TRUE) ||
			(ICON == TRUE && i == IC))
			{
				if (spec[i]->type == EX ||
					spec[i]->type == SURF || spec[i]->type == SURF_PSI)
					continue;	
				IPRSNT[i] = TRUE;	
				s_list.push_back(i);	
				if (i < count_s)
				{
					cation_list.push_back(i);
				}
				if (i >= count_s && i < 2*count_s)
				{
					neutral_list.push_back(i);
				}
				if (i >= 2*count_s)
				{
					anion_list.push_back(i);
				}
				if (i < count_s || i >= 2*count_s)
				{
					ion_list.push_back(i);
				}
				if (spec[i]->lm > log_min)
				{
					M[i] = under(spec[i]->lm);
				}
			}
		}
	}	
	if (ICON == TRUE)
	{
		IPRSNT[IC] = TRUE;
	}
	for (int i = 0; i < count_pitz_param; i++)
	{
		/*
		TYPE_B0, TYPE_B1, TYPE_B2, TYPE_C0, TYPE_THETA, TYPE_LAMDA, TYPE_ZETA,
		TYPE_PSI, TYPE_ETHETA, TYPE_ALPHAS, TYPE_MU, TYPE_ETA, TYPE_Other,
		TYPE_SIT_EPSILON, TYPE_SIT_EPSILON_MU
		*/
		int i0 = pitz_params[i]->ispec[0];
		int i1 = pitz_params[i]->ispec[1];
		if (IPRSNT[i0] == FALSE || IPRSNT[i1] == FALSE) continue;
		int i2 = pitz_params[i]->ispec[2];
		if (pitz_params[i]->type == TYPE_PSI ||
			pitz_params[i]->type == TYPE_ZETA ||
			pitz_params[i]->type == TYPE_MU ||
			pitz_params[i]->type == TYPE_ETA)
		{
			if (IPRSNT[i2] == FALSE)
				continue;
		}
		param_list.push_back(i);
	}
}