#include "Phreeqc.h"
#include "phqalloc.h"

#include "Utils.h"
#include "NameDouble.h"
#include "PBasic.h"
#include "Exchange.h"
#include "GasPhase.h"
#include "PPassemblage.h"
#include "SSassemblage.h"
#include "cxxKinetics.h"
#include "Solution.h"
#include "Parser.h"
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
activity(const char *species_name)
/* ---------------------------------------------------------------------- */
{
	struct species *s_ptr;
	LDBLE a;

	s_ptr = s_search(species_name);
	if (s_ptr == s_h2o)
	{
		a = pow((LDBLE) 10., s_h2o->la);
	}
	else if (s_ptr == s_eminus)
	{
		a = pow((LDBLE) 10., s_eminus->la);
	}
	else if (s_ptr == NULL || s_ptr->in == FALSE)
	{
		a = 1e-99;
	}
	else
	{
		a = pow((LDBLE) 10., s_ptr->lm + s_ptr->lg);
	}
	return (a);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
activity_coefficient(const char *species_name)
/* ---------------------------------------------------------------------- */
{
	struct species *s_ptr;
	LDBLE g;

	s_ptr = s_search(species_name);
	if (s_ptr != NULL && s_ptr->in != FALSE && ((s_ptr->type < EMINUS) || (s_ptr->type == EX) || (s_ptr->type == SURF)))
	{
		g = pow((LDBLE) 10., s_ptr->lg);
	}
	else
	{
		g = 0;
	}
	return (g);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
log_activity_coefficient(const char *species_name)
/* ---------------------------------------------------------------------- */
{
	struct species *s_ptr;
	LDBLE g;

	s_ptr = s_search(species_name);
	if (s_ptr != NULL && s_ptr->in != FALSE && ((s_ptr->type < EMINUS) || (s_ptr->type == EX) || (s_ptr->type == SURF)))
	{
		g = s_ptr->lg;
	}
	else
	{
		g = 0;
	}
	return (g);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
aqueous_vm(const char *species_name)
/* ---------------------------------------------------------------------- */
{
	struct species *s_ptr;
	LDBLE g;

	s_ptr = s_search(species_name);
	if (s_ptr != NULL && s_ptr->in != FALSE && s_ptr->type < EMINUS)
	{
		g = s_ptr->logk[vm_tc];
	}
	else
	{
		g = 0;
	}
	return (g);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
sa_declercq(double sa_type, double Sa, double d, double m, double m0, double gfw)
/* ---------------------------------------------------------------------- */
{
	if (sa_type == 0) 
	{
		// surface-area-calculation-Fixed_Surface
		return Sa;
	}
	else if (sa_type == 1)
		// surface-area-calculation-Square
	{
		double mass0 = m0 * gfw;
		double V0 = mass0 / d;
		double St0 = mass0 * Sa;        // total surface
		double a0 = pow(V0, 1.0/3.0);   // side length
		double Sp0 = 6.0 * a0*a0;       // surface particle
		double np = St0 / Sp0;          // number of particles
	    double RATS = Sa / St0;
		double mass = m * gfw;
		double V = mass / d;
		double a = pow(V, 1.0/3.0); 
		double St = 6.0 * a*a*np;
		return St * RATS;               // total current surface
	}
	else if (sa_type == 2)
	{
		//double pi = 3.14159265359;
		double mass0 = m0 * gfw;
		double V0 = mass0 / d;                         // volume
		double St0 = mass0 * Sa;                       // total surface
		double a0 = pow(3.0 * V0/(4.0 * pi), 1.0/3.0); // ((3*V0)/(4 * 3.14159265359))^(1/3)  
		double Sp0 = (4.0 * pi) * a0 * a0;             // surface particle
		double np = St0 / Sp0;                         // number of particles
		double RATS = Sa / St0;
 
		double mass = m * gfw;
		double V = mass / d;
		double a = pow(3.0 * V/(4.0 * pi), 1.0/3.0);  //((3*V)/(4 * 3.14159265359))^(1/3)
		double St = 4.0 * pi * a * a * np;
		return St * RATS;                             // total current surface
	}
	error_string = sformatf( "Unknown surface area type in SA_DECLERCQ %d.", (int) sa_type);
	error_msg(error_string, CONTINUE);
	input_error++;
	return (MISSING);
 
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
diff_c(const char *species_name)
/* ---------------------------------------------------------------------- */
{
	struct species *s_ptr;
	LDBLE g;

	s_ptr = s_search(species_name);
	if (s_ptr != NULL && s_ptr->in != FALSE && s_ptr->type < EMINUS)
	{
		g = s_ptr->dw;
	}
	else
	{
		g = 0;
	}
	return (g);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
calc_SC(void)
/* ---------------------------------------------------------------------- */
{
	int i;
	LDBLE lm, a, l_z, Dw, SC, ff;

	SC = 0;
	for (i = 0; i < count_species_list; i++)
	{
		if (species_list[i].s->type == EX)
			continue;
		if (species_list[i].s->type == SURF)
			continue;
		if (i > 0
			&& strcmp(species_list[i].s->name,
					  species_list[i - 1].s->name) == 0)
			continue;
		if (species_list[i].s == s_h2o)
			continue;
		if ((Dw = species_list[i].s->dw) == 0)
			continue;
		if ((l_z = fabs(species_list[i].s->z)) == 0)
			continue;

		lm = species_list[i].s->lm;
		if (lm > -9)
		{
/*
      if (l_z < 1.5) {
	ff = (mu_x < 0.36 ? 0.6 :
	sqrt(mu_x));
      }
      else {
	ff = (mu_x < pow(0.4*l_z, 2.0) ? 0.4 :
	sqrt(mu_x) / l_z);
      }
*/
			ff = (mu_x < .36 * l_z ? 0.6 / sqrt(l_z) : sqrt(mu_x) / l_z);

			ff *= species_list[i].s->lg;
			if (ff > 0) ff = 0;
			a = under(lm + ff);
			SC += a * l_z * l_z * Dw;
		}
	}
	SC *= 1e7 * F_C_MOL * F_C_MOL / (R_KJ_DEG_MOL * 298160.0);
/* correct for temperature dependency...
       SC_T = SC_298 * (Dw_T / T) * (298 / Dw_298) and
	     Dw_T = Dw_298 * (T / 298) * (viscos_298 / viscos_T) give:
	 SC_T = SC_298 * (viscos_298 / viscos_T)
 */
	SC *= 0.88862 / viscosity();

	return (SC);
}

/* VP: Density Start */
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
calc_dens(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Calculates density based on the formulas and parameters from Millero,
 *   2000.
 *
 *   Millero, F.J. (2000) The equation of state of lakes. Aquatic geochemistry
 *   volume 6, pages 1-17
 */
	int i;
	LDBLE M_T, /*rho_new,*/ gfw;
	/* 2 options: original VP, assign the volumes of species with zero molar volume to their master species,
	              but this doubles counts of complexes with -Vm defined. And, cation-OH and H-anion
				  complexes are counted once. Also, must add H+ and OH-... */
	//struct species *s_ptr;

	//V_solutes = M_T = 0.0;
	//for (i = 0; i < count_species_list; i++)
	//{
	//	if (species_list[i].s->type != AQ && species_list[i].s->type != HPLUS)
	//	  continue;

	//	//if (species_list[i].master_s->secondary != NULL)
	//	//	gfw = species_list[i].master_s->secondary->gfw;
	//	//else
	//	//	gfw = species_list[i].master_s->primary->gfw;

	//	/* OH-... */
	//	if (!strcmp(species_list[i].s->name, "OH-"))
	//		s_ptr = s_search("OH-");
	//	else if (species_list[i].s->logk[vm_tc] == 0)
	//	{
	//		s_ptr = species_list[i].master_s;
	//		if (s_ptr->secondary)
	//			gfw = s_ptr->secondary->gfw;
	//		else
	//			gfw = s_ptr->primary->gfw;
	//	}
	//	else
	//	{
	//		s_ptr = species_list[i].s;
	//		compute_gfw(species_list[i].s->name, &gfw);
	//	}

	//	/* Special case: CO3-2 species */
	//	if (strcmp(s_ptr->name, "CO3-2") == 0)
	//	{
	//		if (strstr(species_list[i].s->name, "HCO3") != NULL)
	//		{
	//			s_ptr = s_search("HCO3-");
	//		} else if (strstr(species_list[i].s->name, "CO3") != NULL)
	//		{
	//			compute_gfw("CO3-2", &gfw);
	//		}
	//	}
	//	if (!gfw || !strcmp(species_list[i].s->name, "CO2")) /* CO2, H+ and OH- */
	//		compute_gfw(species_list[i].s->name, &gfw);

	//	M_T += (species_list[i].s->moles * gfw);
	//	V_solutes += species_list[i].s->moles * s_ptr->logk[vm_tc];
	//}

	/* 2nd option, use species_x, vm = 0 for complexes with undefined volume... */
	V_solutes = M_T = 0.0;
	for (i = 0; i < count_s_x; i++)
	{
		if (s_x[i]->type != AQ && s_x[i]->type != HPLUS)
		  continue;

		//compute_gfw(s_x[i]->name, &gfw);
		gfw = s_x[i]->gfw;
		M_T += s_x[i]->moles * gfw;
		V_solutes += s_x[i]->moles * s_x[i]->logk[vm_tc];
	}
	/* If pure water then return rho_0 */
	if (M_T == 0)
		return rho_0;
	else
		return rho_0 * (1e3 + M_T / mass_water_aq_x) / (rho_0 * V_solutes / mass_water_aq_x + 1e3);

	//M_T /= 1e3;
	//solution_mass =  mass_water_aq_x + M_T;
	//V_solutes = M_T - rho_0 * V_solutes / 1e3;

	//rho_new = halve(f_rho, 0.5, 2.0, 1e-7);
	//if (!PHR_ISFINITE(rho_new) || rho_new > 1.99999) rho_new = 1.99999;

	//return rho_new;
}
/* VP: Density End */
/* DP: Function for interval halving */

LDBLE Phreeqc::
f_rho(LDBLE rho_old, void *cookie)
/* ---------------------------------------------------------------------- */
{
	LDBLE rho = 1.0;
	Phreeqc * pThis;

	pThis = (Phreeqc *) cookie;

	pThis->solution_volume = pThis->solution_mass / rho_old;
	if (pThis->solution_volume != 0)
	{
		rho = pThis->V_solutes / pThis->solution_volume;
	}
	rho = rho + pThis->rho_0;
	return (rho - rho_old);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
calc_solution_volume(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Calculates solution volume based on sum of mass of element plus density
 */
	LDBLE total_mass = 0;
	LDBLE gfw;
	//compute_gfw("H", &gfw);
	gfw = s_hplus->primary->gfw;
	total_mass = total_h_x * gfw;
	//compute_gfw("O", &gfw);
	gfw = s_h2o->primary->gfw;
	total_mass += total_o_x * gfw;

	for (int i = 0; i < count_master; i++)
	{
		if (master[i]->s->type != AQ) continue;
		struct master *master_ptr = master[i];
		if (master_ptr->primary == TRUE && strcmp(master_ptr->elt->name, "Alkalinity"))
		{
			total_mass += master_ptr->total_primary * master_ptr->elt->gfw; 
		}
	}
	LDBLE rho = calc_dens();
	LDBLE vol = 1e-3 * total_mass / rho;
	return (vol);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
calc_logk_n(const char *name)
/* ---------------------------------------------------------------------- */
{
	char token[MAX_LENGTH];
	int i;
	LDBLE lk;
	struct logk *logk_ptr;
	LDBLE l_logk[MAX_LOG_K_INDICES];
	struct name_coef add_logk;

	for (i = 0; i < MAX_LOG_K_INDICES; i++)
	{
		l_logk[i] = 0.0;
	}
	strcpy(token, name);
	logk_ptr = logk_search(token);
	if (logk_ptr != NULL)
	{
		add_logk.name = token;
		add_logk.coef = 1.0;
		add_other_logk(l_logk, 1, &add_logk);
	lk = k_calc(l_logk, tk_x, patm_x * PASCAL_PER_ATM);
		return (lk);
	}
	return (-999.99);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
calc_logk_p(const char *name)
/* ---------------------------------------------------------------------- */
{
	int i, j;
	char token[MAX_LENGTH];
	struct phase *phase_ptr;
	LDBLE lk=-999.9;
	LDBLE l_logk[MAX_LOG_K_INDICES];

	strcpy(token, name);
	phase_ptr = phase_bsearch(token, &j, FALSE);

	if (phase_ptr != NULL)
	{		
		struct reaction *reaction_ptr;
		if (phase_ptr->replaced)
			reaction_ptr = phase_ptr->rxn_s;
		else
			reaction_ptr = phase_ptr->rxn;
		/*
		*   Print saturation index
		*/
		reaction_ptr->logk[delta_v] = calc_delta_v(reaction_ptr, true) -
			phase_ptr->logk[vm0];
		if (reaction_ptr->logk[delta_v])
			mu_terms_in_logk = true;
		for (i = 0; i < MAX_LOG_K_INDICES; i++)
		{
			l_logk[i] = 0.0;
		}
		//lk = k_calc(reaction_ptr->logk, tk_x, patm_x * PASCAL_PER_ATM);
		select_log_k_expression(reaction_ptr->logk, l_logk);
		add_other_logk(l_logk, phase_ptr->count_add_logk, phase_ptr->add_logk); 
		lk = k_calc(l_logk, tk_x, patm_x * PASCAL_PER_ATM);
	}
	return (lk);
}
#ifdef SKIP
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
calc_logk_s(const char *name)
/* ---------------------------------------------------------------------- */
{
	int i;
	char token[MAX_LENGTH];
	struct species *s_ptr;
	LDBLE lk, l_logk[MAX_LOG_K_INDICES];

	strcpy(token, name);
	s_ptr = s_search(token);
	if (s_ptr != NULL)
	{
		for (i = 0; i < MAX_LOG_K_INDICES; i++)
		{
			l_logk[i] = 0.0;
		}
		//if (s_ptr->moles)
			//select_log_k_expression(s_ptr->rxn_x->logk, l_logk);
		    select_log_k_expression(s_ptr->rxn->logk, l_logk);
		//{
			// perhaps calculate species' delta_v if absent?
		//	select_log_k_expression(s_ptr->rxn_s->logk, l_logk);
		//}
		add_other_logk(l_logk, s_ptr->count_add_logk, s_ptr->add_logk);
		lk = k_calc(l_logk, tk_x, patm_x * PASCAL_PER_ATM);
		return (lk);
	}
	return (-999.99);
}
#endif
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
calc_logk_s(const char *name)
/* ---------------------------------------------------------------------- */
{
	int i;
	char token[MAX_LENGTH];
	struct species *s_ptr;
	LDBLE lk, l_logk[MAX_LOG_K_INDICES];

	strcpy(token, name);
	s_ptr = s_search(token);
	if (s_ptr != NULL)
	{
		//if (s_ptr->logk[vm_tc])
		/* calculate delta_v for the reaction... */
			s_ptr->logk[delta_v] = calc_delta_v(s_ptr->rxn, false);
		for (i = 0; i < MAX_LOG_K_INDICES; i++)
		{
			l_logk[i] = 0.0;
		}
		select_log_k_expression(s_ptr->logk, l_logk);
		mu_terms_in_logk = true;
		add_other_logk(l_logk, s_ptr->count_add_logk, s_ptr->add_logk);
		lk = k_calc(l_logk, tk_x, patm_x * PASCAL_PER_ATM);
		return (lk);
	}
	return (-999.99);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
calc_surface_charge(const char *surface_name)
/* ---------------------------------------------------------------------- */
{
	char token[MAX_LENGTH], token1[MAX_LENGTH];
	char *ptr;
	int i, j, k;
	LDBLE charge;
	struct rxn_token_temp *token_ptr;
	struct master *master_ptr;

	/*
	 *   Go through species, sum charge
	 */
	charge = 0;
	for (k = 0; k < count_s_x; k++)
	{
		if (s_x[k]->type != SURF)
			continue;
		/*
		 *   Match surface_name
		 */
		count_trxn = 0;
		trxn_add(s_x[k]->rxn_s, 1.0, FALSE);	/* rxn_s is set in tidy_model */
		for (i = 1; i < count_trxn; i++)
		{
			token_ptr = &(trxn.token[i]);
			if (token_ptr->s->type != SURF)
				continue;
			master_ptr = trxn.token[i].s->primary;
			strcpy(token, master_ptr->elt->name);
			replace("_", " ", token);
			ptr = token;
			copy_token(token1, &ptr, &j);
			if (strcmp(surface_name, token1) == 0)
			{
				charge += s_x[k]->moles * s_x[k]->z;
			}
		}
	}
	return (charge);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
diff_layer_total(const char *total_name, const char *surface_name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in DDL layer
 */
	cxxSurfaceCharge *surface_charge_ptr1;
	std::string name, token, surface_name_local;
	struct master *master_ptr;

	LDBLE mass_water_surface;
	LDBLE molality, moles_excess, moles_surface, charge;

	if (use.Get_surface_ptr() == NULL || (dl_type_x == cxxSurface::NO_DL &&
									strcmp_nocase("psi", total_name) != 0 &&
									strcmp_nocase("psi1", total_name) != 0 &&
									strcmp_nocase("psi2", total_name) != 0 &&
									strcmp_nocase("charge", total_name) != 0
									&& strcmp_nocase("charge1",
													 total_name) != 0
									&& strcmp_nocase("charge2",
													 total_name) != 0
									&& strcmp_nocase("sigma",
													 total_name) != 0
									&& strcmp_nocase("sigma1",
													 total_name) != 0
									&& strcmp_nocase("sigma2",
													 total_name) != 0))
		return (0);

/*
 *   Find surface...
 */
	int j;
	for (j = 0; j < count_unknowns; j++)
	{
		if (use.Get_surface_ptr()->Get_type() == cxxSurface::DDL || use.Get_surface_ptr()->Get_type() == cxxSurface::CCM)
		{
			if (x[j]->type != SURFACE_CB)
				continue;
			name = x[j]->master[0]->elt->name;
			Utilities::replace("_psi", "", name);
		}
		else if (use.Get_surface_ptr()->Get_type() == cxxSurface::CD_MUSIC)
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
			token =  x[j]->master[0]->elt->name;
			Utilities::replace("_", " ", token);
			std::string::iterator b = token.begin();
			std::string::iterator e = token.end();
			CParser::copy_token(name, b, e);
		}
		if (surface_name != NULL)
		{
			if (strcmp(name.c_str(), surface_name) == 0)
				break;
		}
		else
		{
			break;
		}
	}
	if (j >= count_unknowns)
		return (0);
	surface_name_local = name;
	/*
	 *   Psi, charge, sigma
	 */
	if (strcmp_nocase("psi", total_name) == 0)
	{
		if (use.Get_surface_ptr()->Get_type() == cxxSurface::DDL || use.Get_surface_ptr()->Get_type() == cxxSurface::CCM)
		{
			return ((LDBLE) (x[j]->master[0]->s->la * 2 * R_KJ_DEG_MOL *
							 tk_x * LOG_10 / F_KJ_V_EQ));
		}
		else if (use.Get_surface_ptr()->Get_type() == cxxSurface::CD_MUSIC)
		{
			master_ptr = surface_get_psi_master(surface_name, SURF_PSI);
			if (master_ptr != NULL)
			{
				return ((LDBLE)
						(-master_ptr->s->la * R_KJ_DEG_MOL * tk_x * LOG_10 /
						 F_KJ_V_EQ));
			}
			else
			{
				return (0.0);
			}
		}
		else
		{
			return (0);
		}
	}
	else if (strcmp_nocase("psi1", total_name) == 0)
	{
		master_ptr = surface_get_psi_master(surface_name, SURF_PSI1);
		if (master_ptr != NULL)
		{
			return ((LDBLE)
					(-master_ptr->s->la * R_KJ_DEG_MOL * tk_x * LOG_10 /
					 F_KJ_V_EQ));
		}
		else
		{
			return (0.0);
		}
	}
	else if (strcmp_nocase("psi2", total_name) == 0)
	{
		master_ptr = surface_get_psi_master(surface_name, SURF_PSI2);
		if (master_ptr != NULL)
		{
			return ((LDBLE)
					(-master_ptr->s->la * R_KJ_DEG_MOL * tk_x * LOG_10 /
					 F_KJ_V_EQ));
		}
		else
		{
			return (0.0);
		}
	}
	else if (strcmp_nocase("charge", total_name) == 0)
	{
		if ((use.Get_surface_ptr()->Get_type() == cxxSurface::DDL || use.Get_surface_ptr()->Get_type() == cxxSurface::CCM) && dl_type_x == cxxSurface::NO_DL)
		{
			return ((LDBLE) (x[j]->f));
		}
		else if (use.Get_surface_ptr()->Get_type() == cxxSurface::CD_MUSIC)
		{
			cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[j]->surface_charge);
			return ((charge_ptr->Get_sigma0() *
					(charge_ptr->Get_specific_area() *
					 charge_ptr->Get_grams()) / F_C_MOL));
		}
		else
		{
			return (calc_surface_charge(surface_name_local.c_str()));
		}
	}
	else if (strcmp_nocase("charge1", total_name) == 0)
	{
		if (use.Get_surface_ptr()->Get_type() == cxxSurface::CD_MUSIC)
		{
			cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[j]->surface_charge);
			return ((charge_ptr->Get_sigma1() *
					(charge_ptr->Get_specific_area() *
					 charge_ptr->Get_grams()) / F_C_MOL));
		}
		else
		{
			return (0);
		}
	}
	else if (strcmp_nocase("charge2", total_name) == 0)
	{
		if (use.Get_surface_ptr()->Get_type() == cxxSurface::CD_MUSIC)
		{
			cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[j]->surface_charge);
			return ((charge_ptr->Get_sigma2() *
					(charge_ptr->Get_specific_area() *
					 charge_ptr->Get_grams()) / F_C_MOL));
		}
		else
		{
			return (0);
		}
	}
	else if (strcmp_nocase("sigma", total_name) == 0)
	{
		if (use.Get_surface_ptr()->Get_type() == cxxSurface::DDL || use.Get_surface_ptr()->Get_type() == cxxSurface::CCM)
		{
			cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[j]->surface_charge);
			if (dl_type_x != cxxSurface::NO_DL)
			{
				charge = calc_surface_charge(surface_name_local.c_str());
			}
			else
			{
				charge = x[j]->f;
			}
			if ((charge_ptr->Get_specific_area() *
				 charge_ptr->Get_grams()) > 0)
			{
				return ((charge * F_C_MOL /
						 (charge_ptr->Get_specific_area() *
						  charge_ptr->Get_grams())));
			}
			else
			{
				return (0);
			}
		}
		else if (use.Get_surface_ptr()->Get_type() == cxxSurface::CD_MUSIC)
		{
			cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[j]->surface_charge);
			return ((LDBLE) (charge_ptr->Get_sigma0()));
		}
		else
		{
			return (0);
		}
	}
	else if (strcmp_nocase("sigma1", total_name) == 0)
	{
		if (use.Get_surface_ptr()->Get_type() == cxxSurface::CD_MUSIC)
		{
			cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[j]->surface_charge);
			return ((LDBLE) (charge_ptr->Get_sigma1()));
		}
		else
		{
			return (0);
		}
	}
	else if (strcmp_nocase("sigma2", total_name) == 0)
	{
		if (use.Get_surface_ptr()->Get_type() == cxxSurface::CD_MUSIC)
		{
			cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[j]->surface_charge);
			return ((LDBLE) (charge_ptr->Get_sigma2()));
		}
		else
		{
			return (0);
		}
	}
	else if (strcmp_nocase("water", total_name) == 0)
	{
		if (dl_type_x != cxxSurface::NO_DL)
		{
			cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[j]->surface_charge);
			return (charge_ptr->Get_mass_water());
		}
		else
		{
			return (0);
		}
	}
/*
 *   find total moles of each element in diffuse layer...
 */
	surface_charge_ptr1 = use.Get_surface_ptr()->Find_charge(x[j]->surface_charge);
	if (surface_charge_ptr1)
	{
		mass_water_surface = surface_charge_ptr1->Get_mass_water();
		count_elts = 0;
		paren_count = 0;
		for (j = 0; j < count_s_x; j++)
		{
			if (s_x[j]->type > HPLUS)
				continue;
			molality = under(s_x[j]->lm);
			LDBLE g = surface_charge_ptr1->Get_g_map()[s_x[j]->z].Get_g();

			moles_excess = mass_water_aq_x * molality * (g * s_x[j]->erm_ddl +
										  mass_water_surface /
										  mass_water_aq_x * (s_x[j]->erm_ddl - 1));
			moles_surface = mass_water_surface * molality + moles_excess;
/*
 *   Accumulate elements in diffuse layer
 */
			add_elt_list(s_x[j]->next_elt, moles_surface);
		}
		if (count_elts > 0)
		{
			qsort(elt_list, (size_t) count_elts,
				  (size_t) sizeof(struct elt_list), Phreeqc:: elt_list_compare);
			elt_list_combine();
		}
/*
 *   Return totals
 */
		for (j = 0; j < count_elts; j++)
		{
			if (strcmp(elt_list[j].elt->name, total_name) == 0)
			{
				return ((LDBLE) elt_list[j].coef);
			}
		}
	}
	return (0);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
equi_phase(const char *phase_name)
/* ---------------------------------------------------------------------- */
{
	int j;

	if (use.Get_pp_assemblage_in() == FALSE || use.Get_pp_assemblage_ptr() == NULL)
		return (0);
	for (j = 0; j < count_unknowns; j++)
	{
		if (x[j]->type != PP)
			continue;
		if (strcmp_nocase(x[j]->pp_assemblage_comp_name, phase_name) == 0)
		{
			break;
		}
	}
/*
 *   Print pure phase assemblage data
 */
	cxxPPassemblage * pp_assemblage_ptr = use.Get_pp_assemblage_ptr();
	if (j == count_unknowns)
	{
		/* if not an unknown */
		std::map<std::string, cxxPPassemblageComp>::iterator it;
		it =  pp_assemblage_ptr->Get_pp_assemblage_comps().begin();
		for ( ; it != pp_assemblage_ptr->Get_pp_assemblage_comps().end(); it++)
		{
			if (strcmp_nocase
				(it->second.Get_name().c_str(), phase_name) == 0)
			{
				return (it->second.Get_moles());
			}
		}
	}
	else
	{
		/* if an unknown */
		if (x[j]->moles < 0.0)
			x[j]->moles = 0.0;
		return (x[j]->moles);
	}
	return (0);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
equi_phase_delta(const char *phase_name)
/* ---------------------------------------------------------------------- */
{
	int j;

	if (use.Get_pp_assemblage_in() == FALSE || use.Get_pp_assemblage_ptr() == NULL)
		return (0);
	for (j = 0; j < count_unknowns; j++)
	{
		if (x[j]->type != PP)
			continue;
		if (strcmp_nocase(x[j]->pp_assemblage_comp_name, phase_name) == 0)
		{
			break;
		}
	}
/*
 *   Print pure phase assemblage data
 */
	cxxPPassemblage * pp_assemblage_ptr = use.Get_pp_assemblage_ptr();
	if (j == count_unknowns)
	{
		/* if not an unknown */
		std::map<std::string, cxxPPassemblageComp>::iterator it;
		it =  pp_assemblage_ptr->Get_pp_assemblage_comps().begin();
		for ( ; it != pp_assemblage_ptr->Get_pp_assemblage_comps().end(); it++)
		{
			if (strcmp_nocase
				(it->second.Get_name().c_str(), phase_name) == 0)
			{
				cxxPPassemblageComp * comp_ptr = &(it->second);
				if (state != TRANSPORT && state != PHAST)
				{
					//LDBLE moles = it->second.Get_moles();
					//LDBLE delta_moles = moles - comp_ptr->Get_moles() -
					//	comp_ptr->Get_delta();
					return(0);
				}
				else
				{
					LDBLE moles = it->second.Get_moles();
					LDBLE delta_moles =	moles - comp_ptr->Get_initial_moles();
					return(delta_moles);
				}
			}
		}
	}
	else
	{
		//cxxPPassemblageComp * comp_ptr = pp_assemblage_ptr->Find(x[j]->pp_assemblage_comp_name);
		cxxPPassemblageComp * comp_ptr = (cxxPPassemblageComp *) x[j]->pp_assemblage_comp_ptr;
		if (state != TRANSPORT && state != PHAST)
		{
			LDBLE delta_moles =
				x[j]->moles - comp_ptr->Get_moles() -
				comp_ptr->Get_delta();
			return(delta_moles);
		}
		else
		{
			LDBLE delta_moles =
				x[j]->moles - comp_ptr->Get_initial_moles();
			return(delta_moles);
		}
	}
	return (0);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
equivalent_fraction(const char *name, LDBLE *eq, std::string &elt_name)
/* ---------------------------------------------------------------------- */
{
	struct species *s_ptr = s_search(name);
	*eq = 0;
	elt_name.clear();
	LDBLE f = 0;
	if (s_ptr != NULL && (s_ptr->type == EX || s_ptr->type == SURF))
	{
		*eq = s_ptr->equiv;
		struct elt_list *next_elt;
		LDBLE tot=0.0;
		for (next_elt = s_ptr->next_elt; next_elt->elt != NULL;	next_elt++)
		{
			if (next_elt->elt->master->s->type == SURF ||
				next_elt->elt->master->s->type == EX)
			{
				tot = total_mole(next_elt->elt->name);
				elt_name = next_elt->elt->name;
			}
		}
		if (s_ptr->in == TRUE && tot > 0.0)
		{
			f = s_ptr->moles * s_ptr->equiv / tot; 
		}
	}
	return f;
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
find_gas_comp(const char *gas_comp_name)
/* ---------------------------------------------------------------------- */
{
	int i;

	if (use.Get_gas_phase_in() == FALSE || use.Get_gas_phase_ptr() == NULL)
		return (0);
	cxxGasPhase * gas_phase_ptr = use.Get_gas_phase_ptr();
	for (size_t j = 0; j < gas_phase_ptr->Get_gas_comps().size(); j++)
	{
		if (strcmp_nocase(gas_phase_ptr->Get_gas_comps()[j].Get_phase_name().c_str(), gas_comp_name) == 0)
		{
			struct phase *phase_ptr = phase_bsearch(gas_comp_name, &i, false);
			if (phase_ptr)
			{
				return (phase_ptr->moles_x);
			}
		}
	}
	return (0);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
find_gas_p(void)
/* ---------------------------------------------------------------------- */
{
	if (use.Get_gas_phase_in() == FALSE || use.Get_gas_phase_ptr() == NULL)
		return (0);
	cxxGasPhase *gas_phase_ptr = use.Get_gas_phase_ptr();
	if (gas_phase_ptr->Get_type() == cxxGasPhase::GP_PRESSURE)
	{
		if (gas_unknown == NULL)
			return (0);
		if (gas_unknown->moles < 1e-12)
			return (0);
	}
	return (gas_phase_ptr->Get_total_p());
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
find_gas_vm(void)
/* ---------------------------------------------------------------------- */
{
	if (use.Get_gas_phase_in() == FALSE || use.Get_gas_phase_ptr() == NULL)
		return (0);
	cxxGasPhase *gas_phase_ptr = use.Get_gas_phase_ptr();
	if (gas_phase_ptr->Get_type() == cxxGasPhase::GP_PRESSURE)
	{
		if (gas_unknown == NULL)
			return (0);
		if (gas_unknown->moles < 1e-12)
			return (0);
		gas_phase_ptr->Set_total_moles(gas_unknown->moles);
		gas_phase_ptr->Set_volume(gas_phase_ptr->Get_total_moles() * R_LITER_ATM * tk_x /
			gas_phase_ptr->Get_total_p());
		if (gas_phase_ptr->Get_v_m() >= 0.01)
			gas_phase_ptr->Set_volume(gas_phase_ptr->Get_v_m() * gas_unknown->moles);
	}
	return gas_phase_ptr->Get_volume() / gas_phase_ptr->Get_total_moles();
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
find_misc1(const char *ss_name)
/* ---------------------------------------------------------------------- */
{
	if (use.Get_ss_assemblage_in() == FALSE || use.Get_ss_assemblage_ptr() == NULL)
		return (0.0);
	std::vector<cxxSS *> ss_ptrs = use.Get_ss_assemblage_ptr()->Vectorize();
	for (size_t j = 0; j < ss_ptrs.size(); j++)
	{
		cxxSS *ss_ptr = ss_ptrs[j];
		if (strcmp_nocase(ss_ptr->Get_name().c_str(), ss_name) == 0)
		{
			if (ss_ptr->Get_miscibility())
			{
				return (ss_ptr->Get_xb1());
			}
			else
			{
				return (1.0);
			}
		}
	}
	return (0);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
find_misc2(const char *ss_name)
/* ---------------------------------------------------------------------- */
{
	if (use.Get_ss_assemblage_in() == FALSE || use.Get_ss_assemblage_ptr() == NULL)
		return (0.0);
	std::vector<cxxSS *> ss_ptrs = use.Get_ss_assemblage_ptr()->Vectorize();
	for (size_t j = 0; j < ss_ptrs.size(); j++)
	{
		cxxSS *ss_ptr = ss_ptrs[j];
		if (strcmp_nocase(ss_ptr->Get_name().c_str(), ss_name) == 0)
		{
			if (ss_ptr->Get_miscibility())
			{
				return (ss_ptr->Get_xb2());
			}
			else
			{
				return (1.0);
			}
		}
	}
	return (0);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
find_ss_comp(const char *ss_comp_name)
/* ---------------------------------------------------------------------- */
{
	if (use.Get_ss_assemblage_in() == FALSE || use.Get_ss_assemblage_ptr() == NULL)
		return (0);

	std::vector<cxxSS *> ss_ptrs = use.Get_ss_assemblage_ptr()->Vectorize();
	for (size_t j = 0; j < ss_ptrs.size(); j++)
	{
		cxxSS *ss_ptr = ss_ptrs[j];
		for (size_t i = 0; i < ss_ptr->Get_ss_comps().size(); i++)
		{
			cxxSScomp *comp_ptr = &(ss_ptr->Get_ss_comps()[i]);
			if (strcmp_nocase(comp_ptr->Get_name().c_str(), ss_comp_name) == 0)
			{
				if (ss_ptr->Get_ss_in())
				{
					return (comp_ptr->Get_moles());
				}
				else
				{
					return (0.0);
				}
			}
		}
	}
	return (0);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
get_calculate_value(const char *name)
/* ---------------------------------------------------------------------- */
/*
 *   Gets value from a calclate_value structure
 *	  arguments: name
 *	  return: LDBLE of value
 */
{
	struct calculate_value *calculate_value_ptr;
	calculate_value_ptr = calculate_value_search(name);
	if (calculate_value_ptr == NULL)
	{
		error_string = sformatf( "CALC_VALUE Basic function, %s not found.",
				name);
		error_msg(error_string, CONTINUE);
		input_error++;
		return (MISSING);
	}
	if (name == NULL)
	{
		error_string = sformatf(
				"Definition for calculated value not found, %s", name);
		input_error++;
		error_msg(error_string, CONTINUE);
		return (MISSING);
	}

	char l_command[] = "run";
	PBasic interp(this, this->phrq_io);
	if (calculate_value_ptr->new_def == TRUE)
	{
		if (interp.basic_compile
			(calculate_value_ptr->commands, 
			&calculate_value_ptr->linebase,
			&calculate_value_ptr->varbase,
			&calculate_value_ptr->loopbase) != 0)
		{
			error_string = sformatf(
					"Fatal Basic error in CALCULATE_VALUES %s.",
					calculate_value_ptr->name);
			error_msg(error_string, STOP);
		}
		calculate_value_ptr->new_def = FALSE;
	}

	if (interp.basic_run(l_command, 
		calculate_value_ptr->linebase,
		calculate_value_ptr->varbase, 
		calculate_value_ptr->loopbase) != 0)
	{
		error_string = sformatf( "Fatal Basic error in calculate_value %s.",
				calculate_value_ptr->name);
		error_msg(error_string, STOP);
	}
	if (rate_moles == NAN)
	{
		error_string = sformatf( "Calculated value not SAVEed for %s.",
				calculate_value_ptr->name);
		error_msg(error_string, STOP);
	}
	else
	{
		calculate_value_ptr->calculated = TRUE;
		calculate_value_ptr->value = rate_moles;
	}
	return (calculate_value_ptr->value);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
kinetics_moles(const char *kinetics_name)
/* ---------------------------------------------------------------------- */
{

	if (use.Get_kinetics_in() == FALSE || use.Get_kinetics_ptr() == NULL)
		return (0);
	for (size_t i = 0; i < use.Get_kinetics_ptr()->Get_kinetics_comps().size(); i++)
	{
		cxxKineticsComp *kinetics_comp_ptr = &(use.Get_kinetics_ptr()->Get_kinetics_comps()[i]);
		if (strcmp_nocase
			(kinetics_comp_ptr->Get_rate_name().c_str(), kinetics_name) == 0)
		{
			return (kinetics_comp_ptr->Get_m());
		}
	}

	error_string = sformatf( "No data for rate %s in KINETICS keyword.",
			kinetics_name);
	warning_msg(error_string);
	return (0);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
kinetics_moles_delta(const char *kinetics_name)
/* ---------------------------------------------------------------------- */
{

	if (use.Get_kinetics_in() == FALSE || use.Get_kinetics_ptr() == NULL)
		return (0);
	for (size_t i = 0; i < use.Get_kinetics_ptr()->Get_kinetics_comps().size(); i++)
	{
		cxxKineticsComp *kinetics_comp_ptr = &(use.Get_kinetics_ptr()->Get_kinetics_comps()[i]);
		if (strcmp_nocase
			(kinetics_comp_ptr->Get_rate_name().c_str(), kinetics_name) == 0)
		{
			//return (kinetics_comp_ptr->Get_m());

			if (state != TRANSPORT && state != PHAST)
			{
				//LDBLE moles = kinetics_comp_ptr->Get_m();
				LDBLE delta_moles = - kinetics_comp_ptr->Get_moles();
				return delta_moles;
			}
			else
			{
				//moles =  kinetics_comp_ptr->Get_m();
				LDBLE delta_moles =
						kinetics_comp_ptr->Get_m() -
						kinetics_comp_ptr->Get_initial_moles();
				return delta_moles;
			}
		}
	}

	//error_string = sformatf( "No data for rate %s in KINETICS keyword.",
	//		kinetics_name);
	//warning_msg(error_string);
	return (0);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
log_activity(const char *species_name)
/* ---------------------------------------------------------------------- */
{
	struct species *s_ptr;
	LDBLE la;

	s_ptr = s_search(species_name);

	if (s_ptr == s_eminus)
	{
		la = s_eminus->la;
	}
	else if (s_ptr == NULL || s_ptr->in == FALSE)
	{
		la = -99.99;
	}
	else if (s_ptr == s_h2o)
	{
		la = s_h2o->la;
	}
	else
	{
		la = s_ptr->lm + s_ptr->lg;
	}
	return (la);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
log_molality(const char *species_name)
/* ---------------------------------------------------------------------- */
{
	struct species *s_ptr;
	LDBLE lm;

	s_ptr = s_search(species_name);

	if (s_ptr == s_eminus)
	{
		lm = -99.99;
	}
	else if (s_ptr == NULL || s_ptr->in == FALSE)
	{
		lm = -99.99;
	}
	else if (s_ptr == s_h2o)
	{
		lm = log10(s_ptr->moles / mass_water_aq_x);
	}
	else
	{
		lm = s_ptr->lm;
	}
	return (lm);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
molality(const char *species_name)
/* ---------------------------------------------------------------------- */
{
	struct species *s_ptr;
	LDBLE m;

	s_ptr = s_search(species_name);
	if (s_ptr == NULL || s_ptr == s_eminus || s_ptr->in == FALSE)
	{
		m = 1e-99;
	}
	else
	{
		/* m = pow(10., s_ptr->lm); */
		m = s_ptr->moles / mass_water_aq_x;
	}
	return (m);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
pr_pressure(const char *phase_name)
/* ---------------------------------------------------------------------- */
{
	struct phase *phase_ptr;
	int l;

	phase_ptr = phase_bsearch(phase_name, &l, FALSE);
	if (phase_ptr == NULL)
	{
		error_string = sformatf( "Gas %s, not found.", phase_name);
		warning_msg(error_string);
		return (1e-99);
	}
	else if (phase_ptr->in != FALSE && phase_ptr->pr_in)
	{
		return phase_ptr->pr_p;
	}
	return (0.0);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
pressure(void)
/* ---------------------------------------------------------------------- */
{
	return (patm_x);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
pr_phi(const char *phase_name)
/* ---------------------------------------------------------------------- */
{
	struct phase *phase_ptr;
	int l;

	phase_ptr = phase_bsearch(phase_name, &l, FALSE);
	if (phase_ptr == NULL)
	{
		error_string = sformatf( "Gas %s, not found.", phase_name);
		warning_msg(error_string);
		return (1e-99);
	}
	else if (phase_ptr->in != FALSE && phase_ptr->pr_in)
	{
		return phase_ptr->pr_phi;
	}
	return (1.0);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
saturation_ratio(const char *phase_name)
/* ---------------------------------------------------------------------- */
{
	struct rxn_token *rxn_ptr;
	struct phase *phase_ptr;
	int l;
	LDBLE si, iap;

	iap = 0.0;
	phase_ptr = phase_bsearch(phase_name, &l, FALSE);
	if (phase_ptr == NULL)
	{
		error_string = sformatf( "Mineral %s, not found.", phase_name);
		warning_msg(error_string);
		return (1e-99);
	}
	else if (phase_ptr->in != FALSE)
	{
		for (rxn_ptr = phase_ptr->rxn_x->token + 1; rxn_ptr->s != NULL;
			 rxn_ptr++)
		{
			iap += rxn_ptr->s->la * rxn_ptr->coef;
		}
		si = iap - phase_ptr->lk;
		return (pow((LDBLE) 10.0, si));
	}
	return (0.0);

}

/* ---------------------------------------------------------------------- */
int Phreeqc::
saturation_index(const char *phase_name, LDBLE * iap, LDBLE * si)
/* ---------------------------------------------------------------------- */
{
	struct rxn_token *rxn_ptr;
	struct phase *phase_ptr;
	int l;

	*si = -99.99;
	*iap = 0.0;
	phase_ptr = phase_bsearch(phase_name, &l, FALSE);
	if (phase_ptr == NULL)
	{
		error_string = sformatf( "Mineral %s, not found.", phase_name);
		warning_msg(error_string);
		*si = -99;
	}
	else if (phase_ptr->in != FALSE)
	{
		for (rxn_ptr = phase_ptr->rxn_x->token + 1; rxn_ptr->s != NULL;
			 rxn_ptr++)
		{
			*iap += rxn_ptr->s->la * rxn_ptr->coef;
		}
		*si = *iap - phase_ptr->lk;
	}
	else
	{
		return (ERROR);
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
sum_match_gases(const char *mytemplate, const char *name)
/* ---------------------------------------------------------------------- */
{
	int i;
	LDBLE tot;
	struct elt_list *next_elt;

	if (use.Get_gas_phase_in() == FALSE || use.Get_gas_phase_ptr() == NULL)
		return (0);
	cxxGasPhase *gas_phase_ptr = use.Get_gas_phase_ptr();
	tot = 0;
	for (size_t j = 0; j < gas_phase_ptr->Get_gas_comps().size(); j++)
	{
		struct phase * phase_ptr = phase_bsearch(gas_phase_ptr->Get_gas_comps()[j].Get_phase_name().c_str(),
			&i, FALSE);
		if (match_elts_in_species(phase_ptr->formula, mytemplate) == TRUE)
		{
			if (name == NULL)
			{
				tot += phase_ptr->moles_x;
			}
			else
			{
				for (next_elt = phase_ptr->next_elt;
					 next_elt->elt != NULL; next_elt++)
				{
					if (strcmp(next_elt->elt->name, name) == 0)
					{
						tot += next_elt->coef * phase_ptr->moles_x;
						break;
					}
				}
			}
		}
	}
	return (tot);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
sum_match_species(const char *mytemplate, const char *name)
/* ---------------------------------------------------------------------- */
{
	int i;
	LDBLE tot;
	struct elt_list *next_elt;

	count_elts = 0;
	paren_count = 0;
	tot = 0;
	if (sum_species_map.find(mytemplate) == sum_species_map.end())
	{
		std::vector<std::string> species_list;
		for (i = 0; i < count_s_x; i++)
		{
			struct species *s_ptr = s_x[i];
			if (match_elts_in_species(s_ptr->name, mytemplate) == TRUE)
			{
				species_list.push_back(s_ptr->name);
			}
		}
		sum_species_map[mytemplate] = species_list;
	}
	std::vector<std::string> &species_list = (sum_species_map.find(mytemplate))->second;
	for (size_t i=0; i < species_list.size(); i++)
	{
		struct species *s_ptr = s_search(species_list[i].c_str());
		if (s_ptr->in == FALSE) continue;
		if (name == NULL)
		{
			tot += s_ptr->moles;
		}
		else
		{
			for (next_elt = s_ptr->next_elt; next_elt->elt != NULL;
					next_elt++)
			{
				if (strcmp(next_elt->elt->name, name) == 0)
				{
					tot += next_elt->coef * s_ptr->moles;
					break;
				}
			}
		}
	}
	return (tot);
}



/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
sum_match_ss(const char *mytemplate, const char *name)
/* ---------------------------------------------------------------------- */
{
	LDBLE tot;
	struct elt_list *next_elt;

	if (use.Get_ss_assemblage_in() == FALSE || use.Get_ss_assemblage_ptr() == NULL)
		return (0);
	tot = 0;
	std::vector<cxxSS *> ss_ptrs = use.Get_ss_assemblage_ptr()->Vectorize();
	for (size_t j = 0; j < ss_ptrs.size(); j++)
	{
		cxxSS *ss_ptr = ss_ptrs[j];
		if (strcmp_nocase(ss_ptr->Get_name().c_str(), mytemplate) == 0)
		{
			if (!ss_ptr->Get_ss_in())
			{
				tot = 0;
				break;
			}
			for (size_t i = 0; i < ss_ptr->Get_ss_comps().size(); i++)
			{
				cxxSScomp *comp_ptr = &(ss_ptr->Get_ss_comps()[i]);
				if (name == NULL)
				{
					tot += comp_ptr->Get_moles();
				}
				else
				{
					int l;
					struct phase *phase_ptr = phase_bsearch(comp_ptr->Get_name().c_str(), &l, FALSE);
					for (next_elt = phase_ptr->next_elt; next_elt->elt != NULL; next_elt++)
					{
						if (strcmp(next_elt->elt->name, name) == 0)
						{
							tot += next_elt->coef *	comp_ptr->Get_moles();
							break;
						}
					}
				}
			}
			break;
		}
	}
	return (tot);
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
list_ss(std::string ss_name, cxxNameDouble &composition)
/* ---------------------------------------------------------------------- */
{
	LDBLE tot = 0;
	composition.clear();
	if (use.Get_ss_assemblage_in() == FALSE || use.Get_ss_assemblage_ptr() == NULL)
		return (0);

	std::vector<cxxSS *> ss_ptrs = use.Get_ss_assemblage_ptr()->Vectorize();
	for (size_t j = 0; j < ss_ptrs.size(); j++)
	{
		cxxSS *ss_ptr = ss_ptrs[j];
		if (strcmp_nocase(ss_ptr->Get_name().c_str(), ss_name.c_str()) == 0)
		{
			for (size_t i = 0; i < ss_ptr->Get_ss_comps().size(); i++)
			{
				cxxSScomp *comp_ptr = &(ss_ptr->Get_ss_comps()[i]);
				composition.add(comp_ptr->Get_name().c_str(), comp_ptr->Get_moles());
				tot += comp_ptr->Get_moles();
			}
			break;
		}
	}
	return (tot);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
match_elts_in_species(const char *name, const char *mytemplate)
/* ---------------------------------------------------------------------- */
{
/*
 *	Makes a list of elements with their coefficients, stores elements
 *	in elt_list at position count_elts.  Global variable count_elts is
 *	updated with each stored element.  Also uses static global variable
 *	paren_count.
 *
 *	Arguments:
 *	   **t_ptr	input, point in token string to start looking
 *				  output, is next position to start looking
 *		 coef	 input, coefficient to multiply subscripts by
 */
	int i, i1, l, case_no, match;
	char c, c1;
	char *ptr, *ptr1;
	LDBLE d;
	char element[MAX_LENGTH];
	char token[MAX_LENGTH], equal_list[MAX_LENGTH]; 
	char token1[MAX_LENGTH], template1[MAX_LENGTH], equal_list1[MAX_LENGTH];
	char str[2];

	strcpy(token, name);
	squeeze_white(token);
	replace("(+", "(", token);
	if (strstr("token", "++") != NULL)
	{
		replace("++++++", "+6", token);
		replace("+++++", "+5", token);
		replace("++++", "+4", token);
		replace("+++", "+3", token);
		replace("++", "+2", token);
	}
	if (strstr("token", "--") != NULL)
	{
		replace("------", "-6", token);
		replace("-----", "-5", token);
		replace("----", "-4", token);
		replace("---", "-3", token);
		replace("--", "-2", token);
	}

	ptr = token;
	std::vector< std::pair<std::string, LDBLE> > match_vector;
	while ((c = *ptr) != '\0')
	{
		c1 = *(ptr + 1);
		str[0] = c;
		str[1] = '\0';
/*
 * New element
 */
		if (isupper((int) c) || (c == 'e' && c1 == '-') || (c == '['))
		{
			/*
			 *   Get new element and subscript
			 */
			if (get_elt(&ptr, element, &l) == ERROR)
			{
				return (ERROR);
			}
			if (get_num(&ptr, &d) == ERROR)
			{
				return (ERROR);
			}
			std::pair<std::string, LDBLE> pr(element, d);
			match_vector.push_back(pr);			
		}
		else
		{
			std::pair<std::string, LDBLE> pr(str, 1.0);
			match_vector.push_back(pr);		
			ptr += 1;
		}
	}
	/*
	 *  Replace elements with first of equivalent elements
	 */
	strcpy(template1, mytemplate);
	squeeze_white(template1);
	ptr = template1;
	while (extract_bracket(&ptr, equal_list) == TRUE)
	{
		replace("{", "", equal_list);
		replace("}", "", equal_list);
		while (replace(",", " ", equal_list) == TRUE);
		ptr1 = equal_list;
		/*
		 *   Get first name in a list from template
		 */
		std::string elt_name;
		if (copy_token(elt_name, &ptr1) == EMPTY)
		{
			error_string = sformatf(
					"Expecting a nonempty list of element names in isotope sum. %s",
					mytemplate);
			error_msg(error_string, CONTINUE);
			return (ERROR);
		}
		std::string replace_name = elt_name;
		/*
		 *   Replace in species all equivalent names from template
		 */
		while (copy_token(elt_name, &ptr1) != EMPTY)
		{
			for (i = 0; i < (int) match_vector.size(); i++)
			{
				if (elt_name == match_vector[i].first)
				{
					match_vector[i].first = replace_name;
				}
			}
		}
	}
	/*
	 *  Combine contiguous elements
	 */
	i1 = 0;		
	for (i = 1; i < (int) match_vector.size(); i++)
	{
		if ((isupper((int) (match_vector[i].first[0])) != FALSE)
		&& (match_vector[i].first == match_vector[i1].first))
		{
			match_vector[i1].second += match_vector[i].second;
		}
		else
		{
			i1++;
			match_vector[i1].first = match_vector[i].first;
			match_vector[i1].second = match_vector[i].second;
		}
	}
	int count_match_tokens = i1 + 1;
	/*
	 * write out string
	 */
	token[0] = '\0';
	for (i = 0; i < count_match_tokens; i++)
	{
		strcat(token, match_vector[i].first.c_str());
		if (match_vector[i].second != 1.0)
		{
			sprintf(token1, "%g", (double) match_vector[i].second);
			strcat(token, token1);
		}
	}
	/*
	 *  Write a template name using first of equivalent elements
	 */
	strcpy(template1, mytemplate);
	squeeze_white(template1);
	ptr = template1;
	while (extract_bracket(&ptr, equal_list) == TRUE)
	{
		strcpy(equal_list1, equal_list);
		replace("{", "", equal_list);
		replace("}", "", equal_list);
		while (replace(",", " ", equal_list) == TRUE);
		ptr1 = equal_list;
		/*
		 *   Get first name in a list
		 */
		std::string elt_name;
		if (copy_token(elt_name, &ptr1) == EMPTY)
		{
			error_string = sformatf(
					"Expecting a nonempty list of element names in isotope sum. %s",
					mytemplate);
			error_msg(error_string, CONTINUE);
			return (ERROR);
		}
		replace(equal_list1, elt_name.c_str(), template1);
		squeeze_white(template1);
		ptr = template1;
	}
	/*
	 *   Compare string
	 */
	/* Cases: 0 exact match
	 *	  1 leading wild card
	 *	  2 trailing wild card
	 *	  3 leading and trailing wild card
	 */
	case_no = 0;
	if (template1[0] == '*')
		case_no = 1;
	l = (int) strlen(template1);
	if (template1[l - 1] == '*')
	{
		if (case_no != 1)
		{
			case_no = 2;
		}
		else
		{
			case_no = 3;
		}
	}
	while (replace("*", "", template1));
	match = FALSE;
	switch (case_no)
	{
	case 0:
		/* exact match */
		if (strcmp(token, template1) == 0)
			match = TRUE;
		break;
	case 1:
		/* leading wild card */
		if ((ptr = strstr(token, template1)) == NULL)
		{
			match = FALSE;
		}
		else
		{
			if (strcmp(ptr, template1) == 0)
				match = TRUE;
		}
		break;
	case 2:
		/* trailing wild card */
		if (strstr(token, template1) == token)
			match = TRUE;
		break;
	case 3:
		/* trailing wild card */
		if (strstr(token, template1) != NULL)
			match = TRUE;
		break;
	}
	return (match);
}
#ifdef SKIP
/* ---------------------------------------------------------------------- */
int Phreeqc::
match_elts_in_species(const char *name, const char *mytemplate)
/* ---------------------------------------------------------------------- */
{
/*
 *	Makes a list of elements with their coefficients, stores elements
 *	in elt_list at position count_elts.  Global variable count_elts is
 *	updated with each stored element.  Also uses static global variable
 *	paren_count.
 *
 *	Arguments:
 *	   **t_ptr	input, point in token string to start looking
 *				  output, is next position to start looking
 *		 coef	 input, coefficient to multiply subscripts by
 */
	int i, i1, l, case_no, match;
	char c, c1;
	char *ptr, *ptr1;
	const char *replace_name, *name_ptr;
	LDBLE d;
	char element[MAX_LENGTH];
	char token[MAX_LENGTH], equal_list[MAX_LENGTH], elt_name[MAX_LENGTH];
	char token1[MAX_LENGTH], template1[MAX_LENGTH], equal_list1[MAX_LENGTH];
	char str[2];

	strcpy(token, name);
	squeeze_white(token);
	while (replace("(+", "(", token));
	while (replace("++++++", "+6", token));
	while (replace("+++++", "+5", token));
	while (replace("++++", "+4", token));
	while (replace("+++", "+3", token));
	while (replace("++", "+2", token));
	while (replace("------", "-6", token));
	while (replace("-----", "-5", token));
	while (replace("----", "-4", token));
	while (replace("---", "-3", token));
	while (replace("--", "-2", token));

	ptr = token;
	count_match_tokens = 0;
	while ((c = *ptr) != '\0')
	{
		c1 = *(ptr + 1);
		str[0] = c;
		str[1] = '\0';
/*
 * New element
 */
		if (isupper((int) c) || (c == 'e' && c1 == '-') || (c == '['))
		{
			/*
			 *   Get new element and subscript
			 */
			if (get_elt(&ptr, element, &l) == ERROR)
			{
				return (ERROR);
			}
			match_tokens[count_match_tokens].name = string_hsave(element);
			if (get_num(&ptr, &d) == ERROR)
			{
				return (ERROR);
			}
			match_tokens[count_match_tokens++].coef = d;
		}
		else
		{
			match_tokens[count_match_tokens].name = string_hsave(str);
			match_tokens[count_match_tokens++].coef = 1.0;
			ptr += 1;
		}
	}
	/*
	 *  Replace elements with first of equivalent elements
	 */
	strcpy(template1, mytemplate);
	squeeze_white(template1);
	ptr = template1;
	while (extract_bracket(&ptr, equal_list) == TRUE)
	{
		replace("{", "", equal_list);
		replace("}", "", equal_list);
		while (replace(",", " ", equal_list) == TRUE);
		ptr1 = equal_list;
		/*
		 *   Get first name in a list from template
		 */
		if (copy_token(elt_name, &ptr1, &l) == EMPTY)
		{
			error_string = sformatf(
					"Expecting a nonempty list of element names in isotope sum. %s",
					mytemplate);
			error_msg(error_string, CONTINUE);
			return (ERROR);
		}
		replace_name = string_hsave(elt_name);
		/*
		 *   Replace in species all equivalent names from template
		 */
		while (copy_token(elt_name, &ptr1, &l) != EMPTY)
		{
			name_ptr = string_hsave(elt_name);
			for (i = 0; i < count_match_tokens; i++)
			{
				if (name_ptr == match_tokens[i].name)
				{
					match_tokens[i].name = replace_name;
				}
			}
		}
	}
	/*
	 *  Combine contiguous elements
	 */
	i1 = 0;
	for (i = 1; i < count_match_tokens; i++)
	{
		if ((isupper((int) (match_tokens[i].name[0])) != FALSE)
			&& (match_tokens[i].name == match_tokens[i1].name))
		{
			match_tokens[i1].coef += match_tokens[i].coef;
		}
		else
		{
			i1++;
			match_tokens[i1].name = match_tokens[i].name;
			match_tokens[i1].coef = match_tokens[i].coef;
		}
	}
	count_match_tokens = i1 + 1;
	/*
	 * write out string
	 */
	token[0] = '\0';
	for (i = 0; i < count_match_tokens; i++)
	{
		strcat(token, match_tokens[i].name);
		if (match_tokens[i].coef != 1.0)
		{
			sprintf(token1, "%g", (double) match_tokens[i].coef);
			strcat(token, token1);
		}
	}
	/*
	 *  Write a template name using first of equivalent elements
	 */
	strcpy(template1, mytemplate);
	squeeze_white(template1);
	ptr = template1;
	while (extract_bracket(&ptr, equal_list) == TRUE)
	{
		strcpy(equal_list1, equal_list);
		replace("{", "", equal_list);
		replace("}", "", equal_list);
		while (replace(",", " ", equal_list) == TRUE);
		ptr1 = equal_list;
		/*
		 *   Get first name in a list
		 */
		if (copy_token(elt_name, &ptr1, &l) == EMPTY)
		{
			error_string = sformatf(
					"Expecting a nonempty list of element names in isotope sum. %s",
					mytemplate);
			error_msg(error_string, CONTINUE);
			return (ERROR);
		}
		replace_name = string_hsave(elt_name);
		replace(equal_list1, replace_name, template1);
		squeeze_white(template1);
		ptr = template1;
	}
	/*
	 *   Compare string
	 */
	/* Cases: 0 exact match
	 *	  1 leading wild card
	 *	  2 trailing wild card
	 *	  3 leading and trailing wild card
	 */
	case_no = 0;
	if (template1[0] == '*')
		case_no = 1;
	l = (int) strlen(template1);
	if (template1[l - 1] == '*')
	{
		if (case_no != 1)
		{
			case_no = 2;
		}
		else
		{
			case_no = 3;
		}
	}
	while (replace("*", "", template1));
	match = FALSE;
	switch (case_no)
	{
	case 0:
		/* exact match */
		if (strcmp(token, template1) == 0)
			match = TRUE;
		break;
	case 1:
		/* leading wild card */
		if ((ptr = strstr(token, template1)) == NULL)
		{
			match = FALSE;
		}
		else
		{
			if (strcmp(ptr, template1) == 0)
				match = TRUE;
		}
		break;
	case 2:
		/* trailing wild card */
		if (strstr(token, template1) == token)
			match = TRUE;
		break;
	case 3:
		/* trailing wild card */
		if (strstr(token, template1) != NULL)
			match = TRUE;
		break;
	}
	return (match);
}
#endif
/* ---------------------------------------------------------------------- */
int Phreeqc::
extract_bracket(char **string, char *bracket_string)
/* ---------------------------------------------------------------------- */
{
	char *ptr, *ptr1;

	if ((ptr = strstr(*string, "{")) == NULL)
		return (FALSE);
	strcpy(bracket_string, ptr);
	if ((ptr1 = strstr(bracket_string, "}")) == NULL)
	{
		error_string = sformatf(
				"No matching bracket (}) in isotope template string %s",
				*string);
		error_msg(error_string, CONTINUE);
		input_error++;
		return (FALSE);
	}
	ptr1[1] = '\0';
	*string = strstr(*string, "}");
	*string += 1;
	return (TRUE);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
surf_total(const char *total_name, const char *surface_name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in LDBLE layer
 */
	int j;

	if (use.Get_surface_ptr() == NULL || surface_name == NULL || total_name == NULL)
		return (0);

	bool redox = false;
	if (strstr(total_name, "(") != NULL)
	{
		redox = true;
	}
	if (!redox)
	{
		if (strcmp(total_name, "H") == 0 || strcmp(total_name, "O") == 0)
		{
			return surf_total_no_redox(total_name, surface_name);
		}
	}	
/*
 *   Find surface...
 */
	for (j = 0; j < count_unknowns; j++)
	{
		if (x[j]->type != SURFACE)
			continue;
		
		std::string token;
		token = x[j]->master[0]->elt->name;
		replace("_", " ", token);
		std::string::iterator b = token.begin();
		std::string::iterator e = token.end();
		std::string name;
		CParser::copy_token(name, b, e);
		if (strcmp(name.c_str(), surface_name) == 0)
				break;
	}
	if (j >= count_unknowns)
		return (0);
/*
 *   find total moles for redox state
 */
	LDBLE t = 0;
	for (j = 0; j < count_s_x; j++)
	{
		if (s_x[j]->type != SURF)
			continue;

		std::string token;
		bool match = false; 

		// find if surface matches 
		for (int i = 0; s_x[j]->next_elt[i].elt != NULL; i++)
		{
			if (s_x[j]->next_elt[i].elt->master->type != SURF) continue;

			//strcpy(token, s_x[j]->next_elt[i].elt->name);
			//replace("_", " ", token);
			//ptr = token;
			//copy_token(name, &ptr, &k);
			token = s_x[j]->next_elt[i].elt->name;
			replace("_", " ", token);
			std::string::iterator b = token.begin();
			std::string::iterator e = token.end();
			std::string name;
			CParser::copy_token(name, b, e);
			if (strcmp(name.c_str(), surface_name) == 0)
			{
				match = true;
				break;
			}
		}
		if (!match) continue;

		// surface matches, now match element or redox state
		struct rxn_token *rxn_ptr;
		for (rxn_ptr = s_x[j]->rxn_s->token + 1; rxn_ptr->s != NULL; rxn_ptr++)
		{
			if (redox && rxn_ptr->s->secondary)
			{
				token = rxn_ptr->s->secondary->elt->name;
			}
			else if (!redox && rxn_ptr->s->secondary)
			{
				token = rxn_ptr->s->secondary->elt->primary->elt->name;
			}
			else if (!redox && rxn_ptr->s->primary)
			{
				token = rxn_ptr->s->primary->elt->name;
			}
			else
			{
				continue;
			}
			if (strcmp(token.c_str(), total_name) == 0)
			{
				t += rxn_ptr->coef * s_x[j]->moles;
				break;
			}
			else
			// sum all sites in case total_name is a surface name without underscore surf ("Hfo_w", "Hfo")
			{
				if (rxn_ptr->s->type == SURF)
				{
					if (token.find("_") != std::string::npos)
					{
						token = token.substr(0, token.find("_"));
					}
					if (strcmp(token.c_str(), total_name) == 0)
					{
						t += rxn_ptr->coef * s_x[j]->moles;
						break;
					}
				}
			}
		}
	}
	return t;
}
#ifdef SKIP
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
surf_total(const char *total_name, const char *surface_name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in LDBLE layer
 */
	int j;

	if (use.Get_surface_ptr() == NULL || surface_name == NULL)
		return (0);

/*
 *   Find surface...
 */
	for (j = 0; j < count_unknowns; j++)
	{
		if (x[j]->type != SURFACE)
			continue;
		
		std::string token;
		token = x[j]->master[0]->elt->name;
		replace("_", " ", token);
		std::string::iterator b = token.begin();
		std::string::iterator e = token.end();
		std::string name;
		CParser::copy_token(name, b, e);
		if (strcmp(name.c_str(), surface_name) == 0)
				break;
	}
	if (j >= count_unknowns)
		return (0);
/*
 *   find total moles for redox state
 */
	LDBLE t = 0;
	bool redox = false;
	if (strstr(total_name, "(") != NULL)
	{
		redox = true;
	}
	for (j = 0; j < count_s_x; j++)
	{
		if (s_x[j]->type != SURF)
			continue;

		std::string token;
		bool match = false; 

		// find if surface matches 
		for (int i = 0; s_x[j]->next_elt[i].elt != NULL; i++)
		{
			if (s_x[j]->next_elt[i].elt->master->type != SURF) continue;

			//strcpy(token, s_x[j]->next_elt[i].elt->name);
			//replace("_", " ", token);
			//ptr = token;
			//copy_token(name, &ptr, &k);
			token = s_x[j]->next_elt[i].elt->name;
			replace("_", " ", token);
			std::string::iterator b = token.begin();
			std::string::iterator e = token.end();
			std::string name;
			CParser::copy_token(name, b, e);
			if (strcmp(name.c_str(), surface_name) == 0)
			{
				match = true;
				break;
			}
		}
		if (!match) continue;

		// surface matches, now match element or redox state
		struct rxn_token *rxn_ptr;
		for (rxn_ptr = s_x[j]->rxn_s->token + 1; rxn_ptr->s != NULL; rxn_ptr++)
		{
			if (redox && rxn_ptr->s->secondary)
			{
				token = rxn_ptr->s->secondary->elt->name;
			}
			else if (!redox && rxn_ptr->s->secondary)
			{
				token = rxn_ptr->s->secondary->elt->primary->elt->name;
			}
			else if (!redox && rxn_ptr->s->primary)
			{
				token = rxn_ptr->s->primary->elt->name;
			}
			else
			{
				continue;
			}
			if (strcmp(token.c_str(), total_name) == 0)
			{
				t += rxn_ptr->coef * s_x[j]->moles;
				break;
			}
			else
			// sum all sites in case total_name is a surface name without underscore surf ("Hfo_w", "Hfo")
			{
				if (rxn_ptr->s->type == SURF)
				{
					if (token.find("_") != std::string::npos)
					{
						token = token.substr(0, token.find("_"));
					}
					if (strcmp(token.c_str(), total_name) == 0)
					{
						t += rxn_ptr->coef * s_x[j]->moles;
						break;
					}
				}
			}
		}
	}
	return t;
}
#endif
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
surf_total_no_redox(const char *total_name, const char *surface_name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in LDBLE layer
 */
	int i, j, k;
	char name[MAX_LENGTH], token[MAX_LENGTH];
	char surface_name_local[MAX_LENGTH];
	char *ptr;

	if (use.Get_surface_ptr() == NULL)
		return (0);

/*
 *   Find surface...
 */
	for (j = 0; j < count_unknowns; j++)
	{
		if (x[j]->type != SURFACE)
			continue;
		strcpy(token, x[j]->master[0]->elt->name);
		replace("_", " ", token);
		ptr = token;
		copy_token(name, &ptr, &k);
		if (surface_name != NULL)
		{
			if (strcmp(name, surface_name) == 0)
				break;
		}
		else
		{
			break;
		}
	}
	if (j >= count_unknowns)
		return (0);
	strcpy(surface_name_local, name);
/*
 *   find total moles of each element in diffuse layer...
 */
	count_elts = 0;
	paren_count = 0;
	for (j = 0; j < count_s_x; j++)
	{
		if (s_x[j]->type != SURF)
			continue;
		for (i = 0; s_x[j]->next_elt[i].elt != NULL; i++)
		{
			if (s_x[j]->next_elt[i].elt->master->type != SURF) continue;

			strcpy(token, s_x[j]->next_elt[i].elt->name);
			replace("_", " ", token);
			ptr = token;
			copy_token(name, &ptr, &k);
			if (strcmp(name, surface_name_local) == 0)
			{
/*
 *   Accumulate elements in diffuse layer
 */
				add_elt_list(s_x[j]->next_elt, s_x[j]->moles);
				//fprintf(stderr, "%15s\t%e\t%s\t%s\n", s_x[j]->name, s_x[j]->moles, name, surface_name_local );
				break;
			}
		}
	}
	if (count_elts > 0)
	{
		qsort(elt_list, (size_t) count_elts,
			  (size_t) sizeof(struct elt_list), elt_list_compare);
		elt_list_combine();
	}
/*
 *   Return totals
 */
	for (j = 0; j < count_elts; j++)
	{
		if (strcmp(elt_list[j].elt->name, total_name) == 0)
		{
			return ((LDBLE) elt_list[j].coef);
		}
	}
	return (0);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
total(const char *total_name)
/* ---------------------------------------------------------------------- */
{
	struct master *master_ptr;
	LDBLE t;
	int i;

	if (strcmp(total_name, "H") == 0)
	{
		return (total_h_x / mass_water_aq_x);
	}
	if (strcmp(total_name, "O") == 0)
	{
		return (total_o_x / mass_water_aq_x);
	}
	master_ptr = master_bsearch(total_name);
	t = 0.0;
	if (master_ptr == NULL)
	{
		if (strcmp_nocase(total_name, "water") == 0)
		{
			return (mass_water_aq_x);
		}
		else if (strcmp_nocase(total_name, "charge") == 0)
		{
			return (cb_x / mass_water_aq_x);
		}
/*
        sprintf (error_string, "Cannot find definition for master species, %s.",
	         total_name);
        warning_msg (error_string);
*/
	}
/*
 *  Primary master species
 */
	else if (master_ptr->primary == TRUE)
	{
		/*
		 *  Not a redox element
		 */
		if (master_ptr->s->secondary == NULL)
		{
			t = master_ptr->total / mass_water_aq_x;
			/*
			 * Redox element, need to sum totals of all redox states
			 */
		}
		else
		{
			t = 0;
			for (i = master_ptr->number + 1;
				 (i < count_master && master[i]->elt->primary == master_ptr);
				 i++)
			{
				t += master[i]->total / mass_water_aq_x;
			}
		}
	}
/*
 *  Secondary master species
 */
	else
	{
		t = master_ptr->total / mass_water_aq_x;
	}
	return (t);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
total_mole(const char *total_name)
/* ---------------------------------------------------------------------- */
{
	struct master *master_ptr;
	LDBLE t;
	int i;

	if (strcmp(total_name, "H") == 0)
	{
		return (total_h_x);
	}
	if (strcmp(total_name, "O") == 0)
	{
		return (total_o_x);
	}
	master_ptr = master_bsearch(total_name);
	t = 0.0;
	if (master_ptr == NULL)
	{
		if (strcmp_nocase(total_name, "water") == 0)
		{
			return (mass_water_aq_x / gfw_water);
		}
		else if (strcmp_nocase(total_name, "charge") == 0)
		{
			return (cb_x);
		}
/*
        sprintf (error_string, "Cannot find definition for master species, %s.",
	         total_name);
        warning_msg (error_string);
*/
	}
/*
 *  Primary master species
 */
	else if (master_ptr->primary == TRUE)
	{
		/*
		 *  Not a redox element
		 */
		if (master_ptr->s->secondary == NULL)
		{
			t = master_ptr->total;
			/*
			 * Redox element, need to sum totals of all redox states
			 */
		}
		else
		{
			t = 0;
			for (i = master_ptr->number + 1;
				 (i < count_master && master[i]->elt->primary == master_ptr);
				 i++)
			{
				t += master[i]->total;
			}
		}
	}
/*
 *  Secondary master species
 */
	else
	{
		t = master_ptr->total;
	}
	return (t);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
get_edl_species(cxxSurfaceCharge & charge_ref)
/* ---------------------------------------------------------------------- */
{

	double mass_water_surface = charge_ref.Get_mass_water();
	space((void **) ((void *) &sys), count_s_x, &max_sys, sizeof(struct system_species));
	count_sys = 0;
	for (int j = 0; j < count_s_x; j++)
	{
		if (s_x[j]->type == H2O)
		{
			sys[count_sys].name = string_duplicate(s_x[j]->name);
			sys[count_sys].moles = mass_water_surface / gfw_water;
			sys_tot += sys[count_sys].moles;
			count_sys++;
		}
		else if (s_x[j]->type < H2O)
		{
			double molality = under(s_x[j]->lm);
			double moles_excess = mass_water_aq_x * molality * charge_ref.Get_g_map()[s_x[j]->z].Get_g();
			double moles_surface = mass_water_surface * molality + moles_excess;
			sys[count_sys].name = string_duplicate(s_x[j]->name);
			sys[count_sys].moles = moles_surface;
			sys_tot += sys[count_sys].moles;
			count_sys++;
#ifdef SKIP
			double g = charge_ref.Get_g_map()[s_x[j]->z].Get_g();
			double moles_excess = mass_water_aq_x * molality * (g * s_x[j]->erm_ddl +
				mass_water_surface /
				mass_water_aq_x * (s_x[j]->erm_ddl - 1));
			double c = (mass_water_surface * molality + moles_excess) / mass_water_surface;
			charge_ref.Get_dl_species_map()[s_x[j]->number] = c;
#endif
		}
		else
		{
			continue;
		}
	}
#ifdef SKIP
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */
	int i;
/*
 *   find total moles in aq, surface, and exchange
 */

	for (i = 0; i < count_s_x; i++)
	{
		//if (s_x[i]->type != AQ)
		if (s_x[i]->type > HPLUS)
			continue;
		sys[count_sys].name = string_duplicate(s_x[i]->name);
		sys[count_sys].moles = s_x[i]->moles;
		sys_tot += sys[count_sys].moles;
		sys[count_sys].type = string_duplicate("aq");
		count_sys++;
		space((void **) ((void *) &sys), count_sys, &max_sys,
			  sizeof(struct system_species));
	}
	
#endif
	return (OK);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
edl_species(const char *surf_name, LDBLE * count, char ***names, LDBLE ** moles, LDBLE * area, LDBLE * thickness)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */
	int i;
	sys_tot = 0;
	count_sys = 0;
	max_sys = 100;
	space((void **) ((void *) &sys), INIT, &max_sys,
		  sizeof(struct system_species));
	if (!(dl_type_x == cxxSurface::NO_DL))
	{
		cxxSurface *surface_ptr = use.Get_surface_ptr();
		for (size_t i = 0; i < surface_ptr->Get_surface_charges().size(); i++)
		{
			cxxSurfaceCharge & charge_ref = surface_ptr->Get_surface_charges()[i];
			if (strcmp(charge_ref.Get_name().c_str(), surf_name) == 0)
			{	
				get_edl_species(charge_ref);
				*area = charge_ref.Get_specific_area() * charge_ref.Get_grams();
				*thickness = surface_ptr->Get_thickness();
				break;
			}
		}
	}
	/*
	 *   Sort system species
	 */
	if (count_sys > 1)
	{
		qsort(sys, (size_t) count_sys,
			  (size_t) sizeof(struct system_species), system_species_compare);
	}
	/*
	 * malloc space
	 */
	*names = (char **) PHRQ_malloc((size_t) (count_sys + 1) * sizeof(char *));
	if (names == NULL)
		malloc_error();
	*moles = (LDBLE *) PHRQ_malloc((size_t) (count_sys + 1) * sizeof(LDBLE));
	if (moles == NULL)
		malloc_error();

	(*names)[0] = NULL;
	(*moles)[0] = 0;
	for (i = 0; i < count_sys; i++)
	{
		(*names)[i + 1] = sys[i].name;
		(*moles)[i + 1] = sys[i].moles;
	}
	*count = (LDBLE) count_sys;

	PHRQ_free(sys);
	return (sys_tot);
}				
				/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
system_total(const char *total_name, LDBLE * count, char ***names,
			 char ***types, LDBLE ** moles)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */
	int i;

	sys_tot = 0;
	count_sys = 0;
	max_sys = 100;
	space((void **) ((void *) &sys), INIT, &max_sys,
		  sizeof(struct system_species));
	if (strcmp_nocase(total_name, "elements") == 0)
	{
		system_total_elements();
	}
	else if (strcmp_nocase(total_name, "phases") == 0)
	{
		system_total_si();
	}
	else if (strcmp_nocase(total_name, "aq") == 0)
	{
		system_total_aq();
	}
	else if (strcmp_nocase(total_name, "ex") == 0)
	{
		system_total_ex();
	}
	else if (strcmp_nocase(total_name, "surf") == 0)
	{
		system_total_surf();
	}
	else if (strcmp_nocase(total_name, "s_s") == 0)
	{
		system_total_ss();
	}
	else if (strcmp_nocase(total_name, "gas") == 0)
	{
		system_total_gas();
	}
	else if (strcmp_nocase(total_name, "equi") == 0)
	{
		system_total_equi();
	}
	else
	{
		if (strstr(total_name, "(") == NULL)
		{
			system_total_elt(total_name);
		}
		else
		{
			system_total_elt_secondary(total_name);
		}
	}
	/*
	 *   Sort system species
	 */
	if (count_sys > 1)
	{
		qsort(sys, (size_t) count_sys,
			  (size_t) sizeof(struct system_species), system_species_compare);
	}
	/*
	 * malloc space
	 */
	*names = (char **) PHRQ_malloc((size_t) (count_sys + 1) * sizeof(char *));
	if (names == NULL)
		malloc_error();
	*types = (char **) PHRQ_malloc((size_t) (count_sys + 1) * sizeof(char *));
	if (types == NULL)
		malloc_error();
	*moles = (LDBLE *) PHRQ_malloc((size_t) (count_sys + 1) * sizeof(LDBLE));
	if (moles == NULL)
		malloc_error();

	(*names)[0] = NULL;
	(*types)[0] = NULL;
	(*moles)[0] = 0;
	for (i = 0; i < count_sys; i++)
	{
		(*names)[i + 1] = sys[i].name;
		(*types)[i + 1] = sys[i].type;
		(*moles)[i + 1] = sys[i].moles;
	}
	*count = (LDBLE) count_sys;
	if (strcmp_nocase(total_name, "elements") == 0)
	{
		sys_tot = 0;;
		for (i = 0; i < count_sys; i++)
		{
			if (strcmp(sys[i].type, "dis") == 0 &&
				strstr(sys[i].name, "(") == NULL &&
				strcmp(sys[i].name, "H") != 0
				&& strcmp(sys[i].name, "O") != 0)
			{
				sys_tot += sys[i].moles;
			}
		}
	}
	PHRQ_free(sys);
	return (sys_tot);
}

/* ---------------------------------------------------------------------- */
std::string Phreeqc::
phase_formula(std::string phase_name, cxxNameDouble &stoichiometry)
/* ---------------------------------------------------------------------- */
{
/*
 *   Returns formula of mineral
 *   Also returns arrays of elements and stoichiometry in elts_arg and coef_arg
 */
	stoichiometry.clear();
	std::string formula;

	int j;
	struct phase *phase_ptr = phase_bsearch(phase_name.c_str(), &j, FALSE);
	if (phase_ptr != NULL)
	{
		formula.append(phase_ptr->formula);
		cxxNameDouble nd(phase_ptr->next_elt);
		stoichiometry = nd;
	}

	return (formula);
}
/* ---------------------------------------------------------------------- */
std::string Phreeqc::
species_formula(std::string phase_name, cxxNameDouble &stoichiometry)
/* ---------------------------------------------------------------------- */
{
/*
 *   Returns formula of mineral
 *   Also returns arrays of elements and stoichiometry in elts_arg and coef_arg
 */
	stoichiometry.clear();
	std::string formula;
	formula = "none";
	struct species *s_ptr = s_search(phase_name.c_str());
	if (s_ptr != NULL)
	{
		cxxNameDouble nd(s_ptr->next_elt);
		stoichiometry = nd;
		stoichiometry["charge"] = s_ptr->z;
		if (s_ptr->type == EX)
		{
			formula = "ex";
		}
		else if (s_ptr->type == SURF)
		{
			formula = "surf";
		}
		else
		{
			formula = "aq";
		}
	}
	return (formula);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
system_total_elements(void)
/* ---------------------------------------------------------------------- */
{
	int i, j;
	LDBLE t;
	char name[MAX_LENGTH];
	struct master *master_ptr;

	/*
	 * Include H and O
	 */
	sys[count_sys].name = string_duplicate("H");
	sys[count_sys].moles = total_h_x;
	sys_tot += sys[count_sys].moles;
	sys[count_sys].type = string_duplicate("dis");
	count_sys++;
	space((void **) ((void *) &sys), count_sys, &max_sys,
		  sizeof(struct system_species));
	sys[count_sys].name = string_duplicate("O");
	sys[count_sys].moles = total_o_x;
	sys_tot += sys[count_sys].moles;
	sys[count_sys].type = string_duplicate("dis");
	count_sys++;
	space((void **) ((void *) &sys), count_sys, &max_sys,
		  sizeof(struct system_species));
	/*
	 * Include H(1) and O(-2)
	 */
	sys[count_sys].name = string_duplicate("H(1)");
	sys[count_sys].moles = solution_sum_secondary("H(1)");
	sys_tot += sys[count_sys].moles;
	sys[count_sys].type = string_duplicate("dis");
	count_sys++;
	space((void **) ((void *) &sys), count_sys, &max_sys,
		  sizeof(struct system_species));
	sys[count_sys].name = string_duplicate("O(-2)");
	sys[count_sys].moles = solution_sum_secondary("O(-2)");
	sys_tot += sys[count_sys].moles;
	sys[count_sys].type = string_duplicate("dis");
	count_sys++;
	space((void **) ((void *) &sys), count_sys, &max_sys,
		  sizeof(struct system_species));

	for (i = 0; i < count_master; i++)
	{
		master_ptr = master[i];
		if (master_ptr->primary == TRUE && master_ptr->total_primary <= 0)
			continue;
		if (master_ptr->in == FALSE
			&& (master_ptr->primary == FALSE
				|| master_ptr->total_primary == 0))
			continue;
		/*
		 *  H and O
		 */
		if (master_ptr->s == s_hplus)
		{
			continue;
		}
		else if (master_ptr->s == s_h2o)
		{
			continue;
		}
		if (master_ptr->primary == TRUE)
		{
			if (master_ptr->total_primary > 0)
			{
				t = master_ptr->total_primary;
				/*
				 *  Not a redox element
				 */
			}
			else if (master_ptr->s->secondary == NULL)
			{
				t = master_ptr->total;
				/*
				 * Redox element, need to sum totals of all redox states
				 */
			}
			else
			{
				t = 0;
				for (j = master_ptr->number + 1;
					 master[j]->elt->primary == master_ptr; j++)
				{
					t += master[j]->total;
				}
			}
			/*
			 *  Secondary master species
			 */
		}
		else
		{
			t = master_ptr->total;
		}
		strcpy(name, master[i]->elt->name);
		sys[count_sys].name = string_duplicate(name);
		sys[count_sys].moles = t;
		sys_tot += sys[count_sys].moles;
		if (master[i]->s->type <= SOLID)
		{
			sys[count_sys].type = string_duplicate("dis");
		}
		else if (master[i]->s->type == EX)
		{
			sys[count_sys].type = string_duplicate("ex");
		}
		else if (master[i]->s->type == SURF || master[i]->s->type == SURF_PSI)
		{
			sys[count_sys].type = string_duplicate("surf");
		}
		count_sys++;
		space((void **) ((void *) &sys), count_sys, &max_sys,
			  sizeof(struct system_species));

	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
system_total_si(void)
/* ---------------------------------------------------------------------- */
{
	int i;
	LDBLE si, iap;
	struct rxn_token *rxn_ptr;
	char name[MAX_LENGTH];

	sys_tot = -999.9;
	for (i = 0; i < count_phases; i++)
	{
		if (phases[i]->in == FALSE || phases[i]->type != SOLID)
			continue;
/*
 *   Print saturation index
 */
		iap = 0.0;
		for (rxn_ptr = phases[i]->rxn_x->token + 1; rxn_ptr->s != NULL;
			 rxn_ptr++)
		{
			iap += rxn_ptr->s->la * rxn_ptr->coef;
		}
		si = -phases[i]->lk + iap;
		strcpy(name, phases[i]->name);
		sys[count_sys].name = string_duplicate(name);
		sys[count_sys].moles = si;
		if (si > sys_tot)
			sys_tot = si;
		sys[count_sys].type = string_duplicate("phase");
		count_sys++;
		space((void **) ((void *) &sys), count_sys, &max_sys,
			  sizeof(struct system_species));
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
system_total_aq(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */
	int i;
/*
 *   find total moles in aq, surface, and exchange
 */
	for (i = 0; i < count_s_x; i++)
	{
		//if (s_x[i]->type != AQ)
		if (s_x[i]->type > HPLUS)
			continue;
		sys[count_sys].name = string_duplicate(s_x[i]->name);
		sys[count_sys].moles = s_x[i]->moles;
		sys_tot += sys[count_sys].moles;
		sys[count_sys].type = string_duplicate("aq");
		count_sys++;
		space((void **) ((void *) &sys), count_sys, &max_sys,
			  sizeof(struct system_species));
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
system_total_ex(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */
	int i;
/*
 *   find total moles in aq, surface, and exchange
 */
	for (i = 0; i < count_s_x; i++)
	{
		if (s_x[i]->type != EX)
			continue;
		if (s_x[i]->primary != NULL)
			continue;
		sys[count_sys].name = string_duplicate(s_x[i]->name);
		sys[count_sys].moles = s_x[i]->moles;
		sys_tot += sys[count_sys].moles;
		sys[count_sys].type = string_duplicate("ex");
		count_sys++;
		space((void **) ((void *) &sys), count_sys, &max_sys,
			  sizeof(struct system_species));
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
system_total_surf(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */
	int i;
/*
 *   find total moles in aq, surface, and exchange
 */
	for (i = 0; i < count_s_x; i++)
	{
		if (s_x[i]->type != SURF)
			continue;
		sys[count_sys].name = string_duplicate(s_x[i]->name);
		sys[count_sys].moles = s_x[i]->moles;
		sys_tot += sys[count_sys].moles;
		sys[count_sys].type = string_duplicate("surf");
		count_sys++;
		space((void **) ((void *) &sys), count_sys, &max_sys,
			  sizeof(struct system_species));
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
system_total_gas(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */
	int i;

/*
 *   find total in gas phase
 */
	if (use.Get_gas_phase_ptr() == NULL)
		return (OK);
	cxxGasPhase *gas_phase_ptr = use.Get_gas_phase_ptr();
	for (size_t j = 0; j < gas_phase_ptr->Get_gas_comps().size(); j++)
	{
		struct phase *phase_ptr = phase_bsearch(gas_phase_ptr->Get_gas_comps()[j].Get_phase_name().c_str(), 
			&i, FALSE);
		assert(phase_ptr);
		sys[count_sys].name = string_duplicate(phase_ptr->name);
		sys[count_sys].moles = phase_ptr->moles_x;
		sys_tot += sys[count_sys].moles;
		sys[count_sys].type = string_duplicate("gas");
		count_sys++;
		space((void **) ((void *) &sys), count_sys, &max_sys,
			  sizeof(struct system_species));
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
system_total_equi(void)
/* ---------------------------------------------------------------------- */
{
/*
 *  Equilibrium phases
 */
	if (use.Get_pp_assemblage_ptr() == NULL)
		return (OK);
	std::map <std::string, cxxPPassemblageComp > comps = use.Get_pp_assemblage_ptr()->Get_pp_assemblage_comps();
	std::map <std::string, cxxPPassemblageComp >::iterator it = comps.begin();
	for  ( ; it != comps.end(); it++)
	{
			cxxPPassemblageComp *comp_ptr = &(it->second);
			int l;
			struct phase *phase_ptr = phase_bsearch(comp_ptr->Get_name().c_str(), &l, FALSE);
			sys[count_sys].name = string_duplicate(phase_ptr->name);
			sys[count_sys].moles = comp_ptr->Get_moles();
			sys_tot += sys[count_sys].moles;
			sys[count_sys].type = string_duplicate("equi");
			count_sys++;
			space((void **) ((void *) &sys), count_sys, &max_sys,
				  sizeof(struct system_species));
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
system_total_ss(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */

/*
 *  Solid solutions
 */
	if (use.Get_ss_assemblage_ptr() == NULL)
		return (OK);
	std::vector<cxxSS *> ss_ptrs = use.Get_ss_assemblage_ptr()->Vectorize();
	for (size_t k = 0; k < ss_ptrs.size(); k++)
	{
		cxxSS *ss_ptr = ss_ptrs[k];
		for (size_t i = 0; i < ss_ptr->Get_ss_comps().size(); i++)
		{
			cxxSScomp *comp_ptr = &(ss_ptr->Get_ss_comps()[i]);
			int l;
			struct phase *phase_ptr = phase_bsearch(comp_ptr->Get_name().c_str(), &l, FALSE);
			sys[count_sys].name = string_duplicate(phase_ptr->name);
			sys[count_sys].moles = comp_ptr->Get_moles();
			sys_tot += sys[count_sys].moles;
			sys[count_sys].type = string_duplicate("s_s");
			count_sys++;
			space((void **) ((void *) &sys), count_sys, &max_sys,
				  sizeof(struct system_species));
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
system_total_elt(const char *total_name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */
	int i, j, k;
	LDBLE molality, moles_excess, moles_surface, mass_water_surface;
	char name[MAX_LENGTH];

/*
 *   find total moles in aq, surface, and exchange
 */
	for (i = 0; i < count_s_x; i++)
	{
		count_elts = 0;
		paren_count = 0;
		add_elt_list(s_x[i]->next_elt, s_x[i]->moles);

		if (count_elts > 0)
		{
			qsort(elt_list, (size_t) count_elts,
				(size_t) sizeof(struct elt_list), elt_list_compare);
			elt_list_combine();
		}
		/*
		 *   Look for element
		 */
		for (j = 0; j < count_elts; j++)
		{
			if (strcmp(elt_list[j].elt->name, total_name) == 0)
			{
				sys[count_sys].name = string_duplicate(s_x[i]->name);
				sys[count_sys].moles = elt_list[j].coef;
				sys_tot += sys[count_sys].moles;
				if (s_x[i]->type == AQ)
				{
					sys[count_sys].type = string_duplicate("aq");
				}
				else if (s_x[i]->type == EX)
				{
					sys[count_sys].type = string_duplicate("ex");
					/* subtract again the dummy moles of primary exchange species... */
					if (s_x[i]->primary != NULL)
					{
						sys_tot -= elt_list[j].coef;
					}
				}
				else if (s_x[i]->type == SURF)
				{
					sys[count_sys].type = string_duplicate("surf");
				}
				else if (s_x[i]->type == HPLUS)
				{
					sys[count_sys].type = string_duplicate("aq");
					/* sys[count_sys].moles = total_h_x; */
				}
				else if (s_x[i]->type == H2O)
				{
					sys[count_sys].type = string_duplicate("aq");
					/* sys[count_sys].moles = total_o_x; */
				}
				else
				{
					error_msg("System_total", STOP);
				}
				count_sys++;
				space((void **) ((void *) &sys), count_sys, &max_sys,
					  sizeof(struct system_species));
				break;
			}
		}
	}
	if (use.Get_surface_ptr() != NULL && dl_type_x != cxxSurface::NO_DL)
	{
		/*
		 *   Find position of component in surface charge data
		 */
		i = -1;
		for (k = 0; k < count_unknowns; k++)
		{
			if (x[k]->type != SURFACE_CB)
				continue;
			cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[k]->surface_charge);
			i++;
			/*
			 *   Loop through all surface components, calculate each H2O surface (diffuse layer),
			 *   H2O aq, and H2O bulk (diffuse layers plus aqueous).
			 */
			mass_water_surface = charge_ptr->Get_mass_water();
			count_elts = 0;
			paren_count = 0;
			for (j = 0; j < count_s_x; j++)
			{
				if (s_x[j]->type > HPLUS)
					continue;
				molality = under(s_x[j]->lm);
				moles_excess =
					mass_water_aq_x * molality *
					(charge_ptr->Get_g_map()[s_x[j]->z].Get_g() * s_x[j]->erm_ddl +
					 mass_water_surface / mass_water_aq_x * (s_x[j]->erm_ddl -
															 1));
				moles_surface = mass_water_surface * molality + moles_excess;
				/*
				 *   Accumulate elements in diffuse layer
				 */
				add_elt_list(s_x[j]->next_elt, moles_surface);
			}
			if (count_elts > 0)
			{
				qsort(elt_list, (size_t) count_elts,
					  (size_t) sizeof(struct elt_list), elt_list_compare);
				elt_list_combine();
			}
			/*
			 *   Print totals
			 */
			for (j = 0; j < count_elts; j++)
			{
				if (strcmp(elt_list[j].elt->name, total_name) == 0)
				{
					strcpy(name, x[k]->master[0]->elt->name);
					replace("_psi", "", name);
					sys[count_sys].name = string_duplicate(name);
					sys[count_sys].moles = elt_list[j].coef;
					sys_tot += sys[count_sys].moles;
					sys[count_sys].type = string_duplicate("diff");
					count_sys++;
					space((void **) ((void *) &sys), count_sys, &max_sys,
						  sizeof(struct system_species));
					break;
				}
			}
		}
	}
/*
 *   find total moles in mineral phases
 */
	if (use.Get_pp_assemblage_in() == TRUE && use.Get_pp_assemblage_ptr() != NULL)
	{
		for (i = 0; i < count_unknowns; i++)
		{
			if (x[i]->type != PP)
				continue;
			//std::map<std::string, cxxPPassemblageComp>::iterator it;
			//it =  pp_assemblage_ptr->Get_pp_assemblage_comps().find(x[i]->pp_assemblage_comp_name);
			cxxPPassemblageComp * comp_ptr = (cxxPPassemblageComp * ) x[i]->pp_assemblage_comp_ptr;
			//if (it->second.Get_add_formula().size() > 0)
			if (comp_ptr->Get_add_formula().size() > 0)
				continue;
			count_elts = 0;
			paren_count = 0;
			int j;
			//struct phase * phase_ptr = phase_bsearch(x[i]->pp_assemblage_comp_name, &j, FALSE);
			struct phase * phase_ptr = x[i]->phase;
			add_elt_list(phase_ptr->next_elt, x[i]->moles);
			if (count_elts > 0)
			{
				qsort(elt_list, (size_t) count_elts,
					  (size_t) sizeof(struct elt_list), elt_list_compare);
				elt_list_combine();
			}
			for (j = 0; j < count_elts; j++)
			{
				if (strcmp(elt_list[j].elt->name, total_name) == 0)
				{
					sys[count_sys].name =
						string_duplicate(phase_ptr->name);
					sys[count_sys].moles = elt_list[j].coef;
					sys_tot += sys[count_sys].moles;
					sys[count_sys].type = string_duplicate("equi");
					count_sys++;
					space((void **) ((void *) &sys), count_sys, &max_sys,
						  sizeof(struct system_species));
					break;
				}
			}
		}
	}
/*
 *  Solid solutions
 */
	if (use.Get_ss_assemblage_ptr() != NULL)
	{
		std::vector<cxxSS *> ss_ptrs = use.Get_ss_assemblage_ptr()->Vectorize();
		for (size_t k = 0; k < ss_ptrs.size(); k++)
		{
			cxxSS *ss_ptr = ss_ptrs[k];
			if (ss_ptr->Get_ss_in())
			{
				for (size_t i = 0; i < ss_ptr->Get_ss_comps().size(); i++)
				{
					cxxSScomp *comp_ptr = &(ss_ptr->Get_ss_comps()[i]);
					int l;
					struct phase *phase_ptr = phase_bsearch(comp_ptr->Get_name().c_str(), &l, FALSE);
					count_elts = 0;
					paren_count = 0;
					add_elt_list(phase_ptr->next_elt,
								 comp_ptr->Get_moles());
					if (count_elts > 0)
					{
						qsort(elt_list, (size_t) count_elts,
							  (size_t) sizeof(struct elt_list),
							  elt_list_compare);
						elt_list_combine();
					}
					for (j = 0; j < count_elts; j++)
					{
						if (strcmp(elt_list[j].elt->name, total_name) == 0)
						{
							sys[count_sys].name =
								string_duplicate(phase_ptr->name);
							sys[count_sys].moles = elt_list[j].coef;
							sys_tot += sys[count_sys].moles;
							sys[count_sys].type = string_duplicate("s_s");
							count_sys++;
							space((void **) ((void *) &sys), count_sys,
								  &max_sys, sizeof(struct system_species));
							break;
						}
					}
				}
			}
		}
	}
/*
 *   find total in gas phase
 */
	if (use.Get_gas_phase_ptr() != NULL)
	{
		cxxGasPhase *gas_phase_ptr = use.Get_gas_phase_ptr();
		for (size_t i = 0; i < gas_phase_ptr->Get_gas_comps().size(); i++)
		{
			struct phase *phase_ptr = 
				phase_bsearch(gas_phase_ptr->Get_gas_comps()[i].Get_phase_name().c_str(), &k, FALSE);
			assert(phase_ptr);
			if (phase_ptr->in == TRUE)
			{
				count_elts = 0;
				paren_count = 0;
				add_elt_list(phase_ptr->next_elt, phase_ptr->moles_x);
				if (count_elts > 0)
				{
					qsort(elt_list, (size_t) count_elts,
						  (size_t) sizeof(struct elt_list), elt_list_compare);
					elt_list_combine();
				}
				/*
				 *   Look for element
				 */
				for (j = 0; j < count_elts; j++)
				{
					if (strcmp(elt_list[j].elt->name, total_name) == 0)
					{
						sys[count_sys].name = string_duplicate(phase_ptr->name);
						sys[count_sys].moles = elt_list[j].coef;
						sys_tot += sys[count_sys].moles;
						sys[count_sys].type = string_duplicate("gas");
						count_sys++;
						space((void **) ((void *) &sys), count_sys, &max_sys,
							  sizeof(struct system_species));
						break;
					}
				}
			}
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
system_total_elt_secondary(const char *total_name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */
	int i, j, k, l;
	LDBLE molality, moles_excess, moles_surface, mass_water_surface, sum,
		coef;
	char name[MAX_LENGTH];
/*
 *   find total moles in aq, surface, and exchange
 */
	for (i = 0; i < count_s_x; i++)
	{
		count_elts = 0;
		paren_count = 0;
		if (s_x[i]->next_secondary != NULL)
		{
			add_elt_list(s_x[i]->next_secondary, s_x[i]->moles);
		}
		else
		{
			add_elt_list(s_x[i]->next_sys_total, s_x[i]->moles);
		}
		if (count_elts > 0)
		{
			qsort(elt_list, (size_t) count_elts,
				  (size_t) sizeof(struct elt_list), elt_list_compare);
			elt_list_combine();
		}
		/*
		 *   Look for element
		 */
		for (j = 0; j < count_elts; j++)
		{
			if (strcmp(elt_list[j].elt->name, total_name) == 0)
			{
				sys[count_sys].name = string_duplicate(s_x[i]->name);
				sys[count_sys].moles = elt_list[j].coef;
				sys_tot += sys[count_sys].moles;
				if (s_x[i]->type == AQ)
				{
					sys[count_sys].type = string_duplicate("aq");
				}
				else if (s_x[i]->type == EX)
				{
					sys[count_sys].type = string_duplicate("ex");
				}
				else if (s_x[i]->type == SURF)
				{
					sys[count_sys].type = string_duplicate("surf");
				}
				else if (s_x[i]->type == HPLUS)
				{
					sys[count_sys].type = string_duplicate("aq");
					/* sys[count_sys].moles = total_h_x; */
				}
				else if (s_x[i]->type == H2O)
				{
					sys[count_sys].type = string_duplicate("aq");
					/* sys[count_sys].moles = total_o_x; */
				}
				else
				{
					error_msg("System_total", STOP);
				}
				count_sys++;
				space((void **) ((void *) &sys), count_sys, &max_sys,
					  sizeof(struct system_species));
				break;
			}
		}
	}
	if (use.Get_surface_ptr() != NULL && dl_type_x != cxxSurface::NO_DL)
	{
		/*
		 *   Find position of component in surface charge data
		 */
		i = -1;
		for (k = 0; k < count_unknowns; k++)
		{
			if (x[k]->type != SURFACE_CB)
				continue;
			cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[k]->surface_charge);
			i++;
			/*
			 *   Loop through all surface components, calculate each H2O surface (diffuse layer),
			 *   H2O aq, and H2O bulk (diffuse layers plus aqueous).
			 */
			mass_water_surface = charge_ptr->Get_mass_water();
			sum = 0;
			for (j = 0; j < count_s_x; j++)
			{
				count_elts = 0;
				paren_count = 0;
				if (s_x[i]->next_secondary != NULL)
				{
					add_elt_list(s_x[i]->next_secondary, 1);
				}
				else
				{
					add_elt_list(s_x[i]->next_sys_total, 1);
				}
				for (l = 0; l < count_elts; l++)
				{
					if (strcmp(elt_list[l].elt->name, total_name) == 0)
					{
						coef = elt_list[l].coef;
						if (s_x[j]->type > H2O)
							continue;
						molality = under(s_x[j]->lm);
						moles_excess =
							mass_water_aq_x * molality *
							charge_ptr->Get_g_map()[s_x[j]->z].Get_g();
						moles_surface =
							mass_water_surface * molality + moles_excess;
						sum += moles_surface * coef;
						break;
					}
				}
				if (l >= count_elts)
					continue;
				strcpy(name, x[k]->master[0]->elt->name);
				replace("_psi", "", name);
				sys[count_sys].name = string_duplicate(name);
				sys[count_sys].moles = sum;
				sys_tot += sys[count_sys].moles;
				sys[count_sys].type = string_duplicate("diff");
				count_sys++;
				space((void **) ((void *) &sys), count_sys, &max_sys,
					  sizeof(struct system_species));
				break;
			}
		}
	}
/*
 *   find total moles in mineral phases
 */
	if (use.Get_pp_assemblage_in() == TRUE && use.Get_pp_assemblage_ptr() != NULL)
	{
		for (i = 0; i < count_unknowns; i++)
		{
			if (x[i]->type != PP)
				continue;
			//std::map<std::string, cxxPPassemblageComp>::iterator it;
			//it =  pp_assemblage_ptr->Get_pp_assemblage_comps().find(x[i]->pp_assemblage_comp_name);
			cxxPPassemblageComp * comp_ptr = (cxxPPassemblageComp * ) x[i]->pp_assemblage_comp_ptr;
			//if (it->second.Get_add_formula().size() > 0)
			if (comp_ptr->Get_add_formula().size() > 0)
				continue;
			count_elts = 0;
			paren_count = 0;
			int j;
			//struct phase * phase_ptr = phase_bsearch(x[i]->pp_assemblage_comp_name, &j, FALSE);
			struct phase * phase_ptr = x[i]->phase;
			add_elt_list(phase_ptr->next_sys_total,	 x[i]->moles);
			if (count_elts > 0)
			{
				qsort(elt_list, (size_t) count_elts,
					  (size_t) sizeof(struct elt_list), elt_list_compare);
				elt_list_combine();
			}
			for (j = 0; j < count_elts; j++)
			{
				if (strcmp(elt_list[j].elt->name, total_name) == 0)
				{
					sys[count_sys].name =
						string_duplicate(phase_ptr->name);
					sys[count_sys].moles = elt_list[j].coef;
					sys_tot += sys[count_sys].moles;
					sys[count_sys].type = string_duplicate("equi");
					count_sys++;
					space((void **) ((void *) &sys), count_sys, &max_sys,
						  sizeof(struct system_species));
					break;
				}
			}
		}
	}
/*
 *  Solid solutions
 */
	if (use.Get_ss_assemblage_ptr() != NULL)
	{
		std::vector<cxxSS *> ss_ptrs = use.Get_ss_assemblage_ptr()->Vectorize();
		for (size_t i = 0; i < ss_ptrs.size(); i++)
		{
			cxxSS *ss_ptr = ss_ptrs[i];
			if (ss_ptr->Get_ss_in())
			{
				for (size_t k = 0; k < ss_ptr->Get_ss_comps().size(); k++)
				{
					cxxSScomp *comp_ptr = &(ss_ptr->Get_ss_comps()[k]);
					int l;
					struct phase *phase_ptr = phase_bsearch(comp_ptr->Get_name().c_str(), &l, FALSE);
					count_elts = 0;
					paren_count = 0;
					add_elt_list(phase_ptr->next_sys_total,
								 comp_ptr->Get_moles());
					if (count_elts > 0)
					{
						qsort(elt_list, (size_t) count_elts,
							  (size_t) sizeof(struct elt_list),
							  elt_list_compare);
						elt_list_combine();
					}
					for (j = 0; j < count_elts; j++)
					{
						if (strcmp(elt_list[j].elt->name, total_name) == 0)
						{
							sys[count_sys].name =
								string_duplicate(phase_ptr->name);
							sys[count_sys].moles = elt_list[j].coef;
							sys_tot += sys[count_sys].moles;
							sys[count_sys].type = string_duplicate("s_s");
							count_sys++;
							space((void **) ((void *) &sys), count_sys,
								  &max_sys, sizeof(struct system_species));
							break;
						}
					}
				}
			}
		}
	}
/*
 *   find total in gas phase
 */
	if (use.Get_gas_phase_ptr() != NULL)
	{
		cxxGasPhase *gas_phase_ptr = use.Get_gas_phase_ptr();
		for (size_t j = 0; j < gas_phase_ptr->Get_gas_comps().size(); j++)	
		{
			struct phase *phase_ptr = 
				phase_bsearch(gas_phase_ptr->Get_gas_comps()[j].Get_phase_name().c_str(), &i, FALSE);
			assert(phase_ptr);
			if (phase_ptr->in == TRUE)
			{
				count_elts = 0;
				paren_count = 0;
				add_elt_list(phase_ptr->next_sys_total,
							 phase_ptr->moles_x);

				if (count_elts > 0)
				{
					qsort(elt_list, (size_t) count_elts,
						  (size_t) sizeof(struct elt_list), elt_list_compare);
					elt_list_combine();
				}
				/*
				 *   Look for element
				 */
				for (size_t j1 = 0; j1 < (size_t) count_elts; j1++)
				{
					if (strcmp(elt_list[j1].elt->name, total_name) == 0)
					{
						sys[count_sys].name =
							string_duplicate(phase_ptr->name);
						sys[count_sys].moles = elt_list[j1].coef;
						sys_tot += sys[count_sys].moles;
						sys[count_sys].type = string_duplicate("gas");
						count_sys++;
						space((void **) ((void *) &sys), count_sys, &max_sys,
							  sizeof(struct system_species));
						break;
					}
				}
			}
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
solution_number(void)
/* ---------------------------------------------------------------------- */
{
	Phreeqc * PhreeqcPtr = this;
	int soln_no = -999;
	if (PhreeqcPtr->state == TRANSPORT)
	{
		soln_no = PhreeqcPtr->cell_no;
	}
	else if (PhreeqcPtr->state == PHAST)
	{
		soln_no = PhreeqcPtr->cell_no;
	}
	else if (PhreeqcPtr->state == ADVECTION)
	{
		soln_no = PhreeqcPtr->cell_no;
	}
	else if (PhreeqcPtr->state < REACTION)
	{
		soln_no = PhreeqcPtr->use.Get_solution_ptr()->Get_n_user();
	}
	else
	{
		if (PhreeqcPtr->use.Get_mix_in())
		{
			soln_no = PhreeqcPtr->use.Get_n_mix_user();
		}
		else
		{
			soln_no = PhreeqcPtr->use.Get_n_solution_user();
		}
	}
	return soln_no;
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
solution_sum_secondary(const char *total_name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */
	int i, j;
	LDBLE sum;
/*
 *   find total moles in aq, surface, and exchange
 */
	sum = 0;
	for (i = 0; i < count_s_x; i++)
	{
		if (s_x[i]->type > H2O)
			continue;
		count_elts = 0;
		paren_count = 0;
		if (s_x[i]->next_secondary != NULL)
		{
			add_elt_list(s_x[i]->next_secondary, s_x[i]->moles);
		}
		else
		{
			add_elt_list(s_x[i]->next_sys_total, s_x[i]->moles);
		}
		if (count_elts > 0)
		{
			qsort(elt_list, (size_t) count_elts,
				  (size_t) sizeof(struct elt_list), elt_list_compare);
			elt_list_combine();
		}
		/*
		 *   Look for element
		 */
		for (j = 0; j < count_elts; j++)
		{
			if (strcmp(elt_list[j].elt->name, total_name) == 0)
			{
				sum += elt_list[j].coef;
				break;
			}
		}
	}
	return (sum);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
system_species_compare(const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
	const struct system_species *a, *b;

	a = (const struct system_species *) ptr1;
	b = (const struct system_species *) ptr2;
	if (a->moles < b->moles)
		return (1);
	if (a->moles > b->moles)
		return (-1);
	return (0);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
system_total_solids(cxxExchange *exchange_ptr,
					cxxPPassemblage *pp_assemblage_ptr,
					cxxGasPhase *gas_phase_ptr,
					cxxSSassemblage *ss_assemblage_ptr,
					cxxSurface *surface_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in solid phases
 */
	count_elts = 0;
	paren_count = 0;
/*
 *   find total moles in exchanger
 */
	if (exchange_ptr != NULL)
	{
		for (size_t i = 0; i < exchange_ptr->Get_exchange_comps().size(); i++)
		{
			add_elt_list(exchange_ptr->Get_exchange_comps()[i].Get_totals(), 1.0);
		}
	}
	if (surface_ptr != NULL)
	{
		for (size_t i = 0; i < surface_ptr->Get_surface_comps().size(); i++)
		{
			add_elt_list(surface_ptr->Get_surface_comps()[i].Get_totals(), 1.0);
		}
	}
	if (ss_assemblage_ptr != NULL)
	{
		std::vector<cxxSS *> ss_ptrs = ss_assemblage_ptr->Vectorize();
		for (size_t i = 0; i < ss_ptrs.size(); i++)
		{
			cxxSS *ss_ptr = ss_ptrs[i];
			for (size_t j = 0; j < ss_ptr->Get_ss_comps().size(); j++)
			{
				cxxSScomp *comp_ptr = &(ss_ptr->Get_ss_comps()[j]);
				int l;
				struct phase *phase_ptr = phase_bsearch(comp_ptr->Get_name().c_str(), &l, FALSE);
				add_elt_list(phase_ptr->next_elt,
							 comp_ptr->Get_moles());
			}
		}
	}
	if (gas_phase_ptr != NULL)
	{
		for (size_t j = 0; j < gas_phase_ptr->Get_gas_comps().size(); j++)
		{
			int i;
			struct phase *phase_ptr = 
				phase_bsearch(gas_phase_ptr->Get_gas_comps()[j].Get_phase_name().c_str(), &i, FALSE);
			add_elt_list(phase_ptr->next_elt, gas_phase_ptr->Get_gas_comps()[j].Get_moles());
		}
	}
	if (pp_assemblage_ptr != NULL)
	{
		std::map<std::string, cxxPPassemblageComp>::iterator it;
		it =  pp_assemblage_ptr->Get_pp_assemblage_comps().begin();
		for ( ; it != pp_assemblage_ptr->Get_pp_assemblage_comps().end(); it++)
		{
			int j;
			struct phase * phase_ptr = phase_bsearch(it->first.c_str(), &j, FALSE);
			add_elt_list(phase_ptr->next_elt,
						 it->second.Get_moles());
		}
	}

	if (count_elts > 0)
	{
		qsort(elt_list, (size_t) count_elts,
			  (size_t) sizeof(struct elt_list), elt_list_compare);
		elt_list_combine();
	}
	return (OK);
}

LDBLE Phreeqc::
iso_value(const char *total_name)
{
	int j;
	char token[MAX_LENGTH];
	char my_total_name[MAX_LENGTH];
	strcpy(token, "");
	strcpy(my_total_name, total_name);
	while (replace(" ","_",my_total_name));
	for (j = 0; j < count_isotope_ratio; j++)
	{
		if (isotope_ratio[j]->ratio == MISSING)
			continue;
		if (strcmp(my_total_name, isotope_ratio[j]->name) != 0)
			continue;
		return (isotope_ratio[j]->converted_ratio);
	}
	strcpy(my_total_name, total_name);
	while (replace("[","",my_total_name));
	while (replace("]","",my_total_name));
	strcat(token,"R(");
	strcat(token,my_total_name);
	strcat(token,")");
	for (j = 0; j < count_isotope_ratio; j++)
	{
		if (isotope_ratio[j]->ratio == MISSING)
			continue;
		if (strcmp(token, isotope_ratio[j]->name) != 0)
			continue;
		return (isotope_ratio[j]->converted_ratio);
	}
	return -1000.;
}

char * Phreeqc::
iso_unit(const char *total_name)
{
	int j;
	char token[MAX_LENGTH], unit[MAX_LENGTH];
	struct master_isotope *master_isotope_ptr;
	char my_total_name[MAX_LENGTH];
	strcpy(token, "");
	strcpy(my_total_name, total_name);
	while (replace(" ","_",my_total_name));
	strcpy(unit, "unknown");
	for (j = 0; j < count_isotope_ratio; j++)
	{
		if (isotope_ratio[j]->ratio == MISSING)
			continue;
		if (strcmp(my_total_name, isotope_ratio[j]->name) != 0)
			continue;
		master_isotope_ptr = master_isotope_search(isotope_ratio[j]->isotope_name);
		if (master_isotope_ptr != NULL)
		{
			strcpy(unit, master_isotope_ptr->units);
		}
		return string_duplicate(unit);
	}
	strcpy(my_total_name, total_name);
	while (replace("[","",my_total_name));
	while (replace("]","",my_total_name));
	strcat(token,"R(");
	strcat(token,my_total_name);
	strcat(token,")");
	for (j = 0; j < count_isotope_ratio; j++)
	{
		if (isotope_ratio[j]->ratio == MISSING)
			continue;
		if (strcmp(token, isotope_ratio[j]->name) != 0)
			continue;
		master_isotope_ptr = master_isotope_search(isotope_ratio[j]->isotope_name);
		if (master_isotope_ptr != NULL)
		{
			strcpy(unit, master_isotope_ptr->units);
		}
		return string_duplicate(unit);
	}
	return string_duplicate(unit);
}

int Phreeqc::
basic_compile(char *commands, void **lnbase, void **vbase, void **lpbase)
{
	return this->basic_interpreter->basic_compile(commands, lnbase, vbase, lpbase);
}

int Phreeqc::
basic_run(char *commands, void *lnbase, void *vbase, void *lpbase)
{
	return this->basic_interpreter->basic_run(commands, lnbase, vbase, lpbase);
}

void Phreeqc::
basic_free(void)
{
	delete this->basic_interpreter;
}

#if defined(SWIG) || defined(SWIG_IPHREEQC)

#include "BasicCallback.h"

double Phreeqc::
basic_callback(double x1, double x2, const char * str)
{
	if (this->basicCallback)
	{
		return this->basicCallback->Callback(x1, x2, str);
	}
	return 0.0;
}

#else  /* defined(SWIG) || defined(SWIG_IPHREEQC) */

#ifdef IPHREEQC_NO_FORTRAN_MODULE
double Phreeqc::
basic_callback(double x1, double x2, char * str)
#else
double Phreeqc::
basic_callback(double x1, double x2, const char * str)
#endif
{
	double local_x1 = x1;
	double local_x2 = x2;

	if (basic_callback_ptr != NULL)
	{
		return (*basic_callback_ptr) (x1, x2, (const char *) str, basic_callback_cookie);
	}
	if (basic_fortran_callback_ptr != NULL)
	{
#ifdef IPHREEQC_NO_FORTRAN_MODULE
		return (*basic_fortran_callback_ptr) (&local_x1, &local_x2, str, (int) strlen(str));
#else
		return (*basic_fortran_callback_ptr) (&local_x1, &local_x2, str, (int) strlen(str));
#endif
	}
	return 0;
}

void 
Phreeqc::register_basic_callback(double (*fcn)(double x1, double x2, const char *str, void *cookie), void *cookie1)
{
	this->basic_callback_ptr = fcn;
	this->basic_callback_cookie = cookie1;
}
#ifdef IPHREEQC_NO_FORTRAN_MODULE
void 
Phreeqc::register_fortran_basic_callback(double ( *fcn)(double *x1, double *x2, char *str, size_t l))
{
	this->basic_fortran_callback_ptr = fcn;
}
#else

void 
Phreeqc::register_fortran_basic_callback(double ( *fcn)(double *x1, double *x2, const char *str, int l))
{
	this->basic_fortran_callback_ptr = fcn;
}
#endif

#endif  /* defined(SWIG) || defined(SWIG_IPHREEQC) */
