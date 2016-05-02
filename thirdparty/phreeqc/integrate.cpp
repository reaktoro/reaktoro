#include "Phreeqc.h"
#include "phqalloc.h"
#include "Utils.h"
#include "Solution.h"

#define MAX_QUAD 20
#define K_POLY 5

/* ---------------------------------------------------------------------- */
int Phreeqc::
calc_all_g(void)
/* ---------------------------------------------------------------------- */
{
	int converge, converge1;
	LDBLE new_g, xd1;
	LDBLE epsilon;

	if (use.Get_surface_ptr() == NULL)
		return (OK);
/*
 *   calculate g for each surface
 */
	epsilon = convergence_tolerance;
	if (convergence_tolerance >= 1e-8)
	{
		G_TOL = 1e-9;
	}
	else
	{
		G_TOL = 1e-10;
	}

	converge = TRUE;

	for (int j = 0; j < count_unknowns; j++)
	{
		if (x[j]->type != SURFACE_CB)
			continue;
		if (debug_diffuse_layer == TRUE)
			output_msg(sformatf( "Calc_all_g, X[%d]\n", j));
		cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[j]->surface_charge);
		std::map<LDBLE, cxxSurfDL> temp_g_map;
		cxxSurfDL temp_g;
		charge_ptr->Get_g_map()[0] = temp_g;
		temp_g_map[0] = temp_g;
		xd_global = exp(-2 * x[j]->master[0]->s->la * LOG_10);
		/* alpha = 0.02935 @ 25;        (ee0RT/2)**1/2, (L/mol)**1/2 C / m**2 */
		/* 1000 J/kJ and 1000 L/m**3 */
		//alpha_global =	sqrt(EPSILON * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000.0) * 1000.0 *
		//		 tk_x * 0.5);
		alpha_global =	sqrt(eps_r * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000.0) * 1000.0 *
				 tk_x * 0.5);
/*
 *   calculate g for given surface for each species
 */
		for (int i = 0; i < count_s_x; i++)
		{
			if (s_x[i]->type > HPLUS)
				continue;
			if (temp_g_map.find(s_x[i]->z) != temp_g_map.end())
				continue;
			z_global = s_x[i]->z;
			if (charge_ptr->Get_grams() > 0.0)
			{
				
				if ((use.Get_surface_ptr()->Get_only_counter_ions() == false) ||
					(((x[j]->master[0]->s->la > 0) && (z_global < 0))
					 || ((x[j]->master[0]->s->la < 0) && (z_global > 0))))
				{
					if (xd_global > 0.1)
					{
						new_g = qromb_midpnt(charge_ptr, 1.0, xd_global);
					}
					else if (xd_global > 0.01)
					{
						new_g = qromb_midpnt(charge_ptr, 1.0, 0.1);
						new_g += qromb_midpnt(charge_ptr, 0.1, xd_global);
					}
					else if (xd_global > 0.001)
					{
						new_g = qromb_midpnt(charge_ptr, 1.0, 0.1);
						new_g += qromb_midpnt(charge_ptr, 0.1, 0.01);
						new_g += qromb_midpnt(charge_ptr, 0.01, xd_global);
					}
					else if (xd_global > 0.0001)
					{
						new_g = qromb_midpnt(charge_ptr, 1.0, 0.1);
						new_g += qromb_midpnt(charge_ptr, 0.1, 0.01);
						new_g += qromb_midpnt(charge_ptr, 0.01, .001);
						new_g += qromb_midpnt(charge_ptr, 0.001, xd_global);
					}
					else if (xd_global > 0.00001)
					{
						new_g = qromb_midpnt(charge_ptr, 1.0, 0.1);
						new_g += qromb_midpnt(charge_ptr, 0.1, 0.01);
						new_g += qromb_midpnt(charge_ptr, 0.01, .001);
						new_g += qromb_midpnt(charge_ptr, 0.001, .0001);
						new_g += qromb_midpnt(charge_ptr, 0.0001, xd_global);
					}
					else if (xd_global > 0.000001)
					{
						new_g = qromb_midpnt(charge_ptr, 1.0, 0.1);
						new_g += qromb_midpnt(charge_ptr, 0.1, 0.01);
						new_g += qromb_midpnt(charge_ptr, 0.01, .001);
						new_g += qromb_midpnt(charge_ptr, 0.001, .0001);
						new_g += qromb_midpnt(charge_ptr, 0.0001, .00001);
						new_g += qromb_midpnt(charge_ptr, 0.00001, xd_global);
					}
					else if (xd_global > 0.0000001)
					{
						new_g = qromb_midpnt(charge_ptr, 1.0, 0.1);
						new_g += qromb_midpnt(charge_ptr, 0.1, 0.01);
						new_g += qromb_midpnt(charge_ptr, 0.01, .001);
						new_g += qromb_midpnt(charge_ptr, 0.001, .0001);
						new_g += qromb_midpnt(charge_ptr, 0.0001, .00001);
						new_g += qromb_midpnt(charge_ptr, 0.00001, .000001);
						new_g += qromb_midpnt(charge_ptr, 0.000001, xd_global);
					}
					else if (xd_global > 0.00000001)
					{
						new_g = qromb_midpnt(charge_ptr, 1.0, 0.1);
						new_g += qromb_midpnt(charge_ptr, 0.1, 0.01);
						new_g += qromb_midpnt(charge_ptr, 0.01, .001);
						new_g += qromb_midpnt(charge_ptr, 0.001, .0001);
						new_g += qromb_midpnt(charge_ptr, 0.0001, .00001);
						new_g += qromb_midpnt(charge_ptr, 0.00001, .000001);
						new_g += qromb_midpnt(charge_ptr, 0.000001, .0000001);
						new_g += qromb_midpnt(charge_ptr, 0.0000001, xd_global);
					}
					else
					{
						new_g = qromb_midpnt(charge_ptr, 1.0, 0.1);
						new_g += qromb_midpnt(charge_ptr, 0.1, 0.01);
						new_g += qromb_midpnt(charge_ptr, 0.01, .001);
						new_g += qromb_midpnt(charge_ptr, 0.001, .0001);
						new_g += qromb_midpnt(charge_ptr, 0.0001, .00001);
						new_g += qromb_midpnt(charge_ptr, 0.00001, .000001);
						new_g += qromb_midpnt(charge_ptr, 0.000001, .0000001);
						new_g += qromb_midpnt(charge_ptr, 0.0000001, .00000001);
						new_g += qromb_midpnt(charge_ptr, 0.00000001, xd_global);
					}
				}
				else
				{
					new_g = 0;
				}
			}
			else
			{
				new_g = 0.0;
			}
			if ((use.Get_surface_ptr()->Get_only_counter_ions()) && new_g < 0)
				new_g = 0;
			converge1 = TRUE;
			if (fabs(new_g) >= 1.)
			{
				if (fabs((new_g - charge_ptr->Get_g_map()[z_global].Get_g()) / new_g) > epsilon)
				{
					converge1 = FALSE;
				}
			}
			else
			{
				if (fabs(new_g - charge_ptr->Get_g_map()[z_global].Get_g()) > epsilon)
				{
					converge1 = FALSE;
				}
			}
			if (converge1 == FALSE)
			{
				converge = FALSE;
				if (debug_diffuse_layer == TRUE)
				{
					output_msg(sformatf(
							   "\t%12f\t%12.4e\t%12.4e\t%12.4e\n",
							   (double) z_global,
							   (double) charge_ptr->Get_g_map()[z_global].Get_g(),
							   (double) new_g,
							   (double) (new_g - charge_ptr->Get_g_map()[z_global].Get_g())));
				}
			}
			charge_ptr->Get_g_map()[z_global].Set_g(new_g);
			if (new_g == 0)
			{
				charge_ptr->Get_g_map()[z_global].Set_dg(0.);
			}
			else
			{
				if (charge_ptr->Get_grams() > 0.0)
				{
					LDBLE dg = charge_ptr->Get_grams() *
						charge_ptr->Get_specific_area() * alpha_global *
						g_function(xd_global) / F_C_MOL;
					dg *=
						-2. / (exp(x[j]->master[0]->s->la * LOG_10) *
							   exp(x[j]->master[0]->s->la * LOG_10));
					if ((xd_global - 1) < 0.0)
					{
						dg *= -1.0;
					}

					if (fabs(dg) < 1e-8)
					{
						xd1 = exp(-2 * 1e-3 * LOG_10);
						new_g = qromb_midpnt(charge_ptr, 1.0, xd1);
						dg = new_g / .001;
					}
					charge_ptr->Get_g_map()[z_global].Set_dg(dg);
				}
				else
				{
					charge_ptr->Get_g_map()[z_global].Set_dg(0.0);
				}
			}
			temp_g_map[z_global] = charge_ptr->Get_g_map()[z_global];
		}
		if (debug_diffuse_layer == TRUE)
		{
			output_msg(sformatf("\nSurface component %d: charge,\tg,\tdg/dlny,\txd\n",
					  (int) charge_ptr->Get_g_map().size()));
			std::map<LDBLE, cxxSurfDL>::iterator it;
			for (it = charge_ptr->Get_g_map().begin(); it != charge_ptr->Get_g_map().end(); it++)
			{
				output_msg(sformatf(
						   "\t%12f\t%12.4e\t%12.4e\t%12.4e\n",
						   (double) it->first,
						   (double) it->second.Get_g(),
						   (double) it->second.Get_dg(),
						   (double) xd_global));
			}
		}
	}
	return (converge);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
g_function(LDBLE x_value)
/* ---------------------------------------------------------------------- */
{
	LDBLE sum, return_value, sum1;
	int i;
	LDBLE ln_x_value;

	if (equal(x_value, 1.0, G_TOL * 100) == TRUE)
		return (0.0);
	sum = 0.0;
	ln_x_value = log(x_value);
	
	cxxSurfaceCharge *charge_ptr = &(use.Get_surface_ptr()->Get_surface_charges()[0]);
	std::map<LDBLE, cxxSurfDL>::iterator it = charge_ptr->Get_g_map().begin();
	for ( ; it != charge_ptr->Get_g_map().end(); it++)
	{
		it->second.Set_psi_to_z(exp(ln_x_value * it->first) - 1.0);
	}
	for (i = 0; i < count_s_x; i++)
	{
		if (s_x[i]->type < H2O && s_x[i]->z != 0.0)
		{
			sum += s_x[i]->moles * charge_ptr->Get_g_map()[s_x[i]->z].Get_psi_to_z();
		}
	}
	if (sum < 0.0)
	{
		sum = 0.0;
		sum1 = 0.0;
		output_msg(sformatf(
				   "Species\tmoles\tX**z-1\tsum\tsum charge\n"));
		for (i = 0; i < count_s_x; i++)
		{
			if (s_x[i]->type < H2O && s_x[i]->z != 0.0)
			{
				sum += s_x[i]->moles * (pow(x_value, s_x[i]->z) - 1.0);
				sum1 += s_x[i]->moles * s_x[i]->z;
				output_msg(sformatf( "%s\t%e\t%e\t%e\t%e\n",
						   s_x[i]->name, (double) s_x[i]->moles,
						   (double) (pow((LDBLE) x_value, (LDBLE) s_x[i]->z) -
									 1.0), (double) sum, (double) sum1));
			}
		}
		error_string = sformatf( "Negative sum in g_function, %e\t%e.",
				(double) sum, (double) x_value);
		error_msg(error_string, CONTINUE);
		error_string = sformatf(
				"Solutions must be charge balanced, charge imbalance is %e\n",
				(double) sum1);
		error_msg(error_string, STOP);
	}

	return_value =
		(exp(ln_x_value * z_global) -
		 1) / sqrt((x_value * x_value * mass_water_aq_x * sum));
	return (return_value);
}
/* ---------------------------------------------------------------------- */
void Phreeqc::
polint(LDBLE * xa, LDBLE * ya, int n, LDBLE xv, LDBLE * yv, LDBLE * dy)
/* ---------------------------------------------------------------------- */
{
	int i, m, ns;
	LDBLE den, dif, dift, ho, hp, w;
	LDBLE *c, *d;

	ns = 1;
	dif = fabs(xv - xa[1]);
/*
 *   Malloc work space
 */
	c = (LDBLE *) PHRQ_malloc((size_t) (n + 1) * sizeof(LDBLE));
	if (c == NULL)
		malloc_error();
	d = (LDBLE *) PHRQ_malloc((size_t) (n + 1) * sizeof(LDBLE));
	if (d == NULL)
		malloc_error();



	for (i = 1; i <= n; i++)
	{
		dift = fabs(xv - xa[i]);
		if (dift < dif)
		{
			ns = i;
			dif = dift;
		}
		c[i] = ya[i];
		d[i] = ya[i];
	}

	*yv = ya[ns--];
	for (m = 1; m < n; m++)
	{
		for (i = 1; i <= n - m; i++)
		{
			ho = xa[i] - xv;
			hp = xa[i + m] - xv;
			w = c[i + 1] - d[i];
			if ((den = ho - hp) == 0.0)
			{
				error_msg("In subroutine polint.", STOP);
			}
			den = w / den;
			d[i] = hp * den;
			c[i] = ho * den;
		}
		if (2 * ns < (n - m))
		{
			*dy = c[ns + 1];
		}
		else
		{
			*dy = d[ns--];
		}
		*yv += *dy;

/*		*yv += (*dy = (2 * ns < (n-m) ? c[ns+1] : d[ns--])); */
	}
	c = (LDBLE *) free_check_null(c);
	d = (LDBLE *) free_check_null(d);
	return;
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
midpnt(LDBLE x1, LDBLE x2, int n)
/* ---------------------------------------------------------------------- */
{
	LDBLE xv, tnm, sum, del, ddel;
	int it, j;

	if (n == 1)
	{
		midpoint_sv = (x2 - x1) * g_function(0.5 * (x1 + x2));
		return (midpoint_sv);
	}
	else
	{
		for (it = 1, j = 1; j < n - 1; j++)
			it *= 3;
		tnm = (LDBLE) it;
		del = (x2 - x1) / (3 * tnm);
		ddel = del + del;
		xv = x1 + 0.5 * del;
		sum = 0.0;
		for (j = 1; j <= it; j++)
		{
#if defined(PHREEQCI_GUI)
			PhreeqcIWait(this);
#endif
			sum += g_function(xv);
			xv += ddel;
			sum += g_function(xv);
			xv += del;
		}
		midpoint_sv = (midpoint_sv + (x2 - x1) * sum / tnm) / 3.0;
		return midpoint_sv;
	}
}

/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
qromb_midpnt(cxxSurfaceCharge *charge_ptr, LDBLE x1, LDBLE x2)
/* ---------------------------------------------------------------------- */
{
	LDBLE ss, dss;
	LDBLE sv[MAX_QUAD + 2], h[MAX_QUAD + 2];
	int j;

	h[0] = 1.0;
	sv[0] = midpnt(x1, x2, 1);
	for (j = 1; j < MAX_QUAD; j++)
	{
		sv[j] = midpnt(x1, x2, j + 1);
		h[j] = h[j - 1] / 9.0;

		if (fabs(sv[j] - sv[j - 1]) <= G_TOL * fabs(sv[j]))
		{
			sv[j] *= charge_ptr->Get_grams() * charge_ptr->Get_specific_area() * alpha_global / F_C_MOL;	/* (ee0RT/2)**1/2, (L/mol)**1/2 C / m**2 */
			if ((x2 - 1) < 0.0)
				sv[j] *= -1.0;
			if (debug_diffuse_layer == TRUE)
			{
				output_msg(sformatf(
						   "Iterations in qromb_midpnt: %d\n", j));
			}
			return (sv[j]);
		}

		if (j >= K_POLY - 1)
		{
			polint(&h[j - K_POLY], &sv[j - K_POLY], K_POLY, 0.0, &ss, &dss);
			if (fabs(dss) <= G_TOL * fabs(ss) || fabs(dss) < G_TOL)
			{
				ss *= charge_ptr->Get_grams() * charge_ptr->Get_specific_area() * alpha_global / F_C_MOL;	/* (ee0RT/2)**1/2, (L/mol)**1/2 C / m**2 */
				if ((x2 - 1) < 0.0)
					ss *= -1.0;
				if (debug_diffuse_layer == TRUE)
				{
					output_msg(sformatf(
							   "Iterations in qromb_midpnt: %d\n", j));
				}
				return (ss);
			}
		}

	}
	error_string = sformatf(
			"\nToo many iterations integrating diffuse layer.\n");
	error_msg(error_string, STOP);
	return (-999.9);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
calc_init_g(void)
/* ---------------------------------------------------------------------- */
{
	if (use.Get_surface_ptr() == NULL)
		return (OK);
/*
 *   calculate g for each surface
 */
	for (int j = 0; j < count_unknowns; j++)
	{
		if (x[j]->type != SURFACE_CB)
			continue;
		cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[j]->surface_charge);
		xd_global = exp(-2 * x[j]->master[0]->s->la * LOG_10);
		/* alpha = 0.02935 @ 25;      (ee0RT/2)**1/2, (L/mol)**1/2 C / m**2 */
		/*  second 1000 is liters/m**3 */
		//alpha_global =	sqrt(EPSILON * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000.0) *
		//	1000.0 * tk_x * 0.5);
		alpha_global =	sqrt(eps_r * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000.0) *
			1000.0 * tk_x * 0.5);

		if (charge_ptr->Get_g_map().size() == 0)
		{
			cxxSurfDL temp_g;
			charge_ptr->Get_g_map()[0.0] = temp_g;
		}
/*
 *   calculate g for given surface for each species
 */
		for (int i = 0; i < count_s_x; i++)
		{
			if (s_x[i]->type > HPLUS)
				continue;

			if (charge_ptr->Get_g_map().find(s_x[i]->z) == charge_ptr->Get_g_map().end())
			{
				cxxSurfDL temp_g;
				/* save g for charge */
				if (charge_ptr->Get_grams() > 0.0)
				{
					temp_g.Set_g(2 * alpha_global * sqrt(mu_x) * (pow(xd_global, s_x[i]->z / 2.0) - 1) *
						charge_ptr->Get_grams() *
						charge_ptr->Get_specific_area() / F_C_MOL);
					temp_g.Set_dg(-s_x[i]->z);
					if (use.Get_surface_ptr()->Get_only_counter_ions() &&
						temp_g.Get_g() < 0)
					{
						temp_g.Set_g(0);
						temp_g.Set_dg(0);
					}
				}
				else
				{
					temp_g.Set_g(0);
					temp_g.Set_dg(-s_x[i]->z);
				}
				charge_ptr->Get_g_map()[s_x[i]->z] = temp_g;
			}

			{
				int is = s_x[i]->number;
				assert (is < (int) s_diff_layer.size());
				// species found in diff_layer
				s_diff_layer[is][charge_ptr->Get_name()].Set_g_moles(0);
				s_diff_layer[is][charge_ptr->Get_name()].Set_dg_g_moles(0);
			}
		}
		if (debug_diffuse_layer == TRUE)
		{
			output_msg(sformatf(
					   "\nSurface component %d: charge,\tg,\tdg\n",
					   (int) charge_ptr->Get_g_map().size()));
			std::map<LDBLE, cxxSurfDL>::iterator it;
			for (it = charge_ptr->Get_g_map().begin(); it != charge_ptr->Get_g_map().end(); it++)
			{
				output_msg(sformatf( "\t%12f\t%12.4e\t%12.4e\n",
						   (double) it->first,
						   (double) it->second.Get_g(),
						   (double) it->second.Get_dg()));
			}
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
initial_surface_water(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   In initial surface calculation, need to calculate
 *   mass of water in diffuse layer.
 *   diffuse layer water + aqueous solution water = bulk water.
 *   Ionic strength is fixed, so diffuse-layer water will not change
 */
	LDBLE debye_length, b, r, rd, ddl_limit, rd_limit, fraction, sum_surfs, l_s;
	LDBLE damp_aq;
/*
 *   Debye  length = 1/k = sqrt[eta*eta_zero*R*T/(2*F**2*mu_x*1000)], Dzombak and Morel, p 36
 *
 *   1000 converts kJ to J; 1000 converts Liters to meter**3; debye_length is in meters.
 */
	//debye_length = (EPSILON * EPSILON_ZERO * R_KJ_DEG_MOL * 1000.0 * tk_x)
	//	/ (2. * F_C_MOL * F_C_MOL * mu_x * 1000.);
	debye_length = (eps_r * EPSILON_ZERO * R_KJ_DEG_MOL * 1000.0 * tk_x)
		/ (2. * F_C_MOL * F_C_MOL * mu_x * 1000.);
	debye_length = sqrt(debye_length);

	/*  ddl is at most the fraction ddl_limit of bulk water */
	ddl_limit = use.Get_surface_ptr()->Get_DDL_limit();

/*
 *   Loop through all surface components, calculate each H2O surface (diffuse layer),
 *   H2O aq, and H2O bulk (diffuse layers plus aqueous).
 */

	if (use.Get_surface_ptr()->Get_debye_lengths() > 0)
	{
		sum_surfs = 0.0;
		for (int i = 0; i < count_unknowns; i++)
		{
			if (x[i]->type != SURFACE_CB)
				continue;
			cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[i]->surface_charge);
			sum_surfs +=
				charge_ptr->Get_specific_area() *
				charge_ptr->Get_grams();
		}

		rd = debye_length * use.Get_surface_ptr()->Get_debye_lengths();
		use.Get_surface_ptr()->Set_thickness(rd);

		if (state == INITIAL_SURFACE)
		{
			/* distribute water over DDL (rd) and free pore (r - rd) */
			/* find r: free pore (m3) = pi * (r - rd)^2 * L, where L = A / (2*pi*r),
			   A = sum_surfs = sum of the surface areas */
			b = -2 * (rd + use.Get_solution_ptr()->Get_mass_water() / (1000 * sum_surfs));
			r = 0.5 * (-b + sqrt(b * b - 4 * rd * rd));
			/* DDL (m3) = pi * (r^2 - (r - rd)^2) * L */
			rd_limit = (1 - sqrt(1 - ddl_limit)) * r;
			/* rd should be smaller than r and the ddl limit */
			if (rd > rd_limit)
			{
				mass_water_surfaces_x =
					use.Get_solution_ptr()->Get_mass_water() * ddl_limit / (1 -
																ddl_limit);
				r = 0.002 * (mass_water_surfaces_x +
							 use.Get_solution_ptr()->Get_mass_water()) / sum_surfs;
				rd_limit = (1 - sqrt(1 - ddl_limit)) * r;
				rd = rd_limit;
				use.Get_surface_ptr()->Set_thickness(rd);
			}
			else
				mass_water_surfaces_x =
					(r * r / pow(r - rd, 2) -
					 1) * use.Get_solution_ptr()->Get_mass_water();
			for (int i = 0; i < count_unknowns; i++)
			{
				if (x[i]->type != SURFACE_CB)
					continue;
				cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[i]->surface_charge);
				l_s =charge_ptr->Get_specific_area() *
					charge_ptr->Get_grams();
				charge_ptr->Set_mass_water(mass_water_surfaces_x * l_s / sum_surfs);
			}
		}
		else
		{
			r = 0.002 * mass_water_bulk_x / sum_surfs;
			rd_limit = (1 - sqrt(1 - ddl_limit)) * r;
			if (rd > rd_limit)
			{
				rd = rd_limit;
				use.Get_surface_ptr()->Set_thickness(rd);
				fraction = ddl_limit;
			}
			else
				fraction = 1 - pow(r - rd, 2) / (r * r);
			damp_aq = 1.0;
			if (g_iterations > 10)
				damp_aq = 0.2;
			else if (g_iterations > 5)
				damp_aq = 0.5;
			mass_water_surfaces_x = damp_aq * fraction * mass_water_bulk_x +
				(1 - damp_aq) * mass_water_surfaces_x;
			for (int i = 0; i < count_unknowns; i++)
			{
				if (x[i]->type != SURFACE_CB)
					continue;
				cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[i]->surface_charge);
				l_s = charge_ptr->Get_specific_area() *
					charge_ptr->Get_grams();
				charge_ptr->Set_mass_water(mass_water_surfaces_x * l_s / sum_surfs);
			}
		}
	}
	else
	{
		/* take constant thickness of, default 1e-8 m (100 Angstroms) */
		mass_water_surfaces_x = 0.0;
		for (int i = 0; i < count_unknowns; i++)
		{
			if (x[i]->type != SURFACE_CB)
				continue;
			cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[i]->surface_charge);
			charge_ptr->Set_mass_water(charge_ptr->Get_specific_area() *
				charge_ptr->Get_grams() * use.Get_surface_ptr()->Get_thickness() *
				1000);
			mass_water_surfaces_x += charge_ptr->Get_mass_water();
		}
	}

	if (use.Get_surface_ptr()->Get_type() == cxxSurface::CD_MUSIC)
		mass_water_bulk_x = mass_water_aq_x + mass_water_surfaces_x;
	else
	{
		/*  for variable distribution of water over DDL and bulk... */
		if (state > INITIAL_SURFACE)
			mass_water_aq_x = mass_water_bulk_x - mass_water_surfaces_x;
		else
			mass_water_bulk_x = mass_water_aq_x + mass_water_surfaces_x;
	}

	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
sum_diffuse_layer(cxxSurfaceCharge *charge_ptr)
/* ---------------------------------------------------------------------- */
{
	LDBLE mass_water_surface;
	LDBLE molality, moles_excess, moles_surface;

	if (use.Get_surface_ptr() == NULL)
		return (OK);
/*
 *   Find position of component in list of components
 */

/*
 *   Loop through all surface components, calculate each H2O surface (diffuse layer),
 *   H2O aq, and H2O bulk (diffuse layers plus aqueous).
 */
	count_elts = 0;
	paren_count = 0;
	mass_water_surface = charge_ptr->Get_mass_water();
	for (int j = 0; j < count_s_x; j++)
	{
		if (s_x[j]->type > HPLUS)
			continue;
		molality = under(s_x[j]->lm);
		LDBLE g = charge_ptr->Get_g_map()[s_x[j]->z].Get_g();
		if (s_x[j]->erm_ddl != 1)
		{
			LDBLE ratio_aq = mass_water_surface / mass_water_aq_x;
			LDBLE g2 = g / ratio_aq + 1;
			g = ratio_aq * (g2 * s_x[j]->erm_ddl - 1);
		}
		moles_excess = mass_water_aq_x * molality * g;
		moles_surface = mass_water_surface * molality + moles_excess;
/*
 *   Accumulate elements in diffuse layer
 */
		add_elt_list(s_x[j]->next_elt, moles_surface);
	}
	add_elt_list(s_h2o->next_elt, mass_water_surface / gfw_water);

	if (count_elts > 0)
	{
		qsort(elt_list, (size_t) count_elts,
			  (size_t) sizeof(struct elt_list), elt_list_compare);
		elt_list_combine();
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
calc_all_donnan(void)
/* ---------------------------------------------------------------------- */
{
	bool converge; 
	int cd_m;
	LDBLE new_g, f_psi, surf_chrg_eq, psi_avg, f_sinh, A_surf, ratio_aq;
	LDBLE new_g2, f_psi2, surf_chrg_eq2, psi_avg2, dif, var1;

	if (use.Get_surface_ptr() == NULL)
		return (OK);
	//f_sinh = sqrt(8000.0 * EPSILON * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000.0) *
	//		 tk_x * mu_x);
	f_sinh = sqrt(8000.0 * eps_r * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000.0) *
			 tk_x * mu_x);
/*
 *   calculate g for each surface...
 */
	converge = TRUE;
	for (int j = 0; j < count_unknowns; j++)
	{
		if (x[j]->type != SURFACE_CB)
			continue;
		cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[j]->surface_charge);

		if (debug_diffuse_layer == TRUE)
			output_msg(sformatf( "Calc_all_g, X[%d]\n", j));
/*
 *  sum eq of each charge number in solution...
 */
		std::map<LDBLE, LDBLE>::iterator it;
		for (it = charge_group_map.begin(); it != charge_group_map.end(); it++)
		{
			it->second = 0.0;
		}
		charge_group_map.clear();
		for (int i = 0; i < count_s_x; i++)
		{
			if (s_x[i]->type > HPLUS)
				continue;
			charge_group_map[s_x[i]->z] += s_x[i]->z * s_x[i]->moles * s_x[i]->erm_ddl;
		}
		/* find surface charge from potential... */
		A_surf = charge_ptr->Get_specific_area() * charge_ptr->Get_grams();
		if (use.Get_surface_ptr()->Get_type() == cxxSurface::CD_MUSIC)
		{
			f_psi = x[j + 2]->master[0]->s->la * LOG_10;	/* -FPsi/RT */
			f_psi = f_psi / 2;
			cd_m = 1;
		} else
		{
			f_psi = x[j]->master[0]->s->la * LOG_10;
			cd_m = -1;
		}
		surf_chrg_eq = A_surf * f_sinh * sinh(f_psi) / F_C_MOL;
		if (surf_chrg_eq < -5e3)
		{
			surf_chrg_eq = -5e3;
			var1 = surf_chrg_eq / (A_surf * f_sinh / F_C_MOL);
			var1 = (var1 + sqrt(var1 * var1 + 1));
			f_psi = (var1 > 1e-8 ? log(var1) : -18.4);
			surf_chrg_eq = A_surf * f_sinh * sinh(f_psi) / F_C_MOL;
			x[j]->master[0]->s->la = f_psi / LOG_10;
		}
		/* also for the derivative... */
		dif = 1e-5;
		f_psi2 = f_psi + dif;
		surf_chrg_eq2 = A_surf * f_sinh * sinh(f_psi2) / F_C_MOL;

		/* find psi_avg that matches surface charge... */
		psi_avg = calc_psi_avg(charge_ptr, surf_chrg_eq);
		psi_avg2 = calc_psi_avg(charge_ptr, surf_chrg_eq2);

		/*output_msg(sformatf( "psi's  %e %e %e\n", f_psi, psi_avg, surf_chrg_eq)); */

		/* fill in g's */
		ratio_aq = charge_ptr->Get_mass_water() / mass_water_aq_x;

		for (it = charge_group_map.begin(); it != charge_group_map.end(); it++)
		{
			LDBLE z = it->first;
			new_g = ratio_aq * (exp(cd_m * z * psi_avg) - 1);
			if (use.Get_surface_ptr()->Get_only_counter_ions() &&
				((surf_chrg_eq < 0 && z < 0)
				 || (surf_chrg_eq > 0 && z > 0)))
				new_g = -ratio_aq;
			if (new_g <= -ratio_aq)
				new_g = -ratio_aq + G_TOL * 1e-3;
			new_g2 = ratio_aq * (exp(cd_m * z * psi_avg2) - 1);
			if (use.Get_surface_ptr()->Get_only_counter_ions() &&
				((surf_chrg_eq < 0 && z < 0)
				 || (surf_chrg_eq > 0 && z > 0)))
				new_g2 = -ratio_aq;
			if (new_g2 <= -ratio_aq)
				new_g2 = -ratio_aq + G_TOL * 1e-3;
			if (fabs(new_g) >= 1)
			{
				if (fabs((new_g - charge_ptr->Get_g_map()[z].Get_g()) / new_g) > convergence_tolerance)
				{
					converge = FALSE;
				}
			}
			else
			{
				if (fabs(new_g - charge_ptr->Get_g_map()[z].Get_g()) > convergence_tolerance)
				{
					converge = FALSE;
				}
			}
			charge_ptr->Get_g_map()[z].Set_g(new_g);
			if (new_g != 0)
			{
				charge_ptr->Get_g_map()[z].Set_dg((new_g2 - new_g) / dif);
			}
			else
			{
				charge_ptr->Get_g_map()[z].Set_dg(-z);
			}
			/* save g for species */
		}
		if (debug_diffuse_layer == TRUE)
		{
			std::string name =  x[j]->master[0]->elt->name;
			Utilities::replace("_psi", "", name);
			output_msg(sformatf(
					   "\nDonnan all on %s (%d): charge, \tg, \tdg, Psi_surface = %8f V. \n",
					   name.c_str(), (int) charge_ptr->Get_g_map().size(),
					   x[j]->master[0]->s->la * 2 * LOG_10 * R_KJ_DEG_MOL *
					   tk_x / F_KJ_V_EQ));
			for (std::map<LDBLE, cxxSurfDL>::iterator i_it = charge_ptr->Get_g_map().begin();
				i_it != charge_ptr->Get_g_map().end(); i_it++)
			{
				output_msg(sformatf( "\t%12f\t%12.4e\t%12.4e\n",
						   (double) i_it->first,
						   (double) i_it->second.Get_g(),
						   (double) i_it->second.Get_dg()));
			}
		}
	}
	return (converge);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
calc_init_donnan(void)
/* ---------------------------------------------------------------------- */
{
	LDBLE f_psi, surf_chrg_eq, psi_avg, f_sinh, A_surf, ratio_aq;

	if (use.Get_surface_ptr() == NULL)
		return (OK);
	f_sinh =
		//sqrt(8000.0 * EPSILON * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000.0) *
		//	 tk_x * mu_x);
		sqrt(8000.0 * eps_r * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000.0) *
			 tk_x * mu_x);
	if (convergence_tolerance >= 1e-8)
	{
		G_TOL = 1e-9;
	}
	else
	{
		G_TOL = 1e-13;
	}
/*
 *  sum eq of each charge number in solution...
 */
	charge_group_map.clear();
	charge_group_map[0.0] = 0.0;

	for (int i = 0; i < count_s_x; i++)
	{
		if (s_x[i]->type > HPLUS)
			continue;
		if (charge_group_map.find(s_x[i]->z) != charge_group_map.end())
		{
			charge_group_map.find(s_x[i]->z)->second += s_x[i]->z * s_x[i]->moles * s_x[i]->erm_ddl;
		}
		else
		{
			charge_group_map[s_x[i]->z] = s_x[i]->z * s_x[i]->moles * s_x[i]->erm_ddl;
		}
	}
/*
 *   calculate g for each surface...
 */
	for (int j = 0; j < count_unknowns; j++)
	{
		if (x[j]->type != SURFACE_CB)
			continue;
		cxxSurfaceCharge *charge_ptr = use.Get_surface_ptr()->Find_charge(x[j]->surface_charge);
		charge_ptr->Get_g_map().clear();

		/* find surface charge from potential... */
		A_surf = charge_ptr->Get_specific_area() * charge_ptr->Get_grams();
		if (use.Get_surface_ptr()->Get_type() == cxxSurface::CD_MUSIC)
		{
			f_psi = x[j + 2]->master[0]->s->la * LOG_10;	/* -FPsi/RT */
			f_psi = f_psi / 2;
		} else
			f_psi = x[j]->master[0]->s->la * LOG_10;
		surf_chrg_eq = A_surf * f_sinh * sinh(f_psi) / F_C_MOL;

		/* find psi_avg that matches surface charge... */
/*    psi_avg = calc_psi_avg (0);
    appt 7/9/8... may have to change above one */
		psi_avg = calc_psi_avg(charge_ptr, 0 * surf_chrg_eq);

		/* fill in g's */
		ratio_aq = charge_ptr->Get_mass_water() / mass_water_aq_x;

		std::map<LDBLE, LDBLE>::iterator kit;
		for (kit = charge_group_map.begin(); kit != charge_group_map.end(); kit++)
		{
			LDBLE z = kit->first;
			LDBLE eq = kit->second;

			charge_ptr->Get_g_map()[z].Set_g(ratio_aq * (exp(-z * psi_avg) - 1));

			if (use.Get_surface_ptr()->Get_only_counter_ions()
				&& ((surf_chrg_eq < 0 && z < 0)
					|| (surf_chrg_eq > 0 && z > 0)))
				charge_ptr->Get_g_map()[z].Set_g(-ratio_aq);

			if (charge_ptr->Get_g_map()[z].Get_g() != 0)
			{
				charge_ptr->Get_g_map()[z].Set_dg(-A_surf * f_sinh * cosh(f_psi) / 
					(eq * F_C_MOL));
			}
			else
			{
				charge_ptr->Get_g_map()[z].Set_dg(-z);
			}

			/* save g for species */
			for (int i = 0; i < count_s_x; i++)
			{
				int is = s_x[i]->number;
				assert (is < (int) s_diff_layer.size());

				s_diff_layer[is][charge_ptr->Get_name()].Set_g_moles(0.0);
				s_diff_layer[is][charge_ptr->Get_name()].Set_dg_g_moles(0.0);
			}
		}
		if (debug_diffuse_layer == TRUE)
		{
			std::string name = x[j]->master[0]->elt->name;
			Utilities::replace("_psi", "", name);
			output_msg(sformatf(
					   "\nDonnan init on %s : charge, \tg, \tdg, Psi_surface = %8f V. \n",
					   name.c_str(),
					   x[j]->master[0]->s->la * 2 * LOG_10 * R_KJ_DEG_MOL *
					   tk_x / F_KJ_V_EQ));
			for (std::map<LDBLE, cxxSurfDL>::iterator i_it = charge_ptr->Get_g_map().begin();
				i_it != charge_ptr->Get_g_map().end(); i_it++)
			{
				output_msg(sformatf( "\t%12f\t%12.4e\t%12.4e\n",
						   (double) i_it->first,
						   (double) i_it->second.Get_g(),
						   (double) i_it->second.Get_dg()));
			}
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
calc_psi_avg(cxxSurfaceCharge *charge_ptr, LDBLE surf_chrg_eq)
/* ---------------------------------------------------------------------- */
{
/*
 * calculate the average (F * Psi / RT) that lets the DL charge counter the surface charge
 */
	LDBLE fd, fd1, p, temp, ratio_aq;

	ratio_aq = charge_ptr->Get_mass_water() / mass_water_aq_x;
	p = 0;
	if (surf_chrg_eq == 0 || ratio_aq == 0)
		return (0.0);
	else if (surf_chrg_eq < 0)
		p = -0.5 * log(-surf_chrg_eq * ratio_aq / mu_x + 1);
	else if (surf_chrg_eq > 0)
		p = 0.5 * log(surf_chrg_eq * ratio_aq / mu_x + 1);
/*
 * Optimize p in SS{s_x[i]->moles * z_i * g(p)} = -surf_chrg_eq
 *  g(p) = exp(-p * z_i) * ratio_aq
 * Elsewhere in PHREEQC, g is the excess, after subtraction of conc's for p = 0:
 *		      g(p) = (exp(-p *z_i) - 1) * ratio_aq
 */
	int l_iter = 0;
	do
	{
		fd = surf_chrg_eq;
		fd1 = 0.0;
		std::map<LDBLE, LDBLE>::iterator it;
		for (it = charge_group_map.begin(); it != charge_group_map.end(); it++)
		{
			LDBLE z = it->first;
			LDBLE eq = it->second;
			/*  multiply with ratio_aq for multiplier options cp and cm
				in calc_all_donnan (not used now)...  */
			temp = exp(-z * p) * ratio_aq;

			if (use.Get_surface_ptr()->Get_only_counter_ions() &&
				((surf_chrg_eq < 0 && z < 0)
				 || (surf_chrg_eq > 0 && z > 0)))
				temp = 0.0;
			fd += eq * temp;
			fd1 -= z * eq * temp;
		}
		fd /= -fd1;
		p += (fd > 1) ? 1 : ((fd < -1) ? -1 : fd);
		if (fabs(p) < G_TOL)
			p = 0.0;
		l_iter++;
		if (l_iter > 50)
		{
			error_string = sformatf(
					"\nToo many iterations in subroutine calc_psi_avg; surface charge = %12.4e; surface water = %12.4e.\n",
					(double) surf_chrg_eq, (double) charge_ptr->Get_mass_water());
			error_msg(error_string, STOP);
		}
	}
	while (fabs(fd) > 1e-12 && p != 0.0);
	if (debug_diffuse_layer == TRUE)
		output_msg(sformatf(
				   "iter in calc_psi_avg = %d. g(+1) = %8f. surface charge = %12.4e.\n",
				   l_iter, (double) (exp(-p) - 1), (double) surf_chrg_eq));

	return (p);
}
