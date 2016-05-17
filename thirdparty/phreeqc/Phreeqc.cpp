#include "Phreeqc.h"
#include <algorithm>			// std::replace

#include "NameDouble.h"
#include "Solution.h"
#include "Reaction.h"
#include "PPassemblage.h"
#include "Exchange.h"
#include "Surface.h"
#include "GasPhase.h"
#include "SSassemblage.h"
#include "cxxKinetics.h"
#include "phqalloc.h"
#include "PBasic.h"
#include "Temperature.h"
#include "SSassemblage.h"

const struct const_iso Phreeqc::iso_defaults[] = {
	{"13C", -10, 1},
	{"13C(4)", -10, 1},
	{"13C(-4)", -50, 5},
	{"34S", 10, 1},
	{"34S(6)", 10, 1},
	{"34S(-2)", -30, 5},
	{"2H", -28, 1},
	{"2H(1)", -28, 1},
	{"2H(0)", -28, 1},
	{"18O", -5, .1},
	{"18O(-2)", -5, .1},
	{"18O(0)", -5, .1},
	{"87Sr", .71, .01},
	{"11B", 20, 5}
};

const int Phreeqc::count_iso_defaults = (sizeof(iso_defaults) / sizeof(struct const_iso));

Phreeqc::~Phreeqc(void)
{

	clean_up();
	
	PHRQ_free_all();
	if (phrq_io == &ioInstance)
	{
		this->phrq_io->clear_istream();
		this->phrq_io->close_ostreams();
	}
}

void Phreeqc::set_phast(int tf)
{
	this->phast = tf;
}
size_t Phreeqc::list_components(std::list<std::string> &list_c)
/*
 *	 Find all elements in any class definition
 */
{
	cxxNameDouble accumulator;
	//accumulator.add("H", 1);
	//accumulator.add("O", 1);

	// solutions
	{
		std::map<int, cxxSolution>::const_iterator cit = Rxn_solution_map.begin();
		for (; cit !=  Rxn_solution_map.end(); cit++)
		{
			cxxSolution entity(cit->second);
			accumulator.add_extensive(entity.Get_totals(), 1.0);
		}
	}

	// irreversible reactions
	{
		std::map<int, cxxReaction>::const_iterator cit = Rxn_reaction_map.begin();
		for (; cit !=  Rxn_reaction_map.end(); cit++)
		{
			cxxReaction r_ptr(cit->second);
			reaction_calc(&r_ptr);
			accumulator.add_extensive(r_ptr.Get_elementList(), 1.0);
		}
	}

	// pure phases
	{
		std::map<int, cxxPPassemblage>::const_iterator cit = Rxn_pp_assemblage_map.begin();
		for (; cit !=  Rxn_pp_assemblage_map.end(); cit++)
		{
			cxxPPassemblage entity = cit->second;
			entity.totalize(this);
			accumulator.add_extensive(entity.Get_eltList(), 1.0);
		}
	}
	// exchangers
	{
		std::map<int, cxxExchange>::const_iterator cit = Rxn_exchange_map.begin();
		for (; cit !=  Rxn_exchange_map.end(); cit++)
		{
			cxxExchange entity = cit->second;
			entity.totalize();
			accumulator.add_extensive(entity.Get_totals(), 1.0);
		}
	}

	// surfaces
	{
		std::map<int, cxxSurface>::const_iterator cit = Rxn_surface_map.begin();
		for (; cit !=  Rxn_surface_map.end(); cit++)
		{
			cxxSurface entity = cit->second;
			entity.totalize();
			accumulator.add_extensive(entity.Get_totals(), 1.0);
		}
	}
	// gas phases
	{
		std::map<int, cxxGasPhase>::const_iterator cit = Rxn_gas_phase_map.begin();
		for (; cit !=  Rxn_gas_phase_map.end(); cit++)
		{
			cxxGasPhase entity = cit->second;
			entity.totalize(this);
			accumulator.add_extensive(entity.Get_totals(), 1.0);
		}
	}

	// solid-solutions
	{
		std::map<int, cxxSSassemblage>::const_iterator cit = Rxn_ss_assemblage_map.begin();
		for (; cit !=  Rxn_ss_assemblage_map.end(); cit++)
		{
			cxxSSassemblage entity = cit->second;
			entity.totalize(this);
			accumulator.add_extensive(entity.Get_totals(), 1.0);
		}
	}
	// kinetics
	{
		std::map<int, cxxKinetics>::iterator it = Rxn_kinetics_map.begin();
		for (; it !=  Rxn_kinetics_map.end(); it++)
		{
			calc_dummy_kinetic_reaction_tally(&(it->second));
			cxxKinetics entity = it->second;
			accumulator.add_extensive(entity.Get_totals(), 1.0);
		}
	}
	// Put in all primaries
	cxxNameDouble::iterator it;
	for (it = accumulator.begin(); it != accumulator.end(); it++)
	{
		if (it->first == "Charge") continue;
		char string[MAX_LENGTH];
		strcpy(string, it->first.c_str());
		struct master *master_ptr = master_bsearch_primary(string);
		if (master_ptr == NULL) continue;
		if (master_ptr->type != AQ) continue;
		accumulator.add(master_ptr->elt->name, 1);
	}
	// print list
	for (it = accumulator.begin(); it != accumulator.end(); it++)
	{
		struct master *master_ptr = master_bsearch(it->first.c_str());
		if (master_ptr == NULL) continue;
		if (master_ptr->type != AQ) continue;
		if (master_ptr->primary == 0) continue;
		if (it->first == "Charge") continue;
		if (it->first == "O") continue;
		if (it->first == "H") continue;
		list_c.push_back(it->first);
	}
	return(list_c.size());
}
Phreeqc::Phreeqc(PHRQ_io *io)
{
	// phrq_io
	if (io)
	{
		this->phrq_io = io;
	}
	else
	{
		this->phrq_io = &this->ioInstance;
	}
	// auto PHRQ_io ioInstance;

	// initialize data members
	init();

#if defined(SWIG) || defined(SWIG_IPHREEQC)
	basicCallback = NULL;
#endif
}
void Phreeqc::init(void)
{
	same_model                      = FALSE;
	current_tc                      = NAN;
	current_pa                      = NAN;
	current_mu                      = NAN;
	mu_terms_in_logk                = true;

	/* ----------------------------------------------------------------------
	*   STRUCTURES
	* ---------------------------------------------------------------------- */
/*
 *	 last model
 */
	last_model.force_prep           = TRUE;
	last_model.temperature          = -100;
	last_model.pressure             = 0;
	last_model.count_exchange       = -1;
	last_model.exchange             = NULL;
	last_model.count_kinetics       = -1;
	last_model.kinetics             = NULL;
	last_model.count_gas_phase      = -1;
	last_model.gas_phase            = NULL;
	last_model.count_ss_assemblage  = -1;
	last_model.ss_assemblage        = NULL;
	last_model.count_pp_assemblage  = -1;
	last_model.pp_assemblage        = NULL;
	last_model.add_formula          = NULL;
	last_model.si                   = NULL;
	last_model.dl_type              = cxxSurface::NO_DL;
	last_model.surface_type         = cxxSurface::UNKNOWN_DL;
	last_model.only_counter_ions    = FALSE;
	last_model.thickness            = 1e-8;
	last_model.count_surface_comp   = -1;
	last_model.surface_comp         = NULL;
	last_model.count_surface_charge = -1;
	last_model.surface_charge       = NULL;

	current_selected_output         = NULL;
	current_user_punch              = NULL;
	high_precision                  = false;
#ifdef SKIP
	//struct punch punch;
/*
 *	 Initialize punch
 */
	punch.in				= FALSE;
	punch.count_totals		= 0;
	punch.totals			= 0;
	punch.count_molalities	= 0;
	punch.molalities		= 0;
	punch.count_activities	= 0;
	punch.activities		= 0;
	punch.count_pure_phases = 0;
	punch.pure_phases		= 0;
	punch.count_si			= 0;
	punch.si				= 0;
	punch.count_gases		= 0;
	punch.gases				= 0;
	punch.count_s_s			= 0;
	punch.s_s               = 0;
	punch.count_kinetics	= 0;
	punch.kinetics			= 0;
	punch.count_isotopes	= 0;
	punch.isotopes			= 0;
	punch.count_calculate_values = 0;
	punch.calculate_values	= 0;
	punch.inverse			= TRUE;
	punch.sim				= TRUE;
	punch.state				= TRUE;
	punch.soln				= TRUE;
	punch.dist				= TRUE;
	punch.time				= TRUE;
	punch.step				= TRUE;
	punch.rxn				= FALSE;
	punch.temp				= FALSE;
	punch.ph				= TRUE;
	punch.pe				= TRUE;
	punch.alk				= FALSE;
	punch.mu				= FALSE;
	punch.water				= FALSE;
	punch.high_precision	= FALSE;
	punch.user_punch		= TRUE;
	punch.charge_balance	= FALSE;
	punch.percent_error		= FALSE;
#endif

	MIN_LM = -30.0;			    /* minimum log molality allowed before molality set to zero */
	LOG_ZERO_MOLALITY = -30;	/* molalities <= LOG_ZERO_MOLALITY are considered equal to zero */
	MIN_RELATED_LOG_ACTIVITY = -30;
	MIN_TOTAL = 1e-25;
	MIN_TOTAL_SS = MIN_TOTAL/100;
	MIN_RELATED_SURFACE = MIN_TOTAL*100;
	// auto Rxn_temperature_map;
	// auto Rxn_pressure_map;

	/* ----------------------------------------------------------------------
	*   Surface
	* --------------------------------------------------------------------- */
	g_iterations               = -1;
	G_TOL                      = 1e-8;
	// auto Rxn_surface_map;
	// auto charge_group_map;
	change_surf_count          = 0;
	change_surf                = NULL;
	/* ----------------------------------------------------------------------
	*   Exchange
	* ---------------------------------------------------------------------- */
	// auto Rxn_exchange_map;

	/* ----------------------------------------------------------------------
	*   Kinetics
	* ---------------------------------------------------------------------- */
	// auto Rxn_kinetics_map;

	/*----------------------------------------------------------------------
	*   Save
	*---------------------------------------------------------------------- */
	count_save_values          = 0;
	save_values                = NULL;
	save_init(-1);             // set initial save values

	// auto use

	// copier structures
	copy_solution.n_user		= copy_solution.start		= copy_solution.end			= 0;
	copy_solution.count		    = copy_solution.max                                     = 0;
	copy_pp_assemblage.n_user	= copy_pp_assemblage.start	= copy_pp_assemblage.end	= 0;
	copy_pp_assemblage.count	= copy_pp_assemblage.max                                = 0;
	copy_exchange.n_user		= copy_exchange.start		= copy_exchange.end			= 0;
	copy_exchange.count		    = copy_exchange.max                                     = 0;
	copy_surface.n_user			= copy_surface.start		= copy_surface.end			= 0;
	copy_surface.count		    = copy_surface.max                                      = 0;
	copy_ss_assemblage.n_user	= copy_ss_assemblage.start = copy_ss_assemblage.end		= 0;
	copy_ss_assemblage.count	= copy_ss_assemblage.max                                = 0;
	copy_gas_phase.n_user		= copy_gas_phase.start		= copy_gas_phase.end		= 0;
	copy_gas_phase.count		= copy_gas_phase.max                                    = 0;
	copy_kinetics.n_user		= copy_kinetics.start		= copy_kinetics.end			= 0;
	copy_kinetics.count		    = copy_kinetics.max                                     = 0;
	copy_mix.n_user				= copy_mix.start			= copy_mix.end				= 0;
	copy_mix.count		        = copy_mix.max                                          = 0;
	copy_reaction.n_user		= copy_reaction.start		= copy_reaction.end			= 0;
	copy_reaction.count		    = copy_reaction.max                                     = 0;
	copy_temperature.n_user		= copy_temperature.start	= copy_temperature.end		= 0;
	copy_temperature.count		= copy_temperature.max                                  = 0;
	copy_pressure.n_user		= copy_pressure.start		= copy_pressure.end			= 0;
	copy_pressure.count		    = copy_pressure.max                                     = 0;
	/*----------------------------------------------------------------------
	*   Inverse
	*---------------------------------------------------------------------- */
	inverse					= NULL;
	count_inverse			= 0;
	/*----------------------------------------------------------------------
	*   Mix
	*---------------------------------------------------------------------- */
	// auto Rxn_mix_map;
	// auto Dispersion_mix_map;
	// auto Rxn_solution_mix_map;
	// auto Rxn_exchange_mix_map;
	// auto Rxn_gas_phase_mix_map;
	// auto Rxn_kinetics_mix_map;
	// auto Rxn_pp_assemblage_mix_map;
	// auto Rxn_ss_assemblage_mix_map;
	// auto Rxn_surface_mix_map;
	/*----------------------------------------------------------------------
	*   Irreversible reaction
	*---------------------------------------------------------------------- */
	run_cells_one_step = false;
	// auto Rxn_reaction_map;
	/*----------------------------------------------------------------------
	*   Gas phase
	*---------------------------------------------------------------------- */
	// auto Rxn_gas_phase_map;
	/*----------------------------------------------------------------------
	*   Solid solution
	*---------------------------------------------------------------------- */
	// auto Rxn_ss_assemblage_map;
	/*----------------------------------------------------------------------
	*   Pure-phase assemblage
	*---------------------------------------------------------------------- */
	// auto Rxn_pp_assemblage_map;
	/*----------------------------------------------------------------------
	*   Species_list
	*---------------------------------------------------------------------- */
	count_species_list      = 0;
	max_species_list        = 0;
	species_list            = NULL;
	/*----------------------------------------------------------------------
	*   Jacobian and Mass balance lists
	*---------------------------------------------------------------------- */
	count_sum_jacob0        = 0;
	max_sum_jacob0          = 0;
	sum_jacob0              = NULL;	
	count_sum_mb1           = 0;
	max_sum_mb1             = 0;
	sum_mb1                 = NULL;	
	count_sum_jacob1        = 0;
	max_sum_jacob1          = 0;
	sum_jacob1              = NULL;
	count_sum_mb2           = 0;
	max_sum_mb2             = 0;
	sum_mb2                 = NULL;
	count_sum_jacob2        = 0;
	max_sum_jacob2          = 0;
	sum_jacob2              = NULL;
	count_sum_delta         = 0;
	max_sum_delta           = 0;
	sum_delta               = NULL;
	/*----------------------------------------------------------------------
	*   Solution
	*---------------------------------------------------------------------- */
	// auto Rxn_solution_map;
	// auto unnumbered_solutions;
	save_species = false;
	/*----------------------------------------------------------------------
	*   Global solution
	*---------------------------------------------------------------------- */
	title_x                 = NULL;
	new_x                   = FALSE;
	description_x			= NULL;
	tc_x                    = 0;
	tk_x                    = 0;
	patm_x                  = 1;
	last_patm_x             = 1;
	numerical_fixed_volume  = false;
	force_numerical_fixed_volume = false;
	//switch_numerical        = false;
	ph_x                    = 0;
	solution_pe_x           = 0;
	mu_x                    = 0;
	ah2o_x                  = 1.0;
	density_x               = 0;
	total_h_x               = 0;
	total_o_x               = 0;
	cb_x                    = 0;
	total_ions_x            = 0;
	mass_water_aq_x         = 0;
	mass_water_surfaces_x   = 0;
	mass_water_bulk_x       = 0;
	units_x					= NULL;
	// auto pe_x
	// auto isotopes_x
	// auto default_pe_x
	dl_type_x                = cxxSurface::NO_DL;
	total_carbon             = 0;
	total_co2                = 0;
	total_alkalinity         = 0;
	gfw_water                = 0;
	step_x                   = 0;
	kin_time_x               = 0;
	/*----------------------------------------------------------------------
	*   Transport data
	*---------------------------------------------------------------------- */
	count_cells              = 1;
	count_shifts             = 1;
	ishift                   = 1;
	bcon_first = bcon_last   = 3;
	correct_disp             = FALSE;
	tempr                    = 2.0;
	timest                   = 0.0;
	simul_tr                 = 0;
	diffc                    = 0.3e-9;
	heat_diffc               = -0.1;
	cell                     = 0;
	mcd_substeps             = 1.0;
	stag_data                = NULL;
	print_modulus            = 1;
	punch_modulus            = 1;
	dump_in                  = FALSE;
	dump_modulus             = 0;
	transport_warnings       = TRUE;
	cell_data                = NULL;
	old_cells                = 0;
	max_cells                = 0;
	all_cells                = 0;
	multi_Dflag              = FALSE;
	interlayer_Dflag         = FALSE;
	default_Dw               = 0;
	multi_Dpor               = 0;
	interlayer_Dpor          = 0.1;
	multi_Dpor_lim           = 0;
	interlayer_Dpor_lim      = 0;
	multi_Dn                 = 0;
	interlayer_tortf         = 100.0;
	cell_no                  = 0;
	/*----------------------------------------------------------------------
	*   Advection data
	*---------------------------------------------------------------------- */
	count_ad_cells           = 1;
	count_ad_shifts          = 1;
	print_ad_modulus         = 1;
	punch_ad_modulus         = 1;
	advection_punch          = NULL;
	advection_kin_time       = 0.0;
	advection_kin_time_defined = FALSE;
	advection_print          = NULL;
	advection_warnings       = TRUE;
	/*----------------------------------------------------------------------
	*   Tidy data
	*---------------------------------------------------------------------- */
	new_model                = TRUE;
	new_exchange             = FALSE;
	new_pp_assemblage        = FALSE;
	new_surface              = FALSE;
	new_reaction             = FALSE;
	new_temperature          = FALSE;
	new_mix                  = FALSE;
	new_solution             = FALSE;
	new_gas_phase            = FALSE;
	new_inverse              = FALSE;
	new_punch                = FALSE;
	new_ss_assemblage        = FALSE;
	new_kinetics             = FALSE;
	new_copy                 = FALSE;
	new_pitzer               = FALSE;
	/*----------------------------------------------------------------------
	*   Elements
	*---------------------------------------------------------------------- */
	elements                 = NULL;
	count_elements           = 0;
	max_elements             = MAX_ELEMENTS;
	element_h_one            = NULL;
	/*----------------------------------------------------------------------
	*   Element List
	*---------------------------------------------------------------------- */
	elt_list                 = NULL;
	count_elts               = 0;
	max_elts                 = MAX_ELTS;
	/*----------------------------------------------------------------------
	*   Species
	*---------------------------------------------------------------------- */
	logk                     = NULL;
	count_logk               = 0;
	max_logk                 = MAX_S;
	moles_per_kilogram_string= NULL;
	pe_string                = NULL;
	s                        = NULL;
	count_s                  = 0;
	max_s                    = MAX_S;
	// auto s_diff_layer;
	s_x                      = NULL;
	count_s_x                = 0;
	max_s_x                  = 0;
	s_h2o					= NULL;
	s_hplus					= NULL;
	s_h3oplus				= NULL;
	s_eminus				= NULL;
	s_co3					= NULL;
	s_h2					= NULL;
	s_o2					= NULL;
	/*----------------------------------------------------------------------
	*   Phases
	*---------------------------------------------------------------------- */
	phases					= NULL;
	count_phases            = 0;
	max_phases              = MAX_PHASES;
	/*----------------------------------------------------------------------
	*   Master species
	*---------------------------------------------------------------------- */
	master                  = NULL;
	dbg_master              = NULL;
	count_master            = 0;
	max_master              = MAX_MASTER;
	/*----------------------------------------------------------------------
	*   Unknowns
	*---------------------------------------------------------------------- */
	x                       = NULL;
	count_unknowns          = 0;
	max_unknowns            = 0;
	ah2o_unknown            = NULL;
	alkalinity_unknown      = NULL;
	carbon_unknown          = NULL;
	charge_balance_unknown  = NULL;
	exchange_unknown        = NULL;
	mass_hydrogen_unknown   = NULL;
	mass_oxygen_unknown     = NULL;
	mb_unknown              = NULL;
	mu_unknown              = NULL;
	pe_unknown              = NULL;
	ph_unknown              = NULL;
	pure_phase_unknown      = NULL;
	solution_phase_boundary_unknown = NULL;
	surface_unknown         = NULL;
	gas_unknown             = NULL;
	ss_unknown              = NULL;
	// auto gas_unknowns;
	/*----------------------------------------------------------------------
	*   Reaction work space
	*---------------------------------------------------------------------- */
	// struct trxn;	
	trxn.token				= 0;
	for (int i = 0; i < MAX_LOG_K_INDICES; i++)
	{
		trxn.logk[i] = 0;
	}
	for (int i = 0; i < 3; i++)
	{
		trxn.dz[i] = 0;
	}
	count_trxn              = 0;
	max_trxn                = MAX_TRXN;

	mb_unknowns             = NULL;
	count_mb_unknowns       = 0;
	max_mb_unknowns         = MAX_TRXN;
	/* ----------------------------------------------------------------------
	*   Print
	* ---------------------------------------------------------------------- */
	pr.all                  = TRUE;
	pr.initial_solutions    = TRUE;
	pr.initial_exchangers   = TRUE;
	pr.reactions            = TRUE;
	pr.gas_phase            = TRUE;
	pr.ss_assemblage        = TRUE;
	pr.pp_assemblage        = TRUE;
	pr.surface              = TRUE;
	pr.exchange             = TRUE;
	pr.kinetics             = TRUE;
	pr.totals               = TRUE;
	pr.eh                   = TRUE;
	pr.species              = TRUE;
	pr.saturation_indices   = TRUE;
	pr.irrev                = TRUE;
	pr.mix                  = TRUE;
	pr.reaction             = TRUE;
	pr.use                  = TRUE;
	pr.logfile              = FALSE;
	pr.punch                = TRUE;
	pr.status               = TRUE;
	pr.inverse              = TRUE;
	pr.dump                 = TRUE;
	pr.user_print           = TRUE;
	pr.headings             = TRUE;
	pr.user_graph           = TRUE;
	pr.echo_input           = TRUE;
	pr.warnings             = 100;
	pr.initial_isotopes     = TRUE;
	pr.isotope_ratios       = TRUE;
	pr.isotope_alphas       = TRUE;
	pr.hdf                  = FALSE;
	pr.alkalinity           = FALSE;
	status_on               = true;
#ifdef NPP
	status_interval         = 40;
#else
	status_interval         = 250;
#endif
	status_timer            = clock();
	count_warnings          = 0;
	/* ----------------------------------------------------------------------
	*   RATES
	* ---------------------------------------------------------------------- */
	rates                   = NULL;
	count_rates				= 0;
	rate_m					= 0;
	rate_m0					= 0;
	rate_time				= 0;
	rate_kin_time           = 1.0;
	rate_sim_time_start		= 0;
	rate_sim_time_end		= 0;
	rate_sim_time			= 0;
	rate_moles				= 0;
	initial_total_time		= 0;
	// auto rate_p
	count_rate_p            = 0;
	/* ----------------------------------------------------------------------
	*   USER PRINT COMMANDS
	* ---------------------------------------------------------------------- */
	user_print				= NULL;
#ifdef SKIP
	user_punch				= NULL;
	user_punch_headings		= NULL;
	user_punch_count_headings = 0;
#endif
	n_user_punch_index      = 0;
	fpunchf_user_s_warning  = 0;
	fpunchf_user_buffer[0]  = 0;

#if defined PHREEQ98 
	struct rate *user_graph;
	char **user_graph_headings;
	int user_graph_count_headings;
#endif
#if defined MULTICHART
	// auto chart_handler;
	chart_handler.Set_io(phrq_io);
#endif
	/* ----------------------------------------------------------------------
	*   GLOBAL DECLARATIONS
	* ---------------------------------------------------------------------- */
	error_string            = NULL;
	simulation				= 0;
	state                   = INITIALIZE;
	reaction_step           = 0;
	transport_step          = 0;
	transport_start         = 0;
	advection_step          = 0;
	stop_program            = FALSE;
	incremental_reactions   = FALSE;
	count_strings           = 0;
	array					= NULL;
	delta					= NULL;
	residual				= NULL;
	input_error             = 0;
	next_keyword            = Keywords::KEY_NONE;
	parse_error             = 0;
	paren_count             = 0;
	iterations              = 0;
	gamma_iterations        = 0;
	run_reactions_iterations= 0;
	max_line				= MAX_LINE;
	line                    = NULL;
	line_save				= NULL;
	LOG_10                  = log(10.0);
	debug_model             = FALSE;
	debug_prep              = FALSE;
	debug_set               = FALSE;
	debug_diffuse_layer     = FALSE;
	debug_inverse           = FALSE;
#ifdef USE_LONG_DOUBLE
	/* from float.h, sets tolerance for cl1 routine */
	inv_tol_default         = pow((long double) 10, (long double) -LDBL_DIG + 5);
#else
	inv_tol_default         = pow((double) 10, (double) -DBL_DIG + 5);
#endif
	itmax                   = 100;
	max_tries               = 1000;
#ifdef USE_LONG_DOUBLE
	/* from float.h, sets tolerance for cl1 routine */
	ineq_tol                = pow((long double) 10, (long double) -LDBL_DIG);
#else
	ineq_tol                = pow((double) 10, (double) -DBL_DIG);
#endif
	convergence_tolerance   = 1e-8;	
	step_size				= 100.;
	pe_step_size			= 10.;
	step_size_now           = step_size;
	pe_step_size_now        = pe_step_size;
	pp_scale				= 1.0;
	pp_column_scale			= 1.0;
	diagonal_scale			= FALSE;
	mass_water_switch		= FALSE;
	delay_mass_water		= FALSE;
	equi_delay      		= 0;
	dampen_ah2o             = false;
	censor					= 0.0;
	aqueous_only			= 0;
	negative_concentrations = FALSE;
	calculating_deriv		= FALSE;
	numerical_deriv			= FALSE;
	count_total_steps       = 0;
	phast                   = FALSE;
	llnl_temp				= 0;
	llnl_count_temp			= 0;
	llnl_adh				= 0;
	llnl_count_adh			= 0;
	llnl_bdh				= 0;
	llnl_count_bdh			= 0;
	llnl_bdot				= 0;
	llnl_count_bdot			= 0;
	llnl_co2_coefs			= 0;
	llnl_count_co2_coefs	= 0;
	//selected_output_file_name = NULL;
	dump_file_name			= NULL;
	remove_unstable_phases  = FALSE;
	// auto screen_string;
	spread_length           = 10;
	/* ---------------------------------------------------------------------- */
	/*
	*   Hash definitions
	*/
	// auto strings_map;
#ifdef HASH
	// auto strings_hash;
#endif
	elements_hash_table     = NULL;
	species_hash_table      = NULL;
	phases_hash_table       = NULL;
	logk_hash_table         = NULL;
	master_isotope_hash_table = NULL;
	/* ----------------------------------------------------------------------
	*   ISOTOPES
	* ---------------------------------------------------------------------- */
	count_master_isotope	= 0;
	master_isotope			= NULL;
	max_master_isotope		= MAX_ELTS;
	initial_solution_isotopes = FALSE;
	count_calculate_value	= 0;
	calculate_value			= NULL;
	max_calculate_value		= MAX_ELTS;
	calculate_value_hash_table = NULL;	
	count_isotope_ratio		= 0;
	isotope_ratio			= 0;
	max_isotope_ratio		= MAX_ELTS;
	isotope_ratio_hash_table = 0;	
	count_isotope_alpha		= 0;
	isotope_alpha			= 0;
	max_isotope_alpha		= MAX_ELTS;
	isotope_alpha_hash_table = 0;


	phreeqc_mpi_myself		= 0;
	first_read_input		= TRUE;
	user_database			= NULL;
	//have_punch_name			= FALSE;
	print_density		    = 0;
	zeros                   = NULL;	
	zeros_max			    = 1;
	cell_pore_volume	    = 0;
	cell_volume			    = 0;
	cell_porosity		    = 0;
	cell_saturation		    = 0;
	sys                     = NULL;
	count_sys               = 0;
	max_sys                 = 0;
	sys_tot                 = 0;

	V_solutes               = 0.0;
	rho_0                   = 0;
	kappa_0                 = 0.0;
	p_sat                   = 0.0;
	eps_r                   = EPSILON;
	DH_A                    = 0.0;
	DH_B                    = 0.0;
	DH_Av                   = 0.0;
	QBrn                    = 0.0;
	ZBrn                    = 0.0;
	dgdP                    = 0.0;

	need_temp_msg           = 0;
	solution_mass           = 0;
	solution_volume         = 0;
	/* phqalloc.cpp ------------------------------- */
	s_pTail                 = NULL;
	/* Basic */
	basic_interpreter       = NULL;
	basic_callback_ptr      = NULL;
	basic_callback_cookie   = NULL;
	basic_fortran_callback_ptr  = NULL;

	/* cl1.cpp ------------------------------- */
	x_arg                   = NULL; 
	res_arg                 = NULL; 
	scratch                 = NULL;
	x_arg_max               = 0; 
	res_arg_max             = 0; 
	scratch_max             = 0;
#ifdef SKIP
	/* dw.cpp ------------------------------- */
	/* COMMON /QQQQ/ */	
	Q0                      = 0;
	Q5                      = 0;
	GASCON                  = 0.461522e0;
	TZ                      = 647.073e0;
	AA                      = 1.e0;
	Z                       = 0;
	DZ                      = 0;
	Y                       = 0;
	G1                      = 11.e0;
	G2                      = 44.333333333333e0;
	GF                      = 3.5e0;
	B1                      = 0;
	B2                      = 0;
	B1T                     = 0;
	B2T                     = 0;
	B1TT                    = 0;
	B2TT                    = 0;
#endif
	/* gases.cpp ------------------------------- */
	a_aa_sum                = 0;
	b2                      = 0;
	b_sum                   = 0;
	R_TK                    = 0;
	/* input.cpp ------------------------------- */
	check_line_return       = 0;  
	reading_db              = FALSE;
	/* integrate.cpp ------------------------------- */
	midpoint_sv             = 0;
	z_global                = 0;
	xd_global               = 0;
	alpha_global            = 0;
	/* integrate.cpp ------------------------------- */
	max_row_count           = 50;
	max_column_count        = 50;
	carbon                  = FALSE;
	col_name                = NULL;
	row_name                = NULL;
	count_rows              = 0;
	count_optimize          = 0;
	col_phases              = 0;
	col_redox               = 0;
	col_epsilon             = 0;
	col_ph                  = 0;
	col_water               = 0;
	col_isotopes            = 0;
	col_phase_isotopes      = 0;
	row_mb                  = 0;
	row_fract               = 0;
	row_charge              = 0;
	row_carbon              = 0;
	row_isotopes            = 0;
	row_epsilon             = 0;
	row_isotope_epsilon     = 0;
	row_water               = 0;
	inv_zero                = NULL;
	array1                  = 0;
	inv_res                 = NULL;
	inv_delta1              = NULL;
	delta2                  = NULL;
	delta3                  = NULL;
	inv_cu                  = NULL;
	delta_save              = NULL;
	min_delta               = NULL;
	max_delta               = NULL;
	inv_iu                  = NULL;
	inv_is                  = NULL;
	klmd                    = 0;
	nklmd                   = 0;
	n2d                     = 0;
	kode                    = 0;
	iter                    = 0;
	toler                   = 0;
	error                   = 0;
	max_pct                 = 0;
	scaled_error            = 0;
	master_alk              = NULL;
	row_back                = NULL;
	col_back                = NULL;
	good                    = NULL;
	bad                     = NULL;
	minimal                 = NULL;
	max_good                = 0;
	max_bad                 = 0;
	max_minimal             = 0;
	count_good              = 0;
	count_bad               = 0;
	count_minimal           = 0;
	count_calls             = 0;
	soln_bits               = 0;
	phase_bits              = 0;
	current_bits            = 0;
	temp_bits               = 0;
	netpath_file            = NULL;
	count_inverse_models    = 0;
	count_pat_solutions     = 0;
	for (int i = 0; i < 32; i++)
	{
		min_position[i]     = 0;
		max_position[i]     = 0;
		now[i]              = 0;
	}
	/* kinetics.cpp ------------------------------- */
	count_pp = count_pg = count_ss = 0; 
	cvode_kinetics_ptr      = NULL;
	cvode_test              = FALSE;
	cvode_error             = FALSE;
	cvode_n_user            = -99;
	cvode_n_reactions       = -99;
	cvode_step_fraction     = 0.0;
	cvode_rate_sim_time     = 0.0;
	cvode_rate_sim_time_start = 0.0;
	cvode_last_good_time    = 0.0;
	cvode_prev_good_time    = 0.0;
	cvode_last_good_y       = NULL;
	cvode_prev_good_y       = NULL;
	kinetics_machEnv        = NULL;
	kinetics_y              = NULL;
	kinetics_abstol         = NULL;
	kinetics_cvode_mem      = NULL;
	cvode_pp_assemblage_save= NULL;
	cvode_ss_assemblage_save= NULL;
	m_original              = NULL;
	m_temp                  = NULL;
	rk_moles                = NULL;
	set_and_run_attempt     = 0;
	x0_moles                = NULL;
	/* model.cpp ------------------------------- */
	gas_in                  = FALSE;
	min_value               = 1e-10;
	normal                  = NULL;
	ineq_array              = NULL;
	res                     = NULL;
	cu                      = NULL;
	zero                    = NULL;
	delta1                  = NULL;
	iu                      = NULL;
	is                      = NULL;
	back_eq                 = NULL;
	normal_max              = 0;
	ineq_array_max          = 0;
	res_max                 = 0;
	cu_max                  = 0;
	zero_max                = 0;
	delta1_max              = 0;
	iu_max                  = 0;
	is_max                  = 0;
	back_eq_max             = 0;
	/* phrq_io_output.cpp ------------------------------- */
	forward_output_to_log   = 0;
	/* phreeqc_files.cpp ------------------------------- */
		default_data_base       = string_duplicate("phreeqc.dat");
#ifdef PHREEQ98
	int outputlinenr;
	char *LogFileNameC;
	char progress_str[512];
#endif
	/* Pitzer  */	
	pitzer_model			= FALSE;
	sit_model				= FALSE;
	pitzer_pe				= FALSE;
	full_pitzer = FALSE;
	always_full_pitzer = FALSE;
	ICON					= TRUE;
	IC                      = -1;
	COSMOT                  = 0;
	AW                      = 0;
	VP                      = 0;
	DW0                     = 0;
	pitz_params				= NULL;
	count_pitz_param		= 0;
	max_pitz_param			= 100;
	// auto pitz_param_map
	theta_params			= 0;
	count_theta_param		= 0;
	max_theta_param			= 100;
	use_etheta				= TRUE;
	OTEMP					= -100.;
	OPRESS					= -100.;
	A0                      = 0;	
	spec                    = NULL;
	cations                 = NULL;
	anions                  = NULL;
	neutrals                = NULL;
	count_cations           = 0;
	count_anions            = 0;
	count_neutrals          = 0;
	MAXCATIONS              = 0;
	FIRSTANION              = 0;
	MAXNEUTRAL              = 0;
	mcb0                    = NULL;
	mcb1                    = NULL;
	mcc0                    = NULL;
	IPRSNT                  = NULL;
	M                       = NULL;
	LGAMMA                  = NULL;
	for (int i = 0; i < 23; i++)
	{
		BK[i]				= 0.0;
		DK[i]				= 0.0;
	}
#ifdef PHREEQ98
	int connect_simulations, graph_initial_solutions;
	int shifts_as_points;
	int chart_type;
	int ShowChart;
	int RowOffset, ColumnOffset;
#endif
	dummy                   = 0;
	/* print.cpp ------------------------------- */
	sformatf_buffer = (char *) PHRQ_malloc(256 * sizeof(char));
	if (sformatf_buffer == NULL) 
			malloc_error();
		sformatf_buffer_size = 256;
#ifdef PHREEQ98
	int colnr, rownr;
	int graph_initial_solutions;
	int prev_advection_step, prev_transport_step;	/*, prev_reaction_step */
	/* int shifts_as_points; */
	int chart_type;
	int AddSeries;
	int FirstCallToUSER_GRAPH;
#endif
	/* read.cpp */
	prev_next_char          = NULL;
#if defined PHREEQ98 
	int shifts_as_points;
#endif
	/* read_class.cxx */
	// auto dump_info
	// auto delete_info
	// auto run_info
	run_info.Set_io(phrq_io);
	/* readtr.cpp */
	// auto dump_file_name_cpp;
	/* sit.cpp ------------------------------- */
	sit_params              = NULL;
	count_sit_param			= 0;
	max_sit_param			= 100;
	// auto sit_param_map	
	sit_A0                  = 0;
	sit_count_cations       = 0;
	sit_count_anions        = 0;
	sit_count_neutrals      = 0;
	sit_MAXCATIONS          = 0;
	sit_FIRSTANION          = 0;
	sit_MAXNEUTRAL          = 0;
	sit_IPRSNT              = NULL;
	sit_M                   = NULL;
	sit_LGAMMA              = NULL;
	/* tidy.cpp ------------------------------- */
	a0                      = 0;
	a1                      = 0;
	kc                      = 0;
	kb                      = 0;
	/* tally.cpp ------------------------------- */
	t_buffer                = NULL;
	tally_count_component   = 0;
	tally_table             = NULL;
	count_tally_table_columns = 0;
	count_tally_table_rows  = 0;
	/* transport.cpp ------------------------------- */
	sol_D                   = NULL;
	sol_D_dbg               = NULL;
	J_ij                    = NULL;
	J_ij_il                 = NULL;
	J_ij_count_spec         = 0;
	m_s                     = NULL;
	count_m_s               = 0;
	tot1_h                  = 0;
	tot1_o                  = 0;
	tot2_h                  = 0;
	tot2_o                  = 0;
	diffc_max               = 0;
	diffc_tr                = 0;
	J_ij_sum                = 0;
	transp_surf             = FALSE;
	heat_mix_array          = NULL;
	temp1                   = NULL;
	temp2                   = NULL;
	nmix                    = 0;
	heat_nmix               = 0;
	heat_mix_f_imm          = 0;
	heat_mix_f_m            = 0;
	warn_MCD_X              = 0;
	warn_fixed_Surf         = 0;
#ifdef PHREEQ98
	int AutoLoadOutputFile, CreateToC;
	int ProcessMessages, ShowProgress, ShowProgressWindow, ShowChart;
	int outputlinenr;
	int stop_calculations;
	char err_str98[80];
#endif
	/* utilities.cpp ------------------------------- */
	spinner                 = 0;
	// keycount;
	for (int i = 0; i < Keywords::KEY_COUNT_KEYWORDS; i++)
	{
		keycount.push_back(0);
	}

	return;
}
/*-----------------------------------------------------*/
Phreeqc::Phreeqc(const Phreeqc &src)
{
	this->phrq_io = src.phrq_io;
	this->init();
	this->initialize();
	InternalCopy(&src);
}
void
Phreeqc::InternalCopy(const Phreeqc *pSrc)
{
	// phrq_io
	/*
	if (io)
	{
		this->phrq_io = io;
	}
	else
	{
		this->phrq_io = &this->ioInstance;
	}
	*/

	same_model                      = FALSE;
	current_tc                      = pSrc->current_tc;
	current_pa                      = pSrc->current_pa;
	current_mu                      = pSrc->current_mu;
	mu_terms_in_logk                = pSrc->mu_terms_in_logk;

	MIN_LM = pSrc->MIN_LM;			    /* minimum log molality allowed before molality set to zero */
	LOG_ZERO_MOLALITY = pSrc->LOG_ZERO_MOLALITY;	/* molalities <= LOG_ZERO_MOLALITY are considered equal to zero */
	MIN_RELATED_LOG_ACTIVITY = pSrc->MIN_RELATED_LOG_ACTIVITY;
	MIN_TOTAL = pSrc->MIN_TOTAL;
	MIN_TOTAL_SS = pSrc->MIN_TOTAL_SS;
	MIN_RELATED_SURFACE = pSrc->MIN_RELATED_SURFACE;
	/* ----------------------------------------------------------------------
	*   STRUCTURES
	* ---------------------------------------------------------------------- */
/*
 *	 last model
 */
	//-- skip last model, accept init

/*
 *	 Initialize punch
 */
	//-- skip punch, accept init
	high_precision = pSrc->high_precision;

	Rxn_temperature_map = pSrc->Rxn_temperature_map;
	Rxn_pressure_map = pSrc->Rxn_pressure_map;

	/* ----------------------------------------------------------------------
	*   Surface
	* --------------------------------------------------------------------- */
	g_iterations               = -1;
	G_TOL                      = 1e-8;
	Rxn_surface_map = pSrc->Rxn_surface_map;
	// auto charge_group_map;
	/*
	change_surf_count          = 0;
	change_surf                = NULL;
	*/
	/* ----------------------------------------------------------------------
	*   Exchange
	* ---------------------------------------------------------------------- */
	Rxn_exchange_map = pSrc->Rxn_exchange_map;

	/* ----------------------------------------------------------------------
	*   Kinetics
	* ---------------------------------------------------------------------- */
	Rxn_kinetics_map = pSrc->Rxn_kinetics_map;

	/*----------------------------------------------------------------------
	*   Save
	*---------------------------------------------------------------------- */
	count_save_values          = 0;
	/*
	save_values                = NULL;	
	save_init(-1);             // set initial save values
	*/

	// auto use

	// copier structures
	//-- skip copier, accept init

	/*----------------------------------------------------------------------
	*   Inverse
	*---------------------------------------------------------------------- */
	
	/*
	inverse					= NULL;
	*/
	count_inverse			= 0;
	/*----------------------------------------------------------------------
	*   Mix
	*---------------------------------------------------------------------- */
	// Should be empty after each END
	// auto Rxn_mix_map;
	// auto Dispersion_mix_map;
	// auto Rxn_solution_mix_map;
	// auto Rxn_exchange_mix_map;
	// auto Rxn_gas_phase_mix_map;
	// auto Rxn_kinetics_mix_map;
	// auto Rxn_pp_assemblage_mix_map;
	// auto Rxn_ss_assemblage_mix_map;
	// auto Rxn_surface_mix_map;
	/*----------------------------------------------------------------------
	*   Irreversible reaction
	*---------------------------------------------------------------------- */
	Rxn_reaction_map = pSrc->Rxn_reaction_map;
	run_cells_one_step = pSrc->run_cells_one_step;
	/*----------------------------------------------------------------------
	*   Gas phase
	*---------------------------------------------------------------------- */
	Rxn_gas_phase_map = pSrc->Rxn_gas_phase_map;
	/*----------------------------------------------------------------------
	*   Solid solution
	*---------------------------------------------------------------------- */
	Rxn_ss_assemblage_map = pSrc->Rxn_ss_assemblage_map;
	/*----------------------------------------------------------------------
	*   Pure-phase assemblage
	*---------------------------------------------------------------------- */
	Rxn_pp_assemblage_map = pSrc->Rxn_pp_assemblage_map;
	/*----------------------------------------------------------------------
	*   Species_list
	*---------------------------------------------------------------------- */
	/*
	count_species_list      = 0;
	max_species_list        = 0;
	species_list            = NULL;
	*/
	/*----------------------------------------------------------------------
	*   Jacobian and Mass balance lists
	*---------------------------------------------------------------------- */
	/*
	count_sum_jacob0        = 0;
	max_sum_jacob0          = 0;
	sum_jacob0              = NULL;	
	count_sum_mb1           = 0;
	max_sum_mb1             = 0;
	sum_mb1                 = NULL;	
	count_sum_jacob1        = 0;
	max_sum_jacob1          = 0;
	sum_jacob1              = NULL;
	count_sum_mb2           = 0;
	max_sum_mb2             = 0;
	sum_mb2                 = NULL;
	count_sum_jacob2        = 0;
	max_sum_jacob2          = 0;
	sum_jacob2              = NULL;
	count_sum_delta         = 0;
	max_sum_delta           = 0;
	sum_delta               = NULL;
	*/
	/*----------------------------------------------------------------------
	*   Solution
	*---------------------------------------------------------------------- */
	Rxn_solution_map = pSrc->Rxn_solution_map;
	save_species = pSrc->save_species;
	// auto Rxn_solution_map;
	// auto unnumbered_solutions;
	/*----------------------------------------------------------------------
	*   Global solution
	*---------------------------------------------------------------------- */
	/*
	title_x                 = NULL;
	new_x                   = FALSE;
	description_x			= NULL;
	tc_x                    = 0;
	tk_x                    = 0;
	patm_x                  = 1;
	last_patm_x             = 1;
	numerical_fixed_volume  = false;
	force_numerical_fixed_volume = false;
	//switch_numerical        = false;
	ph_x                    = 0;
	solution_pe_x           = 0;
	mu_x                    = 0;
	ah2o_x                  = 1.0;
	density_x               = 0;
	total_h_x               = 0;
	total_o_x               = 0;
	cb_x                    = 0;
	total_ions_x            = 0;
	mass_water_aq_x         = 0;
	mass_water_surfaces_x   = 0;
	mass_water_bulk_x       = 0;
	units_x					= NULL;
	*/
	// auto pe_x
	// auto isotopes_x
	// auto default_pe_x
	/*
	dl_type_x                = cxxSurface::NO_DL;
	total_carbon             = 0;
	total_co2                = 0;
	total_alkalinity         = 0;
	gfw_water                = 0;
	step_x                   = 0;
	kin_time_x               = 0;
	*/
	/*----------------------------------------------------------------------
	*   Transport data
	*---------------------------------------------------------------------- */
	count_cells              = pSrc->count_cells;
	count_shifts             = pSrc->count_shifts;
	ishift                   = pSrc->ishift;
	bcon_first				 = pSrc->bcon_first;
	bcon_last				 = pSrc->bcon_last;
	correct_disp             = pSrc->correct_disp;
	tempr                    = pSrc->tempr;
	timest                   = pSrc->timest;
	simul_tr                 = pSrc->simul_tr;
	diffc                    = pSrc->diffc;
	heat_diffc               = pSrc->heat_diffc;
	cell                     = pSrc->cell;
	mcd_substeps             = pSrc->mcd_substeps;
	/* stag_data */
	memcpy(stag_data, pSrc->stag_data, sizeof(struct stag_data));
	print_modulus            = pSrc->print_modulus;
	punch_modulus            = pSrc->punch_modulus;
	dump_in                  = pSrc->dump_in;
	dump_modulus             = pSrc->dump_modulus;
	transport_warnings       = pSrc->transport_warnings;
	/* cell_data */
	if (count_cells > 0)
	{
		cell_data = (struct cell_data *) free_check_null(cell_data);
		cell_data = (struct cell_data *) PHRQ_malloc((size_t) (count_cells * sizeof(struct cell_data)));
		if (cell_data == NULL) malloc_error();
		memcpy(cell_data, pSrc->cell_data, ((size_t) (count_cells * sizeof(struct cell_data))));
	}
	old_cells = pSrc->old_cells;
	max_cells = pSrc->max_cells;
	all_cells = pSrc->all_cells;
	multi_Dflag              = pSrc->multi_Dflag;
	interlayer_Dflag         = pSrc->interlayer_Dflag;
	default_Dw               = pSrc->default_Dw;
	multi_Dpor               = pSrc->multi_Dpor;
	interlayer_Dpor          = pSrc->interlayer_Dpor;
	multi_Dpor_lim           = pSrc->multi_Dpor_lim;
	interlayer_Dpor_lim      = pSrc->interlayer_Dpor_lim;
	multi_Dn                 = pSrc->multi_Dn;
	interlayer_tortf         = pSrc->interlayer_tortf;
	cell_no                  = pSrc->cell_no;
	/*----------------------------------------------------------------------
	*   Advection data
	*---------------------------------------------------------------------- */
	count_ad_cells           = pSrc->count_ad_cells;
	count_ad_shifts          = pSrc->count_ad_shifts;
	print_ad_modulus         = pSrc->print_ad_modulus;
	punch_ad_modulus         = pSrc->punch_ad_modulus;
	/* advection_punch */
	if (count_ad_cells > 0)
	{
		advection_punch = (int *) free_check_null(advection_punch);
		advection_punch = (int *) PHRQ_malloc((size_t) (count_ad_cells * sizeof(int)));
		if (advection_punch == NULL) malloc_error();
		memcpy(advection_punch, pSrc->advection_punch, (size_t) (count_ad_cells * sizeof(int)));
	}
	/* advection_print */
	if (count_ad_cells > 0)
	{
		advection_print = (int *) free_check_null(advection_print);
		advection_print = (int *) PHRQ_malloc((size_t) (count_ad_cells * sizeof(int)));
		if (advection_print == NULL) malloc_error();
		memcpy(advection_print, pSrc->advection_print, (size_t) (count_ad_cells * sizeof(int)));
	}
	advection_kin_time       = pSrc->advection_kin_time;
	advection_kin_time_defined = pSrc->advection_kin_time_defined;
	advection_warnings       = pSrc->advection_warnings;
	/*----------------------------------------------------------------------
	*   Tidy data
	*---------------------------------------------------------------------- */
	/*
	new_model                = TRUE;
	new_exchange             = FALSE;
	new_pp_assemblage        = FALSE;
	new_surface              = FALSE;
	new_reaction             = FALSE;
	new_temperature          = FALSE;
	new_mix                  = FALSE;
	new_solution             = FALSE;
	new_gas_phase            = FALSE;
	new_inverse              = FALSE;
	new_punch                = FALSE;
	new_ss_assemblage        = FALSE;
	new_kinetics             = FALSE;
	new_copy                 = FALSE;
	new_pitzer               = FALSE;
	*/
	/*----------------------------------------------------------------------
	*   Elements
	*---------------------------------------------------------------------- */
	for (int i = 0; i < pSrc->count_elements; i++)
	{
		string_hsave(pSrc->elements[i]->name);
		struct element *elt_ptr = element_store(pSrc->elements[i]->name);
		elt_ptr->gfw = pSrc->elements[i]->gfw;
	}
	element_h_one = element_store("H(1)");
	/*
	elements                 = NULL;
	count_elements           = 0;
	max_elements             = MAX_ELEMENTS;
	element_h_one            = NULL;
	*/
	/*----------------------------------------------------------------------
	*   Element List
	*---------------------------------------------------------------------- */
	/*
	elt_list                 = NULL;
	count_elts               = 0;
	max_elts                 = MAX_ELTS;
	*/
	/*----------------------------------------------------------------------
	*   Species
	*---------------------------------------------------------------------- */
	/*
	logk                     = NULL;
	count_logk               = 0;
	max_logk                 = MAX_S;
	moles_per_kilogram_string= NULL;
	pe_string                = NULL;
	s                        = NULL;
	count_s                  = 0;
	max_s                    = MAX_S;
	// auto s_diff_layer;
	s_x                      = NULL;
	count_s_x                = 0;
	max_s_x                  = 0;
	s_h2o					= NULL;
	s_hplus					= NULL;
	s_h3oplus				= NULL;
	s_eminus				= NULL;
	s_co3					= NULL;
	s_h2					= NULL;
	s_o2					= NULL;
	*/	
	// logk
	for (int i = 0; i < pSrc->count_logk; i++)
	{
		char * name = string_duplicate(pSrc->logk[i]->name);
		struct logk *logk_ptr = logk_store(name, FALSE);
		free_check_null(name);
		memcpy(logk_ptr, pSrc->logk[i], sizeof(struct logk));
		logk_ptr->name = string_hsave(pSrc->logk[i]->name);
		logk_ptr->add_logk = NULL;
		if (logk_ptr->count_add_logk > 0)
		{
			logk_ptr->add_logk = (struct name_coef *) PHRQ_malloc((size_t) pSrc->logk[i]->count_add_logk * sizeof(struct name_coef));
			if (logk[i]->add_logk == NULL) malloc_error();
			for (int j = 0; j < logk_ptr->count_add_logk; j++)
			{
				logk_ptr->add_logk[j].coef = pSrc->logk[i]->add_logk[j].coef;
				logk_ptr->add_logk[j].name = string_hsave( pSrc->logk[i]->add_logk[j].name);
			}
		}	
	}
	// s, species
	for (int i = 0; i < pSrc->count_s; i++)
	{
		struct species *s_ptr = s_store(pSrc->s[i]->name, pSrc->s[i]->z, FALSE);
		memcpy(s_ptr, pSrc->s[i], sizeof(struct species));
		s_ptr->name = string_hsave(pSrc->s[i]->name);
		// fix up all pointers
		s_ptr->mole_balance = NULL;
		if (pSrc->s[i]->mole_balance != NULL)
		{
			s_ptr->mole_balance = string_hsave(pSrc->s[i]->mole_balance);
		}
		s_ptr->primary = NULL;
		s_ptr->secondary = NULL;
		//add_logk
		s_ptr->add_logk = NULL;
		if (s_ptr->count_add_logk > 0)
		{
			s_ptr->add_logk = (struct name_coef *) PHRQ_malloc((size_t) s_ptr->count_add_logk * sizeof(struct name_coef));
			if (s_ptr->add_logk == NULL) malloc_error();
			for (int j = 0; j < s_ptr->count_add_logk; j++)
			{
				s_ptr->add_logk[j].coef = pSrc->s[i]->add_logk[j].coef;
				s_ptr->add_logk[j].name = string_hsave( pSrc->s[i]->add_logk[j].name);
			}
		}
		//next_elt
		s_ptr->next_elt = NULL;
		if (pSrc->s[i]->next_elt)
		{
			cxxNameDouble next_elt(pSrc->s[i]->next_elt);
			s_ptr->next_elt = NameDouble2elt_list(next_elt);
		}
		//next_secondary
		s_ptr->next_secondary = NULL;
		if (pSrc->s[i]->next_secondary)
		{
			cxxNameDouble next_secondary(pSrc->s[i]->next_secondary);
			s_ptr->next_secondary = NameDouble2elt_list(next_secondary);
		}
		//next_sys_total
		s_ptr->next_sys_total = NULL;
		if (pSrc->s[i]->next_sys_total)
		{
			cxxNameDouble next_sys_total(pSrc->s[i]->next_sys_total);
			s_ptr->next_sys_total = NameDouble2elt_list(next_sys_total);
		}
		//rxn
		s_ptr->rxn = NULL;
		if (pSrc->s[i]->rxn != NULL)
		{
			cxxChemRxn rxn(pSrc->s[i]->rxn);
			s_ptr->rxn = cxxChemRxn2rxn(rxn);
			//s_ptr->rxn = rxn_copy_operator(pSrc->s[i]->rxn);
		}
		//rxn_s	
		s_ptr->rxn_s = NULL;
		if (pSrc->s[i]->rxn_s != NULL)
		{
			cxxChemRxn rxn_s(pSrc->s[i]->rxn_s);
			s_ptr->rxn_s = cxxChemRxn2rxn(rxn_s);
		}
		//rxn_x
		s_ptr->rxn_x = NULL;
		if (pSrc->s[i]->rxn_x != NULL)
		{
			cxxChemRxn rxn_x(pSrc->s[i]->rxn_x);
			s_ptr->rxn_x = cxxChemRxn2rxn(rxn_x);
		}
	}
	s_h2o					= s_search("H2O");
	s_hplus					= s_search("H+");
	s_h3oplus				= s_search("H3O+");
	s_eminus				= s_search("e-");
	s_co3					= s_search("CO3-2");
	s_h2					= s_search("H2");
	s_o2					= s_search("O2");
	/*----------------------------------------------------------------------
	*   Phases
	*---------------------------------------------------------------------- */
	/*
	phases					= NULL;
	count_phases            = 0;
	max_phases              = MAX_PHASES;
	*/
	for (int i = 0; i < pSrc->count_phases; i++)
	{
		struct phase *phase_ptr = phase_store(pSrc->phases[i]->name);
		memcpy(phase_ptr, pSrc->phases[i], sizeof(struct phase));
		// clean up pointers
		phase_ptr->name = string_hsave(pSrc->phases[i]->name);
		phase_ptr->formula = string_hsave(pSrc->phases[i]->formula);
		//add_logk
		phase_ptr->add_logk = NULL;
		if (phase_ptr->count_add_logk > 0)
		{
			phase_ptr->add_logk = (struct name_coef *) PHRQ_malloc((size_t) pSrc->phases[i]->count_add_logk * sizeof(struct name_coef));
			if (phase_ptr->add_logk == NULL) malloc_error();
			for (int j = 0; j < phase_ptr->count_add_logk; j++)
			{
				phase_ptr->add_logk[j].coef = pSrc->phases[i]->add_logk[j].coef;
				phase_ptr->add_logk[j].name = string_hsave( pSrc->phases[i]->add_logk[j].name);
			}
		}
		//next_elt
		phase_ptr->next_elt = NULL;
		if (pSrc->phases[i]->next_elt)
		{
			cxxNameDouble next_elt(pSrc->phases[i]->next_elt);
			phase_ptr->next_elt = NameDouble2elt_list(next_elt);
		}
		//next_sys_total
		phase_ptr->next_sys_total = NULL;
		if (pSrc->phases[i]->next_sys_total)
		{
			cxxNameDouble next_sys_total(pSrc->phases[i]->next_sys_total);
			phase_ptr->next_sys_total = NameDouble2elt_list(next_sys_total);
		}
		//rxn
		phase_ptr->rxn = NULL;
		if (pSrc->phases[i]->rxn != NULL)
		{
			cxxChemRxn rxn(pSrc->phases[i]->rxn);
			phase_ptr->rxn = cxxChemRxn2rxn(rxn);
		}
		//rxn_s
		phase_ptr->rxn_s = NULL;
		if (pSrc->phases[i]->rxn_s != NULL)
		{
			cxxChemRxn rxn_s(pSrc->phases[i]->rxn_s);
			phase_ptr->rxn_s = cxxChemRxn2rxn(rxn_s);
		}
		//rxn_x
		phase_ptr->rxn_x = NULL;
		if (pSrc->phases[i]->rxn_x != NULL)
		{
			cxxChemRxn rxn_x(pSrc->phases[i]->rxn_x);
			phase_ptr->rxn_x = cxxChemRxn2rxn(rxn_x);	
		}
	}
	/*----------------------------------------------------------------------
	*   Master species
	*---------------------------------------------------------------------- */
	/*
	master                  = NULL;
	dbg_master              = NULL;
	count_master            = 0;
	max_master              = MAX_MASTER;
	*/
	count_master = pSrc->count_master;
	max_master = pSrc->max_master;
	master = (struct master **) free_check_null(master);
	master = (struct master **) PHRQ_malloc((size_t) max_master * sizeof(struct master *));
	if (master == NULL) malloc_error();
	dbg_master = master;
	for (int i = 0; i < count_master; i++)
	{
		master[i] = (struct master *) PHRQ_malloc( sizeof(struct master));
		if (master[i] == NULL) malloc_error();
		memcpy(master[i], pSrc->master[i], sizeof(struct master));
		// clean up pointers
		master[i]->gfw_formula = NULL;
		if (pSrc->master[i]->gfw_formula != NULL)
		{
			master[i]->gfw_formula = string_hsave(pSrc->master[i]->gfw_formula);
		}
		master[i]->elt = element_store(pSrc->master[i]->elt->name);
		master[i]->unknown = NULL;
		master[i]->s = s_store(pSrc->master[i]->s->name, pSrc->master[i]->s->z, false);
		//rxn_primary
		master[i]->rxn_primary = NULL;
		if (pSrc->master[i]->rxn_primary != NULL)
		{
			cxxChemRxn rxn_primary(pSrc->master[i]->rxn_primary);
			master[i]->rxn_primary = cxxChemRxn2rxn(rxn_primary);
		}
		//rxn_secondary
		master[i]->rxn_secondary = NULL;
		if (pSrc->master[i]->rxn_secondary != NULL)
		{
			cxxChemRxn rxn_secondary(pSrc->master[i]->rxn_secondary);
			master[i]->rxn_secondary = cxxChemRxn2rxn(rxn_secondary);	
		}
	}
	/*----------------------------------------------------------------------
	*   Unknowns
	*---------------------------------------------------------------------- */
	/*
	x                       = NULL;
	count_unknowns          = 0;
	max_unknowns            = 0;
	ah2o_unknown            = NULL;
	alkalinity_unknown      = NULL;
	carbon_unknown          = NULL;
	charge_balance_unknown  = NULL;
	exchange_unknown        = NULL;
	mass_hydrogen_unknown   = NULL;
	mass_oxygen_unknown     = NULL;
	mb_unknown              = NULL;
	mu_unknown              = NULL;
	pe_unknown              = NULL;
	ph_unknown              = NULL;
	pure_phase_unknown      = NULL;
	solution_phase_boundary_unknown = NULL;
	surface_unknown         = NULL;
	gas_unknown             = NULL;
	ss_unknown              = NULL;
	*/
	// auto gas_unknowns;
	/*----------------------------------------------------------------------
	*   Reaction work space
	*---------------------------------------------------------------------- */
	// struct trxn;	
	/*
	trxn.token				= 0;
	for (int i = 0; i < MAX_LOG_K_INDICES; i++)
	{
		trxn.logk[i] = 0;
	}
	for (int i = 0; i < 3; i++)
	{
		trxn.dz[i] = 0;
	}
	count_trxn              = 0;
	max_trxn                = MAX_TRXN;
	*/
	/*
	mb_unknowns             = NULL;
	count_mb_unknowns       = 0;
	max_mb_unknowns         = MAX_TRXN;
	*/
	/* ----------------------------------------------------------------------
	*   Print
	* ---------------------------------------------------------------------- */
	/*
	pr.all                  = TRUE;
	pr.initial_solutions    = TRUE;
	pr.initial_exchangers   = TRUE;
	pr.reactions            = TRUE;
	pr.gas_phase            = TRUE;
	pr.ss_assemblage        = TRUE;
	pr.pp_assemblage        = TRUE;
	pr.surface              = TRUE;
	pr.exchange             = TRUE;
	pr.kinetics             = TRUE;
	pr.totals               = TRUE;
	pr.eh                   = TRUE;
	pr.species              = TRUE;
	pr.saturation_indices   = TRUE;
	pr.irrev                = TRUE;
	pr.mix                  = TRUE;
	pr.reaction             = TRUE;
	pr.use                  = TRUE;
	pr.logfile              = FALSE;
	pr.punch                = TRUE;
	pr.status               = TRUE;
	pr.inverse              = TRUE;
	pr.dump                 = TRUE;
	pr.user_print           = TRUE;
	pr.headings             = TRUE;
	pr.user_graph           = TRUE;
	pr.echo_input           = TRUE;
	pr.warnings             = 100;
	pr.initial_isotopes     = TRUE;
	pr.isotope_ratios       = TRUE;
	pr.isotope_alphas       = TRUE;
	pr.hdf                  = FALSE;
	pr.alkalinity           = FALSE;
	*/
	pr = pSrc->pr;
	status_on               = pSrc->status_on;
	status_interval         = pSrc->status_interval;
	status_timer            = clock();
	count_warnings          = 0;
	/* ----------------------------------------------------------------------
	*   RATES
	* ---------------------------------------------------------------------- */
	/*
	rates                   = NULL;
	count_rates				= 0;
	rate_m					= 0;
	rate_m0					= 0;
	rate_time				= 0;
	rate_kin_time           = 1.0;
	rate_sim_time_start		= 0;
	rate_sim_time_end		= 0;
	rate_sim_time			= 0;
	rate_moles				= 0;
	initial_total_time		= 0;
	// auto rate_p
	count_rate_p            = 0;
	*/
	rates = (struct rate *) free_check_null(rates);
	count_rates = pSrc->count_rates;
	if (count_rates > 0)
	{
		rates = (struct rate *) PHRQ_malloc((size_t) count_rates * sizeof(struct rate));
		if (rates == NULL) malloc_error();
		for (int i = 0; i < count_rates; i++)
		{
			rates[i].name = string_hsave(pSrc->rates[i].name);
			rates[i].commands = string_duplicate(pSrc->rates[i].commands); 
			rates[i].new_def = TRUE;
			rates[i].linebase = NULL;
			rates[i].varbase = NULL;
			rates[i].loopbase = NULL;
		}
	}
	/* ----------------------------------------------------------------------
	*   USER PRINT COMMANDS
	* ---------------------------------------------------------------------- */
	/*
	user_print				= NULL;
	*/
	{
		user_print->name = NULL;
		user_print->commands = NULL;
		if (pSrc->user_print->commands != NULL)
		{
			user_print->commands = string_duplicate(pSrc->user_print->commands); 
		}
		user_print->new_def = TRUE;
		user_print->linebase = NULL;
		user_print->varbase = NULL;
		user_print->loopbase = NULL;
	}

	// For now, User Punch is not copied
#ifdef SKIP
	/*
		user_punch				= NULL;
	*/
	{
		user_punch->name = NULL;
		user_punch->commands = NULL;
		if (pSrc->user_punch->commands != NULL)
		{
			user_punch->commands = string_duplicate(pSrc->user_punch->commands); 
		}
		user_punch->new_def = TRUE;
		user_punch->linebase = NULL;
		user_punch->varbase = NULL;
		user_punch->loopbase = NULL;
	}	
	/*
	user_punch_headings		= NULL;
	user_punch_count_headings = 0;
	*/
	user_punch_count_headings = pSrc->user_punch_count_headings;
	if (user_punch_count_headings > 0)
	{
		user_punch_headings = (const char **) free_check_null(user_punch_headings);
		user_punch_headings = (const char **) PHRQ_malloc((size_t) user_punch_count_headings * sizeof(char *));
		if (user_punch_headings == NULL) malloc_error();
		for (int i = 0; i < user_punch_count_headings; i++)
		{
			user_punch_headings[i] = string_hsave(pSrc->user_punch_headings[i]);
		}
	}
#endif
	n_user_punch_index      = pSrc->n_user_punch_index;
	fpunchf_user_s_warning  = pSrc->fpunchf_user_s_warning;
	//fpunchf_user_buffer[0]  = 0;

#if defined PHREEQ98 
	struct rate *user_graph;
	char **user_graph_headings;
	int user_graph_count_headings;
#endif
#if defined MULTICHART
	// auto chart_handler;
	chart_handler.Set_io(phrq_io);
#endif
	/* ----------------------------------------------------------------------
	*   GLOBAL DECLARATIONS
	* ---------------------------------------------------------------------- */
	/*
	error_string            = NULL;
	simulation				= 0;
	int state               = INITIALIZE;
	reaction_step           = 0;
	transport_step          = 0;
	transport_start         = 0;
	advection_step          = 0;
	stop_program            = FALSE;
	incremental_reactions   = FALSE;
	count_strings           = 0;
	array					= NULL;
	delta					= NULL;
	residual				= NULL;
	input_error             = 0;
	next_keyword            = Keywords::KEY_NONE;
	parse_error             = 0;
	paren_count             = 0;
	iterations              = 0;
	gamma_iterations        = 0;
	run_reactions_iterations= 0;
	max_line				= MAX_LINE;
	line                    = NULL;
	line_save				= NULL;
	LOG_10                  = log(10.0);
	debug_model             = FALSE;
	debug_prep              = FALSE;
	debug_set               = FALSE;
	debug_diffuse_layer     = FALSE;
	debug_inverse           = FALSE;
	*/
	inv_tol_default         = pSrc->inv_tol_default;
	itmax                   = pSrc->itmax;
	max_tries               = pSrc->max_tries;
	ineq_tol                = pSrc->ineq_tol;
	convergence_tolerance   = pSrc->convergence_tolerance;
	step_size				= pSrc->step_size;
	pe_step_size			= pSrc->pe_step_size;
	step_size_now           = step_size;
	pe_step_size_now        = pe_step_size;
	pp_scale				= pSrc->pp_scale;
	pp_column_scale			= pSrc->pp_column_scale;
	diagonal_scale			= pSrc->diagonal_scale;
	mass_water_switch		= pSrc->mass_water_switch;
	delay_mass_water		= pSrc->delay_mass_water;
	equi_delay      		= pSrc->equi_delay;
	dampen_ah2o             = pSrc->dampen_ah2o;
	censor					= pSrc->censor;
	aqueous_only			= pSrc->aqueous_only;
	negative_concentrations = pSrc->negative_concentrations;
	calculating_deriv		= FALSE;
	numerical_deriv			= FALSE;
	count_total_steps       = 0;
	phast                   = FALSE;
	/*
	llnl_temp				= 0;
	llnl_count_temp			= 0;
	llnl_adh				= 0;
	llnl_count_adh			= 0;
	llnl_bdh				= 0;
	llnl_count_bdh			= 0;
	llnl_bdot				= 0;
	llnl_count_bdot			= 0;
	llnl_co2_coefs			= 0;
	llnl_count_co2_coefs	= 0;
	*/
	llnl_count_temp			= pSrc->llnl_count_temp;
	if (llnl_count_temp > 0)
	{
		llnl_temp = (LDBLE *) free_check_null(llnl_temp);
		llnl_temp = (LDBLE *) PHRQ_malloc((size_t) llnl_count_temp * sizeof(LDBLE));
		if (llnl_temp == NULL) malloc_error();
		memcpy(llnl_temp, pSrc->llnl_temp, (size_t) llnl_count_temp * sizeof(LDBLE));
	}
	llnl_count_adh			= pSrc->llnl_count_adh;
	if (llnl_count_adh > 0)
	{
		llnl_adh = (LDBLE *) free_check_null(llnl_adh);
		llnl_adh = (LDBLE *) PHRQ_malloc((size_t) llnl_count_adh * sizeof(LDBLE));
		if (llnl_adh == NULL) malloc_error();
		memcpy(llnl_adh, pSrc->llnl_adh, (size_t) llnl_count_adh * sizeof(LDBLE));
	}
	llnl_count_bdh			= pSrc->llnl_count_bdh;
	if (llnl_count_bdh > 0)
	{
		llnl_bdh = (LDBLE *) free_check_null(llnl_bdh);
		llnl_bdh = (LDBLE *) PHRQ_malloc((size_t) llnl_count_bdh * sizeof(LDBLE));
		if (llnl_bdh == NULL) malloc_error();
		memcpy(llnl_bdh, pSrc->llnl_bdh, (size_t) llnl_count_bdh * sizeof(LDBLE));
	}
	llnl_count_bdot			= pSrc->llnl_count_bdot;
	if (llnl_count_bdot > 0)
	{
		llnl_bdot = (LDBLE *) free_check_null(llnl_bdot);
		llnl_bdot = (LDBLE *) PHRQ_malloc((size_t) llnl_count_bdot * sizeof(LDBLE));
		if (llnl_bdot == NULL) malloc_error();
		memcpy(llnl_bdot, pSrc->llnl_bdot, (size_t) llnl_count_bdot * sizeof(LDBLE));
	}
	llnl_count_co2_coefs	= pSrc->llnl_count_co2_coefs;
	if (llnl_count_co2_coefs > 0)
	{
		llnl_co2_coefs = (LDBLE *) free_check_null(llnl_co2_coefs);
		llnl_co2_coefs = (LDBLE *) PHRQ_malloc((size_t) llnl_count_co2_coefs * sizeof(LDBLE));
		if (llnl_co2_coefs == NULL) malloc_error();
		memcpy(llnl_co2_coefs, pSrc->llnl_co2_coefs, (size_t) llnl_count_co2_coefs * sizeof(LDBLE));
	}

	// Not implemented for now
	//SelectedOutput_map = pSrc->SelectedOutput_map;
	SelectedOutput_map.clear();

	//selected_output_file_name = NULL;
	//dump_file_name			= NULL;
	//remove_unstable_phases  = FALSE;
	// auto screen_string;
	spread_length           = 10;
	/* ---------------------------------------------------------------------- */
	/*
	*   Hash definitions
	*/
	// auto strings_map;
#ifdef HASH
	// auto strings_hash;
#endif
	/*
	elements_hash_table     = NULL;
	species_hash_table      = NULL;
	phases_hash_table       = NULL;
	logk_hash_table         = NULL;
	master_isotope_hash_table = NULL;
	*/
	/* ----------------------------------------------------------------------
	*   ISOTOPES
	* ---------------------------------------------------------------------- */
	/*
	count_master_isotope	= 0;
	master_isotope			= NULL;
	max_master_isotope		= MAX_ELTS;
	*/
	for (int i = 0; i < pSrc->count_master_isotope; i++)
	{
		struct master_isotope *master_isotope_ptr = master_isotope_store(pSrc->master_isotope[i]->name, FALSE);
		memcpy(master_isotope_ptr, pSrc->master_isotope[i], sizeof(struct master_isotope));
		master_isotope_ptr->name = string_hsave(pSrc->master_isotope[i]->name);
		int n;
		master_isotope_ptr->master = NULL;
		if (pSrc->master_isotope[i]->master)
		{
			char * name = string_duplicate(pSrc->master_isotope[i]->master->elt->name);
			master_isotope_ptr->master = master_search(name, &n);
			free_check_null(name);
		}
		if (master_isotope_ptr->master == NULL)
		{
			//error_msg("Error in copy constructor for master_isotope.", STOP);
		}
		master_isotope_ptr->elt = NULL;
		if (pSrc->master_isotope[i]->elt)
		{
			master_isotope_ptr->elt = element_store(pSrc->master_isotope[i]->elt->name);
		}
		master_isotope_ptr->units = NULL;
		if (pSrc->master_isotope[i]->units)
		{
			master_isotope_ptr->units = string_hsave(pSrc->master_isotope[i]->units);
		}
	}
	initial_solution_isotopes = pSrc->initial_solution_isotopes;
	/*
	count_calculate_value	= 0;
	calculate_value			= NULL;
	max_calculate_value		= MAX_ELTS;
	calculate_value_hash_table = NULL;	
	*/
	for (int i = 0; i < pSrc->count_calculate_value; i++)
	{
		struct calculate_value *calculate_value_ptr = calculate_value_store(pSrc->calculate_value[i]->name, FALSE);
		memcpy(calculate_value_ptr, pSrc->calculate_value[i], sizeof(struct calculate_value));
		calculate_value_ptr->value = pSrc->calculate_value[i]->value;
		calculate_value_ptr->commands = NULL;
		if (pSrc->calculate_value[i]->commands)
		{
			calculate_value_ptr->commands = string_duplicate(pSrc->calculate_value[i]->commands);
		}
		calculate_value_ptr->new_def = TRUE;
		calculate_value_ptr->calculated = FALSE;
		calculate_value_ptr->linebase = NULL;
		calculate_value_ptr->varbase = NULL;
		calculate_value_ptr->loopbase = NULL;
	}
	/*
	count_isotope_ratio		= 0;
	isotope_ratio			= 0;
	max_isotope_ratio		= MAX_ELTS;
	isotope_ratio_hash_table = 0;	
	*/
	for (int i = 0; i < pSrc->count_isotope_ratio; i++)
	{
		struct isotope_ratio *isotope_ratio_ptr = isotope_ratio_store(pSrc->isotope_ratio[i]->name, FALSE);
		isotope_ratio_ptr->name = string_hsave(pSrc->isotope_ratio[i]->name);
		isotope_ratio_ptr->isotope_name = string_hsave(pSrc->isotope_ratio[i]->isotope_name);
		isotope_ratio_ptr->ratio = pSrc->isotope_ratio[i]->ratio;
		isotope_ratio_ptr->converted_ratio = pSrc->isotope_ratio[i]->converted_ratio;
	}
	/*
	count_isotope_alpha		= 0;
	isotope_alpha			= 0;
	max_isotope_alpha		= MAX_ELTS;
	isotope_alpha_hash_table = 0;
	*/
	for (int i = 0; i < pSrc->count_isotope_alpha; i++)
	{
		struct isotope_alpha *isotope_alpha_ptr = isotope_alpha_store(pSrc->isotope_alpha[i]->name, FALSE);
		isotope_alpha_ptr->named_logk = string_hsave(pSrc->isotope_alpha[i]->named_logk);
		isotope_alpha_ptr->value = pSrc->isotope_alpha[i]->value;
	}

	phreeqc_mpi_myself		= 0;
	first_read_input		= TRUE;
	user_database			= string_duplicate(pSrc->user_database);
	//have_punch_name			= pSrc->have_punch_name;
	print_density		    = pSrc->print_density;
#ifdef SKIP
	zeros                   = NULL;	
	zeros_max			    = 1;
	cell_pore_volume	    = 0;
	cell_volume			    = 0;
	cell_porosity		    = 0;
	cell_saturation		    = 0;
	sys                     = NULL;
	count_sys               = 0;
	max_sys                 = 0;
	sys_tot                 = 0;

	V_solutes               = 0.0;
	rho_0                   = 0;
	kappa_0                 = 0.0;
	p_sat                   = 0.0;
	eps_r                   = EPSILON;
	DH_A                    = 0.0;
	DH_B                    = 0.0;
	DH_Av                   = 0.0;
	QBrn                    = 0.0;
	ZBrn                    = 0.0;
	dgdP                    = 0.0;

	need_temp_msg           = 0;
	solution_mass           = 0;
	solution_volume         = 0;
	/* phqalloc.cpp ------------------------------- */
	s_pTail                 = NULL;
	/* Basic */
	basic_interpreter       = NULL;
	/* cl1.cpp ------------------------------- */
	x_arg                   = NULL; 
	res_arg                 = NULL; 
	scratch                 = NULL;
	x_arg_max               = 0; 
	res_arg_max             = 0; 
	scratch_max             = 0;
	/* dw.cpp ------------------------------- */
	/* COMMON /QQQQ/ */	
	Q0                      = 0;
	Q5                      = 0;
	GASCON                  = 0.461522e0;
	TZ                      = 647.073e0;
	AA                      = 1.e0;
	Z                       = 0;
	DZ                      = 0;
	Y                       = 0;
	G1                      = 11.e0;
	G2                      = 44.333333333333e0;
	GF                      = 3.5e0;
	B1                      = 0;
	B2                      = 0;
	B1T                     = 0;
	B2T                     = 0;
	B1TT                    = 0;
	B2TT                    = 0;
	/* gases.cpp ------------------------------- */
	a_aa_sum                = 0;
	b2                      = 0;
	b_sum                   = 0;
	R_TK                    = 0;
	/* input.cpp ------------------------------- */
	check_line_return       = 0;  
	reading_db              = FALSE;
	/* integrate.cpp ------------------------------- */
	midpoint_sv             = 0;
	z_global                = 0;
	xd_global               = 0;
	alpha_global            = 0;
	/* inverse.cpp ------------------------------- */
	max_row_count           = 50;
	max_column_count        = 50;
	carbon                  = FALSE;
	col_name                = NULL;
	row_name                = NULL;
	count_rows              = 0;
	count_optimize          = 0;
	col_phases              = 0;
	col_redox               = 0;
	col_epsilon             = 0;
	col_ph                  = 0;
	col_water               = 0;
	col_isotopes            = 0;
	col_phase_isotopes      = 0;
	row_mb                  = 0;
	row_fract               = 0;
	row_charge              = 0;
	row_carbon              = 0;
	row_isotopes            = 0;
	row_epsilon             = 0;
	row_isotope_epsilon     = 0;
	row_water               = 0;
	inv_zero                = NULL;
	array1                  = 0;
	inv_res                 = NULL;
	inv_delta1              = NULL;
	delta2                  = NULL;
	delta3                  = NULL;
	inv_cu                  = NULL;
	delta_save              = NULL;
	min_delta               = NULL;
	max_delta               = NULL;
	inv_iu                  = NULL;
	inv_is                  = NULL;
	klmd                    = 0;
	nklmd                   = 0;
	n2d                     = 0;
	kode                    = 0;
	iter                    = 0;
	toler                   = 0;
	error                   = 0;
	max_pct                 = 0;
	scaled_error            = 0;
	master_alk              = NULL;
	row_back                = NULL;
	col_back                = NULL;
	good                    = NULL;
	bad                     = NULL;
	minimal                 = NULL;
	max_good                = 0;
	max_bad                 = 0;
	max_minimal             = 0;
	count_good              = 0;
	count_bad               = 0;
	count_minimal           = 0;
	count_calls             = 0;
	soln_bits               = 0;
	phase_bits              = 0;
	current_bits            = 0;
	temp_bits               = 0;
	netpath_file            = NULL;
	count_inverse_models    = 0;
	count_pat_solutions     = 0;
	for (int i = 0; i < 32; i++)
	{
		min_position[i]     = 0;
		max_position[i]     = 0;
		now[i]              = 0;
	}
	/* kinetics.cpp ------------------------------- */
	count_pp = count_pg = count_ss = 0; 
	cvode_kinetics_ptr      = NULL;
	cvode_test              = FALSE;
	cvode_error             = FALSE;
	cvode_n_user            = -99;
	cvode_n_reactions       = -99;
	cvode_step_fraction     = 0.0;
	cvode_rate_sim_time     = 0.0;
	cvode_rate_sim_time_start = 0.0;
	cvode_last_good_time    = 0.0;
	cvode_prev_good_time    = 0.0;
	cvode_last_good_y       = NULL;
	cvode_prev_good_y       = NULL;
	kinetics_machEnv        = NULL;
	kinetics_y              = NULL;
	kinetics_abstol         = NULL;
	kinetics_cvode_mem      = NULL;
	cvode_pp_assemblage_save= NULL;
	cvode_ss_assemblage_save= NULL;
	m_original              = NULL;
	m_temp                  = NULL;
	rk_moles                = NULL;
	set_and_run_attempt     = 0;
	x0_moles                = NULL;
	/* model.cpp ------------------------------- */
	gas_in                  = FALSE;
	min_value               = 1e-10;
	normal                  = NULL;
	ineq_array              = NULL;
	res                     = NULL;
	cu                      = NULL;
	zero                    = NULL;
	delta1                  = NULL;
	iu                      = NULL;
	is                      = NULL;
	back_eq                 = NULL;
	normal_max              = 0;
	ineq_array_max          = 0;
	res_max                 = 0;
	cu_max                  = 0;
	zero_max                = 0;
	delta1_max              = 0;
	iu_max                  = 0;
	is_max                  = 0;
	back_eq_max             = 0;
	/* phrq_io_output.cpp ------------------------------- */
	forward_output_to_log   = 0;
	/* phreeqc_files.cpp ------------------------------- */
	default_data_base       = string_duplicate("phreeqc.dat");
#ifdef PHREEQ98
	int outputlinenr;
	char *LogFileNameC;
	char progress_str[512];
#endif
#endif
	/* Pitzer  */	
	pitzer_model			= pSrc->pitzer_model;
	sit_model				= pSrc->sit_model;
	pitzer_pe				= pSrc->pitzer_pe;

	//full_pitzer             = FALSE;
	//always_full_pitzer      = FALSE;
	//ICON					= TRUE;
	//IC                      = -1;
	//COSMOT                  = 0;
	//AW                      = 0;
	//VP                      = 0;
	//DW0                     = 0;

	ICON					= pSrc->ICON;
	/*
	pitz_params				= NULL;
	count_pitz_param		= 0;
	max_pitz_param			= 100;
	*/
	for (int i = 0; i < pSrc->count_pitz_param; i++)
	{
		pitz_param_store(pSrc->pitz_params[i], true);
	}

	// auto pitz_param_map
	/*
	theta_params			= 0;
	count_theta_param		= 0;
	max_theta_param			= 100;
	use_etheta				= TRUE;
	OTEMP					= -100.0;
	OPRESS					= -100.0;
	A0                      = 0;	
	spec                    = NULL;
	cations                 = NULL;
	anions                  = NULL;
	neutrals                = NULL;
	count_cations           = 0;
	count_anions            = 0;
	count_neutrals          = 0;
	MAXCATIONS              = 0;
	FIRSTANION              = 0;
	MAXNEUTRAL              = 0;
	mcb0                    = NULL;
	mcb1                    = NULL;
	mcc0                    = NULL;
	IPRSNT                  = NULL;
	M                       = NULL;
	LGAMMA                  = NULL;
	for (int i = 0; i < 23; i++)
	{
		BK[i]				= 0.0;
		DK[i]				= 0.0;
	}
	*/

#ifdef PHREEQ98
	int connect_simulations, graph_initial_solutions;
	int shifts_as_points;
	int chart_type;
	int ShowChart;
	int RowOffset, ColumnOffset;
#endif
	dummy                   = 0;
	/* print.cpp ------------------------------- */
	/*
	sformatf_buffer = (char *) PHRQ_malloc(256 * sizeof(char));
	if (sformatf_buffer == NULL) 
		malloc_error();
	sformatf_buffer_size = 256;
	*/
#ifdef PHREEQ98
	int colnr, rownr;
	int graph_initial_solutions;
	int prev_advection_step, prev_transport_step;	/*, prev_reaction_step */
	/* int shifts_as_points; */
	int chart_type;
	int AddSeries;
	int FirstCallToUSER_GRAPH;
#endif
	/* read.cpp */
	prev_next_char          = NULL;
#if defined PHREEQ98 
	int shifts_as_points;
#endif
	/* read_class.cxx */
	// auto dump_info
	// auto delete_info
	// auto run_info
	/*
	run_info.Set_io(phrq_io);
	*/
	/* readtr.cpp */
	// auto dump_file_name_cpp;
	/* sit.cpp ------------------------------- */
/*
	sit_params              = NULL;
	count_sit_param			= 0;
	max_sit_param			= 100;
	// auto sit_param_map	
	sit_A0                  = 0;
	sit_count_cations       = 0;
	sit_count_anions        = 0;
	sit_count_neutrals      = 0;
	sit_MAXCATIONS          = 0;
	sit_FIRSTANION          = 0;
	sit_MAXNEUTRAL          = 0;
	sit_IPRSNT              = NULL;
	sit_M                   = NULL;
	sit_LGAMMA              = NULL;
*/
	
	for (int i = 0; i < pSrc->count_sit_param; i++)
	{
		sit_param_store(pSrc->sit_params[i], true);
	}
	/* tidy.cpp ------------------------------- */
	//a0                      = 0;
	//a1                      = 0;
	//kc                      = 0;
	//kb                      = 0;
	/* tally.cpp ------------------------------- */
	//t_buffer                = NULL;
	//tally_count_component   = 0;
	//tally_table             = NULL;
	//count_tally_table_columns = 0;
	//count_tally_table_rows  = 0;

	/* transport.cpp ------------------------------- */
	/* storage is created and freed in tranport.cpp */
	sol_D                   = NULL;   
	sol_D_dbg               = NULL;
	J_ij                    = NULL;
	J_ij_il                 = NULL;
	J_ij_count_spec         = pSrc->J_ij_count_spec;
	m_s                     = NULL;
	count_m_s               = pSrc->count_m_s;
	tot1_h                  = pSrc->tot1_h;
	tot1_o                  = pSrc->tot1_o;
	tot2_h                  = pSrc->tot2_h;
	tot2_o                  = pSrc->tot2_o;
	diffc_max               = pSrc->diffc_max;
	diffc_tr                = pSrc->diffc_tr;
	J_ij_sum                = pSrc->J_ij_sum;
	transp_surf             = pSrc->transp_surf;
	heat_mix_array          = NULL;
	temp1                   = NULL;
	temp2                   = NULL;
	nmix                    = pSrc->nmix;
	heat_nmix               = pSrc->heat_nmix;
	heat_mix_f_imm          = pSrc->heat_mix_f_imm;
	heat_mix_f_m            = pSrc->heat_mix_f_m;
	warn_MCD_X              = pSrc->warn_MCD_X;
	warn_fixed_Surf         = pSrc->warn_fixed_Surf;

#ifdef PHREEQ98
	int AutoLoadOutputFile, CreateToC;
	int ProcessMessages, ShowProgress, ShowProgressWindow, ShowChart;
	int outputlinenr;
	int stop_calculations;
	char err_str98[80];
#endif
	/* utilities.cpp ------------------------------- */
	//spinner                 = 0;
	//// keycount;
	//for (int i = 0; i < Keywords::KEY_COUNT_KEYWORDS; i++)
	//{
	//	keycount.push_back(0);
	//}

	// make sure new_model gets set
	this->keycount[Keywords::KEY_SOLUTION_SPECIES] = 1;
	this->tidy_model();
	return;
}
// Operator overloaded using a member function
Phreeqc &Phreeqc::operator=(const Phreeqc &rhs) 
{
	if (this == &rhs)      // Same object?
		return *this; 

	// clean up this here
	this->clean_up();

	this->PHRQ_free_all();
	if (this->phrq_io == &this->ioInstance)
	{
		this->phrq_io->clear_istream();
		this->phrq_io->close_ostreams();
	}

	// copy Phreeqc object to this
	this->phrq_io = rhs.phrq_io;
	this->init();
	this->initialize();
	this->InternalCopy(&rhs);
	return *this;
}

int Phreeqc::next_user_number(Keywords::KEYWORDS key)
{
	switch (key)
	{
	case Keywords::KEY_REACTION_TEMPERATURE:
		return Utilities::Rxn_next_user_number(Rxn_temperature_map);
		break;
	case Keywords::KEY_REACTION_PRESSURE:
		return Utilities::Rxn_next_user_number(Rxn_pressure_map);
		break;
	case Keywords::KEY_SURFACE:
		return Utilities::Rxn_next_user_number(Rxn_surface_map);
		break;
	case Keywords::KEY_EXCHANGE:
		return Utilities::Rxn_next_user_number(Rxn_exchange_map);
		break;
	case Keywords::KEY_KINETICS:
		return Utilities::Rxn_next_user_number(Rxn_kinetics_map);
		break;
	case Keywords::KEY_MIX:
		return Utilities::Rxn_next_user_number(Rxn_mix_map);
		break;
	case Keywords::KEY_REACTION:
		return Utilities::Rxn_next_user_number(Rxn_reaction_map);
		break;
	case Keywords::KEY_GAS_PHASE:
		return Utilities::Rxn_next_user_number(Rxn_gas_phase_map);
		break;
	case Keywords::KEY_SOLID_SOLUTIONS:
		return Utilities::Rxn_next_user_number(Rxn_ss_assemblage_map);
		break;
	case Keywords::KEY_EQUILIBRIUM_PHASES:
		return Utilities::Rxn_next_user_number(Rxn_pp_assemblage_map);
		break;
	case Keywords::KEY_SOLUTION:
		return Utilities::Rxn_next_user_number(Rxn_solution_map);
		break;
	default:
		assert(false);
		return -999;
	}
}
