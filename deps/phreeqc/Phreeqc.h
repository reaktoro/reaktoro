#ifndef _INC_PHREEQC_H
#define _INC_PHREEQC_H
#if defined(WIN32)
#include <windows.h>
#endif
#if defined(WIN32_MEMORY_DEBUG)
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif
/* ----------------------------------------------------------------------
*   INCLUDE FILES
* ---------------------------------------------------------------------- */
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#ifdef HASH
#include <hash_map>
#endif
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <errno.h>
#include <float.h>
#include "phrqtype.h"
#include "cvdense.h"	
#include "runner.h"
#include "dumper.h"
#include "PHRQ_io.h"
#include "SelectedOutput.h"
#include "UserPunch.h"
#ifdef MULTICHART
#include "ChartHandler.h"
#endif
#include "Keywords.h"
#include "Pressure.h"
#include "cxxMix.h"
#include "Use.h"
#include "Surface.h"
#ifdef SWIG_SHARED_OBJ
#include "thread.h"
#endif

class cxxNameDouble;
class cxxKinetics;
//class cxxMix;
class cxxKineticsComp;
class cxxExchange;
class cxxExchComp;
class cxxGasPhase;
class cxxTemperature;
class cxxPPassemblage;
class cxxPPassemblageComp;
class cxxReaction;
class cxxSolution;
class cxxISolutionComp;
class cxxSolutionIsotope;
class cxxSSassemblage;
class cxxSS;
class cxxStorageBin;

#include "global_structures.h"
class PBasic;

class Phreeqc
{
public:
	Phreeqc(PHRQ_io *io = NULL);
	Phreeqc(const Phreeqc &src);
	void InternalCopy(const Phreeqc *pSrc);
	Phreeqc &operator=(const Phreeqc &rhs); 
	~Phreeqc(void);

public:
	//
	// Phreeqc class methods
	//

	// advection.cpp -------------------------------
	int advection(void);

	// basicsubs.cpp -------------------------------
	int basic_compile(char *commands, void **lnbase, void **vbase, void **lpbase);
	int basic_run(char *commands, void *lnbase, void *vbase, void *lpbase);
	void basic_free(void);
#ifdef IPHREEQC_NO_FORTRAN_MODULE
	double basic_callback(double x1, double x2, char * str);
#else
	double basic_callback(double x1, double x2, const char * str);
#endif
	void register_basic_callback(double ( *fcn)(double x1, double x2, const char *str, void *cookie), void *cookie1);
#ifdef IPHREEQC_NO_FORTRAN_MODULE
	void register_fortran_basic_callback(double ( *fcn)(double *x1, double *x2, char *str, size_t l));
#else
	void register_fortran_basic_callback(double ( *fcn)(double *x1, double *x2, const char *str, int l));
#endif

	LDBLE activity(const char *species_name);
	LDBLE activity_coefficient(const char *species_name);
	LDBLE log_activity_coefficient(const char *species_name);
	LDBLE aqueous_vm(const char *species_name);
	LDBLE diff_c(const char *species_name);
	LDBLE sa_declercq(double type, double sa, double d, double m, double m0, double gfw);
	LDBLE calc_SC(void);
	/* VP: Density Start */
	LDBLE calc_dens(void);
	/* VP: Density End */
	LDBLE calc_logk_n(const char *name);
	LDBLE calc_logk_p(const char *name);
	LDBLE calc_logk_s(const char *name);
	LDBLE calc_surface_charge(const char *surface_name);
	LDBLE diff_layer_total(const char *total_name, const char *surface_name);
	LDBLE edl_species(const char *surf_name, LDBLE * count, char ***names, LDBLE ** moles, LDBLE * area, LDBLE * thickness);
	int get_edl_species(cxxSurfaceCharge & charge_ref);
	LDBLE equi_phase(const char *phase_name);
	LDBLE equi_phase_delta(const char *phase_name);
	LDBLE equivalent_fraction(const char *name, LDBLE *eq, std::string &elt_name);
	LDBLE find_gas_comp(const char *gas_comp_name);
	LDBLE find_gas_p(void);
	LDBLE find_gas_vm(void);
	LDBLE find_misc1(const char *ss_name);
	LDBLE find_misc2(const char *ss_name);
	LDBLE find_ss_comp(const char *ss_comp_name);
	LDBLE get_calculate_value(const char *name);
	char * iso_unit(const char *total_name);
	LDBLE iso_value(const char *total_name);
	LDBLE kinetics_moles(const char *kinetics_name);
	LDBLE kinetics_moles_delta(const char *kinetics_name);
	LDBLE log_activity(const char *species_name);
	LDBLE log_molality(const char *species_name);
	LDBLE molality(const char *species_name);
	LDBLE pressure(void);
	LDBLE pr_pressure(const char *phase_name);
	LDBLE pr_phi(const char *phase_name);
	LDBLE saturation_ratio(const char *phase_name);
	int saturation_index(const char *phase_name, LDBLE * iap, LDBLE * si);
	int solution_number(void);
	LDBLE solution_sum_secondary(const char *total_name);
	LDBLE sum_match_gases(const char *stemplate, const char *name);
	LDBLE sum_match_species(const char *stemplate, const char *name);
	LDBLE sum_match_ss(const char *stemplate, const char *name);
	int match_elts_in_species(const char *name, const char *stemplate);
	int extract_bracket(char **string, char *bracket_string);
	LDBLE surf_total(const char *total_name, const char *surface_name);
	LDBLE surf_total_no_redox(const char *total_name, const char *surface_name);
	static int system_species_compare(const void *ptr1, const void *ptr2);
	LDBLE system_total(const char *total_name, LDBLE * count, char ***names,
		char ***types, LDBLE ** moles);
	std::string phase_formula(std::string phase_name, cxxNameDouble &stoichiometry);
	std::string species_formula(std::string phase_name, cxxNameDouble &stoichiometry);
	LDBLE list_ss(std::string ss_name, cxxNameDouble &composition);
	int system_total_elements(void);
	int system_total_si(void);
	int system_total_aq(void);
	int system_total_ex(void);
	int system_total_surf(void);
	int system_total_gas(void);
	int system_total_equi(void);
	int system_total_ss(void);
	int system_total_elt(const char *total_name);
	int system_total_elt_secondary(const char *total_name);
	LDBLE total(const char *total_name);
	LDBLE total_mole(const char *total_name);
	int system_total_solids(cxxExchange *exchange_ptr,
		cxxPPassemblage *pp_assemblage_ptr,
		cxxGasPhase *gas_phase_ptr,
		cxxSSassemblage *ss_assemblage_ptr,
		cxxSurface *surface_ptr);

	static LDBLE f_rho(LDBLE rho_old, void *cookie);
	static LDBLE f_Vm(LDBLE v1, void *cookie);
	LDBLE calc_solution_volume(void);

	// chart.cpp
#if defined PHREEQ98 
	void DeleteCurves(void);
	void ExtractCurveInfo(char *line, int curvenr);
	void  GridChar(char *s, char *a);
	void MallocCurves(int nc, int ncxy);
	int OpenCSVFile(char file_name[MAX_LENGTH]);
	void SaveCurvesToFile(char file_name[MAX_LENGTH]);
	void PlotXY(char *x, char *y);
	void ReallocCurves(int new_nc);
	void ReallocCurveXY(int i);
	void SetAxisScale(char *a, int j, char *token, int true_);
	void SetAxisTitles(char *s, int i);
	void SetChartTitle(char *s);
	void start_chart(bool end);
#endif

	// cl1.cpp -------------------------------
	int cl1(int k, int l, int m, int n,
		int nklmd, int n2d,
		LDBLE * q,
		int *kode, LDBLE toler,
		int *iter, LDBLE * x, LDBLE * res, LDBLE * error,
		LDBLE * cu, int *iu, int *s, int check);
	void cl1_space(int check, int n2d, int klm, int nklmd);

	// cl1mp.cpp -------------------------------
	int cl1mp(int k, int l, int m, int n,
		int nklmd, int n2d,
		LDBLE * q_arg,
		int *kode, LDBLE toler,
		int *iter, LDBLE * x_arg, LDBLE * res_arg, LDBLE * error,
		LDBLE * cu_arg, int *iu, int *s, int check, LDBLE censor_arg);

	// class_main.cpp -------------------------------
	int write_banner(void);

	/* default.cpp */
public:
	int close_input_files(void);
	int close_output_files(void);
	static int istream_getc(void *cookie);
	int process_file_names(int argc, char *argv[], std::istream **db_cookie,
		std::istream **input_cookie, int log);

	/* PHRQ_io_output.cpp */
	void screen_msg(const char * str);

	void echo_msg(const char *err_str);
	int warning_msg(const char *err_str);
	void set_forward_output_to_log(int value);
	int get_forward_output_to_log(void);

	// dump_ostream
	bool dump_open(const char *file_name);
	void dump_flush(void);
	void dump_close(void);
	void dump_msg(const char * str);

	// log_ostream
	bool log_open(const char *file_name);
	void log_flush(void);
	void log_close(void);
	void log_msg(const char * str);

	// error_ostream
	bool error_open(const char *file_name);
	void error_flush(void);
	void error_close(void);
	void error_msg(const char * str, bool stop=false);

	// output_ostream
	bool output_open(const char *file_name);
	void output_flush(void);
	void output_close(void);
	void output_msg(const char * str);

	// punch_ostream
	bool punch_open(const char *file_name, int n_user);
	void punch_flush(void);
	void punch_close(void);
	void punch_msg(const char * str);

	void fpunchf_heading(const char *name);
	void fpunchf(const char *name, const char *format, double d);
	void fpunchf(const char *name, const char *format, char * d);
	void fpunchf(const char *name, const char *format, int d);
	void fpunchf_user(int user_index, const char *format, double d);
	void fpunchf_user(int user_index, const char *format, char * d);
	int fpunchf_end_row(const char *format);
#ifdef SKIP
	// dw.cpp -------------------------------
	int BB(LDBLE T);
	LDBLE PS(LDBLE T);
	LDBLE VLEST(LDBLE T);
	int DFIND(LDBLE * DOUT, LDBLE P, LDBLE D, LDBLE T);
	int QQ(LDBLE T, LDBLE D);
	LDBLE BASE(LDBLE D);
#endif
	// input.cpp -------------------------------
	int reading_database(void);
	void set_reading_database(int reading_database);
	int check_line(const char *string, int allow_empty, int allow_eof,
		int allow_keyword, int print);
	int add_char_to_line(int *i, char c);
	int check_line_impl(const char *string, int allow_empty,
		int allow_eof, int allow_keyword, int print);
	int get_line(void);
	int get_logical_line(void *cookie, int *l);
	int read_database(void);
	int run_simulations(void);

	// integrate.cpp -------------------------------
	int calc_all_g(void);
	int calc_init_g(void);
	int initial_surface_water(void);
	int sum_diffuse_layer(cxxSurfaceCharge *surface_charge_ptr1);
	int calc_all_donnan(void);
	int calc_init_donnan(void);
	LDBLE g_function(LDBLE x_value);
	LDBLE midpnt(LDBLE x1, LDBLE x2, int n);
	void polint(LDBLE * xa, LDBLE * ya, int n, LDBLE xv, LDBLE * yv,
		LDBLE * dy);
	LDBLE qromb_midpnt(cxxSurfaceCharge *charge_ptr, LDBLE x1, LDBLE x2);
	LDBLE calc_psi_avg(cxxSurfaceCharge *charge_ptr, LDBLE surf_chrg_eq);
	int calc_all_donnan_music(void);
	int calc_init_donnan_music(void);

	// inverse.cpp -------------------------------
	int inverse_models(void);
	int add_to_file(const char *filename, const char *string);
	int bit_print(unsigned long bits, int l);
	int carbon_derivs(struct inverse *inv_ptr);
	int check_isotopes(struct inverse *inv_ptr);
	int check_solns(struct inverse *inv_ptr);
	int count_isotope_unknowns(struct inverse *inv_ptr,
		struct isotope **isotope_unknowns);
	cxxSolutionIsotope *get_isotope(cxxSolution *solution_ptr, const char *elt);
	LDBLE get_inv_total(cxxSolution *solution_ptr, const char *elt);
	int isotope_balance_equation(struct inverse *inv_ptr, int row, int n);
	int post_mortem(void);
	bool test_cl1_solution(void);
	unsigned long get_bits(unsigned long bits, int position, int number);
	unsigned long minimal_solve(struct inverse *inv_ptr,
		unsigned long minimal_bits);
	void dump_netpath(struct inverse *inv_ptr);
	int dump_netpath_pat(struct inverse *inv_ptr);
	int next_set_phases(struct inverse *inv_ptr, int first_of_model_size,
		int model_size);
	int phase_isotope_inequalities(struct inverse *inv_ptr);
	int print_model(struct inverse *inv_ptr);
	int punch_model_heading(struct inverse *inv_ptr);
	int punch_model(struct inverse *inv_ptr);
	void print_isotope(FILE * netpath_file, cxxSolution *solution_ptr,
		const char *elt, const char *string);
	void print_total(FILE * netpath_file, cxxSolution *solution_ptr,
		const char *elt, const char *string);
	void print_total_multi(FILE * netpath_file, cxxSolution *solution_ptr,
		const char *string, const char *elt0,
		const char *elt1, const char *elt2, const char *elt3,
		const char *elt4);

	void print_total_pat(FILE * netpath_file, const char *elt,
		const char *string);
	int range(struct inverse *inv_ptr, unsigned long cur_bits);
	int save_bad(unsigned long bits);
	int save_good(unsigned long bits);
	int save_minimal(unsigned long bits);
	unsigned long set_bit(unsigned long bits, int position, int value);
	int setup_inverse(struct inverse *inv_ptr);
	int set_initial_solution(int n_user_old, int n_user_new);
	int set_ph_c(struct inverse *inv_ptr,
		int i, cxxSolution *soln_ptr_orig,	int n_user_new,
		LDBLE d_alk, LDBLE ph_factor, LDBLE alk_factor);
	int shrink(struct inverse *inv_ptr, LDBLE * array_in,
		LDBLE * array_out, int *k, int *l, int *m, int *n,
		unsigned long cur_bits, LDBLE * delta_l, int *col_back_l,
		int *row_back_l);
	int solve_inverse(struct inverse *inv_ptr);
	int solve_with_mask(struct inverse *inv_ptr, unsigned long cur_bits);
	int subset_bad(unsigned long bits);
	int subset_minimal(unsigned long bits);
	int superset_minimal(unsigned long bits);
	int write_optimize_names(struct inverse *inv_ptr);

	// isotopes.cpp -------------------------------
	int add_isotopes(cxxSolution &solution_ptr);
	int calculate_values(void);
	int calculate_isotope_moles(struct element *elt_ptr,
		cxxSolution *solution_ptr, LDBLE total_moles);
	LDBLE convert_isotope(struct master_isotope *master_isotope_ptr, LDBLE ratio);
	int from_pcil(struct master_isotope *master_isotope_ptr);
	int from_permil(struct master_isotope *master_isotope_ptr, LDBLE major_total);
	int from_pct(struct master_isotope *master_isotope_ptr, LDBLE major_total);
	int from_tu(struct master_isotope *master_isotope_ptr);
	struct calculate_value *calculate_value_alloc(void);
	int calculate_value_free(struct calculate_value *calculate_value_ptr);
	struct calculate_value *calculate_value_search(const char *name);
	struct calculate_value *calculate_value_store(const char *name,
		int replace_if_found);
	struct isotope_alpha *isotope_alpha_alloc(void);
	struct isotope_alpha *isotope_alpha_search(const char *name);
	struct isotope_alpha *isotope_alpha_store(const char *name,
		int replace_if_found);
	struct isotope_ratio *isotope_ratio_alloc(void);
	struct isotope_ratio *isotope_ratio_search(const char *name);
	struct isotope_ratio *isotope_ratio_store(const char *name,
		int replace_if_found);
	struct master_isotope *master_isotope_store(const char *name,
		int replace_if_found);
	struct master_isotope *master_isotope_alloc(void);
	struct master_isotope *master_isotope_search(const char *name);
	int print_initial_solution_isotopes(void);
	int print_isotope_ratios(void);
	int print_isotope_alphas(void);
	int punch_isotopes(void);
	int punch_calculate_values(void);
	int read_calculate_values(void);
	int read_isotopes(void);
	int read_isotope_ratios(void);
	int read_isotope_alphas(void);
	int calculate_value_init(struct calculate_value *calculate_value_ptr);
	int isotope_alpha_init(struct isotope_alpha *isotope_alpha_ptr);
	int isotope_ratio_init(struct isotope_ratio *isotope_ratio_ptr);
	int master_isotope_init(struct master_isotope *master_isotope_ptr);

	// kinetics.cpp -------------------------------
	void cvode_init(void);
	bool cvode_update_reactants(int i, int nsaver, bool save_it);
	int run_reactions(int i, LDBLE kin_time, int use_mix, LDBLE step_fraction);
	int set_and_run(int i, int use_mix, int use_kinetics, int nsaver,
		LDBLE step_fraction);
	int set_and_run_wrapper(int i, int use_mix, int use_kinetics, int nsaver,
		LDBLE step_fraction);
	int set_advection(int i, int use_mix, int use_kinetics, int nsaver);
	int free_cvode(void);
public:
	static void f(integertype N, realtype t, N_Vector y, N_Vector ydot,
		void *f_data);
	static void Jac(integertype N, DenseMat J, RhsFn f, void *f_data, realtype t,
		N_Vector y, N_Vector fy, N_Vector ewt, realtype h,
		realtype uround, void *jac_data, long int *nfePtr,
		N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

	int calc_final_kinetic_reaction(cxxKinetics *kinetics_ptr);
	int calc_kinetic_reaction(cxxKinetics *kinetics_ptr,
		LDBLE time_step);
	bool limit_rates(cxxKinetics *kinetics_ptr);
	int rk_kinetics(int i, LDBLE kin_time, int use_mix, int nsaver,
		LDBLE step_fraction);
	int set_reaction(int i, int use_mix, int use_kinetics);
	int set_transport(int i, int use_mix, int use_kinetics, int nsaver);
	int store_get_equi_reactants(int k, int kin_end);

	// mainsubs.cpp  -------------------------------
	std::ifstream * open_input_stream(char *query, char *default_name, std::ios_base::openmode mode, bool batch);
	std::ofstream * open_output_stream(char *query, char *default_name, std::ios_base::openmode mode, bool batch);
	int copy_entities(void);
	void do_mixes(void);
	void initialize(void);
	int initial_exchangers(int print);
	int initial_gas_phases(int print);
	int initial_solutions(int print);
	int solution_mix(void);
	int step_save_exch(int n_user);
	int step_save_surf(int n_user);
	int initial_surfaces(int print);
	int reactions(void);
	int saver(void);
	int xsolution_save(int k_user);
	int xexchange_save(int n_user);
	int xgas_save(int n_user);
	int xpp_assemblage_save(int n_user);
	int xss_assemblage_save(int n_user);
	int xsurface_save(int n_user);
	int do_initialize(void);
	int do_status(void);
	void save_init(int i);
	int copy_use(int i);
	int set_use(void);

	// model.cpp -------------------------------
	int check_residuals(void);
	int free_model_allocs(void);
	int ineq(int kode);
	int model(void);
	int jacobian_sums(void);
	int mb_gases(void);
	int mb_ss(void);
	int mb_sums(void);
	int molalities(int allow_overflow);
	int reset(void);
	int residuals(void);
	int set(int initial);
	int sum_species(void);
	int surface_model(void);
	LDBLE ss_root(LDBLE a0, LDBLE a1, LDBLE kc, LDBLE kb, LDBLE xcaq,
		LDBLE xbaq);
	LDBLE ss_halve(LDBLE a0, LDBLE a1, LDBLE x0, LDBLE x1, LDBLE kc,
		LDBLE kb, LDBLE xcaq, LDBLE xbaq);
	LDBLE ss_f(LDBLE xb, LDBLE a0, LDBLE a1, LDBLE kc, LDBLE kb,
		LDBLE xcaq, LDBLE xbaq);
	int numerical_jacobian(void);
	void set_inert_moles(void);
	void unset_inert_moles(void);
#ifdef SLNQ
	int add_trivial_eqns(int rows, int cols, LDBLE * matrix);
	//int slnq(int n, LDBLE * a, LDBLE * delta, int ncols, int print);
#endif
	int calc_gas_pressures(void);
	int calc_fixed_volume_gas_pressures(void);
	int calc_ss_fractions(void);
	int gammas(LDBLE mu);
	int initial_guesses(void);
	int revise_guesses(void);
	int ss_binary(cxxSS *ss_ptr);
	int ss_ideal(cxxSS *ss_ptr);
	void ineq_init(int max_row_count, int max_column_count);

	// parse.cpp -------------------------------
	int check_eqn(int association);
	int get_charge(char *charge, LDBLE * z);
	int get_elt(char **t_ptr, char *element, int *i);
	int get_elts_in_species(char **t_ptr, LDBLE coef);
	int get_num(char **t_ptr, LDBLE * num);
	int get_secondary_in_species(char **t_ptr, LDBLE coef);
	int parse_eq(char *eqn, struct elt_list **elt_ptr, int association);
	int get_coef(LDBLE * coef, char **eqnaddr);
	int get_secondary(char **t_ptr, char *element, int *i);
	int get_species(char **ptr);

	// phqalloc.cpp -------------------------------
public:
#if !defined(NDEBUG)
	void *PHRQ_malloc(size_t, const char *, int);
	void *PHRQ_calloc(size_t, size_t, const char *, int);
	void *PHRQ_realloc(void *, size_t, const char *, int);
#else
	void *PHRQ_malloc(size_t);
	void *PHRQ_calloc(size_t, size_t);
	void *PHRQ_realloc(void *, size_t);
#endif
	void PHRQ_free(void *ptr);
	void PHRQ_free_all(void);

public:

	// pitzer.cpp -------------------------------
	struct pitz_param *pitz_param_read(char *string, int n);
	void pitz_param_store(struct pitz_param *pzp_ptr, bool force_copy);
	void sit_param_store(struct pitz_param *pzp_ptr, bool force_copy);
	struct theta_param *theta_param_search(LDBLE zj, LDBLE zk);
	struct theta_param *theta_param_alloc(void);
	int theta_param_init(struct theta_param *theta_param_ptr);
	void pitzer_make_lists(void);
	int gammas_pz(void);
	int model_pz(void);
	int pitzer(void);
	int pitzer_clean_up(void);
	int pitzer_init(void);
	int pitzer_tidy(void);
	int read_pitzer(void);
	int set_pz(int initial);
	int calc_pitz_param(struct pitz_param *pz_ptr, LDBLE TK, LDBLE TR);
	int check_gammas_pz(void);	
#ifdef SKIP
	LDBLE DC(LDBLE T);
	int DW(LDBLE T);
#endif
	int ISPEC(const char *name);
	LDBLE G(LDBLE Y);
	LDBLE GP(LDBLE Y);
	int ETHETAS(LDBLE ZJ, LDBLE ZK, LDBLE I, LDBLE * etheta,
		LDBLE * ethetap);
	void ETHETA_PARAMS(LDBLE X, LDBLE& JAY, LDBLE& JPRIME );
	//int BDK(LDBLE X);
	int pitzer_initial_guesses(void);
	int pitzer_revise_guesses(void);
	int PTEMP(LDBLE TK);
	//LDBLE JAY(LDBLE X);
	//LDBLE JPRIME(LDBLE Y);
	int jacobian_pz(void);

	// pitzer_structures.cpp -------------------------------
	struct pitz_param *pitz_param_alloc(void);
	int pitz_param_init(struct pitz_param *pitz_param_ptr);
	struct pitz_param *pitz_param_duplicate(struct pitz_param *old_ptr);
	int pitz_param_copy(struct pitz_param *old_ptr,
	struct pitz_param *new_ptr);

	// prep.cpp -------------------------------
	int add_potential_factor(void);
	int add_cd_music_factors(int n);
	int add_surface_charge_balance(void);
	int add_cd_music_charge_balances(int i);
	int build_gas_phase(void);
	int build_fixed_volume_gas(void);
	int build_jacobian_sums(int k);
	int build_mb_sums(void);
	int build_min_exch(void);
	int build_model(void);
	int build_pure_phases(void);
	int build_ss_assemblage(void);
	int build_solution_phase_boundaries(void);
	int build_species_list(int n);
	int build_min_surface(void);
	LDBLE calc_lk_phase(phase * p_ptr, LDBLE TK, LDBLE pa);
	LDBLE calc_delta_v(reaction * r_ptr, bool phase);
	LDBLE calc_PR(std::vector<struct phase *> phase_ptrs, LDBLE P, LDBLE TK, LDBLE V_m);
	LDBLE calc_PR();
	int calc_vm(LDBLE tc, LDBLE pa);
	int change_hydrogen_in_elt_list(LDBLE charge);
	int clear(void);
	//int convert_units(struct solution *solution_ptr);
	int convert_units(cxxSolution *solution_ptr);
	LDBLE f_Vm(LDBLE v1);
	struct unknown *find_surface_charge_unknown(std::string &str_ptr, int plane);
	struct master **get_list_master_ptrs(char *ptr,
	struct master *master_ptr);
	int inout(void);
	int is_special(struct species *spec);
	int mb_for_species_aq(int n);
	int mb_for_species_ex(int n);
	int mb_for_species_surf(int n);
	int quick_setup(void);
	int resetup_master(void);
	int save_model(void);
	int setup_exchange(void);
	int setup_gas_phase(void);
	int setup_fixed_volume_gas(void);
	int setup_master_rxn(struct master **master_ptr_list,
		const std::string &pe_rxn);
	int setup_pure_phases(void);
	int adjust_setup_pure_phases(void);
	int setup_related_surface(void);
	int setup_ss_assemblage(void);
	int setup_solution(void);
	int adjust_setup_solution(void);
	int setup_surface(void);
	int setup_unknowns(void);
	int store_dn(int k, LDBLE * source, int row, LDBLE coef_in,
		LDBLE * gamma_source);
	int store_jacob(LDBLE * source, LDBLE * target, LDBLE coef);
	int store_jacob0(int row, int column, LDBLE coef);
	int store_mb(LDBLE * source, LDBLE * target, LDBLE coef);
	int store_mb_unknowns(struct unknown *unknown_ptr, LDBLE * LDBLE_ptr,
		LDBLE coef, LDBLE * gamma_ptr);
	int store_sum_deltas(LDBLE * source, LDBLE * target, LDBLE coef);
	int tidy_redox(void);
	struct master **unknown_alloc_master(void);
	int write_mb_eqn_x(void);
	int write_mb_for_species_list(int n);
	int write_mass_action_eqn_x(int stop);

	int check_same_model(void);
	int k_temp(LDBLE tc, LDBLE pa);
	LDBLE k_calc(LDBLE * logk, LDBLE tempk, LDBLE presPa);
	int prep(void);
	int reprep(void);
	int rewrite_master_to_secondary(struct master *master_ptr1,
	struct master *master_ptr2);
	int switch_bases(void);
	int write_phase_sys_total(int n);

	// print.cpp -------------------------------
	char *sformatf(const char *format, ...);
	int array_print(LDBLE * array_l, int row_count, int column_count,
		int max_column_count);
	int set_pr_in_false(void);
	int print_all(void);
	int print_exchange(void);
	int print_gas_phase(void);
	int print_master_reactions(void);
	int print_reaction(struct reaction *rxn_ptr);
	int print_species(void);
	int print_surface(void);
	int print_user_print(void);
	int punch_all(void);
	int print_alkalinity(void);
	int print_diffuse_layer(cxxSurfaceCharge *surface_charge_ptr);
	int print_eh(void);
	int print_reaction(void);
	int print_kinetics(void);
	int print_mix(void);
	int print_pp_assemblage(void);
	int print_ss_assemblage(void);
	int print_saturation_indices(void);
	int print_surface_cd_music(void);
	int print_totals(void);
	int print_using(void);
	/*int print_user_print(void);*/
	int punch_gas_phase(void);
	int punch_identifiers(void);
	int punch_kinetics(void);
	int punch_molalities(void);
	int punch_activities(void);
	int punch_pp_assemblage(void);
	int punch_ss_assemblage(void);
	int punch_saturation_indices(void);
	int punch_totals(void);
	int punch_user_punch(void);
#if defined MULTICHART
	int punch_user_graph(void);
#endif

	// read.cpp -------------------------------
	int read_input(void);
	int read_conc(cxxSolution *solution_ptr, int count_mass_balance, char *str);
	int *read_list_ints_range(char **ptr, int *count_ints, int positive,
		int *int_list);
	int read_log_k_only(char *ptr, LDBLE * log_k);
	int read_t_c_only(char *ptr, LDBLE *t_c);
	int read_p_c_only(char *ptr, LDBLE * p_c);
	int read_omega_only(char *ptr, LDBLE *omega);
	int read_number_description(char *ptr, int *n_user, int *n_user_end,
		char **description, int allow_negative=FALSE);
	int check_key(const char *str);
	int check_units(std::string &tot_units, bool alkalinity, bool check_compatibility,
			const char *default_units, bool print);
	int find_option(const char *item, int *n, const char **list, int count_list,
		int exact);
	int get_option(const char **opt_list, int count_opt_list, char **next_char);
	int get_true_false(char *string, int default_value);

	int add_psi_master_species(char *token);
	int read_advection(void);
	int read_analytical_expression_only(char *ptr, LDBLE * log_k);
	/* VP: Density Start */
	int read_millero_abcdef (char *ptr, LDBLE * abcdef);
	/* VP: Density End */
	int read_copy(void);
	int read_debug(void);
	int read_delta_h_only(char *ptr, LDBLE * delta_h,
		DELTA_H_UNIT * units);
	int read_aq_species_vm_parms(char *ptr, LDBLE * delta_v);
	int read_vm_only(char *ptr, LDBLE * delta_v,
		DELTA_V_UNIT * units);
	int read_phase_vm(char *ptr, LDBLE * delta_v,
		DELTA_V_UNIT * units);
	int read_llnl_aqueous_model_parameters(void);
	int read_exchange(void);
	int read_exchange_master_species(void);
	int read_exchange_species(void);
	int read_gas_phase(void);
	int read_incremental_reactions(void);
	int read_inverse(void);
	int read_inv_balances(struct inverse *inverse_ptr, char *next_char);
	int read_inv_isotopes(struct inverse *inverse_ptr, char *ptr);
	int read_inv_phases(struct inverse *inverse_ptr, char *next_char);
	int read_kinetics(void);
	int read_line_doubles(char *next_char, LDBLE ** d, int *count_d,
		int *count_alloc);
	int read_lines_doubles(char *next_char, LDBLE ** d, int *count_d,
		int *count_alloc, const char **opt_list,
		int count_opt_list, int *opt);
	LDBLE *read_list_doubles(char **ptr, int *count_doubles);
	int *read_list_ints(char **ptr, int *count_ints, int positive);
	int *read_list_t_f(char **ptr, int *count_ints);
	int read_master_species(void);
	int read_mix(void);
	int read_entity_mix(std::map<int, cxxMix> &mix_map);
	//int read_solution_mix(void);
	int read_named_logk(void);
	int read_phases(void);
	int read_print(void);
	int read_pp_assemblage(void);
	int read_rates(void);
	int read_reaction(void);
	int read_reaction_reactants(cxxReaction *reaction_ptr);
	int read_reaction_steps(cxxReaction *reaction_ptr);
	int read_solid_solutions(void);
	int read_temperature(void);
	int read_reaction_temps(struct temperature *temperature_ptr);
	int read_reaction_pressure(void);
	int read_reaction_pressure_raw(void);
	int read_save(void);
	int read_selected_output(void);
	int read_solution(void);
	int read_species(void);
	int read_surface(void);
	int read_surface_master_species(void);
	int read_surface_species(void);
	int read_use(void);
	int read_title(void);
	int read_user_print(void);
	int read_user_punch(void);
#if defined PHREEQ98 
	int read_user_graph(void);
#endif
#if defined MULTICHART
	int read_user_graph_handler();
#endif
	int next_keyword_or_option(const char **opt_list, int count_opt_list);
	int cleanup_after_parser(CParser &parser);

	// ReadClass.cxx
	int read_dump(void);
	int read_delete(void);
	int read_run_cells(void);
	int streamify_to_next_keyword(std::istringstream & lines);
	int dump_entities(void);
	int delete_entities(void);
	int run_as_cells(void);
	void dump_ostream(std::ostream& os);

	// readtr.cpp -------------------------------
	int read_transport(void);
	int dump(void);
	int dump_exchange(int k);
	int dump_gas_phase(int k);
	int dump_kinetics(int k);
	int dump_mix(int k);
	int dump_pp_assemblage(int k);
	int dump_reaction(int k);
	int dump_ss_assemblage(int k);
	int dump_solution(int k);
	int dump_surface(int k);
	int dump_cpp(void);
	int read_line_LDBLEs(char *next_char, LDBLE ** d, int *count_d,
		int *count_alloc);

	// sit.cpp -------------------------------
	int gammas_sit(void);
	int model_sit(void);
	int sit(void);
	int sit_clean_up(void);
	int sit_init(void);
	int sit_tidy(void);
	int read_sit(void);
	int set_sit(int initial);
	int calc_sit_param(struct pitz_param *pz_ptr, LDBLE TK, LDBLE TR);
	int check_gammas_sit(void);
	int sit_ISPEC(const char *name);
	/*int DH_AB (LDBLE TK, LDBLE *A, LDBLE *B);*/
	int sit_initial_guesses(void);
	int sit_revise_guesses(void);
	int PTEMP_SIT(LDBLE tk);
	void sit_make_lists(void);
	int jacobian_sit(void);

	// spread.cpp -------------------------------
	int read_solution_spread(void);
	int copy_token_tab(char *token_ptr, char **ptr, int *length);
	int get_option_string(const char **opt_list, int count_opt_list,
		char **next_char);
	int spread_row_free(struct spread_row *spread_row_ptr);
	int spread_row_to_solution(struct spread_row *heading,
		struct spread_row *units,
		struct spread_row *data,
		struct defaults defaults);
	struct spread_row *string_to_spread_row(char *string);
#ifdef PHREEQCI_GUI
	void add_row(struct spread_row *spread_row_ptr);
	void copy_defaults(struct defaults *dest_ptr,
	struct defaults *src_ptr);
	void free_spread(void);
	struct spread_row *copy_row(struct spread_row *spread_row_ptr);
#endif

	// step.cpp -------------------------------
	int step(LDBLE step_fraction);
	int xsolution_zero(void);
	int add_exchange(cxxExchange *exchange_ptr);
	int add_gas_phase(cxxGasPhase *gas_phase_ptr);
	int add_kinetics(cxxKinetics *kinetics_ptr);
	int add_mix(cxxMix * mix_ptr);
	int add_pp_assemblage(cxxPPassemblage *pp_assemblage_ptr);
	int add_reaction(cxxReaction *reaction_ptr, int step_number, LDBLE step_fraction);
	int add_ss_assemblage(cxxSSassemblage *ss_assemblage_ptr);
	int add_solution(cxxSolution *solution_ptr, LDBLE extensive,
		LDBLE intensive);
	int add_surface(cxxSurface *surface_ptr);
	int check_pp_assemblage(cxxPPassemblage *pp_assemblage_ptr);
	int gas_phase_check(cxxGasPhase *gas_phase_ptr);
	int pp_assemblage_check(cxxPPassemblage *pp_assemblage_ptr);
	int reaction_calc(cxxReaction *reaction_ptr);
	int solution_check(void);
	int ss_assemblage_check(cxxSSassemblage *ss_assemblage_ptr);

	// structures.cpp -------------------------------
	int clean_up(void);
	int reinitialize(void);
	int copier_add(struct copier *copier_ptr, int n_user, int start, int end);
	int copier_free(struct copier *copier_ptr);
	int copier_init(struct copier *copier_ptr);
	static int element_compare(const void *ptr1, const void *ptr2);
public:
	struct element *element_store(const char *element);
	int elt_list_combine(void);
	static int elt_list_compare(const void *ptr1, const void *ptr2);
protected:
	struct elt_list *elt_list_dup(struct elt_list *elt_list_ptr_old);
	int elt_list_print(struct elt_list *elt_list_ptr);
	struct elt_list *elt_list_save(void);
	cxxNameDouble elt_list_NameDouble(void);
	struct elt_list * NameDouble2elt_list(const cxxNameDouble &nd);
public:
	enum entity_type get_entity_enum(char *name);
	struct inverse *inverse_alloc(void);
	int inverse_delete(int i);
	static int inverse_isotope_compare(const void *ptr1, const void *ptr2);
	struct inverse *inverse_search(int n_user, int *n);
	int inverse_sort(void);
protected:
	struct logk *logk_alloc(void);
	int logk_copy2orig(struct logk *logk_ptr);
	struct logk *logk_store(char *name, int replace_if_found);
	struct logk *logk_search(const char *name);
	struct master *master_alloc(void);
	static int master_compare(const void *ptr1, const void *ptr2);
	int master_delete(char *ptr);
public:
	struct master *master_bsearch(const char *ptr);
	struct master *master_bsearch_primary(const char *ptr);
	struct master *master_bsearch_secondary(char *ptr);
	struct master *master_search(char *ptr, int *n);
	struct pe_data *pe_data_alloc(void);
public:
	struct pe_data *pe_data_dup(struct pe_data *pe_ptr_old);
	struct pe_data *pe_data_free(struct pe_data *pe_data_ptr);
protected:
	int pe_data_store(struct pe_data **pe, const char *token);
public:
	struct phase *phase_bsearch(const char *ptr, int *j, int print);
protected:
	static int phase_compare(const void *ptr1, const void *ptr2);
	int phase_delete(int i);
	struct phase *phase_store(const char *name);
public:
	struct rate *rate_bsearch(char *ptr, int *j);
	int rate_free(struct rate *rate_ptr);
	struct rate *rate_search(const char *name, int *n);
	int rate_sort(void);
	struct reaction *rxn_alloc(int ntokens);
	struct reaction *rxn_dup(struct reaction *rxn_ptr_old);
	struct reaction * cxxChemRxn2rxn(cxxChemRxn &cr);
	LDBLE rxn_find_coef(struct reaction *r_ptr, const char *str);
	int rxn_free(struct reaction *rxn_ptr);
	int rxn_print(struct reaction *rxn_ptr);
	static int s_compare(const void *ptr1, const void *ptr2);
	int s_delete(int i);
	struct species *s_search(const char *name);
	struct species *s_store(const char *name, LDBLE z, int replace_if_found);
protected:
	struct save_values *save_values_bsearch(struct save_values *k, int *n);
	static int save_values_compare(const void *ptr1, const void *ptr2);
	int save_values_sort(void);
	int save_values_store(struct save_values *s_v);
	static int isotope_compare(const void *ptr1, const void *ptr2);
	static int species_list_compare_alk(const void *ptr1, const void *ptr2);
	static int species_list_compare_master(const void *ptr1, const void *ptr2);
	int species_list_sort(void);
	struct Change_Surf *change_surf_alloc(int count);
public:
	struct master *surface_get_psi_master(const char *name, int plane);
	int system_duplicate(int i, int save_old);
	int trxn_add(struct reaction *r_ptr, LDBLE coef, int combine);
	int trxn_add(cxxChemRxn &r_ptr, LDBLE coef, int combine);
	int trxn_add_phase(struct reaction *r_ptr, LDBLE coef, int combine);
	int trxn_combine(void);
	int trxn_copy(struct reaction *rxn_ptr);
	LDBLE trxn_find_coef(const char *str, int start);
	int trxn_print(void);
	int trxn_reverse_k(void);
	int trxn_sort(void);
	int trxn_swap(const char *token);
	struct unknown *unknown_alloc(void);
	int unknown_delete(int i);
	int unknown_free(struct unknown *unknown_ptr);
	int entity_exists(const char *name, int n_user);
	static int inverse_compare(const void *ptr1, const void *ptr2);
	int inverse_free(struct inverse *inverse_ptr);
	static int kinetics_compare_int(const void *ptr1, const void *ptr2);
	int logk_init(struct logk *logk_ptr);
	static int master_compare_string(const void *ptr1, const void *ptr2);
	int master_free(struct master *master_ptr);
	struct phase *phase_alloc(void);
	static int phase_compare_string(const void *ptr1, const void *ptr2);
	int phase_free(struct phase *phase_ptr);
	int phase_init(struct phase *phase_ptr);
	static int rate_compare(const void *ptr1, const void *ptr2);
	static int rate_compare_string(const void *ptr1, const void *ptr2);
	struct species *s_alloc(void);
	int s_free(struct species *s_ptr);
	int s_init(struct species *s_ptr);
	static int ss_assemblage_compare_int(const void *ptr1, const void *ptr2);
	static int solution_compare(const void *ptr1, const void *ptr2);
	static int solution_compare_int(const void *ptr1, const void *ptr2);
	static int species_list_compare(const void *ptr1, const void *ptr2);
	static int surface_compare_int(const void *ptr1, const void *ptr2);
	static int rxn_token_temp_compare(const void *ptr1, const void *ptr2);
	int trxn_multiply(LDBLE coef);

	struct elt_list * cxxNameDouble2elt_list(const cxxNameDouble * nd);
	struct name_coef * cxxNameDouble2name_coef(const cxxNameDouble * nd);
	struct master_activity * cxxNameDouble2master_activity(const cxxNameDouble * nd);
	struct master * cxxNameDouble2surface_master(const cxxNameDouble * totals);

	void Use2cxxStorageBin(cxxStorageBin & sb);
	void phreeqc2cxxStorageBin(cxxStorageBin & sb);
	void phreeqc2cxxStorageBin(cxxStorageBin & sb, int n);
	void cxxStorageBin2phreeqc(cxxStorageBin & sb, int n);
	void cxxStorageBin2phreeqc(cxxStorageBin & sb);

	/* tally.cpp */
	void add_all_components_tally(void);
	int build_tally_table(void);
	int calc_dummy_kinetic_reaction_tally(cxxKinetics *kinetics_ptr);
	int diff_tally_table(void);
	int extend_tally_table(void);
	int free_tally_table(void);
	int fill_tally_table(int *n_user, int index_conservative, int n_buffer);
	int get_tally_table_rows_columns(int *rows, int *columns);
	int get_tally_table_column_heading(int column, int *type, char *string);
	int get_tally_table_row_heading(int column, char *string);
	int store_tally_table(LDBLE * array, int row_dim, int col_dim,
		LDBLE fill_factor);
	int zero_tally_table(void);
	int elt_list_to_tally_table(struct tally_buffer *buffer_ptr);
	int master_to_tally_table(struct tally_buffer *buffer_ptr);
	int get_all_components(void);
	int print_tally_table(void);
	int set_reaction_moles(int n_user, LDBLE moles);
	int set_reaction_temperature(int n_user, LDBLE tc);
	int set_kinetics_time(int n_user, LDBLE step);

	// tidy.cpp -------------------------------
	int add_other_logk(LDBLE * source_k, int count_add_logk,
	struct name_coef *add_logk);
	int add_logks(struct logk *logk_ptr, int repeats);
	LDBLE halve(LDBLE f(LDBLE x, void *), LDBLE x0, LDBLE x1, LDBLE tol);
	int replace_solids_gases(void);
	int ss_prep(LDBLE t, cxxSS *ss_ptr, int print);
	int select_log_k_expression(LDBLE * source_k, LDBLE * target_k);
	int slnq(int n, LDBLE * a, LDBLE * delta, int ncols, int print);
public:
	int tidy_punch(void);
	int tidy_model(void);
	int check_species_input(void);
	LDBLE coef_in_master(struct master *master_ptr);
	int phase_rxn_to_trxn(struct phase *phase_ptr,
	struct reaction *rxn_ptr);
	int reset_last_model(void);
	int rewrite_eqn_to_primary(void);
	int rewrite_eqn_to_secondary(void);
	int species_rxn_to_trxn(struct species *s_ptr);
	int tidy_logk(void);
	int tidy_exchange(void);
	int tidy_min_exchange(void);
	int tidy_kin_exchange(void);
	int tidy_gas_phase(void);
	int tidy_inverse(void);
	int tidy_isotopes(void);
	int tidy_isotope_ratios(void);
	int tidy_isotope_alphas(void);
	int tidy_kin_surface(void);
	int tidy_master_isotope(void);
	int tidy_min_surface(void);
	int tidy_phases(void);
	int tidy_pp_assemblage(void);
	int tidy_solutions(void);
	int tidy_ss_assemblage(void);
	int tidy_species(void);
	int tidy_surface(void);
	int scan(LDBLE f(LDBLE x, void *), LDBLE * xx0, LDBLE * xx1);
	static LDBLE f_spinodal(LDBLE x, void *);
	int solve_misc(LDBLE * xxc1, LDBLE * xxc2, LDBLE tol);
	int ss_calc_a0_a1(cxxSS *ss_ptr);

	// transport.cpp -------------------------------
	int transport(void);
	int set_initial_moles(int i);
	cxxSurface sum_surface_comp(cxxSurface *source1, LDBLE f1,
		cxxSurface *source2, std::string charge_name, LDBLE f2,
		LDBLE new_Dw);
	int reformat_surf(const char *comp_name, LDBLE fraction, const char *new_comp_name,
		LDBLE new_Dw, int cell);
	LDBLE viscosity(void);
	int multi_D(LDBLE DDt, int mobile_cell, int stagnant);
	int find_J(int icell, int jcell, LDBLE mixf, LDBLE DDt, int stagnant);
	int fill_spec(int cell_no);
	int fill_m_s(struct J_ij *J_ij, int J_ij_count_spec);
	static int sort_species_name(const void *ptr1, const void *ptr2);
	int disp_surf(LDBLE stagkin_time);
	int diff_stag_surf(int mobile_cell);
	int check_surfaces(cxxSurface *surface_ptr1, cxxSurface *surface_ptr2);
	cxxSurface mobile_surface_copy(cxxSurface *surface_old_ptr,
		int n_user_new,
		bool move_old);
	int init_mix(void);
	int init_heat_mix(int nmix);
	int heat_mix(int heat_nmix);
	int mix_stag(int i, LDBLE stagkin_time, int punch,
		LDBLE step_fraction_kin);

	// utilities.cpp -------------------------------
public:
	int add_elt_list(struct elt_list *elt_list_ptr, LDBLE coef);
	int add_elt_list_multi_surf(struct elt_list *elt_list_ptr, LDBLE coef, struct element *surf_elt_ptr);
	int add_elt_list(const cxxNameDouble & nd, LDBLE coef);
protected:
	int backspace_screen(int spaces);
	LDBLE calc_alk(struct reaction *rxn_ptr);
public:
	LDBLE calc_rho_0(LDBLE tc, LDBLE pa);
	LDBLE calc_dielectrics(LDBLE tc, LDBLE pa);
	int compute_gfw(const char *string, LDBLE * gfw);
#if defined PHREEQ98 
	int copy_title(char *token_ptr, char **ptr, int *length);
#endif
	int copy_token(char *token_ptr, char **ptr, int *length);
	int copy_token(std::string &token, char **ptr);
	int dup_print(const char *ptr, int emphasis);
	int equal(LDBLE a, LDBLE b, LDBLE eps);
public:
	void *free_check_null(void *ptr);
protected:
	void free_hash_strings(HashTable * Table);
	int get_token(char **eqnaddr, char *string, LDBLE * z, int *l);
	int hcreate_multi(unsigned Count, HashTable ** HashTable_ptr);
	void hdestroy_multi(HashTable * HashTable_ptr);
	ENTRY *hsearch_multi(HashTable * Table, ENTRY item, ACTION action);
	int islegit(const char c);
public:
	void malloc_error(void);
protected:
	int parse_couple(char *token);
	int print_centered(const char *string);
public:
	static int replace(const char *str1, const char *str2, char *str);
	static bool replace(const char *str1, const char *str2, std::string & str);
	static int strcmp_nocase(const char *str1, const char *str2);
	static int strcmp_nocase_arg1(const char *str1, const char *str2);
protected:
	void space(void **ptr, int i, int *max, int struct_size);
	void squeeze_white(char *s_l);
	int status(int count, const char *str, bool kinetics = false);
	void str_tolower(char *str);
	void str_toupper(char *str);
public:
#if !defined(NDEBUG) && defined(WIN32_MEMORY_DEBUG)
	char *_string_duplicate(const char *token, const char *szFileName, int nLine);
#else
	char *string_duplicate(const char *token);
#endif
	const char *string_hsave(const char *str);
	void strings_map_clear();
#ifdef HASH
	void strings_hash_clear();
#endif
protected:
	char *string_pad(const char *str, int i);
	int string_trim(char *str);
	int string_trim_right(char *str);
	int string_trim_left(char *str);
	static LDBLE under(LDBLE xval);
	void zero_double(LDBLE * target, int n);
	int get_input_errors(void);
#ifdef PHREEQ98
	void AddToCEntry(char *a, int l, int i);
	void ApplicationProcessMessages(void);
	int copy_title(char *token_ptr, char **ptr, int *length);
	extern int clean_up_null(void);
#endif
	int isamong(char c, const char *s_l);
	Address Hash_multi(HashTable * Table, const char *Key);
	void ExpandTable_multi(HashTable * Table);
public:
	int main_method(int argc, char *argv[]);
	void set_phast(int);
	int next_user_number(Keywords::KEYWORDS key);
	size_t list_components(std::list<std::string> &list_c);
	PHRQ_io * Get_phrq_io(void) {return this->phrq_io;}
	void Set_run_cells_one_step(const bool tf) {this->run_cells_one_step = tf;}


	std::map<int, cxxSolution> & Get_Rxn_solution_map() {return this->Rxn_solution_map;}
	std::map<int, cxxExchange> & Get_Rxn_exchange_map() {return this->Rxn_exchange_map;}
	std::map<int, cxxGasPhase> & Get_Rxn_gas_phase_map() {return this->Rxn_gas_phase_map;}
	std::map<int, cxxKinetics> & Get_Rxn_kinetics_map() {return this->Rxn_kinetics_map;}
	std::map<int, cxxPPassemblage> & Get_Rxn_pp_assemblage_map() {return this->Rxn_pp_assemblage_map;}
	std::map<int, cxxSSassemblage> & Get_Rxn_ss_assemblage_map() {return this->Rxn_ss_assemblage_map;}
	std::map<int, cxxSurface> & Get_Rxn_surface_map() {return this->Rxn_surface_map;}
	std::map<int, cxxTemperature> & Get_Rxn_temperature_map() {return this->Rxn_temperature_map;}
	std::map<int, cxxPressure> & Get_Rxn_pressure_map() {return this->Rxn_pressure_map;}


protected:
	void init(void);

	//
	//Data members
	//
protected:
	PHRQ_io *phrq_io;
	PHRQ_io ioInstance;
	int same_model;

	LDBLE current_tc;
	LDBLE current_pa;
	LDBLE current_mu;
	bool mu_terms_in_logk;

	/* ----------------------------------------------------------------------
	*   STRUCTURES
	* ---------------------------------------------------------------------- */

	struct model last_model;
	//struct punch punch;
	bool high_precision;

	/* ----------------------------------------------------------------------
	*   Temperatures
	* ---------------------------------------------------------------------- */

	std::map<int, cxxTemperature> Rxn_temperature_map;

	/* ----------------------------------------------------------------------
	*   Pressures
	* ---------------------------------------------------------------------- */
	std::map<int, cxxPressure> Rxn_pressure_map;

	/* ----------------------------------------------------------------------
	*   Surface
	* --------------------------------------------------------------------- */

	int g_iterations;
	LDBLE G_TOL;
	std::map <int, cxxSurface> Rxn_surface_map;
	std::map <LDBLE, LDBLE> charge_group_map;
	int change_surf_count;
	struct Change_Surf *change_surf;

	/* ----------------------------------------------------------------------
	*   Exchange
	* ---------------------------------------------------------------------- */
	std::map<int, cxxExchange> Rxn_exchange_map;

	/* ----------------------------------------------------------------------
	*   Kinetics
	* ---------------------------------------------------------------------- */
	std::map<int, cxxKinetics> Rxn_kinetics_map;
	bool use_kinetics_limiter;

	/*----------------------------------------------------------------------
	*   Save
	*---------------------------------------------------------------------- */
	int count_save_values;
	struct save_values *save_values;
	struct save save;

	/*----------------------------------------------------------------------
	*   Use
	*---------------------------------------------------------------------- */
	cxxUse use;

	/*----------------------------------------------------------------------
	*   Copy
	*---------------------------------------------------------------------- */
	struct copier copy_solution;
	struct copier copy_pp_assemblage;
	struct copier copy_exchange;
	struct copier copy_surface;
	struct copier copy_ss_assemblage;
	struct copier copy_gas_phase;
	struct copier copy_kinetics;
	struct copier copy_mix;
	struct copier copy_reaction;
	struct copier copy_temperature;
	struct copier copy_pressure;

	/*----------------------------------------------------------------------
	*   Inverse
	*---------------------------------------------------------------------- */

	struct inverse *inverse;
	int count_inverse;

	/*----------------------------------------------------------------------
	*   Mix
	*---------------------------------------------------------------------- */
	std::map<int, cxxMix> Rxn_mix_map;
	std::map<int, cxxMix> Dispersion_mix_map;
	std::map<int, cxxMix> Rxn_solution_mix_map;
	std::map<int, cxxMix> Rxn_exchange_mix_map;
	std::map<int, cxxMix> Rxn_gas_phase_mix_map;
	std::map<int, cxxMix> Rxn_kinetics_mix_map;
	std::map<int, cxxMix> Rxn_pp_assemblage_mix_map;
	std::map<int, cxxMix> Rxn_ss_assemblage_mix_map;
	std::map<int, cxxMix> Rxn_surface_mix_map;
	/*
	* List new definitions
	*/
	std::set<int> Rxn_new_exchange;
	std::set<int> Rxn_new_gas_phase;
	std::set<int> Rxn_new_kinetics;     // not used
	std::set<int> Rxn_new_mix;          // not used
	std::set<int> Rxn_new_pp_assemblage;
	std::set<int> Rxn_new_pressure;     // not used
	std::set<int> Rxn_new_reaction;     // not used
	std::set<int> Rxn_new_solution;
	std::set<int> Rxn_new_ss_assemblage;
	std::set<int> Rxn_new_surface;
	std::set<int> Rxn_new_temperature;  // not used
	/*----------------------------------------------------------------------
	*   Irreversible reaction
	*---------------------------------------------------------------------- */
	std::map<int, cxxReaction> Rxn_reaction_map;

	/*----------------------------------------------------------------------
	*   Gas phase
	*---------------------------------------------------------------------- */
	std::map<int, cxxGasPhase> Rxn_gas_phase_map;

	/*----------------------------------------------------------------------
	*   Solid solution
	*---------------------------------------------------------------------- */
	std::map<int, cxxSSassemblage> Rxn_ss_assemblage_map;

	/*----------------------------------------------------------------------
	*   Pure-phase assemblage
	*---------------------------------------------------------------------- */
	std::map<int, cxxPPassemblage> Rxn_pp_assemblage_map;

	/*----------------------------------------------------------------------
	*   Species_list
	*---------------------------------------------------------------------- */
	int count_species_list;
	int max_species_list;
	struct species_list *species_list;

	/*----------------------------------------------------------------------
	*   Jacobian and Mass balance lists
	*---------------------------------------------------------------------- */

	int count_sum_jacob0;	/* number of elements in sum_jacob0 */
	int max_sum_jacob0;	/* calculated maximum number of elements in sum_jacob0 */
	struct list0 *sum_jacob0;	/* array of pointers to targets and coefficients for array */

	int count_sum_mb1;		/* number of elements in sum_mb1 */
	int max_sum_mb1;		/* calculated maximum number of elements in sum_mb1 */
	struct list1 *sum_mb1;	/* array of pointers to sources and targets for mass
							balance summations with coef = 1.0 */
	int count_sum_jacob1;	/* number of elements in sum_jacob1 */
	int max_sum_jacob1;	/* calculated maximum number of elements in sum_jacob1 */
	struct list1 *sum_jacob1;	/* array of pointers to sources and targets for array
								equations with coef = 1.0 */
	int count_sum_mb2;		/* number of elements in sum_mb2 */
	int max_sum_mb2;		/* calculated maximum number of elements in sum_mb2 */
	struct list2 *sum_mb2;	/* array of coefficients and pointers to sources and
							targets for mass balance summations with coef != 1.0 */
	int count_sum_jacob2;	/* number of elements in sum_jacob2 */
	int max_sum_jacob2;	/* calculated maximum number of elements in sum_jacob2 */
	struct list2 *sum_jacob2;	/* array of coefficients and pointers to sources and
								targets, coef != 1.0 */
	int count_sum_delta;	/* number of elements in sum_delta */
	int max_sum_delta;		/* calculated maximum number of elements in sum_delta */
	struct list2 *sum_delta;	/* array of pointers to sources, targets and coefficients for
								summing deltas for mass balance equations */
	/*----------------------------------------------------------------------
	*   Solution
	*---------------------------------------------------------------------- */
	std::map<int, cxxSolution> Rxn_solution_map;
	std::vector<cxxSolution> unnumbered_solutions;
	bool save_species;

	/*----------------------------------------------------------------------
	*   Global solution
	*---------------------------------------------------------------------- */
	char *title_x;
	int new_x;
	char *description_x;
	LDBLE tc_x;
	LDBLE tk_x;
	LDBLE patm_x;
	LDBLE last_patm_x;
	bool numerical_fixed_volume;
	bool force_numerical_fixed_volume;
	//bool switch_numerical;
	LDBLE ph_x;
	LDBLE solution_pe_x;
	LDBLE mu_x;
	LDBLE ah2o_x;
	LDBLE density_x;
	LDBLE total_h_x;
	LDBLE total_o_x;
	LDBLE cb_x;
	LDBLE total_ions_x;
	LDBLE mass_water_aq_x;
	LDBLE mass_water_surfaces_x;
	LDBLE mass_water_bulk_x;
	char *units_x;
	std::map < std::string, cxxChemRxn > pe_x;
	std::map<std::string, cxxSolutionIsotope> isotopes_x;
	std::string default_pe_x;
	cxxSurface::DIFFUSE_LAYER_TYPE dl_type_x;
	LDBLE total_carbon;
	LDBLE total_co2;
	LDBLE total_alkalinity;
	LDBLE gfw_water;
	LDBLE step_x;
	LDBLE kin_time_x;

	/*----------------------------------------------------------------------
	*   Transport data
	*---------------------------------------------------------------------- */
	int count_cells;
	int count_shifts;
	int ishift;
	int bcon_first;
	int bcon_last;
	int correct_disp;
	LDBLE tempr;
	LDBLE timest;
	int simul_tr;
	LDBLE diffc;
	LDBLE heat_diffc;
	int cell;
	LDBLE mcd_substeps;
	struct stag_data *stag_data;
	int print_modulus;
	int punch_modulus;
	int dump_in;
	int dump_modulus;
	int transport_warnings;
	struct cell_data *cell_data;
	int old_cells, max_cells, all_cells;
	int multi_Dflag;		/* signals calc'n of multicomponent diffusion */
	int interlayer_Dflag;	/* multicomponent diffusion and diffusion through interlayer porosity */
	LDBLE default_Dw;		/* default species diffusion coefficient in water at 25oC, m2/s */
	LDBLE multi_Dpor;		/* uniform porosity of free porewater in solid medium */
	LDBLE interlayer_Dpor;	/* uniform porosity of interlayer space of montmorillonite in solid medium */
	LDBLE multi_Dpor_lim;	/* limiting free porewater porosity where transport stops */
	LDBLE interlayer_Dpor_lim;	/* limiting interlayer porosity where transport stops */
	LDBLE multi_Dn;		/* exponent to calculate pore water diffusion coefficient,
						Dp = Dw * (multi_Dpor)^multi_Dn */
	LDBLE interlayer_tortf;	/* tortuosity_factor in interlayer porosity,
							Dpil = Dw / interlayer_tortf */

	int cell_no, mixrun;
	/*----------------------------------------------------------------------
	*   Advection data
	*---------------------------------------------------------------------- */
	int count_ad_cells;
	int count_ad_shifts;
	int print_ad_modulus;
	int punch_ad_modulus;
	int *advection_punch, *advection_print;
	LDBLE advection_kin_time;
	LDBLE advection_kin_time_defined;
	int advection_warnings;

	/*----------------------------------------------------------------------
	*   Tidy data
	*---------------------------------------------------------------------- */
	int new_model, new_exchange, new_pp_assemblage, new_surface,
		new_reaction, new_temperature, new_mix, new_solution, new_gas_phase,
		new_inverse, new_punch, new_ss_assemblage, new_kinetics, new_copy,
		new_pitzer;

	/*----------------------------------------------------------------------
	*   Elements
	*---------------------------------------------------------------------- */

	struct element **elements;
	int count_elements;
	int max_elements;
	struct element *element_h_one;

	/*----------------------------------------------------------------------
	*   Element List
	*---------------------------------------------------------------------- */

	struct elt_list *elt_list;	/* structure array of working space while reading equations
								names are in "strings", initially in input order */
	int count_elts;		/* number of elements in elt_list = position of next */
	int max_elts;
	/*----------------------------------------------------------------------
	*   Reaction
	*---------------------------------------------------------------------- */
	bool run_cells_one_step;
	/*----------------------------------------------------------------------
	*   Species
	*---------------------------------------------------------------------- */

	struct logk **logk;
	int count_logk;
	int max_logk;

	char *moles_per_kilogram_string;
	char *pe_string;

	struct species **s;
	int count_s;
	int max_s;
	std::vector< std::map < std::string, cxxSpeciesDL > > s_diff_layer;

	struct species **s_x;
	int count_s_x;
	int max_s_x;

	struct species *s_h2o;
	struct species *s_hplus;
	struct species *s_h3oplus;
	struct species *s_eminus;
	struct species *s_co3;
	struct species *s_h2;
	struct species *s_o2;

	/*----------------------------------------------------------------------
	*   Phases
	*---------------------------------------------------------------------- */
	struct phase **phases;
	int count_phases;
	int max_phases;

	/*----------------------------------------------------------------------
	*   Master species
	*---------------------------------------------------------------------- */
	struct master **master;	/* structure array of master species */
	struct master **dbg_master;
	int count_master;
	int max_master;

	/*----------------------------------------------------------------------
	*   Unknowns
	*---------------------------------------------------------------------- */

	struct unknown **x;
	int count_unknowns;
	int max_unknowns;

	struct unknown *ah2o_unknown;
	struct unknown *alkalinity_unknown;
	struct unknown *carbon_unknown;
	struct unknown *charge_balance_unknown;
	struct unknown *exchange_unknown;
	struct unknown *mass_hydrogen_unknown;
	struct unknown *mass_oxygen_unknown;
	struct unknown *mb_unknown;
	struct unknown *mu_unknown;
	struct unknown *pe_unknown;
	struct unknown *ph_unknown;
	struct unknown *pure_phase_unknown;
	struct unknown *solution_phase_boundary_unknown;
	struct unknown *surface_unknown;
	struct unknown *gas_unknown;
	struct unknown *ss_unknown;
	std::vector<struct unknown *> gas_unknowns;

	/*----------------------------------------------------------------------
	*   Reaction work space
	*---------------------------------------------------------------------- */
	struct reaction_temp trxn;	/* structure array of working space while reading equations
								species names are in "temp_strings" */
	int count_trxn;		/* number of reactants in trxn = position of next */
	int max_trxn;

	struct unknown_list *mb_unknowns;
	int count_mb_unknowns;
	int max_mb_unknowns;

	/* ----------------------------------------------------------------------
	*   Print
	* ---------------------------------------------------------------------- */
	struct prints pr;
	bool status_on;
	clock_t status_interval;
	clock_t status_timer;
	std::string status_string;
	int count_warnings;

	/* ----------------------------------------------------------------------
	*   RATES
	* ---------------------------------------------------------------------- */
	struct rate *rates;
	int count_rates;
	LDBLE rate_m, rate_m0, rate_time, rate_kin_time, rate_sim_time_start,
		rate_sim_time_end, rate_sim_time, rate_moles, initial_total_time;
	std::vector<LDBLE> rate_p;
	int count_rate_p;

	/* ----------------------------------------------------------------------
	*   USER PRINT COMMANDS
	* ---------------------------------------------------------------------- */
	struct rate *user_print;
	//struct rate *user_punch;
	//const char **user_punch_headings;
	//int user_punch_count_headings;
	int n_user_punch_index;

	int fpunchf_user_s_warning;
	char fpunchf_user_buffer[80];

#if defined PHREEQ98 
	struct rate *user_graph;
	char **user_graph_headings;
	int user_graph_count_headings;
#endif
#if defined MULTICHART
	ChartHandler chart_handler;
public:
	ChartHandler& Get_chart_handler(void)
	{
		return chart_handler;
	}
	const ChartHandler& Get_chart_handler(void)const
	{
		return chart_handler;
	}
protected:
#endif

	/* ----------------------------------------------------------------------
	*   GLOBAL DECLARATIONS
	* ---------------------------------------------------------------------- */
	const char * error_string;
	int simulation;
	int state;
	int reaction_step;
	int transport_step;
	int transport_start;
	int advection_step;
	int stop_program;
	int incremental_reactions;

	double MIN_LM;
	double LOG_ZERO_MOLALITY;
	double MIN_TOTAL;
	double MIN_TOTAL_SS;
	double MIN_RELATED_SURFACE;
	double MIN_RELATED_LOG_ACTIVITY;

	int count_strings;
	int max_strings;

	LDBLE *array;
	LDBLE *delta;
	LDBLE *residual;

	int input_error;

	Keywords::KEYWORDS next_keyword;
	int parse_error;
	int paren_count;
	int iterations;
	int gamma_iterations;
	int run_reactions_iterations;

	int max_line;
	char *line;
	char *line_save;

	LDBLE LOG_10;

	int debug_model;
	int debug_prep;
	int debug_set;
	int debug_diffuse_layer;
	int debug_inverse;

	LDBLE inv_tol_default;
	int itmax;
	int max_tries;
	LDBLE ineq_tol;
	LDBLE convergence_tolerance;
	LDBLE step_size;
	LDBLE pe_step_size;
	LDBLE step_size_now;
	LDBLE pe_step_size_now;
	LDBLE pp_scale;
	LDBLE pp_column_scale;
	int diagonal_scale;	/* 0 not used, 1 used */
	int mass_water_switch;
	int delay_mass_water;
	int equi_delay;
	bool dampen_ah2o;
	LDBLE censor;
	int aqueous_only;
	int negative_concentrations;
	int calculating_deriv;
	int numerical_deriv;

	int count_total_steps;
	int phast;
	LDBLE *llnl_temp, *llnl_adh, *llnl_bdh, *llnl_bdot, *llnl_co2_coefs;
	int llnl_count_temp, llnl_count_adh, llnl_count_bdh, llnl_count_bdot,
		llnl_count_co2_coefs;

	//char *selected_output_file_name;
	std::map<int, SelectedOutput> SelectedOutput_map;
	SelectedOutput * current_selected_output;
	
	std::map <int, UserPunch> UserPunch_map;
	UserPunch * current_user_punch;

	char *dump_file_name;
	int remove_unstable_phases;
	std::string screen_string;
#ifdef PHREEQCI_GUI
	struct spread_sheet g_spread_sheet;
#endif
	int spread_length;

	/* ---------------------------------------------------------------------- */
	/*
	*   Hash definitions
	*/

	std::map<std::string, std::string *> strings_map;
#ifdef HASH
	std::hash_map<std::string, std::string *> strings_hash;
#endif
	HashTable *elements_hash_table;
	HashTable *species_hash_table;
	HashTable *phases_hash_table;
	HashTable *logk_hash_table;
	HashTable *master_isotope_hash_table;

#if defined(PHREEQCI_GUI)
#include "../../phreeqci_gui.h"
#endif /* defined(PHREEQCI_GUI) */
	/* ----------------------------------------------------------------------
	*   ISOTOPES
	* ---------------------------------------------------------------------- */
	//struct name_coef match_tokens[50];
	//int count_match_tokens;
	int count_master_isotope;
	struct master_isotope **master_isotope;
	int max_master_isotope;
	int initial_solution_isotopes;
	int count_calculate_value;
	struct calculate_value **calculate_value;
	int max_calculate_value;
	HashTable *calculate_value_hash_table;
	int count_isotope_ratio;
	struct isotope_ratio **isotope_ratio;
	int max_isotope_ratio;
	HashTable *isotope_ratio_hash_table;
	int count_isotope_alpha;
	struct isotope_alpha **isotope_alpha;
	int max_isotope_alpha;
	HashTable *isotope_alpha_hash_table;
	int phreeqc_mpi_myself;
	int first_read_input;
	char *user_database;

	//int have_punch_name;
	/* VP: Density Start */
	int print_density;
	/* VP: Density End */

	LDBLE *zeros;
	int zeros_max;

	LDBLE cell_pore_volume;
	LDBLE cell_porosity;
	LDBLE cell_volume;
	LDBLE cell_saturation;
	struct system_species *sys;
	int count_sys, max_sys;
	LDBLE sys_tot;

	LDBLE V_solutes, rho_0, rho_0_sat, kappa_0, p_sat/*, ah2o_x0*/;
	LDBLE eps_r; // relative dielectric permittivity
	LDBLE DH_A, DH_B, DH_Av; // Debye-Hueckel A, B and Av
	LDBLE QBrn; // Born function d(ln(eps_r))/dP / eps_r * 41.84004, for supcrt calc'n of molal volume
	LDBLE ZBrn; // Born function (-1/eps_r + 1) * 41.84004, for supcrt calc'n of molal volume
	LDBLE dgdP; // dg / dP, pressure derivative of g-function, for supcrt calc'n of molal volume

	int need_temp_msg;
	LDBLE solution_mass, solution_volume;

	/* phqalloc.cpp ------------------------------- */
	PHRQMemHeader *s_pTail;

	/* Basic */
	PBasic * basic_interpreter;
	double (*basic_callback_ptr) (double x1, double x2, const char *str, void *cookie);
	void *basic_callback_cookie;
#ifdef IPHREEQC_NO_FORTRAN_MODULE
	double (*basic_fortran_callback_ptr) (double *x1, double *x2, char *str, size_t l);
#else
	double (*basic_fortran_callback_ptr) (double *x1, double *x2, const char *str, int l);
#endif
#if defined(SWIG) || defined(SWIG_IPHREEQC)
	class BasicCallback *basicCallback;
    void SetCallback(BasicCallback *cb) { basicCallback = cb; }
#endif

	/* cl1.cpp ------------------------------- */
	LDBLE *x_arg, *res_arg, *scratch;
	int x_arg_max, res_arg_max, scratch_max;
#ifdef SKIP
	/* dw.cpp ------------------------------- */
	/* COMMON /QQQQ/ */
	LDBLE Q0, Q5;
	LDBLE GASCON, TZ, AA;
	LDBLE Z, DZ, Y;
	LDBLE G1, G2, GF;
	LDBLE B1, B2, B1T, B2T, B1TT, B2TT;
#endif
	/* gases.cpp ------------------------------- */
	LDBLE a_aa_sum, b2, b_sum, R_TK;

	/* input.cpp ------------------------------- */
	int check_line_return;  
	int reading_db;

	/* integrate.cpp ------------------------------- */
	LDBLE midpoint_sv;
	LDBLE z_global, xd_global, alpha_global;

	/* inverse.cpp ------------------------------- */
	int max_row_count, max_column_count;
	int carbon;
	const char **col_name, **row_name;
	int count_rows, count_optimize;
	int col_phases, col_redox, col_epsilon, col_ph, col_water,
		col_isotopes, col_phase_isotopes;
	int row_mb, row_fract, row_charge, row_carbon, row_isotopes,
		row_epsilon, row_isotope_epsilon, row_water;
	LDBLE *inv_zero, *array1, *inv_res, *inv_delta1, *delta2, *delta3, *inv_cu,
		*delta_save;
	LDBLE *min_delta, *max_delta;
	int *inv_iu, *inv_is;
	int klmd, nklmd, n2d, kode, iter;
	LDBLE toler, error, max_pct, scaled_error;
	struct master *master_alk;
	int *row_back, *col_back;
	unsigned long *good, *bad, *minimal;
	int max_good, max_bad, max_minimal;
	int count_good, count_bad, count_minimal, count_calls;
	unsigned long soln_bits, phase_bits, current_bits, temp_bits;
	FILE *netpath_file;
	int count_inverse_models, count_pat_solutions;
	int min_position[32], max_position[32], now[32];
	std::vector <std::string> inverse_heading_names;

	/* kinetics.cpp ------------------------------- */
public:
	int count_pp, count_pg, count_ss;
	void *cvode_kinetics_ptr;
	int cvode_test;
	int cvode_error;
	int cvode_n_user;
	int cvode_n_reactions;
	realtype cvode_step_fraction;
	realtype cvode_rate_sim_time;
	realtype cvode_rate_sim_time_start;
	realtype cvode_last_good_time;
	realtype cvode_prev_good_time;
	N_Vector cvode_last_good_y;
	N_Vector cvode_prev_good_y;
	M_Env kinetics_machEnv;
	N_Vector kinetics_y, kinetics_abstol;
	void *kinetics_cvode_mem;
	cxxSSassemblage *cvode_ss_assemblage_save;
	cxxPPassemblage *cvode_pp_assemblage_save;
protected:
	LDBLE *m_original;
	LDBLE *m_temp;
	LDBLE *rk_moles;
	int set_and_run_attempt;
	LDBLE *x0_moles;

	/* model.cpp ------------------------------- */
	int gas_in;
	LDBLE min_value;
	LDBLE *normal, *ineq_array, *res, *cu, *zero, *delta1;
	int *iu, *is, *back_eq;
	int normal_max, ineq_array_max, res_max, cu_max, zero_max,
		delta1_max, iu_max, is_max, back_eq_max;

	/* phrq_io_output.cpp ------------------------------- */
	int forward_output_to_log;

	/* phreeqc_files.cpp ------------------------------- */
	char *default_data_base;
#ifdef PHREEQ98
	int outputlinenr;
	char *LogFileNameC;
	char progress_str[512];
#endif

	/* Pitzer  */
	int pitzer_model, sit_model, pitzer_pe;
	int full_pitzer, always_full_pitzer, ICON, IC;
	LDBLE COSMOT;
	LDBLE AW;
	LDBLE VP, DW0;
	struct pitz_param **pitz_params;
	int count_pitz_param, max_pitz_param;
	std::map< std::string, size_t > pitz_param_map;
	struct theta_param **theta_params;
	int count_theta_param, max_theta_param;
	int use_etheta;
	LDBLE OTEMP, OPRESS;
	LDBLE A0;
	struct species **spec, **cations, **anions, **neutrals;
	int count_cations, count_anions, count_neutrals;
	int MAXCATIONS, FIRSTANION, MAXNEUTRAL;
	struct pitz_param *mcb0, *mcb1, *mcc0;
	int *IPRSNT;
	LDBLE *M, *LGAMMA;
	LDBLE BK[23], DK[23];

#ifdef PHREEQ98
	int connect_simulations, graph_initial_solutions;
	int shifts_as_points;
	int chart_type;
	int ShowChart;
	int RowOffset, ColumnOffset;
#endif
	LDBLE dummy;

	/* print.cpp ------------------------------- */
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
	char *prev_next_char;
#if defined PHREEQ98 
	int shifts_as_points;
#endif

	/* read_class.cxx */
	dumper dump_info;
	StorageBinList delete_info;
	runner run_info;
	char * sformatf_buffer;
	size_t sformatf_buffer_size;

	/* readtr.cpp */
	std::string dump_file_name_cpp;

	/* sit.cpp ------------------------------- */
	struct pitz_param **sit_params;
	int count_sit_param, max_sit_param;
	std::map< std::string, size_t > sit_param_map;
	LDBLE sit_A0;
	int sit_count_cations, sit_count_anions, sit_count_neutrals;
	int sit_MAXCATIONS, sit_FIRSTANION, sit_MAXNEUTRAL;
	int *sit_IPRSNT;
	LDBLE *sit_M, *sit_LGAMMA;
	std::vector<int> s_list, cation_list, neutral_list, anion_list, ion_list, param_list;

	/* tidy.cpp ------------------------------- */
	LDBLE a0, a1, kc, kb;

	/* tally.cpp ------------------------------- */
	struct tally_buffer *t_buffer;
	int tally_count_component;
	struct tally *tally_table;
	int count_tally_table_columns;
	int count_tally_table_rows;

	/* transport.cpp ------------------------------- */
	struct sol_D *sol_D;
	struct sol_D *sol_D_dbg;
	struct J_ij *J_ij, *J_ij_il;
	int J_ij_count_spec;

	struct M_S *m_s;
	int count_m_s;
	LDBLE tot1_h, tot1_o, tot2_h, tot2_o;
	LDBLE diffc_max, diffc_tr, J_ij_sum;
	int transp_surf;
	LDBLE *heat_mix_array;
	LDBLE *temp1, *temp2;
	int nmix, heat_nmix;
	LDBLE heat_mix_f_imm, heat_mix_f_m;
	int warn_MCD_X, warn_fixed_Surf;

#ifdef PHREEQ98
	int AutoLoadOutputFile, CreateToC;
	int ProcessMessages, ShowProgress, ShowProgressWindow, ShowChart;
	int outputlinenr;
	int stop_calculations;
	char err_str98[80];
#endif
	/* utilities.cpp ------------------------------- */
	int spinner;
	std::map<std::string, double> gfw_map;
	std::map<const char *, int> rates_map;

	/* new after release of Version 3 */
	std::map<std::string, std::vector < std::string> > sum_species_map; 
	std::map<std::string, std::vector < std::string> > sum_species_map_db;

	friend class PBasic;
	friend class ChartObject;
	friend class IPhreeqc;
	friend class TestIPhreeqc;
	friend class TestSelectedOutput;
	friend class IPhreeqcMMS;
	friend class IPhreeqcPhast;
	friend class PhreeqcRM;

	std::vector<int> keycount;  // used to mark keywords that have been read 

public:
	static const struct const_iso iso_defaults[];
	static const int count_iso_defaults;
};
#endif /* _INC_PHREEQC_H */

#ifndef _INC_ISFINITE_H
#define _INC_ISFINITE_H
	/*********************************
	isfinite handling
	(Note: Should NOT be guarded)
	**********************************/

#if defined (PHREEQ98) || defined (_MSC_VER)
#  define HAVE_FINITE
#  define finite _finite
#else  /*defined (PHREEQ98) || defined (_MSC_VER)*/
#  if defined(DJGPP)
#    define HAVE_FINITE
#  endif
#endif /*defined (PHREEQ98) || defined (_MSC_VER)*/

#if defined(HAVE_ISFINITE)
#  define PHR_ISFINITE(x) isfinite(x)
#elif defined(HAVE_FINITE)
#  define PHR_ISFINITE(x) finite(x)
#elif defined(HAVE_ISNAN)
#  define PHR_ISFINITE(x) ( ((x) == 0.0) || ((!isnan(x)) && ((x) != (2.0 * (x)))) )
#else
#  define PHR_ISFINITE(x) ( ((x) == 0.0) || (((x) == (x)) && ((x) != (2.0 * (x)))) )
#endif
#endif // _INC_ISFINITE_H

#ifndef _INC_UTILITIES_NAMESPACE_H
#define _INC_UTILITIES_NAMESPACE_H
namespace Utilities
{
	LDBLE get_nan(void);

	// operations on maps of entities (Solution, Exchange, ...)
	template < typename T >
	void Rxn_dump_raw(const T & b, std::ostream & s_oss, unsigned int indent)
	{
		typename T::const_iterator it;
		for (it = b.begin(); it != b.end(); ++it)
		{
			// Adding logic to dump only non-negative entities
			if (it->second.Get_n_user() >= 0)
			{
				it->second.dump_raw(s_oss, indent);
			}
		}
		return;
	}

	template < typename T >
	T * Rxn_find(std::map < int, T > &b, int i)
	{
		if (b.find(i) != b.end())
		{
			return (&(b.find(i)->second));
		}
		else
		{
			return (NULL);
		}
	}

	template < typename T >
	int Rxn_next_user_number(std::map < int, T > &b)
	{
		int ret = 0;
		if (b.size() != 0)
		{
			ret = b.rbegin()->first + 1;
		}
		return ret;
	}

	template < typename T >
	T * Rxn_copy(std::map < int, T > &b, int i, int j)
	{
		typename std::map < int, T >::iterator it;
		it = b.find(i);
		if (it != b.end())
		{
			b[j] = it->second;
			it = b.find(j);
			it->second.Set_n_user(j);
			it->second.Set_n_user_end(j);
			return &(it->second);
		}
		else
		{
			return (NULL);
		}
	}

	template < typename T >
	void Rxn_copies(std::map < int, T > &b, int n_user, int n_user_end)
	{
		if (n_user_end <= n_user) return;
		typename std::map < int, T >::iterator it;
		it = b.find(n_user);
		if (it != b.end())
		{
			for (int j = n_user + 1; j <= n_user_end; j++)
			{
				b[j] = it->second;
				it = b.find(j);
				it->second.Set_n_user(j);
				it->second.Set_n_user_end(j);
			}
		}
	}
	template < typename T >
	int Rxn_read_raw(std::map < int, T > &m, std::set < int > &s, Phreeqc * phreeqc_cookie)
	{
		typename std::map < int, T >::iterator it;
		assert(!phreeqc_cookie->reading_database());

		T entity(phreeqc_cookie->Get_phrq_io());

		CParser parser(phreeqc_cookie->Get_phrq_io());
		entity.read_raw(parser);

		// Store
		if (entity.Get_base_error_count() == 0)
		{
			m[entity.Get_n_user()] = entity;
		}

		// Make copies if necessary
		Utilities::Rxn_copies(m, entity.Get_n_user(), entity.Get_n_user_end());
		for (int i = entity.Get_n_user(); i <= entity.Get_n_user_end(); i++)
		{
			s.insert(i);
		}
		return phreeqc_cookie->cleanup_after_parser(parser);
	}

#ifdef SKIP
	template < typename T >
	int Rxn_read_modify(std::map < int, T > &m, std::set < int > &s, Phreeqc * phreeqc_cookie)
	{
		typename std::map < int, T >::iterator it;
		
		CParser parser(phreeqc_cookie->Get_phrq_io());

		std::string key_name;
		std::string::iterator b = parser.line().begin();
		std::string::iterator e = parser.line().end();
		CParser::copy_token(key_name, b, e);

		cxxNumKeyword nk;
		nk.read_number_description(parser);
		T * entity_ptr = Utilities::Rxn_find(m, nk.Get_n_user());
		if (!entity_ptr)
		{
			std::ostringstream errstr;
			errstr <<  "Could not find " << key_name << " " << nk.Get_n_user() << " to modify.\n";
			phreeqc_cookie->error_msg(errstr.str().c_str(), PHRQ_io::OT_STOP);
		}

		entity_ptr->read_raw(parser, false);
		entity_ptr->Set_n_user(nk.Get_n_user());
		entity_ptr->Set_n_user_end(nk.Get_n_user_end());
		entity_ptr->Set_description(nk.Get_description());
		s.insert(entity_ptr->Get_n_user());

		return phreeqc_cookie->cleanup_after_parser(parser);
	}
#endif

	template < typename T >
	int Rxn_read_modify(std::map < int, T > &m, std::set < int > &s, Phreeqc * phreeqc_cookie)
	{
		typename std::map < int, T >::iterator it;
		
		CParser parser(phreeqc_cookie->Get_phrq_io());

		std::string key_name;
		std::string::iterator b = parser.line().begin();
		std::string::iterator e = parser.line().end();
		CParser::copy_token(key_name, b, e);

		cxxNumKeyword nk;
		nk.read_number_description(parser);
		T * entity_ptr = Utilities::Rxn_find(m, nk.Get_n_user());
		if (!entity_ptr)
		{
			std::ostringstream errstr;
			errstr <<  "Could not find " << key_name << " " << nk.Get_n_user() << ", ignoring modify data.\n";
			phreeqc_cookie->warning_msg(errstr.str().c_str());
			//phreeqc_cookie->error_msg(errstr.str().c_str(), PHRQ_io::OT_STOP);

			// Don't throw, read data into dummy entity, then ignore
			T entity;
			entity_ptr = &entity;
			entity_ptr->read_raw(parser, false);
			return phreeqc_cookie->cleanup_after_parser(parser);
		}

		entity_ptr->read_raw(parser, false);
		entity_ptr->Set_n_user(nk.Get_n_user());
		entity_ptr->Set_n_user_end(nk.Get_n_user_end());
		entity_ptr->Set_description(nk.Get_description());
		s.insert(entity_ptr->Get_n_user());

		return phreeqc_cookie->cleanup_after_parser(parser);
	}

	template < typename T >
	void Rxn_mix(std::map <int, cxxMix> &mix_map, std::map < int, T > &entity_map, Phreeqc * phreeqc_cookie)
	{
		std::map<int, cxxMix>::iterator mix_it;
		for (mix_it = mix_map.begin(); mix_it != mix_map.end(); mix_it++)
		{
			T entity(entity_map, mix_it->second, mix_it->second.Get_n_user(), phreeqc_cookie->Get_phrq_io());
			entity_map[mix_it->second.Get_n_user()] = entity;
			Utilities::Rxn_copies(entity_map, mix_it->second.Get_n_user(), mix_it->second.Get_n_user_end());
		}
		mix_map.clear();
	}

} // namespace Utilities


#if defined(PHREEQCI_GUI)
void PhreeqcIWait(Phreeqc *phreeqc);
#endif

#if !defined(NDEBUG) && defined(WIN32_MEMORY_DEBUG)
#define   string_duplicate(s)             _string_duplicate(s, __FILE__, __LINE__)
#endif
#if defined(_DEBUG)
	char * _string_duplicate(const char *token, const char *szFileName, int nLine);
#endif

#endif //_INC_UTILITIES_NAMESPACE_H