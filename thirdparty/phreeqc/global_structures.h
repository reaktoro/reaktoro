#ifndef _INC_GLOBAL_STRUCTURES_H
#define _INC_GLOBAL_STRUCTURES_H
#include "Surface.h"
/* ----------------------------------------------------------------------
 *   #define DEFINITIONS
 * ---------------------------------------------------------------------- */
#ifndef NAN
#   define NAN -99999999
#endif
#define MISSING -9999.999            
#include "NA.h"   /* NA = not available */

#define F_C_MOL 96493.5			/* C/mol or joule/volt-eq */
#define F_KJ_V_EQ  96.4935		/* kJ/volt-eq */
#define F_KCAL_V_EQ 23.0623		/* kcal/volt-eq */
#define R_LITER_ATM 0.0820597	/* L-atm/deg-mol */
#define R_KCAL_DEG_MOL 0.00198726	/* kcal/deg-mol */
#define R_KJ_DEG_MOL 0.00831470	/* kJ/deg-mol */
#define EPSILON 78.5			/* dialectric constant, dimensionless. Is calculated as eps_r(P, T) in calc_dielectrics. Update the code?? */
#define EPSILON_ZERO 8.854e-12	/* permittivity of free space, C/V-m = C**2/m-J */
#define JOULES_PER_CALORIE 4.1840
#define PASCAL_PER_ATM 1.01325E5 /* conversion from atm to Pa */
#define AVOGADRO 6.02252e23		/* atoms / mole */
#define pi 3.14159265358979
#define AH2O_FACTOR 0.017

#define TRUE 1
#define FALSE 0
#define OK 1
#define ERROR 0
#define STOP 1
#define CONTINUE 0

#define OPTION_EOF -1
#define OPTION_KEYWORD -2
#define OPTION_ERROR -3
#define OPTION_DEFAULT -4
#define OPT_1 -5

#define DISP 2
#define STAG 3
#define NOMIX 4

#define CONVERGED 2
#define MASS_BALANCE 3

#define REWRITE 2
#define INIT -1

/* check_line values, plus EMPTY, EOF, OK */
#define KEYWORD 3

/* copy_token values */
#define EMPTY 2
#define UPPER 4
#define LOWER 5
#define DIGIT 6
#define UNKNOWN 7
#define OPTION 8

/* species types */
#define AQ 0
#define HPLUS 1
#define H2O 2
#define EMINUS 3
#define SOLID 4
#define EX 5
#define SURF 6
#define SURF_PSI 7
#define SURF_PSI1 8
#define SURF_PSI2 9

/* unknown types */
#define MB 10
#define ALK 11
#define CB 12
#define SOLUTION_PHASE_BOUNDARY 13
#define MU 14
#define AH2O 15
#define MH 16
#define MH2O 17
#define PP 18
#define EXCH 19
#define SURFACE 20
#define SURFACE_CB 21
#define SURFACE_CB1 22
#define SURFACE_CB2 23
#define GAS_MOLES 24
#define SS_MOLES 25
#define PITZER_GAMMA 26
#define SLACK 28
/* state */
#define INITIALIZE	       0
#define INITIAL_SOLUTION   1
#define INITIAL_EXCHANGE   2
#define INITIAL_SURFACE 3
#define INITIAL_GAS_PHASE  4
#define REACTION		   5
#define INVERSE		 6
#define ADVECTION		 7
#define TRANSPORT		 8
#define PHAST		     9

/* constaints in mass balance */
#define EITHER 0
#define DISSOLVE 1
#define PRECIPITATE -1

/* gas phase type */
#define PRESSURE 1
#define VOLUME 2

#define MAX_PP_ASSEMBLAGE 10	/* default estimate of the number of phase assemblages */
#define MAX_ADD_EQUATIONS 20	/* maximum number of equations added together to reduce eqn to
								   master species */
#define MAX_ELEMENTS 50			/* default estimate of the number of elements */
#define MAX_LENGTH 256			/* maximum number of characters component name */
#define MAX_LINE 4096			/* estimate of maximum line length */
#define MAX_MASS_BALANCE 10		/* initial guess of number mass balance equations for a solution */
#define MAX_MASTER 50			/* default estimate of the number of master species */
#define MAX_ELTS 15				/* default estimate for maximum number of times elements occur in
								   an equation */
#define MAX_PHASES 500			/* initial guess of number of phases defined */
#define MAX_S 500				/* default estimate for maximum number of species in aqueous model */
#define MAX_SUM_JACOB0 50		/* list used to calculate jacobian */
#define MAX_SUM_JACOB1 500		/* list used to calculate jacobian */
#define MAX_SUM_JACOB2 500		/* list used to calculate jacobian */
#define MAX_SUM_MB 500			/* list used to calculate mass balance sums */
#define MAX_TRXN 16				/* default estimate for maximum number of components in an eqn */
#define MAX_UNKNOWNS 15			/* default estimate for maximum number of unknowns in model */
#define TOL 1e-9				/* tolerance for comparisons of double numbers */
#define MAX_LM 3.0				/* maximum log molality allowed in intermediate iterations */
#define MAX_M 1000.0
#ifdef USE_DECIMAL128
//#define MIN_LM -80.0			/* minimum log molality allowed before molality set to zero */
//#define LOG_ZERO_MOLALITY -80	/* molalities <= LOG_ZERO_MOLALITY are considered equal to zero */
//#define MIN_TOTAL 1e-60
//#define MIN_TOTAL_SS MIN_TOTAL/100
//#define MIN_RELATED_SURFACE MIN_TOTAL*100
//#define MIN_RELATED_LOG_ACTIVITY -60
#else
//#define MIN_LM -30.0			/* minimum log molality allowed before molality set to zero */
//#define LOG_ZERO_MOLALITY -30	/* molalities <= LOG_ZERO_MOLALITY are considered equal to zero */
//#define MIN_TOTAL 1e-25
//#define MIN_TOTAL_SS MIN_TOTAL/100
//#define MIN_RELATED_SURFACE MIN_TOTAL*100
//#define MIN_RELATED_LOG_ACTIVITY -30
#endif
#define REF_PRES_PASCAL 1.01325E5   /* Reference pressure: 1 atm */
/*
 *   Hash definitions
 */
# define SegmentSize		    256
# define SegmentSizeShift	  8	/* log2(SegmentSize) */
# define DirectorySize	    256
# define DirectorySizeShift      8	/* log2(DirectorySize)  */
# define Prime1			  37
# define Prime2			  1048583
# define DefaultMaxLoadFactor   5
//
// Typedefs and structure definitions
//
typedef enum { kcal, cal, kjoules, joules } DELTA_H_UNIT;
typedef enum { cm3_per_mol, dm3_per_mol, m3_per_mol } DELTA_V_UNIT;
enum entity_type
{ Solution, Reaction, Exchange, Surface, Gas_phase, Pure_phase, Ss_phase,
	Kinetics, Mix, Temperature, Pressure, UnKnown
};

typedef enum {
	logK_T0,
	delta_h,
	T_A1,
	T_A2,
	T_A3,
	T_A4,
	T_A5,
	T_A6,
	delta_v,	/* set in calc_delta_v: calculated molar volume-change of the reaction */
	vm_tc,		/* set in calc_vm: calculated molal volume of the species at tc */
	vm0,		/* read: molar volume of a phase */
	vma1, vma2, vma3, vma4, /* read: a1..a4 from supcrt, see calc_vm */
	wref,       /* from supcrt */
	b_Av,		/* b in z^2 * A_v * log(1 + b * I^0.5) / (2 * b) */
	vmi1, vmi2, vmi3, vmi4, /* ionic strength terms: (i1 + i2/(TK - 228) + i3 * (TK - 228) ) * I^i4 */
	MAX_LOG_K_INDICES	/* Keep this definition at the end of the enum */
} LOG_K_INDICES;
/* HSEARCH(3C) */
typedef struct entry
{
	const char *key;
	void *data;
} ENTRY;
typedef enum
{ FIND, ENTER } ACTION;

/* TSEARCH(3C) */
typedef enum
{ preorder, postorder, endorder, leaf } VISIT;

typedef struct Element
{
	/*
	 ** The user only sees the first two fields,
	 ** as we pretend to pass back only a pointer to ENTRY.
	 ** {S}he doesn`t know what else is in here.
	 */
	const char *Key;
	char *Data;
	struct Element *Next;		/* secret from user    */
} Element, *Segment;

typedef struct
{
	short p;					/* Next bucket to be split      */
	short maxp;					/* upper bound on p during expansion */
	long KeyCount;				/* current # keys       */
	short SegmentCount;			/* current # segments   */
	short MinLoadFactor;
	short MaxLoadFactor;
	Segment *Directory[DirectorySize];
} HashTable;

typedef unsigned long Address;

typedef struct PHRQMemHeader
{
	struct PHRQMemHeader *pNext;	/* memory allocated just after this one */
	struct PHRQMemHeader *pPrev;	/* memory allocated just prior to this one */
	size_t size;				/* memory request + sizeof(PHRQMemHeader) */
#if !defined(NDEBUG)
	char *szFileName;			/* file name */
	int nLine;					/* line number */
	int dummy;					/* alignment */
#endif
} PHRQMemHeader;

struct model
{
	int force_prep;
	LDBLE temperature;
	int count_exchange;
	struct master **exchange;

	int count_kinetics;
	struct kinetics *kinetics;

	int count_gas_phase;
	struct phase **gas_phase;

	int count_ss_assemblage;
	const char **ss_assemblage;

	int count_pp_assemblage;
	struct phase **pp_assemblage;
	const char **add_formula;
	LDBLE *si;

	cxxSurface::DIFFUSE_LAYER_TYPE dl_type;
	cxxSurface::SURFACE_TYPE surface_type;
	int only_counter_ions;

	LDBLE thickness;
	int count_surface_comp;
	const char **surface_comp;
	int count_surface_charge;
	const char **surface_charge;
	LDBLE pressure;
	bool numerical_fixed_volume;
};



struct name_master
{
	const char *name;
	struct master *master;
};
struct name_species
{
	const char *name;
	struct species *s;
};
struct name_phase
{
	const char *name;
	struct phase *phase;
};
#ifdef SKIP
struct punch
{
	int in;
	int new_def;
	struct name_master *totals;
	int count_totals;
	struct name_species *molalities;
	int count_molalities;
	struct name_species *activities;
	int count_activities;
	struct name_phase *pure_phases;
	int count_pure_phases;
	struct name_phase *si;
	int count_si;
	struct name_phase *gases;
	int count_gases;
	struct name_phase *s_s;
	int count_s_s;
	struct name_phase *kinetics;
	int count_kinetics;
	struct name_master *isotopes;
	int count_isotopes;
	struct name_master *calculate_values;
	int count_calculate_values;
	int inverse;
	int sim;
	int state;
	int soln;
	int dist;
	int time;
	int step;
	int rxn;
	int temp;
	int ph;
	int pe;
	int alk;
	int mu;
	int water;
	int high_precision;
	int user_punch;
	int charge_balance;
	int percent_error;
};
#endif
struct Change_Surf
{
	const char *comp_name;
	LDBLE fraction;
	const char *new_comp_name;
	LDBLE new_Dw;
	int cell_no;
	int next;
};

struct Charge_Group
{
	LDBLE z;
	LDBLE eq;
};

/*----------------------------------------------------------------------
 *   Save
 *---------------------------------------------------------------------- */
struct save_values
{
	LDBLE value;
	int count_subscripts;
	int *subscripts;
};

struct save
{
	int solution;
	int n_solution_user;
	int n_solution_user_end;
	int mix;
	int n_mix_user;
	int n_mix_user_end;
	int reaction;
	int n_reaction_user;
	int n_reaction_user_end;
	int pp_assemblage;
	int n_pp_assemblage_user;
	int n_pp_assemblage_user_end;
	int exchange;
	int n_exchange_user;
	int n_exchange_user_end;
	int kinetics;
	int n_kinetics_user;
	int n_kinetics_user_end;
	int surface;
	int n_surface_user;
	int n_surface_user_end;
	int gas_phase;
	int n_gas_phase_user;
	int n_gas_phase_user_end;
	int ss_assemblage;
	int n_ss_assemblage_user;
	int n_ss_assemblage_user_end;
};

/*----------------------------------------------------------------------
 *   Copy
 *---------------------------------------------------------------------- */
struct copier
{
	int count;
	int max;
	int *n_user;
	int *start;
	int *end;
};

/*----------------------------------------------------------------------
 *   Inverse
 *---------------------------------------------------------------------- */
struct inverse
{
	int n_user;
	char *description;
	int new_def;
	int minimal;
	int range;
	int mp;
	LDBLE mp_censor;
	LDBLE range_max;
	LDBLE tolerance;
	LDBLE mp_tolerance;
	int count_uncertainties;
	LDBLE *uncertainties;
	int count_ph_uncertainties;
	LDBLE *ph_uncertainties;
	LDBLE water_uncertainty;
	int mineral_water;
	int carbon;
	LDBLE *dalk_dph;
	LDBLE *dalk_dc;
	int count_solns;
	int *solns;
	int count_force_solns;
	int *force_solns;
	int count_elts;
	struct inv_elts *elts;
	int count_phases;
	struct inv_phases *phases;
	int count_master_list;
	struct master **master_list;
	int count_redox_rxns;
	int count_isotopes;
	struct inv_isotope *isotopes;
	int count_i_u;
	struct inv_isotope *i_u;
	int count_isotope_unknowns;
	struct isotope *isotope_unknowns;
	const char *netpath;
	const char *pat;
};
struct inv_elts
{
	const char *name;
	struct master *master;
	int row;
	int count_uncertainties;
	LDBLE *uncertainties;
};
struct inv_isotope
{
	const char *isotope_name;
	LDBLE isotope_number;
	const char *elt_name;
	int count_uncertainties;
	LDBLE *uncertainties;
};
struct inv_phases
{
	const char *name;
	struct phase *phase;
	int column;
	int constraint;
	int force;
	int count_isotopes;
	struct isotope *isotopes;
};
struct name_coef
{
	const char *name;
	LDBLE coef;
};
/*----------------------------------------------------------------------
 *   Species_list
 *---------------------------------------------------------------------- */
struct species_list
{
	struct species *master_s;
	struct species *s;
	LDBLE coef;
};

/*----------------------------------------------------------------------
 *   Jacobian and Mass balance lists
 *---------------------------------------------------------------------- */
struct list0
{
	LDBLE *target;
	LDBLE coef;
};
struct list1
{
	LDBLE *source;
	LDBLE *target;
};
struct list2
{
	LDBLE *source;
	LDBLE *target;
	LDBLE coef;
};

 struct isotope
 {
 	LDBLE isotope_number;
 	const char *elt_name;
 	const char *isotope_name;
 	LDBLE total;
 	LDBLE ratio;
 	LDBLE ratio_uncertainty;
 	LDBLE x_ratio_uncertainty;
 	struct master *master;
 	struct master *primary;
 	LDBLE coef;					/* coefficient of element in phase */
};
struct iso
{
	const char *name;
	LDBLE value;
	LDBLE uncertainty;
};
/*----------------------------------------------------------------------
 *   Transport data
 *---------------------------------------------------------------------- */
struct stag_data
{
	int count_stag;
	LDBLE exch_f;
	LDBLE th_m;
	LDBLE th_im;
};
struct cell_data
{
	LDBLE length;
	LDBLE mid_cell_x;
	LDBLE disp;
	LDBLE temp;
	LDBLE por;					/* free (uncharged) porewater porosities */
	LDBLE por_il;				/* interlayer water porosities */
	int punch;
	int print;
};

/*----------------------------------------------------------------------
 *   Keywords
 *---------------------------------------------------------------------- */
 struct key
 {
 	char *name;
 	int keycount;
 };
 struct const_key
 {
 	const char *name;
 	int keycount;
};

/*----------------------------------------------------------------------
 *   Elements
 *---------------------------------------------------------------------- */
struct element
{
	const char *name;					/* element name */
	/*    int in; */
	struct master *master;
	struct master *primary;
	LDBLE gfw;
};
/*----------------------------------------------------------------------
 *   Element List
 *---------------------------------------------------------------------- */
struct elt_list
{								/* list of name and number of elements in an equation */
	struct element *elt;		/* pointer to element structure */
	LDBLE coef;					/* number of element e's in eqn */
};
/*----------------------------------------------------------------------
 *   Reaction
 *---------------------------------------------------------------------- */
struct reaction
{
	LDBLE logk[MAX_LOG_K_INDICES];
	LDBLE dz[3];
	struct rxn_token *token;
};
struct rxn_token
{
	struct species *s;
	LDBLE coef;
	const char *name;
};
class cxxChemRxn
{
public:
	cxxChemRxn(void)
	{
		for (size_t i = 0; i < MAX_LOG_K_INDICES; i++)
		{
			logk[i] = 0.0;
		}
		for (size_t i = 0; i < 3; i++)
		{
			dz[i] =0.0;
		}
	}
	cxxChemRxn(struct reaction *rxn)
	{
		for (size_t i = 0; i < MAX_LOG_K_INDICES; i++)
		{
			logk[i] = rxn->logk[i];
		}
		for (size_t i = 0; i < 3; i++)
		{
			dz[i] = rxn->dz[i];
		}
		struct rxn_token *next_token;
		next_token = rxn->token;
		this->tokens.push_back(*next_token++);
		while (next_token->s != NULL || next_token->name != NULL)
		{
			this->tokens.push_back(*next_token++);
		}
	}
	~cxxChemRxn(void) {}
	LDBLE *Get_logk(void) {return this->logk;}
	void   Set_logk(LDBLE *d)
	{
		for (size_t i = 0; i < MAX_LOG_K_INDICES; i++)
		{
			logk[i] = d[i];
		}

	}
	LDBLE *Get_dz(void) {return this->dz;}
	void   Set_dz(LDBLE *d)
	{
		for (size_t i = 0; i < 3; i++)
		{
			dz[i] = d[i];
		}
	}
	std::vector<struct rxn_token> &Get_tokens(void) {return this->tokens;}
	void Set_tokens(const std::vector<struct rxn_token> &t) {this->tokens = t;}

protected:
	LDBLE logk[MAX_LOG_K_INDICES];
	LDBLE dz[3];
	std::vector<struct rxn_token> tokens;
};
/*----------------------------------------------------------------------
 *   Species
 *---------------------------------------------------------------------- */
struct species
{								/* all data pertinent to an aqueous species */
	const char *name;					/* name of species */
	const char *mole_balance;			/* formula for mole balance */
	int in;						/* species used in model if TRUE */
	int number;
	struct master *primary;		/* points to master species list, NULL if not primary master */
	struct master *secondary;	/* points to master species list, NULL if not secondary master */
	LDBLE gfw;					/* gram formula wt of species */
	LDBLE z;					/* charge of species */
	LDBLE dw;					/* tracer diffusion coefficient in water at 25oC, m2/s */
	LDBLE erm_ddl;				/* enrichment factor in DDL */
	LDBLE equiv;				/* equivalents in exchange species */
	LDBLE alk;					/* alkalinity of species, used for cec in exchange */
	LDBLE carbon;				/* stoichiometric coefficient of carbon in species */
	LDBLE co2;					/* stoichiometric coefficient of C(4) in species */
	LDBLE h;					/* stoichiometric coefficient of H in species */
	LDBLE o;					/* stoichiometric coefficient of O in species */
	LDBLE dha, dhb, a_f;		/* WATEQ Debye Huckel a and b-dot; active_fraction coef for exchange species */
	LDBLE lk;					/* log10 k at working temperature */
	LDBLE logk[MAX_LOG_K_INDICES];				/* log kt0, delh, 6 coefficients analytical expression + volume terms */
/* VP: Density Start */
	LDBLE millero[7];		    /* regression coefficients to calculate temperature dependent phi_0 and b_v of Millero density model */
	/* VP: Density End */
	DELTA_H_UNIT original_units;	/* enum with original delta H units */
	int count_add_logk;
	struct name_coef *add_logk;
	LDBLE lg;					/* log10 activity coefficient, gamma */
	LDBLE lg_pitzer;			/* log10 activity coefficient, from pitzer calculation */
	LDBLE lm;					/* log10 molality */
	LDBLE la;					/* log10 activity */
	LDBLE dg;					/* gamma term for jacobian */
	LDBLE dg_total_g;
	LDBLE moles;				/* moles in solution; moles/mass_water = molality */
	int type;					/* flag indicating presence in model and types of equations */
	int gflag;					/* flag for preferred activity coef eqn */
	int exch_gflag;				/* flag for preferred activity coef eqn */
	struct elt_list *next_elt;	/* pointer to next element */
	struct elt_list *next_secondary;
	struct elt_list *next_sys_total;
	int check_equation;			/* switch to check equation for charge and element balance */
	struct reaction *rxn;		/* pointer to data base reaction */
	struct reaction *rxn_s;		/* pointer to reaction converted to secondary and primary
								   master species */
	struct reaction *rxn_x;		/* reaction to be used in model */
	LDBLE tot_g_moles;			/* (1 + sum(g)) * moles */
	LDBLE tot_dh2o_moles;		/* sum(moles*g*Ws/Waq) */
	LDBLE cd_music[5];
	LDBLE dz[3];
	DELTA_V_UNIT original_deltav_units;
};
struct logk
{								/* Named log K's */
	const char *name;					/* name of species */
	LDBLE lk;					/* log10 k at working temperature */
	LDBLE log_k[MAX_LOG_K_INDICES];				/* log kt0, delh, 6 coefficients analalytical expression */
	DELTA_H_UNIT original_units;	/* enum with original delta H units */
	int count_add_logk;
	int done;
	struct name_coef *add_logk;
	LDBLE log_k_original[MAX_LOG_K_INDICES];	/* log kt0, delh, 5 coefficients analalytical expression */
	DELTA_V_UNIT original_deltav_units;
};

/*----------------------------------------------------------------------
 *   Phases
 *---------------------------------------------------------------------- */
struct phase
{								/* all data pertinent to a pure solid phase */
	const char *name;					/* name of species */
	const char *formula;				/* chemical formula */
	int in;						/* species used in model if TRUE */
	LDBLE lk;					/* log10 k at working temperature */
	LDBLE logk[MAX_LOG_K_INDICES];				/* log kt0, delh, 6 coefficients analalytical expression */
	DELTA_H_UNIT original_units;	/* enum with original delta H units */
	DELTA_V_UNIT original_deltav_units;
	int count_add_logk;
	struct name_coef *add_logk;
	LDBLE moles_x;
	LDBLE delta_max;
	LDBLE p_soln_x;
	LDBLE fraction_x;
	LDBLE log10_lambda, log10_fraction_x;
	LDBLE dn, dnb, dnc;
	LDBLE gn, gntot;
	LDBLE gn_n, gntot_n;
	LDBLE t_c, p_c, omega; /* gas: critical TK, critical P(atm), Pitzer acentric coeff */
	LDBLE pr_a, pr_b, pr_alpha;		/* Peng-Robinson parm's */
	LDBLE pr_tk, pr_p;			/* Temperature (K), Pressure (atm) */
	LDBLE pr_phi;				/* fugacity coefficient (-) */
	LDBLE pr_aa_sum2;			/* for calculating multicomponent phi */
	LDBLE delta_v[9];			/* delta_v[0] = [1] + [2]*T + [3]/T + [4]*log10(T) + [5]/T^2 + [6]*T^2 + [7]*P */
	LDBLE pr_si_f;				/* si adapter: log10(phi) - delta_v[0] * (P - 1) /RT */
	bool pr_in;					/* Peng-Robinson in the calc's, or not */

	int type;					/* flag indicating presence in model and types of equations */
	struct elt_list *next_elt;	/* pointer to list of elements in phase */
	struct elt_list *next_sys_total;
	int check_equation;			/* switch to check equation for charge and element balance */
	struct reaction *rxn;		/* pointer to data base reaction */
	struct reaction *rxn_s;		/* pointer to reaction converted to secondary and primary
								   master species */
	struct reaction *rxn_x;		/* reaction to be used in model */
	int replaced;               /* equation contains solids or gases */
	int in_system;
};
/*----------------------------------------------------------------------
 *   Master species
 *---------------------------------------------------------------------- */
 struct master
 {								/* list of name and number of elements in an equation */
 	int in;						/* TRUE if in model, FALSE if out, REWRITE if other mb eq */
 	int number;					/* sequence number in list of masters */
 	int last_model;				/* saved to determine if model has changed */
 	int type;					/* AQ or EX */
 	int primary;				/* TRUE if master species is primary */
 	LDBLE coef;					/* coefficient of element in master species */
 	LDBLE total;				/* total concentration for element or valence state */
 	LDBLE isotope_ratio;
 	LDBLE isotope_ratio_uncertainty;
 	int isotope;
 	LDBLE total_primary;
 	/*    LDBLE la;  */ /* initial guess of master species log activity */
 	struct element *elt;		/* element structure */
 	LDBLE alk;					/* alkalinity of species */
 	LDBLE gfw;					/* default gfw for species */
 	const char *gfw_formula;			/* formula from which to calcuate gfw */
 	struct unknown *unknown;	/* pointer to unknown structure */
 	struct species *s;			/* pointer to species structure */
 	struct reaction *rxn_primary;	/* reaction writes master species in terms of primary
 									   master species */
 	struct reaction *rxn_secondary;	/* reaction writes master species in terms of secondary
 									   master species */
	const char * pe_rxn;
 	int minor_isotope;
};
/*----------------------------------------------------------------------
 *   Unknowns
 *---------------------------------------------------------------------- */
struct unknown
{
	int type;
	LDBLE moles;
	LDBLE ln_moles;
	LDBLE f;
	LDBLE sum;
	LDBLE delta;
	LDBLE la;
	int number;
	const char *description;
	struct master **master;
	struct phase *phase;
	LDBLE si;
	int n_gas_phase_user;
	struct species *s;
	const char * exch_comp;
	const char *pp_assemblage_comp_name;
	void *pp_assemblage_comp_ptr; 
	const char * ss_name;
	void *ss_ptr; 
	const char * ss_comp_name;
	void *ss_comp_ptr;
	int ss_comp_number;
	int ss_in;
	const char *surface_comp;
	const char *surface_charge;
	LDBLE related_moles;
	struct unknown *potential_unknown, *potential_unknown1,
		*potential_unknown2;
	int count_comp_unknowns;
	struct unknown **comp_unknowns;	/* list for CD_MUSIC of comps that contribute to 0 plane mass-balance term */
	struct unknown *phase_unknown;
	LDBLE mass_water;
	int dissolve_only;
	LDBLE inert_moles;
	LDBLE V_m;
	LDBLE pressure;
	int mb_number;
	int iteration;
};

/*----------------------------------------------------------------------
 *   Reaction work space
 *---------------------------------------------------------------------- */
struct reaction_temp
{
	LDBLE logk[MAX_LOG_K_INDICES];
	LDBLE dz[3];
	struct rxn_token_temp *token;
};
struct rxn_token_temp
{								/* data for equations, aq. species or minerals */
	const char *name;					/* pointer to a species name (formula) */
	LDBLE z;					/* charge on species */
	struct species *s;
	struct unknown *unknown;
	LDBLE coef;					/* coefficient of species name */
};
struct unknown_list
{
	struct unknown *unknown;
	LDBLE *source;
	LDBLE *gamma_source;
	/*    int row; */
	/*    int col; */
	LDBLE coef;
};
/* ----------------------------------------------------------------------
 *   Print
 * ---------------------------------------------------------------------- */
struct prints
{
	int all;
	int initial_solutions;
	int initial_exchangers;
	int reactions;
	int gas_phase;
	int ss_assemblage;
	int pp_assemblage;
	int surface;
	int exchange;
	int kinetics;
	int totals;
	int eh;
	int species;
	int saturation_indices;
	int irrev;
	int mix;
	int reaction;
	int use;
	int logfile;
	int punch;
	int status;
	int inverse;
	int dump;
	int user_print;
	int headings;
	int user_graph;
	int echo_input;
	int warnings;
	int initial_isotopes;
	int isotope_ratios;
	int isotope_alphas;
	int hdf;
	int alkalinity;
};
/* ----------------------------------------------------------------------
 *   RATES
 * ---------------------------------------------------------------------- */
struct rate
{
	const char *name;
	char *commands;
	int new_def;
	void *linebase;
	void *varbase;
	void *loopbase;
};
/* ----------------------------------------------------------------------
 *   GLOBAL DECLARATIONS
 * ---------------------------------------------------------------------- */
struct spread_row
{
	int count;
	int empty, string, number;
	char **char_vector;
	LDBLE *d_vector;
	int *type_vector;
};
struct defaults
{
	LDBLE temp;
	LDBLE density;
	const char *units;
	const char *redox;
	LDBLE ph;
	LDBLE pe;
	LDBLE water;
	int count_iso;
	struct iso *iso;
	LDBLE pressure;	/* pressure in atm */
};
struct spread_sheet
{
	struct spread_row *heading;
	struct spread_row *units;
	int count_rows;
	struct spread_row **rows;
	struct defaults defaults;
};
/* ----------------------------------------------------------------------
 *   ISOTOPES
 * ---------------------------------------------------------------------- */
struct master_isotope
{
	const char *name;
	struct master *master;
	struct element *elt;
	const char *units;
	LDBLE standard;
	LDBLE ratio;
	LDBLE moles;
	int total_is_major;
	int minor_isotope;
};
struct calculate_value
{
	const char *name;
	LDBLE value;
	char *commands;
	int new_def;
	int calculated;
	void *linebase;
	void *varbase;
	void *loopbase;
};
struct isotope_ratio
{
	const char *name;
	const char *isotope_name;
	LDBLE ratio;
	LDBLE converted_ratio;
};
struct isotope_alpha
{
	const char *name;
	const char *named_logk;
	LDBLE value;
};
struct system_species
{
	char *name;
	char *type;
	LDBLE moles;
};

/* tally.c ------------------------------- */
struct tally_buffer
{
	const char *name;
	struct master *master;
	LDBLE moles;
	LDBLE gfw;
};
struct tally
{
	const char *name;
	enum entity_type type;
	const char *add_formula;
	LDBLE moles;
	struct elt_list *formula;
	/*
	 * first total is initial
	 * second total is final
	 * third total is difference (final - initial)
	 */
	struct tally_buffer *total[3];
};

/* transport.c ------------------------------- */
struct spec
{
	const char *name;					/* name of species */
	const char *aq_name;				/* name of aqueous species in EX species */
	int type;					/* type: AQ or EX */
	LDBLE a;					/* activity */
	LDBLE lm;					/* log(concentration) */
	LDBLE lg;					/* log(gamma) */
	LDBLE c;					/* concentration for AQ, equivalent fraction for EX */
	LDBLE z;					/* charge number */
	LDBLE Dwt;					/* temperature corrected free water diffusion coefficient, m2/s */
	LDBLE erm_ddl;				/* enrichment factor in ddl */
};
struct sol_D
{
	int count_spec;				/* number of aqueous + exchange species */
	int count_exch_spec;		/* number of exchange species */
	LDBLE exch_total, x_max;	/* total moles of X-, max X- in transport step in sol_D[1] */
	struct spec *spec;
};
struct J_ij
{
	const char *name;
	LDBLE tot1, tot2;
};
struct M_S
{
	const char *name;
	LDBLE tot1, tot2;
};
// Pitzer definitions
typedef enum
{ TYPE_B0, TYPE_B1, TYPE_B2, TYPE_C0, TYPE_THETA, TYPE_LAMDA, TYPE_ZETA,
  TYPE_PSI, TYPE_ETHETA, TYPE_ALPHAS, TYPE_MU, TYPE_ETA, TYPE_Other,
  TYPE_SIT_EPSILON, TYPE_SIT_EPSILON_MU
} pitz_param_type;

struct pitz_param
{
	const char *species[3];
	int ispec[3];
	pitz_param_type type;
	LDBLE p;
	union
	{
		LDBLE b0;
		LDBLE b1;
		LDBLE b2;
		LDBLE c0;
		LDBLE theta;
		LDBLE lamda;
		LDBLE zeta;
		LDBLE psi;
		LDBLE alphas;
		LDBLE mu;
		LDBLE eta;
	  LDBLE eps;
	  LDBLE eps1;
	} U;
	LDBLE a[6];
	LDBLE alpha;
	LDBLE os_coef;
	LDBLE ln_coef[3];
	struct theta_param *thetas;
};

struct theta_param
{
	LDBLE zj;
	LDBLE zk;
	LDBLE etheta;
	LDBLE ethetap;
};

struct const_iso
{
	const char *name;
	LDBLE value;
	LDBLE uncertainty;
};

#endif /* _INC_GLOBAL_STRUCTURES_H  */

