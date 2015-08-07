//-------------------------------------------------------------------
// $Id: s_kinmet.h 903 2013-11-14 12:12:55Z kulik $
//
// Declaration of TKinMet class for kinetics/metastability models
//
// Copyright (C) 2012-2013  D.Kulik, B.Thien, G.Kosakowski
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//------------------------------------------------------------------
//
#ifndef S_KINMET_H
#define S_KINMET_H

#include <vector>

const int MAXDCNAME_ = 16, MAXPHNAME_ = 16, MAXSYMB_ = 4;   // see also v_mod.h

enum kinmet_controls {   /// Re-declared codes to control kinetic rate models (see also m_phase.h)

    KM_UNDEF_ = 'N',      /// not defined, no account for
KinProCode_ = 2,
    KM_PRO_MWR_ = 'M',     /// Kinetics of generic dissolution/precipitation (no uptake, ionex, adsorption)
    KM_PRO_UPT_ = 'U',     /// Kinetics of uptake/entrapment (of minor/trace element) into solid solution
    KM_PRO_IEX_ = 'X',     /// Kinetics of ion exchange (clays, C-S-H, zeolites, ...)
    KM_PRO_ADS_ = 'A',     /// Kinetics of adsorption (on MWI), redox
    KM_PRO_NUPR_ = 'P',    /// Kinetics of nucleation and precipitation
KinModCode_ = 3,
    KM_MOD_TST_ = 'T',     /// Generic TST dissolution/precipitation model following Shott ea 2012
    KM_MOD_PAL_ = 'P',     /// Dissolution/precipitation model of the form (Palandri 2004)
    KM_MOD_WOL_ = 'W',     /// Linear-rate growth/dissolution model (e.g. calcite, Wolthers 2012)
    KM_MOD_NUGR_ = 'U',    /// Mineral nucleation and growth model with nuclei/particle size distr. (TBD)
KinSorpCode_ = 4,
    KM_UPT_ENTRAP_ = 'E',  /// Unified entrapment model (Thien,Kulik,Curti 2013)
    KM_UPT_UPDP_ = 'M',    /// DePaolo (2011) uptake kinetics model TBD
    KM_UPT_SEMO_ = 'G',    /// Growth (surface) entrapment model (Watson 2004) TBD
    KM_IEX_FAST_ = 'F',    /// Fast ion exchange kinetics (e.g. montmorillonite, CSH) TBD
    KM_IEX_SLOW_ = 'L',    /// Slow ion exchange kinetics (e.g. illite, zeolites)  TBD
    KM_ADS_INHIB_ = 'I',   /// Adsorption inhibition TBD
    KM_NUCL_SSMP_  = 'P',  /// Solid solution nucleation model (Prieto 2013) TBD
KinLinkCode_ = 5,
    KM_LNK_SURF_ = 'S',    ///  Link to (fraction of) solid substrate surface area (default)
    KM_LNK_PVOL_ = 'P',    ///  Link to (fraction of) solid substrate (pore) volume
    KM_LNK_MASS_ = 'M',    ///  Link to (fraction of) solid substrate mass
KinSizedCode_ = 6,  // Codes for dependencies of the sphericity shape factor on system variables
    KM_SIZED_ETM_ = 'T',   ///  Empirical f(time) cubic polynomial f = a + bt +ct^2 + dt^3 (default)
    KM_SIZED_ESI_ = 'S',   ///  Empirical f(lgSI) cubic polynomial f = a + bt +ct^2 + dt^3
    KM_SIZED_ESA_ = 'A',   ///  Empirical f(sarea-change) cubic polynomial f = a + bt +ct^2 + dt^3
    KM_SIZED_EVOL_ = 'V',  ///  Empirical f(volume-change) cubic polynomial f = a + bt +ct^2 + dt^3
    KM_SIZED_MASS_ = 'M',  ///  Empirical f(mass-change) cubic polynomial f = a + bt +ct^2 + dt^3
    KM_SIZED_MOL_ = 'X',  ///  Empirical f(amount-change) cubic polynomial f = a + bt +ct^2 + dt^3
    KM_SIZED_UNI_ = 'U',   ///  Uniform particle/pore size distribution (reserved)
    KM_SIZED_BIN_ = 'B',   ///  Binodal particle/pore size distribution (reserved)
    KM_SIZED_FUN_ = 'F',   ///  Empirical particle/pore size distribution function (reserved)
KinResCode_ = 7,     // Units and type of rate constants
    KM_RES_SURF_N_ = 'A',   /// surface-scaled rate constant (k in mol/m2/s), default
    KM_RES_SURF_M_ = 'M',   /// surface-scaled rate constant (k in kg/m2/s)
    KM_RES_PVS_N_  = 'V',   /// pore-volume-scaled rate constant (k in mol/m3/s)
    KM_RES_PVS_M_  = 'W',   /// pore-volume-scaled rate constant (k in kg/m3/s)
    KM_RES_ABS_N_  = 'F',   /// absolute (unscaled) rate constant (k in mol/s)
    KM_RES_ABS_M_  = 'G',   /// absolute (unscaled) rate constant (k in kg/s)
    KM_LIN_RATE_   = 'L',   /// linear growth/dissolution rate constant (v in m/s)

// Minor end member codes for uptake kinetics models
    DC_SOL_MINOR_ = 'J',
    DC_SOL_MINDEP_ = 'F'
};

enum affin_term_op_codes {
    ATOP_CLASSIC_ = 0,       /// classic TST affinity term (see .../Doc/Descriptions/KinetParams.pdf)
    ATOP_CLASSIC_REV_ = 1,   /// classic TST affinity term, reversed
    ATOP_SCHOTT_ = 2,        /// Schott et al. 2012 fig. 1e
    ATOP_HELLMANN_ = 3,      /// Hellmann Tisserand 2006 eq 9
    ATOP_TENG1_ = 4,         /// Teng et al. 2000, eq 13
    ATOP_TENG2_ = 5,         /// Teng et al. 2000, Fig 6
    ATOP_FRITZ_ = 6          /// Fritz et al. 2009, eq 6 nucleation and growth

};

// This data structure should be passed in order to create an instance of TKinMet derived class for a phase
struct KinMetData {

    char  KinProCod_;   /// Code of the kinetic process (derived TKinMet class), see enum kinmet_controls
    char  KinModCod_;   /// Type code of the kinetic/metastability model, see enum kinmet_controls
    char  KinSorpCod_;  /// Type code of sorption kinetics model (solution/sorption phases only), see enum kinmet_controls
    char  KinLinkCod_;  /// Type code of metastability links of this phase to other phases, see enum kinmet_controls
    char  KinSizedCod_; /// Type of particle/pore size distribution and A_s correction, see enum kinmet_controls
    char  KinResCod_;   /// Reserved model control code
    char  PhasNam_[MAXPHNAME_+1];      /// Phase name (for specific built-in models)

    long int NComp_;   /// Number of components in the phase (nDC in m_phase.h, L1[k] in MULTI)
    long int nlPh_;    /// Number of linked phases (cf. lPh), default 0
    long int nlPc_;    /// TKinMet, TSorpMod: number of parameters per linked phase, default 0.

    long int nPRk_;  /// number of "parallel reactions" that affect amount constraints for k-th phase (1, 2, 3, ...), 1 by default
    long int nSkr_;  /// number of (aqueous or gaseous or surface) species from other reacting phases involved, 0 by default
    long int nrpC_;  /// number of parameter (coefficients) involved in 'parallel reaction' terms (0 or 12 + 3res.)
    long int naptC_; /// number of parameter (coefficients) per species involved in 'activity product' terms (0 or 1)
    long int nAscC_; /// number of parameter coefficients in specific surface area correction equation ( 0 to 5 )
    long int nFaceC_; /// number of (separately considered) crystal faces or surface patches ( 1 to 4 )
//    long int numpC_; /// number of uptake model parameter coefficients (per end member)
//    long int iRes4_;  // reserved

    double T_k_;     /// Temperature, K (initial)
    double P_bar_;   /// Pressure, bar (initial)
    double kTau_;    /// current time, s (initial)
    double kdT_;     /// current time step (initial)
  //
    double IS_;      /// Effective molal ionic strength of aqueous electrolyte
    double pH_;      /// pH of aqueous solution
    double pe_;      /// pe of aqueous solution
    double Eh_;      /// Eh of aqueous solution, V
  //
    double nPh_;     /// current amount of this phase, mol (read-only)
    double mPh_;     /// current mass of this phase, kg (read-only)
    double vPh_;     /// current volume of this phase, m3 (read-only)
    double sAPh_;    /// current surface of this phase, m2
    double LaPh_;    /// phase stability index (log scale)
    double OmPh_;    /// phase stability index (activity scale) 10^LaPh_
double sFact_;   /// phase surface area - volume shape factor ( >= 4.836 for spheres )
    double sSA_;    /// Specific surface area of the phase, m2/kg, default: 0.
    double sgw_;    /// Standard mean surface energy of solid-aqueous interface, J/m2
    double sgg_;    /// Standard mean surface energy of gas-aqueous interface, J/m2
    double rX0_;    /// Mean radius r0 for (spherical or cylindrical) particles, m (reserved)
    double hX0_;    /// Mean thickness h0 for cylindrical or 0 for spherical particles, m (reserved)
    double sVp_;    /// Specific pore volume of phase, m3/g (default: 0)
    double sGP_;    /// surface free energy of the phase, J (YOF*PhM)
    double nPul_;   /// upper restriction to this phase amount, mol (calculated here)
    double nPll_;   /// lower restriction to this phase amount, mol (calculated here)

    double *arlPhc_;   /// TsolMod, TKinMet, TSorpMod: pointer to input array of phase link parameters [nlPh*nlPc]
    double *arfeSAr_;  /// Pointer to input fractions of surface area of the solid related to different parallel reactions [nPRk] read-only
    double *arrpCon_;  /// Pointer to input array of kinetic rate constants for faces and 'parallel reactions' [nPRk*nrpC] read-only
    double *arapCon_;  /// Pointer to array of parameters per species involved in 'activity product' terms [nPRk * nSkr*naptC] read-only
    double *arAscp_;   /// Pointer to array of parameter coefficients of equation for correction of specific surface area [nAscC] read-only
    // new:new: array of nucleation model parameters (A.Testino?)

    char  (*SM_)[MAXDCNAME_];  /// pointer to the classifier of DCs involved in sorption phase [NComp] read-only
    char  *arDCC_;       /// pointer to the classifier of DCs involved in the phase [NComp] read-only

    long int (*arPhXC_)[2];  /// TKinMet: linked phase indexes and linkage type codes [nlPh][2]

    long int (*arocPRk_)[2]; /// pointer to operation codes for kinetic parallel reaction affinity terms [nPRk] read-only
    long int *arxSKr_;  /// pointer to input array of DC indexes used in activity products [nSKr_] read-only

    double *arym_;    /// Pointer to molalities of all species in MULTI (provided), read-only
    double *arla_;    /// Pointer to lg activities of all species in MULTI (provided), read-only

double *arxp_;   /// Pointer to amounts of all phases in MULTI (provided), read-only   mol
double *armp_;   /// Pointer to masses of all phases in MULTI (provided), read-only    g
double *arvp_;   /// Pointer to volumes of all phases in MULTI (provided), read-only   cm3
double *arasp_;  /// Pointer to (current) specific surface areas of all phases in MULTI (provided), read-only

    double *arnx_;    /// Pointer to mole amounts of phase components (provided) [NComp] read-only

    double *arnxul_;  /// Vector of upper kinetic restrictions to nx, moles [NComp]  (DUL) direct access output
    double *arnxll_;  /// Vector of lower kinetic restrictions to nx, moles [NComp]  (DLL) direct access output

    double *arWx_;    /// Species (end member) mole fractions ->NComp
    double *arVol_;   /// molar volumes of end-members (species) cm3/mol ->NSpecies

};

// Class describing a 'parallel reaction' region kinetic rate law data (Schott ea 2012 general form)
struct TKinReact
{
    // Input data
     long int xPR;   /// index of this parallel reaction
     long int nSa;   // the same as nSKr - tot. number of species used in activity products

     long int ocPRk[2]; /// operation code for this kinetic parallel reaction affinity term
     long int *xSKr;  /// pointer to input array of DC indexes used in activity products [nSKr] (copy)

     double feSAr;   /// input fraction of surface area of the solid related to this parallel reaction
     double *rpCon;  /// pointer to input array of kinetic rate constants for this 'parallel reaction' [nrpC]
     double **apCon; /// pointer input array of parameters per species involved in 'activity product' terms [nSkr][naptC]

     // work data: unpacked rpCon[nrpC]
     double ko,  /// kod net rate constant at standard temperature (mol/m2/s)
                ///    if k > 0 it is used for dissolution only; if k < 0 it is used for precipitation only
            Ko,  /// kop gross rate constant at standard temperature (mol/m2/s) (if K != 0 it is used in any case)
                /// k or K convention: negative sign: precipitation; positive sign: dissolution
            Ap,  /// Arrhenius parameter
            Ea,  /// activation energy at st.temperature J/mol
            bI,  /// Ionic strength power coefficient (default 0)
            bpH, /// pH power coefficient (default 0)
            bpe, /// pe power coefficient (default 0)
            bEh, /// Eh power coefficient (default 0)
            pPR, /// Reaction order power coefficient for far-from-equilibrium cases
            qPR, /// reaction order power coefficient in the affinity term (default 1)
            mPR, /// second reaction order parameter in the affinity term (default 0)
            uPR, /// parameter constant in the affinity term
            OmEff, /// Effective saturation index for a simple nucleation model (default 1)
            nucRes; // reserved
//            Omg; /// Input stability index non-log (d-less)

     // Results of rate term calculation
     double
        arf,  // Arrhenius factor (temperature correction on kappa)
        cat,  // catalytic product term (f(prod(a))
        aft,  // affinity term (f(Omega))

        k,   // dissolution rate constant (corrected for T) in mol/m2/s
        K,   // precipitation rate constant (corrected for T) in mol/m2/s
        rPR,   // rate for this region (output) in mol/s
        rmol   // rate for the whole face (output) in mol/s
//        velo,   // velocity of face growth (positive) or dissolution (negative) m/s
        ;
};

class TKinMet  // Base class for MWR kinetics and metastability models
{
    protected:
    char  KinProCode;   /// Code of the kinetic process (derived TKinMet class), see enum kinmet_controls
    char  KinModCode;   /// Type code of the kinetic/metastability model, see enum kinmet_controls
    char  KinSorpCode;   /// Type code of sorption kinetics model (solution/sorption phases only), see enum kinmet_controls
    char  KinLinkCode;   /// Type code of metastability links of this phase to other phases, see enum kinmet_controls
    char  KinSizedCode; /// Type of particle/pore size distribution and A_s correction, see enum kinmet_controls
    char  KinResCode;   /// Reserved model control code
    char  PhasName[MAXPHNAME_+1];      /// Phase name (for specific built-in models)

    long int NComp;   /// Number of components in the phase (nDC in m_phase.h, L1[k] in MULTI)
    long int nlPh;    /// Number of linked phases (cf. arPhXC), default 0
    long int nlPc;    /// TKinMet, TSorpMod: number of parameters per linked phase (<=7), default 0.

    long int nPRk;     /// number of 'parallel reactions' that affect amount constraints for k-th phase (1, 2, 3, ...), 1 by default
    long int nSkr;     /// number of (aqueous or gaseous) species from other reacting phases involved, 0 by default
    long int nrpC;     /// number of parameter (coefficients) involved in 'parallel reaction' terms (0 or 12 + 3res.)
    long int naptC;    /// number of parameter (coefficients) per species involved in 'activity product' terms (0 or 1)
    long int nAscC;    /// number of parameter coefficients in specific surface area correction equation ( 0 to 5 )
    long int nFaceC;   /// number of (separately considered) crystal faces or surface patches ( 1 to 4 )
//    long int numpC;   /// number of sorption/uptake model parameter coefficients (per end member)
//    long int iRes4;   // reserved
    double OmgTol;      /// Tolerance for checking dissolution or precipitation cases (default 1e-6)
    double R_CONST;     /// Gas constant, 8.31451 J/K/mol
    double T_k;         /// Temperature, K
    double P_bar;       /// Pressure, bar
    double kTau;        /// current time, s
    double kdT;         /// current time increment, s
    double sSAi;        /// initial specific surface area of the phase (m2/kg)
    double nPhi;        /// initial amount of the phase, mol
    double mPhi;        /// initial mass of the phase, kg
    double vPhi;        /// initial volume of the phase, m3
double sFacti;      /// inital sphericity factor (0 < sFacti < 1), negative for pores
double sFact;       /// current sphericity factor (0 < sFact < 1)
                    // all three for the built-in correction of specific surface area and kTot-vTot transformation
    // These values will be corrected after GEM converged at each time step
    double IS;          /// Effective molal ionic strength of aqueous electrolyte
    double pH;          /// pH of aqueous solution
    double pe;          /// pe of aqueous solution
    double Eh;          /// Eh of aqueous solution, V
    double nPh;     /// current amount of this phase, mol
    double mPh;     /// current mass of this phase, kg
    double vPh;     /// current volume of this phase, m3
    double sAPh;    /// current surface area of this phase, m2
    double LaPh;    /// phase stability index (log scale)
    double OmPh;    /// phase stability index (activity scale) 10^LaPh

    // These values may be corrected inside of TKinMet class instance over time steps
    double sSA;    /// Specific surface area of the phase, m2/kg, default: 0.
    double sgw;    /// Standard mean surface energy of solid-aqueous interface, J/m2
    double sgg;    /// Standard mean surface energy of gas-aqueous interface, J/m2
    double rX0;    /// Mean radius r0 for (spherical or cylindrical) particles, m (reserved)
    double hX0;    /// Mean thickness h0 for cylindrical or 0 for spherical particles, m (reserved)
    double sVp;    /// Specific pore volume of phase, m3/kg (default: 0)
    double sGP;    /// surface free energy of the phase, J (YOF*PhM)
    double nPul;   /// upper restriction to this phase amount, mol (output)
    double nPll;   /// lower restriction to this phase amount, mol (output)

    double **arlPhc; /// TsolMod, TKinMet, TSorpMod: pointer to input array of phase link parameters [nlPh*nlPc]
    double *arfeSAr;   /// input fractions of surface area of the solid related to different parallel reactions [nPRk] - input
    double **arrpCon;  /// input array of kinetic rate constants for faces and 'parallel reactions' [nPRk*nrpC]
    double ***arapCon; /// input array of parameters per species involved in 'activity product' terms [nPRk * nSkr*naptC]
    double *arAscp;  /// input array of coefficients of the sFact function for correction of specific surface area [nAscC]
    // new:new: array of nucleation model parameters (A.Testino?)

    char  (*SM)[MAXDCNAME_];  /// pointer to the list of DC names in the phase [NComp] read-only
    char  *arDCC;       /// pointer to the classifier of DCs involved in the phase [NComp] read-only

    long int (*arPhXC)[2];  /// TKinMet: linked phase indexes and linkage type codes [nlPh][2]

    long int (*arocPRk)[2]; /// input operation codes for kinetic parallel reaction affinity terms [nPRk][2]
    long int *arxSKr;  /// pointer to input array of DC indexes used in activity products [nSKr]

    double *arym;    /// Pointer to molalities of all species in MULTI (provided), read-only
    double *arla;    /// Pointer to lg activities of all species in MULTI (provided), read-only

double *arxp;   /// Pointer to amounts of all phases in MULTI (provided), read-only   mol
double *armp;   /// Pointer to masses of all phases in MULTI (provided), read-only    g
double *arvp;   /// Pointer to volumes of all phases in MULTI (provided), read-only   cm3
double *arasp;  /// Pointer to (current) specific surface areas of all phases in MULTI (provided), read-only

    double *arnx;   /// Pointer to mole amounts of phase components (provided) [NComp] read-only mol

    double *arnxul; /// Vector of upper kinetic restrictions to nx, moles [NComp]  (DUL) direct access output
    double *arnxll; /// Vector of lower kinetic restrictions to nx, moles [NComp]  (DLL) direct access output

    double *arWx;   /// Species (end member) mole fractions ->NComp
    double *arVol;  /// molar volumes of end-members (species) cm3/mol ->NComp

// Work data and kinetic law calculation results
    TKinReact *arPRt; /// work array of parameters and results for 'parallel reaction' terms [nPRk]

    double *spcfu;  /// work array of coefficients for splitting nPul and nPll into nxul and nxll [NComp]
    double *spcfl;  /// work array of coefficients for splitting nPul and nPll into nxul and nxll [NComp]

    double kTot;    /// Total specific MWR (precipitation/dissolution) rate Rn (mol/m2/s) - output
    double gTot;    /// Total specific MWR rate RG (kg/m2/s)
    double rTot;    /// Current total MWR rate (mol/s) Rn*sAPh
    double vTot;    /// Total orthogonal surface propagation velocity RL (m/s) - output

//  Work data for SSA correction over time step
    double pi;      /// Pi = 3.14159265358979;
    double Np;      /// estimated number of particles
    double Rho;     /// current density of this phase, kg/m3
    double d32;     /// Sauter equivalent diameter of particles, m
    double d43;     /// De Brouckere equivalent diameter, m
    double dsv;     /// Surface-volume equivalent diameter, m
    double dsvi;    /// Initial surface-volume equivalent diameter, m
    double Vp;      /// equivalent volume of one particle, m3
    double Ap;      /// surface area of one partcle, m2
    double Mp;      /// estimated mass of one partcle, kg
    double sSAV;    /// specific surface per unit volume in m2/m3
    double sSAVcor; /// corrected specific surface per unit volume in m2/m3
    double sSAcor;  /// Corrected specific surface area (m2/kg) - output
    double sAph_c;  /// Corrected surface area of the phase (m2/kg)
    double kdT_c;   /// new: modified (suggested) time increment, s

//  Work data for totals for linked phases
    double sSAlp;   /// specific surface per unit mass in m2/kg
    double sSAVlp;  /// specific surface per unit volume in m2/m3
    double sAPhlp;  /// total surface area, m2
    double mPhlp;   /// total mass, kg
    double vPhlp;   /// total volume, m3
    double Rholp;   /// total density, kg/m3
    double nPhlp;   /// total moles
    double nPhlpi;  /// initial total moles
    double mPhlpi;  /// initial total mass, kg
    double vPhlpi;  /// initial total volume, m3

    // SS dissolution

    // SS precipitation

    // (SS) nucleation

    // functions for allocation and initialization of kinetic rate tables
    void alloc_kinrtabs();
    long int init_kinrtabs( double *p_arlPhc, double *p_arrpCon,  double *p_arapCon );
    void free_kinrtabs();
    // functions for allocation and initialization of the  TKinReact array for parallel reactions data
    void alloc_arPRt();
    void init_arPRt();
    void free_arPRt();

    // Calculation of total properties of linked phases
    bool linked_phases_properties( bool if_init );
    // Calculates mean equivalent properties of particles or pores
    bool particles_pores_properties( bool if_init );

  public:
    // Generic constructor
    TKinMet( const KinMetData *kmd );

    // Destructor
    virtual ~TKinMet();

    virtual   // Calculates temperature/pressure corrections to kinetic rate constants
    bool PTparam( const double TK, const double P );

    virtual  // Calculates phase dissolution/precipitation/nucleation rates
    bool RateMod();

    virtual  // Calculates initial phase dissolution/precipitation/nucleation rates
    bool RateInit();

    virtual  // Calculates (splits) rates to lower/upper metastability constraints
    bool SplitMod();

    virtual  // Calculates (splits) rates to lower/upper metastability constraints
    bool SplitInit();

    virtual // sets the new metastability constraints if the time step is accepted
    bool SetMetCon();

    // Sets new specific surface area of the phase As and other properties;
    // also sets 'parallel reactions' area fractions
    // returns false if these parameters in TKinMet instance did not change; true if they did.
    //
    virtual
    bool UpdateFSA(const double pAsk, const double pXFk, const double pFWGTk, const double pFVOLk,
                    const double pLgOm, const double p_fFact, const double pYOFk,
                    const double pICa, const double ppHa, const double ppea, const double pEha );

    // Returns (modified) specific surface area of the phase, metastability constraints (via parameters),
    // and sends internally back (modified) 'parallel reactions' area fractions.
    virtual
    double GetModFSA(double& p_fFact, double& prTot, double& pkTot, double& pvTot,
                     double& pPULk, double& pPLLk);

    // Updates temperature to T_K and pressure to P_BAR;
    // calculates Arrhenius factors and temperature-corrected rate constants in all PR regions.
    virtual
    long int UpdatePT(const double T_K, const double P_BAR );

    // sets new time Tau and time step dTau
    // returns false if neither kTau nor kdT changed; true otherwise
    virtual
    bool UpdateTime( const double Tau, const double dTau );

    // Returns (modified) time step value that satisfies back-coupling criterion
    //    TBD
    virtual
    double GetMdTime( )
    {
      return kdT_c;
    }

    // Checks dimensions in order to re-allocate class instance, if necessary
    virtual
    bool testSizes( const KinMetData *kmd );

    // calculates the rate constant for r-th parallel reaction
    virtual
    double PRrateCon(TKinReact &rk, const long int r );

    // calculates corrected specific surface area
    virtual
    double CorrSpecSurfArea(const double sFratio, const bool toinit );

};


// - - - - - - - - - - - - - - - - - - - - - - - - - - - derived classes - - - - - - - - - - - - - -

class TMWReaKin: public TKinMet  // Generic MWR kinetics models no nucleation/uptake/sorption TBD
{
    private:


    // SS dissolution

    // SS precipitation

    // (SS) nucleation


    // specific stuff for uptake kinetics (TBD)
            // internal functions
        //        void alloc_internal();
        //        void free_internal();

     public:
            // Constructor
            TMWReaKin( const KinMetData *kmd /*, specific params */ ):TKinMet( kmd )
            {
            }

            // Destructor
            ~TMWReaKin(){}

            // Initializes uptake rates
            bool SSReaKinInit();

            // Calculates uptake rates
            bool SSReaKinMod();

};

class TUptakeKin: public TKinMet  // SS uptake kinetics models Thien,Kulik,Curti 2013
{
    private:
    // specific stuff for uptake kinetics
    long int numpC;    /// number of sorption/uptake model parameter coefficients (per end member)
    long int nElm;     /// number of independent components in IPM work data structure
    long int iRes4;    // reserved

    long int *arxTrDC; /// pointer to input array of aq DC indexes for end-members [NComp]
    long int *arxICu;  /// pointer to input array of aq IC indexes for end members [NComp]

    double **arUmpCon; /// input array of uptake model coefficients [NComp*numpC] read-only
    double *arElm;     /// pointer to total molalities of elements (IC) dissolved in aqueous phase [nElem]

    double *Rdj;       /// pointer to vector of output Rd values (for t-th time step) [NComp], direct access
    double *Dfj;       /// pointer to vector of output Delta(j,rest) values (for t-th time step) [NComp], direct access

    // Uptake model output

    // internal functions
    void alloc_upttabs();
    long int init_upttabs( double *p_arUmpCon );
    void free_upttabs();

    public:

    // Constructor
    TUptakeKin(KinMetData *kmd, long int p_numpC, long int p_nElm, double *p_arUmpCon, long int *p_arxICu,
               double *p_arElm, double *p_Rdj, double *p_Dfj /*, specific params */ );

    // Destructor
    ~TUptakeKin();

    // Calculates temperature/pressure corrections to kinetic rate constants
    bool UptKinPTparam( const double TK, const double P  );

    // Initializes uptake rates
    bool UptakeInit();

    // Calculates uptake rates
    bool UptakeMod();

};


class TIonExKin: public TKinMet  // Ion exchange (layered) kinetics model TBD
{
    private:

    // specific stuff for uptake kinetics

    // internal functions
    //        void alloc_internal();
    //        void free_internal();

    public:

    // Constructor
    TIonExKin( KinMetData *kmd /*, specific params */ ):TKinMet( kmd )
    {

    }

    // Destructor
    ~TIonExKin(){}

    // Calculates ion exchange rates
    long int IonExRatesMod();


};


class TAdsorpKin: public TKinMet  // Adsorption (layered) kinetics model TBD
{
    private:

    // specific stuff for adsorption kinetics

    // internal functions
    //        void alloc_internal();
    //        void free_internal();

    public:

    // Constructor
    TAdsorpKin( KinMetData *kmd /*, specific params */ ):TKinMet( kmd )
    {

    }

    // Destructor
    ~TAdsorpKin(){}

    // Calculates uptake rates
    long int AdsorpRatesMod();


};


class TNucleKin: public TKinMet  // Mineral nucleation/growth kinetics models TBD
{
    private:

    // specific stuff for uptake kinetics
  // new:new: array of nucleation model parameters (A.Testino?)

    // internal functions
    //        void alloc_internal();
    //        void free_internal();

    public:

    // Constructor
    TNucleKin( KinMetData *kmd /*, specific params */ ):TKinMet( kmd )
    {

    }

    // Destructor
    ~TNucleKin(){}

    // Calculates uptake rates
    long int NucleGrowthMod();

};

// More derived classes can be added here, if necessary


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif // S_KINMET_H
