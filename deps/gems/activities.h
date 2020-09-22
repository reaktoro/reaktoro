//-------------------------------------------------------------------
// $Id: activities.h 879 2013-10-10 14:47:33Z kulik $
/// \file activities.h
/// Contains definition of the ACTIVITY structure - module for computing
/// activities for GEM IPM 3 and new speciation algorithms.
//
/// \struct ACTIVITY activity.h
/// Defines the structure of node-dependent data for
/// exchange of activity terms between the node instance and the GEM algorithms.
/// Used in TNode class.
//
// Copyright (c) 2014- by D.Kulik, A.Leal
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
//-------------------------------------------------------------------
//
#ifndef _Activity_H_
#define _Activity_H_

#include "datach.h"
#include "databr.h"
#include "s_solmod.h"
#include "s_sorpmod.h"
#include "m_const.h"

/*
const double R_CONSTANT = 8.31451,
              NA_CONSTANT = 6.0221367e23,
                F_CONSTANT = 96485.309,
                  e_CONSTANT = 1.60217733e-19,
                    k_CONSTANT = 1.380658e-23,
// Conversion factors
                      cal_to_J = 4.184,
                        C_to_K = 273.15,
                          lg_to_ln = 2.302585093,
                            ln_to_lg = 0.434294481,
                             H2O_mol_to_kg = 55.50837344,
                               Min_phys_amount = 1.66e-24;
*/
typedef struct  /// ACTIVITY - data structure for computing DC activities in phases
{               /// DATACH indexation throughout, no I/O file exchange

//      Usage of this variable                                              node GEM-in GEM-out node
   double
// \section Chemical scalar variables
    TK,     ///< Node temperature T (Kelvin)                     	         +      +      -     -
    P, 	    ///< Node Pressure P (Pa)                         	             +      +      -     -
    RT,     // RT product
    ln5551, // ln(H2O_mol_to_kg) = log(55.50837344)
//
    IC,     ///< Effective aqueous ionic strength (molal)                    -      -      +     +
    pH,     ///< pH of aqueous solution in the activity scale (-log10 molal) -      -      +     +
    pe,     ///< pe of aqueous solution in the activity scale (-log10 molal) -      -      +     +
    Eh;     ///< Eh of aqueous solution (V)                                  -      -      +     +

    double Pc, Tc;  // T in K, P in bar
   // Here TSolMod, TSorpMod, TKinMet parameters & coefficients moved from TMULTI MULTI
   //   DK, AL September 2014
      long int
      // TSolMod stuff
      *LsMod, ///< Number of interaction parameters. Max parameter order (cols in IPx),
                ///< and number of coefficients per parameter in PMc table [3*nPS]
      *LsMdc, ///<  for multi-site models: [3*nPS] - number of nonid. params per component;
               /// number of sublattices nS; number of moieties nM
      *LsMdc2, ///<  new: [3*nPS] - number of DQF coeffs; reciprocal coeffs per end member;
              /// reserved
      *IPx,   ///< Collected indexation table for interaction parameters of non-ideal solutions
                ///< ->LsMod[k,0] x LsMod[k,1]   over nPS

      *LsPhl,  ///< new: Number of phase links; number of link parameters; [nPH][2]
      (*PhLin)[2];  ///< new: indexes of linked phases and link type codes (sum 2*LsPhl[k][0] over nPH)

      /* TSolMod !! arrays and counters to be added (for mixed-solvent electrolyte phase) TW

      ncsolv, /// TW new: number of solvent parameter coefficients (columns in solvc array)
      nsolv,  /// TW new: number of solvent interaction parameters (rows in solvc array)
      *ixsolv, /// new: array of indexes of solvent interaction parameters [nsolv*2]
      *solvc, /// TW new: array of solvent interaction parameters [ncsolv*nsolv]

      ncdiel, /// TW new: number of dielectric constant coefficients (colums in dielc array)
      ndiel,  /// TW new: number of dielectric constant parameters (rows in dielc array)
      *ixdiel /// new: array of indexes of dielectric interaction parameters [ndiel*2]
      *dielc, /// TW new: array of dielectric constant parameters [ncdiel*ndiel]

      ndh,    /// TW new: number of generic DH coefficients (rows in dhc array)
      *dhc,   /// TW new: array of generic DH parameters [ndh]
      */
      double
        denW[5],   ///< Density of water, first T, second T, first P, second P derivative for Tc,Pc
        denWg[5],  ///< Density of steam for Tc,Pc
        epsW[5],   ///< Diel. constant of H2O(l)for Tc,Pc
        epsWg[5];  ///< Diel. constant of steam for Tc,Pc
      double
        *PMc,    ///< Collected interaction parameter coefficients for the (built-in) non-ideal mixing models -> LsMod[k,0] x LsMod[k,2]
        *DMc,    ///< Non-ideality coefficients f(TPX) for DC -> L1[k] x LsMdc[k][0]
        *MoiSN,  ///< End member moiety- site multiplicity number tables ->  L1[k] x LsMdc[k][1] x LsMdc[k][2]
        *SitFr;  ///< Tables of sublattice site fractions for moieties -> LsMdc[k][1] x LsMdc[k][2]
      double
      *lPhc,  ///< new: Collected array of phase link parameters (sum(LsPhl[k][1] over Fi)
      *DQFc;  ///< new: Collected array of DQF parameters for DCs in phases -> L1[k] x LsMdc2[k][0]
    //  *rcpc,  ///< new: Collected array of reciprocal parameters for DCs in phases -> L1[k] x LsMdc2[k][1]


      // TSorpMod stuff
      long int
      *LsESmo, ///< new: number of EIL model layers; EIL params per layer; CD coefs per DC; reserved  [nPS][4]
      *LsISmo, ///< new: number of surface sites; isotherm coeffs per site; isotherm coeffs per DC; max.denticity of DC [nPS][4]
      *xSMd;   ///< new: denticity of surface species per surface site (site allocation) (-> L1[k]*LsISmo[k][3]] )
      long int  (*SATX)[4]; ///< Setup of surface sites and species (will be applied separately within each sorption phase) [Lads]
                 /// link indexes to surface type [XL_ST]; sorbent em [XL_EM]; surf.site [XL-SI] and EDL plane [XL_SP]
       double
      *EImc,  ///< new: Collected EIL model coefficients k -> += LsESmo[k][0]*LsESmo[k][1]
      *mCDc,  ///< new: Collected CD EIL model coefficients per DC k -> += L1[k]*LsESmo[k][2]
      *IsoPc, ///< new: Collected isotherm coefficients per DC k -> += L1[k]*LsISmo[k][2];
      *IsoSc; ///< new: Collected isotherm coeffs per site k -> += LsISmo[k][0]*LsISmo[k][1];
      // TSorpMod & TKinMet stuff
      double
      *SorMc; ///< new: Phase-related kinetics and sorption model parameters: [nPS][16]
              ///< in the same order as from Asur until fRes2 in TPhase

//      Usage of this variable                                                                     node GEM-in GEM-out node
     double *H0;     ///< DC molar enthalpies, reserved [L]
//     double *A0;     ///< DC molar Helmholtz energies, reserved [L]
     double *S0;     ///< DC molar entropies, reserved [L]
     double *Cp0;    ///< DC molar heat capacity, reserved [L]

     double *tpp_G; ///< input Gibbs energy of species in J/mol [nDC]
     double *G0;   ///< Input normalized g0_j(T,P) for DC at unified standard scale[nDC]               +      +      -     -
     double *G;    ///< Normalized DC energy function c(j), incl. primal mu mole/mole [nDC]            -      -      +     +

     double *lnAct;  ///< ln of DC activities (mu - mu0) [nDC]                                         -      -      +     +
     double *lnGam;  ///< ln of DC activity coefficients in unified (mole-fraction) scale [nDC]        -      -      +     +
     double *lnGmo;   ///< Copy of lnGam from previous IPM iteration (reserved)
    // TSolMod stuff (detailed output on partial energies of mixing)   lnGam components
    double *lnDQFt; ///< new: DQF terms adding to overall activity coefficients [nDCs]               -      -      +     +
    double *lnRcpt; ///< new: reciprocal terms adding to overall activity coefficients [nDCs]        -      -      +     +
    double *lnExet; ///< new: excess energy terms adding to overall activity coefficients [nDCs]     -      -      +     +
    double *lnCnft; ///< new: configurational terms adding to overall activity [nDCs]                -      -      +     +
    // Takeover from MULTI
    double *Gamma;   ///< DC activity coefficients in molal or other phase-specific scale [0:L-1]
    double *lnGmf;   ///< ln of initial DC activity coefficients for correcting G0 [0:L-1]
    double *lnGmM;   ///< ln of DC pure gas fugacity (or metastability) coeffs or DDF correction [0:L-1]

    // TSorpMod stuff (TBD)
    double *lnScalT;  ///< new: Surface/volume scaling activity correction terms [nDCs]              -      -      +     +
    double *lnSACT;   ///< new: ln isotherm-specific SACT for surface species [nDCs]                 -      -      +     +
    double *lnGammF;  ///< new: Frumkin or BET non-electrostatic activity coefficients [nDCs]        -      -      +     +
    double *CTerms;   ///< new: Coulombic terms (electrostatic activity coefficients) [nDCs]         -      -      +     +

    double (*VPh)[MIXPHPROPS],     ///< Volume properties for mixed phases [FIs]                     -      -      +     +
           (*GPh)[MIXPHPROPS],     ///< Gibbs energy properties for mixed phases [FIs]               -      -      +     +
           (*HPh)[MIXPHPROPS],     ///< Enthalpy properties for mixed phases [FIs]                   -      -      +     +
           (*SPh)[MIXPHPROPS],     ///< Entropy properties for mixed phases [FIs]                    -      -      +     +
           (*CPh)[MIXPHPROPS],     ///< Heat capacity Cp properties for mixed phases [FIs]           -      -      +     +
           (*APh)[MIXPHPROPS],     ///< Helmholtz energy properties for mixed phases [FIs]           -      -      +     +
           (*UPh)[MIXPHPROPS];     ///< Internal energy properties for mixed phases [FIs]            -      -      +     +
            // MIXPHPROPS: see m_const.h

// Coding - takeover from MULTI
    char (*sMod)[8];   ///< new: Codes for built-in mixing models of multicomponent phases [FIs]
    char  (*dcMod)[6];   ///< Codes for PT corrections for dependent component data [L]
    char  (*SM)[MAXDCNAME];  ///< List of DC names in the system [L]
    char  (*SF)[MAXPHNAME+MAXSYMB];  ///< List of phase names in the system [FI]
    char *ICC;   ///< Classifier of IC { e o h a z v i <int> } [N]
    char *DCC;   ///< Classifier of DC { TESKWL GVCHNI JMFD QPR <0-9>  AB  XYZ O } [L]
    char *PHC;   ///< Classifier of phases { a g f p m l x d h } [FI]
    char *DCCW;  ///< internal DC class codes for calculation of primal chemical potentials [L]

     double  *MM;     //< DC molar masses, g/mol [L]
     double *FWGT;    //< phase (carrier) masses, g                [0:FI-1]
     double *YOF;     //< Surface free energy parameter for phases (J/g) (to accomodate for variable phase composition) [FI]

// Work part (to compute concentrations and chem.potentials), takeover from MULTI
    long int
    N,      // Number of ICs
        NR,       	///< NR - dimensions of R matrix
    L,      // Number of DCs
    Ls,     // Total number of DCs in phases-solutions
    LO,     // LO -   index of water-solvent in DC list
    FI,     // Number of phases
    FIs,    // Number of phases-solutions,
    K2,
    *L1;    //< l_a vector - number of DCs included into each phase [nPH]; copy of nDCinPH

    double
    *A,     ///< DC stoichiometry matrix A composed of a_ji [0:N-1][0:L-1]
    *Awt,   ///< IC atomic (molar) mass, g/mole [0:N-1]
    *XF,    ///< Output total number of moles of phases Xa[0:FI-1]
    *YF,    ///< Approximation of X_a in the next IPM iteration [0:FI-1]
    *XFA,   ///< Quantity of carrier in asymmetric phases Xwa, moles [FIs]
    *YFA,   ///< Approximation of XFA in the next IPM iteration [0:FIs-1]
    *X,     ///< DC quantities at eqstate x_j, moles - primal IPM solution [L]
    *Y,     ///< Copy of x_j from previous IPM iteration [0:L-1]
    *Fx,    ///< Dual DC chemical potentials defined via u_i and a_ji [L]
    *Wx,    ///< Mole fractions Wx of DC in multi-component phases [L]
    *F,     ///< Primal DC chemical potentials defined via g0_j, Wx_j and lnGam_j[L]
    *F0,    ///< Excess Gibbs energies for (metastable) DC, mole/mole [L]
    *U,     ///< IC chemical potentials u_i (mole/mole) - dual IPM solution [N]
    *Falp;  ///< Phase stability index (PC==2) [FI]

    double *GamFs;   ///< Copy of activity coefficients Gamma [L]
    double *fDQF,    ///< Increments to molar G0 values of DCs from pure gas fugacities or DQF terms, normalized [L]
    *MU,      ///< mu_j values of differences between dual and primal DC chem.potentials [L]
    *EMU,     ///< Exponents of DC increment to F_a criterion for phase [L]
    *NMU,     ///< DC increments to F_a criterion for phase [L]
    *Y_la,    ///< log activity of DC in phases (mju-mji0) [0:L-1]
    *Y_m,     ///< Molalities of aqueous species and sorbates [0:Ls-1]
    *Pparc,   ///< Partial pressures or fugacities of pure DC, bar (Pc by default) [0:L-1]
    *EZ,      ///< Formula charge of DC in multi-component phases [0:Ls-1]
    *Vol,     ///< DC molar volumes, cm3/mol [L]
    *FVOL,    ///< phase volumes, cm3 comment corrected DK 04.08.2009  [0:FI-1]
    *BF,    ///< Output bulk compositions of multicomponent phases bf_ai[FIs][N]
    *BFC;   ///< Total output bulk composition of all solid phases [1][N]

    double  ///< Numerical tolerances (scaled), taken over from MULTI
    XwMinM, ///< Cutoff mole amount for elimination of water-solvent { 1e-13 }
    ScMinM, ///< Cutoff mole amount for elimination of solid sorbent { 1e-13 }
    DcMinM, ///< Cutoff mole amount for elimination of solution- or surface species { 1e-30 }
    PhMinM, ///< Cutoff mole amount for elimination of non-electrolyte condensed phase { 1e-23 }
            ///< insertion values (re-scaled to system size)
//    DFYwM,  ///< Insertion mole amount for water-solvent { 1e-6 }
//    DFYaqM, ///< Insertion mole amount for aqueous and surface species { 1e-6 }
//    DFYidM, ///< Insertion mole amount for ideal solution components { 1e-6 }
//    DFYrM,  ///< Insertion mole amount for major solution components (incl. sorbent) { 1e-6 }
//    DFYhM,  ///< Insertion mole amount for minor solution components { 1e-6 }
//    DFYcM,  ///< Insertion mole amount for single-component phase { 1e-6 }
//    DFYsM,  ///< Insertion mole amount used in PhaseSelect() for a condensed phase component  { 1e-7 }
      SizeFactor, ///< factor for re-scaling the cutoffs/insertions to the system size
//    TMols,      ///< Input total moles in b vector before rescaling
//    SMols,      ///< Standart total moles (upscaled) {1000}
//    MBX,       ///< Total mass of the system, kg
//    FX,    	    ///< Current Gibbs potential of the system in IPM, moles
//    IC,         ///< Effective molal ionic strength of aqueous electrolyte
//    pH,         ///< pH of aqueous solution
//    pe,         ///< pe of aqueous solution
//    Eh,         ///< Eh of aqueous solution, V
//    DHBM,       ///< balance (relative) precision criterion
     DSM,        ///< min amount of phase DS
//    GWAT,       ///< used in ipm_gamma()
//    YMET,       ///< reserved
//    PCI,        ///< Current value of Dikin criterion of IPM convergence DK>=DX
//    DXM,        ///< IPM convergence criterion threshold DX (1e-5)
//    lnP,        ///< log Ptotal
//    RT,         ///< RT: 8.31451*T (J/mole/K)
//    FRT,        ///< F/RT, F - Faraday constant = 96485.309 C/mol
//    Yw,         ///< Current number of moles of solvent in aqueous phase
//    ln5551,     ///< ln(55.50837344)
    ICmin,      ///< tolerance for minimum ionic strength to calculate aqueous activity models
    aqsTail,    ///< v_j asymmetry correction factor for aqueous species
    lowPosNum,  ///< Minimum mole amount considered in GEM calculations (MinPhysAmount = 1.66e-24)
    logXw,      ///< work variable
    logYFk,     ///< work variable
    YFk;        ///< Current number of moles in a multicomponent phase
}
ACTIVITY;

// class TNode;

class TActivity
{
    ACTIVITY act;

    // Internal arrays for the performance optimization  (since version 2.0.0) <- MULTI
       long int sizeN; /*, sizeL, sizeAN;*/
       double *AA;
       double *BB;
       long int *arrL;
       long int *arrAN;

    long int sizeFIs;     ///< current size of phSolMod
    TSolMod* (*phSolMod); ///< size current FIs - number of multicomponent phases

    // new - allocation of TsorpMod
    long int sizeFIa;       ///< current size of phSorpMod
    TSorpMod* (*phSorpMod); ///< size current FIa - number of adsorption phases

    bool load;
    long int IT, IIM;  // current number of iterations of GEM algorithm (used in smoothing)
    double AG, DGC,    // smoothing
    FitVar[5];  ///< Internal. FitVar[0] is total mass (g) of solids in the system (sum over the BFC array)
                ///<      FitVar[1], [2] reserved
                ///<       FitVar[4] is the AG smoothing parameter;
                ///<       FitVar[3] is the actual smoothing coefficient
protected:

    void set_def( );
    void Alloc_TSolMod( long int newFIs );
    void Free_TSolMod( );
    void Alloc_TSorpMod( long int newFIs );
    void Free_TSorpMod();

    void Alloc_A_B( long int newN );
    void Free_A_B();
    void Build_compressed_xAN();
    void Free_compressed_xAN();
    void Free_internal();

public:
    void Alloc_internal();

    long int CheckMassBalanceResiduals(double *Y );
    double ConvertGj_toUniformStandardState( double g0, long int j, long int k );
    double PhaseSpecificGamma( long int j, long int jb, long int je, long int k, long int DirFlag = 0L ); // Added 26.06.08
    void SetSmoothingFactor( long int mode ); // new smoothing function (3 variants)

    // ipm_chemical.cpp
//        void XmaxSAT_IPM2();
    //    void XmaxSAT_IPM2_reset();
        double DC_DualChemicalPotential( double U[], double AL[], long int N, long int j );
//        void Set_DC_limits( long int Mode );
        void TotalPhasesAmounts( double X[], double XF[], double XFA[] );
        double DC_PrimalChemicalPotentialUpdate( long int j, long int k );
        double  DC_PrimalChemicalPotential( double G,  double logY,  double logYF,
                               double asTail,  double logYw,  char DCCW );
        void PrimalChemicalPotentials( double F[], double Y[],
                                      double YF[], double YFA[] );
        double KarpovCriterionDC( double *dNuG, double logYF, double asTail,
                     double logYw, double Wx,  char DCCW );
        void KarpovsPhaseStabilityCriteria();
        void  StabilityIndexes( );   // added 01.05.2010 DK
        double DC_GibbsEnergyContribution(   double G,  double x,  double logXF,
                                 double logXw,  char DCCW );
        double GX( double LM  );

        void ConvertDCC();
        long int  getXvolume();

    // ipm_chemical2.cpp
        void GasParcP();
        void phase_bcs( long int N, long int M, long int jb, double *A, double X[], double BF[] );
        void phase_bfc( long int k, long int jj );
        double bfc_mass( void );
        void CalculateConcentrationsInPhase( double X[], double XF[], double XFA[],
             double Factor, double MMC, double Dsur, long int jb, long int je, long int k );
        void CalculateConcentrations( double X[], double XF[], double XFA[]);
//        long int GouyChapman(  long int jb, long int je, long int k );

//  Surface activity coefficient terms
//        long int SurfaceActivityCoeff( long int jb, long int je, long int jpb, long int jdb, long int k );
//    void SurfaceActivityTerm( long int jb, long int je, long int k );  // Obsolete / deleted

//  ipm_chemical3.cpp
//        void IS_EtaCalc();
//        void pm_GC_ods_link( long int k, long int jb, long int jpb, long int jdb, long int ipb );
//        double SmoothingFactor( );
//        void SetSmoothingFactor( long int mode ); // new smoothing function (3 variants)
    // Main call for calculation of activity coefficients on IPM iterations
        long int CalculateActivityCoefficients( long int LinkMode );
    // Built-in activity coefficient models
    // Generic solution model calls
        void SolModCreate( long int jb, long int jmb, long int jsb, long int jpb, long int jdb,
                           long int k, long int ipb, char ModCode, char MixCode,
                           /* long int jphl, long int jlphc, */ long int jdqfc, long int jrcpc );
        void SolModParPT( long int k, char ModCode );
        void SolModActCoeff( long int k, char ModCode );
        void SolModExcessProp( long int k, char ModCode );
        void SolModIdealProp ( long int jb, long int k, char ModCode );
        void SolModStandProp ( long int jb, long int k, char ModCode );
        void SolModDarkenProp ( long int jb, long int k, char ModCode );

        /// This allocation is used only in standalone GEMS3K
        TActivity( DATACH *csd, DATABR *sbc )
        {
             sizeN = 0;
             AA = 0;
             BB = 0;
             arrL = 0;
             arrAN = 0;

             sizeFIs = 0;
             phSolMod = 0;
             sizeFIa = 0;
             phSorpMod = 0;
//             ICmin = 0.0001;
//             sizeFI = 0;
//             phKinMet = 0;
             act.N = csd->nIC;      // Number of ICs
                 act.NR = act.N;       	///< NR - dimensions of R matrix
             act.L = csd->nDC;      // Number of DCs
             act.Ls = csd->nDCs;     // Total number of DCs in phases-solutions
//             act.LO;     // LO -   index of water-solvent in DC list
             act.FI = csd->nPH;     // Number of phases
             act.FIs = csd->nPS;    // Number of phases-solutions,
             set_def();
             act.lnAct = new double[act.L];
             act.tpp_G = new double[act.L];
             load = false;
             Alloc_TSolMod( csd->nPS );
             sizeFIs = csd->nPS;
//             Alloc_TSorpMod( na->CSD.nPS );
//             sizeFIa = na->CSD.nPS;
        }

    ~TActivity( )
    {          // destructor
         if(act.lnAct) delete[] act.lnAct;
         if(act.tpp_G) delete[] act.tpp_G;
         Free_TSolMod( );
         Free_internal( );
//        Free_TSorpMod( );
    }

    // Returns pointer to ACTIVITY work data structure
    ACTIVITY* GetActivityDataPtr( void )
    {
        return &this->act;
    }
    // Generic access methods
    void setTemperature(double T); // set temperature (in units of K)
    void setPressure(double P); // set pressure (in units of Pa)
    void updateStandardGibbsEnergies(); // compute standard Gibbs energies (so far only P,T interpolation)
    void updateStandardVolumes();
    void updateStandardEnthalpies();
    void updateStandardEntropies();
//    void updateThermoData();

    void setSpeciation(const double* n); // set speciation (in units of moles)
    void updateConcentrations(); // compute concentrations in all phases
    void updateActivityCoefficients(); // compute activity coefficients
    void updateChemicalPotentials(); // compute primal chemical potentials
    void updateActivities(); // compute primal activities
    void updateChemicalData();
};


#endif

// -----------------------------------------------------------------------------
// end of _Activity_h

