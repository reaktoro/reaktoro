//-------------------------------------------------------------------
// $Id: ms_multi.h 986 2014-08-31 16:06:28Z kulik $
//
/// \file ms_multi.h
/// Declaration of TMulti class, configuration, and related functions
/// based on the IPM work data structure MULTI
//
/// \struct MULTI ms_multi.h
/// Contains chemical thermodynamic work data for GEM IPM-3 algorithm
//
// Copyright (c) 1995-2013 S.Dmytriyeva, D.Kulik, T.Wagner
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
#ifndef _ms_multi_h_
#define _ms_multi_h_

#ifndef IPMGEMPLUGIN

#include "m_param.h"
#include "v_ipnc.h"

// Internal subroutine for ET_translate() to process Phase scripts
typedef int (tget_ndx)( int nI, int nO, int Xplace );

#else

#include <ctime>
#include "m_const.h"
#include "datach.h"

#endif
// TSolMod header
#include "s_solmod.h"
// new: TsorpMod and TKinMet
#include "s_sorpmod.h"
#include "s_kinmet.h"

typedef struct
{  // MULTI is base structure to Project (local values)
  char
    stkey[EQ_RKLEN+1],   ///< Record key identifying IPM minimization problem
    // NV_[MAXNV], nulch, nulch1, ///< Variant Nr for fixed b,P,T,V; index in a megasystem
    PunE,         ///< Units of energy  { j;  J c C N reserved }
    PunV,         ///< Units of volume  { j;  c L a reserved }
    PunP,         ///< Units of pressure  { b;  B p P A reserved }
    PunT;         ///< Units of temperature  { C; K F reserved }
  long int
    N,        	///< N - number of IC in IPM problem
    NR,       	///< NR - dimensions of R matrix
    L,        	///< L -   number of DC in IPM problem
    Ls,       	///< Ls -   total number of DC in multi-component phases
    LO,       	///< LO -   index of water-solvent in IPM DC list
    PG,       	///< PG -   number of DC in gas phase
    PSOL,     	///< PSOL - number of DC in liquid hydrocarbon phase
    Lads,     	///< Total number of DC in sorption phases included into this system.
    FI,       	///< FI -   number of phases in IPM problem
    FIs,      	///< FIs -   number of multicomponent phases
    FIa,      	///< FIa -   number of sorption phases
    FI1,     ///< FI1 -   number of phases present in eqstate
    FI1s,    ///< FI1s -   number of multicomponent phases present in eqstate
    FI1a,    ///< FI1a -   number of sorption phases present in eqstate
    IT,      ///< It - number of completed IPM iterations
    E,       ///< PE - flag of electroneutrality constraint { 0 1 }
    PD,      ///< PD - mode of calling CalculateActivityCoefficients() { 0 1 2 3 4 }
    PV,      ///< Flag for the volume balance constraint (on Vol IC) - for indifferent equilibria at P_Sat { 0 1 }
    PLIM,    ///< PU - flag of activation of DC/phase restrictions { 0 1 }
    Ec,     ///< CalculateActivityCoefficients() return code: 0 (OK) or 1 (error)
    K2,     ///< Number of IPM loops performed ( >1 up to 6 because of PSSC() )
    PZ,     ///< Indicator of PSSC() status (since r1594): 0 untouched, 1 phase(s) inserted
            ///< 2 insertion done after 5 major IPM loops
    pNP,    ///< Mode of FIA selection: 0-automatic-LP AIA, 1-smart SIA, -1-user's choice
    pESU,   ///< Unpack old eqstate from EQSTAT record?  0-no 1-yes
    pIPN,   ///< State of IPN-arrays:  0-create; 1-available; -1 remake
    pBAL,   ///< State of reloading CSD:  1- BAL only; 0-whole CSD
    tMin,   ///< Type of thermodynamic potential to minimize
    pTPD,   ///< State of reloading thermod data: 0-all  -1-full from database   1-new system 2-no
    pULR,   ///< Start recalc kinetic constraints (0-do not, 1-do )internal
pKMM, ///< new: State of KinMet arrays: 0-create; 1-available; -1 remake
    ITaia,  ///< Number of IPM iterations completed in AIA mode (renamed from pRR1)
    FIat,   ///< max. number of surface site types
    MK,     ///< IPM return code: 0 - continue;  1 - converged
    W1,     ///< Indicator ofSpeciationCleanup() status (since r1594) 0 untouched, -1 phase(s) removed, 1 some DCs inserted
    is,     ///< is - index of IC for IPN equations ( CalculateActivityCoefficients() )
    js,     ///< js - index of DC for IPN equations ( CalculateActivityCoefficients() )
    next,   ///< for IPN equations (is it really necessary? TW please check!
    sitNcat,    //< Can be re-used
    sitNan      // Can be re-used
    ;
  double
    TC,  	///< Temperature T, min. (0,2000 C)
    TCc, 	///< Temperature T, max. (0,2000 C)
    T,   	///< T, min. K
    Tc,   	///< T, max. K
    P,      ///< Pressure P, min(0,10000 bar)
    Pc,   	///< Pressure P, max.(0,10000 bar)
    VX_,    ///< V(X) - volume of the system, min., cm3
    VXc,    ///< V(X) - volume of the system, max., cm3
    GX_,    ///< Gibbs potential of the system G(X), min. (J)
    GXc,    ///< Gibbs potential of the system G(X), max. (J)
    AX_,    ///< Helmholtz potential of the system F(X)
    AXc,    ///<  reserved
    UX_,  	///< Internal energy of the system U(X)
    UXc,  	///<  reserved
    HX_,    ///< Total enthalpy of the system H(X)
    HXc, 	///<  reserved
    SX_,    ///< Total entropy of the system S(X)
    SXc,   ///<	 reserved
    CpX_,  ///< reserved
    CpXc,  ///< reserved
    CvX_,  ///< reserved
    CvXc,  ///< reserved
  // TKinMet stuff
  kTau,  ///< current time, s (kinetics)        -->> databr  Tm
  kdT,   ///< current time step, s (kinetics)   -->> databr  dt

    TMols,      ///< Input total moles in b vector before rescaling
    SMols,      ///< Standart total moles (upscaled) {1000}
    MBX,        ///< Total mass of the system, kg
    FX,    	    ///< Current Gibbs potential of the system in IPM, moles
    IC,         ///< Effective molal ionic strength of aqueous electrolyte
    pH,         ///< pH of aqueous solution
    pe,         ///< pe of aqueous solution
    Eh,         ///< Eh of aqueous solution, V
    DHBM,       ///< balance (relative) precision criterion
    DSM,        ///< min amount of phase DS
    GWAT,       ///< used in ipm_gamma()
    YMET,       ///< reserved
    PCI,        ///< Current value of Dikin criterion of IPM convergence DK>=DX
    DXM,        ///< IPM convergence criterion threshold DX (1e-5)
    lnP,        ///< log Ptotal
    RT,         ///< RT: 8.31451*T (J/mole/K)
    FRT,        ///< F/RT, F - Faraday constant = 96485.309 C/mol
    Yw,         ///< Current number of moles of solvent in aqueous phase
    ln5551,     ///< ln(55.50837344)
    aqsTail,    ///< v_j asymmetry correction factor for aqueous species
    lowPosNum,  ///< Minimum mole amount considered in GEM calculations (MinPhysAmount = 1.66e-24)
    logXw,      ///< work variable
    logYFk,     ///< work variable
    YFk,        ///< Current number of moles in a multicomponent phase
    FitVar[5];  ///< Internal. FitVar[0] is total mass (g) of solids in the system (sum over the BFC array)
                ///<      FitVar[1], [2] reserved
                ///<       FitVar[4] is the AG smoothing parameter;
                ///<       FitVar[3] is the actual smoothing coefficient
  long int
  *L1,    ///< l_a vector - number of DCs included into each phase [Fi]
// from here moved to --> datach.h
  // TSolMod stuff
  *LsMod, ///< Number of interaction parameters. Max parameter order (cols in IPx),
            ///< and number of coefficients per parameter in PMc table [3*FIs]
  *LsMdc, ///<  for multi-site models: [3*FIs] - number of nonid. params per component;
           /// number of sublattices nS; number of moieties nM
  *LsMdc2, ///<  new: [3*FIs] - number of DQF coeffs; reciprocal coeffs per end member;
          /// reserved
  *IPx,   ///< Collected indexation table for interaction parameters of non-ideal solutions
            ///< ->LsMod[k,0] x LsMod[k,1]   over FIs

  *LsPhl,  ///< new: Number of phase links; number of link parameters; [Fi][2]
  (*PhLin)[2];  ///< new: indexes of linked phases and link type codes (sum 2*LsPhl[k][0] over Fi)

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

  // TSorpMod stuff
  long int
  *LsESmo, ///< new: number of EIL model layers; EIL params per layer; CD coefs per DC; reserved  [Fis][4]
  *LsISmo, ///< new: number of surface sites; isotherm coeffs per site; isotherm coeffs per DC; max.denticity of DC [Fis][4]
  *xSMd;   ///< new: denticity of surface species per surface site (site allocation) (-> L1[k]*LsISmo[k][3]] )
  long int  (*SATX)[4]; ///< Setup of surface sites and species (will be applied separately within each sorption phase) [Lads]
             /// link indexes to surface type [XL_ST]; sorbent em [XL_EM]; surf.site [XL-SI] and EDL plane [XL_SP]
  // TKinMet stuff
  long int
  *LsKin,  ///< new: number of parallel reactions nPRk[k]; number of species in activity products nSkr[k];
           /// number of parameter coeffs in parallel reaction term nrpC[k]; number of parameters
           /// per species in activity products naptC[k]; nAscC number of parameter coefficients in As correction;
           /// nFaces[k] number of (separately considered) crystal faces or surface patches ( 1 to 4 ) [Fi][6]
  *LsUpt,  ///< new: number of uptake kinetics model parameters (coefficients) numpC[k]; reserved [Fis][2]

  *xSKrC,  ///< new: Collected array of aq/gas/sorption species indexes used in activity products (-> += LsKin[k][1])
  (*ocPRkC)[2], ///< new: Collected array of operation codes for kinetic parallel reaction terms (-> += LsKin[k][0])
                /// and indexes of faces (surface patches)
*xICuC;  ///< new: Collected array of IC species indexes used in partition (fractionation) coefficients  ->L1[k]   TBD

  double
   // TSolMod stuff
    *PMc,    ///< Collected interaction parameter coefficients for the (built-in) non-ideal mixing models -> LsMod[k,0] x LsMod[k,2]
    *DMc,    ///< Non-ideality coefficients f(TPX) for DC -> L1[k] x LsMdc[k][0]
    *MoiSN,  ///< End member moiety- site multiplicity number tables ->  L1[k] x LsMdc[k][1] x LsMdc[k][2]
    *SitFr;  ///< Tables of sublattice site fractions for moieties -> LsMdc[k][1] x LsMdc[k][2]

  double
    denW[5],   ///< Density of water, first T, second T, first P, second P derivative for Tc,Pc
    denWg[5],  ///< Density of steam for Tc,Pc
    epsW[5],   ///< Diel. constant of H2O(l)for Tc,Pc
    epsWg[5];  ///< Diel. constant of steam for Tc,Pc
// until here moved to --> datach.h

  // Stoichiometry basis
  double
    *A,   ///< DC stoichiometry matrix A composed of a_ji [0:N-1][0:L-1]
    *Awt,    ///< IC atomic (molar) mass, g/mole [0:N-1]

 // Reconsider usage
     *Wb,     ///< Relative Born factors (HKF, reserved) [0:Ls-1]
     *Wabs,   ///< Absolute Born factors (HKF, reserved) [0:Ls-1]
     *Rion,   ///< Ionic or solvation radii, A (reserved) [0:Ls-1]
     *HYM,    ///< reserved
     *ENT,    ///< reserved no object

     *H0,     ///< DC molar enthalpies, reserved [L]
     *A0,     ///< DC molar Helmholtz energies, reserved [L]
     *U0,     ///< DC molar internal energies, reserved [L]
     *S0,     ///< DC molar entropies, reserved [L]
     *Cp0,    ///< DC molar heat capacity, reserved [L]
     *Cv0,    ///< DC molar Cv, reserved [L]

    *VL,      ///< ln mole fraction of end members in phases-solutions
  // Old sorption stuff
    *Xcond,   ///< conductivity of phase carrier, sm/m2   [0:FI-1], reserved
    *Xeps,    ///< diel.permeability of phase carrier (solvent) [0:FI-1], reserved
    *Aalp,    ///< Full vector of specific surface areas of phases (m2/g) [0:FI-1]
    *Sigw,    ///< Specific surface free energy for phase-water interface (J/m2)   [0:FI-1]
    *Sigg,  	///< Specific surface free energy for phase-gas interface (J/m2) (not yet used)  [0:FI-1], reserved
// from here moved to --> datach.h
  // TSolMod stuff
  *lPhc,  ///< new: Collected array of phase link parameters (sum(LsPhl[k][1] over Fi)
  *DQFc,  ///< new: Collected array of DQF parameters for DCs in phases -> L1[k] x LsMdc2[k][0]
//  *rcpc,  ///< new: Collected array of reciprocal parameters for DCs in phases -> L1[k] x LsMdc2[k][1]

  // TSorpMod & TKinMet stuff
  *SorMc, ///< new: Phase-related kinetics and sorption model parameters: [Fis][16]
          ///< in the same order as from Asur until fRes2 in TPhase
  // TSorpMod stuff
  *EImc,  ///< new: Collected EIL model coefficients k -> += LsESmo[k][0]*LsESmo[k][1]
  *mCDc,  ///< new: Collected CD EIL model coefficients per DC k -> += L1[k]*LsESmo[k][2]
  *IsoPc, ///< new: Collected isotherm coefficients per DC k -> += L1[k]*LsISmo[k][2];
  *IsoSc, ///< new: Collected isotherm coeffs per site k -> += LsISmo[k][0]*LsISmo[k][1];
  // TKinMet stuff
  *feSArC, ///< new: Collected array of fractions of surface area related to parallel reactions k-> += LsKin[k][0]
  *rpConC,  ///< new: Collected array of kinetic rate constants k-> += LsKin[k][0]*LsKin[k][2];
  *apConC,  ///< new:!! Collected array of parameters per species involved in activity product terms
          ///  k-> += LsKin[k][0]*LsKin[k][1]*LsKin[k][3];
  *AscpC,   /// new: parameter coefficients of equation for correction of specific surface area k-> += LsKin[k][4]
  *UMpcC  ///< new: Collected array of uptake model coefficients k-> += L1[k]*LsUpt[k][0];
      ;
  // until here moved to --> datach.h

  //  Data for old surface comlexation and sorption models (new variant [Kulik,2006])
  double  (*Xr0h0)[2];   ///< mean r & h of particles (- pores), nm  [0:FI-1][2], reserved
  double  (*Nfsp)[MST];  ///< Fractions of the sorbent specific surface area allocated to surface types  [FIs][FIat]
  double  (*MASDT)[MST]; ///< Total maximum site  density per surface type (mkmol/g)  [FIs][FIat]
  double  (*XcapF)[MST]; ///< Capacitance density of Ba EDL layer F/m2 [FIs][FIat]
  double  (*XcapA)[MST]; ///< Capacitance density of 0 EDL layer, F/m2 [FIs][FIat]
  double  (*XcapB)[MST]; ///< Capacitance density of B EDL layer, F/m2 [FIs][FIat]
  double  (*XcapD)[MST]; ///< Eff. cap. density of diffuse layer, F/m2 [FIs][FIat]
  double  (*XdlA)[MST];  ///< Effective thickness of A EDL layer, nm [FIs][FIat], reserved
  double  (*XdlB)[MST];  ///< Effective thickness of B EDL layer, nm [FIs][FIat], reserved
  double  (*XdlD)[MST];  ///< Effective thickness of diffuse layer, nm [FIs][FIat], reserved
  double  (*XlamA)[MST]; ///< Factor of EDL discretness  A < 1 [FIs][FIat], reserved
  double  (*Xetaf)[MST]; ///< Density of permanent surface type charge (mkeq/m2) for each surface type on sorption phases [FIs][FIat]
  double  (*MASDJ)[DFCN];  ///< Parameters of surface species in surface complexation models [Lads][DFCN]
                           // Contents defined in the enum below this structure
// Other data
  double
    *XFs,    ///< Current quantities of phases X_a at IPM iterations [0:FI-1]
    *Falps,  ///< Current Karpov criteria of phase stability  F_a [0:FI-1]
    *Fug,    ///< Demo partial fugacities of gases [0:PG-1]
    *Fug_l,  ///< Demo log partial fugacities of gases [0:PG-1]
    *Ppg_l,  ///< Demo log partial pressures of gases [0:PG-1]

    *DUL,     ///< VG Vector of upper kinetic restrictions to x_j, moles [L]
    *DLL,     ///< NG Vector of lower kinetic restrictions to x_j, moles [L]
    *fDQF,    ///< Increments to molar G0 values of DCs from pure gas fugacities or DQF terms, normalized [L]
// TKinMet stuff (old DODs, new contents )
    *PUL,  ///< Vector of upper restrictions to multicomponent phases amounts [FIs]
    *PLL,  ///< Vector of lower restrictions to multicomponent phases amounts [FIs]
*PfFact, /// new: phase surface area - volume shape factor (taken from TKinMet or set from TNode) [FI]
*PrT,    /// new: Total MWR rate (mol/s) for phases - TKinMet output [FI]
*PkT,    /// new: Total specific MWR rate (mol/m2/s) for phases - TKinMet output [FI]
*PvT,    /// new: Total one-dimensional MWR surface propagation velocity (m/s) - TKinMet output [FI]
//  potentially can be extended to all solution phases?
*emRd,   /// new: output Rd values (partition coefficients) for end members (in uptake kinetics model) [Ls]
*emDf,   /// new: output Df values (fractionation coeffs.) for end members (in uptake kinetics model) [Ls]
//
    *YOF,     ///< Surface free energy parameter for phases (J/g) (to accomodate for variable phase composition) [FI]
    *Vol,     ///< DC molar volumes, cm3/mol [L]
    *MM,      ///< DC molar masses, g/mol [L]
    *Pparc,   ///< Partial pressures or fugacities of pure DC, bar (Pc by default) [0:L-1]
    *Y_m,     ///< Molalities of aqueous species and sorbates [0:Ls-1]
    *Y_la,    ///< log activity of DC in multi-component phases (mju-mji0) [0:L-1]
    *Y_w,     ///< Mass concentrations of DC in multi-component phases,%(ppm)[Ls]
    *Gamma,   ///< DC activity coefficients in molal or other phase-specific scale [0:L-1]
    *lnGmf,   ///< ln of initial DC activity coefficients for correcting G0 [0:L-1]
    *lnGmM,   ///< ln of DC pure gas fugacity (or metastability) coeffs or DDF correction [0:L-1]
    *EZ,      ///< Formula charge of DC in multi-component phases [0:Ls-1]
    *FVOL,    ///< phase volumes, cm3 comment corrected DK 04.08.2009  [0:FI-1]
    *FWGT,    ///< phase (carrier) masses, g                [0:FI-1]
//
    *G,       ///< Normalized DC energy function c(j), mole/mole [0:L-1]            --> activities.h
    *G0,      ///< Input normalized g0_j(T,P) for DC at unified standard scale[L]   --> activities.h
    *lnGam,   ///< ln of DC activity coefficients in unified (mole-fraction) scale [0:L-1] --> activities.h
    *lnGmo;   ///< Copy of lnGam from previous IPM iteration (reserved)
  double  (*lnSAC)[4]; ///< former lnSAT ln surface activity coeff and Coulomb's term  [Lads][4]

  // TSolMod stuff (detailed output on partial energies of mixing)   --> activities.h
  double *lnDQFt; ///< new: DQF terms adding to overall activity coefficients [Ls_]
  double *lnRcpt; ///< new: reciprocal terms adding to overall activity coefficients [Ls_]
  double *lnExet; ///< new: excess energy terms adding to overall activity coefficients [Ls_]
  double *lnCnft; ///< new: configurational terms adding to overall activity [Ls_]
  // TSorpMod stuff
  double *lnScalT;  ///< new: Surface/volume scaling activity correction terms [Ls_]
  double *lnSACT;   ///< new: ln isotherm-specific SACT for surface species [Ls_]
  double *lnGammF;  ///< new: Frumkin or BET non-electrostatic activity coefficients [Ls_]
  double *CTerms;   ///< new: Coulombic terms (electrostatic activity coefficients) [Ls_]

  double  *B,  ///< Input bulk chem. compos. of the system - b vector, moles of IC[N]
    *U,  ///< IC chemical potentials u_i (mole/mole) - dual IPM solution [N]
    *U_r,  ///< IC chemical potentials u_i (J/mole) [0:N-1]
    *C,    ///< Calculated IC mass-balance deviations (moles) [0:N-1]
    *IC_m, ///< Total IC molalities in aqueous phase (excl.solvent) [0:N-1]
    *IC_lm,	///< log total IC molalities in aqueous phase [0:N-1]
    *IC_wm,	///< Total dissolved IC concentrations in g/kg_soln [0:N-1]
    *BF,    ///< Output bulk compositions of multicomponent phases bf_ai[FIs][N]
    *BFC,   ///< Total output bulk composition of all solid phases [1][N]
    *XF,    ///< Output total number of moles of phases Xa[0:FI-1]
    *YF,    ///< Approximation of X_a in the next IPM iteration [0:FI-1]
    *XFA,   ///< Quantity of carrier in asymmetric phases Xwa, moles [FIs]
    *YFA,   ///< Approximation of XFA in the next IPM iteration [0:FIs-1]
    *Falp;  ///< Karpov phase stability criteria F_a [0:FI-1] or phase stability index (PC==2)

  double (*VPh)[MIXPHPROPS],     ///< Volume properties for mixed phases [FIs]
         (*GPh)[MIXPHPROPS],     ///< Gibbs energy properties for mixed phases [FIs]
         (*HPh)[MIXPHPROPS],     ///< Enthalpy properties for mixed phases [FIs]
         (*SPh)[MIXPHPROPS],     ///< Entropy properties for mixed phases [FIs]
         (*CPh)[MIXPHPROPS],     ///< Heat capacity Cp properties for mixed phases [FIs]
         (*APh)[MIXPHPROPS],     ///< Helmholtz energy properties for mixed phases [FIs]
         (*UPh)[MIXPHPROPS];     ///< Internal energy properties for mixed phases [FIs]

// old sorption models - EDL models (data for electrostatic activity coefficients)
   double (*XetaA)[MST]; ///< Total EDL charge on A (0) EDL plane, moles [FIs][FIat]
   double (*XetaB)[MST]; ///< Total charge of surface species on B (1) EDL plane, moles[FIs][FIat]
   double (*XetaD)[MST]; ///< Total charge of surface species on D (2) EDL plane, moles[FIs][FIat]
   double (*XpsiA)[MST]; ///< Relative potential at A (0) EDL plane,V [FIs][FIat]
   double (*XpsiB)[MST]; ///< Relative potential at B (1) EDL plane,V [FIs][FIat]
   double (*XpsiD)[MST]; ///< Relative potential at D (2) plane,V [FIs][FIat]
   double (*XFTS)[MST];  ///< Total number of moles of surface DC at surface type [FIs][FIat]
//
   double *X,  ///< DC quantities at eqstate x_j, moles - primal IPM solution [L]
    *Y,   ///< Copy of x_j from previous IPM iteration [0:L-1]
    *XY,  ///< Copy of x_j from previous loop of Selekt2() [0:L-1]
    *Qp,  ///< Work variables related to non-ideal phases FIs*(QPSIZE=180)
    *Qd,  ///< Work variables related to DC in non-ideal phases FIs*(QDSIZE=60)
    *MU,  ///< mu_j values of differences between dual and primal DC chem.potentials [L]
    *EMU, ///< Exponents of DC increment to F_a criterion for phase [L]
    *NMU, ///< DC increments to F_a criterion for phase [L]
    *W,   ///< Weight multipliers for DC (incl restrictions) in IPM [L]
    *Fx,  ///< Dual DC chemical potentials defined via u_i and a_ji [L]
    *Wx,  ///< Mole fractions Wx of DC in multi-component phases [L]
    *F,   /// <Primal DC chemical potentials defined via g0_j, Wx_j and lnGam_j[L]
    *F0;  ///< Excess Gibbs energies for (metastable) DC, mole/mole [L]
// Old sorption models
   double (*D)[MST];  ///< Reserved; new work array for calc. surface act.coeff.
// Name lists
  char (*sMod)[8];   ///< new: Codes for built-in mixing models of multicomponent phases [FIs]
  char (*kMod)[6];  ///< new: Codes for built-in kinetic models [Fi]
  char  (*dcMod)[6];   ///< Codes for PT corrections for dependent component data [L]
  char  (*SB)[MAXICNAME+MAXSYMB]; ///< List of IC names in the system [N]
  char  (*SB1)[MAXICNAME]; ///< List of IC names in the system [N]
  char  (*SM)[MAXDCNAME];  ///< List of DC names in the system [L]
  char  (*SF)[MAXPHNAME+MAXSYMB];  ///< List of phase names in the system [FI]
  char  (*SM2)[MAXDCNAME];  ///< List of multicomp. phase DC names in the system [Ls]
  char  (*SM3)[MAXDCNAME];  ///< List of adsorption DC names in the system [Lads]
  char  *DCC3;   ///< Classifier of DCs involved in sorption phases [Lads]
  char  (*SF2)[MAXPHNAME+MAXSYMB]; ///< List of multicomp. phase names in the syst [FIs]
  char  (*SFs)[MAXPHNAME+MAXSYMB]; ///< List of phases currently present in non-zero quantities [FI]
  char  *pbuf, 	///< Text buffer for table printouts
// Class codes
    *RLC,   ///< Code of metastability constraints for DCs [L] enum DC_LIMITS
    *RSC,   ///< Units of metastability/kinetic constraints for DCs  [L]
    *RFLC,  ///< Classifier of restriction types for XF_a [FIs]
    *RFSC,  ///< Classifier of restriction scales for XF_a [FIs]
    *ICC,   ///< Classifier of IC { e o h a z v i <int> } [N]
    *DCC,   ///< Classifier of DC { TESKWL GVCHNI JMFD QPR <0-9>  AB  XYZ O } [L]
    *PHC;   ///< Classifier of phases { a g f p m l x d h } [FI]
  char  (*SCM)[MST]; ///< Classifier of built-in electrostatic models applied to surface types in sorption phases [FIs][FIat]
  char  *SATT,  ///< Classifier of applied SACT equations (isotherm corrections) [Lads]
    *DCCW;  ///< internal DC class codes [L]
  // TSorpMod stuff
  char *IsoCt; ///< new: Collected isotherm and SATC codes for surface site types k -> += 2*LsISmo[k][0]

//  SolutionData *asd; ///< Array of data structures to pass info to TSolMod [FIs]

  long int ITF,       ///< Number of completed IA EFD iterations
         ITG,         ///< Number of completed GEM IPM iterations
         ITau,    /// new: Time iteration for TKinMet class calculations
         IRes1;
  clock_t t_start, t_end;
  double t_elap_sec;  ///< work variables for determining IPM calculation time
#ifdef IPMGEMPLUGIN
  double *Guns;     ///<  mu.L work vector of uncertainty space increments to tp->G + sy->GEX
  double *Vuns;     ///<  mu.L work vector of uncertainty space increments to tp->Vm
  double *tpp_G;    ///< Partial molar(molal) Gibbs energy g(TP) (always), J/mole
  double *tpp_S;    ///< Partial molar(molal) entropy s(TP), J/mole/K
  double *tpp_Vm;   ///< Partial molar(molal) volume Vm(TP) (always), J/bar
#endif

  // additional arrays for internal calculation in ipm_main
  long int
  *mui,   ///< IC indices in RMULTS IC list [N]
  *muk,   ///< Phase indices in RMULTS phase list [FI]
  *muj;   ///< DC indices in RMULTS DC list [L]
  double *XU;      ///< dual-thermo calculation of DC amount X(j) from A matrix and u vector [L]
  double (*Uc)[2]; ///< Internal copy of IC chemical potentials u_i (mole/mole) at r-1 and r-2 [N][2]
  double *Uefd;    ///< Internal copy of IC chemical potentials u_i (mole/mole) - EFD function [N]
  char errorCode[100]; ///<  code of error in IPM      (Ec number of error)
  char errorBuf[1024]; ///< description of error in IPM
  double logCDvalues[5]; ///< Collection of lg Dikin crit. values for the new smoothing equation
double *GamFs;   ///< Copy of activity coefficients Gamma before the first enter in PhaseSelection() [L] new

  double // Iterators for MTP interpolation (do not load/unload for IPM)
  Pai[4],    ///< Pressure P, bar: start, end, increment for MTP array in DataCH , Ptol
  Tai[4],    ///< Temperature T, C: start, end, increment for MTP array in DataCH , Ttol
  Fdev1[2],  ///< Function1 and target deviations for  minimization of thermodynamic potentials
  Fdev2[2];  ///< Function2 and target deviations for  minimization of thermodynamic potentials

  // Experimental: modified cutoff and insertion values (DK 28.04.2010)
  double
// cutoffs (rescaled to system size)
  XwMinM, ///< Cutoff mole amount for elimination of water-solvent { 1e-13 }
  ScMinM, ///< Cutoff mole amount for elimination of solid sorbent { 1e-13 }
  DcMinM, ///< Cutoff mole amount for elimination of solution- or surface species { 1e-30 }
  PhMinM, ///< Cutoff mole amount for elimination of non-electrolyte condensed phase { 1e-23 }
          ///< insertion values (re-scaled to system size)
  DFYwM,  ///< Insertion mole amount for water-solvent { 1e-6 }
  DFYaqM, ///< Insertion mole amount for aqueous and surface species { 1e-6 }
  DFYidM, ///< Insertion mole amount for ideal solution components { 1e-6 }
  DFYrM,  ///< Insertion mole amount for major solution components (incl. sorbent) { 1e-6 }
  DFYhM,  ///< Insertion mole amount for minor solution components { 1e-6 }
  DFYcM,  ///< Insertion mole amount for single-component phase { 1e-6 }
  DFYsM,  ///< Insertion mole amount used in PhaseSelect() for a condensed phase component  { 1e-7 }
  SizeFactor; ///< factor for re-scaling the cutoffs/insertions to the system size
}
MULTI;

/// Indexation in a row of the pmp->SATX[][] array.
/// [0] - max site density in mkmol/(g sorbent);
/// [1] - species charge allocated to 0 plane;
/// [2] - surface species charge allocated to beta -or third plane;
/// [3] - Frumkin interaction parameter;
/// [4] species denticity or coordination number;
/// [5]  - reserved parameter (e.g. species charge on 3rd EIL plane)
enum IndexationSATX {
   XL_ST = 0, XL_EM = 1, XL_SI = 2, XL_SP = 3
};

struct SPP_SETTING;
class TProfil;
class TNode;

// Data of MULTI
class TMulti
#ifndef IPMGEMPLUGIN
  : public TSubModule
#endif
{
    MULTI pm;
    MULTI *pmp;
    SPP_SETTING *paTProfil;

// Internal arrays for the performance optimization  (since version 2.0.0)
   long int sizeN; /*, sizeL, sizeAN;*/
   double *AA;
   double *BB;
   long int *arrL;
   long int *arrAN;

   void Alloc_A_B( long int newN );
   void Free_A_B();
   void Build_compressed_xAN();
   void Free_compressed_xAN();
   void Free_internal();

#ifndef IPMGEMPLUGIN
// These pointers and methods are only used in GEMS-PSI
    //SYSTEM *syp;
    //MTPARM *tpp;
    //RMULTS *mup;

    void MultiSystemInit();
    void SystemToLookup();

    void multi_sys_dc();
    void multi_sys_ph();
    void ph_sur_param( int k, int kk );
    void ph_surtype_assign( int k, int kk, int jb, int je,
                            int car_l[], int car_c, int Cjs );
    void sm_text_analyze( int nph, int Type, int JB, int JE, int jb, int je );
    void SolModLoad();
// new TKinMet stuff
void KinMetModLoad();
    bool CompressPhaseIpxt( int kPH );
    gstring PressSolMod( int nP );
    char *ExtractEG( char *Etext, int jp, int *EGlen, int Nes );
    int find_icnum( char *name, int LNmode );
    int find_dcnum( char *name, int jb, int je, int LNmode, char *stmt  );
    int find_phnum( char *name, int LNmode );
    int find_acnum( char *name, int LNmode );
    int find_phnum_multi( const char *name);
    int find_dcnum_multi( const char *name);
    int find_icnum_multi( const char *name);
    const char* GetHtml();

#else

   bool load; // internal value

#endif

protected:
// From here move to activities.h or node.h
   long int sizeFIs;     ///< current size of phSolMod
   TSolMod* (*phSolMod); ///< size current FIs - number of multicomponent phases

   void Alloc_TSolMod( long int newFIs );
   void Free_TSolMod();

   // new - allocation of TsorpMod and TKinMet
   long int sizeFIa;       ///< current size of phSorpMod
   TSorpMod* (*phSorpMod); ///< size current FIa - number of adsorption phases

   void Alloc_TSorpMod( long int newFIa );
   void Free_TSorpMod();

   long int sizeFI;      ///< current size of phKinMet
   TKinMet* (*phKinMet); ///< size current FI -   number of phases

   void Alloc_TKinMet( long int newFI );
   void Free_TKinMet();
// until here move to activities.h or node.h

// Added for implementation of divergence detection in dual solution 06.05.2011 DK
   long int nNu;  ///< number of ICs in the system
   long int cnr;  ///< current IPM iteration
   long int nCNud; ///< number of IC names for divergent dual chemical potentials
   double *U_mean; ///< Cumulative mean dual solution approximation [nNu]
   double *U_M2;   ///< Cumulative sum of squares [nNu]
   double *U_CVo;  ///< Cumulative Coefficient of Variation for dual solution approximation r-1 [nNu]
   double *U_CV;   ///< Cumulative Coefficient of Variation for r-th dual solution approximation [nNu]
   long int *ICNud; ///< List of IC indexes for divergent dual chemical potentials [nNu]
   void Alloc_uDD( long int newN );
   void Free_uDD();
   void Reset_uDD(long int cr, bool trace = false );
   void Increment_uDD( long int r, bool trace = false );
   long int Check_uDD( long int mode, double DivTol, bool trace = false );

#ifdef IPMGEMPLUGIN

   char PAalp_; ///< Flag for using (+) or ignoring (-) specific surface areas of phases
   char PSigm_; ///< Flag for using (+) or ignoring (-) specific surface free energies

#endif

   void setErrorMessage( long int num, const char *code, const char * msg);
   void addErrorMessage( const char * msg);

// ipm_chemical.cpp
    void XmaxSAT_IPM2();
//    void XmaxSAT_IPM2_reset();
    double DC_DualChemicalPotential( double U[], double AL[], long int N, long int j );
    void Set_DC_limits( long int Mode );
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
    long int GouyChapman(  long int jb, long int je, long int k );

//  Surface activity coefficient terms
    long int SurfaceActivityCoeff( long int jb, long int je, long int jpb, long int jdb, long int k );
//    void SurfaceActivityTerm( long int jb, long int je, long int k );  // Obsolete / deleted

// ipm_chemical3.cpp
    void IS_EtaCalc();
    void pm_GC_ods_link( long int k, long int jb, long int jpb, long int jdb, long int ipb );
    double SmoothingFactor( );
    void SetSmoothingFactor( long int mode ); // new smoothing function (3 variants)
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

    // Specific phase property calculation functions  // obsolete (29.11.10 TW)
    // void IdealGas( long int jb, long int k, double *Zid );
    // void IdealOneSite( long int jb, long int k, double *Zid );
    // void IdealMultiSite( long int jb, long int k, double *Zid );

// New stuff for TKinMet class implementation
long int CalculateKinMet( long int LinkMode  );
void KM_Create(long int jb, long int k, long int kc, long int kp, long int kf,
            long int ka, long int ks, long int kd, long ku, long ki, const char *kmod,
            long jphl, long jlphc );
void KM_ParPT( long int k, const char *kMod );
void KM_InitTime( long int k, const char *kMod );
void KM_UpdateTime( long int k, const char *kMod );
void KM_UpdateFSA(long jb, long int k, const char *kMod );
void KM_ReturnFSA(long int k, const char *kMod );
void KM_CalcRates( long int k, const char *kMod );
void KM_InitRates( long int k, const char *kMod );
void KM_CalcSplit( long int jb, long int k, const char *kMod );
void KM_InitSplit( long int jb, long int k, const char *kMod );
void KM_CalcUptake( long int jb, long int k, const char *kMod );
void KM_InitUptake( long int jb, long int k, const char *kMod );
void KM_SetAMRs( long int jb, long int k, const char *kMod );

// ipm_main.cpp - numerical part of GEM IPM-2
    void GEM_IPM( long int rLoop );
    long int MassBalanceRefinement( long int WhereCalledFrom );
    long int InteriorPointsMethod( long int &status, long int rLoop );
    void AutoInitialApproximation( );

// ipm_main.cpp - miscellaneous fuctions of GEM IPM-2
   void MassBalanceResiduals( long int N, long int L, double *A, double *Y,
                               double *B, double *C );
//   long int CheckMassBalanceResiduals(double *Y );
   double OptimizeStepSize( double LM );
   void DC_ZeroOff( long int jStart, long int jEnd, long int k=-1L );
   void DC_RaiseZeroedOff( long int jStart, long int jEnd, long int k=-1L );
   double RaiseDC_Value( const long int j );
   //   void LagrangeMultiplier();
   long int MetastabilityLagrangeMultiplier();
   void WeightMultipliers( bool square );
   long int MakeAndSolveSystemOfLinearEquations( long int N, bool initAppr );
   double DikinsCriterion(  long int N, bool initAppr );
   double StepSizeEstimate( bool initAppr );
   void Restore_Y_YF_Vectors();

   double RescaleToSize( bool standard_size ); // replaced calcSfactor() 30.08.2009 DK
   long int SpeciationCleanup( double AmountThreshold, double ChemPotDiffCutoff ); // added 25.03.10 DK
   long int PhaseSelectionSpeciationCleanup( long int &k_miss, long int &k_unst, long int rLoop );
   long int PhaseSelect( long int &k_miss, long int &k_unst, long int rLoop );
   bool GEM_IPM_InitialApproximation();

   // IPM_SIMPLEX.CPP Simplex method modified with two-sided constraints (Karpov ea 1997)
    void SolveSimplex(long int M, long int N, long int T, double GZ, double EPS,
                 double *UND, double *UP, double *B, double *U,
                 double *AA, long int *STR, long int *NMB );
    void SPOS( double *P, long int STR[],long int NMB[],long int J,long int M,double AA[]);
    void START( long int T,long int *ITER,long int M,long int N,long int NMB[],
                double GZ,double EPS,long int STR[],long int *BASE,
                double B[],double UND[],double UP[],double AA[],double *A,
                double *Q );
    void NEW(long int *OPT,long int N,long int M,double EPS,double *LEVEL,long int *J0,
             long int *Z,long int STR[], long int NMB[], double UP[],
             double AA[], double *A);
    void WORK(double GZ,double EPS,long int *I0, long int *J0,long int *Z,long int *ITER,
              long int M, long int STR[],long int NMB[],double AA[],
              long int BASE[],long int *UNO,double UP[],double *A,double Q[]);
    void FIN(double EPS,long int M,long int N,long int STR[],long int NMB[],
             long int BASE[],double UND[],double UP[],double U[],
             double AA[],double *A,double Q[],long int *ITER);
    void GibbsEnergyMinimization();
    double SystemTotalMolesIC( );
    void ScaleSystemToInternal(  double ScFact );
    void RescaleSystemFromInternal(  double ScFact );
    void MultiConstInit(); // from MultiRemake
    void GEM_IPM_Init();

public:
    TNode *node;

    void set_def( long int i=0);

#ifndef IPMGEMPLUGIN

    static TMulti* sm;

// This is used only in GEM-Selektor
    TIArray<IPNCalc> qEp;
    TIArray<IPNCalc> qEd;

    TMulti( int nrt );
    ~TMulti()
    {
//      Free_TSolMod();     // Added 06.05.2011 DK - possible bugfix
       Free_internal();
       Free_uDD();
    };

    void ods_link( int i=0);
    void dyn_set( int i=0);
    void dyn_kill( int i=0);
    void dyn_new( int i=0);
//    void set_def( int i=0);
//    void sit_dyn_new();

    // ms_muleq.cpp
    void packData();
    void packData( TCIntArray PHon, TCIntArray DCon );
    void setSizes();
    void loadData( bool newRec );
    void unpackData();

    void MultiKeyInit( const char*key );
    void EqstatExpand( const char *key,  bool calcActivityModels, bool calcKineticModels );
    void ET_translate( int nOet, int nOpex, int JB, int JE, int jb, int je,
     tget_ndx *get_ndx = 0 );
    void getNamesList( int nO, TCStringArray& lst );

   class UserCancelException {};

   /// connection to UnSpace
   double pb_GX( double *Gxx  );

#else
/// This allocation is used only in standalone GEMS3K
   TMulti( TNode* na_ )
   {
     pmp = &pm;
     node = na_; // parent
     sizeN = 0;
     AA = 0;
     BB = 0;
     arrL = 0;
     arrAN = 0;

 U_mean = 0;
 U_M2 = 0;
 U_CVo = 0;
 U_CV = 0;
 ICNud = 0;

     sizeFIs = 0;
     phSolMod = 0;
     sizeFIa = 0;
     phSorpMod = 0;
     sizeFI = 0;
     phKinMet = 0;

     pmp->Guns = 0;
     pmp->Vuns = 0;
     pmp->tpp_G = 0;
     pmp->tpp_S = 0;
     pmp->tpp_Vm = 0;

     load = false;
   }

    ~TMulti()
    {  multi_free(); }

    void multi_realloc( char PAalp, char PSigm );
    void multi_free();

    void CheckMtparam(); // Test load thermodynamic data before

    void set_load (bool what); // DM 20.05.2013 - Ensures the re-reading of the system properties into GEM IMP data structure

#endif

    MULTI* GetPM()
    { return &pm; }

    void setPa( TProfil *prof);
    long int testMulti( );

    const char* GetName() const
    {  return "Multi";  }

    // Internal functions for SCMs
//       void getLsModsum( long int& LsModSum, long int& LsIPxSum );
//       void getLsMdcsum( long int& LsMdcSum,long int& LsMsnSum,long int& LsSitSum );

    //connection to mass transport
    void to_file( GemDataStream& ff );
    void to_text_file( const char *path, bool append=false  );
    void from_file( GemDataStream& ff );
    void to_text_file_gemipm( const char *path, bool addMui,
    		bool with_comments = true, bool brief_mode = false );
    void from_text_file_gemipm( const char *path,  DATACH  *dCH );

    // New functions for TSolMod, TKinMet and TSorpMod parameter arrays
       void getLsModsum( long int& LsModSum, long int& LsIPxSum );
       void getLsMdcsum( long int& LsMdcSum,long int& LsMsnSum,long int& LsSitSum );
       /// Get dimensions from LsPhl array
       void getLsPhlsum( long int& PhLinSum,long int& lPhcSum );
       /// Get dimensions from LsMdc2 array
       void getLsMdc2sum( long int& DQFcSum,long int& rcpcSum );
       /// Get dimensions from LsISmo array
       void getLsISmosum( long int& IsoCtSum,long int& IsoScSum, long int& IsoPcSum,long int& xSMdSum );
       /// Get dimensions from LsESmo array
       void getLsESmosum( long int& EImcSum,long int& mCDcSum );
       /// Get dimensions from LsKin array
       void getLsKinsum( long int& xSKrCSum,long int& ocPRkC_feSArC_Sum,
                         long int& rpConCSum,long int& apConCSum, long int& AscpCSum );
       /// Get dimensions from LsUpot array
       void getLsUptsum( long int& UMpcSum );

    // EXTERNAL FUNCTIONS
    // MultiCalc
    void Alloc_internal();
double CalculateEquilibriumState( long int typeMin, long int& NumIterFIA, long int& NumIterIPM );
    void InitalizeGEM_IPM_Data();
    void DC_LoadThermodynamicData( TNode* aNa = 0 );
    //DM 25.02.2014
    void Access_GEM_IMP_init();
    int get_sizeFIs () {return sizeFIs;}
    // acces for node class
    TSolMod * pTSolMod (int xPH);


    long int CheckMassBalanceResiduals(double *Y );
    double ConvertGj_toUniformStandardState( double g0, long int j, long int k );
    double PhaseSpecificGamma( long int j, long int jb, long int je, long int k, long int DirFlag = 0L ); // Added 26.06.08

    double HelmholtzEnergy( double x );
    double InternalEnergy( double TC, double P );

};

// ???? syp->PGmax
typedef enum {  // Symbols of thermodynamic potential to minimize
    G_TP    =  'G',   // Gibbs energy minimization G(T,P)
    A_TV    =  'A',   // Helmholts energy minimization A(T,V)
    U_SV    =  'U',   // isochoric-isentropicor internal energy at isochoric conditions U(S,V)
    H_PS    =  'H',   // isobaric-isentropic or enthalpy H(P,S)
    _S_PH   =  '1',   // negative entropy at isobaric conditions and fixed enthalpy -S(P,H)
    _S_UV   =  '2'    // negative entropy at isochoric conditions and fixed internal energy -S(P,H)

} THERM_POTENTIALS;

typedef enum {  // Symbols of thermodynamic potential to minimize
    G_TP_    =  0,   // Gibbs energy minimization G(T,P)
    A_TV_    =  1,   // Helmholts energy minimization A(T,V)
    U_SV_    =  2,   // isochoric-isentropicor internal energy at isochoric conditions U(S,V)
    H_PS_    =  3,   // isobaric-isentropic or enthalpy H(P,S)
    _S_PH_   =  4,   // negative entropy at isobaric conditions and fixed enthalpy -S(P,H)
    _S_UV_   =  5    // negative entropy at isochoric conditions and fixed internal energy -S(P,H)

} NUM_POTENTIALS;

// kg44: this is not correctly implemented!! 
//double A_P( double x, double x2=0.);
//double U_TP( double TC, double P);


typedef enum {  // Field index into outField structure
  f_pa_PE = 0,  f_PV,  f_PSOL,  f_PAalp,  f_PSigm,
  f_Lads,  f_FIa,  f_FIat
} MULTI_STATIC_FIELDS;

typedef enum {  // Field index into outField structure
  f_sMod = 0,  f_LsMod,  f_LsMdc,  f_B,  f_DCCW,
  f_Pparc,  f_fDQF,  f_lnGmf,  f_RLC,  f_RSC,
  f_DLL,  f_DUL,  f_Aalp,  f_Sigw,  f_Sigg,
  f_YOF,  f_Nfsp,  f_MASDT,  f_C1,  f_C2,
  f_C3,  f_pCh,  f_SATX,  f_MASDJ,  f_SCM,
  f_SACT,  f_DCads,
// static
  f_pa_DB,  f_pa_DHB,  f_pa_EPS,  f_pa_DK,  f_pa_DF,
  f_pa_DP,  f_pa_IIM,  f_pa_PD,  f_pa_PRD,  f_pa_AG,
  f_pa_DGC,  f_pa_PSM,  f_pa_GAR,  f_pa_GAH,  f_pa_DS,
  f_pa_XwMin,  f_pa_ScMin,  f_pa_DcMin,  f_pa_PhMin,  f_pa_ICmin,
  f_pa_PC,  f_pa_DFM,  f_pa_DFYw,  f_pa_DFYaq,  f_pa_DFYid,
  f_pa_DFYr,  f_pa_DFYh,  f_pa_DFYc,  f_pa_DFYs,  f_pa_DW,
  f_pa_DT,  f_pa_GAS,  f_pa_DG,  f_pa_DNS,  f_pa_IEPS,
  f_pKin,  f_pa_DKIN,  f_mui,  f_muk,  f_muj,
  f_pa_PLLG,  f_tMin,  f_dcMod,
//new
    f_kMod, f_LsKin, f_LsUpt, f_xICuC, f_PfFact,
    f_LsESmo, f_LsISmo, f_SorMc, f_LsMdc2, f_LsPhl

} MULTI_DYNAMIC_FIELDS;



#endif   //_ms_multi_h

