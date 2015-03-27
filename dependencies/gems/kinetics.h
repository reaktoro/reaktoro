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
#ifndef _Kinetics_H_
#define _Kinetics_H_

#include "datach.h"
#include "databr.h"
#include "m_const.h"

typedef struct  /// KINETICS - data structure for computing DC AMRs in phases
{               /// DATACH indexation throughout, no I/O file exchange

//      Usage of this variable                                              node GEM-in GEM-out node
   double
// \section Chemical scalar variables
    T,     ///< Node temperature T (Kelvin)                     	         +      +      -     -
    P, 	    ///< Node Pressure P (Pa)                         	             +      +      -     -
    kTau,   // time -< MULTI
    kdT,    // time step - <- MULTI
    Tc,    // temperature K
    Pc;    // pressure bar

    char (*kMod)[6];  ///< new: Codes for built-in kinetic models [nPH]; k is the phase index

    // TKinMet stuff  <- MULTI
    long int
    *LsKin,  ///< new: number of parallel reactions nPRk[k]; number of species in activity products nSkr[k];
             /// number of parameter coeffs in parallel reaction term nrpC[k]; number of parameters
             /// per species in activity products naptC[k]; nAscC number of parameter coefficients in As correction;
             /// nFaces[k] number of (separately considered) crystal faces or surface patches ( 1 to 4 ) [nPH][6]
    *LsUpt,  ///< new: number of uptake kinetics model parameters (coefficients) numpC[k]; reserved [nPS][2]

    *xSKrC,  ///< new: Collected array of aq/gas/sorption species indexes used in activity products (-> += LsKin[k][1])
    (*ocPRkC)[2], ///< new: Collected array of operation codes for kinetic parallel reaction terms (-> += LsKin[k][0])
                  /// and indexes of faces (surface patches)
    *xICuC;  ///< new: Collected array of IC species indexes used in partition (fractionation) coefficients  ->L1[k]   TBD
    double
    *feSArC, ///< new: Collected array of fractions of surface area related to parallel reactions k-> += LsKin[k][0]
    *rpConC,  ///< new: Collected array of kinetic rate constants k-> += LsKin[k][0]*LsKin[k][2];
    *apConC,  ///< new:!! Collected array of parameters per species involved in activity product terms
            ///  k-> += LsKin[k][0]*LsKin[k][1]*LsKin[k][3];
    *AscpC,   /// new: parameter coefficients of equation for correction of specific surface area k-> += LsKin[k][4]
    *UMpcC;  ///< new: Collected array of uptake model coefficients k-> += L1[k]*LsUpt[k][0];

    // work arrays and counters from MULTI - intermediate, to be replaced by TNode data only
    long int
    N,      // Number of ICs
    L,      // Number of DCs
    Ls,     // Total number of DCs in phases-solutions
    FI,     // Number of phases
    FIs,    // Number of phases-solutions
    *L1,    //< l_a vector - number of DCs included into each phase [nPH]; copy of nDCinPH
    *LsPhl,  //< new: Number of phase links; number of link parameters; [nPH][2]
    (*PhLin)[2];  //< new: indexes of linked phases and link type codes (sum 2*LsPhl[k][0] over nPH)

    double
    *XF,    //< Output total number of moles of phases Xa[0:FI-1]
    *YF,    //< Approximation of X_a in the next IPM iteration [0:FI-1]
    *XFA,   //< Quantity of carrier in asymmetric phases Xwa, moles [FIs]
    *YFA,   //< Approximation of XFA in the next IPM iteration [0:FIs-1]
    *X,     //< DC quantities at eqstate x_j, moles - primal IPM solution [L]
    *Y,     //< Copy of x_j from previous IPM iteration [0:L-1]
    *Fx,    //< Dual DC chemical potentials defined via u_i and a_ji [L]
    *Wx,    //< Mole fractions Wx of DC in multi-component phases [L]
    *F,     //<Primal DC chemical potentials defined via g0_j, Wx_j and lnGam_j[L]
    *F0,    //< Excess Gibbs energies for (metastable) DC, mole/mole [L]
    *Falp,  //< phase stability index (PC==2) [FI]

    *Aalp,  //< Full vector of specific surface areas of phases (m2/g) [0:FI-1]
    *Sigw,  //< Specific surface free energy for phase-water interface (J/m2)   [0:FI-1]
    *Sigg,  //< Specific surface free energy for phase-gas interface (J/m2) (not yet used)  [0:FI-1], reserved
    *lPhc,  //< new: Collected array of phase link parameters (sum(LsPhl[k][1] over Fi)
    (*Xr0h0)[2],  //< mean r & h of particles (- pores), nm  [0:FI-1][2], reserved

    *XFs,    //< Current quantities of phases X_a at IPM iterations [0:FI-1]
    *Falps,  //< Current Karpov criteria of phase stability  F_a [0:FI-1]
    *DUL,    //< VG Vector of upper kinetic restrictions to x_j, moles [L]
    *DLL,    //< NG Vector of lower kinetic restrictions to x_j, moles [L]
    *PUL,    //< Vector of upper restrictions to multicomponent phases amounts [FIs]
    *PLL,    //< Vector of lower restrictions to multicomponent phases amounts [FIs]
    *PfFact, // new: phase surface area - volume shape factor (taken from TKinMet or set from TNode) [FI]
    *PrT,    // new: Total MWR rate (mol/s) for phases - TKinMet output [FI]
    *PkT,    // new: Total specific MWR rate (mol/m2/s) for phases - TKinMet output [FI]
    *PvT,    // new: Total one-dimensional MWR surface propagation velocity (m/s) - TKinMet output [FI]
    //  potentially can be extended to all solution phases?
    *emRd,   // new: output Rd values (partition coefficients) for end members (in uptake kinetics model) [Ls]
    *emDf,   // new: output Df values (fractionation coeffs.) for end members (in uptake kinetics model) [Ls]
    //
    *YOF,     //< Surface free energy parameter for phases (J/g) (to accomodate for variable phase composit
    *FVOL,    //< phase volumes, cm3 comment corrected DK 04.08.2009  [0:FI-1]
    *FWGT,    //< phase (carrier) masses, g                [0:FI-1]
    *Vol,     //< DC molar volumes, cm3/mol [L]
    *MM,      //< DC molar masses, g/mol [L]
    *Y_m,     //< Molalities of aqueous species and sorbates [0:Ls-1]
    *Y_la,    //< log activity of DC in multi-component phases (mju-mji0) [0:L-1]
    *IC_m;    //< Total IC molalities in aqueous phase (excl.solvent) [0:N-1]
    char
    *PHC,   //< Classifier of phases { a g f p m l x d h } [FI]
    *DCC,   //< Classifier of DC { TESKWL GVCHNI JMFD QPR <0-9>  AB  XYZ O } [L]
    *RLC,   //< Code of metastability constraints for DCs [L] enum DC_LIMITS
    *RSC,   //< Units of metastability/kinetic constraints for DCs  [L]
    *RFLC,  //< Classifier of restriction types for XF_a [FIs]
    *RFSC;  //< Classifier of restriction scales for XF_a [FIs]
    char  (*SM)[MAXDCNAME];  //< List of DC names in the system [L]
    char  (*SF)[MAXPHNAME+MAXSYMB];  //< List of phase names in the system [FI]
    double
    IC,     //< Effective aqueous ionic strength (molal)                    -      -      +     +
    pH,     //< pH of aqueous solution in the activity scale (-log10 molal) -      -      +     +
    pe,     //< pe of aqueous solution in the activity scale (-log10 molal) -      -      +     +
    Eh;     //< Eh of aqueous solution (V)                                  -      -      +     +

}
KINETICS;

typedef KINETICS*  KINETICSPTR;

class TKinetics
{
    KINETICS kin;
    // How to link to KINETICS?

    long int sizeFI;      ///< current size of phKinMet
    TKinMet* (*phKinMet); ///< size current FI -   number of phases
    bool load;

public:

    /// This allocation is used only in standalone GEMS3K
    TKinetics( DATACH *csd, DATABR *sbc )
    {

            sizeFI = 0;
            phKinMet = 0;

         load = false;
         Alloc_TKinMet( csd->nPH );
         sizeFI = csd->nPH;
    }

~TKinetics( )
{          // destructor
    Free_TKinMet( );
}

    void set_def( void );
    void Alloc_TKinMet( long int newFI );
    void Free_TKinMet();

    void Set_DC_limits( long int Mode );

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


};

#endif

// -----------------------------------------------------------------------------
// end of _Activity_h

