//-------------------------------------------------------------------
// $Id: databr.h 879 2013-10-10 14:47:33Z dmitrieva $
/// \file databr.h
/// Contains definition of the DATABR structure - data bridge between
/// GEMS3K and another code.
//
/// \struct DATABR databr.h
/// DataBRidge defines the structure of node-dependent data for
/// exchange between the coupled GEM IPM and FMT code parts.
/// DATABR structure is used in TNode and TNodeArray classes.
//
// Copyright (c) 2003-2011 by D.Kulik, S.Dmytriyeva, F.Enzmann, W.Pfingsten
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
#ifndef _DataBr_H_
#define _DataBr_H_



typedef struct  /// DATABR - template node data bridge structure
{
   long int
     NodeHandle,    ///< Node identification handle, not used in calculaions on TNode level
     NodeTypeHY,    ///< Node type code (hydraulic), not used on TNode level see NODETYPE
     NodeTypeMT,    ///< Node type (mass transport), not used on TNode level see NODETYPE
     NodeStatusFMT, ///< Node status code in FMT part, not used on TNode level see NODECODEFMT
     NodeStatusCH,  ///< Node status code in GEM (input and output) see NODECODECH
     IterDone;      ///< Number of iterations performed by GEM IPM in the last run

/*  these important data array dimensions are provided in the DATACH structure
   long int
    nICb,       ///< Number of Independent Components kept in the DATABR memory structure (<= nIC)
    nDCb,      	///< Number of Dependent Components kept in the DATABR memory structure (<=nDC)
    nPHb,     	///< Number of Phases to be kept in the DATABR structure (<= nPH)
    nPSb,       ///< Number of Phases-solutions (multicomponent phases) to be kept in the DATABR memory structure (<= nPS)
*/
//      Usage of this variable (DB - data bridge)                           MT-DB DB-GEM GEM-DB DB-MT
   double
// \section Chemical scalar variables
    TK,     ///< Node temperature T (Kelvin)                     	         +      +      -     -
    P, 	    ///< Node Pressure P (Pa)                         	             +      +      -     -
    Vs,     ///< Volume V of reactive subsystem  (m3)                       (+)    (+)     +     +
    Vi,     ///< Volume of inert subsystem (m3)          	                 +      -      -     +
    Ms,     ///< Mass of reactive subsystem (kg)         	                 +     (+)     -     -
    Mi,     ///< Mass of inert subsystem (kg)             	                 +      -      -     +

    Gs,     ///< Total Gibbs energy of the reactive subsystem (J/RT) (norm)  -      -      +     +
    Hs, 	///< Total enthalpy of reactive subsystem (J) (reserved)         -      -      +     +
    Hi,     ///< Total enthalpy of inert subsystem (J) (reserved)            +      -      -     +

    IC,     ///< Effective aqueous ionic strength (molal)                    -      -      +     +
    pH,     ///< pH of aqueous solution in the activity scale (-log10 molal) -      -      +     +
    pe,     ///< pe of aqueous solution in the activity scale (-log10 molal) -      -      +     +
    Eh,     ///< Eh of aqueous solution (V)                                  -      -      +     +
    Tm,     ///< Actual total simulation time (s)                            +      +      -     -
    dt      ///< Actual time step (s) - needed for TKinMet, can change!      +      +     (+)   (+)
#ifdef NODEARRAYLEVEL
      ,
// \section  FMT variables (units or dimensionsless) - to be used for storing them
//  at the nodearray level, normally not used in the single-node FMT-GEM coupling
    Dif,    ///< General diffusivity of disolved matter (m2/s)
    Vt,		///< Total volume of the node (m3) (Vs + Vi)
    vp,		///< Advection velocity (in pores)  (m/s)
    eps,	///< Effective (actual) porosity normalized to 1
    Km,		///< Actual permeability (m2)
    Kf,		///< Actual Darcy`s constant (m2/s)
    S,		///< Specific storage coefficient, dimensionless (default 1.0)
            ///<  if not 1.0 then can be used as mass scaling factor relative to Ms for the
            ///<  bulk composition/speciation of reactive sub-system in this node
    Tr,     ///< transmissivity m2/s
    h,		///< Actual hydraulic head (hydraulic potential) (m)
    rho,	///< Actual carrier density for density-driven flow (kg/m3)
    al,		///< Specific longitudinal dispersivity of porous media (m)
    at,		///< Specific transversal dispersivity of porous media (m)
    av,		///< Specific vertical dispersivity of porous media (m)
    hDl,	///< Hydraulic longitudinal dispersivity (m2/s)
    hDt,	///< Hydraulic transversal dispersivity (m2/s)
    hDv,	///< Hydraulic vertical dispersivity (m2/s)
    nto	    ///< Tortuosity factor (dimensionless)
#endif
   ;
// \section Data arrays - dimensions nICb, nDCb, nPHb, nPSb see in the DATACH structure
// exchange of values occurs through lists of indices, e.g. xIC, xDC, xPH from DATACH

//      Usage of this variable (DB = data bridge)                                               MT-DB DB-GEM GEM-DB DB-MT
   double
// IC (stoichiometry units)
    *bIC,  ///< Bulk composition of (reactive part of) the system (moles)[nICb]                      +      +      -     -
    *rMB,  ///< Mass balance residuals (moles) [nICb]                                                -      -      +     +
    *uIC,  ///< Chemical potentials of ICs (dual GEM solution) normalized scale (mol/mol)[nICb]      -      -      +     +
// DC (species) in reactive subsystem
    *xDC,  ///< Speciation - amounts of DCs in equilibrium state - primal GEM solution(moles)[nDCb] (+)    (+)     +     +
    *gam,  ///< Activity coefficients of DCs in their respective phases  [nDCb]                     (+)    (+)     +     +
// Metastability/kinetic controls
    *dul,  ///< Upper additional metastability restrictions AMR on amounts of DCs (moles) [nDCb]     +      +      -     -
    *dll,  ///< Lower AMR on amounts of DCs (moles) [nDCb]                                           +      +      -     -
// Phases in reactive subsystem
    *aPH,  ///< Specific surface areas of phases (m2/kg) [nPHb], can change in TKinMet               +      +     (+)   (+)
    *xPH,  ///< Amounts of phases in equilibrium state (moles) [nPHb]                                -      -      +     +
    *vPS,  ///< Volumes of multicomponent phases (m3)  [nPSb]                                        -      -      +     +
    *mPS,  ///< Masses of multicomponent phases (kg)    [nPSb]                                       -      -      +     +
    *bPS,  ///< Bulk elemental compositions of multicomponent phases (moles) [nPSb][nICb]            -      -      +     +
    *xPA,  ///< Amount of carrier (sorbent or solvent) in multicomponent phases [nPSb]               -      -      +     +
           ///< (+) can be used as input in "smart initial approximation" mode of GEM IPM-2 algorithm
    *bSP,  ///< Output bulk composition of the equilibrium solid part of the system, moles   [nICb]  -      -      +     +
   // Metastability/kinetic controls on phases-solutions (added in devPhase branch)
    *amru, ///< Upper AMRs on masses of DCs (kg) [nPSb]                                              +      +      -     -
    *amrl; ///< Lower AMRs on masses of DCs (kg) [nPSb]                                              +      +      -     -
}
DATABR;

typedef DATABR*  DATABRPTR;

/// \enum NODECODECH NodeStatus codes with respect to GEMIPM calculations
/*typedef*/ enum NODECODECH {
 NO_GEM_SOLVER= 0,   ///< No GEM re-calculation needed for this node
 NEED_GEM_AIA = 1,   ///< Need GEM calculation with LPP (automatic) initial approximation (AIA)
 OK_GEM_AIA   = 2,   ///< OK after GEM calculation with LPP AIA
 BAD_GEM_AIA  = 3,   ///< Bad (not fully trustful) result after GEM calculation with LPP AIA
 ERR_GEM_AIA  = 4,   ///< Failure (no result) in GEM calculation with LPP AIA
 NEED_GEM_SIA = 5,   ///< Need GEM calculation with no-LPP (smart) IA, SIA
                     ///<   using the previous speciation (full DATABR lists only)
 OK_GEM_SIA   = 6,   ///< OK after GEM calculation with SIA
 BAD_GEM_SIA  = 7,   ///< Bad (not fully trustful) result after GEM calculation with SIA
 ERR_GEM_SIA  = 8,   ///< Failure (no result) in GEM calculation with SIA
 T_ERROR_GEM  = 9    ///< Terminal error has occurred in GEMS3K (e.g. memory corruption). Restart is required.
} /*NODECODECH*/;


// \typedef NODECODEFMT Node status codes set by the FMT (FluidMassTransport) part
typedef enum {
 No_nodearray  = -1, ///< Indicates that no node transport properties are present in this DATABR and DBR file
 No_transport  = 0,  ///< Chemical calculations only, no transport coupled
 Initial_RUN   = 1,
 OK_Hydraulic  = 2,
 BAD_Hydraulic = 3,  ///< insufficient convergence
 OK_Transport  = 4,
 BAD_Transport = 5,  ///< insufficient convergence
 NEED_RecalcMT = 6,
 OK_MassBal    = 7,
 OK_RecalcPar  = 8,
 Bad_Recalc    = 9,
 Bad_Time_Step = 10  ///< converged but TKinMet finds that time step must be less (returns in dt)

} NODECODEFMT; /// Node status codes set by the FMT (FluidMassTransport) part

typedef enum {  /// Node type codes controlling hydraulic/mass-transport behavior
  normal       = 0, ///< normal node
// boundary condition node
  NBC1source   = 1, ///< Dirichlet source ( constant concentration )
  NBC1sink    = -1, ///< Dirichlet sink
  NBC2source   = 2, ///< Neumann source ( constant gradient )
  NBC2sink    = -2, ///< Neumann sink
  NBC3source   = 3, ///< Cauchy source ( constant flux )
  NBC3sink    = -3, ///< Cauchy sink
  INIT_FUNK    = 4  ///< functional conditions (e.g. input time-depended functions)
} NODETYPE;

typedef enum {  /// Node type codes controlling Indexation Code
   undefi  = 0,    ///<
   nICbi   = -1,   ///< index of Independent Components kept in the DBR file and DATABR memory structure
   nDCbi   = -2,   ///< index of Dependent Components kept in the DBR file and DATABR memory structure
   nPHbi   = -3,   ///< index of Phases to be kept in the DBR file and DATABR structure
   nPSbi   = -4,   ///< index of Phases-solutions (multicomponent phases) to be kept in the DBR file and DATABR memory structure
   nPSbnICbi = -5  ///<  index of multicomponent phases
 } NODEINEX;

typedef enum {  /// Field index into outField structure
f_NodeHandle = 0,f_NodeTypeHY,f_NodeTypeMT,f_NodeStatusFMT,f_NodeStatusCH,
f_IterDone, f_TK, f_P, f_Vs,f_Vi,
f_Ms, f_Mi, f_Hs, f_Hi, f_Gs,
f_IS, f_pH, f_pe, f_Eh,
f_Tm, f_dt,
//#ifdef NODEARRAYLEVEL
f_Dif,f_Vt, f_vp, f_eps,
f_Km, f_Kf, f_S,  f_Tr, f_h,
f_rho,f_al, f_at, f_av, f_hDl,
f_hDt, f_hDv, f_nto,
//#endif
    // dynamic arrays (52-38=14+2new)
f_bIC, f_rMB, f_uIC, f_xDC, f_gam,
f_dll, f_dul, f_aPH, f_xPH, f_vPS,
f_mPS, f_bPS, f_xPA, f_bSP,
f_amru, f_amrl,
// only for VTK format output
f_mPH, f_vPH, f_m_t, f_con, f_mju, f_lga

} DATABR_FIELDS;


#endif

// -----------------------------------------------------------------------------
// end of _DataBr_h

