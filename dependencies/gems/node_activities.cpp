//-------------------------------------------------------------------
// $Id: node_activities.cpp 799 2013-03-17 12:33:51Z kulik $
//
/// \file node_activities.cpp
/// Implementation of chemistry-specific functions for kinetics & metastability
/// for the IPM convex programming Gibbs energy minimization algorithm
//
// Copyright (c) 2014  D.Kulik, A.Leal
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

#include <cmath>
#include<iomanip>

#include "node.h"
#include "m_param.h"
#include "activities.h"


// Generic access methods that use the new TActivity class
// set temperature (in units of K)
void TNode::setTemperature(const double T)
{
    if( T > 0. && T < 1e9 )
    {
        CNode->TK = T;
        ACTIVITY* ap = atp->GetActivityDataPtr();
        ap->TK = T;
        ap->RT = T*R_CONSTANT;
    }
}

// set pressure (in units of Pa)
void TNode::setPressure(const double P)
{
    if( P >= 0. && P < 1e12 )
    {
        CNode->P = P;
        ACTIVITY* ap = atp->GetActivityDataPtr();
        ap->P = P/bar_to_Pa;
    };
}

// compute standard Gibbs energies (so far only from P,T interpolation)
//   (in J/mol)
void TNode::updateStandardGibbsEnergies()
{
    bool norm = false;
    ACTIVITY* ap = atp->GetActivityDataPtr();
    for( long int j=0; j<CSD->nDC; j++ )
    {
        ap->tpp_G[j] = this->DC_G0( j, CNode->P, CNode->TK, norm );
    }
}

void TNode::updateStandardVolumes()
{
    ACTIVITY* ap = atp->GetActivityDataPtr();
    for( long int j=0; j<CSD->nDC; j++ )
    {
        ap->Vol[j] = this->DC_V0( j, CNode->P, CNode->TK ) * bar_to_Pa *10.;
    }
}

void TNode::updateStandardEnthalpies()
{
    ACTIVITY* ap = atp->GetActivityDataPtr();
    for( long int j=0; j<CSD->nDC; j++ )
    {
        ap->H0[j] = this->DC_H0( j, CNode->P, CNode->TK );
    }
}

void TNode::updateStandardEntropies()
{
    ACTIVITY* ap = atp->GetActivityDataPtr();
    for( long int j=0; j<CSD->nDC; j++ )
    {
        ap->S0[j] = this->DC_S0( j, CNode->P, CNode->TK );
    }
}

void TNode::updateStandardHeatCapacities()
{
    ACTIVITY* ap = atp->GetActivityDataPtr();
    for( long int j=0; j<CSD->nDC; j++ )
    {
        ap->Cp0[j] = this->DC_Cp0( j, CNode->P, CNode->TK );
    }
}

void TNode::updateThermoData()
{
    bool norm = true;
    ACTIVITY* ap = atp->GetActivityDataPtr();
    for( long int j=0; j<CSD->nDC; j++ )
    {
        ap->tpp_G[j] = this->DC_G0( j, CNode->P, CNode->TK, norm );
        ap->Vol[j] = this->DC_V0( j, CNode->P, CNode->TK )* 10. * bar_to_Pa;
        ap->H0[j] = this->DC_H0( j, CNode->P, CNode->TK );
        ap->S0[j] = this->DC_S0( j, CNode->P, CNode->TK );
        ap->Cp0[j] = this->DC_Cp0( j, CNode->P, CNode->TK );
    }
}

// set speciation (in units of moles) - scaled down?
// then calculates amounts of phases
void TNode::setSpeciation( const double* n )
{
    ACTIVITY* ap = atp->GetActivityDataPtr();
    for( long int j=0; j<CSD->nDC; j++ )
    {
        ap->Y[j] = n[j] * ap->SizeFactor;
        ap->X[j] = n[j] * ap->SizeFactor;
//        ap.XY[j] *= ScFact;
//        ap.XU[j] *= ScFact;
    }
    atp->TotalPhasesAmounts( ap->X, ap->XF, ap->XFA );
    atp->TotalPhasesAmounts( ap->Y, ap->YF, ap->YFA );
}

// compute concentrations in all phases
void TNode::updateConcentrations()
{
    ACTIVITY* ap = atp->GetActivityDataPtr();
    atp->CalculateConcentrations( ap->X, ap->XF, ap->XFA);
}

//initialize models of mixing (TSolMod class instances) before GEM run
// also initializes G0 values
void TNode::initActivityCoefficients()
{
    long int k, j, jb, je, retCode;

    je = 0;

    ACTIVITY* ap = atp->GetActivityDataPtr();

    for( k=0; k< ap->FI; k++ )
    {
      jb = je;
      je += ap->L1[k];
      // load t/d data from DC - to be extended for DCH->H0, DCH->S0, DCH->Cp0, DCH->DD
      // depending on the presence of these arrays in DATACH and Multi structures
       for( j=jb; j<je; j++ )
       {
         ap->G0[j] = atp->ConvertGj_toUniformStandardState( ap->tpp_G[j], j, k );
       }
    }
    retCode = atp->CalculateActivityCoefficients( LINK_TP_MODE );
    //   if(retCode)
           // Errors
}

// compute activity coefficients on GEM iterations
void TNode::updateActivityCoefficients()
{
   long int retCode;
//   ACTIVITY* ap = atp->GetActivityDataPtr();
   retCode = atp->CalculateActivityCoefficients( LINK_UX_MODE );
//   if(retCode)
       // Errors
}

// compute integral phase properties after GEM run (TBD)
void TNode::getIntegralPhaseProperties()
{
    long int retCode;
//    ACTIVITY* ap = atp->GetActivityDataPtr();
    retCode = atp->CalculateActivityCoefficients( LINK_PP_MODE );
    atp->StabilityIndexes( );
    //   if(retCode)
           // Errors
}

// compute primal chemical potentials
void TNode::updateChemicalPotentials()
{
    ACTIVITY* ap = atp->GetActivityDataPtr();
    atp->PrimalChemicalPotentials( ap->F, ap->Y, ap->YF, ap->YFA );
}

// gets the total Gibbs energy of the system (at GEM iteration)
double TNode::updateTotalGibbsEnergy()
{
   double Gn;
//   ACTIVITY* ap = atp->GetActivityDataPtr();
   Gn = atp->GX( 0 );
   return Gn;
}

// compute primal activities
void TNode::updateActivities()
{
    ACTIVITY* ap = atp->GetActivityDataPtr();
    for( long int j=0; j<CSD->nDC; j++ )
    {
        ap->lnAct[j] = ap->F[j] - ap->tpp_G[j]/ap->RT;
    }
}

void TNode::updateChemicalData()
{
    ;
}

// Initializes and fills out the TActivity class instance
// (so far, by copying the scalars and pointers from MULTI;
// at the next stage, the data should be directly read into this instance instead of MULTI;
// finally, the TActivity functionality should be used in GEM IPM)
//
void
TNode::InitCopyActivities( DATACH* csd, MULTI* mp, DATABR* cnd )
{

   atp = new TActivity( csd, cnd );

   ACTIVITY* ap = atp->GetActivityDataPtr();

   // \section Chemical scalar variables
      ap->TK = cnd->TK;     //< Node temperature T (Kelvin)                     	         +      +      -     -
      ap->P = cnd->P; 	    //< Node Pressure P (Pa)                         	             +      +      -     -
      ap->RT = mp->RT;     // RT product
      ap->ln5551 = mp->ln5551; // ln(H2O_mol_to_kg) = log(55.50837344)
   //
       ap->IC = mp->IC;     //< Effective aqueous ionic strength (molal)                    -      -      +     +
       ap->pH = mp->pH;     //< pH of aqueous solution in the activity scale (-log10 molal) -      -      +     +
       ap->pe = mp->pe;     //< pe of aqueous solution in the activity scale (-log10 molal) -      -      +     +
       ap->Eh = mp->Eh;     //< Eh of aqueous solution (V)                                  -      -      +     +

       ap->Pc = mp->Pc;
       ap->Tc = mp->Tc;  // T in K, P in bar
      // Here TSolMod, TSorpMod, TKinMet parameters & coefficients moved from TMULTI MULTI
      //   DK, AL September 2014
         // TSolMod stuff
       ap->LsMod = mp->LsMod; //< Number of interaction parameters. Max parameter order (cols in IPx),
                   //< and number of coefficients per parameter in PMc table [3*nPS]
       ap->LsMdc = mp->LsMdc; //<  for multi-site models: [3*nPS] - number of nonid. params per component;
                  // number of sublattices nS; number of moieties nM
       ap->LsMdc2 = mp->LsMdc2; //<  new: [3*nPS] - number of DQF coeffs; reciprocal coeffs per end member;
                 // reserved
       ap->IPx =  mp->IPx;   //< Collected indexation table for interaction parameters of non-ideal solutions
                   //< ->LsMod[k,0] x LsMod[k,1]   over nPS
       ap->LsPhl = mp->LsPhl;  //< new: Number of phase links; number of link parameters; [nPH][2]
       ap->PhLin = mp->PhLin;  //< new: indexes of linked phases and link type codes (sum 2*LsPhl[k][0] over nPH)

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
        for( long int i=0; i<5; i++ )
        {
          ap->denW[i] = mp->denW[i];     //< Density of water, first T, second T, first P, second P derivative for Tc,Pc
          ap->denWg[i] = mp->denWg[i];   //< Density of steam for Tc,Pc
          ap->epsW[i] = mp->epsW[i];     //< Diel. constant of H2O(l)for Tc,Pc
          ap->epsWg[i] = mp->epsWg[i];   //< Diel. constant of steam for Tc,Pc
        }
        ap->PMc = mp->PMc;      //< Collected interaction parameter coefficients for the (built-in) non-ideal mixing models -> LsMod[k,0] x LsMod[k,2]
        ap->DMc = mp->DMc;      //< Non-ideality coefficients f(TPX) for DC -> L1[k] x LsMdc[k][0]
        ap->MoiSN = mp->MoiSN;  //< End member moiety- site multiplicity number tables ->  L1[k] x LsMdc[k][1] x LsMdc[k][2]
        ap->SitFr = mp->SitFr;  //< Tables of sublattice site fractions for moieties -> LsMdc[k][1] x LsMdc[k][2]
        ap->lPhc = mp->lPhc;    //< new: Collected array of phase link parameters (sum(LsPhl[k][1] over Fi)
        ap->DQFc = mp->DQFc;    //< new: Collected array of DQF parameters for DCs in phases -> L1[k] x LsMdc2[k][0]
       //  *rcpc,  //< new: Collected array of reciprocal parameters for DCs in phases -> L1[k] x LsMdc2[k][1]

        // TSorpMod stuff
        ap->LsESmo = mp->LsESmo; //< new: number of EIL model layers; EIL params per layer; CD coefs per DC; reserved  [nPS][4]
        ap->LsISmo = mp->LsISmo; //< new: number of surface sites; isotherm coeffs per site; isotherm coeffs per DC; max.denticity of DC [nPS][4]
        ap->xSMd = mp->xSMd;     //< new: denticity of surface species per surface site (site allocation) (-> L1[k]*LsISmo[k][3]] )
        ap->SATX = mp->SATX;     //< Setup of surface sites and species (will be applied separately within each sorption phase) [Lads]
                    // link indexes to surface type [XL_ST]; sorbent em [XL_EM]; surf.site [XL-SI] and EDL plane [XL_SP]
        ap->EImc = mp->EImc;     //< new: Collected EIL model coefficients k -> += LsESmo[k][0]*LsESmo[k][1]
        ap->mCDc = mp->mCDc;     //< new: Collected CD EIL model coefficients per DC k -> += L1[k]*LsESmo[k][2]
        ap->IsoPc = mp->IsoPc;   //< new: Collected isotherm coefficients per DC k -> += L1[k]*LsISmo[k][2];
        ap->IsoSc = mp->IsoSc;   //< new: Collected isotherm coeffs per site k -> += LsISmo[k][0]*LsISmo[k][1];
         // TSorpMod & TKinMet stuff
        ap->SorMc = mp->SorMc;   //< new: Phase-related kinetics and sorption model parameters: [nPS][16]
                  //< in the same order as from Asur until fRes2 in TPhase

        ap->H0 = mp->H0;     //< DC molar enthalpies, reserved [L]
   //     double *A0;        //< DC molar Helmholtz energies, reserved [L]
        ap->S0 = mp->S0;     //< DC molar entropies, reserved [L]
        ap->Cp0 = mp->Cp0;   //< DC molar heat capacity, reserved [L]
//      Usage of this variable                                                                     node GEM-in GEM-out node
       ap->G0 = mp->G0;   //< Input normalized g0_j(T,P) for DC at unified standard scale[nDC]                  +      +      -     -
       ap->G = mp->G;     //< Normalized DC energy function c(j), incl. primal mu mole/mole [nDC]               -      -      +     +

//      ap->lnAct = mp->lnAct;   //< ln of DC activities (mu - mu0) [nDC]                                       -      -      +     +
       ap->lnGam = mp->lnGam;   //< ln of DC activity coefficients in unified (mole-fraction) scale [nDC]       -      -      +     +
       ap->lnGmo = mp->lnGmo;   //< Copy of lnGam from previous IPM iteration (reserved)
       // TSolMod stuff (detailed output on partial energies of mixing)   lnGam components
       ap->lnDQFt = mp->lnDQFt;  //< new: DQF terms adding to overall activity coefficients [nDCs]              -      -      +     +
       ap->lnRcpt = mp->lnRcpt; //< new: reciprocal terms adding to overall activity coefficients [nDCs]        -      -      +     +
       ap->lnExet = mp->lnExet; //< new: excess energy terms adding to overall activity coefficients [nDCs]     -      -      +     +
       ap->lnCnft = mp->lnCnft; //< new: configurational terms adding to overall activity [nDCs]                -      -      +     +
       // Takeover from MULTI
       ap->Gamma = mp->Gamma;   //< DC activity coefficients in molal or other phase-specific scale [0:L-1]
       ap->lnGmf = mp->lnGmf;   //< ln of initial DC activity coefficients for correcting G0 [0:L-1]
       ap->lnGmM = mp->lnGmM;   //< ln of DC pure gas fugacity (or metastability) coeffs or DDF correction [0:L-1]

       // TSorpMod stuff (TBD)
       ap->lnScalT = mp->lnScalT; //< new: Surface/volume scaling activity correction terms [nDCs]              -      -      +     +
       ap->lnSACT = mp->lnSACT;   //< new: ln isotherm-specific SACT for surface species [nDCs]                 -      -      +     +
       ap->lnGammF = mp->lnGammF; //< new: Frumkin or BET non-electrostatic activity coefficients [nDCs]        -      -      +     +
       ap->CTerms = mp->CTerms;   //< new: Coulombic terms (electrostatic activity coefficients) [nDCs]         -      -      +     +

       ap->VPh = mp->VPh;     //< Volume properties for mixed phases [FIs]                     -      -      +     +
       ap->GPh = mp->GPh;     //< Gibbs energy properties for mixed phases [FIs]               -      -      +     +
       ap->HPh = mp->HPh;     //< Enthalpy properties for mixed phases [FIs]                   -      -      +     +
       ap->SPh = mp->SPh;     //< Entropy properties for mixed phases [FIs]                    -      -      +     +
       ap->CPh = mp->CPh;     //< Heat capacity Cp properties for mixed phases [FIs]           -      -      +     +
       ap->APh = mp->APh;     //< Helmholtz energy properties for mixed phases [FIs]           -      -      +     +
       ap->UPh = mp->UPh;     //< Internal energy properties for mixed phases [FIs]            -      -      +     +

   // Coding - takeover from MULTI
       ap->sMod = mp->sMod;   //< new: Codes for built-in mixing models of multicomponent phases [FIs]
       ap->dcMod = mp->dcMod; //< Codes for PT corrections for dependent component data [L]
       ap->SM = mp->SM;       //< List of DC names in the system [L]
       ap->SF = mp->SF;       //< List of phase names in the system [FI]
       ap->ICC = mp->ICC;     //< Classifier of IC { e o h a z v i <int> } [N]
       ap->DCC = mp->DCC;     //< Classifier of DC { TESKWL GVCHNI JMFD QPR <0-9>  AB  XYZ O } [L]
       ap->PHC = mp->PHC;     //< Classifier of phases { a g f p m l x d h } [FI]
       ap->DCCW = mp->DCCW;   //< internal DC class codes for calculation of primal chemical potentials [L]

        ap->MM = mp->MM;      //< DC molar masses, g/mol [L]
        ap->FWGT = mp->FWGT;  //< phase (carrier) masses, g                [0:FI-1]
        ap->YOF = mp->YOF;    //< Surface free energy parameter for phases (J/g) (to accomodate for variable phase composition) [FI]

   // Work part (to compute concentrations and chem.potentials), takeover from MULTI
       ap->N = mp->N;       // Number of ICs
           ap->NR = mp->NR;       	//< NR - dimensions of R matrix
       ap->L = mp->L;       // Number of DCs
       ap->Ls = mp->Ls;     // Total number of DCs in phases-solutions
       ap->LO = mp->LO;     // LO -   index of water-solvent in DC list
       ap->FI = mp->FI;     // Number of phases
       ap->FIs = mp->FIs;   // Number of phases-solutions,
       ap->K2 = mp->K2;
       ap->L1 = mp->L1;    //< l_a vector - number of DCs included into each phase [nPH]; copy of nDCinPH

       ap->A = mp->A;      //< DC stoichiometry matrix A composed of a_ji [0:N-1][0:L-1]
       ap->Awt = mp->Awt;  //< IC atomic (molar) mass, g/mole [0:N-1]
       ap->XF = mp->XF;    //< Output total number of moles of phases Xa[0:FI-1]
       ap->YF = mp->YF;    //< Approximation of X_a in the next IPM iteration [0:FI-1]
       ap->XFA = mp->XFA;  //< Quantity of carrier in asymmetric phases Xwa, moles [FIs]
       ap->YFA = mp->YFA;  //< Approximation of XFA in the next IPM iteration [0:FIs-1]
       ap->X = mp->X;      //< DC quantities at eqstate x_j, moles - primal IPM solution [L]
       ap->Y = mp->Y;      //< Copy of x_j from previous IPM iteration [0:L-1]
       ap->Fx = mp->Fx;    //< Dual DC chemical potentials defined via u_i and a_ji [L]
       ap->Wx = mp->Wx;    //< Mole fractions Wx of DC in multi-component phases [L]
       ap->F = mp->F;      //< Primal DC chemical potentials defined via g0_j, Wx_j and lnGam_j[L]
       ap->F0 = mp->F0;    //< Excess Gibbs energies for (metastable) DC, mole/mole [L]
       ap->U = mp->U;      //< IC chemical potentials u_i (mole/mole) - dual IPM solution [N]
       ap->Falp = mp->Falp;  //< Phase stability index (PC==2) [FI]

       ap->GamFs = mp->GamFs;  //< Copy of activity coefficients Gamma [L]
       ap->fDQF = mp->fDQF;    //< Increments to molar G0 values of DCs from pure gas fugacities or DQF terms, normalized [L]
       ap->MU = mp->MU;        //< mu_j values of differences between dual and primal DC chem.potentials [L]
       ap->EMU = mp->EMU;      //< Exponents of DC increment to F_a criterion for phase [L]
       ap->NMU = mp->NMU;      //< DC increments to F_a criterion for phase [L]
       ap->Y_la = mp->Y_la;    //< log activity of DC in phases (mju-mji0) [0:L-1]
       ap->Y_m = mp->Y_m;      //< Molalities of aqueous species and sorbates [0:Ls-1]
       ap->Pparc = mp->Pparc;  //< Partial pressures or fugacities of pure DC, bar (Pc by default) [0:L-1]
       ap->EZ = mp->EZ;        //< Formula charge of DC in multi-component phases [0:Ls-1]
       ap->Vol = mp->Vol;      //< DC molar volumes, cm3/mol [L]
       ap->FVOL = mp->FVOL;    //< phase volumes, cm3 comment corrected DK 04.08.2009  [0:FI-1]
       ap->BF = mp->BF;        //< Output bulk compositions of multicomponent phases bf_ai[FIs][N]
       ap->BFC = mp->BFC;      //< Total output bulk composition of all solid phases [1][N]

       //< Numerical tolerances (scaled), taken over from MULTI
      ap->XwMinM = mp->XwMinM; //< Cutoff mole amount for elimination of water-solvent { 1e-13 }
      ap->ScMinM = mp->ScMinM; //< Cutoff mole amount for elimination of solid sorbent { 1e-13 }
      ap->DcMinM = mp->DcMinM; //< Cutoff mole amount for elimination of solution- or surface species { 1e-30 }
      ap->PhMinM = mp->PhMinM; //< Cutoff mole amount for elimination of non-electrolyte condensed phase { 1e-23 }
       //< insertion values (re-scaled to system size)
   //    DFYwM,  ///< Insertion mole amount for water-solvent { 1e-6 }
   //    DFYaqM, ///< Insertion mole amount for aqueous and surface species { 1e-6 }
   //    DFYidM, ///< Insertion mole amount for ideal solution components { 1e-6 }
   //    DFYrM,  ///< Insertion mole amount for major solution components (incl. sorbent) { 1e-6 }
   //    DFYhM,  ///< Insertion mole amount for minor solution components { 1e-6 }
   //    DFYcM,  ///< Insertion mole amount for single-component phase { 1e-6 }
   //    DFYsM,  ///< Insertion mole amount used in PhaseSelect() for a condensed phase component  { 1e-7 }
      ap->SizeFactor = mp->SizeFactor;  //< factor for re-scaling the cutoffs/insertions to the system size
   //    TMols,      ///< Input total moles in b vector before rescaling
   //    SMols,      ///< Standart total moles (upscaled) {1000}
   //    MBX,        ///< Total mass of the system, kg
   //    FX,    	    ///< Current Gibbs potential of the system in IPM, moles
   //    IC,         ///< Effective molal ionic strength of aqueous electrolyte
   //    pH,         ///< pH of aqueous solution
   //    pe,         ///< pe of aqueous solution
   //    Eh,         ///< Eh of aqueous solution, V
   //    DHBM,       ///< balance (relative) precision criterion
        ap->DSM = mp->DSM;  //< min amount of phase DS
   //    GWAT,       ///< used in ipm_gamma()
   //    YMET,       ///< reserved
   //    PCI,        ///< Current value of Dikin criterion of IPM convergence DK>=DX
   //    DXM,        ///< IPM convergence criterion threshold DX (1e-5)
   //    lnP,        ///< log Ptotal
   //    RT,         ///< RT: 8.31451*T (J/mole/K)
   //    FRT,        ///< F/RT, F - Faraday constant = 96485.309 C/mol
   //    Yw,         ///< Current number of moles of solvent in aqueous phase
   //    ln5551,     ///< ln(55.50837344)
       ap->ICmin = 1e-4; // mp->ICmin;  //< tolerance for minimum ionic strength to calculate aqueous activity models
       ap->aqsTail = mp->aqsTail; //< v_j asymmetry correction factor for aqueous species
       ap->lowPosNum = mp->lowPosNum,  //< Minimum mole amount considered in GEM calculations (MinPhysAmount = 1.66e-24)
       ap->logXw = mp->logXw;      //< work variable
       ap->logYFk = mp->logYFk;    //< work variable
       ap->YFk = mp->YFk;        //< Current number of moles in a multicomponent phase

       atp->Alloc_internal();

}


//--------------------- End of node_activities.cpp ---------------------------


