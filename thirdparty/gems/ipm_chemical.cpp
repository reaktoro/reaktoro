//-------------------------------------------------------------------
// $Id: ipm_chemical.cpp 771 2012-12-13 13:07:43Z kulik $
//
/// \file ipm_chemical.cpp
/// Implementation of chemistry-specific functions (concentrations,
/// activity coefficients, adsorption models etc.)
/// for the IPM convex programming Gibbs energy minimization algorithm
//
// Copyright (c) 1992-2012  D.Kulik, S.Dmytriyeva, K.Chudnenko
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

#include "m_param.h"
#ifndef IPMGEMPLUGIN
#include "service.h"
#include "stepwise.h"
#endif

// #define GEMITERTRACE

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of max.moles of surface species for SACT stabilization
///  to improve IPM-2 convergence at high SACT values.
///  xj0 values are placed as upper kinetic constraints
//
void TMulti::XmaxSAT_IPM2()
{
    long int i, j, ja, k, jb, je=0, ist=0, Cj, iSite[6];
    double XS0, xj0, XVk, XSk, XSkC, xj, Mm, rIEPS, xjn;

  if(!pm.DUL )   // not possible to install upper kinetic constraints!
      return;

  for( k=0; k<pm.FIs; k++ )
  { // loop over phases
     jb = je;
     je += pm.L1[k];
     if( pm.PHC[k] != PH_SORPTION )
          continue;

    if( pm.XFA[k] < pm.DSM ) // No sorbent retained by the IPM
        continue;
    if( pm.XF[k]-pm.XFA[k] < min( pm.lowPosNum, pm.DcMinM ) ) // may need qd_real in subtraction!
        continue;  // No surface species

    for(i=0; i<6; i++)
        iSite[i] = -1;

    // Extraction of site indices
    for( j=jb; j<je; j++ )
    {
        ja = j - ( pm.Ls - pm.Lads );
        if( pm.SATT[ja] != SAT_SOLV )
        {
            if( pm.DCC[j] == DC_PEL_CARRIER || pm.DCC[j] == DC_SUR_MINAL ||
                    pm.DCC[j] == DC_SUR_CARRIER ) continue;
/*!!!!!!*/  ist = pm.SATX[ja][XL_ST]; // / MSPN; MSPN = 2 - number of EDL planes
            continue;
        }
/*!!!!!!*/  ist = pm.SATX[ja][XL_ST]; //  / MSPN;
        iSite[ist] = j;     // To be checked !!!
    }

    for( j=jb; j<je; j++ )
    { // Loop over DCs
        if( pm.X[j] <= min( pm.lowPosNum, pm.DcMinM ) )
            continue;  // This surface DC has been killed by the IPM
        rIEPS = paTProfil->p.IEPS;
        ja = j - ( pm.Ls - pm.Lads );

        switch( pm.DCC[j] )  // code of species class
        {
        default: // pm.lnGam[j] = 0.0;
            continue;
        case DC_SSC_A0:
        case DC_SSC_A1:
        case DC_SSC_A2:
        case DC_SSC_A3:
        case DC_SSC_A4:
        case DC_WSC_A0:
        case DC_WSC_A1:
        case DC_WSC_A2:
        case DC_WSC_A3:
        case DC_WSC_A4:
        case DC_SUR_GROUP:
        case DC_IEWC_B:
        case DC_SUR_COMPLEX:
        case DC_SUR_IPAIR:
        case DC_IESC_A:
            // Calculate ist - index of surface type
/*!!!!!!*/  ist = pm.SATX[ja][XL_ST]; // / MSPN;
            // Cj - index of carrier DC
            Cj = pm.SATX[ja][XL_EM];
            if( Cj < 0 )
            {  // Assigned to the whole sorbent
                XVk = pm.XFA[k];
                Mm = pm.FWGT[k] / XVk;
            }
            else
            { // Assigned to one of the sorbent end-members
                XVk = pm.X[Cj];
                if( XVk < pm.DSM/10.0 )
                    continue; // This end-member is zeroed off by IPM
                Mm = pm.MM[Cj] * XVk/pm.XFA[k];  // mol.mass
            }
            XSk = pm.XFTS[k][ist]; // Tot.moles of sorbates on surf.type
            xj = pm.X[j];  // Current moles of this surf.species
//            a=1.0;  Frumkin factor - reserved for extension to FFG isotherm
            switch( pm.SATT[ja] )
            {
            case SAT_COMP: // Competitive surface species on a surface type
                if( iSite[ist] < 0 )
                    xjn = 0.0;
                else xjn = pm.X[iSite[ist]]; // neutral site does not compete!
                XS0 = max(pm.MASDT[k][ist], pm.MASDJ[ja][PI_DEN]);
                XS0 = XS0 * XVk * Mm / 1e6 * pm.Nfsp[k][ist];
                            // expected total in moles
                XSkC = XSk - xjn - xj; // occupied by the competing species
                if( XSkC < 0.0 )
                    XSkC = rIEPS;
                xj0 = XS0 - XSkC;    // expected moles of this sorbate
                if( xj0 > pm.lnSAC[ja][3] )
                    xj0 = pm.lnSAC[ja][3];
                if( xj0 < rIEPS )
                   xj0 = rIEPS;  // ensuring that it will not zero off
                pm.DUL[j] = xj0;
/*                if( pm.W1 != 1 && pm.IT > 0 && fabs( (pm.DUL[j] - oDUL)/pm.DUL[j] ) > 0.1 )
                {
cout << "XmaxSAT_IPM2 Comp. IT= " << pm.IT << " j= " << j << " oDUL=" << oDUL << " DUL=" << pm.DUL[j] << endl;
                }
*/                break;
            case SAT_L_COMP:
            case SAT_QCA_NCOMP:
            case SAT_QCA1_NCOMP:
            case SAT_QCA2_NCOMP:
            case SAT_QCA3_NCOMP:
            case SAT_QCA4_NCOMP:
            case SAT_BET_NCOMP:
            case SAT_FRUM_COMP:
            case SAT_FRUM_NCOMP:
            case SAT_PIVO_NCOMP:
            case SAT_VIR_NCOMP:
            case SAT_NCOMP: // Non-competitive surface species
                 xj0 = fabs( pm.MASDJ[ja][PI_DEN] ) * XVk * Mm / 1e6
                      * pm.Nfsp[k][ist]; // in moles
                 pm.DUL[j] = xj0 - rIEPS;
                 // Compare with old DUL from previous iteration!
/*                if( pm.W1 != 1 && pm.IT > 0 && fabs( (pm.DUL[j] - oDUL)/pm.DUL[j] ) > 0.1 )
                {
cout << "XmaxSAT_IPM2 Ncomp IT= " << pm.IT << " j= " << j << " oDUL=" << oDUL << " DUL=" << pm.DUL[j] << endl;
                }
*/                break;

            case SAT_SOLV:  // Neutral surface site (e.g. >O0.5H@ group)
                rIEPS = paTProfil->p.IEPS;
                XS0 = (max( pm.MASDT[k][ist], pm.MASDJ[ja][PI_DEN] ));
                XS0 = XS0 * XVk * Mm / 1e6 * pm.Nfsp[k][ist]; // in moles

                pm.DUL[j] =  XS0 - rIEPS;
                if( pm.DUL[j] <= rIEPS )
                   pm.DUL[j] = rIEPS;
                break;
// New methods added by KD 13.04.04
            case SAT_INDEF: // No SAT calculation
            default:        // pm.lnGam[j] = 0.0;
                break;
            }
        }
     }  // j
  } // k
}

/* // clearing pm.DUL constraints!
void TMulti::XmaxSAT_IPM2_reset()
{
    long int j, ja, k, jb, je=0;

  if(!pm.DUL )   // no upper kinetic constraints!
      return;

  for( k=0; k<pm.FIs; k++ )
  { // loop on phases
     jb = je;
     je += pm.L1[k];
     if( pm.PHC[k] != PH_SORPTION )
          continue;

    for( j=jb; j<je; j++ )
    { // Loop for DC
      ja = j - ( pm.Ls - pm.Lads );
      if( pm.lnSAC && ja >= 0 && ja < pm.Lads )
          pm.DUL[j] = pm.lnSAC[ja][3];  // temp. storing initial DUL constr.
    }  // j
  } // k
}
*/
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculating value of dual chemical potential of j-th dependent component
///     performance optimized version  (February 2007)
double TMulti::DC_DualChemicalPotential( double U[], double AL[], long int N, long int j )
{
   long int i, ii;
   double Nu = 0.0;
   for( i = arrL[j]; i < arrL[j+1]; i++ )
   {  ii = arrAN[i];
      if( ii >= N )
       continue;
       Nu += U[ii]*(AL[ii]);
   }
   return Nu;
}
//    for( long int i=0; i<N; i++ )  // obsolete straightforward loop
//         Nu += U[i]*AL[i];
//
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
///  This procedure sets kinetic constraints according to a given
///  concentration units.
//  Needs much more work, elaboration, and performance optimization
//
void TMulti::Set_DC_limits( long int Mode )
{
    double XFL, XFU, XFS=0., XFM, MWXW, MXV, XL=0., XU=0.;
    long int jb, je, j,k, MpL;
    char tbuf[150];

    if( !pm.PLIM )
        return;  // no metastability limits to be set
// ???????????????????????????????????????
    CalculateConcentrations( pm.X, pm.XF, pm.XFA );

    for(k=0; k<pm.FI; k++)
        XFS+=pm.XF[k];  // calculate sum of moles in all phases

    jb=0;
    for( k=0; k<pm.FI; k++ )
    { // cycle over phases
        je=jb+pm.L1[k];
//        XFM=0.;
        MWXW =0.;
        MXV = 0.;
        XFL = 0.;
        XFU = 1e6;
        if( Mode && pm.XF[k] < pm.DSM )
            goto NEXT_PHASE;
        XFM = pm.FWGT[k]; // Mass of a phase
        if( Mode )
        {
            MWXW = XFM/pm.XF[k];         // current molar mass of phase
            MXV = pm.FVOL[k]/pm.XF[k]; // current molar volume of phase
        }
        // Check codes for phase DC
        MpL=0;
        for( j=jb; j<je; j++ )
            if( pm.RLC[j] != NO_LIM )
                MpL = 1;

if( k < pm.FIs )
{					// Temporary workaround - DK  13.12.2007
        if( pm.RFLC[k] == NO_LIM && !MpL )
        { // check type restrictions on phase
            goto NEXT_PHASE;
        }
        switch( pm.RFSC[k] )
        { // check scale restrictions on phase in all system
        case QUAN_MOL:
            XFL = Mode? pm.XF[k]: pm.PLL[k];
            XFU = Mode? pm.XF[k]: pm.PUL[k];
            break;
        case CON_MOLAL:
            XFL = Mode? pm.XF[k]: pm.PLL[k]*pm.GWAT/H2O_mol_to_kg;
            XFU = Mode? pm.XF[k]: pm.PUL[k]*pm.GWAT/H2O_mol_to_kg;
            break;
        case CON_MOLFR:
            XFL = Mode? pm.XF[k]: pm.PLL[k]*XFS;
            XFU = Mode? pm.XF[k]: pm.PUL[k]*XFS;
            break;
        case CON_WTFR:   if(MWXW < 1e-15) break;  // Temp.fix
            XFL = Mode? pm.XF[k]: pm.PLL[k]*pm.MBX/MWXW;
            XFU = Mode? pm.XF[k]: pm.PUL[k]*pm.MBX/MWXW;
            break;
        case CON_VOLFR:   if(MXV < 1e-15) break; // Temp.fix
            XFL = Mode? pm.XF[k]: pm.PLL[k]*pm.VXc/MXV;
            XFU = Mode? pm.XF[k]: pm.PUL[k]*pm.VXc/MXV;
            break;
        default:
            ; // do more?
        }
//        if( pm.RFLC[k] == NO_LIM )
//        {                            Temporary!
            XFL = 0.0;
            XFU = 1e6;
//        }
}
        for( j=jb; j<je; j++ )
        { // loop over DCs
            if( pm.RLC[j] == NO_LIM )
                continue;

            if( Mode )
            {
                XU = pm.DUL[j];
                XL = pm.DLL[j];
            }
            else
                switch( pm.RSC[j] ) // get initial limits
                {
                case QUAN_MOL:
                    XU = pm.DUL[j];
                    XL = pm.DLL[j];
                    break;
                case CON_MOLAL:
                    XU = pm.DUL[j]*pm.GWAT/H2O_mol_to_kg;
                    XL = pm.DLL[j]*pm.GWAT/H2O_mol_to_kg;
                    break;
                case CON_MOLFR:
                    XU = pm.DUL[j]*XFU;
                    XL = pm.DLL[j]*XFL;
                    break;
                case CON_WTFR:
//Ask DK! 20/04/2002
#ifndef IPMGEMPLUGIN
                    XU = pm.DUL[j]*XFU*MWXW /
         TProfil::pm->MolWeight(pm.N, pm.Awt, pm.A+j*pm.N );
                    XL = pm.DLL[j]*XFL*MWXW /
         TProfil::pm->MolWeight(pm.N, pm.Awt, pm.A+j*pm.N );

#endif
                    break;
                case CON_VOLFR:
                    XU = pm.DUL[j]*XFU*MXV/ pm.Vol[j];
                    XL = pm.DLL[j]*XFL*MXV/ pm.Vol[j];
                    break;
                default:
                    ; // do more
                }
            // check combine
            if( XU < 0.0 ) XU = 0.0;
            if( XU > 1e6 ) XU = 1e6;
            if( XL < 0.0 ) XL = 0.0;
            if( XL > 1e6 ) XL = 1e6;
            if( XU > XFU )
            {
 //               JJ = j;
//                KK = k;
                sprintf( tbuf, "Inconsistent upper DC metastability limits j=%ld k=%ld XU=%g XFU=%g",
                         j, k, XU, XFU );
                Error( "E11IPM: Set_DC_limits(): ",tbuf );
//                XU = XFU; // - pm.lowPosNum;
            }
            if( XL < XFL )
            {
//                JJ = j;
//                KK = k;
                sprintf( tbuf, "Inconsistent lower DC metastability limits j=%ld k=%ld XL=%g XFL=%g",
                         j, k, XL, XFL );
                Error( "E12IPM: Set_DC_limits(): ",tbuf );
//                XL = XFL; // - pm.lowPosNum;
            }
            pm.DUL[j]=XU;
            pm.DLL[j]=XL;
        }   // j
NEXT_PHASE:
        jb = je;
    }  // k
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculating total amounts of phases
//
void TMulti::TotalPhasesAmounts( double X[], double XF[], double XFA[] )
{
    long int jj, j, i, k;
    double XFw, XFs, x;

    j=0;
    for( k=0; k< pm.FI; k++ )
    { // cycle by phases
        i=j+pm.L1[k];
        XFw = 0.0;
        XFs=0.0; // calculating mole amount of carrier (solvent/sorbent)
        for(jj=j; jj<i; jj++)
        {
            x = X[jj];
            if( pm.DCCW[jj] == DC_ASYM_CARRIER && pm.FIs )
                XFw += x;
            else XFs += x;
        }
        XF[k] = XFw + XFs;
        if( k < pm.FIs )
            XFA[k] = XFw;
        j=i;
    }

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Corrections to primal chemical potentials F0[j]
///  of j-th species in k-th phase among IPM main iterations.
///  Returns double value of corrected chem. potential.
///  If error, returns +7777777 J/mole.
//  Last modif. 05 Jan 2000 by DK to include BSM EDL model.
double TMulti::DC_PrimalChemicalPotentialUpdate( long int j, long int k )
{
    long int ja=0, ist, isp, jc=-1;
    double F0=0.0, Fold, dF0, Mk=0.0, Ez, psiA, psiB, CD0, CDb, ObS;
    double FactSur, FactSurT;
    SPP_SETTING *pa = paTProfil;

    Fold = pm.F0[j];
    if( pm.FIat > 0 && j < pm.Ls && j >= pm.Ls - pm.Lads )
    {
        ja = j - ( pm.Ls - pm.Lads );
        jc = pm.SATX[ja][XL_EM];
// Workaround 08.02.2011 by DK to prevent crash on inconsistent sorption DC
        if( jc >= pm.L1[k] )
            jc = -1;
    }
    if( k < pm.FIs && pm.XFA[k] > 1e-12)
    {
           if( jc < 0 ) // phase (carrier) molar mass g/mkmol
              Mk = pm.FWGT[k]/pm.XFA[k]*1e-6;
           else Mk = pm.MM[jc]*(pm.X[jc]/pm.XFA[k])*1e-6;
        // DC carrier molar mass g/mkmol
    }
    switch( pm.DCC[j] )
    { // analyse species class code
    case DC_SCP_CONDEN:
    case DC_AQ_PROTON:
    case DC_AQ_ELECTRON:
    case DC_AQ_SPECIES:
case DC_AQ_SURCOMP:
    case DC_GAS_COMP:
    case DC_GAS_H2O:
    case DC_GAS_CO2:
    case DC_GAS_H2:
    case DC_GAS_N2:
    case DC_AQ_SOLVENT:
    case DC_AQ_SOLVCOM:
    case DC_SOL_IDEAL:
    case DC_SOL_MINOR:
    case DC_SOL_MAJOR:
    case DC_SOL_MINDEP:
    case DC_SOL_MAJDEP:
        F0 = pm.lnGmM[j];
        break;
        // adsorption
    case DC_SSC_A0:
    case DC_SSC_A1:
    case DC_SSC_A2:
    case DC_SSC_A3:
    case DC_SSC_A4:
    case DC_WSC_A0:
    case DC_WSC_A1:
    case DC_WSC_A2:
    case DC_WSC_A3:
    case DC_WSC_A4:
    case DC_SUR_GROUP:
    case DC_SUR_COMPLEX:
    case DC_SUR_IPAIR:
    case DC_IESC_A:
    case DC_IEWC_B:
        if( !pm.Lads || !pm.FIat )   // Foolproof - to prevent crash DK 22.12.2009
            break;
        F0 = pm.lnGmM[j]; /* + pm.lnGam[j]; */
        // get ist - index of surface type and isp - index of surface plane
/*!!!!!*/  ist = pm.SATX[ja][XL_ST];  // / MSPN;
/*!!!!!*/  isp = pm.SATX[ja][XL_SP]; // % MSPN;
        CD0 = pm.MASDJ[ja][PI_CD0];  // species charge that goes into 0 plane
        CDb = pm.MASDJ[ja][PI_CDB];  // species charge that goes into B plane
        ObS = pm.MASDJ[ja][PI_DEN];  // obsolete - the sign for outer-sphere charge
        if( ObS >= 0.0 )
            ObS = 1.0;
        else ObS = -1.0;
        psiA = pm.XpsiA[k][ist];
        psiB = pm.XpsiB[k][ist];
        Ez = double(pm.EZ[j]);
        if( !isp )  // This is the 0 (A) plane species
        {
            if( fabs( CD0 ) > 1e-20 )  // Doubtful...
                Ez = CD0;
            F0 += psiA * Ez * pm.FRT;
        }
        else  // This is B plane
        {
            if( pm.SCM[k][ist] == SC_MTL || pm.SCM[k][ist] == SC_MXC )
            { // Modified TL: Robertson, 1997; also XTLM Kulik 2002
                  if( fabs( CDb ) > 1e-20 )  // Doubtful...
                      Ez = CDb;
                  F0 += psiB * Ez * pm.FRT;
            }
            if( pm.SCM[k][ist] == SC_TLM )
            {
// New CD version of TLM  added 25.10.2004
               if( fabs( CD0 ) > 1e-20 && fabs( CDb ) > 1e-20 )
                  F0 += ( psiA*CD0 + psiB*CDb )* pm.FRT;
             // see also Table 4 in Zachara & Westall, 1999
             // Old version:  TLM Hayes & Leckie, 1987 uses the sign indicator at density
                else {
                  if( ObS < 0 )
                  {
                      Ez -= 1.0;
                      F0 += ( psiA + Ez * psiB )* pm.FRT;
                  }
                  else
                  {
                      Ez += 1.0;
                      F0 += ( Ez * psiB - psiA )* pm.FRT;
                  }
               }
            }
            else if( pm.SCM[k][ist] == SC_BSM || pm.SCM[k][ist] == SC_CCM )
            { // Basic Stern model, Christl & Kretzschmar, 1999
            // New CD version of TLM  added 25.10.2004
               if( fabs( CD0 ) > 1e-20 && fabs( CDb ) > 1e-20 )
                  F0 += ( psiA*CD0 + psiB*CDb )* pm.FRT;
                else {
                  if( ObS < 0 )
                  {
                      Ez -= 1.0;
                      F0 += ( psiA + Ez * psiB )* pm.FRT;
                  }
                  else
                  {
                      Ez += 1.0;
                      F0 += ( Ez * psiB - psiA )* pm.FRT;
                  }
               }
            }
        }
        if( Mk > 1e-9 )  // Mk is carrier molar mass in g/mkmol
        {   // Correction for standard density, surface area and surface type fraction
                FactSur = Mk * (pm.Aalp[k]) * pa->p.DNS*1.66054;
        	    // FactSur is adsorbed mole amount at st. surf. density per mole of solid carrier
                FactSurT = FactSur * (pm.Nfsp[k][ist]);
                if( pm.SCM[k][ist] == SC_MXC || pm.SCM[k][ist] == SC_NNE ||
                    pm.SCM[k][ist] == SC_IEV )
                                                // F0 -= log( Mk * (pm.Nfsp[k][ist]) *
                                                // (pm.Aalp[k]) * pa->p.DNS*1.66054 );
                  F0 -= log( FactSurT );
            else  F0 -= log( FactSurT );
                                // F0 -= log( Mk * (pm.Nfsp[k][ist]) *
                                // (pm.Aalp[k]) * pa->p.DNS*1.66054 );
            F0 -= FactSur / ( 1. + FactSur );
                                // F0 -= (pm.Aalp[k])*Mk*pa->p.DNS*1.66054 /
                                // ( 1.0 + (pm.Aalp[k])*Mk*pa->p.DNS*1.66054 );
        }
        break;
    case DC_PEL_CARRIER:
    case DC_SUR_MINAL:  // constant charge of carrier - not completed
    case DC_SUR_CARRIER: // Mk is carrier molar mass in g/mkmol
        if( !pm.Lads || !pm.FIat )   // Foolproof - to prevent crash DK 22.12.2009
            break;
        FactSur = Mk * (pm.Aalp[k]) * pa->p.DNS*1.66054;
        F0 -= FactSur / ( 1. + FactSur );
        F0 += FactSur;
        break;
    }
    F0 += pm.lnGam[j];

    if( k >= pm.FIs )
        return F0;
    // Smoothing procedure for highly non-ideal systems
    if( pm.sMod[k][SGM_MODE] != SM_IDEAL )  // check this condition for sublattice SS models!
            // || pm.sMod[k][SCM_TYPE] != SC_NNE )  // changed, 14.07.2009 (TW)
    {
        double SmoSensT = 1e-5;   // to be adjusted
        dF0 = F0 - Fold;
        if( pm.X[j] > min( pm.lowPosNum, pm.DcMinM ) && fabs( dF0 ) >= SmoSensT )
       	    F0 = Fold + dF0 * SmoothingFactor();    // Changed 18.06.2008 DK
    }
    return F0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of DC primal chemical potential F.
/// From moles of DC Y[], total moles of phase YF[] and DC partial
/// molar Gibbs energy gT (obtained from pm.G[]) which includes
/// activity coefficient terms.
/// On error returns F = +7777777.
double TMulti::DC_PrimalChemicalPotential(
    double G,      ///< gT0+gEx
    double logY,   ///< ln x
    double logYF,  ///< ln Xa
    double asTail, ///< asymmetry non-log term or 0 for symmetric phases
    double logYw,  ///< ln Xv
    char DCCW      ///< generalized species class code
)
{
    double F;
    switch( DCCW )
    {
    case DC_SINGLE:
        F = G;
        break;
    case DC_ASYM_SPECIES:
        F = G + logY - logYw + asTail;
        break;
    case DC_ASYM_CARRIER:
        F = G + logY - logYF + asTail + 1.0 -
            1.0/(1.0 - asTail);
        break;
    case DC_SYMMETRIC:
        F = G + logY - logYF;
        break;
    default:
        F = 7777777.;
    }
    return F;
}


// Kernel functions of IPM - rewritten by DK for adsorption
//
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// VJ - Update of primal chemical potentials
//
void
TMulti::PrimalChemicalPotentials( double F[], double Y[], double YF[], double YFA[] )
{
    long int i,j,k;
    double NonLogTerm=0., v, Yf; // v is debug variable

    for( j=0; j<pm.L; j++)
       F[j] =0;

    j=0;
    for( k=0; k<pm.FI; k++ )
    { // loop over phases
        i=j+pm.L1[k];
//        if( YF[k] <= pm.lowPosNum*100. || ( pm.PHC[k] == PH_AQUEL &&
//        ( YF[k] <= TProfil::pm->pa.p.XwMin || Y[pm.LO] <= pm.lowPosNum*1e3 )))
        if( pm.L1[k] == 1L && YF[k] < pm.PhMinM )
        	goto NEXT_PHASE;
        if( YF[k] <= pm.DSM || ( pm.PHC[k] == PH_AQUEL &&
            ( YF[k] <= pm.DSM || Y[pm.LO] <= pm.XwMinM )))
            goto NEXT_PHASE;

        pm.YFk = 0.0;
        Yf= YF[k]; // calculate number of moles of carrier
        if( pm.FIs && k<pm.FIs )
            pm.YFk = YFA[k];
        if( Yf >= 1e6 )
        {                 // error - will result in zerodivide!
           gstring pbuf(pm.SF[k],0,20);
           char buf[200];
           sprintf( buf, "Broken phase amount from primal approximation: Phase %s  Yf= %lg", pbuf.c_str(), Yf );
           Error( "E13IPM: PrimalChemicalPotentials():", buf);
//           Yf = pm.YFk;
        }
//        if( pm.YFk > pm.lowPosNum*10. )
        if( ( pm.PHC[k] == PH_AQUEL && pm.YFk >= pm.XwMinM )
                        || ( pm.PHC[k] == PH_SORPTION && pm.YFk >= pm.ScMinM )
                        || ( pm.PHC[k] == PH_POLYEL && pm.YFk >= pm.ScMinM ) )
        {
            pm.logXw = log(pm.YFk);
            NonLogTerm = 1.- pm.YFk / Yf;
#ifdef NOMUPNONLOGTERM
NonLogTerm = 0.0;
#endif
        }
        if( pm.L1[k] > 1 )
        {
            pm.logYFk = log( Yf );
        }
        if( pm.PHC[k] == PH_AQUEL )
        {    // ln moles of solvent in aqueous phase
            pm.Yw = pm.YFk;
            pm.aqsTail = NonLogTerm;
        }
        for( ; j<i; j++ )
        { //  cycle by DC
            if( Y[j] < min( pm.DcMinM, pm.lowPosNum ))
                continue;  // exception by minimum DC quantity
                           // calculate chemical potential of j-th DC
            v = DC_PrimalChemicalPotential( pm.G[j], log(Y[j]), pm.logYFk,
                              NonLogTerm, pm.logXw, pm.DCCW[j] );
            F[j] = v;
       }   // j
NEXT_PHASE:
        j = i;
    }  // k
    if( pm.Yw >= pm.DSM ) // pm.lowPosNum*1e3 )
        pm.logXw = log(pm.Yw);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of a species contribution to the total Gibbs energy G(X)
/// of the system. On error returns +7777777.
//
double TMulti::DC_GibbsEnergyContribution(
    double G,      ///< gT0+gEx
    double x,      ///< x - mole amount of species
    double logXF,  ///< ln Xa - mole amount of phase
    double logXw,  ///< ln Xv - mole amount of the solvent/sorbent
    char DCCW      /// generalized species class code
)
{
    double Gi;

    switch( DCCW )
    {
    case DC_ASYM_SPECIES:
        Gi = x * ( G + log(x) - logXw );
        break;
    case DC_ASYM_CARRIER:
    case DC_SYMMETRIC:
        Gi = x * ( G + log(x) - logXF );
        break;
    case DC_SINGLE:
        Gi = G * x;
        break;
    default:
        Gi = 7777777.;
    }
    return Gi;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of the total Gibbs energy of the system G(X)
/// and copying of Y, YF vectors into X,XF, respectively.
///  Parameter LM is the IPM step size for calculation of new
///  quantities of all species (vector X[]) using the direction
///  of descent (MU[] vector). If LM == 0, this function
///  just copies vector Y[] into X[].
///  \return of G(X) in moles.
//
double TMulti::GX( double LM  )
{
    long int i, j, k;
    double x, XF, XFw, FX, Gi; // debug variable
//    double const1= pm.lowPosNum*10.,
//           const2 = pm.lowPosNum*1000.;

    if( LM < pm.lowPosNum )     // copy vector Y into X
        for(i=0;i<pm.L;i++)
            pm.X[i]=pm.Y[i];
    else  // calculate new values of X
        for(i=0;i<pm.L;i++ )
        {  // gradient vector pm.MU - the direction of descent!
            pm.X[i]=pm.Y[i]+LM*pm.MU[i];
//            if( pm.X[i] <  pm.lowPosNum )   // this is the Ls set cutoff !!!!!!!!!!
            if( pm.X[i] <  pm.DcMinM )
                pm.X[i]=0.;
        }
    // calculate new total quantities of phases
    TotalPhasesAmounts( pm.X, pm.XF, pm.XFA );

    // calculating G(X)
    FX=0.;
    j=0;
    for( k=0; k<pm.FI; k++ )
    { // loop for phases
        i=j+pm.L1[k];
        pm.logXw = -101.;
        XFw = 0.0;  // calculating mole amount of the solvent/sorbent
        if( pm.FIs && k<pm.FIs )
            XFw = pm.XFA[k];
 //       if( XFw > const1 )
        if( ( pm.PHC[k] == PH_AQUEL && XFw >= pm.XwMinM )
                        || ( pm.PHC[k] == PH_SORPTION && XFw >= pm.ScMinM )
                        || ( pm.PHC[k] == PH_POLYEL && XFw >= pm.ScMinM ) )
             pm.logXw = log( XFw );
        /*   */
        XF = pm.XF[k];
        if( !(pm.FIs && k < pm.FIs) )
        {
                if( XF < pm.PhMinM )
        		goto NEXT_PHASE;
        }
        else if( XF < pm.DSM && pm.logXw < -100. )
        	goto NEXT_PHASE;

        pm.logYFk = log( XF );

        for( ; j<i; j++ )
        { // DCs (species)
            x = pm.X[j];
            if( x < pm.DcMinM )
                continue;
            // calculating increment of G(x)
            // Gi = DC_GibbsEnergyContribution( pm.G[j], x, pm.logYFk, pm.logXw,
            //                     pm.DCCW[j] );
            // call replaced here by inline variant for higher performance
            switch( pm.DCCW[j] )
            {
             case DC_ASYM_SPECIES:
                    Gi = x * ( pm.G[j] + log(x) - pm.logXw );
                    break;
            case DC_ASYM_CARRIER:
            case DC_SYMMETRIC:
                   Gi = x * ( pm.G[j] + log(x) - pm.logYFk );
                   break;
            case DC_SINGLE:
                   Gi = pm.G[j] * x;
                   break;
           default:
                    Gi = 7777777.;
           }
          FX += Gi;
        }   // j
NEXT_PHASE:
        j = i;
    }  // k
//    cout << "GX  " << setprecision(16) << scientific <<  FX << endl;
    return(FX);
}

#ifndef IPMGEMPLUGIN
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Variant of GX() function for use in the UnSpace module (non-optimized)
// Should not be called from within GEMIPM!
//
double TMulti::pb_GX( double *Gxx  )
{
    long int i, j, k;
    double Gi, x, XF, XFw, FX;
    SPP_SETTING *pa = paTProfil;

    // calculating G(X)
    FX=0.;
    j=0;
    for( k=0; k<pm.FI; k++ )
    { // phase loop
        i=j+pm.L1[k];
        pm.logXw = -101.;
        XFw = 0.0;  // calculating mole amount of the solvent/sorbent
        if( pm.FIs && k<pm.FIs )
            XFw = pm.XFA[k];
 //       if( XFw > const1 )
        if( ( pm.PHC[k] == PH_AQUEL && XFw >= pa->p.XwMin )
                        || ( pm.PHC[k] == PH_SORPTION && XFw >= pa->p.ScMin )
                        || ( pm.PHC[k] == PH_POLYEL && XFw >= pa->p.ScMin ) )
             pm.logXw = log( XFw );
        /*   */
        XF = pm.XF[k];
        if( !(pm.FIs && k < pm.FIs) )
        {
                if( XF < pa->p.PhMin )
        		goto NEXT_PHASE;
        }
        else if( XF < pa->p.DS && pm.logXw < 100. )
        	goto NEXT_PHASE;
        pm.logYFk = log( XF );

        for( ; j<i; j++ )
        { // DC loop
            x = pm.X[j];
            if( x < pa->p.DcMin )
                continue;
            // calculating DC increment to G(x)
            Gi = DC_GibbsEnergyContribution( Gxx[j], x, pm.logYFk, pm.logXw,
                                 pm.DCCW[j] );
            FX += Gi;
        }   // j
NEXT_PHASE:
        j = i;
    }  // k
    return(FX);
}
#endif


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Conversion of g(T,P) value for DCs into the uniform cj scale.
/// \param k - index of phase, \param j - index DC in phase
/// \return if error code, returns 777777777.
//
double TMulti:: ConvertGj_toUniformStandardState( double g0, long int j, long int k )
{
    double G, YOF=0;

    G = g0/pm.RT;
    if( pm.YOF )
        YOF = pm.YOF[k];     // should be already normalized (J/mol/RT)
    // Calculation of standard concentration scaling terms
    switch( pm.DCC[j] )
    { // Aqueous electrolyte
    case DC_AQ_PROTON:
    case DC_AQ_ELECTRON:
    case DC_AQ_SPECIES:
case DC_AQ_SURCOMP:
        G += pm.ln5551;
        // calculate molar mass of solvent
    case DC_AQ_SOLVCOM:
    case DC_AQ_SOLVENT:
#ifndef IPMGEMPLUGIN
        if( TSyst::sm->GetSY()->PYOF != S_OFF )
#endif
          if( YOF != 0.0 )
        	G += YOF;  // In GEMS3K, YOF[k] is the only way to influence G directly

        break;
    case DC_GAS_COMP: // gases except H2O and CO2
    case DC_GAS_H2O: // index to switch off?
    case DC_GAS_CO2:
    case DC_GAS_H2:
    case DC_GAS_N2:
    case DC_SOL_IDEAL:
    case DC_SOL_MINOR:
    case DC_SOL_MAJOR: // changed by DK on 4.12.2006
    case DC_SOL_MINDEP:
    case DC_SOL_MAJDEP:
        if( pm.PHC[k] == PH_GASMIX || pm.PHC[k] == PH_FLUID
            || pm.PHC[k] == PH_PLASMA )
        {
//        if( pm.Pparc[j] != 1.0 && pm.Pparc[j] > 1e-30 )
//           G += log( pm.Pparc[j] ); // log partial pressure/fugacity
//        else
               G += log( pm.Pc ); // log general pressure (changed 04.12.2006)
        }
        // non-electrolyte condensed mixtures
    case DC_SCP_CONDEN: // single-component phase
    case DC_SUR_MINAL:
    case DC_SUR_CARRIER:
    case DC_PEL_CARRIER:
#ifndef IPMGEMPLUGIN
        if( TSyst::sm->GetSY()->PYOF != S_OFF )
#endif
          if( YOF != 0.0 )
        	 G += YOF;
        break;
        // Sorption phases
    case DC_SSC_A0:
    case DC_SSC_A1:
    case DC_SSC_A2:
    case DC_SSC_A3:
    case DC_SSC_A4:
    case DC_WSC_A0:
    case DC_WSC_A1:
    case DC_WSC_A2:
    case DC_WSC_A3:
    case DC_WSC_A4:
    case DC_SUR_GROUP:
    case DC_SUR_COMPLEX:
    case DC_SUR_IPAIR:
    case DC_IESC_A:
    case DC_IEWC_B:
        G += pm.ln5551;
        break;
    default: // error - returning 7777777
        return 7777777.;
    }
    return G;
}

/// Converting DC class codes into generic internal codes of IPM
//
void TMulti::ConvertDCC()
{
    long int i, j, k, iRet=0;
    char DCCW;

    j=0;
    for( k=0; k< pm.FI; k++ )
    { // phase loop
        i=j+pm.L1[k];
        if( pm.L1[k] == 1 )
        {
            pm.DCCW[j] = DC_SINGLE;
            goto NEXT_PHASE;
        }
        for( ; j<i; j++ )
        { // DC loop
            switch( pm.DCC[j] ) // select v_j expression
            {
            case DC_SCP_CONDEN:
                DCCW = DC_SINGLE;
                break;
            case DC_GAS_COMP:
            case DC_GAS_H2O:
            case DC_GAS_CO2:
            case DC_GAS_H2:
            case DC_GAS_N2:
            case DC_SOL_IDEAL:
            case DC_SOL_MINOR:
            case DC_SOL_MAJOR:
            case DC_SOL_MINDEP:
            case DC_SOL_MAJDEP:
                DCCW = DC_SYMMETRIC;
                break;
            case DC_AQ_PROTON:
            case DC_AQ_ELECTRON:
            case DC_AQ_SPECIES:
            case DC_AQ_SURCOMP:
                DCCW = DC_ASYM_SPECIES;
                break;
            case DC_AQ_SOLVCOM:
            case DC_AQ_SOLVENT:
                DCCW = DC_ASYM_CARRIER;
                break;
            case DC_IESC_A:
            case DC_IEWC_B:
                DCCW = DC_ASYM_SPECIES;
                break;
                // Remapping
            case DC_SUR_GROUP:
            case DC_SUR_COMPLEX:
                DCCW = DC_ASYM_SPECIES;
                pm.DCC[j] = DC_SSC_A0;
                break;
            case DC_SUR_IPAIR:
                DCCW = DC_ASYM_SPECIES;
                pm.DCC[j] = DC_WSC_A0;
                break;
            case DC_SUR_MINAL:
            case DC_SUR_CARRIER:
            case DC_PEL_CARRIER:
                 DCCW = DC_ASYM_CARRIER;
                break;
            default:
                if( isdigit( pm.DCC[j] ))
                {
                    if( pm.PHC[k] == PH_SORPTION ||
                            pm.PHC[k] == PH_POLYEL )
                    {
                        DCCW = DC_ASYM_SPECIES;
                        break;
                    }
                }
                DCCW = DC_SINGLE;
                iRet++;  // error the class code
            }
            pm.DCCW[j] = DCCW;
        }   // j
NEXT_PHASE:
        j = i;
    }  // k
    ErrorIf( iRet>0, "E19IPM: ConvertDCC()", "Invalid DC class code. Memory corruption?");
}

/// Get the index of volume IC ("Vv") for the volume balance constraint
long int TMulti::getXvolume()
{
 long int ii, ret = -1;
 for( ii = pm.N-1; ii>=0; ii--)
 {
  if( pm.ICC[ii] == IC_VOLUME )
  { ret = ii; break; }
 }
 return ret;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of Karpov stability criteria for a DC
// Modified for kinetic constraints 05.11.2007, 16.04.2012 by DK
//
double TMulti::KarpovCriterionDC(
    double *dNuG,  // logarithmic Nu[j]-c[j] difference - is modified here
    double logYF,  // ln Xa   (Xa is mole amount of the whole phase)
    double asTail, // asymmetry correction (0 for symmetric phases)
    double logYw,  // ln Xw   (Xw is mole amount of solvent)
    double Wx,     // primal mole fraction of this DC
    char DCCW      // Generic class code of DC
)
{
    double Fj = 0.0;  // output phase stability criterion

    if( logYF > -76. && Wx > 1e-33 )    // Check thresholds!
        switch( DCCW ) // expressions for fj
        {
        default: // error code would be needed here !!!
            *dNuG = -228.;
        case DC_SINGLE:
//            Wx = 1.0;
        case DC_SYMMETRIC:
            break;
        case DC_ASYM_SPECIES:
            *dNuG += logYw - logYF - asTail;
            break;
        case DC_ASYM_CARRIER:
            *dNuG += 1.0/(1.0 - asTail) - asTail - 1.0;
        }    
    if( *dNuG < 13.8155 && *dNuG > -228. )
        Fj = exp( *dNuG );  // dual estimate of mole fraction
    if( *dNuG >= 13.8155 )
        *dNuG = 13.8155;
    if( *dNuG <= -228. )
        *dNuG = -228.;
    Fj -= Wx;     // If primal Wx = 0 then this DC is not in L_S set
    return Fj;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of Karpov stability criteria for all phases
// Modifications for metastability-controlled components: DK 16.04.2012
//
void TMulti::KarpovsPhaseStabilityCriteria()
{
    bool KinConstr, fRestore;
    long int k, j, ii;
    double *EMU,*NMU, YF, Nu, dNuG, Wx, Yj, Fj, sumWx, NonLogTerm = 0.;
    SPP_SETTING *pa = paTProfil;

    EMU = pm.EMU;
    NMU = pm.NMU;
    for(ii=0; ii<pm.L; ii++ )
        EMU[ii] = NMU[ii]=0.0;
    j=0;
    pm.YMET = 0.0;  // check setup!
    for( k=0; k<pm.FI; k++ )
    { // phases
        ii=j+pm.L1[k];
        pm.Falp[k] = pm.YMET; // metastability parameter
        pm.logXw = -76.;
        pm.logYFk = -76.;

        pm.YFk = 0.0;
        YF= pm.YF[k]; // moles of carrier
        if( pm.FIs && k<pm.FIs )
            pm.YFk = pm.YFA[k];
//        if( pm.PHC[k] == PH_AQUEL ) {
            if( pm.YFk > 1e-33 )   // amount of phase or carrier cannot be less than 1e-33 mol!
            {
               pm.logXw = log(pm.YFk);
               NonLogTerm = 1.- pm.YFk / YF;
#ifdef NOMUPNONLOGTERM
NonLogTerm = 0.0;
#endif
            }
            else {
               pm.logXw = -76.;
               NonLogTerm = 0.0;
            }
//        }
        if( pm.L1[k] > 1 && YF > 1e-33 )
            pm.logYFk = log( YF );
        else pm.logYFk = -76.;
        if( pm.PHC[k] == PH_AQUEL) // number of moles of solvent
        {
            pm.Yw = pm.YFk;
            pm.aqsTail = NonLogTerm;
        }
// The code below streamlined on 16.04.2012 by DK
        sumWx = 0.0; fRestore = false;
        for( ; j<ii; j++ )
        {
            KinConstr = false;
            Nu = DC_DualChemicalPotential( pm.U, pm.A+j*pm.N, pm.NR, j );
            dNuG = Nu - pm.G[j]; // this is -s_j (6pot paper 1)
            Wx = 0.0;
            Yj = pm.Y[j];
if( YF > pm.DSM )  // if-else rearranged 24.08.2012 DK
{
    if( Yj > pm.DcMinM )
        Wx = Yj / YF; // calculating primal mole fraction of DC
}
else  { // if phase is removed then a saved copy of activity coefficients is used
        //   dNuG += pm.fDQF[j];   // bugfix 28.08.2012 DK
            if( pm.K2 && pm.lnGam[j] == 0. && pm.GamFs[j] != 0 )
                dNuG -= pm.GamFs[j];
            if( pm.GamFs[j] != 0 )
            {
                fRestore = true;
                Wx = pow( 10., (pm.Y_la[j] - pm.GamFs[j]/lg_to_ln ));
            }
            if( pm.L1[k] > 1 && pm.sMod )
            {
                if( pm.sMod[k][SPHAS_TYP] == SM_IDEAL )
                {
                   Wx = pow( 10., pm.Y_la[j] );
                   fRestore = true; // can always insert a simple ideal solution phase
                }
            }
      }
            if( ( pm.DUL[j] < 1e6 && Yj >= ( pm.DUL[j] - pa->p.DKIN ) )
                || ( pm.DLL[j] > 0 && Yj <= ( pm.DLL[j] + pa->p.DKIN ) ) )
                KinConstr = true; // DC with the amount lying on the non-trivial kinetic constraint
            // calculating Karpov stability criteria for DCs
            Fj = KarpovCriterionDC( &dNuG, pm.logYFk, NonLogTerm,
                     pm.logXw, Wx, pm.DCCW[j] );
            NMU[j] = dNuG;  // dNuG is stored for all DCs, not only those in L_S set
            EMU[j] = Fj + Wx; // fix dual mole fraction estimate DK 13.04.2012
            if( KinConstr == false )
                sumWx += EMU[j];
            else
                sumWx += Wx; // under active kinetic constraint, only primal mole fraction is valid!
        } // j
        j = ii;
//        if( pm.L1[k] == 1 )
        pm.Falp[k] = sumWx - 1.;  // generalized Karpov critetion with metastable components
        if( fabs( pm.Falp[k] ) < pm.DSM )
            pm.Falp[k] = 0.;
if( pm.L1[k] > 1 && YF <= pm.DSM && fRestore == false )
    pm.Falp[k] = -1.;   // provisional - set to
    }  // k
}

//===================================================================
/// Speciation cleanup subroutine if CleanupStatus is not zero
/// \param  AmountCorrectionThreshold - the maximum DC amount correction that can be cleaned ( 1e-5 )
/// \param  MjuDiffCutoff - normalized chem.pot. difference threshold (dMu = ln a - ln a,dual)
//
/// \return 0 if no subsequent refinement of mass balance is needed;
///         1 if species amounts were cleaned up to more than requested overall mass balance accuracy
///        -1 if cleanup has been done and the degeneration of the chemical system occurred
///         2 solution is seriously distorted and full PhaseSelect3() loop is necessary
///
long int TMulti::SpeciationCleanup( double AmountCorrectionThreshold, double MjuDiffCutoff )
{
    long int NeedToImproveMassBalance = 0, L1k, L1kZeroDCs, k, j, jb = 0;
    double MjuPrimal, MjuDual, MjuDiff, Yj, YjDiff=0., YjCleaned;
    double CutoffDistortionMBR = 0.1 * pm.DHBM;
    bool KinConstr, Degenerated = false;
    SPP_SETTING *pa = paTProfil;

    PrimalChemicalPotentials( pm.F, pm.Y, pm.YF, pm.YFA );
    jb=0;
    for(k=0;k<pm.FI;k++)
    {
       if( ( pm.YF[k] >= pm.DcMinM ) ) // Only in phase present in mass balance!
       {                            // (acc. to definition of the L_S set)
            L1k = pm.L1[k]; // Number of components in the phase
            L1kZeroDCs = 0;
            for(j=jb; j<jb+L1k; j++)
            {
               Yj = YjCleaned = pm.Y[j];
               KinConstr = false;
                // Detecting the DC having the non-trivial kinetic constraint
               // Fixing a very small component constrained from below
               if( pm.DUL[j] < 1e6 && Yj >= pm.DUL[j] ) // Yj >= ( pm.DUL[j] - pa->p.DKIN ) )
               {
                   pm.Y[j] = pm.DUL[j];
                   KinConstr = true;
               }
               if( pm.DLL[j] > 0 && Yj <= pm.DLL[j] )   // Yj <=( pm.DLL[j] + pa->p.DKIN ) )
               { // Fixing a small component constrained from above
                   pm.Y[j] = pm.DLL[j];
                   KinConstr = true;
               }
               if( KinConstr == true )
               {
                   YjDiff = fabs( pm.Y[j] -Yj );
                   if( YjDiff > CutoffDistortionMBR )
                       NeedToImproveMassBalance = 1;
                   if( YjDiff > AmountCorrectionThreshold )
                       NeedToImproveMassBalance = 2;
                   continue;   // skipping DC if it has active non-trivial metastability constraint
               }
//               else
               if( Yj >= pm.DcMinM )
               {   // we check in the Ls set only, except metastability constraints
                  if( ( pm.DUL[j] < 1e6 && Yj > pm.DUL[j] - pa->p.DKIN )
                      || ( pm.DLL[j] > 0 && Yj < pm.DLL[j] + pa->p.DKIN ) )
                      continue; // we don't clean xDC sitting too close to metast. constraints
                  MjuPrimal = pm.F[j];   // normalized
                  MjuDual = pm.Fx[j];    // /pmp->RT; Fixed DK 12.03.2012
                  MjuDiff = MjuPrimal - MjuDual;
                  if( fabs( MjuDiff ) > MjuDiffCutoff )
                  {
                      if( L1k == 1 && MjuDiff > 0. )
                      {  // Pure phase
                         YjCleaned = 0.0; // Cleaning out a "phantom" pure phase
                      }
                      if(L1k > 1)
                      {  // Component of a solution phase
                         // The species is present in a larger or smaller amount than necessary
                           YjCleaned = Yj / exp( MjuDiff );
                      }
                  }
               }
               else if( !(pm.PHC[k] == PH_SORPTION || pm.PHC[k] == PH_POLYEL)  )
               {  // component has been zeroed off - primal chemical potential not defined
                    if( L1k == 1 )
                       continue; // L1kZeroDCs++; No way to improve for pure substance
                    else  // possibly lost DC in a solution phase (only estimated mole fraction available)
                       YjCleaned = pm.EMU[j]*pm.YF[k]; // does not work for adsorption yet
               }

               YjDiff = YjCleaned - Yj;
               if( fabs( YjDiff ) > CutoffDistortionMBR )
               {
                   NeedToImproveMassBalance = 1;
                   if( fabs( YjDiff ) > AmountCorrectionThreshold )
                   {   // Correction was too large
                       NeedToImproveMassBalance = 2;
                       // Temporary: only correction of the size of threshold
                       if( YjDiff > 0. )
                           pm.Y[j] += AmountCorrectionThreshold;
                       else
                           pm.Y[j] -= AmountCorrectionThreshold;
                       }
                       else {  // Reasonable correction
                           pm.Y[j] = YjCleaned;
                       }
                   }
               else {  // Correction that does not affect the mass balance
                          pm.Y[j] = YjCleaned;
               }
               if( pm.Y[j] < pm.DcMinM )
               {  // Corrected amount is too small - DC amount is zeroed off
                   pm.Y[j] = 0.;
//                   DCremoved++;
                   L1kZeroDCs++;
               }
            }  // for j
            if(( pm.L1[k] - L1kZeroDCs <= 1 && k < pm.FIs )
                ||( pm.L1[k] - L1kZeroDCs == 0 && k >= pm.FIs ))
            {   Degenerated = true;
                NeedToImproveMassBalance = 1;
            }
       }
       jb+=pm.L1[k];
   }
   if( NeedToImproveMassBalance )
   { // diagnostics to be implemented
      if( Degenerated && NeedToImproveMassBalance == 1 )
          NeedToImproveMassBalance = -1;
      if( Degenerated && NeedToImproveMassBalance == 2 )
      {  // Diagnostic output here
          NeedToImproveMassBalance = -2;
      }
   }
   return NeedToImproveMassBalance;
}

//====================================================================================
/// New simplified PSSC() algorithm   DK 01.05.2010.
/// PhaseSelection() part only looks for phases to be inserted, also checks if some
/// solution phases are unstable. Removal of unstable phases is done afterwards in
/// SpeciationCleanup subroutine if CleanupStatus is not zero.
/// As phase stability criterion, uses (log) phase stability (saturation) index
/// computed from DualTh activities of components and activity coefficients
/// \return 1L if Ok; 0 if one more IPM loop should be done;
///          -1L if 3 loops did not fix the problem
//
long int TMulti::PhaseSelectionSpeciationCleanup( long int &kfr, long int &kur, long int CleanupStatus )
{
    double logSI, PhaseAmount = 0., AmThExp, AmountThreshold = 0.;
    double YFcleaned, Yj, YjDiff, YjCleaned=0., MjuPrimal, MjuDual, MjuDiff;
    double CutoffDistortionMBR = 0.1 * pm.DHBM;
    bool KinConstrDC, KinConstrPh;
    bool MassBalanceViolation = false;
    bool NeedToImproveMassBalance = false;
    long int L1k, L1kZeroDCs, k, j, jb = 0, status,
        DCinserted = 0, DCremoved = 0, PHinserted = 0, PHremoved = 0;
    double MjuDiffCutoff = 1e-3; // InsValue;
    SPP_SETTING *pa = paTProfil;
    if( pa->p.GAS > 1e-6 )
         MjuDiffCutoff = pa->p.GAS;
    AmThExp = (double)abs( pa->p.PRD );
    if( AmThExp && AmThExp < 4.)
    {
        AmThExp = 4.;
    }
    AmountThreshold = pow(10.,-AmThExp);

    kfr = -1; kur = -1;
    if( !pm.K2 )
        for( j=0; j<pm.L; j++ )
        {                               // only after the first GEM run
            pm.GamFs[j]=pm.Gamma[j];    // storing a copy of the activity coefficients
        }
    for( j=0; j<pm.L; j++ )
        pm.XY[j]=pm.Y[j];    // Storing a copy of the new speciation vector

    PrimalChemicalPotentials( pm.F, pm.Y, pm.YF, pm.YFA );
    StabilityIndexes( ); // Calculation of phase stability criteria
    (pm.K2)++;

    for(k=0;k<pm.FI;k++)
    {
       L1k = pm.L1[k]; // Number of components in the phase
       KinConstrPh = false;
  /*
       if( pm.PHC[k] == PH_SORPTION || pm.PHC[k] == PH_POLYEL )
       {
           KinConstrDC = false;
           bool AllSC_Suppr = true;
           // Examining if surface species all have been suppressed to zero on this sorption phase
           switch( pm.DCC[j] )
           {
              case DC_SUR_GROUP:
              case DC_SSC_A0: case DC_SSC_A1: case DC_SSC_A2: case DC_SSC_A3: case DC_SSC_A4:
              case DC_WSC_A0: case DC_WSC_A1: case DC_WSC_A2: case DC_WSC_A3: case DC_WSC_A4:
              case DC_SUR_COMPLEX: case DC_SUR_IPAIR: case DC_IESC_A: case DC_IEWC_B:
                   ln_ax_dual -= lnFmol;
                   if( pm.DUL[j] > 0. )
                      AllSC_Suppr = false;

                  break;
              case DC_PEL_CARRIER: case DC_SUR_MINAL: case DC_SUR_CARRIER: // sorbent
                   // Detecting the DC having the non-trivial kinetic constraint
                   if( pm.DUL[j] < 1e6 ) // && Yj >= ( pm.DUL[j] - pa->p.DKIN ) )
                      KinConstrDC = true;
                   if( pm.DLL[j] > 0.0 ) // && Yj <= ( pm.DLL[j] + pa->p.DKIN ) )
                      KinConstrDC = true;
                   break;
              default:
                   break; // error in DC class code
           }
           || KinConstrPh == true )
           goto NextPhase;  // Temporary workaround
       }
*/
 //      else
 //      { // Non - sorption phases
          for(j=jb; j<jb+L1k; j++)
          {  // Checking if a DC in phase is under kinetic control
//          Yj = pm.Y[j];
            KinConstrDC = false;
            // Detecting the DC having the non-trivial kinetic constraint
            if( pm.DUL[j] < 1e6 ) // && Yj >= ( pm.DUL[j] - pa->p.DKIN ) )
               KinConstrDC = true;
            if( pm.DLL[j] > 0.0 ) // && Yj <= ( pm.DLL[j] + pa->p.DKIN ) )
               KinConstrDC = true;
            if( KinConstrDC == true ) // Bug fixed 21.07.2010 DK
               KinConstrPh = true;
          } //  j
 //      }
       if( pm.PHC[k] == PH_SORPTION || pm.PHC[k] == PH_POLYEL
               || KinConstrPh == true )
           goto NextPhase;  // Temporary workaround

       PhaseAmount = pm.XF[k];
       logSI = pm.Falp[k];
       if( logSI > -pa->p.DFM && logSI < pa->p.DF && PhaseAmount < pm.DSM )
       {  // Phase is stable and present in zero or less than DS amount - zeroing off
          bool RemFlagDC = false;
          for(j=jb; j<jb+L1k; j++)
          {
             if( pm.Y[j] )
             {
                DCremoved++; RemFlagDC = true;
             }
             pm.Y[j] = 0.;
          }
          if( RemFlagDC == true )
             PHremoved++;
          goto NextPhase;
       }
       if( logSI >= pa->p.DF )  // 2 - INSERTION CASE
       {  // this phase is stable or over-stable
           if( PhaseAmount < pm.DSM ) // pm.DFYsM )
           {  // phase appears to be lost - insertion of all components of the phase
               if( L1k > 1 )
                  DC_RaiseZeroedOff( jb, jb+L1k, k );
               else
                  pm.Y[jb] = pm.DFYsM; // Spec. value for pure phase insertion
               DCinserted += L1k;
               PHinserted++;
               kfr = k;
               MassBalanceViolation = true;
           } // otherwise (if present), the phase is cleaned up
           goto NextPhase;
       }
       if( logSI <= -pa->p.DFM )  // 3 - ELIMINATION CASE
       {
//         bool Incomplete = false;
          if( PhaseAmount >= pm.DcMinM )
          {  // this phase is present - checking elimination if unstable
             kur = k;
             if( PhaseAmount <= AmountThreshold )
             {   // can be zeroed off
                for(j=jb; j<jb+L1k; j++)
                {
                    if( pm.Y[j] >= pm.DcMinM )
                        DCremoved++;
                    pm.Y[j] = 0.;
                }
             }
             else { // Phase amount too high - elimination may break the mass balance
                MassBalanceViolation = true;
                for(j=jb; j<jb+L1k; j++)
                {
                    if( pm.Y[j] >= pm.DcMinM )
                        DCremoved++;
                    pm.Y[j] = 0.;
                }
             }
             PHremoved++;
          }
          goto NextPhase;
       }
     NextPhase: jb+=pm.L1[k];
   } // k
   // First loop over phases finished

   // Speciation Cleanup mode (PRD: amount threshold exponent)
   if( CleanupStatus && pa->p.PRD ) // && MassBalanceViolation == false ) temporarily disabled DK 13.04.2012
   {  //  not done if insertion or elimination of some phases requires another IPM loop
     jb = 0;
     for(k=0;k<pm.FI;k++)  // Speciation Cleanup loop on phases
     {
       L1k = pm.L1[k]; // Number of components in the phase
       KinConstrPh = false;
       YFcleaned = 0.0;
       for(j=jb; j<jb+L1k; j++)
       {  // Checking if a DC in phase is under kinetic control
          Yj = pm.Y[j];
          KinConstrDC = false;
          // Detecting the DC having the non-trivial kinetic constraint
          if( pm.DUL[j] < 1e6 ) // && Yj >= ( pm.DUL[j] - pa->p.DKIN ) )
          {
              if( Yj >= pm.DUL[j] )
                  // Fixing a very small component constrained from above
                  pm.Y[j] = pm.DUL[j];
              KinConstrDC = true;
          }
          if( pm.DLL[j] >= 0 ) // && Yj <= ( pm.DLL[j] + pa->p.DKIN ) )
          {
              if( Yj < pm.DLL[j] )
              // Fixing a small component constrained from below
                 pm.Y[j] = pm.DLL[j];
              KinConstrDC = true;
          }
          if( KinConstrDC == true )
          {
             YjDiff = fabs( pm.Y[j] -Yj );
             if( YjDiff > CutoffDistortionMBR )
                  NeedToImproveMassBalance = true;
             if( YjDiff > AmountThreshold )
                  MassBalanceViolation = true;
             KinConstrPh = true;
          }
          YFcleaned += pm.Y[j];  // Added by DK 27.08.2012
       } //  j
       if( YFcleaned == 0. )     // No cleanup in the eliminated phase!
          goto NextPhaseC;       // Added by DK 27.08.2012
       PhaseAmount = pm.XF[k];                           // bugfix 04.04.2011 DK
       logSI = pm.Falp[k];  // phase stability criterion
       if( !KinConstrPh && ( logSI >= pa->p.DF || logSI <= -pa->p.DFM ) // DK 13.04.2012
           && !(( pm.PHC[k] == PH_SORPTION || pm.PHC[k] == PH_POLYEL) && PhaseAmount < pm.DSM ))
           goto NextPhaseC;  // Already done in previous part
//       PhaseAmount = pm.XF[k];                         // bugfix 04.04.2011 DK
       if( /* logSI > -pa->p.DFM && */ PhaseAmount >= pm.DSM )
       { // Cleaning up a phase which is present in mass balance
          bool Degenerated = false;
          L1kZeroDCs = 0;
          for(j=jb; j<jb+L1k; j++)
          {
            Yj = YjCleaned = pm.Y[j];
            if( ( pm.DUL[j] < 1e6 && Yj == pm.DUL[j] )
                 || ( pm.DLL[j] > 0 && Yj == pm.DLL[j] ) )
                 continue; // metastability-controlled DCs were already fixed
            if( Yj >= pm.DcMinM )
            {  // we check here in the Ls set only, except metastability constraints
               if( ( pm.DUL[j] < 1e6 && Yj > pm.DUL[j] - pa->p.DKIN )
                  || ( pm.DLL[j] > 0 && Yj < pm.DLL[j] + pa->p.DKIN ) )
                  continue; // we don't clean xDC sitting close to metast. constraints
               MjuPrimal = pm.F[j];   // normalized
               MjuDual = pm.Fx[j];    // /pmp->RT; Fixed DK 12.03.2012
               MjuDiff = MjuPrimal - MjuDual;
               if( fabs( MjuDiff ) > MjuDiffCutoff )
               {
                  YjCleaned = Yj / exp( MjuDiff ); // also applies to a DC in a solution phase
                  if( L1k == 1 )
                  {  // Pure phase
                      if( logSI <= -0.4343*MjuDiffCutoff && YjCleaned < AmountThreshold )
                          YjCleaned = 0.;
                      if( logSI >= pa->p.DF && YjCleaned < pm.DFYsM )
                      {   // over-stable phase in too small amount - insertion and next IPM loop (experimental)
                          YjCleaned = pm.DFYsM;
                          kfr = k;
                          MassBalanceViolation = true;
                      }
                  }
               }
             }
             else if( !(pm.PHC[k] == PH_SORPTION || pm.PHC[k] == PH_POLYEL)  )
             {  // component has been zeroed off - primal chemical potential not defined
                  if( L1k == 1 )
                     continue; // L1kZeroDCs++; No way to improve for pure substance
                  else  // possibly lost DC in a solution phase (only estimated mole fraction available)
                     YjCleaned = pm.EMU[j]*pm.YF[k]; // does not work for adsorption yet
             }
             YjDiff = YjCleaned - Yj;
             if( fabs( YjDiff ) > CutoffDistortionMBR )
             {
                    NeedToImproveMassBalance = true;
                    if( fabs( YjDiff ) > AmountThreshold )
                    {   // Correction was too large - next IPM loop required
                       kur = k;
                       MassBalanceViolation = true;
                      // Provisional: only correction up to the amount threshold
                       if( YjDiff > 0. )
                          pm.Y[j] += AmountThreshold;
                       else
                          pm.Y[j] -= AmountThreshold;
                    }
                    else {  // Reasonable correction - no balance violation expected
                      pm.Y[j] = YjCleaned;
                    }
              }
              else {  // Correction that does not affect the mass balance at all
                  pm.Y[j] = YjCleaned;
              }
              if( pm.Y[j] < pm.DcMinM )
              {  // Corrected amount is too small - DC amount is zeroed off
                  pm.Y[j] = 0.;
                  DCremoved++;
                  L1kZeroDCs++;
              }

        }  // for j
        if( L1k - L1kZeroDCs <= 1 && L1k > 1 )
        {
           if( L1k - L1kZeroDCs )
              Degenerated = true;
           else
              PHremoved++;
           NeedToImproveMassBalance = true;
        }
        if( L1k - L1kZeroDCs == 0 && L1k == 1 )
        {
           PHremoved++;
           NeedToImproveMassBalance = true;
        }
//        goto NextPhaseC;
     }
     NextPhaseC: jb+=pm.L1[k];
   } // k
}

#ifndef IPMGEMPLUGIN
    STEP_POINT("PSSC()");
#ifndef Use_mt_mode
        pVisor->Update(false);  // "PhaseSelectionSpeciationCleanup()"
#endif
#endif
    // Analysis of phase selection and cleanup status
// PZ    // Indicator of PhaseSelection() status (since r1594):
//            0 untouched, 1 phase(s) inserted, 2 insertion done after 5 major IPM loops
// W1     // Indicator of SpeciationCleanup() status (since r1594) 0 untouched,
//           -1 phase(s) removed, 1 some DCs inserted
// K2     // Number of IPM loops performed ( >1 up to 3 because of PhaseSelection() )
    status = 1L;
    CleanupStatus = 0;
    if( !PHinserted )
    {  // No phases were inserted back to mass balance - only cleanup
       status = 1L;
       if( DCinserted )
       {  // some DC in multicomponent phases were restored
          CleanupStatus = 1L;
       }
       else if( NeedToImproveMassBalance )
       {
          CleanupStatus = -1L;
       }
       if( MassBalanceViolation )
       {
           if( pm.K2 < 6 )
              status = 0L;  // attempt will be done to improve by doing one more IPM loop
           else
              status = -1L;  // five loops done, violation persistent - bail out
       }
       if( PHremoved || DCremoved )
           CleanupStatus = -1L;
    }
    else { // phases were inserted - go to another IPM loop
       if( pm.K2 < 6 )
          status = 0L;
       else
          status = -1L;
    }
    if( status == -1L )
    {  // changes in Y vector are not accepted - restore (too many IPM loops)
       for(j=0;j<pm.L;j++)
          pm.Y[j]=pm.XY[j];
    }
// cout << "CleanupStatus= " << CleanupStatus << endl;
    return status;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of (logarithmic) stability indexes logSI for all phases
//
void TMulti::StabilityIndexes( void )
{
    long int L1k, k, j, jb = 0;
    double ln_ax_dual, gamma_primal, x_estimate, StabIndex, logSI;
    double lnFmol = log( H2O_mol_to_kg );  // may not work with mixed-solvent electrolyte
    double lnPc = 0., Xw = 1., lnXw = 0., lnFugPur=0., YFk;
    bool fRestore; char sModPT = SM_UNDEF;

    if( pm.Pc > 1e-29 )
       lnPc = log( pm.Pc );
    if( pm.PHC[0] == PH_AQUEL && pm.YFA[0] >= pm.XwMinM  ) // number of moles of solvent
    {
        Xw = pm.YFA[0] / pm.YF[0];
        lnXw = log( Xw );
    }
    jb=0;
    for(k=0;k<pm.FI;k++)
    {
       L1k = pm.L1[k]; // Number of components in the phase
       YFk = pm.YF[k];
       StabIndex = 0.;
       fRestore = false;
       for(j=jb; j<jb+L1k; j++)
       {  // calculation for all components in all phases
          gamma_primal = pm.Gamma[j];  // primal (external) activity coefficient
if( YFk <= pm.DSM )
{
    if( !pm.K2 && gamma_primal != 1.0 ) // can insert because gamma is available
       fRestore = true;
    if( pm.K2 && pm.GamFs[j] != 1.0 && pm.Gamma[j] == 1.0 )
    {
       gamma_primal = pm.GamFs[j];  // taking saved gamma if the phase was removed
       fRestore = true;
    }
    if( L1k > 1 && pm.sMod )
        sModPT =  pm.sMod[k][SPHAS_TYP];
    if( sModPT == SM_IDEAL )
       fRestore = true; // can always insert a simple ideal solution phase
}
else fRestore = true;
          if( gamma_primal < 1e-33 || gamma_primal > 1e33 )
              gamma_primal = 1.;

          ln_ax_dual = lg_to_ln * pm.Y_la[j];  // DualTh activity
          if( ln_ax_dual < -777. )
              ln_ax_dual = -777.;
          lnFugPur = pm.fDQF[j];  // Pure gas fugacity or end-member DQF parameter

          switch( pm.DCC[j] ) // choice of corrections for estimated mole fractions
          {
             case DC_AQ_ELECTRON: case DC_AQ_PROTON:  case DC_AQ_SPECIES: case DC_AQ_SURCOMP:
                  ln_ax_dual -= lnFmol - lnXw;
                  break;
             case DC_AQ_SOLVENT: case DC_AQ_SOLVCOM:
                  break;
             case DC_GAS_COMP: case DC_GAS_H2O:  case DC_GAS_CO2: case DC_GAS_H2: case DC_GAS_N2:
                  ln_ax_dual -= lnPc;
                  ln_ax_dual -= lnFugPur;
                  break;
             case DC_SCP_CONDEN: case DC_SOL_IDEAL: case DC_SOL_MINOR: case DC_SOL_MAJOR:
             case DC_SOL_MINDEP: case DC_SOL_MAJDEP:
                  ln_ax_dual -= lnFugPur;
                  break;
             case DC_SUR_GROUP:
                  ln_ax_dual -= lnFmol;   // maybe more correction is needed for surface species
                  break;
             case DC_SSC_A0: case DC_SSC_A1: case DC_SSC_A2: case DC_SSC_A3: case DC_SSC_A4:
             case DC_WSC_A0: case DC_WSC_A1: case DC_WSC_A2: case DC_WSC_A3: case DC_WSC_A4:
             case DC_SUR_COMPLEX: case DC_SUR_IPAIR: case DC_IESC_A: case DC_IEWC_B:
                  ln_ax_dual -= lnFmol;
  //                gamma_primal = exp( pm.F0[j] );
                 break;
             case DC_PEL_CARRIER: case DC_SUR_MINAL: case DC_SUR_CARRIER: // sorbent
                  ln_ax_dual -= lnFugPur;
  //                gamma_primal = exp( pm.F0[j] );
                  break;
             default:
                  break; // error in DC class code
          }
          x_estimate = exp( ln_ax_dual )/ gamma_primal;   // estimate of DC concentration
          StabIndex += x_estimate;  // Increment to stability index
          pm.NMU[j] = log( x_estimate );  // may be used for something more constructive
          pm.EMU[j] = x_estimate;         // stored the estimated mole fraction of phase component
       }  // for j
       logSI = log10( StabIndex );
       if( fabs( logSI ) < log10( pm.DSM ) )
           logSI = 0.;
       pm.Falp[k] = logSI; // NormDoubleRound( logSI, 3 );
       if( L1k > 1 && fRestore == false && YFk < pm.DSM )
           pm.Falp[k] = -1.; // provisional - to indicate impossibility to restore
       jb += pm.L1[k];
    }  // for k
}

//=================================================================== old ===========================
/// Checking Karpov phase stability criteria Fa for phases and DCs.
///  Using Selekt2() algorithm by Karpov & Chudnenko (1989)
///  modified by DK in 1995 and in 2007, 2012
///  RaiseStatus: if 1 then zeroed-off DCs are raised to constant in solution phases present
///  in eq state, otherwise (0) this is skipped assuming that SpeciationCleanup will be done.
///  \return 0, if some phases were inserted and a new IPM loop is needed
///             (up to 3 loops possible);
///           1, if the IPM solution is final and consistent, no phases were inserted
///          -1, if the IPM solution is inconsistent after 3 Selekt2() loops
///  In this case, the index of most problematic phase is passed through kfr or
///  kur parameter (parameter value -1 means that no problematic phases were found)
//
long int TMulti::PhaseSelect( long int &kfr, long int &kur, long int RaiseStatus ) //  rLoop )
{
    long int k, j, jb, kf, ku;
    double F1, F2, *F0; // , sfactor;
    SPP_SETTING *pa = paTProfil;
    int rLoop = -1;
//    sfactor = calcSfactor();
    if( !pm.K2 )
        for( j=0; j<pm.L; j++ )
        {                               // only after the first GEM run
            pm.GamFs[j]=pm.lnGam[j];    // Storing a copy of the ln activity coefficients term
        }

    KarpovsPhaseStabilityCriteria( );  // calculation of Karpov phase stability criteria (in pm.Falp)
    F0 = pm.Falp;

    (pm.K2)++;
    kf = -1; ku = -1;  // Index for phase diagnostics
    F1 = pa->p.DF;  // Fixed 29.10.2007  DK
    F2 = -pa->p.DFM;  // Meaning of DFM changed 02.11.2007

    for(k=0;k<pm.FI;k++)
    {
        if( F0[k] > F1 && pm.YF[k] < pm.DSM ) //  pa->p.DS )
        {            // stable phase not in mass balance - to be inserted
            F1=F0[k];
            kf=k;
        }
        if( F0[k] < F2 && pm.YF[k] >= pm.DSM ) // pa->p.DS )  // Fixed 2.11.2007
        {            // unstable phase in mass balance - to be excluded
            F2=F0[k];
            ku=k;
        }
    }
kfr = kf;
kur = ku;
    if( kfr < 0 && kur < 0 )
    {    // No phases to insert/exclude or no Fa distortions found
          // Successful end of iterations of SELEKT2()
        return 1L;
    }

    if( (F2 < -pa->p.DFM ) && ( ku >= 0 ) )
    {
        if( pm.K2 > 4 ) // Three Selekt2() loops have already been done!
                return -1L;   // Persistent presence of unstable phase(s) - bad system!

        // Excluding problematic phases
        do
        {  // excluding all phases with  F2 < DF*sfactor
            for( jb=0, k=0; k < ku; k++ )
                 jb += pm.L1[k];

            DC_ZeroOff( jb, jb+pm.L1[ku], ku ); // Zeroing the phase off
            pm.FI1--;
            // find a new phase to exclude, if any exists
            F2= -pa->p.DFM;
            ku = -1;
            for( k=0; k<pm.FI; k++ )
                if( F0[k] < F2 && pm.YF[k] >= pm.DSM )
                {
                    F2=F0[k];
                    ku=k;
                }
        }
        while( ( F2 <= -pa->p.DFM ) && ( ku >= 0 ) );
    } //if ku

    // Inserting problematic phases
    if( F1 > pa->p.DF && kf >= 0 )
    {
        if( pm.K2 > 4 )
           return -1L;   // Persistent absence of stable phase(s) - bad system!

        // There is a phase for which DF*sfactor threshold is exceeded
        do
        {   // insert this phase and set Y[j] for its components
            // with account for asymmetry and non-ideality
            for( jb=0, k=0; k < kf; k++ )
                 jb += pm.L1[k];

            DC_RaiseZeroedOff( jb, jb+pm.L1[kf], kf );

            pm.FI1++;  // check phase rule

            if( pm.FI1 >= pm.NR+1 )
               break;   // No more phases can be inserted

            // find a new phase to insert, if any exists
            F1= pm.lowPosNum; // was about 1e-16
            kf = -1;
            for( k=0; k<pm.FI; k++ )
                if( F0[k] > F1 && pm.YF[k] < pm.DSM )
                {
                    F1=F0[k];
                    kf=k;
                }
        }
        while( F1 > pa->p.DF && kf >= 0 );
        // end of insertion cycle
    } // if kf changed SD 03/02/2009

    if( RaiseStatus )
    {    // Raise zeros in DC amounts in phases-solutions in the case if some phases
         // were inserted or excluded      - this option is probably not needed at all!
         //        double RaiseZeroVal = pm.DHBM;  // Added 29.10.07  by DK
        jb=0;
        for(k=0;k<pm.FIs;k++)
        {
            if( ( pm.YF[k] >= pm.DSM ) || ( pm.pNP && rLoop < 0 ) ) // Only in phase present in mass balance!
            {                            // (acc. to definition of L_S set) PIA only if initial!
                 pm.YF[k]=0.;
                 for(j=jb;j<jb+pm.L1[k];j++)
                 {
                    if( pm.Y[j] < min( pm.lowPosNum, pm.DcMinM ) )  // fixed 30.08.2009
                        pm.Y[j] = RaiseDC_Value( j ); // bugfix 29.10.07
                    pm.YF[k] += pm.Y[j]; // calculate new amounts of phases
                 }
            }
            jb+=pm.L1[k];
        }
    }
#ifndef IPMGEMPLUGIN
    STEP_POINT("Selekt2 procedure");
#ifndef Use_mt_mode
        pVisor->Update(false);  // "PhaseSelection"
#endif
#endif
    if( pm.K2 > 1 )
    { // more then the first step - but the IPM solution has not improved
//      double RaiseZeroVal = pm.DHBM*0.1;   // experimental
       for(j=0;j<pm.L;j++)
          if( fabs(pm.Y[j]- pm.XY[j]) > RaiseDC_Value( j ) ) //
               goto S6;
       // pm.PZ=2; // No significant change has been done by Selekt2()
       return 1L;
    }
S6: // copy of X vector has been changed by Selekt2() algorithm - store
    for(j=0;j<pm.L;j++)
        pm.XY[j]=pm.Y[j];

    return 0L;  // Another loop is needed
}

/// New function to improve on raising zero values in PhaseSelect() and after SolveSimplex()
double TMulti::RaiseDC_Value( const long int j )
{
        double RaiseZeroVal = pm.DFYsM;

        switch(pm.DCC[j] )
        {
    case DC_AQ_PROTON:
    case DC_AQ_ELECTRON:
    case DC_AQ_SPECIES:
    case DC_AQ_SURCOMP:	RaiseZeroVal = pm.DFYaqM;
                                                break;
    case DC_SOL_IDEAL:
    case DC_GAS_COMP:
    case DC_GAS_H2O:
    case DC_GAS_CO2:
    case DC_GAS_H2:
    case DC_GAS_N2:	RaiseZeroVal = pm.DFYidM;
                                                break;
    case DC_AQ_SOLVENT:
    case DC_AQ_SOLVCOM: RaiseZeroVal = pm.DFYwM;
                                                break;
    case DC_SOL_MINOR: case DC_SOL_MINDEP: RaiseZeroVal = pm.DFYhM;
                                                break;
    case DC_SOL_MAJOR: case DC_SOL_MAJDEP: RaiseZeroVal = pm.DFYrM;
                                                break;
        // adsorption
    case DC_SSC_A0:    case DC_SSC_A1:    case DC_SSC_A2:    case DC_SSC_A3:    case DC_SSC_A4: // obsolete
    case DC_WSC_A0:    case DC_WSC_A1:    case DC_WSC_A2:    case DC_WSC_A3:    case DC_WSC_A4: // obsolete
    case DC_SUR_GROUP:
    case DC_SUR_COMPLEX:
    case DC_SUR_IPAIR:
    case DC_IESC_A:
    case DC_IEWC_B:		RaiseZeroVal = pm.DFYaqM;
                                                break;
    case DC_PEL_CARRIER:
    case DC_SUR_MINAL:
    case DC_SUR_CARRIER: RaiseZeroVal = pm.DFYrM;
                                                 break;
    case DC_SCP_CONDEN:	 RaiseZeroVal = pm.DFYcM;
                                                 break;
    default:
                 break;
        }
        return RaiseZeroVal;
}

//--------------------- End of ipm_chemical.cpp ---------------------------
