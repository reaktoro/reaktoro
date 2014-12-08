//-------------------------------------------------------------------
// $Id: ipm_chemical2.cpp 845 2013-06-20 15:58:57Z kulik $
//
/// \file ipm_chemical2.cpp
/// Implementation of chemistry-specific functions (concentrations,
/// activity coefficients, adsorption models etc.)
/// for the IPM convex programming Gibbs energy minimization algorithm
//
// Copyright (c) 1992-2012  D.Kulik, S.Dmitrieva, K.Chudnenko
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
#include "m_param.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculating bulk stoichiometry of a multicomponent phase
//
void TMulti::phase_bcs( long int N, long int M, long int jb, double *A, double X[], double BF[] )
{
    long int ii, i, j;
    double Xx;

    if( !A || !X || !BF )
        return;
    for( i=0; i<N; i++ )
          BF[i] = 0.;

    for( j=0; j<M; j++ )
    {
        Xx = X[j];
        if( fabs( Xx ) < 1e-16 )  // was 1e-12
            continue;
        for( ii=arrL[j+jb]; ii<arrL[j+jb+1]; ii++ )
        {  i = arrAN[ii];
            BF[i] += A[i+j*N] * Xx;
        }
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Adds phase to total bulk stoichiometry of all solid phases in the system
// Done on request by TW in November 2006
//
void TMulti::phase_bfc( long int k, long int jj )
{
    long int ii, i, j;
    double Xx;

    if( pm.PHC[k] == PH_AQUEL || pm.PHC[k] == PH_GASMIX ||
        pm.PHC[k] == PH_FLUID || pm.PHC[k] == PH_PLASMA ||
        pm.PHC[k] == PH_SIMELT || pm.PHC[k] == PH_LIQUID )
        return;
    for( j=0; j<pm.L1[k]; j++ )
    {
        Xx = pm.X[j+jj];
        if( fabs( Xx ) < 1e-12 )
            continue;
        for( ii=arrL[j+jj]; ii<arrL[j+jj+1]; ii++ )
        {  i = arrAN[ii];
           pm.BFC[i] += pm.A[i+(jj+j)*pm.N] * Xx;
        }
    }
}

/// Returns mass of all solid phases in grams (from the BFC vector)
double TMulti::bfc_mass( void )
{
   double TotalMass = 0.;
   for(long int i = 0; i<pm.N; i++ )
     TotalMass += pm.BFC[i]*pm.Awt[i];
   return TotalMass;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#define  a(j,i) ((*(pm.A+(i)+(j)*pm.N)))
/// Calculation of dual chemical potentials, activities, and primal
/// concentrations for DCs (indexed jb to je) in a k-th phase.
// Input arrays X, XF, XFA,  input factors: Factor, MMC
//
//  Do we need this all in GEMIPM ?
//
void TMulti::CalculateConcentrationsInPhase( double X[], double XF[], double XFA[],
              double Factor, double MMC, double /*Dsur*/, long int jb, long int je, long int k)
{
    long int j, ii, i;
    double Muj, /* DsurT=0.0,*/ SPmol, lnFmol=4.016535;
//    SPP_SETTING *pa = &TProfil::pm->pa;

    if( pm.PHC[0] == PH_AQUEL )
    {  // mole fraction to molality conversion
        if( !k ) lnFmol = log(1000./MMC);  // aq species
        else lnFmol = log( H2O_mol_to_kg ); // 4.016535; 	   // other species
    }

    for( j=jb; j<je; j++ )
    { // loop over DC - with important bugfixes from 02.04.2003
        Muj = DC_DualChemicalPotential( pm.U, pm.A+j*pm.N, pm.NR, j );
        pm.Fx[j] = Muj; // *pmp->RT; Fixed DK 12.03.2012

//        if( X[j] <= pm.lowPosNum )
        if( X[j] <= pm.DcMinM )
        { // zeroing off
            pm.Wx[j] = 0.0;
            pm.VL[j] = log( pm.lowPosNum );
            pm.Y_w[j] = 0.0;
            pm.lnGam[j] = 0.0;
            if( pm.PHC[0] == PH_AQUEL )
               pm.Y_m[j] = 0.0;
            switch( pm.DCC[j] ) // choice of expressions
            {                      // since 10.03.2008, changed the concept of DualTh activity
               case DC_SCP_CONDEN:
                    pm.Y_la[j] = ln_to_lg * ( Muj - pm.G0[j] );
                    break;
               case DC_AQ_ELECTRON: case DC_AQ_PROTON:  case DC_AQ_SPECIES: case DC_AQ_SURCOMP:
                    pm.Y_la[j] = ln_to_lg*(Muj - pm.G0[j] + lnFmol );
                    break;
               case DC_AQ_SOLVENT: case DC_AQ_SOLVCOM:
                    pm.Y_la[j] = ln_to_lg* (Muj - pm.G0[j] );
                    break;
               case DC_GAS_COMP: case DC_GAS_H2O:  case DC_GAS_CO2:   // gases
               case DC_GAS_H2: case DC_GAS_N2:
                    pm.Y_la[j] = ln_to_lg * ( Muj - pm.G0[j] );
                    if( pm.Pc > 1e-29 )
                        pm.Y_la[j] += log10( pm.Pc );
                    break;
               case DC_SOL_IDEAL: case DC_SOL_MINOR: case DC_SOL_MAJOR: case DC_SOL_MINDEP: case DC_SOL_MAJDEP:
                    pm.Y_la[j] = ln_to_lg * ( Muj - pm.G0[j] );
                    break;
               case DC_SUR_GROUP:
                    pm.Y_la[j] = ln_to_lg * ( Muj - pm.G0[j] ); // + lnFmol ); corr. 06.10.10 DK
                    break;
               case DC_SSC_A0: case DC_SSC_A1: case DC_SSC_A2: case DC_SSC_A3:
               case DC_SSC_A4: case DC_WSC_A0: case DC_WSC_A1: case DC_WSC_A2:
               case DC_WSC_A3: case DC_WSC_A4: case DC_SUR_COMPLEX:
               case DC_SUR_IPAIR: case DC_IESC_A: case DC_IEWC_B:
                    pm.Y_la[j] = ln_to_lg * ( Muj - pm.G0[j] ); // + lnFmol ); corr. 06.10.10 DK
                    break;
               case DC_PEL_CARRIER: case DC_SUR_MINAL: case DC_SUR_CARRIER: // sorbent
                    pm.Y_la[j] = ln_to_lg * ( Muj - pm.G0[j] );
                    break;
               default:
                    break; // error in DC class code
            }
            continue;
        }
        // calculation of the mole fraction
        pm.Wx[j] = X[j]/XF[k];
        if( X[j] > min( pm.lowPosNum, pm.DcMinM ) )
            pm.VL[j] = log( pm.Wx[j] );     // this is used nowhere except in some scripts. Remove?
        else pm.VL[j] = log( pm.lowPosNum );   // debugging 29.11.05 KD
        pm.Y_la[j] = 0.0;
        switch( pm.DCC[j] ) // choice of expressions
        {
        case DC_SCP_CONDEN:
            pm.Wx[j] = 1;
            pm.VL[j] = 0.0;
            if( pm.LO )
            {   //  bugfix DK 08.02.10
                pm.Y_m[j] = X[j]*Factor; // molality
            }
            pm.Y_w[j] = // mass % of the system
                          1e2 * X[j] * pm.MM[j] / pm.MBX;
            pm.Y_la[j] = ln_to_lg * ( Muj - pm.G0[j] );
            pm.FVOL[k] += pm.Vol[j]*X[j];
            break;
        case DC_AQ_ELECTRON:
            pm.Y_m[j] = 0.0;
            pm.Y_la[j] = 0.0 - pm.pe;
            pm.Y_w[j] = 0.0;
            break;
        case DC_AQ_PROTON:  // in molal scale!
            pm.pH = -ln_to_lg*(Muj-pm.G0[j] + lnFmol );
        case DC_AQ_SPECIES:
            pm.IC += 0.5* X[j]*Factor *(pm.EZ[j]*pm.EZ[j]); // increment to effective molal ionic strength
        case DC_AQ_SURCOMP:
            SPmol = X[j]*Factor;  // molality
 //           pm.IC += 0.5* SPmol *(pm.EZ[j]*pm.EZ[j]); // Bugfix DK 21.10.2011 'K' species don't count here!
            pm.FVOL[k] += pm.Vol[j]*X[j]; // fixed 04.02.03 KD
            pm.Y_m[j] = SPmol;
            pm.Y_la[j] = ln_to_lg*(Muj - pm.G0[j] + lnFmol );
            pm.Y_w[j] = 1e6 * X[j] * pm.MM[j] / pm.FWGT[k];
            if(  pm.DCC[j] != DC_AQ_SURCOMP )
            {  // Bugfix 21.10.2011  DK - excluding 'K' species from total IC molality
               //  Optimized for performance - calculation inline
               for( i=arrL[j]; i<arrL[j+1]; i++ )
               {  ii = arrAN[i];
                  if( ii>= pm.NR )
                    continue;
                  pm.IC_m[ii] += SPmol* a(j,ii);  // total aqueous molality
                  pm.IC_wm[ii] += X[j]* a(j,ii);  // total aqueous mass concentration
               } 
            }
            break;
        case DC_AQ_SOLVENT: // mole fractions of solvent
        case DC_AQ_SOLVCOM:
//            pm.Y_m[j] = X[j]/XFA[k];  Replaced by DK 12.03.2012
            pm.Y_m[j] = H2O_mol_to_kg;
            pm.Y_w[j] = 1e3*X[j]*pm.MM[j]/pm.FWGT[k];
            pm.FVOL[k] += pm.Vol[j]*X[j];
            pm.Y_la[j] = ln_to_lg* (Muj - pm.G0[j] );
            break;
        case DC_GAS_COMP:
        case DC_GAS_H2O:
        case DC_GAS_CO2:   // gases
        case DC_GAS_H2:
        case DC_GAS_N2:
            pm.FVOL[k] += pm.Vol[j]*X[j];
            pm.Y_la[j] = ln_to_lg * ( Muj - pm.G0[j] );
            if( pm.Pc > 1e-9 )
                pm.Y_la[j] += log10( pm.Pc );
            break;
        case DC_SOL_IDEAL:
        case DC_SOL_MINOR:   //solution end member
        case DC_SOL_MAJOR:
        case DC_SOL_MINDEP:
        case DC_SOL_MAJDEP:
            pm.FVOL[k] += pm.Vol[j]*X[j];
            pm.Y_la[j] = ln_to_lg * ( Muj - pm.G0[j] );
            if( pm.LO )
            {   //  bugfix DK 16.04.2012
                pm.Y_m[j] = X[j]*Factor; // molality
            }
            break;
        case DC_SUR_GROUP: // adsorption:
            pm.Y_m[j] = X[j]*Factor; // molality
            pm.Y_w[j] =  // mg/g sorbent
                1e3 * X[j] * pm.MM[j] / (MMC*XFA[k]);
            pm.Y_la[j] = ln_to_lg * ( Muj - pm.G0[j] ); // + lnFmol ); corr. 06.10.10
            pm.FVOL[k] += pm.Vol[j]*X[j]; // fixed 11.03.2008 KD
            break;
        case DC_SSC_A0:
        case DC_SSC_A1:
        case DC_SSC_A2:
        case DC_SSC_A3:
        case DC_SSC_A4:
        case DC_WSC_A0:
        case DC_WSC_A1:
        case DC_WSC_A2:
        case DC_WSC_A3:
        case DC_WSC_A4:  // case DC_SUR_GROUP:
        case DC_SUR_COMPLEX:
        case DC_SUR_IPAIR:
        case DC_IESC_A:
        case DC_IEWC_B:
            pm.Y_m[j] = X[j]*Factor; // molality
            pm.Y_w[j] =  // mg/g sorbent
                1e3 * X[j] * pm.MM[j] / (MMC*XFA[k]);
//           DsurT = MMC * pm.Aalp[k] * pa->p.DNS*1.66054e-6;
            pm.Y_la[j] = ln_to_lg * ( Muj - pm.G0[j] ); // + lnFmol - pm.GEX[j] + Dsur + DsurT/( 1.0+DsurT )
            pm.FVOL[k] += pm.Vol[j]*X[j]; // fixed 11.03.2008   KD
            break;
        case DC_PEL_CARRIER:
        case DC_SUR_MINAL:
        case DC_SUR_CARRIER: // sorbent
            pm.Y_m[j] = X[j]*Factor; // molality
            pm.Y_w[j] = 0.0;
            if( pm.YF[0] >= pm.DSM )
              pm.Y_w[j] = // mg of sorbent per kg aq solution
                1e6 * X[j] * pm.MM[j] / pm.FWGT[0];
            pm.Y_la[j] = ln_to_lg * ( Muj - pm.G0[j] ); // - pm.GEX[j] + Dsur - 1. + 1./(1.+Dsur) - DsurT + DsurT/(1+DsurT)
            pm.FVOL[k] += pm.Vol[j]*X[j];
            break;
        default:
            break; // error in DC class code
        }
        ;
    }   // j
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of derived values (concentrations etc.) on IPM iteration
///  from X,XF, and XFA vectors. Also calculates pH, pe, Eh
// This function has to be rewritten using new set of built-in
// chemical functions.
//
void TMulti::CalculateConcentrations( double X[], double XF[], double XFA[])
{
    long int k, ii, i, j, ist, jj, jja;
    double Factor=0.0, Dsur=0.0, MMC=0.0;
    SPP_SETTING *pa = paTProfil;

//    if( pm.Ls < 2 || !pm.FIs )  Temporary disabled  09.03.2010 DK
//        return;

    for( i=0; i< pm.N; i++ )
     pm.BFC[i] = 0.;         // cleanup for phase_bfc() calculation

    for( j=0; j<pm.Ls; j++ )
    {
        pm.Wx[j] = 0.;
        pm.VL[j] = 0.;
    }
    j=0;
    pm.VXc = 0.0;
    for( k=0; k<pm.FI; k++ )
    { // cycle by phases
        i=j+pm.L1[k];
        pm.FWGT[k] = 0.0;
        pm.FVOL[k] = 0.0;
        //   Dsur = 0.0;

        if( XF[k] > pm.DSM && !( pm.PHC[k] == PH_SORPTION && XFA[k] <= pa->p.ScMin ))
           phase_bfc( k, j );

        if( k >= pm.FIs || pm.L1[k] == 1 )
        { // this is a single- component phase
            pm.Wx[j] = 1.0; // SD 04/05/2010
            if( XF[k] < pm.DSM )
            {   // This phase to be zeroed off
                if( pm.LO )
                    pm.Y_m[j] = 0.0;
                pm.Wx[j] = 0.0;
                pm.Y_w[j] = 0.0;
                pm.Fx[j] = DC_DualChemicalPotential( pm.U, pm.A+j*pm.N, pm.NR, j );
                pm.Y_la[j] = ln_to_lg * ( pm.Fx[j] - pm.G0[j] ); // -pm.GEX[j]
//                pm.Fx[j] *= pm.RT;     // el-chem potential
                goto NEXT_PHASE;
            }
            pm.VL[j] = 0.0;
            if( pm.LO && XFA[0] > 0 )
                pm.Y_m[j] = X[j] * 1000./18.01528/XFA[0]; // molality
            pm.Y_w[j] = // mass % in the system
                1e2 * X[j] * pm.MM[j] / pm.MBX;
            pm.Fx[j] = DC_DualChemicalPotential( pm.U, pm.A+j*pm.N, pm.NR, j );
            pm.Y_la[j] = ln_to_lg * ( pm.Fx[j] - pm.G0[j] ); // - pm.GEX[j]
//            pm.Fx[j] *= pm.RT;     // el-chem potential
            pm.FWGT[k] += X[j] * pm.MM[j];
            pm.FVOL[k] += X[j] * pm.Vol[j];
            goto NEXT_PHASE;
        }
        if( pm.LO && !k )
        {    // aqueous phase present
            pm.IC=0.;                   // fix 20.06.13 DK
            for( ii=0; ii<pm.N; ii++ )
            {
                pm.IC_m[ii] = 0.0;
                pm.IC_lm[ii] = 0.0;
                pm.IC_wm[ii] = 0.0;
            }
        }
        if( XF[k] <= pm.DSM ||
            (pm.PHC[k] == PH_AQUEL && ( XFA[k] <= pm.XwMinM || XF[k] <= pm.DSM ) )
                || ( pm.PHC[k] == PH_SORPTION && XFA[k] <= pm.ScMinM ))
        {
            for( jj=0; jj<pm.N; jj++)
             pm.BF[k*pm.N+jj] = 0.;

            for(jj=j; jj<i; jj++)   // Loop added 10.03.01  KD (GTDEMO)
            {
                pm.Wx[j] = 0.0;
                if( pm.LO )
                    pm.Y_m[jj] = 0.0;
                pm.Y_w[jj] = 0.0;
                pm.Fx[jj] = DC_DualChemicalPotential( pm.U, pm.A+jj*pm.N, pm.NR, jj );
                pm.Y_la[jj] = ln_to_lg * ( pm.Fx[jj] - pm.G0[jj] );
                if(pm.PHC[k] == PH_AQUEL ) // || pm.PHC[k] == PH_SORPTION ) corr. 06.10.10 DK
                   pm.Y_la[jj] += 1.74438;
                if(pm.PHC[k] == PH_GASMIX || pm.PHC[k] == PH_PLASMA )
                   pm.Y_la[jj] += log10( pm.Pc );
//                pm.Fx[jj] *= pm.RT;     // el-chem potential
                pm.lnGam[jj] = 0.0;
            }
            goto NEXT_PHASE;
        }
        // calculate bulk stoichiometry of a multicomponent phase
        phase_bcs( pm.N, pm.L1[k], j, pm.A+j*pm.N, pm.X+j, pm.BF+k*pm.N );

        switch( pm.PHC[k] )
        {
        case PH_AQUEL:
            MMC = 0.0; // molar mass of carrier
//            Dsur = XFA[k]/XF[k] - 1.0; // Asymm.corr. - aq only!
//            if( XFA[k] > pm.lowPosNum )
            if( XFA[k] > pm.XwMinM )
            {
                for(jj=j; jj<i; jj++)
                    if( pm.DCC[jj] == DC_AQ_SOLVENT ||
                            pm.DCC[jj] == DC_AQ_SOLVCOM )
                        MMC += pm.MM[jj]*X[jj]/XFA[k];
            }
            else MMC=18.01528; // Assuming water-solvent
//            if( (XFA[k] > pm.lowPosNum) && (MMC > pm.lowPosNum) )
            if( (XFA[k] > pm.XwMinM) && (MMC > pm.lowPosNum) )
                Factor = 1000./MMC/XFA[k]; // molality
            else Factor = 0.0;
//            pm.IC=0.;
            pm.pe = ln_to_lg* pm.U[pm.N-1];
            pm.Eh = 0.000086 * pm.U[pm.N-1] * pm.T;
        case PH_GASMIX:
        case PH_FLUID:
        case PH_PLASMA:
        case PH_SIMELT:
        case PH_HCARBL:
        case PH_SINCOND:
        case PH_SINDIS:
        case PH_LIQUID:
            pm.YFk = XF[k];
            for(jj=j; jj<i; jj++)
            {
                if( X[jj] > pm.DcMinM)      // fixed 30.08.2009 DK
                    pm.FWGT[k] += X[jj]*pm.MM[jj];
            }
            break;
        case PH_POLYEL:
        case PH_SORPTION: // only sorbent end-members!
            pm.YFk = XFA[k];
            MMC=0.0;

            for( ist=0; ist<pm.FIat; ist++ )
                pm.XFTS[k][ist] = 0.0;
//           if( XFA[k] < pm.lowPosNum ) XFA[k] = pm.lowPosNum;
            if( XFA[k] < pm.ScMinM ) XFA[k] = pm.ScMinM;
            for( jj=j; jj<i; jj++ )
            {
               jja = jj - ( pm.Ls - pm.Lads );
                if( pm.DCC[jj] == DC_SUR_CARRIER ||
                        pm.DCC[jj] == DC_SUR_MINAL ||
                        pm.DCC[jj] == DC_PEL_CARRIER )
                {
                    MMC += pm.MM[jj]*X[jj]/XFA[k];
                    // Only sorbent mass
                    pm.FWGT[k] += X[jj]*pm.MM[jj];
                }
                else
                {
/*!!!!!*/           ist = pm.SATX[jja][XL_ST];
                    pm.XFTS[k][ist] += X[jj];
                }
            }
            pm.logYFk = log(pm.YFk);
//            Dsur = XFA[k]/XF[k] - 1.0;  // Also for sorption phases
//            if( Dsur <= -1.0 ) Dsur = -0.999999;
            break;
        default:
             return; // Phase class code error!
        }
        // calculation of species concentrations in k-th phase
        CalculateConcentrationsInPhase( X, XF, XFA, Factor, MMC, Dsur, j, i, k );

NEXT_PHASE:
        pm.VXc += pm.FVOL[k];
        if( pm.PHC[k] == PH_AQUEL && XF[k] > pm.DSM && XFA[k] > pm.XwMinM )
            for( ii=0; ii<pm.NR; ii++ )
            {
               if( pm.LO  )
               { if( pm.IC_m[ii] >= pa->p.DB )
                    pm.IC_lm[ii] = ln_to_lg*log( pm.IC_m[ii] );
                else
                    pm.IC_lm[ii] = 0;
                if( pm.FWGT[k] >= pa->p.DB )
                    pm.IC_wm[ii] *= pm.Awt[ii]*1000./pm.FWGT[k];
                else
                    pm.IC_wm[ii] = 0;
               }
            }
        j = i;
    }  // k
}

//--------------------------------------------------------------------------------
/// Calculation of surface charge densities on multi-surface sorption phase
void TMulti::IS_EtaCalc()
{
    long int k, i, ist, isp, j=0, ja;
    double XetaS=0., XetaW=0.,  Ez, CD0, CDb;
//    SPP_SETTING *pa = &TProfil::pm->pa;

    for( k=0; k<pm.FIs; k++ )
    { // loop over phases
        i=j+pm.L1[k];
        if( pm.FIat > 0 )
            for( ist=0; ist<pm.FIat; ist++ )
            {
                pm.XetaA[k][ist] = 0.0;
                pm.XetaB[k][ist] = 0.0;
                pm.XetaD[k][ist] = 0.0;     // added 12.09.05  KD
            }

        if( pm.XF[k] <= pm.DSM ||
                (pm.PHC[k] == PH_AQUEL && ( pm.X[pm.LO] <= pm.XwMinM //  pa->p.XwMin
                 || pm.XF[k] <= pm.DHBM ) )
             || (pm.PHC[k] == PH_SORPTION && pm.XF[k] <= pm.ScMinM ) ) //  pa->p.ScMin) )
            goto NEXT_PHASE;

        switch( pm.PHC[k] )
        {  // initialization according to the phase class
        case PH_AQUEL:  // aqueous solution
            pm.Yw = pm.XFA[k];
            XetaW = 0.0;
            break;
        case PH_PLASMA:
        case PH_SIMELT:
            XetaS = 0.0;
            break;
        case PH_POLYEL:
        case PH_SORPTION: // reserved
            break;
        default:
            break;
        }
        for( ; j<i; j++ )
        { // loop over DC for calculating total phase charge
            if( pm.X[j] <= pm.lowPosNum*100. )
                continue; // Skipping too low concentrations
            ja = j - ( pm.Ls - pm.Lads );

            switch( pm.DCC[j] ) // select expressions for species classes
            {
            case DC_AQ_ELECTRON:    case DC_AQ_PROTON:    case DC_AQ_SPECIES:  case DC_AQ_SURCOMP:
                XetaW += pm.X[j]*pm.EZ[j];
            case DC_AQ_SOLVCOM:    case DC_AQ_SOLVENT:
                break;
            case DC_PEL_CARRIER:  case DC_SUR_MINAL:
            case DC_SUR_CARRIER: // charge of carrier: ???????
                                 // pm.XetaA[k] += pm.X[j]*pm.EZ[j];
                break;
                // surface species
            case DC_SSC_A0: case DC_SSC_A1: case DC_SSC_A2:  case DC_SSC_A3:
            case DC_SSC_A4: case DC_WSC_A0: case DC_WSC_A1:  case DC_WSC_A2:
            case DC_WSC_A3: case DC_WSC_A4:
            case DC_SUR_GROUP: case DC_SUR_COMPLEX: case DC_SUR_IPAIR:
            case DC_IESC_A:
            case DC_IEWC_B: // Get ist - index of surface type
                            // and  isp - index of surface plane
                ist = pm.SATX[ja][XL_ST]; // / MSPN;
                isp = pm.SATX[ja][XL_SP]; // % MSPN;
                        // isp  index of outer surface charge allocation  (new)
                // Getting charge distribution information
                CD0 = pm.MASDJ[ja][PI_CD0];
                    // species charge that goes into 0 plane
                CDb = pm.MASDJ[ja][PI_CDB];
          // species charge that goes into 1, 2 or 3 plane according to isp value
                Ez = pm.EZ[j];  // take formula charge as default
                if( !isp )
                { // This is the 0 (A) plane only - no charge distribution!
                    if( fabs( CD0 ) > 1e-20 ) // Only if 0-plane charge is given in the table
                       Ez = CD0;
                    pm.XetaA[k][ist] += pm.X[j]*Ez;
                }
                else
                { // The charge distribution (CD) is specified
                    if( pm.SCM[k][ist] == SC_MTL )
                    {   // Modified TL: Robertson, 1997; also XTLM Kulik 2002
//                        if( fabs( CDb ) > 1e-20 )  // Doubtful...
//                           Ez = CDb;
                        pm.XetaB[k][ist] += pm.X[j]*CDb;
                    }
                    else if( pm.SCM[k][ist] == SC_TLM )
                    {
// New CD version of TLM Hayes & Leckie, 1987  added 25.10.2004
                        pm.XetaB[k][ist] += pm.X[j] * CDb;
                        pm.XetaA[k][ist] += pm.X[j] * CD0;
                    }
                    else if( pm.SCM[k][ist] == SC_3LM )
                    {
// CD 3-layer model (Hiemstra e.a. 1996) added 12.09.2005 by KD
                        if( isp == 1 )
                            pm.XetaB[k][ist] += pm.X[j] * CDb;
                        if( isp == 2 )
                            pm.XetaD[k][ist] += pm.X[j] * CDb;
                        pm.XetaA[k][ist] += pm.X[j] * CD0;
                    }
                    else if( pm.SCM[k][ist] == SC_BSM )
                    { // Basic Stern model Christl & Kretzschmar, 1999
// New CD version of BSM  added 25.10.2004
                        pm.XetaB[k][ist] += pm.X[j] * CDb;
                        pm.XetaA[k][ist] += pm.X[j] * CD0;
                    }
                    else if( pm.SCM[k][ist] == SC_MXC )
                    { // BSM for ion exchange on perm.charge surface
                        if( fabs( CDb ) > 1e-20 )  // Doubtful...
                           Ez = CDb;
                        pm.XetaB[k][ist] += pm.X[j]*Ez;
                        pm.XetaA[k][ist] += pm.X[j]*CD0;  // added for testing
                    }
                    else if( pm.SCM[k][ist] == SC_CCM )
                    { // Added 25.07.03 to implement the extended CCM Nilsson ea 1996
// New CD version of BSM  added 25.10.2004
                           pm.XetaB[k][ist] += pm.X[j] * CDb;
                           pm.XetaA[k][ist] += pm.X[j] * CD0;
                    }
                 //    case DC_SUR_DL_ION:  XetaS += pm.X[j]*pm.EZ[j];
                }
                break;
            default:
                XetaS += pm.X[j]*pm.EZ[j];
                break;
            }
        }   // j
        // compare pm.Xetaf[k]+pm.XetaA[k]+pm.XetaB[k] and XetaS
        // Test XetaW
NEXT_PHASE:
        j = i;
        if( pm.LO && !k && pm.FIat > 0 )
        {
            pm.XetaA[k][0] = XetaW;
            pm.XetaB[k][0] = XetaW;
            pm.XetaD[k][0] = XetaW;
        }
        if( (pm.PHC[k] == PH_PLASMA || pm.PHC[k] == PH_SIMELT)
                && pm.FIat)
            pm.XetaA[k][0] = XetaS;
    }  // k
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of the surface potential pm[q].XpsiD[k] on diffuse.
/// layer plane on k-th sorption phase. From total charge pm.Xeta[k]
/// ( in moles ) using Gouy-Chapman equation.
/// Strictly valid at PSI < 30 mV. Modified by DAK 5 Jan 2000
/// to add a Basic Stern EDL model.
///    Added 13.03.2008 by DK: returns int value showing (if not 0)
///    that some extreme values were reached for charge densities or
///    electric potentials (for detecting bad PIA), 0 otherwise
//
long int
TMulti::GouyChapman(  long int, long int, long int k )
{
    long int ist, status=0;
    double SigA=0., SigD=0., SigB=0., SigDDL=0.,
      XetaA[MST], XetaB[MST], XetaD[MST], f1, f3, A, Sig, F2RT, I, Cap;
    if( pm.XF[k] < pm.ScMinM ) // TProfil::pm->pa.p.ScMin )
        return status; // no sorbent

    // sorbent mass in grams
    pm.YFk = pm.FWGT[k];
    if( pm.XF[k] < pm.DSM )
       pm.YFk = pm.lowPosNum;

    for( ist=0; ist<pm.FIat; ist++ )  // loop over surface types
    {
        double PsiD=0.0, PSIo=0.0, PsiA=0.0, PsiB=0.0; // Cleanup 05.12.2009 DK
        double ConvFactor = 1.;

        XetaA[ist] = XetaB[ist] = XetaD[ist] = 0.0;
        if( pm.SCM[k][ist] == SC_NOT_USED || pm.Nfsp[k][ist] < 1e-9  )
            continue;
        ConvFactor = F_CONSTANT / pm.YFk / pm.Aalp[k] / pm.Nfsp[k][ist];
        // Calculation of charge densities (now limited to total charges > max. balance residual)
        if( fabs( pm.XetaA[k][ist]) > pm.DHBM ) // pm.lowPosNum*100. )
            XetaA[ist] = pm.XetaA[k][ist] * ConvFactor; // in C/m2
        if( fabs( pm.XetaB[k][ist]) > pm.DHBM ) // pm.lowPosNum*100. )   // moles
            XetaB[ist] = pm.XetaB[k][ist] * ConvFactor; // C/m2
        if( fabs( pm.XetaD[k][ist]) > pm.DHBM ) // pm.lowPosNum*100. ) // moles
            XetaD[ist] = pm.XetaD[k][ist] * ConvFactor; // C/m2

        // Limit maximum charge densities to prevent divergence
        if( fabs(XetaA[ist]) > 1.4 )
        {
// cout << "EDL charge density A " << XetaA[ist] << " truncated to +- 0.7 C/m2" <<
//        "  IT= " << pm.IT << " k= " << k << " ist= " << ist << endl;
            XetaA[ist] = XetaA[ist] < 0.0 ? -1.4: 1.4;
            status = 60;
        }
        if( fabs(XetaB[ist]) > 2.0 )
        {
// cout << "EDL charge density B " << XetaB[ist] << " truncated to +- 1.7 C/m2" <<
//        "  IT= " << pm.IT << " k= " << k << " ist= " << ist << endl;
            XetaB[ist] = XetaB[ist] < 0.0 ? -2.0: 2.0;
            status = 61;
        }
        if( fabs(XetaD[ist]) > 1.4 )
        {
// cout << "EDL charge density D " << XetaD[ist] << " truncated to +- 0.7 C/m2" <<
//        "  IT= " << pm.IT << " k= " << k << " ist= " << ist << endl;
            XetaD[ist] = XetaD[ist] < 0.0 ? -1.4: 1.4;
            status = 62;
        }

        SigA = 0.;  SigB = 0.;   SigD = 0.;  SigDDL = 0.; // Cleanup 07.12.2009 DK
        pm.XcapD[k][ist] = 0.0;  pm.XpsiD[k][ist] = 0.0;
        if( fabs( XetaA[ist] ) < pm.DHBM  && // pm.lowPosNum*1e6 &&
               fabs( XetaB[ist] ) < pm.DHBM && // pm.lowPosNum*1e6 &&
               fabs( XetaD[ist] ) < pm.DHBM ) // pm.lowPosNum*1e6 )
            goto GEMU_CALC;  // skipping at near-zero charge
        // calculating charge density at diffuse layer
        switch( pm.SCM[k][ist] )
        {
        case SC_CCM:  // Constant-Capacitance Model Schindler, extension Nilsson
            SigA = pm.Xetaf[k][ist] + XetaA[ist];
            SigDDL = -SigA - XetaB[ist];
            SigB = XetaB[ist];
            break;
        case SC_DDLM: // Generalized Double Layer Model [Dzombak and Morel, 1990]
            SigA = pm.Xetaf[k][ist] + XetaA[ist];
            SigDDL = -SigA;
            SigB = 0.0;
            break;
        case SC_TLM:  // Triple-Layer Model [Hayes and Leckie, 1987]
            SigA = pm.Xetaf[k][ist] + XetaA[ist];
            SigB = XetaB[ist];
            SigDDL = -SigA - XetaB[ist];
            break;
        case SC_MTL:  // Modified Triple-Layer Model [Robertson, 1997]
            SigA = pm.Xetaf[k][ist] + XetaA[ist];
            SigB = XetaB[ist];
            SigDDL = -SigA - XetaB[ist];
            break;
        case SC_BSM: // Basic Stern model: [Christl and Kretzschmar, 1999]
            SigA = pm.Xetaf[k][ist] + XetaA[ist];
            SigB = XetaB[ist];
            SigDDL = -SigA - XetaB[ist];
            break;
        case SC_3LM: // Three-layer model: Hiemstra ea 1996; Tadanier & Eick 2002
            SigA = pm.Xetaf[k][ist] + XetaA[ist];
            SigB = XetaB[ist];
            SigD = XetaD[ist];
            SigDDL = -SigA - SigB -SigD;
            break;
        case SC_MXC:  // BSM for Ion-Exchange on permanent-charge surface
            SigA = pm.Xetaf[k][ist] + XetaA[ist];
            SigB = XetaB[ist];
            SigDDL = -SigA - XetaB[ist];
            break;
        case SC_NNE:  // Non-Electrostatic
            break;
        default:
            continue;
        }
//        if( fabs( SigD ) > 1 )
//        {
//cout << "EDL charge density D " << SigD << " truncated to +- 1 C/m2" <<
//        "  IT= " << pm.IT << " k= " << k << " ist= " << ist << endl;
//            SigD = SigD < 0.0 ? -1.0: 1.0;
//        }
        // Gouy-Chapman equation
        // parameters of diffuse layer using [Damaskin, 1987,p.192-195]
        A = 1e-9;
        F2RT = pm.FRT / 2.;
        Sig = SigDDL; //  - XetaW[ist] ;
        I=pm.IC;
        if( I > 1e-7 )
            // Aq solution density Ro included acc. to [Machesky ea., 1999]
            // Only for the basic Stern model for now
        {
            double Ro = 1.0;
            if( pm.SCM[k][ist] == SC_BSM && pm.FVOL[0] > 1e-16 )
                Ro = pm.FWGT[0] / pm.FVOL[0];
            A = sqrt( 2000. * 8.854e-12 * pm.epsW[0] * pm.RT * I * Ro );
        }
        Cap = F2RT * sqrt( 4.*A*A + Sig*Sig );

        // SD: workaround because of problems with log argument
        f3 =  sqrt( 1.+Sig*Sig/(4.*A*A) ) - Sig/(2.*A);
// std::cout<< f1  << ' ' << f3 << endl;
        if( f3 < 1 )
        {
            f1 = exp( -3. * F2RT );
            if( f3<f1) f3 = f1;
        }
        else
        {
            f1 = exp( 3. * F2RT );
            if( f3>f1 ) f3 = f1;
        }
        PSIo = log(f3)/F2RT;
//          PSIo = log( sqrt( 1.+Sig*Sig/(4.*A*A) ) - Sig/(2.*A) ) / F2RT;
//          Cap0 = fabs(Sig/PSIo);
//          Del = A*1e9/(2.*I*F)/cosh(PSIo*F2RT);
//          pm.XdlD[k] = Del;
        pm.XcapD[k][ist] = Cap;
        pm.XpsiD[k][ist] = PSIo;
        PsiD = PSIo;
        // Truncating diffuse plane potential to avoid divergence
        if( fabs( PsiD ) > 0.4 )
        {
// cout << "All EDL models: PsiD = " << PsiD << " truncated to +- 0.4 V" <<
//      "  IT= " << pm.IT << " k= " << k << " ist= " << ist << endl;
            PsiD = PsiD<0? -0.4: 0.4;
            status = 63;
        }
GEMU_CALC:
        // calculating potentials at EDL planes
        switch( pm.SCM[k][ist] )
        {
        case SC_DDLM: // Generalized Double Layer Model [Dzombak & Morel 1990]
            pm.XpsiA[k][ist] = PsiD;
            pm.XpsiB[k][ist] = PsiD;
            break;
        case SC_CCM:  // Constant-Capacitance Model Schindler, ext. Nilsson
            if( pm.XcapB[k][ist] > 0.001 )
            {  // Classic CCM Schindler with inner-sphere species only
               PsiA = SigA / pm.XcapA[k][ist];
               if( fabs( PsiA ) > 0.7 ) // truncated 0-plane potential
               {
                   PsiA = PsiA<0? -0.7: 0.7;
                   status = 64;
               }
               pm.XpsiA[k][ist] = PsiA;
            }
            else { // Extended CCM model [Nilsson ea 1996] as TLM with PsiD = 0
               PsiB = - SigB / pm.XcapB[k][ist];
               if( fabs( PsiB ) > 0.3 )  // truncated B-plane potential
               {
                   PsiB = PsiB<0? -0.3: 0.3;
                   status = 65;
               }
               PsiA = PsiB + SigA / pm.XcapA[k][ist];
               if( fabs( PsiA ) > 0.7 )
               {
                  PsiA = PsiA<0? -0.7: 0.7;
                  status = 66;
               }
               pm.XpsiA[k][ist] = PsiA;
               pm.XpsiB[k][ist] = PsiB;
            }
            break;
        case SC_MTL:  // Modified Triple Layer Model for X- Robertson | Kulik
// PsiD = 0.0; // test
            PsiB = PsiD - SigDDL / pm.XcapB[k][ist];
            if( fabs( PsiB ) > 0.6)  // truncated B-plane potential
            {
// cout << "EDL (MTL) PsiB = " << PsiB << " truncated to +- 0.6 V" <<
//      "  IT= " << pm.IT << " k= " << k << " ist= " << ist << endl;
                PsiB = PsiB<0? -0.6: 0.6;
                status = 67;
            }
            PsiA = PsiB + SigA / pm.XcapA[k][ist];
            if( fabs( PsiA ) > 1.1 )  // truncated 0 plane potential
            {
// cout << "EDL (MTL) PsiA = " << PsiA << " truncated to +- 1.1 V" <<
//      "  IT= " << pm.IT << " k= " << k << " ist= " << ist << endl;
                PsiA = PsiA<0? -1.1: 1.1;
                status = 68;
            }
            pm.XpsiA[k][ist] = PsiA;
            pm.XpsiB[k][ist] = PsiB;
            break;
        case SC_TLM:  // Triple-Layer Model   [Hayes 1987]
            PsiB = PsiD - SigDDL / pm.XcapB[k][ist];
            if( fabs( PsiB ) > 0.6 )  // // truncated B-plane potential
            {
// cout << "EDL (TLM) PsiB = " << PsiB << " truncated to +- 0.6 V" <<
//      "  IT= " << pm.IT << " k= " << k << " ist= " << ist << endl;
                PsiB = PsiB<0? -0.6: 0.6;
                status = 69;
            }
            PsiA = PsiB + SigA / pm.XcapA[k][ist];
            if( fabs( PsiA ) > 1.1 )  // truncated 0-plane potential
            {
// cout << "EDL (TLM) PsiA = " << PsiA << " truncated to +- 1.1 V" <<
//      "  IT= " << k << " k= " << k << " ist= " << ist << endl;
                PsiA = PsiA<0? -1.1: 1.1;
                status = 70;
            }
            pm.XpsiA[k][ist] = PsiA;
            pm.XpsiB[k][ist] = PsiB;
            break;
        case SC_3LM: // Three-Layer Model [Hiemstra & van Riemsdijk 1996]
//            PsiB = PsiD + SigD / pm.XcapB[k][ist];
// cout << "EDL (3LM) PsiB(D) = " << PsiB << "  IT= " << pm.IT << " k= "
// << k << " ist= " << ist << endl;
            PsiB = PsiD + ( SigA + SigB ) / pm.XcapB[k][ist];  // Compare!
// cout << "EDL (3LM) PsiB(AB) = " << PsiB << "  IT= " << pm.IT << " k= "
// << k << " ist= " << ist << endl;
            if( fabs( PsiB ) > 0.6 )  // truncated B-plane potential
            {
// cout << "EDL (3LM) PsiB = " << PsiB << " truncated to +- 0.6 V" <<
//      "  IT= " << pm.IT << " k= " << k << " ist= " << ist << endl;
                PsiB = PsiB<0? -0.6: 0.6;
                status = 71;
            }
            PsiA = PsiB + SigA / pm.XcapA[k][ist];
            if( fabs( PsiA ) > 1.1 )   // truncated 0-plane potential
            {
// cout << "EDL (3LM) PsiA = " << PsiA << " truncated to +- 1.1 V" <<
//      "  IT= " << pm.IT << " k= " << k << " ist= " << ist << endl;
                PsiA = PsiA<0? -1.1: 1.1;
                status = 72;
            }
            pm.XpsiA[k][ist] = PsiA;
            pm.XpsiB[k][ist] = PsiB;
            break;
        case SC_BSM: /* Basic Stern model, Christl & Kretzschmar, 1999 */
            PsiB = PsiD;
            if( fabs( PsiB ) > 0.6 )  // truncated B-plane potential
            {
//cout << "EDL (BSM) PsiB = " << PsiB << " truncated to +- 0.6 V" <<
//      "  IT= " << pm.IT << " k= " << k << " ist= " << ist << endl;
                PsiB = PsiB<0? -0.6: 0.6;
                status = 73;
            }
            PsiA = PsiB + SigA / pm.XcapA[k][ist];
            if( fabs( PsiA ) > 1.1 )  // truncated 0-plane potential
            {
//cout << "EDL (BSM) PsiA = " << PsiA << " truncated to +- 1.1 V" <<
//      "  IT= " << pm.IT << " k= " << k << " ist= " << ist << endl;
                PsiA = PsiA<0? -1.1: 1.1;
                status = 74;
            }
            pm.XpsiA[k][ist] = PsiA;
            pm.XpsiB[k][ist] = PsiB;
            break;
        case SC_MXC:  // BSM for permanent charge surfaces
            PsiB = PsiD;
            if( fabs( PsiB ) > 0.6 )  // truncated B-plane potential
            {
// cout << "EDL (MXC) PsiB = " << PsiB << " truncated to +- 0.6 V" <<
//      "  IT= " << pm.IT << " k= " << k << " ist= " << ist << endl;
                PsiB = PsiB<0? -0.6: 0.6;
                status = 75;
            }
            PsiA = PsiB + SigA / pm.XcapA[k][ist];
            if( fabs( PsiA ) > 1.1 ) // truncated 0-plane potential
            {
// cout << "EDL (MXC) PsiA = " << PsiA << " truncated to +- 1.1 V" <<
//      "  IT= " << pm.IT << " k= " << k << " ist= " << ist << endl;
                PsiA = PsiA<0? -1.1: 1.1;
                status = 76;
            }
            pm.XpsiA[k][ist] = PsiA;
            pm.XpsiB[k][ist] = PsiB;
            break;
        case SC_NNE:  // Non-Electrostatic
            pm.XpsiA[k][ist] = 0.0;
            pm.XpsiB[k][ist] = 0.0;
            break;
        default:
            break;
        }
    }  // ist   end of loop over surface types
    return status;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
///  Calculation of new surface activity coefficient terms SACT (Kulik, 2004)
//
///  Revised by KD in April 2004 (PSI) to introduce new activity
///  coefficient terms SACT rigorously derived from Langmuir and QCA
///  isotherms (Kulik 2006, Radiochimica Acta).
//
///  SACT are placed into pm.lnGam[j], as other activity coefficients except
///  relative surface potentials (Coulombic terms) kept separately.
//
///  Old (obsolete) SAT calculations (Kulik 2000, 2002) are retained.
//
///  pm.lnSAC[*][3] vector is now used to keep original DUL[j] to restore
///  them after IPM-2 refinements for surface complexes.
///    Added 13.03.2008 by DK: returns int value showing (if true)
///    that some extreme values were obtained for some SACTs,
///    0 otherwise (for detecting bad PIA)
//
long int
TMulti::SurfaceActivityCoeff( long int jb, long int je, long int, long int, long int k )
{
	long int status = 0;
        long int i, ii, j, ja, ist=0, iss, dent, Cj, iSite[MST];
    double XS0,  xj0, XVk, XSk, XSkC, xj, Mm, rIEPS, ISAT, XSs,
           SATst, xjn, q1, q2, aF, cN, eF, lnGamjo, lnDiff, lnFactor;
    SPP_SETTING *pa = paTProfil;

    if( pm.XF[k] <= pm.DSM ) // No sorbent retained by the IPM - phase killed
        return status;
    if( pm.XFA[k] <=  pa->p.ScMin )  // elimination of sorption phase
        return status;  // No surface species left

    for(i=0; i<MST; i++)
    {
        iSite[i] = -1;
        for( ii=0; ii<MST; ii++ )
           pm.D[i][ii] = 0.0;  // cleaning the matrix for sites totals
    }
    // Extraction of site indices for the neutral >OH group
    for( j=jb; j<je; j++ )
    {
        ja = j - ( pm.Ls - pm.Lads );
        if( pm.SATT[ja] == SAT_SOLV )
        {
           ist = pm.SATX[ja][XL_ST]; // / MSPN;
           iSite[ist] = j;
        }
        // Counting current sites totals
        if( pm.DCC[j] == DC_SUR_CARRIER ||
            pm.DCC[j] == DC_SUR_MINAL ||
            pm.DCC[j] == DC_PEL_CARRIER ||
            pm.SATT[ja] == SAT_SOLV )
            continue;
        // Calculating ist (index of surface type)
        ist = pm.SATX[ja][XL_ST];  // / MSPN;
        if( ist < 0 || ist >= MST )
            ist = 0;  // default: zero surface type
        // Calculate iss - index of site on surface type
        iss = pm.SATX[ja][XL_SI];
        if( iss < 0 || iss >= MST )
            iss = 0;  // default: zero site is the weekest and the most abundant one
        pm.D[iss][ist] += pm.X[j]; // adding to total amount on a site
    }

    for( j=jb; j<je; j++ )
    { // Main loop for DCs - surface complexes
        lnGamjo = pm.lnGmo[j];             // bugfix 16.03.2008 DK
        if( pm.X[j] < min( pm.DcMinM, pm.lowPosNum ) )
            continue;  // This surface DC has been killed by IPM
//        OSAT = pm.lnGmo[j]; // added 6.07.01 by KDA
        ja = j - ( pm.Ls - pm.Lads );
        rIEPS =  pa->p.IEPS;   // default 1e-3 (for old SAT - 1e-9)
//        dent = 1;  // default - monodentate
        switch( pm.DCC[j] )  // code of species class
        {
        default: // pm.lnGam[j] = 0.0;  not a surface species
            continue;
        case DC_SSC_A0: case DC_SSC_A1: case DC_SSC_A2: case DC_SSC_A3: case DC_SSC_A4:
        case DC_WSC_A0: case DC_WSC_A1: case DC_WSC_A2: case DC_WSC_A3: case DC_WSC_A4:
        case DC_SUR_GROUP: case DC_IEWC_B: case DC_SUR_COMPLEX: case DC_SUR_IPAIR:
        case DC_IESC_A:
            // Calculating ist (index of surface type)
            ist = pm.SATX[ja][XL_ST]; // / MSPN;
            // Calculating iss - index of site on surf.type
            iss = pm.SATX[ja][XL_SI];
            if( iss < 0 || iss >= MST )
              iss = 0;  // zero site is the weekest and the most abundant one
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
                if( XVk < pm.ScMinM ) // pm.DSM*0.1 )
                    continue;    // This end-member is zeroed off by IPM
                Mm = pm.MM[Cj] * XVk/pm.XFA[k];  // mol.mass
            }
            XSk = pm.XFTS[k][ist];  // Total moles of sorbates on surface type
            XSs = pm.D[iss][ist];   // Total moles of SC on site type
            xj = pm.X[j];           // Current moles of this surface species

            // Extracting isotherm parameters
            if( pm.MASDJ )
            {
               cN = pm.MASDJ[ja][PI_P2];  // Frumkin/Pivovarov water coord. number
               if( cN > 0 )
                   dent = cN;   // denticity for L and QCA isoterms
               else dent = 1;   // Cleanup DK 07.12.2009
               aF = pm.MASDJ[ja][PI_P1];  // Frumkin lateral interaction energy term
            //   bet = pm.MASDJ[ja][PI_P3];   // BET beta parameter (reserved)
            }
            else {  // defaults
               cN = 0.0; aF = 0.0; dent = 1; // bet = 1.0;
            }
            switch( pm.SATT[ja] ) // selection of the SACT model
            {
            case SAT_L_COMP: // Competitive monodentate Langmuir on a surface and site type
                XSkC = XSs / XVk / Mm * 1e6 /pm.Nfsp[k][ist]/
                      pm.Aalp[k]/1.66054;  // per nm2
                XS0 = (fabs(pm.MASDJ[ja][PI_DEN])/pm.Aalp[k]/1.66054);
                        // max. density per nm2
                if( pa->p.PC <= 2 )
                    rIEPS = pa->p.IEPS * XS0;   // relative IEPS
                if( XSkC < 0.0 )
                    XSkC = 0.0;
                if( XSkC >= XS0 )               // Setting limits
                    XSkC = XS0 - 2.0 * rIEPS;
                q1 = XS0 - XSkC;
                if( (pa->p.PC == 3 && !pm.W1) || pa->p.PC != 3 )
                {
                  q2 = rIEPS * XS0;
                  if( q1 > q2 )
                      q2 = q1;
                }
                else {
                   q2 = q1;
                   if( q2 <= 1e-22 )
                       q2 = 1e-22;
                }
                ISAT = log( XS0 ) - log( q2 );
                pm.lnGam[j] = ISAT;
                pm.lnSAC[ja][0] = ISAT;
                break;
              // (Non)competitive QCA-L for 1 to 4 dentate species
//            case SAT_QCA4_NCOMP: dent++;     // code '4'
//            case SAT_QCA3_NCOMP: dent++;     // code '3'
//            case SAT_QCA2_NCOMP:             // code '2'
//            case SAT_QCA1_NCOMP:             // code '1'
            case SAT_QCA_NCOMP:  // dent++;   bidentate is default for QCA
                xj0 =
                 (fabs(pm.MASDJ[ja][PI_DEN])/pm.Aalp[k]/1.66054);
                                             // Max site density per nm2
                xj = XSs / XVk / Mm / pm.Nfsp[k][ist] * 1e6     // xj
                     /pm.Aalp[k]/1.66054; // Density per nm2 on site type iss
                if( pa->p.PC <= 2 )
                    rIEPS = pa->p.IEPS * xj0; // relative IEPS
                if(xj >= xj0/(double)dent)
                     xj = xj0/(double)dent - rIEPS;  // upper limit
//                ISAT = 0.0;
                q2 = xj0 - xj*dent;  // Computing differences in QCA gamma
                q1 = xj0 - xj;
                if( q1 < 1e-22 )
                    q1 = 1e-22;
                if( q2 < 1e-22 )
                    q2 = 1e-22;
                ISAT = log(xj0) + log(q1)*(dent-1) - log(q2)*dent;
                pm.lnGam[j] = ISAT;
                pm.lnSAC[ja][0] = ISAT;
                break;
            case SAT_FRUM_COMP: // Frumkin (FFG) isotherm
                dent = 1; // monodentate for now
                xj0 =
                 (fabs(pm.MASDJ[ja][PI_DEN])/pm.Aalp[k]/1.66054);
                                             // Max site density per nm2
                xj = XSs / XVk / Mm / pm.Nfsp[k][ist] * 1e6  //  xj
                     /pm.Aalp[k]/1.66054; // Current density per nm2
                if( pa->p.PC <= 2 )
                    rIEPS = pa->p.IEPS * xj0; // relative IEPS
                if(xj >= xj0/dent)
                     xj = xj0/dent - rIEPS;  // upper limit
                q2 = xj0 - xj*dent;  // Computing differences in QCA gamma
                q1 = xj0 - xj;
                if( q1 < 1e-22 )
                    q1 = 1e-22;
                if( q2 < 1e-22 )
                    q2 = 1e-22;
                ISAT = log(xj0) + log(q1)*(dent-1) - log(q2)*dent;
                // Calculation of the Frumkin exponential term
                if( fabs (aF) < 1e-9 || fabs (cN) < 1e-9 )
                   eF = 0.0;
                else {  // changed from cN to -cN on 16.09.2005 by KD
                   eF = -cN * aF * xj*dent / xj0 ;  // Fi = Fi'/(kT) Bockris p.938
                }
                pm.lnGam[j] = ISAT + eF;
                pm.lnSAC[ja][0] = ISAT;
                pm.lnSAC[ja][1] = eF;
                break;
            case SAT_FRUM_NCOMP: // (Non)competitive Frumkin isotherm
                                 // for permanent charge surfaces
                XSkC = xj / XVk / Mm / pm.Nfsp[k][ist] * 1e6
                       / pm.Aalp[k]/1.66054;  // per nm2
                XS0 = (pm.MASDJ[ja][PI_DEN]/pm.Aalp[k]/1.66054);
                         // max.dens.per nm2
                if( pa->p.PC <= 2 )
                    rIEPS = pa->p.IEPS * XS0;  // relative IEPS
                if( XSkC < 0.0 )
                    XSkC = 0.0;
                if( XSkC >= XS0 )  // Limits
                    XSkC = XS0 - 2.0 * rIEPS;
                q1 = XS0 - XSkC;
                if(( pa->p.PC == 3 && !pm.W1) || pa->p.PC != 3 )
                {
                  q2 = rIEPS * XS0;
                  if( q1 > q2 )
                      q2 = q1;
                }
                else {
                   q2 = q1;
                   if( q2 <= 1e-22 )
                       q2 = 1e-22;
                }
                ISAT = log( XS0 ) - log( q2 );
                // Calculation of the Frumkin exponential term (competitive)
                if( fabs (aF) < 1e-9 || fabs (cN) < 1e-9 )
                   eF = 0.0;
                else { // changed from cN to -cN on 16.09.2005 by KD
                   eF = -cN * aF * XSkC / XS0;  // Fi = Fi'/(kT) Bockris p.938
                }
                pm.lnGam[j] = ISAT + eF;
                pm.lnSAC[ja][0] = ISAT;
                pm.lnSAC[ja][1] = eF;
                break;
            case SAT_PIVO_NCOMP: // (Non)competitive Pivovarov isotherm
                dent = 1; // monodentate for now
                xj0 =
                 (fabs(pm.MASDJ[ja][PI_DEN])/pm.Aalp[k]/1.66054);
                                             // Max site density per nm2
                xj = XSs / XVk / Mm / pm.Nfsp[k][ist] * 1e6
                     /pm.Aalp[k]/1.66054; // Current density per nm2
                if( pa->p.PC <= 2 )
                    rIEPS = pa->p.IEPS * xj0; // relative IEPS
                if(xj >= xj0/dent)
                     xj = xj0/dent - rIEPS;  // upper limit
                q2 = xj0 - xj*dent;  // Computing differences in gamma
                q1 = xj0 - xj;
                if( q1 < 1e-22 )
                    q1 = 1e-22;
                if( q2 < 1e-22 )
                    q2 = 1e-22;
                ISAT = log(xj0) + log(q1)*(dent-1) - log(q2)*dent;
               // Calculation of the Frumkin exponential term
                if( fabs (aF) < 1e-9 || fabs (cN) < 1e-9 )
                   eF = 0.0;
                else {
                   double pivovar;
                   eF = cN * aF * xj / xj0 ;  // Fi = Fi'/(kT) Bockris p.938
            // Calculation of the Pivovarov 98 exponential correction term
                   pivovar = xj0 / ( xj0 + xj * ( cN -1 ));
                   eF *= pivovar;
                }
                pm.lnGam[j] = ISAT + eF;
                pm.lnSAC[ja][0] = ISAT;
                pm.lnSAC[ja][1] = eF;
                break;
            case SAT_VIR_NCOMP: // Non-Competitive virial isotherm
                dent = 1; // monodentate for now
                xj0 =
                 (fabs(pm.MASDJ[ja][PI_DEN])/pm.Aalp[k]/1.66054);
                                             // Max site density per nm2
                xj = XSs / XVk / Mm / pm.Nfsp[k][ist] * 1e6
                     /pm.Aalp[k]/1.66054; // Current density per nm2
                if( pa->p.PC <= 2 )
                    rIEPS = pa->p.IEPS * xj0; // relative IEPS
                if(xj >= xj0/dent)
                     xj = xj0/dent - rIEPS;  // upper limit
                ISAT = 0.0;
                if( fabs (aF) < 1e-9 || fabs (cN) < 1e-9 )
                   eF = 0.0;
                else {   // changed from cN to -cN on 16.09.2005 by KD
                   eF = -cN * aF * xj / xj0 ;  // Fi = Fi'/(kT) Bockris p.938
                }
                pm.lnGam[j] = ISAT + eF;
                pm.lnSAC[ja][0] = ISAT;
                pm.lnSAC[ja][1] = eF;
                break;
            case SAT_BET_NCOMP: // Non-competitive BET for surface precipitation
                ISAT = 0.0;
// To be completed
//
//
                pm.lnGam[j] = ISAT;
                pm.lnSAC[ja][0] = ISAT;
                break;
            case SAT_INDEF: // No SAT calculation whatsoever
                pm.lnGam[j] = 0.0;
                pm.lnSAC[ja][0] = 0;
                break;
            default:        // pm.lnGam[j] = 0.0;
                break;
//  Obsolete old SAT calculations (retained from previous versions)
            case SAT_COMP: // Competitive SAT (obsolete) on a surface type
                if( iSite[ist] < 0 )
                    xjn = 0.0;
                else xjn = pm.X[iSite[ist]]; // neutral site does not compete!
                XS0 = pm.MASDT[k][ist] * XVk * Mm / 1e6
                      * pm.Nfsp[k][ist]; // expected total in moles
                if( pa->p.PC <= 2 )
                    rIEPS = pa->p.IEPS * XS0;  // relative IEPS
                XSkC = XSk - xjn - xj; // occupied by the competing species;
                                     // this sorbate cannot compete to itself
                if( XSkC < 0.0 )
                    XSkC = 0.0;
                if( XSkC >= XS0 )  // Limits
                    XSkC = XS0 - 2.0 * rIEPS;
                xj0 = XS0 - XSkC;    // expected moles of this sorbate
                if(xj >= xj0)
                       xj = xj0 - rIEPS; // Limits
                if( xj * 2 <= xj0 )
                    ISAT = 0.0;      // ideal case
                else
                {
                   q1 = xj0 - xj;
                   q2 = rIEPS * XS0;
                   if( (pa->p.PC == 3 && !pm.W1) || pa->p.PC != 3 )
                   {
                      if( q1 > q2 )
                        q2 = q1;
                   }
                   else {
                      q2 = q1;
                      if( q2 <= 1e-33 )
                         q2 = 1e-33;
                   }
                   ISAT = log( xj ) - log( q2 );
                }
                pm.lnGam[j] = ISAT;
                pm.lnSAC[ja][0] = ISAT;
                break;
            case SAT_NCOMP: // Non-competitive truncated Langmuir SAT (obsolete)
                // rIEPS = pa->p.IEPS * 2;
                xj0 = fabs( pm.MASDJ[ja][PI_DEN] ) * XVk * Mm / 1e6
                      * pm.Nfsp[k][ist]; // in moles
                if( pa->p.PC <= 2 )
                    rIEPS = pa->p.IEPS * xj0;  // relative IEPS
                if(xj >= xj0)
                     xj = xj0 - rIEPS;  // upper limit
                if( xj * 2.0 <= xj0 )   // Linear adsorption - to improve !
                    ISAT = 0.0;
                else
                {
                    q1 = xj0 - xj;      // limits: rIEPS to 0.5*xj0
                    q2 = xj0 * rIEPS;
                    if( pa->p.PC == 3 && pm.W1 )
                       ISAT = log( xj ) - log( q1 );
                    else {
                       if( q1 > q2 )
                          ISAT = log( xj ) - log( q1 );
                       else
                          ISAT = log( xj ) - log( q2 );
                    }
                 }
                pm.lnGam[j] = ISAT;
                pm.lnSAC[ja][0] = ISAT;
                break;
            case SAT_SOLV:  // Neutral surface site (e.g. >O0.5H@ group)
                            // applies to the whole surface type!
                XSs = 0.0;  // calc total moles on all sites on surface type
                for( i=0; i<MST; i++ )
                   XSs += pm.D[i][ist];
                XSkC = XSs / XVk / Mm * 1e6  // total non-solvent surf.species
                   /pm.Nfsp[k][ist]/ pm.Aalp[k]/1.66054;  // per nm2
                XS0 = (max( pm.MASDT[k][ist], pm.MASDJ[ja][PI_DEN] ));
                SATst = pa->p.DNS*1.66054*pm.Aalp[k]/XS0;
                XS0 = XS0 / pm.Aalp[k]/1.66054;
                if( pa->p.PC <= 2 )
                    rIEPS = pa->p.IEPS * XS0;  // relative IEPS
                if( XSkC < 0.0 )
                    XSkC = 0.0;
                if( XSkC >= XS0 )  // Limits
                    XSkC = XS0 - 2.0 * rIEPS;
                q1 = XS0 - XSkC;
                if( (pa->p.PC == 3 && !pm.W1) || pa->p.PC != 3 )
                {
                  q2 = rIEPS * XS0;
                  if( q1 > q2 )
                      q2 = q1;
                }
                else {
                   q2 = q1;
                   if( q2 <= 1e-22 )
                       q2 = 1e-22;
                }
                ISAT = log( XS0 ) - log( q2 );
                pm.lnGam[j] = ISAT + log( SATst );
                pm.lnSAC[ja][0] = ISAT;
                pm.lnSAC[ja][1] = log(SATst);
                break;
            }
        }

        if( lnGamjo > pm.DHBM )
        {                               // Workaround DK 07.12.2009
            lnFactor = 0.2;
            lnDiff = pm.lnGam[j] - lnGamjo;
            if( fabs( lnDiff ) > fabs( lnFactor ) ) // e times
            {
                if( fabs( lnDiff ) > 6.907755 )  // 1000 times
                    status = 101;   // the SACT has changed too much;
                            // the threshold pa->p.IEPS needs adjustment!
                // Smoothing (provisional)
                pm.lnGam[j] = lnDiff > 0? lnGamjo + lnDiff - lnFactor:
                                lnGamjo + lnDiff + lnFactor;
            }
        }
    }  // j
   return status;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculating demo partial pressures of gases (works only in GEMS-PSI)
//
void TMulti::GasParcP()
{
#ifndef IPMGEMPLUGIN

        long int k,  i, jj=0;
    long int jb, je, j;

    if( !pm.PG )
        return;

  char (*SMbuf)[MAXDCNAME] =
      (char (*)[MAXDCNAME])aObj[ o_w_tprn].Alloc( pm.PG, 1, MAXDCNAME );
  pm.Fug = (double *)aObj[ o_wd_fug].Alloc( pm.PG, 1, D_ );
  pm.Fug_l = (double *)aObj[ o_wd_fugl].Alloc( pm.PG, 1, D_ );
  pm.Ppg_l = (double *)aObj[ o_wd_ppgl].Alloc( pm.PG, 1, D_ );

    for( k=0, je=0; k<pm.FIs; k++ ) // phase
    {
        jb = je;
        je = jb+pm.L1[k];
        if( pm.PHC[k] == PH_GASMIX || pm.PHC[k] == PH_PLASMA
           || pm.PHC[k] == PH_FLUID )
        {
            for( j=jb; j<je; j++,jj++ )
            {  // fixed 02.03.98 DK

                copyValues(SMbuf[jj], pm.SM[j], MAXDCNAME );
                pm.Fug_l[jj] = -(pm.G0[j] + pm.fDQF[j]);
                if( pm.Pc > 1e-9 )
                    pm.Fug_l[jj] += log(pm.Pc);
                for( i=0; i<pm.N; i++ )
                    pm.Fug_l[jj] += *(pm.A+j*pm.N+i) * pm.U[i];
                if( pm.Fug_l[jj] > -37. && pm.Fug_l[jj] < 16. )
                    pm.Fug[jj] = exp( pm.Fug_l[jj] );
                else  pm.Fug[jj] = 0.0;
                // Partial pressure
                pm.Ppg_l[jj] = pm.Fug_l[jj] - pm.lnGam[j];
                pm.Fug_l[jj] *= .43429448;
                pm.Ppg_l[jj] *= .43429448;
            }
            // break;
        }
    }
#endif
}

//--------------------- End of ipm_chemical2.cpp ---------------------------
