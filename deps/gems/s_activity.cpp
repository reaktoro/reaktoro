//-------------------------------------------------------------------
// $Id: s_activity.cpp 771 2012-12-13 13:07:43Z kulik $
//
/// \file s_activity.cpp
/// Implementation of chemistry-specific functions (concentrations,
/// activity coefficients, adsorption models etc.)
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

#include "m_param.h"
#include "activities.h"

//#ifndef IPMGEMPLUGIN
//#include "service.h"
//#include "stepwise.h"
//#endif


void TActivity::setTemperature(double T)
{
    ;
}

// set pressure (in units of Pa)
void TActivity::setPressure(double P)
{
    ;
}

// compute standard Gibbs energies (so far only P,T interpolation)
void TActivity::updateStandardGibbsEnergies()
{
    ;
}

void TActivity::updateStandardVolumes()
{
    ;
}

void TActivity::updateStandardEnthalpies()
{
    ;
}

void TActivity::updateStandardEntropies()
{
    ;
}

//void TActivity::updateThermoData()
//{
//    ;
//}

// set speciation (in units of moles)
void TActivity::setSpeciation( const double* n )
{
 ;

}

// compute concentrations in all phases
void TActivity::updateConcentrations()
{
    ;
}

// compute activity coefficients
void TActivity::updateActivityCoefficients()
{
;
}

// compute primal chemical potentials
void TActivity::updateChemicalPotentials()
{
;
}

// compute primal activities
void TActivity::updateActivities()
{
    ;
}

void TActivity::updateChemicalData()
{
    ;
}




//-------------------------------------------------------------
// setting defaults in TActivity class data
void TActivity::set_def( )
{
    ;
}
// #define GEMITERTRACE

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Internal memory allocation for IPM performance optimization
/// (since version 2.2.0)
//
void TActivity::Alloc_A_B( long int newN )
{
  if( AA && BB && (newN == sizeN) )
    return;
  Free_A_B();
  AA = new  double[newN*newN];
  BB = new  double[newN];
  sizeN = newN;
}

void TActivity::Free_A_B()
{
  if( AA  )
    { delete[] AA; AA = 0; }
  if( BB )
    { delete[] BB; BB = 0; }
  sizeN = 0;
}

#define  a(j,i) ((*(act.A+(i)+(j)*act.N)))
/// Building an index list of non-zero elements of the matrix act.A
void TActivity::Build_compressed_xAN()
{
 long int ii, jj, k;

 // Calculate number of non-zero elements in A matrix
 k = 0;
 for( jj=0; jj<act.L; jj++ )
   for( ii=0; ii<act.N; ii++ )
     if( fabs( a(jj,ii) ) > 1e-12 )
       k++;

   // Free old memory allocation
    Free_compressed_xAN();

   // Allocate memory
   arrL = new long int[act.L+1];
   arrAN = new long int[k];

   // Set indexes in the index arrays
   k = 0;
   for( jj=0; jj<act.L; jj++ )
   { arrL[jj] = k;
     for( ii=0; ii<act.N; ii++ )
       if( fabs( a(jj,ii) ) > 1e-12 )
       {
        arrAN[k] = ii;
        k++;
       }
   }
   arrL[jj] = k;
}
#undef a

void TActivity::Free_compressed_xAN()
{
  if( arrL  )
    { delete[] arrL; arrL = 0;  }
  if( arrAN )
    { delete[] arrAN; arrAN = 0;  }
}

void TActivity::Free_internal()
{
  Free_compressed_xAN();
  Free_A_B();
 }

/// Internal memory allocation for IPM performance optimization
void TActivity::Alloc_internal()
{
// optimization 08/02/2007
 Alloc_A_B( act.N );
 Build_compressed_xAN();
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculating value of dual chemical potential of j-th dependent component
///     performance optimized version  (February 2007)
double TActivity::DC_DualChemicalPotential( double U[], double AL[], long int N, long int j )
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

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculating total amounts of phases
//
void TActivity::TotalPhasesAmounts( double X[], double XF[], double XFA[] )
{
    long int jj, j, i, k;
    double XFw, XFs, x;

    j=0;
    for( k=0; k< act.FI; k++ )
    { // cycle by phases
        i=j+act.L1[k];
        XFw = 0.0;
        XFs=0.0; // calculating mole amount of carrier (solvent/sorbent)
        for(jj=j; jj<i; jj++)
        {
            x = X[jj];
            if( act.DCCW[jj] == DC_ASYM_CARRIER && act.FIs )
                XFw += x;
            else XFs += x;
        }
        XF[k] = XFw + XFs;
        if( k < act.FIs )
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
double TActivity::DC_PrimalChemicalPotentialUpdate( long int j, long int k )
{
    long int ja=0, ist, isp, jc=-1;
    double F0=0.0, Fold, dF0, Mk=0.0, Ez, psiA, psiB, CD0, CDb, ObS;
    double FactSur, FactSurT;
//    SPP_SETTING *pa = paTProfil;

    Fold = act.F0[j];    // No old sorption models in this implementation!
//    if( act.FIat > 0 && j < act.Ls && j >= act.Ls - act.Lads )
//    {
//        ja = j - ( act.Ls - act.Lads );
//        jc = act.SATX[ja][XL_EM];
// Workaround 08.02.2011 by DK to prevent crash on inconsistent sorption DC
//        if( jc >= act.L1[k] )
//            jc = -1;
//    }
    if( k < act.FIs && act.XFA[k] > 1e-12)
    {
           if( jc < 0 ) // phase (carrier) molar mass g/mkmol
              Mk = act.FWGT[k]/act.XFA[k]*1e-6;
           else Mk = act.MM[jc]*(act.X[jc]/act.XFA[k])*1e-6;
        // DC carrier molar mass g/mkmol
    }
    switch( act.DCC[j] )
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
        F0 = act.lnGmM[j];
        break;
/*
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
        if( !act.Lads || !act.FIat )   // Foolproof - to prevent crash DK 22.12.2009
            break;
        F0 = act.lnGmM[j]; / * + act.lnGam[j]; * /
        // get ist - index of surface type and isp - index of surface plane
/ *!!!!!* /  ist = act.SATX[ja][XL_ST];  // / MSPN;
/ *!!!!!* /  isp = act.SATX[ja][XL_SP]; // % MSPN;
        CD0 = act.MASDJ[ja][PI_CD0];  // species charge that goes into 0 plane
        CDb = act.MASDJ[ja][PI_CDB];  // species charge that goes into B plane
        ObS = act.MASDJ[ja][PI_DEN];  // obsolete - the sign for outer-sphere charge
        if( ObS >= 0.0 )
            ObS = 1.0;
        else ObS = -1.0;
        psiA = act.XpsiA[k][ist];
        psiB = act.XpsiB[k][ist];
        Ez = double(act.EZ[j]);
        if( !isp )  // This is the 0 (A) plane species
        {
            if( fabs( CD0 ) > 1e-20 )  // Doubtful...
                Ez = CD0;
            F0 += psiA * Ez * act.FRT;
        }
        else  // This is B plane
        {
            if( act.SCM[k][ist] == SC_MTL || act.SCM[k][ist] == SC_MXC )
            { // Modified TL: Robertson, 1997; also XTLM Kulik 2002
                  if( fabs( CDb ) > 1e-20 )  // Doubtful...
                      Ez = CDb;
                  F0 += psiB * Ez * act.FRT;
            }
            if( act.SCM[k][ist] == SC_TLM )
            {
// New CD version of TLM  added 25.10.2004
               if( fabs( CD0 ) > 1e-20 && fabs( CDb ) > 1e-20 )
                  F0 += ( psiA*CD0 + psiB*CDb )* act.FRT;
             // see also Table 4 in Zachara & Westall, 1999
             // Old version:  TLM Hayes & Leckie, 1987 uses the sign indicator at density
                else {
                  if( ObS < 0 )
                  {
                      Ez -= 1.0;
                      F0 += ( psiA + Ez * psiB )* act.FRT;
                  }
                  else
                  {
                      Ez += 1.0;
                      F0 += ( Ez * psiB - psiA )* act.FRT;
                  }
               }
            }
            else if( act.SCM[k][ist] == SC_BSM || act.SCM[k][ist] == SC_CCM )
            { // Basic Stern model, Christl & Kretzschmar, 1999
            // New CD version of TLM  added 25.10.2004
               if( fabs( CD0 ) > 1e-20 && fabs( CDb ) > 1e-20 )
                  F0 += ( psiA*CD0 + psiB*CDb )* act.FRT;
                else {
                  if( ObS < 0 )
                  {
                      Ez -= 1.0;
                      F0 += ( psiA + Ez * psiB )* act.FRT;
                  }
                  else
                  {
                      Ez += 1.0;
                      F0 += ( Ez * psiB - psiA )* act.FRT;
                  }
               }
            }
        }
        if( Mk > 1e-9 )  // Mk is carrier molar mass in g/mkmol
        {   // Correction for standard density, surface area and surface type fraction
                FactSur = Mk * (act.Aalp[k]) * pa->p.DNS*1.66054;
        	    // FactSur is adsorbed mole amount at st. surf. density per mole of solid carrier
                FactSurT = FactSur * (act.Nfsp[k][ist]);
                if( act.SCM[k][ist] == SC_MXC || act.SCM[k][ist] == SC_NNE ||
                    act.SCM[k][ist] == SC_IEV )
                                                // F0 -= log( Mk * (act.Nfsp[k][ist]) *
                                                // (act.Aalp[k]) * pa->p.DNS*1.66054 );
                  F0 -= log( FactSurT );
            else  F0 -= log( FactSurT );
                                // F0 -= log( Mk * (act.Nfsp[k][ist]) *
                                // (act.Aalp[k]) * pa->p.DNS*1.66054 );
            F0 -= FactSur / ( 1. + FactSur );
                                // F0 -= (act.Aalp[k])*Mk*pa->p.DNS*1.66054 /
                                // ( 1.0 + (act.Aalp[k])*Mk*pa->p.DNS*1.66054 );
        }
        break;
    case DC_PEL_CARRIER:
    case DC_SUR_MINAL:  // constant charge of carrier - not completed
    case DC_SUR_CARRIER: // Mk is carrier molar mass in g/mkmol
        if( !act.Lads || !act.FIat )   // Foolproof - to prevent crash DK 22.12.2009
            break;
        FactSur = Mk * (act.Aalp[k]) * pa->p.DNS*1.66054;
        F0 -= FactSur / ( 1. + FactSur );
        F0 += FactSur;
        break;
*/  default:
        break;
    }
    F0 += act.lnGam[j];
// No smoothing factor either (internal in GEMIPM3)
//    if( k >= act.FIs )
        return F0;
    // Smoothing procedure for highly non-ideal systems
//    if( act.sMod[k][SGM_MODE] != SM_IDEAL )  // check this condition for sublattice SS models!
            // || act.sMod[k][SCM_TYPE] != SC_NNE )  // changed, 14.07.2009 (TW)
//    {
//        double SmoSensT = 1e-5;   // to be adjusted
//        dF0 = F0 - Fold;
//        if( act.X[j] > min( act.lowPosNum, act.DcMinM ) && fabs( dF0 ) >= SmoSensT )
//       	    F0 = Fold + dF0 * SmoothingFactor();    // Changed 18.06.2008 DK
//    }
//    return F0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of DC primal chemical potential F.
/// From moles of DC Y[], total moles of phase YF[] and DC partial
/// molar Gibbs energy gT (obtained from act.G[]) which includes
/// activity coefficient terms.
/// On error returns F = +7777777.
double TActivity::DC_PrimalChemicalPotential(
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
TActivity::PrimalChemicalPotentials( double F[], double Y[], double YF[], double YFA[] )
{
    long int i,j,k;
    double NonLogTerm=0., v, Yw, Yf, YFk, logXw, logYFk, aqsTail; // v is debug variable

    for( j=0; j<act.L; j++)
       F[j] =0;

    j=0;
    for( k=0; k<act.FI; k++ )
    { // loop over phases
        i=j+act.L1[k];
        if( act.L1[k] == 1L && YF[k] < act.PhMinM )
        	goto NEXT_PHASE;
        if( YF[k] <= act.DSM || ( act.PHC[k] == PH_AQUEL &&
            ( YF[k] <= act.DSM || Y[act.LO] <= act.XwMinM )))
            goto NEXT_PHASE;

        YFk = 0.0;
        Yf= YF[k]; // calculate number of moles of carrier
        if( act.FIs && k<act.FIs )
            YFk = YFA[k];
        if( Yf >= 1e6 )
        {                 // error - will result in zerodivide!
           string pbuf(act.SF[k],0,20);
           char buf[200];
           sprintf( buf, "Broken phase amount from primal approximation: Phase %s  Yf= %lg", pbuf.c_str(), Yf );
           Error( "E13IPM: PrimalChemicalPotentials():", buf);
//           Yf = act.YFk;
        }
        if( ( act.PHC[k] == PH_AQUEL && YFk >= act.XwMinM )
                        || ( act.PHC[k] == PH_SORPTION && YFk >= act.ScMinM )
                        || ( act.PHC[k] == PH_POLYEL && YFk >= act.ScMinM ) )
        {
            logXw = log(YFk);
            NonLogTerm = 1.- YFk / Yf;
#ifdef NOMUPNONLOGTERM
NonLogTerm = 0.0;
#endif
        }
        if( act.L1[k] > 1 )
        {
            logYFk = log( Yf );
        }
        if( act.PHC[k] == PH_AQUEL )
        {    // ln moles of solvent in aqueous phase
            Yw = YFk;
            aqsTail = NonLogTerm;
        }
        for( ; j<i; j++ )
        { //  cycle by DC
            if( Y[j] < act.DcMinM )
                continue;  // exception by minimum DC quantity
                           // calculate chemical potential of j-th DC
            v = DC_PrimalChemicalPotential( act.G[j], log(Y[j]), logYFk,
                              NonLogTerm, logXw, act.DCCW[j] );
            F[j] = v;
       }   // j
NEXT_PHASE:
        j = i;
    }  // k
//    if( Yw >= DSM ) // act.lowPosNum*1e3 )
//       logXw = log(Yw);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of a species contribution to the total Gibbs energy G(X)
/// of the system. On error returns +7777777.
//
double TActivity::DC_GibbsEnergyContribution(
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
double TActivity::GX( double LM  )
{
    long int i, j, k;
    double x, XF, XFw, FX, Gi; // debug variable
//    double const1= act.lowPosNum*10.,
//           const2 = act.lowPosNum*1000.;
/*
    if( LM < act.lowPosNum )     // copy vector Y into X
        for(i=0;i<act.L;i++)
            act.X[i]=act.Y[i];
    else  // calculate new values of X
        for(i=0;i<act.L;i++ )
        {  // gradient vector act.MU - the direction of descent!
            act.X[i]=act.Y[i]+LM*act.MU[i];
//            if( act.X[i] <  act.lowPosNum )   // this is the Ls set cutoff !!!!!!!!!!
            if( act.X[i] <  act.DcMinM )
                act.X[i]=0.;
        }
    // calculate new total quantities of phases
    TotalPhasesAmounts( act.X, act.XF, act.XFA );
*/
    // calculating G(X)
    FX=0.;
    j=0;
    for( k=0; k<act.FI; k++ )
    { // loop for phases
        i=j+act.L1[k];
        act.logXw = -101.;
        XFw = 0.0;  // calculating mole amount of the solvent/sorbent
        if( act.FIs && k<act.FIs )
            XFw = act.XFA[k];
 //       if( XFw > const1 )
        if( ( act.PHC[k] == PH_AQUEL && XFw >= act.XwMinM )
                        || ( act.PHC[k] == PH_SORPTION && XFw >= act.ScMinM )
                        || ( act.PHC[k] == PH_POLYEL && XFw >= act.ScMinM ) )
             act.logXw = log( XFw );

        XF = act.XF[k];
        if( !(act.FIs && k < act.FIs) )
        {
                if( XF < act.PhMinM )
        		goto NEXT_PHASE;
        }
        else if( XF < act.DSM && act.logXw < -100. )
        	goto NEXT_PHASE;

        act.logYFk = log( XF );

        for( ; j<i; j++ )
        { // DCs (species)
            x = act.X[j];
            if( x < act.DcMinM )
                continue;
            // calculating increment of G(x)
            // Gi = DC_GibbsEnergyContribution( act.G[j], x, act.logYFk, act.logXw,
            //                     act.DCCW[j] );
            // call replaced here by inline variant for higher performance
            switch( act.DCCW[j] )
            {
             case DC_ASYM_SPECIES:
                    Gi = x * ( act.G[j] + log(x) - act.logXw );
                    break;
            case DC_ASYM_CARRIER:
            case DC_SYMMETRIC:
                   Gi = x * ( act.G[j] + log(x) - act.logYFk );
                   break;
            case DC_SINGLE:
                   Gi = act.G[j] * x;
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

//#ifndef IPMGEMPLUGIN
/*
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Variant of GX() function for use in the UnSpace module (non-optimized)
// Should not be called from within GEMIPM!
//
double TActivity::pb_GX( double *Gxx  )
{
    long int i, j, k;
    double Gi, x, XF, XFw, FX;
    SPP_SETTING *pa = paTProfil;

    // calculating G(X)
    FX=0.;
    j=0;
    for( k=0; k<act.FI; k++ )
    { // phase loop
        i=j+act.L1[k];
        act.logXw = -101.;
        XFw = 0.0;  // calculating mole amount of the solvent/sorbent
        if( act.FIs && k<act.FIs )
            XFw = act.XFA[k];
 //       if( XFw > const1 )
        if( ( act.PHC[k] == PH_AQUEL && XFw >= pa->p.XwMin )
                        || ( act.PHC[k] == PH_SORPTION && XFw >= pa->p.ScMin )
                        || ( act.PHC[k] == PH_POLYEL && XFw >= pa->p.ScMin ) )
             act.logXw = log( XFw );
        / *   * /
        XF = act.XF[k];
        if( !(act.FIs && k < act.FIs) )
        {
                if( XF < pa->p.PhMin )
        		goto NEXT_PHASE;
        }
        else if( XF < pa->p.DS && act.logXw < 100. )
        	goto NEXT_PHASE;
        act.logYFk = log( XF );

        for( ; j<i; j++ )
        { // DC loop
            x = act.X[j];
            if( x < pa->p.DcMin )
                continue;
            // calculating DC increment to G(x)
            Gi = DC_GibbsEnergyContribution( Gxx[j], x, act.logYFk, act.logXw,
                                 act.DCCW[j] );
            FX += Gi;
        }   // j
NEXT_PHASE:
        j = i;
    }  // k
    return(FX);
}
//#endif
*/
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Conversion of g(T,P) value for DCs into the uniform cj scale.
/// \param k - index of phase, \param j - index DC in phase
/// \return if error code, returns 777777777.
//
double TActivity:: ConvertGj_toUniformStandardState( double g0, long int j, long int k )
{
    double G, YOF=0;

    G = g0/act.RT;
    if( act.YOF )
        YOF = act.YOF[k];     // should be already normalized (J/mol/RT)
    // Calculation of standard concentration scaling terms
    switch( act.DCC[j] )
    { // Aqueous electrolyte
    case DC_AQ_PROTON:
    case DC_AQ_ELECTRON:
    case DC_AQ_SPECIES:
case DC_AQ_SURCOMP:
        G += act.ln5551;
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
        if( act.PHC[k] == PH_GASMIX || act.PHC[k] == PH_FLUID
            || act.PHC[k] == PH_PLASMA )
        {
//        if( act.Pparc[j] != 1.0 && act.Pparc[j] > 1e-30 )
//           G += log( act.Pparc[j] ); // log partial pressure/fugacity
//        else
               G += log( act.Pc ); // log general pressure (changed 04.12.2006)
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
        G += act.ln5551;
        break;
    default: // error - returning 7777777
        return 7777777.;
    }
    return G;
}

/// Converting DC class codes into generic internal codes of IPM
//
void TActivity::ConvertDCC()
{
    long int i, j, k, iRet=0;
    char DCCW;

    j=0;
    for( k=0; k< act.FI; k++ )
    { // phase loop
        i=j+act.L1[k];
        if( act.L1[k] == 1 )
        {
            act.DCCW[j] = DC_SINGLE;
            goto NEXT_PHASE;
        }
        for( ; j<i; j++ )
        { // DC loop
            switch( act.DCC[j] ) // select v_j expression
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
                act.DCC[j] = DC_SSC_A0;
                break;
            case DC_SUR_IPAIR:
                DCCW = DC_ASYM_SPECIES;
                act.DCC[j] = DC_WSC_A0;
                break;
            case DC_SUR_MINAL:
            case DC_SUR_CARRIER:
            case DC_PEL_CARRIER:
                 DCCW = DC_ASYM_CARRIER;
                break;
            default:
                if( isdigit( act.DCC[j] ))
                {
                    if( act.PHC[k] == PH_SORPTION ||
                            act.PHC[k] == PH_POLYEL )
                    {
                        DCCW = DC_ASYM_SPECIES;
                        break;
                    }
                }
                DCCW = DC_SINGLE;
                iRet++;  // error the class code
            }
            act.DCCW[j] = DCCW;
        }   // j
NEXT_PHASE:
        j = i;
    }  // k
    ErrorIf( iRet>0, "E19IPM: ConvertDCC()", "Invalid DC class code. Memory corruption?");
}

/// Get the index of volume IC ("Vv") for the volume balance constraint
long int TActivity::getXvolume()
{
 long int ii, ret = -1;
 for( ii = act.N-1; ii>=0; ii--)
 {
  if( act.ICC[ii] == IC_VOLUME )
  { ret = ii; break; }
 }
 return ret;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// Calculation of (logarithmic) stability indexes logSI for all phases
//
void TActivity::StabilityIndexes( void )
{
    long int L1k, k, j, jb = 0;
    double ln_ax_dual, gamma_primal, x_estimate, StabIndex, logSI;
    double lnFmol = log( H2O_mol_to_kg );  // may not work with mixed-solvent electrolyte
    double lnPc = 0., Xw = 1., lnXw = 0., lnFugPur=0., YFk;
    bool fRestore; char sModPT = SM_UNDEF;

    if( act.Pc > 1e-29 )
       lnPc = log( act.Pc );
    if( act.PHC[0] == PH_AQUEL && act.YFA[0] >= act.XwMinM  ) // number of moles of solvent
    {
        Xw = act.YFA[0] / act.YF[0];
        lnXw = log( Xw );
    }
    jb=0;
    for(k=0;k<act.FI;k++)
    {
       L1k = act.L1[k]; // Number of components in the phase
       YFk = act.YF[k];
       StabIndex = 0.;
       fRestore = false;
       for(j=jb; j<jb+L1k; j++)
       {  // calculation for all components in all phases
          gamma_primal = act.Gamma[j];  // primal (external) activity coefficient
if( YFk <= act.DSM )
{
    if( !act.K2 && gamma_primal != 1.0 ) // can insert because gamma is available
       fRestore = true;
    if( act.K2 && act.GamFs[j] != 1.0 && act.Gamma[j] == 1.0 )
    {
       gamma_primal = act.GamFs[j];  // taking saved gamma if the phase was removed
       fRestore = true;
    }
    if( L1k > 1 && act.sMod )
        sModPT =  act.sMod[k][SPHAS_TYP];
    if( sModPT == SM_IDEAL )
       fRestore = true; // can always insert a simple ideal solution phase
}
else fRestore = true;
          if( gamma_primal < 1e-33 || gamma_primal > 1e33 )
              gamma_primal = 1.;

          ln_ax_dual = lg_to_ln * act.Y_la[j];  // DualTh activity
          if( ln_ax_dual < -777. )
              ln_ax_dual = -777.;
          lnFugPur = act.fDQF[j];  // Pure gas fugacity or end-member DQF parameter

          switch( act.DCC[j] ) // choice of corrections for estimated mole fractions
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
  //                gamma_primal = exp( act.F0[j] );
                 break;
             case DC_PEL_CARRIER: case DC_SUR_MINAL: case DC_SUR_CARRIER: // sorbent
                  ln_ax_dual -= lnFugPur;
  //                gamma_primal = exp( act.F0[j] );
                  break;
             default:
                  break; // error in DC class code
          }
          x_estimate = exp( ln_ax_dual )/ gamma_primal;   // estimate of DC concentration
          StabIndex += x_estimate;  // Increment to stability index
          act.NMU[j] = log( x_estimate );  // may be used for something more constructive
          act.EMU[j] = x_estimate;         // stored the estimated mole fraction of phase component
       }  // for j
       logSI = log10( StabIndex );
       if( fabs( logSI ) < log10( act.DSM ) )
           logSI = 0.;
       act.Falp[k] = logSI; // NormDoubleRound( logSI, 3 );
       if( L1k > 1 && fRestore == false && YFk < act.DSM )
           act.Falp[k] = -1.; // provisional - to indicate impossibility to restore
       jb += act.L1[k];
    }  // for k
}


//--------------------- End of s_activity.cpp ---------------------------
