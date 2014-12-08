//-------------------------------------------------------------------
// $Id: s_kinmet.cpp 903 2013-11-14 12:12:55Z kulik $
//
// Stub implementation of TKinMet class for further development
//
// Copyright (C) 2013  D.Kulik, B.Thien, S.Dmitrieva
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
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
using namespace std;
#include "verror.h"
#include "s_kinmet.h"

//=============================================================================================
// TKinMet base class constructor (simulation of phase metastability and MWR kinetics)
// (c) DK March 2013
//
//=============================================================================================
/// generic constructor
//
TKinMet::TKinMet( const KinMetData *kmd ):
    KinProCode(kmd->KinProCod_),    KinModCode(kmd->KinModCod_),   KinSorpCode(kmd->KinSorpCod_),
    KinLinkCode(kmd->KinLinkCod_),  KinSizedCode(kmd->KinSizedCod_),  KinResCode(kmd->KinResCod_),
    NComp(kmd->NComp_), nlPh(kmd->nlPh_), nlPc(kmd->nlPc_), nPRk(kmd->nPRk_), nSkr(kmd->nSkr_),
    nrpC(kmd->nrpC_), naptC(kmd->naptC_), nAscC(kmd->nAscC_), // numpC(kmd->numpC_), iRes4(kmd->iRes4_),
    R_CONST(8.31451), T_k(kmd->T_k_), P_bar(kmd->P_bar_), kTau(kmd->kTau_), kdT(kmd->kdT_),
    sFact(kmd->sFact_), sSAi(kmd->sSA_), nPhi(kmd->nPh_), mPhi(kmd->mPh_), vPhi(kmd->vPh_),
    IS(kmd->IS_), pH(kmd->pH_),  pe(kmd->pe_),  Eh(kmd->Eh_),
    sAPh(kmd->sAPh_), LaPh(kmd->LaPh_), OmPh(kmd->OmPh_),
    sgw(kmd->sgw_),  sgg(kmd->sgg_),  rX0(kmd->rX0_),  hX0(kmd->hX0_),
    sVp(kmd->sVp_), sGP(kmd->sGP_), nPul(kmd->nPul_), nPll(kmd->nPll_)
{
    memcpy( PhasName, kmd->PhasNam_, MAXPHNAME_+1);
    // pointer assignments
    arfeSAr = kmd->arfeSAr_;   // [nPRk]
    arAscp = kmd->arAscp_;     // [nAscC]
    SM = kmd->SM_;             // [NComp]
    arDCC = kmd->arDCC_;       // [NComp]
    arPhXC = kmd->arPhXC_;     // [nlPh]
    arocPRk = kmd->arocPRk_;   // [nPRk]
    arxSKr = kmd->arxSKr_;     // [nSKr]
    arym = kmd->arym_;         // L
    arla = kmd->arla_;         // L

    arxp = kmd->arxp_;     // FI
    armp = kmd->armp_;     // FI
    arvp = kmd->arvp_;     // FI
    arasp = kmd->arasp_;   // FI

    arnx = kmd->arnx_;         // [NComp]
    arnxul = kmd->arnxul_;     // [NComp]  (DUL)
    arnxll = kmd->arnxll_;     // [NComp]  (DLL)
    arWx = kmd->arWx_;         // NComp
    arVol = kmd->arVol_;       // NComp

// 2D and 3D arrays
    arlPhc = NULL;
    arrpCon = NULL;
    arapCon = NULL;
    spcfu = NULL;
    spcfl = NULL;
//    arUmpCon = NULL;
    alloc_kinrtabs( );
    init_kinrtabs( kmd->arlPhc_, kmd->arrpCon_,  kmd->arapCon_ );

    /// allocation of work array of parameters and results for 'parallel reactions'
    arPRt = NULL;
    if(nPRk > 0 )
       arPRt = new TKinReact[nPRk];
//        alloc_arPRt();
    init_arPRt( );   // load data for parallel reactions

    // Work data and kinetic law calculation results
    pi = 3.14159265358979;
    kTot = 0.;   // Total rate constant (mol per m2 per s)
    rTot = 0.;   // Current total MWR rate (mol/s)
    gTot = 0;    // Total rate constant (kg per m2 per s)
    vTot = 0.;   // Total surface propagation velocity (m/s)

    nPh = nPhi;
    mPh = mPhi;
    vPh = vPhi;
    sSA = sSAi;
    sSAcor = sSAi;  // Initialized corrected specific surface area (m2/kg)
    sAph_c = sAPh = mPh*sSA;  // Initialized corrected surface area of the phase (m2/kg)
    sSAV = sSAVcor = sAPh/vPh;
    // Calculation of total properties of linked phases
    if( nlPh )
       linked_phases_properties( true );
    else {
        sSAlp = 0.;    sSAVlp = 0.;  sAPhlp = 0.;   // zero totals for linked phases
        mPhlp = 0.;    vPhlp = 0.;   Rholp =0.;  nPhlp = 0.;
    }
    // Calculates mean equivalent properties of particles or pores
    particles_pores_properties( true );

    kdT_c = kdT;    // Initialized corrected time step (s)
    T_k = 0.; // To trigger P-T recalculation after constructing the class instance
    OmgTol = 1e-6;  // default tolerance for checking dissolution or precipitation cases
}

/// Destructor
TKinMet::~TKinMet()
{
  free_kinrtabs();
  if( arPRt )
      delete[] arPRt;
}

// Checks dimensions in order to re-allocate the class instance, if necessary
bool
TKinMet::testSizes( const KinMetData *kmd )
{
    bool status;
    status = (KinProCode == kmd->KinProCod_) && (KinModCode == kmd->KinModCod_) && (KinSorpCode == kmd->KinSorpCod_)
     && (KinLinkCode == kmd->KinLinkCod_) &&  (KinSizedCode == kmd->KinSizedCod_) && (KinResCode == kmd->KinResCod_)
     && (NComp == kmd->NComp_) && (nlPh == kmd->nlPh_) && (nlPc == kmd->nlPc_) && (nPRk == kmd->nPRk_)
               && (nSkr == kmd->nSkr_);
    return status;
}

/// allocates memory for TKinMet data
void
TKinMet::alloc_kinrtabs()
{
   long int j, s, lp, pr;

   if( nlPh && nlPc )
   {
     arlPhc = new double *[nlPh];
     for( lp=0; lp<nlPh; lp++)
     {
          arlPhc[lp] = new double[nlPc];
     }
   }

   if( nPRk && nrpC )
   {
     arrpCon = new double *[nPRk];
     for( pr=0; pr<nPRk; pr++)
     {
          arrpCon[pr] = new double[nrpC];
     }
   }

   if( nPRk && nSkr && naptC )
   {
        arapCon = new double **[nPRk];
        for(j=0; j<nPRk; j++)
        {
            arapCon[j]   = new double *[nSkr];
            for(s=0; s<nSkr; s++)
            {
                arapCon[j][s] = new double [naptC];
            }
        }
    }
    spcfu = new double[NComp];
    spcfl = new double[NComp];
}

/// returns 0 if o.k. or some arrays were not allocated.
/// returns -1 if error was encountered.
//
long int
TKinMet::init_kinrtabs( double *p_arlPhc, double *p_arrpCon,  double *p_arapCon  )
{
    long int j, i, s, lp, pr;

    if( nlPh && nlPc && arlPhc )
    {
        for( lp=0; lp<nlPh; lp++)
            for( i=0; i<nlPc; i++)
                arlPhc[lp][i] = p_arlPhc[nlPc*lp+i];
    }
    if( arrpCon )
    {
        for( pr=0; pr<nPRk; pr++)
            for( i=0; i<nrpC; i++)
                arrpCon[pr][i] = p_arrpCon[nrpC*pr+i];
    }
    if( nPRk && nSkr && naptC )
    {
        for( j=0; j<nPRk; j++)
            for( s=0; s<nSkr; s++)
                for( i=0; i<naptC; i++)
                    arapCon[j][s][i] = p_arapCon[naptC*nSkr*j+naptC*s+i];
    }
    for( j=0; j<NComp; j++ )
    {
        spcfu[j] = arWx[j];
        spcfl[j] = arWx[j];
    }
    return 0;
}

/// frees memory for TKinMet tables
//
void
TKinMet::free_kinrtabs()
{
    long int j, s, lp, pr;

    if( arlPhc && nlPh )
    {
      for( lp=0; lp<nlPh; lp++)
      {
           delete[]arlPhc[lp];
      }
      delete[]arlPhc;
    }
    if( arrpCon && nPRk )
    {
      for( pr=0; pr<nPRk; pr++)
      {
           delete[]arrpCon[pr];
      }
      delete[]arrpCon;
    }
    if( nPRk && nSkr && naptC )
    {
         for(j=0; j<nPRk; j++)
         {
             for(s=0; s<nSkr; s++)
             {
                 delete[]arapCon[j][s];
             }
         }
         for(j=0; j<nPRk; j++)
         {
             delete[]arapCon[j];
         }
         delete[]arapCon;
    }
    delete[]spcfu;
    delete[]spcfl;
}

// creates the TKinReact array
//void alloc_arPRt()
//{
//    if( nPRk > 0 )
//    {
//        arPRt = new TKinReact[nPRk];
//    }
//}

/// Initializes TKinReact array and loads parameters into it
void
TKinMet::init_arPRt()
{
    long int xj;

    if( nPRk > 0 )
    {
        for(xj=0; xj<nPRk; xj++)
        {
            arPRt[xj].xPR = xj;   /// index of this parallel reaction
            arPRt[xj].nSa = nSkr; // number of species involved in parallel reactions
            arPRt[xj].ocPRk[0] = arocPRk[xj][0]; /// operation code for this kinetic parallel reaction affinity term
            arPRt[xj].ocPRk[1] = arocPRk[xj][1]; /// index of particle face (surface patch)
            arPRt[xj].xSKr = arxSKr;
            arPRt[xj].feSAr = arfeSAr[xj];
            arPRt[xj].rpCon = arrpCon[xj];
            if( nSkr && naptC )
                arPRt[xj].apCon = arapCon[xj];
            else arPRt[xj].apCon = NULL;
            // work data: unpacked rpCon[nrpC]
            if( nrpC >=4 )
            {
                arPRt[xj].ko = arPRt[xj].rpCon[0];  /// rate constant at standard temperature (mol/m2/s)
                arPRt[xj].Ko = arPRt[xj].rpCon[1];  /// rate constant at standard temperature (mol/m2/s)
                arPRt[xj].Ap = arPRt[xj].rpCon[2];  /// Arrhenius parameter
                arPRt[xj].Ea = arPRt[xj].rpCon[3];  /// activation energy at st.temperature J/mol
            }
            else {
                arPRt[xj].ko = arPRt[xj].Ko = arPRt[xj].Ap = arPRt[xj].Ea = 0.0;
            }
            if( nrpC >=8 )
            {
                arPRt[xj].bI = arPRt[xj].rpCon[4];
                arPRt[xj].bpH = arPRt[xj].rpCon[5];
                arPRt[xj].bpe = arPRt[xj].rpCon[6];
                arPRt[xj].bEh = arPRt[xj].rpCon[7];
            }
            else {
               arPRt[xj].bI = arPRt[xj].bpH = arPRt[xj].bpe = arPRt[xj].bEh = 0.0;
            }
            if( nrpC >=12 )
            {
                arPRt[xj].pPR = arPRt[xj].rpCon[8];
                arPRt[xj].qPR = arPRt[xj].rpCon[9];
                arPRt[xj].mPR = arPRt[xj].rpCon[10];
                arPRt[xj].uPR = arPRt[xj].rpCon[11];
            }
            else {
                arPRt[xj].pPR = arPRt[xj].qPR = arPRt[xj].mPR = arPRt[xj].uPR = 0.0;
            }
            if( nrpC > 12 )
            {
                arPRt[xj].OmEff = arPRt[xj].rpCon[12];
                arPRt[xj].nucRes = arPRt[xj].rpCon[13];
            }
            else {
                arPRt[xj].OmEff = 1.;
            }
//            arPRt[xj].Omg = OmPh; /// Input stability index non-log (d-less)

    // Results of rate term calculation
            arPRt[xj].arf = 1.;  // Arrhenius factor (temperature correction on kappa)
            arPRt[xj].cat = 1.;  // catalytic product term (f(prod(a))
            arPRt[xj].aft = 0.;  // affinity term (f(Omega))

            arPRt[xj].k = arPRt[xj].ko;   // rate constant (involving all corrections) in mol/m2/s
            arPRt[xj].K = arPRt[xj].Ko;
            // check direction - for precipitation, arPRt[xj].kPR = arPRt[xj].kop;
            arPRt[xj].rPR = 0.;   // rate for this region (output) in mol/s
            arPRt[xj].rmol = 0.;   // rate for the whole face (output) in mol/s
//        arPRt[xj].velo,   // velocity of face growth (positive) or dissolution (negative) m/s

        } // xj
    }
}

// frees memory for TKinReact array
//void free_arPRt()
//{
//    if( arPRt )
//        delete[]arPRt;
//}


// Calculation of total properties of linked phases
// and initial settings of surface area for this phase (linked to other phases)
bool
TKinMet::linked_phases_properties( bool if_init = false )
{
    // this is a phase linked to one or more other phases!
          long int klp, k, xlc, lpcode, i;
          double lc[8], xpk, mpk, vpk, aspk;

          sSAlp = 0.;    sSAVlp = 0.;  sAPhlp = 0.;   // reset totals for linked phases
          mPhlp = 0.;    vPhlp = 0.;   Rholp =0.;  nPhlp = 0.;

          // Checking if there is a phase linkage
          for( klp=0; klp<nlPh; klp++ )
          {
              k = arPhXC[klp][0];      // index of the phase in MULTI
              lpcode = arPhXC[klp][1]; // phase linkage code
              if(nlPc)
              {
                  for( i=0; i<8; i++ )
                      lc[i] =0.;
                  for( xlc=0; xlc<nlPc; xlc++ )
                      lc[xlc] = arlPhc[klp][xlc];
              }
              // getting properties of the k-th phase
              xpk = arxp[k];  // amount, mol
              mpk = armp[k]/1e3;  // mass, g to kg
              vpk = arvp[k]/1e6;  // volume, cm3 to m3
              aspk = arasp[k]*1e3; // sp.surf.area, m2/g to m2/kg
             // Adding up corrections
              switch( lpcode )   // selecting linkage type
              {
                  default:
                  case 0:     // link to phase amount
                      break;
                  case 1:     // link to (particle/pore) surface area
                  nPhlp += xpk;
                  mPhlp += mpk;
                  vPhlp += vpk;
                  sAPhlp += mpk * aspk * lc[0]; // coeff. lc[0] is fraction of substrate surface area
                                                // lc[1] for correction of shape factor - TBD
                      break;
                  case 2:     // link to (particle/pore) volume
                      break;
                  case 3:     // link to phase mass
                      break;
              }
          }
          //  Work data for totals for linked phases
          sSAlp = sAPhlp/mPhlp;   // specific surface per unit mass in m2/kg
          Rholp = mPhlp/vPhlp;    // total density, kg/m3
          sSAVlp = sAPhlp/vPhlp;  // specific surface area per unit volume in m2
          if( if_init && lpcode == 1 )
          {    // initializing surface area of the composite phase (overgrowth)
               sAPh = sAph_c = sAPhlp;  //   give overgrowth phase (part of) the total surface of substrate
               sSAcor = sSA = sAPh/(mPh+mPhlp);
               sSAVcor = sSAV = sAPh/(vPh+vPhlp);
               nPhlpi = nPhlp;
               mPhlpi = mPhlp;
               vPhlpi = vPhlp;
          }
          sSA = sAPh/(mPh+mPhlp);    // current specific surface of composite per unit mass in m2/kg
          sSAV = sAPh/(vPh+vPhlp);   // current specific surface of composite per unit volume in m2/m3
          return false;
}

// Calculates mean equivalent properties of particles or pores
bool
TKinMet::particles_pores_properties( bool if_init = false )
{
    sFact = 1.;
    if( if_init )  // initialization only
    {
       if( sFact < 1e-6 || fabs(sFact) >= 1 )
       {
           // Setting the sphericity factor (positive for particles or negative for pores)
           sFact = 1.;
           if( nAscC > 0 )
               sFact = arAscp[0];
           sFacti = sFact;
       }
    }
    Rho = (mPh+mPhlp)/(vPh+vPhlp);      // current average density of particles, kg/m3
//    sSAV = sAPh/(vPh+vPhlp);            // current specific surface per unit volume in m2/m3
//    sSA = sAPh/(mPh+mPhlp);             // current specific surface per unit mass in m2/kg
    d32 = 6./sSAV;                      // current Sauter equivalent diameter of particles, m
    dsv = d32/sFact;                    // current surface-volume equivalent diameter, m
    Vp = pi/6.*dsv*dsv*dsv;             // current equivalent volume of one particle, m3
    Mp = Vp*Rho;                        // current average mass of one partcle, kg
    if( if_init )  // initialization only
    {
        Np = (vPh+vPhlp)/Vp;       // so far, estimated number of particles remains constant
//        sSAVcor = sSAV;
        d43 = 0.;
        dsvi = dsv;
    }
    Ap = sAPh/Np;     // current surface area of one partcle, m2

    return false;
}

//=============================================================================================
// TKinReact base class constructor (to keep data for parallel reaction regions)
// (c) DK March 2013
//
//=============================================================================================
// Generic constructor
//TKinReact::TKinReact( )
//{
//
//}

// Sets new specific surface area and other properties of the phase;
// also updates 'parallel reactions' area fractions; if there is a link to other phase(s),
// makes respective corrections using the properties of other phase(s).
// Returns false if these parameters in TKinMet instance did not change; true if they did.
//
bool
TKinMet::UpdateFSA( const double pAsk, const double pXFk, const double pFWGTk, const double pFVOLk,
                    const double pLaPh, const double p_sFact, const double pYOFk,
                    const double pICa, const double ppHa, const double ppea, const double pEha )
{
    long int i;
    bool status = false;
    if( sSA != pAsk || nPh != pXFk || mPh != pFWGTk || vPh != pFVOLk || p_sFact != sFact
        || LaPh != pLaPh || IS != pICa || pH != ppHa || pe != ppea || Eh != pEha || sGP != pYOFk*mPh )
        status = true;
    sSA = pAsk*1000.; // from m2/g to m2/kg
    nPh = pXFk;
    mPh = pFWGTk/1e3;  // from g to kg
// cout << " !!! mPh: " << mPh << endl;
    vPh = pFVOLk/1e6;  // from cm3 to m3
sFacti = p_sFact;
    LaPh = pLaPh;
    sGP = pYOFk*pFWGTk;
    IS = pICa;
    pH = ppHa;
    pe = ppea;
    Eh = pEha;

    sAPh = sSA*(mPh+mPhlp);     // current surface area of the (composite) phase, m2
    sSAV = sAPh/(vPh+vPhlp);    // current specific surface in m2/m3
    OmPh = pow( 10., LaPh );    // phase stability index (activity scale) 10^LaPh_

    if( nlPh )    // get current properties of linked phases
    {
        linked_phases_properties( false );
    }
    // get current properties of particles/pores in this phase
    particles_pores_properties( false );

    for( i = 0; i < nPRk; i++ )
    {
       if( arPRt[i].feSAr != arfeSAr[i] )
       {    // copying (externally modified) surface area fractions for parallel reactions
           status = true;
           arfeSAr[i] = arPRt[i].feSAr;
       }
    }
    // Modification of the shape factor as given in the shape factor correction function
    {
        double x=0., a=0., b=0., c=0., d=0.;
        switch( KinSizedCode )
        {
           default:
           case KM_UNDEF_: // 'N' - no modification
                sFact = sFacti;
                return status;
           case KM_SIZED_ETM_:  // 'T', Empirical f(time) polynomial  (default)
                x = kTau;
                break;
           case KM_SIZED_ESI_:  // 'S', Empirical f(lgSI) polynomial
                x = LaPh;
                break;
           case KM_SIZED_ESA_:  // 'A', Empirical f(sarea-change) polynomial
                x =  ( sSA*(mPh+mPhlp) - sSAi*(mPhi+mPhlpi) ) / ( sSAi*(mPhi+mPhlpi ));
                break;
           case KM_SIZED_EVOL_: // 'V', Empirical f(volume-change) polynomial
                x =  ( vPh+vPhlp - vPhi-vPhlpi ) / (vPhi+vPhlpi);
                break;
           case KM_SIZED_MASS_: // 'M', Empirical f(mass-change) polynomial
                x =  ( mPh+mPhlp - mPhi-mPhlpi ) / (mPhi+mPhlpi);
                break;
           case KM_SIZED_MOL_:  // 'X', Empirical f(amount-change) polynomial
                x =  ( nPh+nPhlp - nPhi-nPhlpi ) / (nPhi+nPhlpi);
                break;
           case KM_SIZED_UNI_:  // 'U', Uniform particle/pore size distribution, f(d_sv)
                x = ( dsv - dsvi ) / dsvi;
                break;
           case KM_SIZED_BIN_:  // 'B', Binodal particle/pore size distribution (reserved)
                return status;
           case KM_SIZED_FUN_:  // 'F', Empirical particle/pore size distribution function (reserved)
                return status;
        }
        if( nAscC >= 1 )
            a = arAscp[0];
        if( nAscC >= 2 )
            b = arAscp[1];
        if( nAscC >= 3 )
            c = arAscp[2];
        if( nAscC >= 4 )
            d = arAscp[3];

        sFact = a + b*x +c*x*x + d*x*x*x;

        if( sFact <= 0. )
            sFact = 0.001;
        if( sFact > 1. )
            sFact = 1.;
    }

    return status;
}

// Returns (modified) specific surface area of the phase;
//  and modifies 'parallel reactions' area fractions
//  (provisionally also modified phase amount constraints)
//
double
TKinMet::GetModFSA( double& p_sFact,  double& prTot,
                    double& pkTot, double& pvTot, double& pPULk, double& pPLLk )
{
    long int i;
    for( i = 0; i < nPRk; i++ )
    {  // gets back (to Multi) modified surface area fractions
       arPRt[i].feSAr = arfeSAr[i];
    }
    p_sFact = sFact;
    pPULk = nPul;
    pPLLk = nPll;
    pkTot = kTot;
    prTot = rTot;
    pvTot = vTot;

    return sSAcor/1e3;
}

// Updates temperature to T_K and pressure to P_BAR;
// calculates Arrhenius factors and temperature-corrected rate constants in all PR regions.
//
long int
TKinMet::UpdatePT ( const double T_K, const double P_BAR )
{
    double Arf=1., kk=0., KK=0.;
    double RT = R_CONST*T_K;
    long int i;

    T_k = T_K; P_bar = P_BAR;

    for( i=0; i<nPRk; i++ )
    {
        Arf = arPRt[i].Ap * exp(-(arPRt[i].Ea)/RT);
        if( arPRt[i].Ap != 0.0 )
            arPRt[i].arf = Arf;
        else
            Arf = 1.0;
        kk = arPRt[i].ko * Arf;
        KK = arPRt[i].Ko * Arf;
        arPRt[i].k = kk;
        arPRt[i].K = KK;
    }
    return 0;
}

// sets new time Tau and time step dTau
// returns false if neither kTau nor kdT changed; true otherwise
//
bool
TKinMet::UpdateTime( const double Tau, const double dTau )
{
    bool status = false;
    if( Tau != kTau || dTau != kdT )
        status = true;
    kTau = Tau;
    kdT = dTau;
    return status;
}

// calculates the rate per m2 for r-th (xPR-th) parallel reaction
double
TKinMet::PRrateCon( TKinReact &rk, const long int r )
{
   long int xj, j, atopc, facex;
   double atp, ajp, aj, kr, bc;

//cout << "kTau: " << kTau << " k: " << rk.k << " K: " << rk.K << " Omg: " << OmPh <<
//        " nPh: " << nPh << " mPh: " << mPh << endl;

if( rk.xPR != r )     // index of this parallel reaction
        cout << rk.xPR << " <-> " << r << " mismatch" << endl;
   atopc = rk.ocPRk[0]; // operation code for this kinetic parallel reaction affinity term
   facex = rk.ocPRk[1]; // particle face index

   // activity (catalysis) product term (f(prod(a))
   rk.cat = 1.;
   if( nSkr && naptC && rk.apCon )
   {
        for( xj=0; xj < rk.nSa; xj++ )
        {
            j = rk.xSKr[xj];
            bc = rk.apCon[xj][0];
            if( bc )
            {
                aj = pow( 10., arla[j] );  // bugfix 4.10.2013 DK
                ajp = pow( aj, bc );
                rk.cat *= ajp;             // may need extension in future
            }
        }
   }
   if( rk.pPR )
       rk.cat = pow( rk.cat, rk.pPR );
   if( rk.bI )
       rk.cat *= pow( IS, rk.bI );
   if( rk.bpH )
       rk.cat *= pow( pH, rk.bpH );
   if( rk.bpe )
       rk.cat *= pow( pe, rk.bpe );
   if( rk.bEh )
       rk.cat *= pow( Eh, rk.bEh );

   // affinity term (f(Omega))
   rk.aft = 0.; atp = 0.;
   switch( atopc )    // selecting a proper affinity term
   {
    default:
    case ATOP_CLASSIC_: // = 0,     classic TST affinity term (see .../Doc/Descriptions/Kinetics-Uptake.doc)
       rk.aft = - rk.uPR + 1. - pow( OmPh, rk.qPR );
       if( rk.mPR )
           rk.aft = pow( rk.aft, rk.mPR );
       break;
    case ATOP_CLASSIC_REV_: // = 1, classic TST affinity term, reversed
       rk.aft = pow( OmPh, rk.qPR ) - 1. - rk.uPR;
       if( rk.mPR )
           rk.aft = pow( rk.aft, rk.mPR );
       rk.aft *= -1.;
       break;
    case ATOP_SCHOTT_: // = 2,      Schott et al. 2012 fig. 1e
       if( OmPh )
           rk.aft = exp( -rk.uPR/OmPh );
       else
           rk.aft = 0.;
       break;
    case ATOP_HELLMANN_: // = 3,    Hellmann Tisserand 2006 eq 9
       if( OmPh )
       {
          atp = pow( rk.qPR*log( OmPh ), rk.uPR );
          rk.aft = 1 - exp( -atp );
       }
       else {
          atp = 0.; rk.aft = 0.;
       }
       break;
    case ATOP_TENG1_: // = 4,       Teng et al. 2000, eq 13
       if( OmPh )
           atp = log( OmPh );
       else
           atp = 0.;
       rk.aft = rk.uPR * ( OmPh - 1. ) * atp;
       break;
    case ATOP_TENG2_: // = 5,       Teng et al. 2000, Fig 6
       rk.aft = pow( OmPh, rk.mPR );
       break;
    case ATOP_FRITZ_: // = 6        Fritz et al. 2009, eq 6 nucleation and growth
       atp = OmPh - rk.OmEff;
       rk.aft = pow( atp, rk.mPR );
//   nucRes
       break;
   }
   // Calculating rate for this region (output) in mol/m2/s
//   rk.rPR = 0.;
   if( LaPh < -OmgTol ) // dissolution  (needs more flexible check based on Fa stability criterion!
   {
       if( rk.k > 0 ) // k    net dissolution rate constant (corrected for T) in mol/m2/s
           rk.rPR = rk.k * rk.cat * rk.aft;
       else if ( rk.K != 0.0 ) // K  gross rate constant (corrected for T) in mol/m2/s
           rk.rPR = rk.K * rk.cat * rk.aft;
       rk.rPR *= rk.feSAr; // Theta: fraction of surface area of the solid related to this parallel reaction
   }
   if( LaPh > OmgTol ) {
       if( rk.k < 0 ) //   net precipitation rate constant (corrected for T) in mol/m2/s
           rk.rPR = fabs(rk.k) * rk.cat * rk.aft;
       else if ( rk.K != 0.0 ) // K  gross rate constant (corrected for T) in mol/m2/s
           rk.rPR = rk.K * rk.cat * rk.aft;
        rk.rPR *= rk.feSAr; // Theta: fraction of surface area of the solid related to this parallel reaction
   }
//   else {  // equilibrium - no rate change
//       ;
//   }
   return rk.rPR;
//   rmol,   // rate for the whole face (output) in mol//m2/s    TBD
//   velo,   // velocity of face growth (positive) or dissolution (negative) nm/s
//       Ap,   /// Arrhenius parameter
//       Ea,   /// activation energy at st.temperature J/mol
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Implementation of TKinMet class (TST kinetics)
//
bool
TKinMet::PTparam( const double TK, const double P )
{
    bool iRet = false;
    if( TK < 273. || TK > 5273. || P < 0. || P > 1e6 )
        iRet = true;  // error
    if( fabs( TK - T_k ) > 0.1 || fabs( P - P_bar ) > 1e-5 )
    {
       iRet = UpdatePT( TK, P );
    }
    return iRet;
}

// Returns a primitive correction for specific surface area, also for the linked phase.
//  Returns the sSAc value for rate calculation and sets fully-corrected sSAcor value.
//     more sophisticated functions to be added here
double
TKinMet::CorrSpecSurfArea( const double sFratio, const bool toinit = false )
{
    double sSAc = sSA, sSAVc = sSAV, r_sSAc;

    // more sophisticated functions to be called here
    if( !toinit )
    {
      if( LaPh < -OmgTol ) // dissolution  (needs a more flexible check based on Fa stability criterion!
      {
          sSAVc = sFratio * sSAV * dsv  / ( dsv - 2.* fabs(vTot) * kdT );
          sSAc = sSAVc / Rho;
      }
      else if( LaPh > OmgTol ) {  // precipitation
          sSAVc = sFratio * sSAV * dsv  / ( dsv + 2.* fabs(vTot) * kdT );
          sSAc = sSAVc / Rho;
      }
      else {  // equilibrium
          sSAVc = sSAV;
          sSAc = sSAcor;   // no change in sSAcor
      }
    }
    //        r_sSAc = (sSAc + sSAcor)/2.;  // Suggested by A.Denisov (PSI ENE) 04.06.2013
    sSAcor = sSAc;
    sAph_c = sSAcor * ( mPh+mPhlp );   // corrected surface area of the phase
//    sAph_c = sAPh;  // leave main surface area yet not corrected
    return sSAc;
}
/*
    { // This phase is not linked to other phases
        // Old SSA correction
        if( LaPh < -OmgTol ) // dissolution  (needs more flexible check based on Fa stability criterion!
        {
            sSAc = sFratio * sSAi * pow( mPh/mPhi, 1./3. );
        }
        else if( LaPh > OmgTol ) {  // precipitation
            sSAc = sFratio * sSAi * pow( mPh/mPhi, -1./3. );
        }
        else {  // equilibrium
            sSAc = sSAcor;   // no change in sSAcor
        }
//        r_sSAc = (sSAc + sSAcor)/2.;  // Suggested by A.Denisov (PSI ENE) 04.06.2013
        sSAcor = sSAc;
        sAPh = sSAcor * mPh;   // corrected surface of the phase.
        sAph_c = sAPh;
    }
*/

bool
TKinMet::RateInit( )
{   
    double RT = R_CONST*T_k, kPR, sSAcr;
    long int r;

    kTot = 0.;
    for( r=0; r<nPRk; r++ )
    {
        kPR = PRrateCon( arPRt[r], r ); // adds the rate constant (mol/m2/s) for r-th parallel reaction
        kTot += kPR;
    }
    nPhi = nPh;         // initial amount of this phase
    mPhi = mPh;         // initial mass
    vPhi = vPh;         // initial volume
    sSAi = sSA;         // initial specific surface area
    sFacti = sFact;     // initial shape factor
//    sSAcr = CorrSpecSurfArea( 1.0, true );
//    sAPh = sSAcr * mPh;   // initial surface of the phase
//    sAph_c = sAPh;
    gTot = kTot * mPh/nPh;  // initial rate in kg/m2/s
    rTot = kTot * sAPh; // overall initial rate (mol/s)
//    dnPh = -kdT * rTot; // overall initial change (moles)
//    dnPh = 0.; // overall initial change (moles)
//    nPh += dnPh;  // New amount of the phase (this operator is doubtful...)
    vTot =  kTot * vPh/nPh;  // orthogonal mean growth/dissolution velocity in m/s
    sSAcr = CorrSpecSurfArea( 1.0, true );
//   Here possibly additional corrections of dissolution/precipitation rates

// cout << " init 0  kTot: " << kTot << " vTot: " << vTot << " rTot: " << rTot << " sSAcor: " << sSAcor << " sAPh: " << sAPh << " nPhi: " << nPhi << endl;
 // Check initial rates and limit them to realistic values

    return false;
}

bool
TKinMet::RateMod( )
{
    double RT = R_CONST*T_k, sFratio = 1., sSAcr;
    long int r;

    kTot = 0.;   // overall specific rate (mol/m2/s)
    for( r=0; r<nPRk; r++ )
    {
       kTot += PRrateCon( arPRt[r], r ); // adds the specific rate (mol/m2/s) for r-th parallel reaction
    }
    // Calculation of shape factor ratio
    // Convert into area-mass shape factor
    // TBD
    if( fabs(sFact) > 0 )
        sFratio = sFacti/sFact;
//    sSAcr = CorrSpecSurfArea( sFratio, false );
//    sAPh = sSAcr * mPh;   // corrected surface of the phase.
//    sAph_c = sAPh;
    gTot = kTot * mPh/nPh;  // rate in kg/m2/s
    rTot = kTot * sAPh; // overall rate for the phase (mol/s)
    vTot = kTot * vPh/nPh;  // linear growth/dissolution velocity in m/s - see eq (2.11)
// cout << " t: " << kTau << " kTot: " << kTot << " vTot: " << vTot << " rTot: " << rTot << " sSAcor: " << sSAcor << " sAPh: " << sAPh << " nPh: " << nPh << endl;
    sSAcr = CorrSpecSurfArea( sFratio, false );
//   Here possibly additional corrections of dissolution/precipitation rates

    return false;
}

// Sets new metastability constraints for the whole phase  based on updated kinetic rates
//
bool
TKinMet::SetMetCon( )
{
   double dnPh = 0., dn_dt = 0.;

   dnPh = -kdT * rTot; // change in phase amount over dt
   dn_dt = -rTot;  // minus total rate mol/s

   // First calculate phase constraints
   if( LaPh < -OmgTol ) // dissolution  (needs more flexible check based on Fa stability criterion!
   {
       nPll = nPh + dnPh;
       if( nPul < nPll + fabs(dnPh) )
           nPul = nPll + fabs(dnPh);    // ensuring slackness
   }
   else if( LaPh > OmgTol ) {  // precipitation rate constant (corrected for T) in mol/m2/s
       nPul = nPh + dnPh;
       if( nPll > nPul - fabs(dnPh) )
           nPll = nPul - fabs(dnPh);    // ensuring slackness
    }
//    else {  // equilibrium  - needs to be checked!
//       nPul = nPh + fabs( dnPh );
//       nPll = nPh - fabs( dnPh );
//    }

    if( NComp > 1 )
       return false;  // DC constraints will be set in SplitMod() or SplitInit()

    // setting DC metastability constraints for a single-component phase
    arnxul[0] = nPul;
    arnxll[0] = nPll;
    return false;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Sets initial metastability constraints on end members of the (solid solution) phase
//
bool
TKinMet::SplitInit( )
{
    double dnPh, dnDC, nDC;
    long int j;

    if( LaPh < -OmgTol ) // dissolution
    {
        dnPh = nPll - nPh;
        for(j=0; j<NComp; j++)
        {
            nDC = arnx[j];
            dnDC = dnPh*spcfl[j];
            arnxll[j] = nDC + dnDC;
            if( arDCC[j] != DC_SOL_MINOR_ && arDCC[j] != DC_SOL_MINDEP_ )
            {    // not a minor/trace element
                if( arnxul[j] < arnxll[j] )
                    arnxul[j] = nDC + (nPul-nPh)*spcfu[j];
            }
            else {  // minor/trace element
                arnxul[j] = arnxll[j];
            }
        }
    }
    else if( LaPh > OmgTol )
    {  // precipitation
        dnPh = nPul - nPh;
        for(j=0; j<NComp; j++)
        {
            nDC = arnx[j];
            dnDC = dnPh*spcfu[j];
            arnxul[j] = nDC + dnDC;
            if( arDCC[j] != DC_SOL_MINOR_ && arDCC[j] != DC_SOL_MINDEP_ )
            {    // not a minor/trace element
                if( arnxll[j] > arnxul[j] )
                    arnxll[j] = nDC + (nPll-nPh)*spcfl[j];
            }
            else {  // minor/trace element
                arnxll[j] = arnxul[j];
            }
        }
    }
    else {  // equilibrium
        double delOmg = pow( 10., OmgTol );
        for(j=0; j<NComp; j++)
        {
            nDC = arnx[j];
            arnxul[j] = nDC + delOmg;
            arnxll[j] = nDC - delOmg;
        }
    }
    return false;
}

// Sets current metastability constraints on end members of the (solid solution) phase
//
bool
TKinMet::SplitMod( )
{
    double dnPh, dnDC, nDC;
    long int j;

    if( LaPh < -OmgTol ) // dissolution
    {
        dnPh = nPll - nPh;
        for(j=0; j<NComp; j++)
        {
            nDC = arnx[j];
            dnDC = dnPh*spcfl[j];
            arnxll[j] = nDC + dnDC;
            if( arDCC[j] != DC_SOL_MINOR_ && arDCC[j] != DC_SOL_MINDEP_ )
            {    // not a minor/trace element
                if( arnxul[j] < arnxll[j] )
//                    arnxul[j] = nPul*spcfu[j];
                    arnxul[j] = nDC + (nPul-nPh)*spcfu[j];
            }
            else {  // minor/trace element
                arnxul[j] = arnxll[j];
            }
        }
    }
    else if( LaPh > OmgTol )
    {  // precipitation
        dnPh = nPul - nPh;
        for(j=0; j<NComp; j++)
        {
            nDC = arnx[j];
            dnDC = dnPh*spcfu[j];
            arnxul[j] = nDC + dnDC;
            if( arDCC[j] != DC_SOL_MINOR_ && arDCC[j] != DC_SOL_MINDEP_ )
            {    // not a minor/trace element
                if( arnxll[j] > arnxul[j] )
//                    arnxll[j] = nPll*spcfl[j];
                    arnxll[j] = nDC + (nPll-nPh)*spcfl[j];
            }
            else {  // minor/trace element
                arnxll[j] = arnxul[j];
            }
        }
    }
    else {  // equilibrium
        double delOmg = pow( 10., OmgTol );
        for(j=0; j<NComp; j++)
        {
            nDC = arnx[j];
            arnxul[j] = nDC + delOmg;
            arnxll[j] = nDC - delOmg;
        }
    }
    return false;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Implementation of TMWReaKin class
//
// Initializes uptake rates
bool
TMWReaKin::SSReaKinInit()
{
    return false;
}

// Calculates uptake rates
bool
TMWReaKin::SSReaKinMod()
{
    return false;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Implementation of TUptakeKin class (uptake kinetics)
//
TUptakeKin::TUptakeKin(KinMetData *kmd, long int p_numpC, long int p_nElm, double *p_arUmpCon ,
            long int *p_arxICu, double *p_arElm, double *p_Rdj, double *p_Dfj ):TKinMet( kmd )
{
    numpC = p_numpC;
nElm = p_nElm;
    arUmpCon = NULL;

arxICu = p_arxICu; // Added 14.06.13 DK
arElm = p_arElm;
Rdj = p_Rdj;  // direct access
Dfj = p_Dfj;  // direct access

    alloc_upttabs();
    init_upttabs( p_arUmpCon );
};

/// Destructor
TUptakeKin::~TUptakeKin()
{
  free_upttabs();
}

/// allocates memory for TUptakeKin data
void
TUptakeKin::alloc_upttabs()
{
   long int j;

   if( NComp && numpC )
   {
     arUmpCon = new double *[NComp];
     for( j=0; j<NComp; j++)
     {
          arUmpCon[j] = new double[numpC];
     }
   }

}

/// returns 0 if o.k. or some arrays were not allocated.
/// returns -1 if error was encountered.
//
long int
TUptakeKin::init_upttabs( double *p_arUmpCon  )
{
    long int j, i;

    if( arUmpCon ) {
        for( j=0; j<NComp; j++)
            for( i=0; i<numpC; i++)
                arUmpCon[j][i] = p_arUmpCon[numpC*j+i];
    }
    return 0;
}

/// frees memory for TKinMet tables
//
void
TUptakeKin::free_upttabs()
{
    long int j;

    if( arUmpCon )
    {
      for( j=0; j<NComp; j++)
      {
           delete[]arUmpCon[j];
      }
      delete[]arUmpCon;
    }
}

bool
TUptakeKin::UptKinPTparam( const double TK, const double P )
{
    // No T,P corrections so far ...
    return false;
}

// Initializes uptake rates by setting spcfu and spcfl vectors
//  proportional to SS mole fractions
bool
TUptakeKin::UptakeInit()
{
    long int i, j;
    double molSum = 0., Rd_rest;

    for( j=0; j<NComp; j++ )
    {
        i = arxICu[j];
        molSum += arElm[i];
    }

    for( j=0; j<NComp; j++ )
    {
       i = arxICu[j];
       spcfu[j] = arWx[j];
       spcfl[j] = arWx[j];
     Rdj[j] = arWx[j]/arElm[i];
     Rd_rest = (1.-arWx[j])/(molSum-arElm[i]);
     Dfj[j] = Rdj[j]/Rd_rest;
    }
    return false;
}

// Calculates uptake rates
bool
TUptakeKin::UptakeMod()
{
    long int j, i;

    switch( KinSorpCode )
    {
        case  KM_UPT_ENTRAP_: //  = 'E',  //	Unified entrapment model (Thien,Kulik,Curti 2013)
        {
            double FTr, DelTr0, Ds, Dl, l, m, xtTr, xtHc, CF, Rd_rest;
            double DelTr, Vml, molSum=0., molMinSum=0., molMajSum=0., spMinSum=0, spMajSum=0;

// Calculating the sums of tot.diss.molal. for all elements relevant to major and minor endmembers
            for( j=0; j<NComp; j++ )
            {
                i = arxICu[j];
                if( arDCC[j] != DC_SOL_MINOR_ && arDCC[j] != DC_SOL_MINDEP_ )
                {    // not a minor/trace element
                    molMajSum += arElm[i];
                }
                else {
                    molMinSum += arElm[i];
                }
            }
            molSum = molMajSum + molMinSum;
//            if(spMajSum < 1e-9)
//               Error - no host components left!
// Calculating the fractionation and splitting coefficients for Tr end members
            for( j=0; j<NComp; j++ )
            {
                 i = arxICu[j];
                if( arDCC[j] != DC_SOL_MINOR_ && arDCC[j] != DC_SOL_MINDEP_ )
                {
                    spMajSum += arWx[j];
                    continue; // not a minor/trace element
                }
                // Minor/trace component
                FTr =    arUmpCon[j][0];  // d-less
                DelTr0 = arUmpCon[j][1];  // eq tr fract.coeff.
                Ds =     arUmpCon[j][2];  // in nm2/s
                Dl =     arUmpCon[j][3];  // in nm2/s
                l =      arUmpCon[j][4];  // in nm
                m =      arUmpCon[j][5];  // d-less

                // Calculate orthogonal linear rate
                Vml = -vTot * m * l * 1e9;     // here vTot is in m/s and Vml is in nm2/s
                DelTr = DelTr0 * ( Ds + Vml ) / ( Ds + Vml/FTr );    // Effective fractionation coeff.
                xtTr = DelTr * arElm[i]/molMajSum;     // Frac.coeff. defined rel to sum of major EMs!
                spcfu[j] = xtTr;
                spcfl[j] = spcfu[j];
                spMinSum += xtTr;
            }
// Correcting splitting coeffs of major EMs for changed sum of split.coeffs. for Tr EMs
            spMajSum = 1.- spMinSum;    // bugfix  14.11.2013 DK
            CF = (1.-spMinSum)/spMajSum;
            for( j=0; j<NComp; j++ )
            {
                i = arxICu[j];
                if( arDCC[j] != DC_SOL_MINOR_ && arDCC[j] != DC_SOL_MINDEP_ )
                {    // not a minor/trace element
                    xtHc = spcfu[j]*CF; //  xtHc = arWx[j]*CF;
                    spcfu[j] = xtHc;
                    spcfl[j] = spcfu[j];   
                }
                // Rd and Df calculation (on the bulk basis)
                Rdj[j] = spcfu[j]/arElm[i];
                Rd_rest = (1-spcfu[j])/(molSum-arElm[i]);
                Dfj[j] = Rdj[j]/Rd_rest;
                // Other ways of Dfj calculation are possible!
            }
            break;
        }
        case KM_UPT_UPDP_:    // = 'M',   //	DePaolo (2011) uptake kinetics model
            break;
        case KM_UPT_SEMO_:    // = 'G',   //  Growth (surface) entrapment model (Watson 2004)
            break;
        case KM_IEX_FAST_:    // = 'F',   //  Fast ion exchange kinetics (e.g. montmorillonite, CSH)
            break;
        case KM_IEX_SLOW_:    // = 'L',   //  Slow ion exchange kinetics (e.g. illite, zeolites)
            break;
        default:
            break;
    }

    return false;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


/*
// -----------------------------------------------------------------------------
// Implementation of derived class main functionality
//
bool PTparam()
{
    return 0;
};

bool RateMod()
{
    return 0;
};

bool SplitMod()
{
    return 0;
};

bool SplitInit()
{
    return 0;
};

bool UptakeMod()
{
    return 0;
};

bool UptakeInit()
{
    return 0;
};

bool SetMetCon()
{
    return 0;
};


*/


//--------------------- End of s_kinmet.cpp ---------------------------

