//-------------------------------------------------------------------
// $Id: ipm_chemical4.cpp 799 2013-03-17 12:33:51Z kulik $
//
/// \file ipm_chemical4.cpp
/// Implementation of chemistry-specific functions for kinetics & metastability
/// for the IPM convex programming Gibbs energy minimization algorithm
//
// Copyright (c) 2013  D.Kulik, S.Dmitrieva, B.Thien
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

//-  static double ICold=0.;
/// \return status code (0 if o.k., non-zero values if there were problems
///     with kinetic/metastability models)
long int
TMulti::CalculateKinMet( long int LinkMode  )
{
   long int k, j, jb, je=0, kf, kfe=0, kp, kpe=0, ka, kae=0, ks, kse=0,
            kc, kd, kce=0, kde=0, ku, kue=0, ki, kie=0, jphl=0, jlphc=0;

   SPP_SETTING *pa = paTProfil;
   char *kMod;

   for( k=0; k<pm.FI; k++ )
   { // loop on solution phases
      jb = je;
      je += pm.L1[k];
      kMod = pmp->kMod[k];
      if( kMod[0] == KM_UNDEF )
          continue;  // skip TKinMet for this phase
      kc = kce;
      kp = kpe;
      kf = kfe;
      ka = kae;
      ks = kse;
      kd = kde;
      ku = kue;
      ki = kie;

//      for( j=jb; j<je; j++ )
//cout << "LM: " << LinkMode << " k: " << k << " dul: " << pm.DUL[j] << " dll: " << pm.DLL[j] << endl;

   // Creating TKinMet instances for phases and passing data, if needed
   switch( LinkMode )
   {
     case LINK_TP_MODE:  // Re-create TKinMet class instances and initialize them
     {
//           for( j= jb; j<je; j++ )
//           {
//                Here cleaning, if necessary
//           }
         switch( pm.PHC[k] )
         {
           //   case PH_AQUEL:
              case PH_LIQUID: case PH_SINCOND: case PH_SINDIS: // case PH_HCARBL:
           // case PH_SIMELT: case PH_GASMIX: case PH_PLASMA: case PH_FLUID:
                KM_Create( jb, k, kc, kp, kf, ka, ks, kd, ku, ki, kMod, jphl, jlphc );
                // Correction of parameters for initial T,P
                KM_ParPT( k, kMod );
                // Reset and initialize time
                KM_InitTime( k, kMod );
                KM_ReturnFSA( k, kMod );
                break;
            default:
                break;
         }
         break;
   }

   case LINK_IN_MODE:  // Initial state calculation of rates etc.
   {
    //
         switch( pm.PHC[k] )
         {
         //   case PH_AQUEL:
             case PH_LIQUID: case PH_SINCOND: case PH_SINDIS: // case PH_HCARBL:
         // case PH_SIMELT: case PH_GASMIX: case PH_PLASMA: case PH_FLUID:
             // Correction for T,P
                KM_ParPT( k, kMod );
                KM_InitTime( k, kMod );
                KM_UpdateFSA( jb, k, kMod );
                KM_InitRates( k, kMod );
                KM_SetAMRs( jb, k, kMod );
                if( k < pm.FIs )
                {
                    KM_InitUptake( jb, k, kMod );
                    KM_InitSplit( jb, k, kMod );
                }
                KM_ReturnFSA( k, kMod );
                break;
            default:
                break;
        }
        break;
    }
   case LINK_PP_MODE:  // Calculation of kinetics and metast. constraints at time step
   {
        switch( pm.PHC[k] )
        {
        //   case PH_AQUEL:
            case PH_LIQUID: case PH_SINCOND: case PH_SINDIS: // case PH_HCARBL:
        // case PH_SIMELT: case PH_GASMIX: case PH_PLASMA: case PH_FLUID:
        // Correction for T,P
                KM_ParPT( k, kMod );
                KM_UpdateTime( k, kMod );
                KM_UpdateFSA( jb, k, kMod );
                KM_CalcRates( k, kMod );
                KM_SetAMRs( jb, k, kMod );
                if( k < pm.FIs )
                {
                    KM_CalcUptake( jb, k, kMod );
                    KM_CalcSplit( jb, k, kMod );
                }
                KM_ReturnFSA( k, kMod );
                break;
            default:
                break;
        }
        break;
     }
     default: ;
   }

   jphl  += pm.LsPhl[k*2];
   jlphc += pm.LsPhl[k*2]*pm.LsPhl[k*2+1];
   kfe += pm.LsKin[k*6];
   kpe += pm.LsKin[k*6];
   kce += pm.LsKin[k*6]*pm.LsKin[k*6+2];
   kae += pm.LsKin[k*6]*pm.LsKin[k*6+1]*pm.LsKin[k*6+3];
   kse += pm.LsKin[k*6+4];  // bugfix 14.12.13 was .LsMod
   kde += pm.LsKin[k*6+1];
   if( k < pm.FIs)  
       kue += pm.LsUpt[k*2]*pmp->L1[k];  // bugfix 10.06.13
   if( kMod[0] == KM_PRO_UPT )
       kie += pm.N;
 } // k
   return 0;
}


//--------------------------------------------------------------------------------
/// Wrapper functions for creating kinetics and metastability models for phases
/// using the TKinMet class.
//
void TMulti::KM_Create( long int jb, long int k, long int kc, long int kp,
                           long int kf, long int ka, long int ks, long int kd, long int ku, long int ki,
                           const char *kmod, long int jphl, long int jlphc )
{
    double *aZ, *aM;
    KinMetData kmd;
    char KinProCode;

    kmd.KinProCod_   = kmod[0];  /// Code of the kinetic process (derived TKinMet class), see enum kinmet_controls
    kmd.KinModCod_   = kmod[1];  /// Type code of the kinetic/metastability model, see enum kinmet_controls
    kmd.KinSorpCod_  = kmod[2];  /// Type code of sorption kinetics model (solution/sorption phases only), see enum kinmet_controls
    kmd.KinLinkCod_  = kmod[3];  /// Type code of metastability links of this phase to other phases, see enum kinmet_controls
    kmd.KinSizedCod_ = kmod[4];  /// Type of particle/pore size distribution and A_s correction, see enum kinmet_controls
    kmd.KinResCod_   = kmod[5];  /// Reserved model control code

    memcpy(kmd.PhasNam_, pm.SF[k]+MAXSYMB_, MAXPHNAME_); /// Phase name (for specific built-in models)
    kmd.PhasNam_[MAXPHNAME]='\0';

    kmd.NComp_ = pm.L1[k];      /// Number of components in the phase (nDC in m_phase.h, L1[k] in MULTI)
    kmd.nlPh_ = pm.LsPhl[k*2];  /// Number of linked phases (cf. lPh), default 0
    kmd.nlPc_ = pm.LsPhl[k*2+1]; /// TKinMet, TSorpMod: number of parameters per linked phase, default 0.

    kmd.nPRk_ = pm.LsKin[k*6];    /// number of «parallel reactions» that affect amount constraints for k-th phase
    kmd.nSkr_ = pm.LsKin[k*6+1];  /// number of (aqueous or gaseous or surface) species from other reacting phases involved
    kmd.nrpC_ = pm.LsKin[k*6+2];  /// number of parameter (coefficients) involved in 'parallel reaction' terms (0 or 12)
    kmd.naptC_ = pm.LsKin[k*6+3]; /// number of parameter (coefficients) per species involved in 'activity product' terms (0 or 1)
    kmd.nAscC_ = pm.LsKin[k*6+4]; /// number of parameter coefficients in specific surface area correction equation ( 0 to 5 )
    kmd.nFaceC_ = pm.LsKin[k*6+5]; /// number of (separately considered) crystal faces or surface patches ( 1 to 4 )

    if( phKinMet[k])
        if(  phKinMet[k]->testSizes( &kmd ) )
        {
            phKinMet[k]->UpdatePT( pm.Tc, pm.Pc );
            phKinMet[k]->UpdateTime( pm.kTau, pm.kdT );
            phKinMet[k]->UpdateFSA( pm.Aalp[k], pm.XF[k], pm.FWGT[k], pm.FVOL[k], pm.Falp[k],
                                    pm.PfFact[k], pm.YOF[k], pm.IC, pm.pH, pm.pe, pm.Eh );
                return; // using old allocation and setup of the kinetics model
        }

    kmd.T_k_ = pm.Tc;     /// Temperature, K (initial)
    kmd.P_bar_ = pm.Pc;   /// Pressure, bar (initial)
    kmd.kTau_ = pm.kTau;  /// current time, s (initial)
    kmd.kdT_ = pm.kdT;    /// current time step (initial)
  //
    kmd.IS_ = pm.IC;      /// Effective molal ionic strength of aqueous electrolyte
    kmd.pH_ = pm.pH;      /// pH of aqueous solution
    kmd.pe_ = pm.pe;      /// pe of aqueous solution
    kmd.Eh_ = pm.Eh;      /// Eh of aqueous solution, V
  //
    kmd.nPh_ = pm.XF[k];   /// current amount of this phase, mol (read-only)
    kmd.mPh_ = pm.FWGT[k]/1e3; /// current mass of this phase, g to kg (read-only)
    kmd.vPh_ = pm.FVOL[k]/1e6; /// current volume of this phase, cm3 to m3 (read-only)
    kmd.sAPh_ = pm.Aalp[k]*pm.FWGT[k];  /// current surface area of this phase, m2
    kmd.LaPh_ = pm.Falp[k];    /// phase stability index (log scale)
    kmd.OmPh_ = pow( 10., pm.Falp[k] );  /// phase stability index (activity scale) 10^LaPh_
  //
    kmd.sFact_= pm.PfFact[k];  /// Shape (sphericity) factor
    kmd.sSA_ = pm.Aalp[k]*1e3; /// Specific surface area of the phase, m2/kg, default: 0.
    kmd.sgw_ = pm.Sigw[k];    /// Standard mean surface energy of solid-aqueous interface, J/m2
    kmd.sgg_ = pm.Sigg[k];    /// Standard mean surface energy of gas-aqueous interface, J/m2
    kmd.rX0_ = pm.Xr0h0[k][0]/1e9;  /// Mean radius r0 for (spherical or cylindrical) particles, n (reserved)
    kmd.hX0_ = pm.Xr0h0[k][1]/1e9;  /// Mean thickness h0 for cylindrical or 0 for spherical particles, n (reserved)
    kmd.sVp_ = 0.;    /// Specific pore volume of phase, m3/g (default: 0)
    kmd.sGP_ = pm.YOF[k]*pm.FWGT[k];    /// surface free energy of the phase, J (YOF*PhM)
    if( k < pm.FIs)
    {
        kmd.nPul_ = pm.PUL[k];   /// upper restriction to this phase amount, mol (calculated here)
        kmd.nPll_ = pm.PLL[k];   /// lower restriction to this phase amount, mol (calculated here)
    }
    else {
       kmd.nPul_ = 1e6;
       kmd.nPll_ = 0.;
    }
    kmd.arlPhc_ = pm.lPhc+jlphc;   /// Pointer to input array of phase link parameters [nlPh*nlPc]
    kmd.arfeSAr_ = pm.feSArC+kf;  /// Pointer to fractions of surface area related to different parallel reactions [nPRk]
    kmd.arrpCon_ = pm.rpConC+kc;  /// Pointer to input array of kinetic rate constants for 'parallel reactions' [nPRk*nrpC]
    kmd.arapCon_ = pm.apConC+ka;  /// Pointer to array of parameters per species involved in 'activity product' terms [nPRk * nSkr*naptC] read-only
    kmd.arAscp_ = pm.AscpC+ks;    /// Pointer to array of parameter coefficients of equation for correction of A_s [nAscC]

    kmd.SM_ = pm.SM+jb;           /// pointer to list of DC names involved in the phase [NComp]
    kmd.arDCC_ = pm.DCC+jb;       /// pointer to the classifier of DCs involved in the phase [NComp]
    kmd.arPhXC_ = pm.PhLin+jphl;  /// TSolMod, TKinMet: Phase linkage type codes [nlPh] { TBA  }

    kmd.arocPRk_ = pm.ocPRkC+kp;  /// operation codes for kinetic parallel reaction affinity terms [nPRk]
    kmd.arxSKr_= pm.xSKrC+kd;     /// pointer to input array of DC indexes used in activity products [nSKr_]

    kmd.arym_ = pm.Y_m;    /// Pointer to molalities of all species in MULTI (provided),
    kmd.arla_ = pm.Y_la;   /// Pointer to lg activities of all species in MULTI (provided)

    kmd.arxp_ = pm.XF;    /// Pointer to amounts of all phases in MULTI (provided)
    kmd.armp_ = pm.FWGT;  /// Pointer to masses of all phases in MULTI (provided)
    kmd.arvp_ = pm.FVOL;  /// Pointer to volumes of all phases in MULTI (provided)
    kmd.arasp_ = pm.Aalp;  /// Pointer to (current) specific surface areas of all phases in MULTI (provided)

    kmd.arnx_ = pm.X + jb; /// Pointer to mole amounts of phase components (provided) [NComp]

    kmd.arnxul_ = pm.DUL+jb;    /// Vector of upper kinetic restrictions to nx, moles [NComp]  (DUL) direct access output
    kmd.arnxll_ = pm.DLL+jb;    /// Vector of lower kinetic restrictions to nx, moles [NComp]  (DLL) direct access output

    kmd.arWx_ = pm.Wx + jb;     /// Species (end member) mole fractions ->NComp
    kmd.arVol_ = pm.Vol + jb;   /// molar volumes of end-members (species) cm3/mol ->NSpecies

// More stuff here, if needed


    KinProCode = kmd.KinProCod_;
    TKinMet* myKM = NULL;

   // creating instances of derived classes from the TKinMet base class
    switch( KinProCode )
    {
        case KM_PRO_MWR_:  // Kinetics of generic dissolution/precipitation (no uptake, ionex, adsorption)
        {
                TMWReaKin* myPT = new TMWReaKin( &kmd );
//                myPT->GetPhaseName( pm.SF[k] );
                myKM = (TKinMet*)myPT;
                break;
        }
        case KM_PRO_UPT_:  // Kinetics of uptake/entrapment (of minor/trace element) into solid solution
        {
             if( k < pm.FIs)
             {
                TUptakeKin* myPT = new TUptakeKin( &kmd, pm.LsUpt[k*2], pm.N, pm.UMpcC+ku, pm.xICuC+ki,
                        pm.IC_m, pm.emRd+jb, pm.emDf+jb );
                myKM = (TKinMet*)myPT;
             }
             break;
        }
        case KM_PRO_IEX_:  // Kinetics of ion exchange (clays, C-S-H, zeolites, ...)
        {
                TIonExKin* myPT = new TIonExKin( &kmd );
                myKM = (TKinMet*)myPT;
                break;
        }
        case KM_PRO_ADS_:  // Kinetics of adsorption (on MWI), redox
        {
                TAdsorpKin* myPT = new TAdsorpKin( &kmd );
                myKM = (TKinMet*)myPT;
                break;
        }
        case KM_PRO_NUPR_:  // Kinetics of nucleation followed by precipitation
        {
           // new:new: array of nucleation model parameters here (A.Testino?)
                TNucleKin* myPT = new TNucleKin( &kmd );
                myKM = (TKinMet*)myPT;
                break;
        }

        // case KM_USERDEF:
        default:
            myKM = NULL;
        	break;
    }
    if(phKinMet[k])
            delete phKinMet[k];
    phKinMet[k] = myKM; // set up new pointer for the kinetics model
    return;
}

/// Wrapper call for calculation of temperature and pressure correction
/// uses TKinMet class
void
TMulti::KM_ParPT( long int k, const char* kMod )
{
    //
    switch( kMod[0] )
    {
        case KM_PRO_MWR_:
        {
            ErrorIf( !phKinMet[k], "KinMetParPT: ","Invalid index of phase");
            TMWReaKin* myKM = (TMWReaKin*)phKinMet[k];
             myKM->PTparam( pm.Tc, pm.Pc );
             break;
        }
        case KM_PRO_UPT_:
        {
            ErrorIf( !phKinMet[k], "KinMetParPT: ","Invalid index of phase");
            TUptakeKin* myKM = (TUptakeKin*)phKinMet[k];
            myKM->PTparam( pm.Tc, pm.Pc );
            myKM->UptKinPTparam( pm.Tc, pm.Pc );
            break;
        }
        case KM_PRO_IEX_: case KM_PRO_ADS_: case KM_PRO_NUPR_:
        default:
              break;
    }
}

/// Wrapper call for initialization of time (step) variables
/// uses TKinMet class
void
TMulti::KM_InitTime( long int k, const char *kMod )
{
    //
    switch( kMod[0] )
    {
        case KM_PRO_MWR_:
        {
            ErrorIf( !phKinMet[k], "KinMetInitTime: ","Invalid index of phase");
            TMWReaKin* myKM = (TMWReaKin*)phKinMet[k];
            myKM->UpdateTime( 0., pm.kdT );
//              myKM->RateInit();
             break;
        }
        case KM_PRO_UPT_:
        {
            ErrorIf( !phKinMet[k], "KinMetInitTime: ","Invalid index of phase");
            TUptakeKin* myKM = (TUptakeKin*)phKinMet[k];
            myKM->UpdateTime( 0., pm.kdT );
            break;
        }
        case KM_PRO_IEX_: case KM_PRO_ADS_: case KM_PRO_NUPR_:
        default:
              break;
    }
}

/// Wrapper call for updating the time (step) variables
/// uses TKinMet class
void
TMulti::KM_UpdateTime( long int k, const char *kMod )
{
    //
    switch( kMod[0] )
    {
        case KM_PRO_MWR_:
        {
             ErrorIf( !phKinMet[k], "KinMetUpdateTime: ","Invalid index of phase");
             TMWReaKin* myKM = (TMWReaKin*)phKinMet[k];
             myKM->UpdateTime( pm.kTau, pm.kdT );
             break;
        }
        case KM_PRO_UPT_:
        {
            ErrorIf( !phKinMet[k], "KinMetUpdateTime: ","Invalid index of phase");
            TUptakeKin* myKM = (TUptakeKin*)phKinMet[k];
            myKM->UpdateTime( pm.kTau, pm.kdT );
            break;
        }
        case KM_PRO_IEX_: case KM_PRO_ADS_: case KM_PRO_NUPR_:
        default:
              break;
    }
}

/// Wrapper call for updating surface area fractions for parallel reactions
/// current properties of the phase (surf.area, amount, mass, volume, log stability index)
///    and current properties of aqueous solution
/// uses TKinMet class
void
TMulti::KM_UpdateFSA( long int jb, long int k, const char *kMod )
{
    double PUL=1e6, PLL=0.;
    if( k < pm.FIs )
    {
        PUL = pm.PUL[k];
        PLL = pm.PLL[k];
    }
    else {
        PUL = pm.DUL[jb];
        PLL = pm.DLL[jb];    // bugfix 9.10.2013 DK
    }
    //
    switch( kMod[0] )
    {
        case KM_PRO_MWR_:
        {
            ErrorIf( !phKinMet[k], "KinMetUpdateFSA: ","Invalid index of phase");
            TMWReaKin* myKM = (TMWReaKin*)phKinMet[k];
            myKM->UpdateFSA( pm.Aalp[k], pm.XF[k], pm.FWGT[k], pm.FVOL[k], pm.Falp[k],
                             pm.PfFact[k], pm.YOF[k], pm.IC, pm.pH, pm.pe, pm.Eh  );
            break;
        }
        case KM_PRO_UPT_:
        {
            ErrorIf( !phKinMet[k], "KinMetUpdateFSA: ","Invalid index of phase");
            TUptakeKin* myKM = (TUptakeKin*)phKinMet[k];
            myKM->UpdateFSA( pm.Aalp[k], pm.XF[k], pm.FWGT[k], pm.FVOL[k], pm.Falp[k],
                             pm.PfFact[k], pm.YOF[k], pm.IC, pm.pH, pm.pe, pm.Eh  );
            break;
        }
        case KM_PRO_IEX_: case KM_PRO_ADS_: case KM_PRO_NUPR_:
        default:
              break;
    }
}

/// Wrapper call for updating surface area of the phase and phase amount metastability constraints
/// uses TKinMet class
void
TMulti::KM_ReturnFSA( long int k, const char *kMod )
{
    double PUL=1e6, PLL=0.;
    //
    switch( kMod[0] )
    {
        case KM_PRO_MWR_:
        {
            ErrorIf( !phKinMet[k], "KinMetGetModFSA: ","Invalid index of phase");
            TMWReaKin* myKM = (TMWReaKin*)phKinMet[k];
            pm.Aalp[k] = myKM->GetModFSA( pm.PfFact[k], pm.PrT[k], pm.PkT[k], pm.PvT[k], PUL, PLL );
            break;
        }
        case KM_PRO_UPT_:
        {
            ErrorIf( !phKinMet[k], "KinMetGetModFSA: ","Invalid index of phase");
            TUptakeKin* myKM = (TUptakeKin*)phKinMet[k];
            pm.Aalp[k] = myKM->GetModFSA( pm.PfFact[k],  pm.PrT[k], pm.PkT[k], pm.PvT[k], PUL, PLL );
            break;
        }
        case KM_PRO_IEX_: case KM_PRO_ADS_: case KM_PRO_NUPR_:
        default:
              break;
    }
    if( k < pm.FIs )
    {
        pm.PUL[k] = PUL;
        pm.PLL[k] = PLL;
    }
}

// Calculation of initial kinetic rates
//
void
TMulti::KM_InitRates( long int k, const char *kMod )
{
    //
    switch( kMod[0] )
    {
        case KM_PRO_MWR_:
        {
            ErrorIf( !phKinMet[k], "KinMetInitRates: ","Invalid index of phase");
            TMWReaKin* myKM = (TMWReaKin*)phKinMet[k];
             myKM->RateInit( );
             break;
        }
        case KM_PRO_UPT_:
        {
            ErrorIf( !phKinMet[k], "KinMetInitRates: ","Invalid index of phase");
            TUptakeKin* myKM = (TUptakeKin*)phKinMet[k];
            myKM->RateInit( );
//            myKM->UptakeInit( );
            break;
        }
        case KM_PRO_IEX_: case KM_PRO_ADS_: case KM_PRO_NUPR_:
        default:
              break;
    }
}


// Calculation of current kinetic rates
//
void
TMulti::KM_CalcRates( long int k, const char *kMod )
{
    //
    switch( kMod[0] )
    {
        case KM_PRO_MWR_:
        {
            ErrorIf( !phKinMet[k], "KinMetCalcRates: ","Invalid index of phase");
            TMWReaKin* myKM = (TMWReaKin*)phKinMet[k];
            myKM->RateMod( );
            break;
        }
        case KM_PRO_UPT_:
        {
            ErrorIf( !phKinMet[k], "KinMetCalcRates: ","Invalid index of phase");
            TUptakeKin* myKM = (TUptakeKin*)phKinMet[k];
            myKM->RateMod( );
//            myKM->UptakeMod( );
            break;
        }
        case KM_PRO_IEX_: case KM_PRO_ADS_: case KM_PRO_NUPR_:
        default:
              break;
    }
}

// Calculation of initial AMR splitting for end members of SS phase
//
void
TMulti::KM_InitSplit( long int jb, long int k, const char *kMod )
{
    //
    switch( kMod[0] )
    {
        case KM_PRO_MWR_:
        {
            ErrorIf( !phKinMet[k], "KinMetInitSplit: ","Invalid index of phase");
            TMWReaKin* myKM = (TMWReaKin*)phKinMet[k];
             myKM->SplitInit( );
             break;
        }
        case KM_PRO_UPT_:
        {
            ErrorIf( !phKinMet[k], "KinMetInitSplit: ","Invalid index of phase");
            TUptakeKin* myKM = (TUptakeKin*)phKinMet[k];
            myKM->SplitInit( );
            break;
        }
        case KM_PRO_IEX_: case KM_PRO_ADS_: case KM_PRO_NUPR_:
        default:
              break;
    }
}

// Calculation of current AMR splitting for end members of SS phase
//
void
TMulti::KM_CalcSplit( long int jb, long int k, const char *kMod )
{
    //
    switch( kMod[0] )
    {
        case KM_PRO_MWR_:
        {
            ErrorIf( !phKinMet[k], "KinMetCalcSplit: ","Invalid index of phase");
            TMWReaKin* myKM = (TMWReaKin*)phKinMet[k];
             myKM->SplitMod( );
             break;
        }
        case KM_PRO_UPT_:
        {
            ErrorIf( !phKinMet[k], "KinMetCalcSplit: ","Invalid index of phase");
            TUptakeKin* myKM = (TUptakeKin*)phKinMet[k];
            myKM->SplitMod( );
            break;
        }
        case KM_PRO_IEX_: case KM_PRO_ADS_: case KM_PRO_NUPR_:
        default:
              break;
    }
}

// Sets new metastability constraints based on updated kinetic rates
//
void
TMulti::KM_SetAMRs( long int jb, long int k,const char *kMod )
{
    //
    switch( kMod[0] )
    {
        case KM_PRO_MWR_:
        {
            ErrorIf( !phKinMet[k], "KinMetCalcSplit: ","Invalid index of phase");
            TMWReaKin* myKM = (TMWReaKin*)phKinMet[k];
             myKM->SetMetCon( );
             break;
        }
        case KM_PRO_UPT_:
        {
            ErrorIf( !phKinMet[k], "KinMetCalcSplit: ","Invalid index of phase");
            TUptakeKin* myKM = (TUptakeKin*)phKinMet[k];
             myKM->SetMetCon( );
             break;
        }
        case KM_PRO_IEX_: case KM_PRO_ADS_: case KM_PRO_NUPR_:
        default:
              break;
    }
}

void
TMulti::KM_CalcUptake( long int jb, long int k, const char *kMod )
{
    //
    switch( kMod[0] )
    {
        case KM_PRO_MWR_:
        {
            ErrorIf( !phKinMet[k], "KinMetCalcUptake: ","Invalid index of phase");
            TMWReaKin* myKM = (TMWReaKin*)phKinMet[k];
            myKM->SSReaKinMod( );
            break;
        }
        case KM_PRO_UPT_:
        {
            ErrorIf( !phKinMet[k], "KinMetCalcUptake: ","Invalid index of phase");
            TUptakeKin* myKM = (TUptakeKin*)phKinMet[k];
            myKM->UptakeMod( );
            break;
        }
        case KM_PRO_IEX_: case KM_PRO_ADS_: case KM_PRO_NUPR_:
        default:
              break;
    }
}


void
TMulti::KM_InitUptake( long int jb, long int k, const char *kMod )
{   
    switch( kMod[0] )
    {      
        case KM_PRO_MWR_:
        {
            ErrorIf( !phKinMet[k], "KinMetCalcUptake: ","Invalid index of phase");
            TMWReaKin* myKM = (TMWReaKin*)phKinMet[k];
            myKM->SSReaKinInit( );
            break;
        }
        case KM_PRO_UPT_:
        {
            ErrorIf( !phKinMet[k], "KinMetInitUptake: ","Invalid index of phase");
            TUptakeKin* myKM = (TUptakeKin*)phKinMet[k];
            myKM->UptakeInit( );
            break;
        }
        case KM_PRO_IEX_: case KM_PRO_ADS_: case KM_PRO_NUPR_:
        default:
              break;
    }
}


// ------------------------------------------------------------------------------------
/// Internal memory allocation for TKinMet performance optimization
/// (since version 3.3.0)
void TMulti::Alloc_TKinMet( long int newFI )
{
  if(  phKinMet && ( newFI == sizeFI) )
    return;

  Free_TKinMet();
  // alloc memory for all multicomponents phases
  phKinMet = new  TKinMet *[newFI];
  sizeFI = newFI;
 for( long int ii=0; ii<newFI; ii++ )
          phKinMet[ii] = 0;
}

void TMulti::Free_TKinMet()
{
  long int kk;

  if( phKinMet )
  {  for(  kk=0; kk<sizeFI; kk++ )
      if( phKinMet[kk] )
           delete phKinMet[kk];

      delete[]  phKinMet;
  }
  phKinMet = 0;
  sizeFI = 0;
}


//--------------------- End of ipm_chemical4.cpp ---------------------------


