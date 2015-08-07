//-------------------------------------------------------------------
// $Id: node_kinetics.cpp 799 2013-03-17 12:33:51Z kulik $
//
/// \file node_kinetics.cpp
/// Implementation of chemistry-specific functions for kinetics & metastability
/// for the IPM convex programming Gibbs energy minimization algorithm
//
// Copyright (c) 2013-2014  D.Kulik, S.Dmitrieva, B.Thien
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
#include "kinetics.h"


void TKinetics::set_def( void )
{

}



// void TNode::setSpeciesUpperAMRs( const double* nu )
// {
// }

// void TNode::setSpeciesLowerAMRs( const double* nl )
// {
// }

// void TNode::setPhasesUpperAMRs( const double* nfu )
// {
// }

// void TNode::setPhasesLowerAMRs( const double* nfl )
// {
// }


// long int TNode::updateKineticsMetastability( long int LinkMode )
// {
// }


/*  Come back to this at a later stage
///  This procedure sets kinetic constraints according to a given
///  concentration units.
//  Needs much more work, elaboration, and performance optimization
//
void TKinetics::Set_DC_limits( long int Mode )
{
    double XFL, XFU, XFS=0., XFM, MWXW, MXV, XL=0., XU=0.;
    long int jb, je, j,k, MpL;
    char tbuf[150];

//    if( !this.PLIM )
//       return;  // no metastability limits to be set
// ???????????????????????????????????????
    CalculateConcentrations( kin.X, kin.XF, kin.XFA );

    for(k=0; k<kin.FI; k++)
        XFS+=kin.XF[k];  // calculate sum of moles in all phases

    jb=0;
    for( k=0; k<kin.FI; k++ )
    { // cycle over phases
        je=jb+kin.L1[k];
//        XFM=0.;
        MWXW =0.;
        MXV = 0.;
        XFL = 0.;
        XFU = 1e6;
        if( Mode && kin.XF[k] < kin.DSM )
            goto NEXT_PHASE;
        XFM = kin.FWGT[k]; // Mass of a phase
        if( Mode )
        {
            MWXW = XFM/kin.XF[k];         // current molar mass of phase
            MXV = kin.FVOL[k]/kin.XF[k]; // current molar volume of phase
        }
        // Check codes for phase DC
        MpL=0;
        for( j=jb; j<je; j++ )
            if( kin.RLC[j] != NO_LIM )
                MpL = 1;

if( k < kin.FIs )
{					// Temporary workaround - DK  13.12.2007
        if( kin.RFLC[k] == NO_LIM && !MpL )
        { // check type restrictions on phase
            goto NEXT_PHASE;
        }
        switch( kin.RFSC[k] )
        { // check scale restrictions on phase in all system
        case QUAN_MOL:
            XFL = Mode? kin.XF[k]: kin.PLL[k];
            XFU = Mode? kin.XF[k]: kin.PUL[k];
            break;
        case CON_MOLAL:
            XFL = Mode? kin.XF[k]: kin.PLL[k]*kin.GWAT/H2O_mol_to_kg;
            XFU = Mode? kin.XF[k]: kin.PUL[k]*kin.GWAT/H2O_mol_to_kg;
            break;
        case CON_MOLFR:
            XFL = Mode? kin.XF[k]: kin.PLL[k]*XFS;
            XFU = Mode? kin.XF[k]: kin.PUL[k]*XFS;
            break;
        case CON_WTFR:   if(MWXW < 1e-15) break;  // Temp.fix
            XFL = Mode? kin.XF[k]: kin.PLL[k]*kin.MBX/MWXW;
            XFU = Mode? kin.XF[k]: kin.PUL[k]*kin.MBX/MWXW;
            break;
        case CON_VOLFR:   if(MXV < 1e-15) break; // Temp.fix
            XFL = Mode? kin.XF[k]: kin.PLL[k]*kin.VXc/MXV;
            XFU = Mode? kin.XF[k]: kin.PUL[k]*kin.VXc/MXV;
            break;
        default:
            ; // do more?
        }
//        if( kin.RFLC[k] == NO_LIM )
//        {                            Temporary!
            XFL = 0.0;
            XFU = 1e6;
//        }
}
        for( j=jb; j<je; j++ )
        { // loop over DCs
            if( kin.RLC[j] == NO_LIM )
                continue;

            if( Mode )
            {
                XU = kin.DUL[j];
                XL = kin.DLL[j];
            }
            else
                switch( kin.RSC[j] ) // get initial limits
                {
                case QUAN_MOL:
                    XU = kin.DUL[j];
                    XL = kin.DLL[j];
                    break;
                case CON_MOLAL:
                    XU = kin.DUL[j]*kin.GWAT/H2O_mol_to_kg;
                    XL = kin.DLL[j]*kin.GWAT/H2O_mol_to_kg;
                    break;
                case CON_MOLFR:
                    XU = kin.DUL[j]*XFU;
                    XL = kin.DLL[j]*XFL;
                    break;
                case CON_WTFR:
//Ask DK! 20/04/2002
#ifndef IPMGEMPLUGIN
                    XU = kin.DUL[j]*XFU*MWXW /
         TProfil::pm->MolWeight(kin.N, kin.Awt, kin.A+j*kin.N );
                    XL = kin.DLL[j]*XFL*MWXW /
         TProfil::pm->MolWeight(kin.N, kin.Awt, kin.A+j*kin.N );

#endif
                    break;
                case CON_VOLFR:
                    XU = kin.DUL[j]*XFU*MXV/ kin.Vol[j];
                    XL = kin.DLL[j]*XFL*MXV/ kin.Vol[j];
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
//                XU = XFU; // - kin.lowPosNum;
            }
            if( XL < XFL )
            {
//                JJ = j;
//                KK = k;
                sprintf( tbuf, "Inconsistent lower DC metastability limits j=%ld k=%ld XL=%g XFL=%g",
                         j, k, XL, XFL );
                Error( "E12IPM: Set_DC_limits(): ",tbuf );
//                XL = XFL; // - kin.lowPosNum;
            }
            kin.DUL[j]=XU;
            kin.DLL[j]=XL;
        }   // j
NEXT_PHASE:
        jb = je;
    }  // k
}
*/

//-  static double ICold=0.;
/// \return status code (0 if o.k., non-zero values if there were problems
///     with kinetic/metastability models)
long int
TKinetics::CalculateKinMet( long int LinkMode  )
{
   long int k, j, jb, je=0, kf, kfe=0, kp, kpe=0, ka, kae=0, ks, kse=0,
            kc, kd, kce=0, kde=0, ku, kue=0, ki, kie=0, jphl=0, jlphc=0;

//   SPP_SETTING *pa = paTProfil;
   char *kMod;

   for( k=0; k< kin.FI; k++ )
   { // loop on solution phases
      jb = je;
      je += kin.L1[k];
      kMod = kin.kMod[k];
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
//cout << "LM: " << LinkMode << " k: " << k << " dul: " << kin.DUL[j] << " dll: " << kin.DLL[j] << endl;

   // Creating TKinMet instances for phases and passing data, if needed
   switch( LinkMode )
   {
     case LINK_TP_MODE:  // Re-create TKinMet class instances and initialize them
     {
//           for( j= jb; j<je; j++ )
//           {
//                Here cleaning, if necessary
//           }
         switch( kin.PHC[k] )
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
         switch( kin.PHC[k] )
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
                if( k < kin.FIs )
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
        switch( kin.PHC[k] )
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
                if( k < kin.FIs )
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

   jphl  += kin.LsPhl[k*2];
   jlphc += kin.LsPhl[k*2]*kin.LsPhl[k*2+1];
   kfe += kin.LsKin[k*6];
   kpe += kin.LsKin[k*6];
   kce += kin.LsKin[k*6]*kin.LsKin[k*6+2];
   kae += kin.LsKin[k*6]*kin.LsKin[k*6+1]*kin.LsKin[k*6+3];
   kse += kin.LsKin[k*6+4];  // bugfix 14.12.13 was .LsMod
   kde += kin.LsKin[k*6+1];
   if( k < kin.FIs)
       kue += kin.LsUpt[k*2]*kin.L1[k];  // bugfix 10.06.13
   if( kMod[0] == KM_PRO_UPT )
       kie += kin.N;
 } // k
   return 0;
}


//--------------------------------------------------------------------------------
/// Wrapper functions for creating kinetics and metastability models for phases
/// using the TKinMet class.
//
void TKinetics::KM_Create( long int jb, long int k, long int kc, long int kp,
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

    memcpy(kmd.PhasNam_, kin.SF[k]+MAXSYMB_, MAXPHNAME_); /// Phase name (for specific built-in models)
    kmd.PhasNam_[MAXPHNAME]='\0';

    kmd.NComp_ = kin.L1[k];      /// Number of components in the phase (nDC in m_phase.h, L1[k] in MULTI)
    kmd.nlPh_ = kin.LsPhl[k*2];  /// Number of linked phases (cf. lPh), default 0
    kmd.nlPc_ = kin.LsPhl[k*2+1]; /// TKinMet, TSorpMod: number of parameters per linked phase, default 0.

    kmd.nPRk_ = kin.LsKin[k*6];    /// number of «parallel reactions» that affect amount constraints for k-th phase
    kmd.nSkr_ = kin.LsKin[k*6+1];  /// number of (aqueous or gaseous or surface) species from other reacting phases involved
    kmd.nrpC_ = kin.LsKin[k*6+2];  /// number of parameter (coefficients) involved in 'parallel reaction' terms (0 or 12)
    kmd.naptC_ = kin.LsKin[k*6+3]; /// number of parameter (coefficients) per species involved in 'activity product' terms (0 or 1)
    kmd.nAscC_ = kin.LsKin[k*6+4]; /// number of parameter coefficients in specific surface area correction equation ( 0 to 5 )
    kmd.nFaceC_ = kin.LsKin[k*6+5]; /// number of (separately considered) crystal faces or surface patches ( 1 to 4 )

    if( phKinMet[k])
        if(  phKinMet[k]->testSizes( &kmd ) )
        {
            phKinMet[k]->UpdatePT( kin.Tc, kin.Pc );
            phKinMet[k]->UpdateTime( kin.kTau, kin.kdT );
            phKinMet[k]->UpdateFSA( kin.Aalp[k], kin.XF[k], kin.FWGT[k], kin.FVOL[k], kin.Falp[k],
                                    kin.PfFact[k], kin.YOF[k], kin.IC, kin.pH, kin.pe, kin.Eh );
                return; // using old allocation and setup of the kinetics model
        }

    kmd.T_k_ = kin.Tc;     /// Temperature, K (initial)
    kmd.P_bar_ = kin.Pc;   /// Pressure, bar (initial)
    kmd.kTau_ = kin.kTau;  /// current time, s (initial)
    kmd.kdT_ = kin.kdT;    /// current time step (initial)
  //
    kmd.IS_ = kin.IC;      /// Effective molal ionic strength of aqueous electrolyte
    kmd.pH_ = kin.pH;      /// pH of aqueous solution
    kmd.pe_ = kin.pe;      /// pe of aqueous solution
    kmd.Eh_ = kin.Eh;      /// Eh of aqueous solution, V
  //
    kmd.nPh_ = kin.XF[k];   /// current amount of this phase, mol (read-only)
    kmd.mPh_ = kin.FWGT[k]/1e3; /// current mass of this phase, g to kg (read-only)
    kmd.vPh_ = kin.FVOL[k]/1e6; /// current volume of this phase, cm3 to m3 (read-only)
    kmd.sAPh_ = kin.Aalp[k]*kin.FWGT[k];  /// current surface area of this phase, m2
    kmd.LaPh_ = kin.Falp[k];    /// phase stability index (log scale)
    kmd.OmPh_ = pow( 10., kin.Falp[k] );  /// phase stability index (activity scale) 10^LaPh_
  //
    kmd.sFact_= kin.PfFact[k];  /// Shape (sphericity) factor
    kmd.sSA_ = kin.Aalp[k]*1e3; /// Specific surface area of the phase, m2/kg, default: 0.
    kmd.sgw_ = kin.Sigw[k];    /// Standard mean surface energy of solid-aqueous interface, J/m2
    kmd.sgg_ = kin.Sigg[k];    /// Standard mean surface energy of gas-aqueous interface, J/m2
    kmd.rX0_ = kin.Xr0h0[k][0]/1e9;  /// Mean radius r0 for (spherical or cylindrical) particles, n (reserved)
    kmd.hX0_ = kin.Xr0h0[k][1]/1e9;  /// Mean thickness h0 for cylindrical or 0 for spherical particles, n (reserved)
    kmd.sVp_ = 0.;    /// Specific pore volume of phase, m3/g (default: 0)
    kmd.sGP_ = kin.YOF[k]*kin.FWGT[k];    /// surface free energy of the phase, J (YOF*PhM)
    if( k < kin.FIs)
    {
        kmd.nPul_ = kin.PUL[k];   /// upper restriction to this phase amount, mol (calculated here)
        kmd.nPll_ = kin.PLL[k];   /// lower restriction to this phase amount, mol (calculated here)
    }
    else {
       kmd.nPul_ = 1e6;
       kmd.nPll_ = 0.;
    }
    kmd.arlPhc_ = kin.lPhc+jlphc;   /// Pointer to input array of phase link parameters [nlPh*nlPc]
    kmd.arfeSAr_ = kin.feSArC+kf;  /// Pointer to fractions of surface area related to different parallel reactions [nPRk]
    kmd.arrpCon_ = kin.rpConC+kc;  /// Pointer to input array of kinetic rate constants for 'parallel reactions' [nPRk*nrpC]
    kmd.arapCon_ = kin.apConC+ka;  /// Pointer to array of parameters per species involved in 'activity product' terms [nPRk * nSkr*naptC] read-only
    kmd.arAscp_ = kin.AscpC+ks;    /// Pointer to array of parameter coefficients of equation for correction of A_s [nAscC]

    kmd.SM_ = kin.SM+jb;           /// pointer to list of DC names involved in the phase [NComp]
    kmd.arDCC_ = kin.DCC+jb;       /// pointer to the classifier of DCs involved in the phase [NComp]
    kmd.arPhXC_ = kin.PhLin+jphl;  /// TSolMod, TKinMet: Phase linkage type codes [nlPh] { TBA  }

    kmd.arocPRk_ = kin.ocPRkC+kp;  /// operation codes for kinetic parallel reaction affinity terms [nPRk]
    kmd.arxSKr_= kin.xSKrC+kd;     /// pointer to input array of DC indexes used in activity products [nSKr_]

    kmd.arym_ = kin.Y_m;    /// Pointer to molalities of all species in MULTI (provided),
    kmd.arla_ = kin.Y_la;   /// Pointer to lg activities of all species in MULTI (provided)

    kmd.arxp_ = kin.XF;    /// Pointer to amounts of all phases in MULTI (provided)
    kmd.armp_ = kin.FWGT;  /// Pointer to masses of all phases in MULTI (provided)
    kmd.arvp_ = kin.FVOL;  /// Pointer to volumes of all phases in MULTI (provided)
    kmd.arasp_ = kin.Aalp;  /// Pointer to (current) specific surface areas of all phases in MULTI (provided)

    kmd.arnx_ = kin.X + jb; /// Pointer to mole amounts of phase components (provided) [NComp]

    kmd.arnxul_ = kin.DUL+jb;    /// Vector of upper kinetic restrictions to nx, moles [NComp]  (DUL) direct access output
    kmd.arnxll_ = kin.DLL+jb;    /// Vector of lower kinetic restrictions to nx, moles [NComp]  (DLL) direct access output

    kmd.arWx_ = kin.Wx + jb;     /// Species (end member) mole fractions ->NComp
    kmd.arVol_ = kin.Vol + jb;   /// molar volumes of end-members (species) cm3/mol ->NSpecies

// More stuff here, if needed


    KinProCode = kmd.KinProCod_;
    TKinMet* myKM = NULL;

   // creating instances of derived classes from the TKinMet base class
    switch( KinProCode )
    {
        case KM_PRO_MWR_:  // Kinetics of generic dissolution/precipitation (no uptake, ionex, adsorption)
        {
                TMWReaKin* myPT = new TMWReaKin( &kmd );
//                myPT->GetPhaseName( kin.SF[k] );
                myKM = (TKinMet*)myPT;
                break;
        }
        case KM_PRO_UPT_:  // Kinetics of uptake/entrapment (of minor/trace element) into solid solution
        {
             if( k < kin.FIs)
             {
                TUptakeKin* myPT = new TUptakeKin( &kmd, kin.LsUpt[k*2], kin.N, kin.UMpcC+ku, kin.xICuC+ki,
                        kin.IC_m, kin.emRd+jb, kin.emDf+jb );
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
TKinetics::KM_ParPT( long int k, const char* kMod )
{
    //
    switch( kMod[0] )
    {
        case KM_PRO_MWR_:
        {
            ErrorIf( !phKinMet[k], "KinMetParPT: ","Invalid index of phase");
            TMWReaKin* myKM = (TMWReaKin*)phKinMet[k];
             myKM->PTparam( kin.Tc, kin.Pc );
             break;
        }
        case KM_PRO_UPT_:
        {
            ErrorIf( !phKinMet[k], "KinMetParPT: ","Invalid index of phase");
            TUptakeKin* myKM = (TUptakeKin*)phKinMet[k];
            myKM->PTparam( kin.Tc, kin.Pc );
            myKM->UptKinPTparam( kin.Tc, kin.Pc );
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
TKinetics::KM_InitTime( long int k, const char *kMod )
{
    //
    switch( kMod[0] )
    {
        case KM_PRO_MWR_:
        {
            ErrorIf( !phKinMet[k], "KinMetInitTime: ","Invalid index of phase");
            TMWReaKin* myKM = (TMWReaKin*)phKinMet[k];
            myKM->UpdateTime( 0., kin.kdT );
//              myKM->RateInit();
             break;
        }
        case KM_PRO_UPT_:
        {
            ErrorIf( !phKinMet[k], "KinMetInitTime: ","Invalid index of phase");
            TUptakeKin* myKM = (TUptakeKin*)phKinMet[k];
            myKM->UpdateTime( 0., kin.kdT );
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
TKinetics::KM_UpdateTime( long int k, const char *kMod )
{
    //
    switch( kMod[0] )
    {
        case KM_PRO_MWR_:
        {
             ErrorIf( !phKinMet[k], "KinMetUpdateTime: ","Invalid index of phase");
             TMWReaKin* myKM = (TMWReaKin*)phKinMet[k];
             myKM->UpdateTime( kin.kTau, kin.kdT );
             break;
        }
        case KM_PRO_UPT_:
        {
            ErrorIf( !phKinMet[k], "KinMetUpdateTime: ","Invalid index of phase");
            TUptakeKin* myKM = (TUptakeKin*)phKinMet[k];
            myKM->UpdateTime( kin.kTau, kin.kdT );
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
TKinetics::KM_UpdateFSA( long int jb, long int k, const char *kMod )
{
    double PUL=1e6, PLL=0.;
    if( k < kin.FIs )
    {
        PUL = kin.PUL[k];
        PLL = kin.PLL[k];
    }
    else {
        PUL = kin.DUL[jb];
        PLL = kin.DLL[jb];    // bugfix 9.10.2013 DK
    }
    //
    switch( kMod[0] )
    {
        case KM_PRO_MWR_:
        {
            ErrorIf( !phKinMet[k], "KinMetUpdateFSA: ","Invalid index of phase");
            TMWReaKin* myKM = (TMWReaKin*)phKinMet[k];
            myKM->UpdateFSA( kin.Aalp[k], kin.XF[k], kin.FWGT[k], kin.FVOL[k], kin.Falp[k],
                             kin.PfFact[k], kin.YOF[k], kin.IC, kin.pH, kin.pe, kin.Eh  );
            break;
        }
        case KM_PRO_UPT_:
        {
            ErrorIf( !phKinMet[k], "KinMetUpdateFSA: ","Invalid index of phase");
            TUptakeKin* myKM = (TUptakeKin*)phKinMet[k];
            myKM->UpdateFSA( kin.Aalp[k], kin.XF[k], kin.FWGT[k], kin.FVOL[k], kin.Falp[k],
                             kin.PfFact[k], kin.YOF[k], kin.IC, kin.pH, kin.pe, kin.Eh  );
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
TKinetics::KM_ReturnFSA( long int k, const char *kMod )
{
    double PUL=1e6, PLL=0.;
    //
    switch( kMod[0] )
    {
        case KM_PRO_MWR_:
        {
            ErrorIf( !phKinMet[k], "KinMetGetModFSA: ","Invalid index of phase");
            TMWReaKin* myKM = (TMWReaKin*)phKinMet[k];
            kin.Aalp[k] = myKM->GetModFSA( kin.PfFact[k], kin.PrT[k], kin.PkT[k], kin.PvT[k], PUL, PLL );
            break;
        }
        case KM_PRO_UPT_:
        {
            ErrorIf( !phKinMet[k], "KinMetGetModFSA: ","Invalid index of phase");
            TUptakeKin* myKM = (TUptakeKin*)phKinMet[k];
            kin.Aalp[k] = myKM->GetModFSA( kin.PfFact[k],  kin.PrT[k], kin.PkT[k], kin.PvT[k], PUL, PLL );
            break;
        }
        case KM_PRO_IEX_: case KM_PRO_ADS_: case KM_PRO_NUPR_:
        default:
              break;
    }
    if( k < kin.FIs )
    {
        kin.PUL[k] = PUL;
        kin.PLL[k] = PLL;
    }
}

// Calculation of initial kinetic rates
//
void
TKinetics::KM_InitRates( long int k, const char *kMod )
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
TKinetics::KM_CalcRates( long int k, const char *kMod )
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
TKinetics::KM_InitSplit( long int jb, long int k, const char *kMod )
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
TKinetics::KM_CalcSplit( long int jb, long int k, const char *kMod )
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
TKinetics::KM_SetAMRs( long int jb, long int k,const char *kMod )
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
TKinetics::KM_CalcUptake( long int jb, long int k, const char *kMod )
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
TKinetics::KM_InitUptake( long int jb, long int k, const char *kMod )
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
void TKinetics::Alloc_TKinMet( long int newFI )
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

void TKinetics::Free_TKinMet()
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


//--------------------- End of node_kinetics.cpp ---------------------------


