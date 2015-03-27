//-------------------------------------------------------------------
// $Id: s_activity3.cpp 986 2014-08-31 16:06:28Z kulik $
//
/// \file ipm_chemical3.cpp
/// Implementation of chemistry-specific functions (concentrations,
/// activity coefficients, chemical potentials, etc.)
/// for the IPM convex programming Gibbs energy minimization algorithm
//
// Copyright (c) 1992-2012  D.Kulik, T.Wagner, S.Dmitrieva, K.Chudnenko
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
// #include "m_param.h"
#include "activities.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
///  Function for converting internal lnGam[j] value into an external (phase-scale-specific)
///      Gamma[j] if DirFlag = 0 or external into internal value if DirFlag = 1.
///  Returns the respectively corrected external gamma activity coefficient or internal lnGam
///  Returns trivial values (lnGam = 0 or Gamma = 1) when the respective component
///    amount is zero (X[j] == 0) (is this a correct policy for zeroed-off components?)
//
//
double
TActivity::PhaseSpecificGamma( long int j, long int jb, long int je, long int k, long int DirFlag )
{
    double NonLogTerm = 0., NonLogTermW = 0., NonLogTermS = 0., MMC = 0.;
//    SPP_SETTING *pa = &TProfil::pm->pa;

    switch( act.PHC[k] )
    {
      case PH_AQUEL:
           if( act.XF[k] && act.XFA[k] )
           {
                NonLogTerm = 1. - act.XFA[k]/act.XF[k];
                NonLogTermW = 2. - act.XFA[k]/act.XF[k] - act.XF[k]/act.XFA[k];
           }
           break;
      case PH_GASMIX:  case PH_FLUID:   case PH_PLASMA:   case PH_SIMELT:
      case PH_HCARBL:  case PH_SINCOND:  case PH_SINDIS:  case PH_LIQUID:
           break;
      case PH_POLYEL:
      case PH_SORPTION: // only sorbent end-members!
           if( act.XF[k] && act.XFA[k] )
           {
              for( long int jj=jb; jj<je; jj++ )
              {
                if( act.DCC[jj] == DC_SUR_CARRIER ||
                    act.DCC[jj] == DC_SUR_MINAL || act.DCC[jj] == DC_PEL_CARRIER )
                    MMC += act.MM[jj]*act.X[jj]/act.XFA[k];
                    // Weighted-average sorbent mole mass
              }
              NonLogTerm = 1. - act.XFA[k]/act.XF[k];  // Also for sorption phases
              NonLogTermS = 2. - act.XFA[k]/act.XF[k] - act.XF[k]/act.XFA[k];
           }
           break;
       default:
          break; // Phase class code error should be generated here!
    }
#ifdef NOMUPNONLOGTERM
NonLogTerm = 0.0;
NonLogTermS = 0.0;
#endif
        if( DirFlag == 0 )
        {	 // Converting lnGam[j] into Gamma[j]
            if( !act.X[j] && !act.XF[k] )   // && !pm->XF[k]  added by DK 13.04.2012
                        return 1.;
            double Gamma = 1.;
            double lnGamS = act.lnGam[j];

            switch( act.DCC[j] )
            { // Aqueous electrolyte
              case DC_AQ_PROTON: case DC_AQ_ELECTRON:  case DC_AQ_SPECIES: case DC_AQ_SURCOMP:
                lnGamS += NonLogTerm;    // Correction by asymmetry term
                break;
                // calculate molar mass of solvent
            case DC_AQ_SOLVCOM:	    case DC_AQ_SOLVENT:
                lnGamS += NonLogTermW;
                break;
            case DC_GAS_COMP: case DC_GAS_H2O:  case DC_GAS_CO2:
            case DC_GAS_H2: case DC_GAS_N2:
                break;
            case DC_SOL_IDEAL:  case DC_SOL_MINOR:  case DC_SOL_MAJOR: case DC_SOL_MINDEP: case DC_SOL_MAJDEP:
            break;
                // non-electrolyte condensed mixtures
            case DC_SCP_CONDEN: case DC_SUR_MINAL:
                break;
            case DC_SUR_CARRIER: case DC_PEL_CARRIER:
                lnGamS += NonLogTermS;
                break;
                // Sorption phases
            case DC_SSC_A0: case DC_SSC_A1: case DC_SSC_A2: case DC_SSC_A3: case DC_SSC_A4:
            case DC_WSC_A0: case DC_WSC_A1: case DC_WSC_A2: case DC_WSC_A3: case DC_WSC_A4:
            case DC_SUR_GROUP: case DC_SUR_COMPLEX: case DC_SUR_IPAIR:  case DC_IESC_A:
            case DC_IEWC_B:
                lnGamS += NonLogTerm;
                break;
            default:
                break;
            }
            Gamma = exp( lnGamS );
            return Gamma;
        }
        else { // Converting Gamma[j] into lnGam[j]
                if( !act.X[j] && !act.XF[k] )   // && !pm->XF[k]  added by DK 13.04.2012
                        return 0.;
                double Gamma = act.Gamma[j];
                double lnGam = 0.0;  // Cleanup by DK 5.12.2009
                if( Gamma != 1.0 && Gamma > act.lowPosNum )
                    lnGam = log( Gamma );
                switch( act.DCC[j] )
        { // Aqueous electrolyte
                   case DC_AQ_PROTON: case DC_AQ_ELECTRON:  case DC_AQ_SPECIES: case DC_AQ_SURCOMP:
                        lnGam -= NonLogTerm;  // Correction by asymmetry term
                        break;
                   case DC_AQ_SOLVCOM:	    case DC_AQ_SOLVENT:
                        lnGam -= NonLogTermW;
                        break;
               case DC_GAS_COMP: case DC_GAS_H2O: case DC_GAS_CO2: case DC_GAS_H2: case DC_GAS_N2:
                                break;
                   case DC_SOL_IDEAL:  case DC_SOL_MINOR:  case DC_SOL_MAJOR: case DC_SOL_MINDEP:
                   case DC_SOL_MAJDEP:
                        break;
               case DC_SCP_CONDEN: case DC_SUR_MINAL:
                            break;
               case DC_SUR_CARRIER: case DC_PEL_CARRIER:
                        lnGam -= NonLogTermS;
                            break;
                                // Sorption phases
               case DC_SSC_A0: case DC_SSC_A1: case DC_SSC_A2: case DC_SSC_A3: case DC_SSC_A4:
                   case DC_WSC_A0: case DC_WSC_A1: case DC_WSC_A2: case DC_WSC_A3: case DC_WSC_A4:
                   case DC_SUR_GROUP: case DC_SUR_COMPLEX: case DC_SUR_IPAIR:  case DC_IESC_A:
                   case DC_IEWC_B:
                        lnGam -= NonLogTerm;
                        break;
                    default:
                        break;
                }
            return lnGam;
        }
}

//--------------------------------------------------------------------------------
static double ICold=0.;
/// Main call point for calculation of DC activity coefficients (lnGam vector)
///    formerly GammaCalc().
/// Controls various built-in models, as well as generic Phase script calculation
/// LinkMode is a parameter indicating the status of Gamma calculations:
/// LINK_TP_MODE - calculation of equations depending on TP only;
/// LINK_UX_MODE - calculation of equations depending on current
///      IPM approximation of the equilibrium state;
/// LINK_PP_MODE - calculation of integral phase properties after GEMIPM has converged
///		needs to be implemented
/// \return status code (0 if o.k., non-zero values if there were problems
///     with surface complexation models)
long int
TActivity::CalculateActivityCoefficients( long int LinkMode  )
{
    long int k, j, jb, je=0, jpb, jdb, ipb,  jpe=0, jde=0, ipe=0;
    //long int  jmb, jme=0, jsb, jse=0;
    //long int jphl=0, jlphc=0, jdqfc=0,  jrcpc=0;
    char *sMod;
    long int statusGam=0, statusGC=0, statusSACT=0, SmMode = 0;
    double LnGam, pmpXFk;
//    SPP_SETTING *pa = paTProfil;

    // calculating concentrations of species in multi-component phases
    switch( LinkMode )
    {
      case LINK_TP_MODE:  // Built-in functions depending on T,P only
      {
        long int  jdqfc=0,  jrcpc=0; // jphl=0, jlphc=0,
        long int  jmb, jme=0, jsb, jse=0;

        FitVar[3] = 1.0;  // resetting the IPM smoothing factor

         for( k=0; k<act.FIs; k++ )
         { // loop on solution phases
            jb = je;
            je += act.L1[k];
            if( act.L1[k] == 1 )
                continue;
            // Indexes for extracting data from IPx, PMc and DMc arrays
            ipb = ipe;
            ipe += act.LsMod[k*3]*act.LsMod[k*3+1];
            jpb = jpe;
            jpe += act.LsMod[k*3]*act.LsMod[k*3+2];
            jdb = jde;
            jde += act.LsMdc[k*3]*act.L1[k];
            sMod = act.sMod[k];

            jmb = jme;
            jme += act.LsMdc[k*3+1]*act.LsMdc[k*3+2]*act.L1[k];
            jsb = jse;
            jse += act.LsMdc[k*3+1]*act.LsMdc[k*3+2];

                    double nxk = 1./act.L1[k];
            for( j= jb; j<je; j++ )
    		{
                if(act.XF[k] < min( act.DSM, act.PhMinM ) ) // act.lowPosNum )   // workaround 10.03.2008 DK
                        act.Wx[j] = nxk;  // need this eventually to avoid problems with zero mole fractions
                act.fDQF[j] =0.0;  // cleaning fDQF in TP mode!
                act.lnGmo[j] = act.lnGam[j]; // saving activity coefficients in TP mode
       	    }
                // if( sMod[SGM_MODE] != SM_STNGAM ) This should not be the case anymore DK 24.11.2010
                // continue;  // The switch below is for built-in functions only!

            // the following section probably needs to be re-written to allow more flexible
            // combinations of fluid models for pure gases with gE mixing models,
            // scheme should probably be the same as in LINK_UX_MODE, 03.06.2008 (TW)
            switch( act.PHC[k] )
            {
                case PH_AQUEL: case PH_LIQUID: case PH_SINCOND: case PH_SINDIS: case PH_HCARBL:
                case PH_SIMELT: case PH_GASMIX: case PH_PLASMA: case PH_FLUID:
                    SolModCreate( jb, jmb, jsb, jpb, jdb, k, ipb,
                      sMod[SPHAS_TYP], sMod[MIX_TYP], /* jphl, jlphc,*/ jdqfc,  jrcpc  );
                    // new solution models (TW, DK 2007)
            	    SolModParPT( k, sMod[SPHAS_TYP] );
            	    break;
              default:
                    break;
            }

            // move pointers
 //           jphl  += (act.LsPhl[k*2]*2);
 //           jlphc += (act.LsPhl[k*2]*act.LsPhl[k*2+1]);
            jdqfc += (act.LsMdc2[k*3]*act.L1[k]);
            jrcpc += (act.LsMdc2[k*3+1]*act.L1[k]);

          } // k
        }
        break;

      case LINK_PP_MODE: // Mode of calculation of integral solution phase properties
        for( k=0; k<act.FIs; k++ )
        { // loop on solution phases
            jb = je;
            je += act.L1[k];
            sMod = act.sMod[k];
                switch( act.PHC[k] )
            {
              case PH_AQUEL: case PH_LIQUID: case PH_SINCOND: case PH_SINDIS: case PH_HCARBL:
              case PH_SIMELT: case PH_GASMIX: case PH_PLASMA: case PH_FLUID:
           	       SolModExcessProp( k, sMod[SPHAS_TYP] ); // extracting integral phase properties
           	       SolModIdealProp( jb, k, sMod[SPHAS_TYP] );
           	       SolModStandProp( jb, k, sMod[SPHAS_TYP] );
           	       SolModDarkenProp( jb, k, sMod[SPHAS_TYP] );
           	       break;
              default:
                       break;
            }
       } // k
       break;

    case LINK_UX_MODE:
    	// Getting actual smoothing parameter
    	SetSmoothingFactor( SmMode );
    	// calculating DC concentrations after this IPM iteration
        CalculateConcentrations( act.X, act.XF, act.XFA );
        // cleaning activity coefficients
        for( j=0; j<act.L; j++ )
        {
            act.lnGam[j] = 0.;
            act.Gamma[j] = 1.;
        }
/*
        if( act.E && act.LO ) // checking electrostatics
        {
          IS_EtaCalc();  //  calculating charges and charge densities
          if( act.FIat > 0 )
             for( k=0; k<act.FIs; k++ )
             {
               if( act.PHC[k] == PH_POLYEL || act.PHC[k] == PH_SORPTION )
               {  long int ist;
                  for( ist=0; ist<act.FIat; ist++ ) // loop over surface types
                  {
                     act.XpsiA[k][ist] = 0.0;        // cleaning Psi before GouyChapman()
                     act.XpsiB[k][ist] = 0.0;
                     act.XpsiD[k][ist] = 0.0;
                  }  // ist
                }
             }  // k
         } // act.E
*/
        break;
    default:
        Error("CalculateActivityCoefficients()","Invalid LinkMode for a built-in solution model");
    }

    jpe=0; jde=0; ipe=0;
    je=0;
    for( k=0; k<act.FI; k++ )
    { // loop on phases
        jb = je;
        je += act.L1[k];
        if( act.L1[k] == 1 )
            goto END_LOOP;
        sMod = act.sMod[k];
            // if( sMod[SGM_MODE] == SM_IDEAL )  Comm.out 29.11.2010 by DK to introduce multi-site ideal models
            // goto END_LOOP;
        pmpXFk = 0.;  // Added 07.01.05 by KD
        for( j = jb; j < je; j++ )
            pmpXFk += act.X[j];
        if( act.XF[k] < act.DSM ) // Bugfix by KD 09.08.2005 (bug report Th.Matschei)
            act.XF[k] = pmpXFk;

        // Indexes for extracting data from IPx, PMc and DMc arrays
        ipb = ipe;                  // added 07.12.2006 by KD
        ipe += act.LsMod[k*3]*act.LsMod[k*3+1];
        jpb = jpe;
        jpe += act.LsMod[k*3]*act.LsMod[k*3+2];  // Changed 07.12.2006  by KD
        jdb = jde;
        jde += act.LsMdc[k*3]*act.L1[k];
        //jmb = jme;
        //jme += act.LsMdc[k*3+1]*act.LsMdc[k*3+2]*act.L1[k];
        //jsb = jse;
        //jse += act.LsMdc[k*3+1]*act.LsMdc[k*3+2];

   if( LinkMode == LINK_UX_MODE && sMod[SGM_MODE] == SM_STNGAM )
   {    // check that SGM_MODE for adsorption or multi-site ideal SS is not SM_IDEAL in Phase records!
        switch( act.PHC[k] )
        {  // calculating activity coefficients using built-in functions
          case PH_AQUEL:   // DH III variant consistent with HKF
             if( pmpXFk > act.DSM && act.X[act.LO] > act.XwMinM && act.IC > act.ICmin )
             {
                switch( sMod[SPHAS_TYP] )
                {
                    case SM_AQDH3: case SM_AQDH2: case SM_AQDH1: case SM_AQDHS: case SM_AQDHH:
                    case SM_AQDAV: case SM_AQSIT: case SM_AQPITZ: case SM_AQEXUQ: case SM_AQELVIS:
						SolModActCoeff( k, sMod[SPHAS_TYP] );
						break;
					default:
						break;
                }
                ICold = act.IC;
             }
             goto END_LOOP;
             break;
          case PH_GASMIX: case PH_PLASMA: case PH_FLUID:
             if( pmpXFk > act.DSM && act.XF[k] > act.PhMinM )
             {
                 if( sMod[SPHAS_TYP] == SM_CGFLUID )  // CG EoS fluid model
                     SolModActCoeff( k, sMod[SPHAS_TYP] );
                 if( sMod[SPHAS_TYP] == SM_PRFLUID )  // PRSV EoS fluid model
                     SolModActCoeff( k, sMod[SPHAS_TYP] );
                 if( sMod[SPHAS_TYP] == SM_SRFLUID )  // SRK EoS fluid model
                     SolModActCoeff( k, sMod[SPHAS_TYP] );
                 if( sMod[SPHAS_TYP] == SM_PR78FL )  // PR78 EoS fluid model
                     SolModActCoeff( k, sMod[SPHAS_TYP] );
                 if( sMod[SPHAS_TYP] == SM_CORKFL )  // CORK EoS fluid model
                     SolModActCoeff( k, sMod[SPHAS_TYP] );
                 if ( sMod[SPHAS_TYP] == SM_STFLUID )  // STP EoS fluid model
                     SolModActCoeff( k, sMod[SPHAS_TYP] );
             }
             goto END_LOOP;
             break;
         case PH_LIQUID: case PH_SIMELT: case PH_SINCOND: case PH_SINDIS: case PH_HCARBL:
             if( pmpXFk > act.DSM )
             {     // solid and liquid mixtures
                switch( sMod[SPHAS_TYP] )
                {
                    case SM_IDEAL:   // Ideal (multi-site) model (DK 29.11.2010)
                    case SM_BERMAN:  // Non-ideal (multi-site) model (DK 07.12.2010)
                    case SM_CEF:     // multi-site non-ideal ss model (CALPHAD) DK 15.08.2014
                    case SM_REDKIS:  // Redlich-Kister model (binary)
                    case SM_MARGB:   // Subregular Margules model (binary)
                    case SM_MARGT:   // Regular Margules model (ternary)
                    case SM_GUGGENM: // Redlich-Kister model (multicomponent), 2007 (TW)
                    case SM_VANLAAR: // VanLaar model (multicomponent), 2007 (TW)
                    case SM_REGULAR: // Regular model (multicomponent), 2007 (TW)
                    case SM_NRTLLIQ: // NRTL model (multicomponent), 03.06.2007 (TW)
                    case SM_WILSLIQ: // Wilson model (multicomponent), 09.06.2007 (TW)
                        SolModActCoeff( k, sMod[SPHAS_TYP] );
                        break;
                    default:
                        break;
                }
             }
             goto END_LOOP;
             break;
/*        case PH_POLYEL:  // PoissonBoltzmann( q, jb, je, k ); break;
        case PH_SORPTION: // electrostatic potenials from Gouy-Chapman eqn
                if( act.PHC[0] == PH_AQUEL && pmpXFk > act.DSM
                && (act.XFA[0] > act.XwMinM && act.XF[0] > act.DSM ))
                {
                    if( act.E )
                    {
                       statusGC = GouyChapman( jb, je, k );
                    // PoissonBoltzmann( q, jb, je, k )
                    }
                    // Calculating surface activity coefficient terms
                    statusSACT = SurfaceActivityCoeff(  jb, je, jpb, jdb, k );
                }
                break;
*/         default:
            goto END_LOOP;
       } // end switch
   }  // end if LinkMode == LINK_UX_MODE

END_LOOP:
        if( LinkMode == LINK_TP_MODE )  // TP mode - added 04.03.2008 by DK
        {
        	for( j=jb; j<je; j++ )
        	{
                   if( act.XF[k] < act.DSM )   // workaround 10.03.2008 DK
                                act.Wx[j] = 0.0;               //
                   LnGam = act.lnGmo[j];
                   act.lnGam[j] = LnGam;
                   if(  fabs( LnGam ) < 84. )
                       act.Gamma[j] = PhaseSpecificGamma( j, jb, je, k, 0 );
                   else act.Gamma[j] = 1.0;
        	}
        }
        else if(LinkMode == LINK_UX_MODE )  // Bugfix! DK 06.04.11
        { // Real mode for activity coefficients
           double lnGamG;
           for( j=jb; j<je; j++ )
           {
             if( act.DCC[j] == DC_AQ_SURCOMP )  // Workaround for aqueous surface complexes DK 22.07.09
                act.lnGam[j] = 0.0;
             lnGamG = PhaseSpecificGamma( j, jb, je, k, 1 );
             LnGam = act.lnGam[j];
             if( fabs( lnGamG ) > 1e-9 )
            	LnGam += lnGamG;
             act.lnGmo[j] = LnGam;
             if( fabs( LnGam ) < 84. )   // before 26.02.08: < 42.
                    act.Gamma[j] = PhaseSpecificGamma( j, jb, je, k, 0 );
             else act.Gamma[j] = 1.0;

             act.F0[j] = DC_PrimalChemicalPotentialUpdate( j, k );
             act.G[j] = act.G0[j] + act.fDQF[j] + act.F0[j];
           }
        }
    }  // k - end loop over phases

    if( statusGC )
        return statusGC;
    if( statusSACT )
    	return statusSACT;
    return statusGam;
}


//--------------------------------------------------------------------------------
/// Wrapper calls for creating multi-component mixing models for phases
/// using  TSolMod class. Now including multi-site ideal and scripted models
//
void TActivity::SolModCreate( long int jb, long int jmb, long int jsb, long int jpb, long int jdb,
                           long int k, long int ipb, char ModCode, char MixCode,
                           /* long int jphl, long int jlphc, */ long int jdqfc, long int  jrcpc)
{
    double *aZ, *aM;//, *aVol;
    //long int *aIPx;
    //char *DCCp;
    SolutionData sd;

    sd.NSpecies = act.L1[k];          // Number of components (end members) in the phase
    sd.NParams = act.LsMod[k*3];      // Number of interaction parameters
    sd.NPcoefs = act.LsMod[k*3+2];    // and number of coefs per parameter in PMc table
    sd.MaxOrder =  act.LsMod[k*3+1];  // max. parameter order (cols in IPx)
    sd.NPperDC = act.LsMdc[k*3];      // Number of non-ideality coeffs per one DC in multicomponent phase
    sd.NSublat = act.LsMdc[k*3+1];    // Number of site types (sublattices) for multi-site SS model
    sd.NMoiet = act.LsMdc[k*3+2];     // Number of moieties for multi-site SS model
    sd.Mod_Code = ModCode;
    sd.Mix_Code = MixCode;

    //new objects to Phase 06/06/12
//    sd.NlPhs = act.LsPhl[k*2];
//    sd.NlPhC = act.LsPhl[k*2+1];
    sd.NDQFpDC = act.LsMdc2[k*3];
//    sd.NrcPpDC = act.LsMdc2[k*3+1];

    if( phSolMod[k])
        if(  phSolMod[k]->testSizes( &sd ) )
    	{
                phSolMod[k]->UpdatePT( act.Tc, act.Pc );
                return; // using old allocation and setup of the solution model
    	}

    // properties generic to all models
    sd.arIPx = act.IPx+ipb;   // Pointer to list of indexes for non-ideal solutions -> NPar x MaxOrd
    sd.arIPc = act.PMc+jpb;   // Interaction parameter coefficients f(TP) -> NPar x NPcoef
    sd.arDCc = act.DMc+jdb;   // End-member parameter coefficients f(TPX) -> NComp x NP_DC
    sd.arWx = act.Wx+jb;       // End member mole fractions
    sd.arlnGam = act.lnGam+jb; // End member ln activity coeffs

    sd.arlnDQFt = act.lnDQFt+jb; // End member ln activity coeffs
    sd.arlnRcpt = act.lnRcpt+jb; // End member ln activity coeffs
    sd.arlnExet = act.lnExet+jb; // End member ln activity coeffs
    sd.arlnCnft = act.lnCnft+jb; // End member ln activity coeffs

    sd.aphVOL = act.FVOL+k;
    sd.DC_Codes = act.DCC+jb;  // pointer to Dcomp class codes (added 02.05.2010 TW)
    sd.arMoiSN = act.MoiSN+jmb;  // Pointer to sublattice-moiety multiplicity array
    sd.arSitFr = act.SitFr+jsb;  // Pointer to sublattice-moiety multiplicity array
    sd.arGEX = act.fDQF+jb;      // DQF parameters or pure-gas fugacities
    sd.arPparc = act.Pparc+jb;
    sd.TP_Code = &act.dcMod[jb];
    sd.T_k = act.Tc;
    sd.P_bar = act.Pc;

    //new objects to Phase 06/06/12
//    sd.arPhLin = act.PhLin+jphl;
//    sd.lPhc = act.lPhc+ jlphc;
    sd.DQFc = act.DQFc+ jdqfc;
//    sd.rcpc = act.rcpc+ jrcpc;
    //sd.arSitFj =

    // specific properties
    aM = act.Y_m+jb;
    aZ = act.EZ+jb;
    sd.arVol = act.Vol+jb;

    TSolMod* mySM = 0;

   // creating instances of subclasses of TSolMod base class
    switch( ModCode )
    {

        case SM_OTHER:  // Hard-coded solid solution models (selected by phase name)
        {
                TModOther* myPT = new TModOther( &sd, act.denW, act.epsW );
                myPT->GetPhaseName( act.SF[k] );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_VANLAAR:  // Van Laar solid solution model (multicomponent)
        {
                TVanLaar* myPT = new TVanLaar( &sd );
                mySM = (TSolMod*)myPT;
                break;
        }
             // break;

        case SM_REGULAR:  // Regular solid solution model (multicomponent)
        {
                TRegular* myPT = new TRegular( &sd );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_GUGGENM:  // Redlich-Kister solid solution model (multicomponent)
        {
                TRedlichKister* myPT = new TRedlichKister( &sd );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_NRTLLIQ:  // NRTL liquid solution model (multicomponent)
        {
                TNRTL* myPT = new TNRTL( &sd );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_WILSLIQ:  // Wilson liquid solution model (multicomponent)
        {
                TWilson* myPT = new TWilson( &sd );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_MARGT:  // Margules ternary (regular) solid solution model
        {
                TMargules* myPT = new TMargules( &sd );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_MARGB:  // Margules binary (subregular) solid solution model
        {
                TSubregular* myPT = new TSubregular( &sd );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_REDKIS:  // Gugenheim binary (REdlich-Kister) solid solution
        {
                TGuggenheim* myPT = new TGuggenheim( &sd );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_AQPITZ:  // Pitzer aqueous electrolyte model (multicomponent)
        {
                TPitzer* myPT = new TPitzer( &sd, aM, aZ, act.denW, act.epsW );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_AQSIT:  // SIT aqueous electrolyte model (multicomponent)
        {
                TSIT* myPT = new TSIT( &sd, aM, aZ, act.denW, act.epsW );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_AQEXUQ:  // EUNIQUAC aqueous electrolyte model (multicomponent)
        {
                TEUNIQUAC* myPT = new TEUNIQUAC( &sd, aM, aZ, act.denW, act.epsW );
                mySM = (TSolMod*)myPT;
                break;
        }
/*
        case SM_AQELVIS:  // ELVIS aqueous electrolyte model (multicomponent)
        {
                TELVIS* myPT = new TELVIS( &sd, aM, aZ, act.denW, act.epsW );
		mySM = (TSolMod*)myPT;
                break;
        }
*/
        case SM_AQDH3:  // extended Debye-Hueckel aqueous electrolyte model (Karpov version)
        {
                TKarpov* myPT = new TKarpov( &sd, aM, aZ, act.denW, act.epsW );
                mySM = (TSolMod*)myPT;
        	break;
        }

        case SM_AQDH2:   // Debye-Hueckel aqueous electrolyte model
        {
                TDebyeHueckel* myPT = new TDebyeHueckel( &sd, aM, aZ, act.denW, act.epsW );
                mySM = (TSolMod*)myPT;
        	break;
        }

        case SM_AQDH1:   // Debye-Hueckel limiting law aqueous electrolyte model
        {
                TLimitingLaw* myPT = new TLimitingLaw( &sd, aM, aZ, act.denW, act.epsW );
                mySM = (TSolMod*)myPT;
        	break;
        }

        case SM_AQDHS:  // extended Debye-Hueckel aqueous electrolyte model (Shvarov version)
        {
                TShvarov* myPT = new TShvarov( &sd, aM, aZ, act.denW, act.epsW );
                mySM = (TSolMod*)myPT;
        	break;
        }

        case SM_AQDHH:  // extended Debye-Hueckel aqueous electrolyte model (Helgeson version)
        {
                THelgeson* myPT = new THelgeson( &sd, aM, aZ, act.denW, act.epsW );
                mySM = (TSolMod*)myPT;
        	break;
        }

        case SM_AQDAV:  // Davies aqueous electrolyte model (in NEA TDB version)
        {
                TDavies* myPT = new TDavies( &sd, aM, aZ, act.denW, act.epsW );
                mySM = (TSolMod*)myPT;
        	break;
        }

        case SM_PRFLUID:  // PRSV fluid mixture (multicomponent)
        {
                TPRSVcalc* myPT = new TPRSVcalc( &sd );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_CGFLUID:  // CG fluid mixture (multicomponent)
        {
                TCGFcalc* myPT = new TCGFcalc( &sd, act.FWGT+k, act.X+jb  );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_SRFLUID:  // SRK fluid mixture (multicomponent)
        {
                TSRKcalc* myPT = new TSRKcalc( &sd );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_PR78FL:  // PR78 fluid mixture (multicomponent)
        {
                TPR78calc* myPT = new TPR78calc( &sd );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_CORKFL:  // CORK fluid mixture (multicomponent)
        {
                TCORKcalc* myPT = new TCORKcalc( &sd );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_STFLUID:  // STP fluid mixture (H2O-CO2)
        {
                TSTPcalc* myPT = new TSTPcalc ( &sd );
                mySM = (TSolMod*)myPT;
                break;
        }

        case SM_BERMAN:  // Non-ideal (multi-site) model
        {
                TBerman* myPT = new TBerman( &sd, act.G0+jb );
                mySM = (TSolMod*)myPT;
                break;
        }
        case SM_CEF:  // Non-ideal (multi-site) model (CALPHAD)
        {
                TCEFmod* myPT = new TCEFmod( &sd, act.G0+jb );
                mySM = (TSolMod*)myPT;
                break;
        }
        case SM_IDEAL:
        {
                TIdeal* myPT = new TIdeal( &sd );
                myPT->GetPhaseName( act.SF[k] );
                mySM = (TSolMod*)myPT;
                break;
        }
        // case SM_USERDEF:
        default:
        	break;
    }
  	if(phSolMod[k])
            delete phSolMod[k];
        phSolMod[k] = mySM; // set up new pointer for the solution model
}

/// Wrapper call for calculation of temperature and pressure correction
/// uses TSolMod class
void TActivity::SolModParPT( long int k, char ModCode )
{
    // Extended constructor to connect to params, coeffs, and mole fractions
    switch( ModCode )
    {  // solid and liquid solutions
        // case SM_USERDEF:
        case SM_IDEAL: case SM_VANLAAR: case SM_REGULAR: case SM_GUGGENM: case SM_NRTLLIQ:
        case SM_WILSLIQ: /* old ss models */ case SM_MARGT: case SM_MARGB: case SM_REDKIS:
        case SM_BERMAN:  case SM_CEF:  // new ss models
        // aqueous DH models
        case SM_AQDH3: case SM_AQDH2: case SM_AQDH1: case SM_AQDHH: case SM_AQDHS: case SM_AQDAV:
        // aqueous SIT models
        case SM_AQPITZ: case SM_AQSIT: case SM_AQEXUQ: case SM_AQELVIS:
        // fluid (gas) models
        case SM_PRFLUID: case SM_CGFLUID: case SM_SRFLUID: case SM_PR78FL: case SM_CORKFL:
        case SM_STFLUID:
        {    ErrorIf( !phSolMod[k], "SolModParPT: ","Invalid index of phase");
              TSolMod* mySM = phSolMod[k];
              mySM->PTparam();
             break;
        }
        default:
              break;
    }
}

/// Wrapper call for calculation of activity coefficients
/// uses TSolMod class
void TActivity::SolModActCoeff( long int k, char ModCode )
{
    switch( ModCode )
    {   // solid and liquid solutions
        // case SM_USERDEF:
        case SM_IDEAL: case SM_VANLAAR: case SM_REGULAR: case SM_GUGGENM: case SM_NRTLLIQ:
        case SM_WILSLIQ: /* old ss models */ case SM_MARGT: case SM_MARGB: case SM_REDKIS:
        case SM_BERMAN:  case SM_CEF:
        // aqueous DH models
        case SM_AQDH3: case SM_AQDH2: case SM_AQDH1: case SM_AQDHH: case SM_AQDHS: case SM_AQDAV:
        // aqueous SIT models
        case SM_AQPITZ: case SM_AQSIT: case SM_AQEXUQ: case SM_AQELVIS:
        // fluid (gas) models
        case SM_PRFLUID: case SM_CGFLUID: case SM_SRFLUID: case SM_PR78FL: case SM_CORKFL:
        case SM_STFLUID:
        {    ErrorIf( !phSolMod[k], "SolModActCoeff: ","Invalid index of phase");
             TSolMod* mySM = phSolMod[k];
             mySM->MixMod();
             break;
        }
        default:
              break;
    }
}

/// Wrapper call for calculation of bulk phase excess properties
/// uses TSolMod class
void TActivity::SolModExcessProp( long int k, char ModCode )
{

    // order of phase properties: G, H, S, CP, V, A, U
    long int j;
    double Gex, Hex, Sex, CPex, Vex, Aex, Uex;
    double zex[7];

    for (j =0; j<7; j++)
    {
        zex[j] = 0.0;
    }
    // insert cases for old solution and activity models
    switch( ModCode )
    {
        // case SM_USERDEF:  case SM_IDEAL:
        case SM_VANLAAR: case SM_REGULAR: case SM_GUGGENM: case SM_NRTLLIQ:
        case SM_WILSLIQ: /* old ss models */ case SM_MARGT: case SM_MARGB: case SM_REDKIS:
        case SM_BERMAN:  case SM_CEF:
        // aqueous DH models
        case SM_AQDH3: case SM_AQDH2: case SM_AQDH1: case SM_AQDHH: case SM_AQDHS: case SM_AQDAV:
        // aqueous SIT models
        case SM_AQPITZ: case SM_AQSIT: case SM_AQEXUQ: case SM_AQELVIS:
        // fluid (gas) models
        case SM_PRFLUID: case SM_CGFLUID: case SM_SRFLUID: case SM_PR78FL: case SM_CORKFL:
        case SM_STFLUID:
        {    ErrorIf( !phSolMod[k], "SolModExcessProp: ","Invalid index of phase");
              TSolMod* mySM = phSolMod[k];
              mySM->ExcessProp( zex );
              break;
        }
        default:
              break;
    }
    // assignments
    Gex = zex[0];
    Hex = zex[1];
    Sex = zex[2];
    CPex = zex[3];
    Vex = zex[4];
    Aex = zex[5];
    Uex = zex[6];
    act.GPh[k][2] = Gex;
    act.HPh[k][2] = Hex;
    act.SPh[k][2] = Sex;
    act.CPh[k][2] = CPex;
    act.VPh[k][2] = Vex;
    act.APh[k][2] = Aex;
    act.UPh[k][2] = Uex;
}


/// Wrapper call for calculation of bulk phase ideal mixing properties
void TActivity::SolModIdealProp( long int jb, long int k, char ModCode )
{
    // order of phase properties: G, H, S, CP, V, A, U
    long int j;
    double Gid, Hid, Sid, CPid, Vid, Aid, Uid;
    double zid[7];

    for (j=0; j<7; j++)
    {
        zid[j] = 0.0;
    }
    // needs to check for DC class state to catch ideal gas phase case
    switch( ModCode )
    {   // check what solution phase (ideal gas?)
        // case SM_USERDEF:
        case SM_IDEAL: case SM_VANLAAR: case SM_REGULAR: case SM_GUGGENM: case SM_NRTLLIQ:
        case SM_WILSLIQ: /* old ss models */ case SM_MARGT: case SM_MARGB: case SM_REDKIS:
        case SM_BERMAN:  case SM_CEF:
        // aqueous DH models
        case SM_AQDH3: case SM_AQDH2: case SM_AQDH1: case SM_AQDHH: case SM_AQDHS: case SM_AQDAV:
        // aqueous SIT models
        case SM_AQPITZ: case SM_AQSIT: case SM_AQEXUQ: case SM_AQELVIS:
        // fluid (gas) models
        case SM_PRFLUID: case SM_CGFLUID: case SM_SRFLUID: case SM_PR78FL: case SM_CORKFL:
        case SM_STFLUID:
        {    ErrorIf( !phSolMod[k], "SolModIdealProp: ","Invalid index of phase");
             TSolMod* mySM = phSolMod[k];
             mySM->IdealProp( zid );
             break;
        }
        default:
            break;
    }
    // assignments
    Gid = zid[0];
    Hid = zid[1];
    Sid = zid[2];
    CPid = zid[3];
    Vid = zid[4];
    Aid = zid[5];
    Uid = zid[6];
    act.GPh[k][1] = Gid;
    act.HPh[k][1] = Hid;
    act.SPh[k][1] = Sid;
    act.CPh[k][1] = CPid;
    act.VPh[k][1] = Vid;
    act.APh[k][1] = Aid;
    act.UPh[k][1] = Uid;
}

/// Wrapper call for retrieving bulk phase Darken quadratic terms
void TActivity::SolModDarkenProp( long int jb, long int k, char ModCode )
{
    // order of phase properties: G, H, S, CP, V, A, U
    long int j;
    double Gdq, Hdq, Sdq, CPdq, Vdq, Adq, Udq;
    double zdq[7];

    for (j=0; j<7; j++)
    {
        zdq[j] = 0.0;
    }

    // data object for derivative properties needs to be added in Multi and DODs in scripts

    // assignments
    Gdq = zdq[0];
    Hdq = zdq[1];
    Sdq = zdq[2];
    CPdq = zdq[3];
    Vdq = zdq[4];
    Adq = zdq[5];
    Udq = zdq[6];
    act.GPh[k][3] = Gdq;
    act.HPh[k][3] = Hdq;
    act.SPh[k][3] = Sdq;
    act.CPh[k][3] = CPdq;
    act.VPh[k][3] = Vdq;
    act.APh[k][3] = Adq;
    act.UPh[k][3] = Udq;

}

/// Wrapper call for retrieving bulk phase standard state terms
void TActivity::SolModStandProp ( long int jb, long int k, char ModCode )
{
    // order of phase properties: G, H, S, CP, V, A, U
    double Gst=0., Hst=0., Sst=0., CPst=0., Vst=0., Ast=0., Ust=0.;

    // add if statement that checks DC class code (aqueous or not)

    // assignments
    act.GPh[k][0] = Gst;
    act.HPh[k][0] = Hst;
    act.SPh[k][0] = Sst;
    act.CPh[k][0] = CPst;
    act.VPh[k][0] = Vst;
    act.APh[k][0] = Ast;
    act.UPh[k][0] = Ust;

}

//-------------------------------------------------------------------------
/// Internal memory allocation for TSolMod performance optimization
/// (since version 2.3.0)
void TActivity::Alloc_TSolMod( long int newFIs )
{
  if(  phSolMod && ( newFIs == sizeFIs) )
    return;

  Free_TSolMod();
  // alloc memory for all multicomponents phases
  phSolMod = new  TSolMod *[newFIs];
  sizeFIs = newFIs;
 for( long int ii=0; ii<newFIs; ii++ )
    	  phSolMod[ii] = 0;
}

void TActivity::Free_TSolMod()
{
  long int kk;

  if( phSolMod )
  {  for(  kk=0; kk<sizeFIs; kk++ )
      if( phSolMod[kk] )
           delete phSolMod[kk];

      delete[]  phSolMod;
  }
  phSolMod = 0;
  sizeFIs = 0;
}

/// Internal memory allocation for TSorpMod performance optimization
/// (since version 3.4.0)
void TActivity::Alloc_TSorpMod( long int newFIs )
{
  if(  phSorpMod && ( newFIs == sizeFIa) )
    return;

  Free_TSorpMod();
  // allocate memory for all multicomponent phases
  phSorpMod = new  TSorpMod *[newFIs];
  sizeFIa = newFIs;
 for( long int ii=0; ii<newFIs; ii++ )
          phSorpMod[ii] = 0;
}

void TActivity::Free_TSorpMod()
{
  long int kk;

  if( phSorpMod )
  {  for(  kk=0; kk<sizeFIa; kk++ )
      if( phSorpMod[kk] )
           delete phSorpMod[kk];

      delete[]  phSorpMod;
  }
  phSorpMod = 0;
  sizeFIa = 0;
}

/*
/// Internal memory allocation for TKinMet performance optimization
/// (since version 3.3.0)
void TActivity::Alloc_TKinMet( long int newFI )
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

void TActivity::Free_TKinMet()
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
*/

/// New correction of smoothing factor for highly non-ideal systems.
// re-written 18.04.2009 DK+TW
/// Smoothing function choice: AG >= 0.0001 and DGC > -0.0001: old f(IT)
///                            AG >= 0.0001 and DGC <= -0.0001: new f(1/IT)
///                            AG <= -0.0001 and DGC <= -0.0001: new f(1/CD)
/// \param mode 0 - taking single log(CD) value for calculation of smoothing factor SF;
///       1, 2, ...  taking log(CD) average from the moving window of length mode
///       (up to 5 consecutive values)
///
void TActivity::SetSmoothingFactor( long int mode )
{
    double TF=1., al, ag, dg, iim, irf; // rg=0.0;
    long int ir; //, Level, itqF, itq;

    ir = IT;
    irf = (double)ir;
    ag = AG; // act.FitVar[4];
    dg = DGC;
    iim = (double)IIM;

    if( dg > -0.0001 && ag >= 0.0001 ) // Smoothing used in the IPM-2 algorithm
    {					// with some improvements
        if(ag>1) ag=1;
        if(ag<0.1) ag=0.1;
        if(dg>0.15) dg=0.15;
        // if( irf > 1000. )
        //	irf = 1000;
        if( dg <= 0.0 )
          TF = ag;
        else
          TF = ag * ( 1 - pow(1-exp(-dg*irf),60.));
        if(TF < 1e-6 )
          TF = 1e-6;
    }
    else if( dg <= -0.0001 && ag >= 0.0001 )
    {
       // New sigmoid smoothing function of 1/IT
        double logr, inv_r = 1., logr_m;
        dg = fabs( dg );
        if( IT )
          inv_r = 1./(double)IT;
        logr = log( inv_r );
        logr_m = log( 1./iim );
        al = dg + ( ag - dg ) / ( 1. + exp( logr_m - logr ) / dg );
        al += exp( log( 1. - ag ) + logr );
        if( al > 1. )
            al = 1.;
        TF = al;
    }
/*
    else if( dg <= -0.0001 && ag <= -0.0001 )
    {
        double dk, cd;   long int i;
        dg = fabs( dg );
        ag = fabs( ag );
        dk = log( DXM );
        // Checking the mode where it is called
        switch( mode )
        {
          default:
          case 0: // MassBalanceRefinement() after SolveSimplex()
                     cd = log( act.PCI );
                 break;
          case 1:
          case 2:
          case 3:
          case 4:
          case 5: // Getting average (log geometric mean) from sampled CD values
                 cd = 0.0;
                 for(i=0; i < mode; i++ )
                         cd += act.logCDvalues[i];
                     cd /= (double)mode; // 5. - bugfix
                 break;
        }
        al = dg + ( ag - dg ) / ( 1. + exp( dk - cd ) / dg );
        al += exp( log( 1. - ag ) + cd );
        if( al > 1. )
            al = 1.;
        TF = al;
    }
*/
//    if( IT )
      FitVar[3] = TF;
}




//--------------------- End of s_activity3.cpp ---------------------------


