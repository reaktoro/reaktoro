//-------------------------------------------------------------------
// $Id: s_solmod.h 725 2012-10-02 15:43:37Z kulik $
//
/// \file s_solmod.h
/// Declarations of TSolMod and derived classes implementing built-in models
/// of mixing in fluid, liquid, aqueous and solid-solution phases

// Copyright (C) 2003-2014  T.Wagner, D.Kulik, S.Dmitrieva, F.Hingerl, S.Churakov
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
//------------------------------------------------------------------
//

#ifndef _s_solmod_h_
#define _s_solmod_h_
#include <cstring>
#include <vector>
#include <iostream>

// re-declaration of enums below required for GEMS3K
// dc_class_codes for fluids will be replaced by tp_codes
enum fluid_mix_rules {  /// codes for mixing rules in EoS models (see m_phase.h)
    MR_UNDEF_ = 'N',
    MR_WAAL_ = 'W',
    MR_CONST_ = 'C',
    MR_TEMP_ = 'T',
    MR_LJ_ = 'J',
    MR_KW1_ = 'K',
    MR_PITZ5_ = '5',
    MR_PITZ6_ = '6',
    MR_PITZ8_ = '8',
    MR_B_RCPT_ = 'R'
};

enum dc_class_codes {  /// codes for fluid types in EoS models (see v_mod.h)
    DC_GAS_H2O_ = 'V',
    DC_GAS_CO2_ = 'C',
    DC_GAS_H2_ = 'H',
    DC_GAS_N2_ = 'N',
    DC_GAS_COMP_ = 'G'
};

enum tp_codes {  /// codes for fluid subroutines in EoS models (see v_mod.h)
    CEM_OFF_ = 'N',
    CEM_GAS_ = 'G',
    CEM_H2O_ = 'V',
    CEM_CO2_ = 'C',
    CEM_CH4_ = 'M',
    CEM_N2_ = 'T',
    CEM_H2_ = 'H',
    CEM_O2_ = 'O',
    CEM_AR_ = 'A',
    CEM_PO_ = 'P',
    CEM_NP_ = 'Q'
};

// ------------------------------------------------------------------

#define MAXPHASENAME 16

/// Base class for subclasses of built-in mixing models.
/// (c) March 2007 DK/TW
struct SolutionData {
    long int NSpecies;  ///< Number of species (end members) in the phase
    long int NParams;   ///< Total number of non-zero interaction parameters
    long int NPcoefs;   ///< Number of coefficients per interaction parameter
    long int MaxOrder;  ///< Maximum order of interaction parameters
    long int NPperDC;   ///< Number of parameters per species (DC)
    long int NSublat;   ///< number of sublattices nS
    long int NMoiet;    ///< number of moieties nM

//    long int NlPhs;     ///< new: Number of linked phases
//    long int NlPhC;     ///< new: Number of linked phase parameter coefficient per link (default 0)
    long int NDQFpDC;   ///< new: Number of DQF parameters per species (end member)
//    long int NrcPpDC;   ///< new: Number of reciprocal parameters per species (end member)

    char Mod_Code;      ///< Code of the mixing model
    char Mix_Code;      ///< Code for specific EoS mixing rule
    char *DC_Codes;     ///< DC class codes for species -> NSpecies
    char (*TP_Code)[6]; ///< Codes for TP correction methods for species ->NSpecies
    long int *arIPx;    ///< Pointer to list of indexes of non-zero interaction parameters

//    long int *arPhLin;  ///< new: indexes of linked phase and link type codes [NlPhs*2] read-only

    double *arIPc;      ///< Table of interaction parameter coefficients
    double *arDCc;      ///< End-member properties coefficients
    double *arMoiSN;    ///< End member moiety- site multiplicity number tables -> NSpecies x NSublat x NMoiet
    double *arSitFr;    ///< Tables of sublattice site fractions for moieties -> NSublat x NMoiet
    double *arSitFj;    ///< new: Table of end member sublattice activity coefficients -> NSpecies x NSublat
    double *arGEX;      ///< Pure-species fugacities, G0 increment terms  -> NSpecies

//    double *lPhc;  ///< new: array of phase link parameters -> NlPhs x NlPhC (read-only)
    double *DQFc;  ///< new: array of DQF parameters for DCs in phases ->  NSpecies x NDQFpDC; (read-only)
//    double *rcpc;  ///< new: array of reciprocal parameters for DCs in phases -> NSpecies x NrcPpDC; (read-only)

    double *arPparc;    ///< Partial pressures -> NSpecies
    double *arWx;       ///< Species (end member) mole fractions ->NSpecies
    double *arlnGam;    ///< Output: activity coefficients of species (end members)   

    // Detailed output on terms of partial end-member properties, allocated in MULTI
    double *arlnDQFt; ///< new: DQF terms adding to overall activity coefficients [Ls_]
    double *arlnRcpt; ///< new: reciprocal terms adding to overall activity coefficients [Ls_]
    double *arlnExet; ///< new: excess energy terms adding to overall activity coefficients [Ls_]
    double *arlnCnft; ///< new: configurational terms adding to overall activity [Ls_]

    double *arVol;      ///< molar volumes of end-members (species) cm3/mol ->NSpecies
    double *aphVOL;     ///< phase volumes, cm3/mol (now obsolete) !!!!!!! check usage!
    double T_k;         ///< Temperature, K (initial)
    double P_bar;       ///< Pressure, bar (initial)
};


class TSolMod
{
	protected:
        char ModCode;   ///< Code of the mixing model
        char MixCode;	///< Code for specific EoS mixing rules
                char *DC_Codes; ///< Class codes of end members (species) ->NComp

        char PhaseName[MAXPHASENAME+1];    ///< Phase name (for specific built-in models)

        long int NComp;   ///< Number of components in the solution phase
        long int NPar;     ///< Number of non-zero interaction parameters
        long int NPcoef;   ///< Number of coeffs per parameter (columns in the aIPc table)
        long int MaxOrd;   ///< max. parameter order (or number of columns in aIPx)
        long int NP_DC;    ///< Number of coeffs per one DC in the phase (columns in aDCc)
        long int NSub;     ///< number of sublattices nS
        long int NMoi;     ///< number of moieties nM

//   long int NlPh;     ///< new: Number of linked phases
//   long int NlPc;     ///< new: Number of linked phase parameter coefficient per link (default 0)
   long int NDQFpc;   ///< new: Number of DQF parameters per species (end member)
//   long int NrcPpc;   ///< new: Number of reciprocal parameters per species (end member)

        //        long int NPTP_DC;  // Number of properties per one DC at T,P of interest (columns in aDC)  !!!! Move to CG EOS subclass
                long int *aIPx;    // Pointer to list of indexes of non-zero interaction parameters
//   long int (*PhLin)[2];  ///< new: indexes of linked phase and link type codes [NlPhs][2] read-only

        double R_CONST; ///< R constant
        double Tk;    	///< Temperature, K
        double Pbar;  	///< Pressure, bar

        double *aIPc;   ///< Table of interaction parameter coefficients
        double *aIP;    ///< Vector of interaction parameters corrected to T,P of interest
        double *aDCc;   ///< End-member properties coefficients
        double *aGEX;   ///< Reciprocal energies, Darken terms, pure fugacities of DC (corrected to TP)
        double *aPparc;  ///< Output partial pressures (activities, fugacities) -> NComp
        double **aDC;   ///< Table of corrected end member properties at T,P of interest  !!!!!! Move to GC EOS subclass!
        double *aMoiSN; ///< End member moiety- site multiplicity number tables -> NComp x NSub x NMoi
        double *aSitFR; ///< Table of sublattice site fractions for moieties -> NSub x NMoi

//    double *lPhcf;  ///< new: array of phase link parameters -> NlPh x NlPc (read-only)
    double *DQFcf;  ///< new: array of DQF parameters for DCs in phases ->  NComp x NDQFpc; (read-only)
//    double *rcpcf;  ///< new: array of reciprocal parameters for DCs in phases -> NComp x NrcPpc; (read-only)

        double *x;      ///< Pointer to mole fractions of end members (provided)
        double *aVol;   ///< molar volumes of species (end members)
        double *phVOL;  ///< phase volume, cm3/mol (now obsolete) !!!!!!!!!!!! Check usage!

        // Results
        // double Gam;   	///< work cell for activity coefficient of end member
        // double lnGamRT;
        // double lnGam;
        double Gex, Hex, Sex, CPex, Vex, Aex, Uex;   ///< molar excess properties of the phase
        double Gid, Hid, Sid, CPid, Vid, Aid, Uid;   ///< molar ideal mixing properties
        double Gdq, Hdq, Sdq, CPdq, Vdq, Adq, Udq;   ///< molar Darken quadratic terms
        double Grs, Hrs, Srs, CPrs, Vrs, Ars, Urs;   ///< molar residual functions (fluids)
        double *lnGamConf, *lnGamRecip, *lnGamEx;    ///< Work arrays for lnGamma components
        double *lnGamma;   ///< Pointer to ln activity coefficients of end members (check that it is collected from three above arrays)

        double **y;       ///< table of moiety site fractions [NSub][NMoi]
        double ***mn;     ///< array of end member moiety-site multiplicity numbers [NComp][NSub][NMoi]
        double *mns;      ///< array of total site multiplicities [NSub]
   double **fjs;     ///< array of site activity coefficients [NComp][NSub]
   double *aSitFj; ///< new: pointer to return table of site activity coefficients NComp x NSub

        // functions for calculation of configurational term for multisite ideal mixing
        void alloc_multisite();
        long int init_multisite();
        void free_multisite();

        /// Functions for calculation of configurational term for multisite ideal mixing
        long int IdealMixing();
        double ideal_conf_entropy();
        void return_sitefr();
        void retrieve_sitefr();


        public:

        /// Generic constructor
        TSolMod( SolutionData *sd );

         /// Generic constructor for DComp/DCthermo
        TSolMod( long int NSpecies,  char Mod_Code,  double T_k, double P_bar );

        /// Destructor
		virtual ~TSolMod();

		virtual long int PureSpecies()
		{
			return 0;
		};

		virtual long int PTparam()
		{
			return 0;
		};

		virtual long int MixMod()
		{
			return 0;
		};

		virtual long int ExcessProp( double *Zex )
		{
			return 0;
		};

		virtual long int IdealProp( double *Zid )
		{
			return 0;
		};

        /// Set new system state
		long int UpdatePT ( double T_k, double P_bar );

        bool testSizes( SolutionData *sd );

        /// Getting phase name
		void GetPhaseName( const char *PhName );

		
        // copy activity coefficients into provided array lngamma
		inline void Get_lnGamma( double* lngamma )		
		{ 
			for( int i=0; i<NComp; i++ )
				lngamma[i] = lnGamma[i]; 
		}


};



/// Subclass for the ideal model (both simple and multi-site)
class TIdeal: public TSolMod
{
            private:

            public:

                    /// Constructor
                    TIdeal( SolutionData *sd );

                    /// Destructor
                    ~TIdeal();

                    /// Calculates T,P corrected interaction parameters
                    long int PTparam();

                    /// Calculates (fictive) activity coefficients
                    long int MixMod();

                    /// Calculates excess properties
                    long int ExcessProp( double *Zex );

                    /// Calculates ideal mixing properties
                    long int IdealProp( double *Zid );

};



/// Churakov & Gottschalk (2003) EOS calculations
/// declaration of EOSPARAM class (used by the TCGFcalc class)
class EOSPARAM
{
	private:
		//static double Told;
		// unsigned long int isize;  // int isize = NComp;
		long int NComp;
		double emix, s3mix;
		double *epspar,*sig3par;
		double *XX;
		double *eps;
		double *eps05;
		double *sigpar;
		double *mpar;
		double *apar;
		double *aredpar;
		double *m2par;
		double **mixpar;

		void allocate();
		void free();

	public:

		double *XX0;

		//EOSPARAM():isize(0),emix(0),s3mix(0),NComp(0){};
		//EOSPARAM(double*data, unsigned nn):isize(0){allocate(nn);init(data,nn);};

		EOSPARAM( double *Xtmp, double *data, long int nn )
			:NComp(nn), emix(0),s3mix(0)
                        { allocate(); init(Xtmp,data,nn); }

		~EOSPARAM()
			{ free(); }

		void init( double*,double *, long int );
                long int NCmp()   {return NComp; }

		double EPS05( long int i)
                        { return eps05[i]; }
		double X( long int i)
                        { return XX[i]; }
		double EPS( long int i)
                        { return eps[i]; }
		double EMIX(void)
                        { return emix; }
		double S3MIX(void)
                        { return s3mix; }

		double MIXS3( long int i, long int j)
		{
			if (i==j) return sig3par[i];
                        if (i<j) return mixpar[i][j];
                            else return mixpar[j][i];
                }

		double MIXES3( long int i, long int j)
		{
			if ( i==j ) return epspar[i];
                        if (i<j) return mixpar[j][i];
                            else return mixpar[i][j];
		};

                double SIG3( long int i){ return sig3par[i]; }
                double M2R( long int i) { return m2par[i]; }
                double A( long int i)   { return apar[i]; }

		long int ParamMix( double *Xin);
};



// -------------------------------------------------------------------------------------
/// Churakov and Gottschalk (2003) EOS calculations.
/// Added 09 May 2003
/// Declaration of a class for CG EOS calculations for fluids
/// Incorporates a C++ program written by Sergey Churakov (CSCS ETHZ)
/// implementing papers by Churakov and Gottschalk (2003a, 2003b)
class TCGFcalc: public TSolMod
{
	private:

                double
                PI_1,    ///< pi
                TWOPI,    ///< 2.*pi
                PISIX,    ///< pi/6.
                TWOPOW1SIX,   ///< 2^(1/6)
                DELTA,
                DELTAMOLLIM,
                R,  NA,  P1,
                PP2, P3, P4,
                P5,  P6, P7,
                P8,  P9, P10,
                AA1, AA2, AA3,
                A4, A5, A6,
                BB1, BB2, BB3,
                B4,  B5,  B6,
                A00, A01, A10,
                A11, A12, A21,
                A22, A23, A31,
                A32, A33, A34;

                //  double PhVol;  // phase volume in cm3
                double *Pparc;     ///< DC partial pressures (pure fugacities)
                double *phWGT;
                double *aX;        ///< DC quantities at eqstate x_j (moles)
                    // double *aGEX;      // Increments to molar G0 values of DCs from pure fugacities
                    // double *aVol;      // DC molar volumes, cm3/mol [L]

                // main work arrays
                EOSPARAM *paar;
                EOSPARAM *paar1;
                double *FugCoefs;
                double *EoSparam;
                double *EoSparam1;
                double (*Cf)[8];   ///< corrected EoS coefficients

                // internal functions
                void alloc_internal();
                void free_internal();
                void set_internal();

                void choose( double *pres, double P,unsigned long int &x1,unsigned long int &x2 );
                double Melt2( double T );
                double Melt( double T );
                void copy( double* sours,double *dest,unsigned long int num );
                void norm( double *X,unsigned long int mNum );
                double RPA( double beta,double nuw );
                double dHS( double beta,double ro );

                inline double fI1_6( double nuw )
                {
                    return (1.+(A4+(A5+A6*nuw)*nuw)*nuw)/
                        ((1.+(AA1+(AA2+AA3*nuw)*nuw)*nuw)*3.);
                }

                inline double fI1_12( double nuw )
                {
                    return (1.+(B4+(B5+B6*nuw)*nuw)*nuw)/
                        ((1.+(BB1+(BB2+BB3*nuw)*nuw)*nuw)*9.);
                }

                inline double fa0( double nuw ,double nu1w2 )
                {
                    return (A00 + A01*nuw)/nu1w2;
                }

                inline double fa1( double nuw ,double nu1w3 )
                {
                    return (A10+(A11+A12*nuw)*nuw)/nu1w3;
                }

                inline double fa2( double nuw ,double nu1w4 )
                {
                    return ((A21+(A22+A23*nuw)*nuw)*nuw)/nu1w4;
                }

                inline double fa3( double nuw ,double nu1w5 )
                {
                    return ((A31+(A32+(A33+A34*nuw)*nuw)*nuw)*nuw)/nu1w5;
                }

                double DIntegral( double T, double ro, unsigned long int IType ); // not used
                double LIntegral( double T, double ro, unsigned long int IType ); // not used
                double KIntegral( double T, double ro, unsigned long int IType ); // not used
                double K23_13( double T, double ro );
                double J6LJ( double T,double ro );
                double FDipPair( double T,double ro,double m2 ); // not used
                double UWCANum( double T,double ro );
                double ZWCANum( double T,double ro );

                double FWCA( double T,double ro );
		double FTOTALMIX( double T_Real,double ro_Real,EOSPARAM* param );
		double UTOTALMIX( double T_Real,double ro_Real,EOSPARAM* param ); // not used
		double ZTOTALMIX( double T_Real,double ro_Real,EOSPARAM* param );
		double PTOTALMIX( double T_Real,double ro_Real,EOSPARAM* param );
		double ROTOTALMIX( double P,double TT,EOSPARAM* param );

		double PRESSURE( double *X, double *param, unsigned long int NN, double ro, double T ); // not used
		double DENSITY( double *X,double *param, unsigned long int NN ,double Pbar, double T );
		long int CGActivCoefRhoT( double *X,double *param, double *act, unsigned long int NN,
				double ro, double T ); // not used

		long int CGActivCoefPT(double *X,double *param,double *act, unsigned long int NN,
				double Pbar, double T, double &roro );

	public:

        /// Constructor
		TCGFcalc( long int NCmp, double Pp, double Tkp );
                TCGFcalc( SolutionData *sd, double *aphWGT, double *arX );

        /// Destructor
		~TCGFcalc();

        /// Calculates of pure species properties (pure fugacities)
		long int PureSpecies( );

        /// Calculates T,P corrected interaction parameters
		long int PTparam();

        /// Calculates activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        ///<  calculates ideal mixing properties
		long int IdealProp( double *Zid );

        /// CGofPureGases, calculates fugacity for 1 species at (X=1)
		long int CGcalcFugPure( double Tmin, float *Cemp, double *FugProps );  // called from DCthermo
		long int CGFugacityPT( double *EoSparam, double *EoSparPT, double &Fugacity,
				double &Volume, double P, double T, double &roro );

        /// Calculates departure functions
		long int CGResidualFunct( double *X, double *param, double *param1, unsigned long int NN,
				double ro, double T );

		double GetDELTA( void )
		{
			return DELTA;
		};
};



// -------------------------------------------------------------------------------------
/// Peng-Robinson-Stryjek-Vera (PRSV) model for fluid mixtures.
/// References: Stryjek and Vera (1986)
/// (c) TW July 2006
class TPRSVcalc: public TSolMod

{
	private:

        double PhVol;   ///< phase volume in cm3
                double *Pparc;  ///< DC partial pressures (pure fugacities)
                    // double *aGEX;   // Increments to molar G0 values of DCs from pure fugacities
                    // double *aVol;   // DC molar volumes, cm3/mol [L]

		// main work arrays
        double (*Eosparm)[6];   ///< EoS parameters
        double (*Pureparm)[4];  ///< Parameters a, b, da/dT, d2a/dT2 for cubic EoS
        double (*Fugpure)[6];   ///< fugacity parameters of pure gas species
        double (*Fugci)[4];     ///< fugacity parameters of species in the mixture

        double **a;		///< arrays of generic parameters
		double **b;
        double **KK;     ///< binary interaction parameter
        double **dKK;    ///< derivative of interaction parameter
        double **d2KK;   ///< second derivative
        double **AA;     ///< binary a terms in the mixture



		// internal functions
		void alloc_internal();
		void free_internal();
		long int AB( double Tcrit, double Pcrit, double omg, double k1, double k2, double k3,
				double &apure, double &bpure, double &da, double &d2a );
		long int FugacityPT( long int i, double *EoSparam );
		long int FugacityPure( long int j ); // Calculates the fugacity of pure species
		long int Cardano( double a2, double a1, double a0, double &z1, double &z2, double &z3 );
		long int MixParam( double &amix, double &bmix );
		long int FugacityMix( double amix, double bmix, double &fugmix, double &zmix, double &vmix );
		long int FugacitySpec( double *fugpure );
		long int ResidualFunct( double *fugpure );
		long int MixingWaals();
		long int MixingConst();
		long int MixingTemp();

	public:

        /// Constructor
		TPRSVcalc( long int NCmp, double Pp, double Tkp );
                TPRSVcalc( SolutionData *sd );

        /// Destructor
		~TPRSVcalc();

        /// Calculates pure species properties (pure fugacities)
		long int PureSpecies();

        /// Calculates T,P corrected interaction parameters
		long int PTparam();

        /// Calculates activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

        /// Calculates pure species properties (called from DCthermo)
		long int PRSVCalcFugPure( double Tmin, float *Cpg, double *FugProps );

};



// -------------------------------------------------------------------------------------
/// Soave-Redlich-Kwong (SRK) model for fluid mixtures.
/// References: Soave (1972); Soave (1993)
/// (c) TW December 2008
class TSRKcalc: public TSolMod

{
	private:

        double PhVol;   ///< phase volume in cm3
                double *Pparc;  ///< DC partial pressures (pure fugacities)
                    // double *aGEX;   // Increments to molar G0 values of DCs from pure fugacities
                    // double *aVol;   // DC molar volumes, cm3/mol [L]

		// main work arrays
        double (*Eosparm)[4];   ///< EoS parameters
        double (*Pureparm)[4];  ///< Parameters a, b, da/dT, d2a/dT2 for cubic EoS
        double (*Fugpure)[6];   ///< Fugacity parameters of pure gas species
        double (*Fugci)[4];     ///< Fugacity parameters of species in the mixture

        double **a;		///< arrays of generic parameters
		double **b;
        double **KK;    ///< binary interaction parameter
        double **dKK;   ///< derivative of interaction parameter
        double **d2KK;  ///< second derivative
        double **AA;    ///< binary a terms in the mixture

		// internal functions
		void alloc_internal();
		void free_internal();
		long int AB( double Tcrit, double Pcrit, double omg, double N,
				double &apure, double &bpure, double &da, double &d2a );
		long int FugacityPT( long int i, double *EoSparam );
		long int FugacityPure( long int j ); // Calculates the fugacity of pure species
		long int Cardano( double a2, double a1, double a0, double &z1, double &z2, double &z3 );
		long int MixParam( double &amix, double &bmix );
		long int FugacityMix( double amix, double bmix, double &fugmix, double &zmix, double &vmix );
		long int FugacitySpec( double *fugpure );
		long int ResidualFunct( double *fugpure );
		long int MixingWaals();
		long int MixingConst();
		long int MixingTemp();

	public:

        /// Constructor
		TSRKcalc( long int NCmp, double Pp, double Tkp );
                TSRKcalc( SolutionData *sd );

        /// Destructor
		~TSRKcalc();

        /// Calculates pure species properties (pure fugacities)
		long int PureSpecies();

        /// Calculates T,P corrected interaction parameters
		long int PTparam();

        /// Calculates activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

        /// Calculates pure species properties (called from DCthermo)
		long int SRKCalcFugPure( double Tmin, float *Cpg, double *FugProps );

};



// -------------------------------------------------------------------------------------
/// Peng-Robinson (PR78) model for fluid mixtures.
/// References: Peng and Robinson (1976); Peng and Robinson (1978)
/// (c) TW July 2009
class TPR78calc: public TSolMod

{
	private:

        double PhVol;   ///< phase volume in cm3
                double *Pparc;  ///< DC partial pressures (pure fugacities)
                    // double *aGEX;   // Increments to molar G0 values of DCs from pure fugacities
                    // double *aVol;   // DC molar volumes, cm3/mol [L]

		// main work arrays
        double (*Eosparm)[4];   ///< EoS parameters
        double (*Pureparm)[4];  ///< Parameters a, b, da/dT, d2a/dT2 for cubic EoS
        double (*Fugpure)[6];   ///< Fugacity parameters of pure gas species
        double (*Fugci)[4];     ///< Fugacity parameters of species in the mixture

        double **a;		///< arrays of generic parameters
		double **b;
        double **KK;    ///< binary interaction parameter
        double **dKK;   ///< derivative of interaction parameter
        double **d2KK;  ///< second derivative
        double **AA;    ///< binary a terms in the mixture

		// internal functions
		void alloc_internal();
		void free_internal();
		long int AB( double Tcrit, double Pcrit, double omg, double N,
				double &apure, double &bpure, double &da, double &d2a );
		long int FugacityPT( long int i, double *EoSparam );
		long int FugacityPure( long int j ); // Calculates the fugacity of pure species
		long int Cardano( double a2, double a1, double a0, double &z1, double &z2, double &z3 );
		long int MixParam( double &amix, double &bmix );
		long int FugacityMix( double amix, double bmix, double &fugmix, double &zmix, double &vmix );
		long int FugacitySpec( double *fugpure );
		long int ResidualFunct( double *fugpure );
		long int MixingWaals();
		long int MixingConst();
		long int MixingTemp();

	public:

        /// Constructor
		TPR78calc( long int NCmp, double Pp, double Tkp );
                TPR78calc( SolutionData *sd );

        /// Destructor
		~TPR78calc();

        /// Calculates pure species properties (pure fugacities)
		long int PureSpecies();

        /// Calculates T,P corrected interaction parameters
		long int PTparam();

        /// Calculates activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

        /// Calculates pure species properties (called from DCthermo)
		long int PR78CalcFugPure( double Tmin, float *Cpg, double *FugProps );

};



// -------------------------------------------------------------------------------------
/// Compensated Redlich-Kwong (CORK) model for fluid mixtures.
/// References: Holland and Powell (1991)
/// (c) TW May 2010
class TCORKcalc: public TSolMod

{
        private:

                // constants and external parameters
                double RR;    ///< gas constant in kbar
                double Pkb;   ///< pressure in kbar
                double PhVol;   ///< phase volume in cm3
                double *Pparc;  ///< DC partial pressures (pure fugacities)
                    // double *aGEX;   // Increments to molar G0 values of DCs from pure fugacities
                    // double *aVol;   // DC molar volumes, cm3/mol [L]

                // internal work data
                double (*Eosparm)[2];   ///< EoS parameters
                double (*Fugpure)[6];   ///< Fugacity parameters of pure gas species
                double (*Fugci)[4];     ///< Fugacity parameters of species in the mixture
                double (*Rho)[11];      ///< density parameters
                char *EosCode;    ///< identifier of EoS routine
                double *phi;
                double *dphi;
                double *d2phi;
                double *dphip;
                double **A;         ///< binary interaction parameters
                double **W;         ///< volume scaled interaction parameters (derivatives)
                double **B;
                double **dB;
                double **d2B;
                double **dBp;

                // internal functions
                void alloc_internal();
                void free_internal();
                long int FugacityPT( long int j, double *EoSparam );
                long int FugacityH2O( long int j );
                long int FugacityCO2( long int j );
                long int FugacityCorresponding( long int j );
                long int VolumeFugacity( long int phState, double pp, double p0, double a, double b, double c,
                        double d, double e, double &vol, double &fc );
                long int Cardano( double cb, double cc, double cd, double &v1, double &v2, double &v3 );
                long int FugacityMix();
                long int ResidualFunct();

        public:

                /// Constructor
                TCORKcalc( long int NCmp, double Pp, double Tkp, char Eos_Code );
                TCORKcalc( SolutionData *sd );

                /// Destructor
                ~TCORKcalc();

                /// Calculates pure species properties (pure fugacities)
                long int PureSpecies();

                /// Calculates T,P corrected interaction parameters
                long int PTparam();

                /// Calculates activity coefficients
                long int MixMod();

                /// Calculates excess properties
                long int ExcessProp( double *Zex );

                /// Calculates ideal mixing properties
                long int IdealProp( double *Zid );

                /// Calculates pure species properties (called from DCthermo)
                long int CORKCalcFugPure( double Tmin, float *Cpg, double *FugProps );

};



// -------------------------------------------------------------------------------------
/// Sterner-Pitzer (STP) model for fluid mixtures.
/// References: Sterner and Pitzer (1994)
/// (c) TW December 2010
class TSTPcalc: public TSolMod

{
        private:

                // constants and external parameters
                double RC, RR, TMIN, TMAX, PMIN, PMAX;
                double Pkbar, Pkb, Pmpa;
                double PhVol;   ///< phase volume in cm3
                double *Pparc;  ///< DC partial pressures (pure fugacities)

                // internal work data
                char *EosCode;
                double *Tc;
                double *Pc;
                double *Psat;
                double *Rhol;
                double *Rhov;
                double *Mw;
                double *Phi;
                double *dPhiD;
                double *dPhiDD;
                double *dPhiT;
                double *dPhiTT;
                double *dPhiDT;
                double *dPhiDDD;
                double *dPhiDDT;
                double *dPhiDTT;
                double (*Fugpure)[7];
                double (*Rho)[11];
                double *phi;
                double *dphi;
                double *d2phi;
                double *dphip;
                double *lng;
                double **cfh;
                double **cfc;
                double **A;
                double **W;
                double **B;
                double **dB;
                double **d2B;
                double **dBp;

                // internal functions
                void alloc_internal();
                void free_internal();
                void set_internal();
                long int UpdateTauP();
                long int FugacityPT( long int j, double *EoSparam );
                long int FugacityH2O( long int j );
                long int FugacityCO2( long int j );
                long int FugacityCorresponding( long int j );
                long int DensityGuess( long int j, double &Delguess );
                long int PsatH2O( long int j );
                long int PsatCO2( long int j );
                long int Pressure( double rho, double &p, double &dpdrho, double **cf );
                long int Helmholtz( long int j, double rho, double **cf );
                long int ResidualFunct();

        public:

                /// Constructor
                TSTPcalc ( long int NCmp, double Pp, double Tkp, char Eos_Code );
                TSTPcalc ( SolutionData *sd );

                /// Destructor
                ~TSTPcalc();

                /// Calculates pure species properties (pure fugacities)
                long int PureSpecies();

                /// Calculates T,P corrected interaction parameters
                long int PTparam();

                /// Calculates activity coefficients
                long int MixMod();

                /// Calculates excess properties
                long int ExcessProp( double *Zex );

                /// Calculates ideal mixing properties
                long int IdealProp( double *Zid );

                /// Calculates pure species properties (called from DCthermo)
                long int STPCalcFugPure( double Tmin, float *Cpg, double *FugProps );

};



// -------------------------------------------------------------------------------------
/// Van Laar model for solid solutions.
/// References:  Holland and Powell (2003)
/// (c) TW March 2007

class TVanLaar: public TSolMod
{
	private:
		double *Wu;
		double *Ws;
		double *Wv;
        double *Wpt;   ///< Interaction coeffs at P-T
        double *Phi;   ///< Mixing terms
        double *PsVol; ///< End member volume parameters

		void alloc_internal();
		void free_internal();

	public:

        /// Constructor
                TVanLaar( SolutionData *sd );

        /// Destructor
		~TVanLaar();

        /// Calculates T,P corrected interaction parameters
		long int PTparam();

        /// Calculates of activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
/// Regular model for multicomponent solid solutions.
/// References: Holland and Powell (1993)
/// (c) TW March 2007
class TRegular: public TSolMod
{
	private:
		double *Wu;
		double *Ws;
		double *Wv;
        double *Wpt;   ///< Interaction coeffs at P-T

		void alloc_internal();
		void free_internal();

	public:

        /// Constructor
                TRegular( SolutionData *sd );

        /// Destructor
		~TRegular();

        /// Calculates T,P corrected interaction parameters
        long int PTparam( );

        /// Calculates of activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
/// Redlich-Kister model for multicomponent solid solutions.
/// References: Hillert (1998)
/// (c) TW March 2007
class TRedlichKister: public TSolMod
{
	private:
		double (*Lu)[4];
		double (*Ls)[4];
		double (*Lcp)[4];
		double (*Lv)[4];
		double (*Lpt)[4];

		void alloc_internal();
		void free_internal();

	public:

        /// Constructor
                TRedlichKister( SolutionData *sd );

        /// Destructor
		~TRedlichKister();

        /// Calculates T,P corrected interaction parameters
		long int PTparam();

        /// Calculates activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
/// Non-random two liquid (NRTL) model for liquid solutions.
/// References: Renon and Prausnitz (1968), Prausnitz et al. (1997)
/// (c) TW June 2008
class TNRTL: public TSolMod
{
	private:
		double **Tau;
		double **dTau;
		double **d2Tau;
		double **Alp;
		double **dAlp;
		double **d2Alp;
		double **G;
		double **dG;
		double **d2G;

		void alloc_internal();
		void free_internal();

	public:

        /// Constructor
                TNRTL( SolutionData *sd );

        /// Destructor
		~TNRTL();

        /// Calculates T,P corrected interaction parameters
		long int PTparam();

        /// Calculates activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
/// Wilson model for liquid solutions.
/// References: Prausnitz et al. (1997)
/// (c) TW June 2008
class TWilson: public TSolMod
{
	private:
		double **Lam;
		double **dLam;
		double **d2Lam;

		void alloc_internal();
		void free_internal();

	public:

        /// Constructor
                TWilson( SolutionData *sd );

        /// Destructor
		~TWilson();

        /// Calculates T,P corrected interaction parameters
		long int PTparam();

        /// Calculates activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
/// Berman model for multi-component sublattice solid solutions.
/// To be extended with reciprocal terms.
/// References: Wood and Nicholls (1978); Berman and Brown (1993)
/// (c) DK/TW December 2010, June 2011
class TBerman: public TSolMod
{
        private:
                long int NrcR;   ///< max. possible number of reciprocal reactions (allocated)
                long int Nrc;    ///< number of reciprocal reactions (actual)
                long int *NmoS;  ///< number of different moieties (in end members) on each sublattice
            long int ***XrcM;  ///< Table of indexes of end members, sublattices and moieties involved in
                               ///< reciprocal reactions [NrecR][4][2], two left and two right side.
                               ///< for each of 4 reaction components: j, mark, // s1, m1, s2, m2.

                double *Wu;    ///< Interaction parameter coefficients a
                double *Ws;    ///< Interaction parameter coefficients b (f(T))
                double *Wv;    ///< Interaction parameter coefficients c (f(P))
                double *Wpt;   ///< Interaction parameters corrected at P-T of interest
            double **fjs;      ///< array of site activity coefficients for end members [NComp][NSub]

                double *Grc;  ///< standard molar reciprocal energies (constant)
                double *oGf;   ///< molar Gibbs energies of end-member compounds
                double *G0f;   ///< standard molar Gibbs energies of end members (constant)
            double *DGrc; ///< molar effects of reciprocal reactions [NrecR]
            double *pyp;  ///< Products of site fractions for end members (CEF mod.) [NComp]
//            double *pyn;  // Products of site fractions for sites not in the end member [NComp]
                void alloc_internal();
                void free_internal();
                long int choose( const long int n, const long int k );
                bool CheckThisReciprocalReaction( const long int r, const long int j, long int *xm );
                long int CollectReciprocalReactions2( void );
//                long int CollectReciprocalReactions3( void );
                long int FindIdenticalSublatticeRow(const long int si, const long int ji, const long jp,
                                                    const long int jb, const long int je );
                                              //      long int &nsx, long int *sx, long int *mx );
                long int ExcessPart();
                               ///< Arrays for ideal conf part must exist in base TSolMod instance
                double PYproduct( const long int j );
                long int em_which(const long int s, const long int m , const long int jb, const long int je);
                long int em_howmany( long int s, long int m );
                double ysigma( const long int j, const long int s );
                double KronDelta( const long int j, const long int s, const long int m );
                double dGref_dysigma(const long int l, const long int s, const long int ex_j );
                double dGref_dysm( const long int s, const long m, const long int ex_j );
                double RefFrameTerm( const long int j, double G_ref );
                long int ReciprocalPart();   ///< Calculation of reciprocal contributions to activity coefficients

        public:

                /// Constructor
                TBerman( SolutionData *sd, double *G0 );

                /// Destructor
                ~TBerman();

                /// Calculates T,P corrected interaction parameters
                long int PTparam();

                /// Calculates activity coefficients
                long int MixMod();

                /// Calculates excess properties
                long int ExcessProp( double *Zex );

                /// Calculates ideal mixing properties
                long int IdealProp( double *Zid );

};


// -------------------------------------------------------------------------------------
/// CEF (Calphad) model for multi-component sublattice solid solutions with reciprocal terms
/// References: Sundman & Agren (1981); Lucas et al. (2006); Hillert (1998).
/// (c) DK/SN since August 2014 (still to change the excess Gibbs energy terms).
class TCEFmod: public TSolMod
{
        private:
                long int *NmoS;  ///< number of different moieties (in end members) on each sublattic

                double *Wu;    ///< Interaction parameter coefficients a
                double *Ws;    ///< Interaction parameter coefficients b (f(T))
                double *Wc;    ///< Interaction parameter coefficients b (f(TlnT))
                double *Wv;    ///< Interaction parameter coefficients c (f(P))
                double *Wpt;   ///< Interaction parameters corrected at P-T of interest
                double **fjs;      ///< array of site activity coefficients for end members [NComp][NSub]

                double *Grc;  ///< standard molar reciprocal energies (constant)
                double *oGf;   ///< molar Gibbs energies of end-member compounds
                double *G0f;   ///< standard molar Gibbs energies of end members (constant)
                double *pyp;  ///< Products of site fractions for end members (CEF mod.) [NComp]
//            double *pyn;  // Products of site fractions for sites not in the end member [NComp]
                void alloc_internal();
                void free_internal();
                long int ExcessPart();
                               ///< Arrays for ideal conf part must exist in base TSolMod instance
                double PYproduct( const long int j );
                long int em_which(const long int s, const long int m , const long int jb, const long int je);
                long int em_howmany( long int s, long int m );
                double ysm( const long int j, const long int s );
                double KronDelta( const long int j, const long int s, const long int m );
                double dGref_dysigma(const long int l, const long int s );
                double dGref_dysm( const long int s, const long m );
                double RefFrameTerm( const long int j, double G_ref );
                long int ReciprocalPart();   ///< Calculation of reciprocal contributions to activity coefficients

        public:

                /// Constructor
                TCEFmod( SolutionData *sd, double *G0 );

                /// Destructor
                ~TCEFmod();

                /// Calculates T,P corrected interaction parameters
                long int PTparam();

                /// Calculates activity coefficients
                long int MixMod();

                /// Calculates excess properties
                long int ExcessProp( double *Zex );

                /// Calculates ideal mixing properties
                long int IdealProp( double *Zid );

};


// -------------------------------------------------------------------------------------
/// SIT model reimplementation for aqueous electrolyte solutions.
/// (c) DK/TW June 2009
class TSIT: public TSolMod
{
	private:

        // data objects copied from MULTI
        double *z;    ///< vector of species charges (for aqueous models)
        double *m;    ///< vector of species molalities (for aqueous models)
        double *RhoW;  ///< water density properties
        double *EpsW;  ///< water dielectrical properties

        // internal work objects
        double I;	///< ionic strength
        double A, dAdT, d2AdT2, dAdP;  ///< A term of DH equation (and derivatives)
        double *LnG;  ///< activity coefficient
        double *dLnGdT;  ///< derivatives
		double *d2LnGdT2;
		double *dLnGdP;
        double **E0;  ///< interaction parameter
		double **E1;
		double **dE0;
		double **dE1;
		double **d2E0;
		double **d2E1;

        // internal functions
		double IonicStrength();
		void alloc_internal();
		void free_internal();

	public:

        /// Constructor
                TSIT( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW );

        /// Destructor
		~TSIT();

        /// Calculates activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

        /// Calculation of internal tables (at each GEM iteration)
		long int PTparam();

};



// -------------------------------------------------------------------------------------
/// Pitzer model, Harvie-Moller-Weare (HMW) version, with explicit temperature dependence.
/// References:
/// (c) SD/FH February 2009
class TPitzer: public TSolMod
{

private:
    long int Nc;	 ///< Number of cations
    long int Na;     ///< Number of anions
    long int Nn;     ///< Number of neutral species
    long int Ns;     ///< Total number of aqueous species (without H2O); index of H2O in aq phase
                     ///< Conversion of species indexes between aq phase and Pitzer parameter tables
    long int *xcx;   ///< list of indexes of Nc cations in aqueous phase
    long int *xax;   ///< list of indexes of Na anions in aq phase
    long int *xnx;   ///< list of indexes of Nn neutral species in aq phase
    double *aZ;    ///< Vector of species charges (for aqueous models)
	double *zc;
	double *za;
    double *aM;    ///< Vector of species molality (for aqueous models)
	double *mc;
	double *ma;
	double *mn;
    double *RhoW;  ///< water density properties
    double *EpsW;  ///< water dielectrical properties

        double Aphi, dAphidT, d2AphidT2, dAphidP;  ///< Computing A-Factor
    double I;  ///< Ionic Strength
    double Is;  ///< Ionic Strength square root
    double Ffac; ///< F-Factor
    double Zfac; ///< Z-Term

    // Input parameter arrays
            //for Gex and activity coefficient calculation
    double **Bet0;     ///< Beta0 table for cation-anion interactions [Nc][Na]
    double **Bet1;	   ///< Beta1 table for cation-anion interactions [Nc][Na]
    double **Bet2;	   ///< Beta2 table for cation-anion interactions [Nc][Na]
    double **Cphi;     ///< Cphi  table for cation-anion interactions [Nc][Na]
    double **Lam;      ///< Lam table for neutral-cation interactions [Nn][Nc]
    double **Lam1;     ///< Lam1 table for neutral-anion interactions [Nn][Na]
    double **Theta;    ///< Theta table for cation-cation interactions [Nc][Nc]
    double **Theta1;   ///< Theta1 table for anion-anion interactions [Na][Na]
    double ***Psi;     ///< Psi array for cation-cation-anion interactions [Nc][Nc][Na]
    double ***Psi1;    ///< Psi1 array for anion-anion-cation interactions [Na][Na][Nc]
    double ***Zeta;    ///< Zeta array for neutral-cation-anion interactions [Nn][Nc][Na]


            // Work parameter arrays
            // double *B1;      /// B' table for cation-anion interactions corrected for IS [Nc][Na]
            // double *B2;      /// B table for cation-anion interactions corrected for IS [Nc][Na]
            // double *B3;      /// B_phi table for cation-anion interactions corrected for IS [Nc][Na]
            // double *Phi1;    /// Phi' table for anion-anion interactions corrected for IS [Na][Na]
            // double *Phi2;    /// Phi table for cation-cation interactions corrected for IS [Nc][Nc]
            // double *Phi3;    /// PhiPhi table for anion-anion interactions corrected for IS [Na][Na]
            // double *C;       /// C table for cation-anion interactions corrected for charge [Nc][Na]
            // double *Etheta;  /// Etheta table for cation-cation interactions [Nc][Nc]
            // double *Ethetap; /// Etheta' table for anion-anion interactions [Na][Na]
            // double bk[21];   /// work space
            // double dk[21];   /// work space

    /// McInnes parameter array and gamma values
	double *McI_PT_array;
	double *GammaMcI;

	enum eTableType
	{
		bet0_ = -10, bet1_ = -11, bet2_ = -12, Cphi_ = -20, Lam_ = -30, Lam1_ = -31,
		Theta_ = -40,  Theta1_ = -41, Psi_ = -50, Psi1_ = -51, Zeta_ = -60
	};

    // internal setup
	void calcSizes();
	void alloc_internal();
	void free_internal();

    /// build conversion of species indexes between aq phase and Pitzer parameter tables
	void setIndexes();
	double setvalue(long int ii, int Gex_or_Sex);

    // internal calculations
    /// Calculation of Etheta and Ethetap values
	void Ecalc( double z, double z1, double I, double DH_term,
					double& Etheta, double& Ethetap );
	inline long int getN() const
	{
		return Nc+Na+Nn;
	}

	double Z_Term( );
	double IonicStr( double& I );
	void getAlp( long int c, long int a, double& alp, double& alp1 );
	double get_g( double x_alp );
	double get_gp( double x_alp );
	double G_ex_par5( long int ii );
	double G_ex_par8( long int ii );
	double S_ex_par5( long int ii );
	double S_ex_par8( long int ii );
	double CP_ex_par5( long int ii );
	double CP_ex_par8( long int ii );
	double F_Factor( double DH_term );
	double lnGammaN( long int N );
	double lnGammaM( long int M, double DH_term );
	double lnGammaX( long int X, double DH_term );
	double lnGammaH2O( double DH_term );

    /// Calc vector of interaction parameters corrected to T,P of interest
	void PTcalc( int Gex_or_Sex );

    /// Calculation KCl activity coefficients for McInnes scaling
	double McInnes_KCl();

	inline long int getIc( long int jj )
    {
		for( long int ic=0; ic<Nc; ic++ )
			if( xcx[ic] == jj )
				return ic;
		return -1;
    }

	inline long int getIa( long int jj )
    {
		for( long int ia=0; ia<Na; ia++ )
			if( xax[ia] == jj )
				return ia;
		return -1;
    }

	inline long int getIn( long int jj )
    {
		for( long int in=0; in<Nn; in++ )
			if( xnx[in] == jj )
				return in;
		return -1;
    }

    inline double p_sum( double* arr, long int *xx, long int Narr )
    {
		double sum_ =0.;
		for( long int i=0; i<Narr; i++ )
          sum_ += arr[xx[i]];
		return sum_;
    }

public:

    /// Constructor
        TPitzer( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW );

    /// Destructor
	~TPitzer();

    /// Calculation of T,P corrected interaction parameters
	long int PTparam();


    long int MixMod();

    /// Calculates activity coefficients
	long int Pitzer_calc_Gamma();
	long int Pitzer_McInnes_KCl();

    /// Calculates excess properties
    long int ExcessProp( double *Zex );

    /// Calculates ideal mixing properties
	long int IdealProp( double *Zid );

	void Pitzer_test_out( const char *path, double Y );

};



// -------------------------------------------------------------------------------------
/// Extended universal quasi-chemical (EUNIQUAC) model for aqueous electrolyte solutions.
/// References: Nicolaisen et al. (1993), Thomsen et al. (1996), Thomsen (2005)
/// (c) TW/FH May 2009
class TEUNIQUAC: public TSolMod
{
	private:

        // data objects copied from MULTI
        double *z;   ///< species charges
        double *m;   ///< species molalities
        double *RhoW;  ///< water density properties
        double *EpsW;  ///< water dielectrical properties

        // internal work objects
        double *R;   ///< volume parameter
        double *Q;   ///< surface parameter
		double *Phi;
		double *Theta;
        double **U;   ///< interaction energies
        double **dU;   ///< first derivative
        double **d2U;   ///< second derivative
		double **Psi;
		double **dPsi;
		double **d2Psi;
        double IS;  ///< ionic strength
        double A, dAdT, d2AdT2, dAdP;  ///< A term of DH equation (and derivatives)

        ///< objects needed for debugging output
		double gammaDH[200];
		double gammaC[200];
		double gammaR[200];

        // internal functions
		void alloc_internal();
		void free_internal();
		long int IonicStrength();

	public:

        /// Constructor
                TEUNIQUAC( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW );

        /// Destructor
		~TEUNIQUAC();

        /// Calculates T,P corrected interaction parameters
		long int PTparam();

        /// Calculates activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

		void Euniquac_test_out( const char *path );

};

// -------------------------------------------------------------------------------------
// ELVIS activity model for aqueous electrolyte solutions
// (c) FFH Aug 2011

class TELVIS: public TSolMod
{
        private:
                // data objects copied from MULTI
                double *z;   							// species charges
                double *m;   							// species molalities
                double *RhoW;  							// water density properties
                double *EpsW;  							// water dielectrical properties
                double aDH;								// averaged ion size term in DH term
                double A, dAdT, d2AdT2, dAdP;  			// A term of DH equation (and derivatives)
                double B, dBdT, d2BdT2, dBdP;  			// B term of DH equation (and derivatives)

                double **beta0;
                double **beta1;
                double **alpha;

                double **coord;         // coordinaiton number parameter

                double **RA;
                double **RC;
                double **QA;
                double **QC;

                double CN;

#ifdef ELVIS_SPEED
#define ELVIS_NCOMP 10
                // internal work objects
                double R[ELVIS_NCOMP];
                double Q[ELVIS_NCOMP];
                double Phi[ELVIS_NCOMP];
                double Theta[ELVIS_NCOMP];

                double EffRad[ELVIS_NCOMP];

                double dRdP[ELVIS_NCOMP];
                double dRdT[ELVIS_NCOMP];
                double d2RdT2[ELVIS_NCOMP];
                double dQdP[ELVIS_NCOMP];
                double dQdT[ELVIS_NCOMP];
                double d2QdT2[ELVIS_NCOMP];

                double WEps[ELVIS_NCOMP][ELVIS_NCOMP];
                double U[ELVIS_NCOMP][ELVIS_NCOMP];
                double dU[ELVIS_NCOMP][ELVIS_NCOMP];
                double d2U[ELVIS_NCOMP][ELVIS_NCOMP];
                double Psi[ELVIS_NCOMP][ELVIS_NCOMP];
                double dPsi[ELVIS_NCOMP][ELVIS_NCOMP];
                double d2Psi[ELVIS_NCOMP][ELVIS_NCOMP];
                double TR[ELVIS_NCOMP][4];

                double U[ELVIS_NCOMP][ELVIS_NCOMP];
                double dUdP[ELVIS_NCOMP][ELVIS_NCOMP];
                double dUdT[ELVIS_NCOMP][ELVIS_NCOMP];
                double d2UdT2[ELVIS_NCOMP][ELVIS_NCOMP];

                double ELVIS_lnGam_DH[ELVIS_NCOMP];
                double ELVIS_lnGam_Born[ELVIS_NCOMP];
                double ELVIS_OsmCoeff_DH[ELVIS_NCOMP];
                double ELVIS_lnGam_UNIQUAC[ELVIS_NCOMP];

#endif

#ifndef ELVIS_SPEED
                // internal work objects
                double *R;   							// volume parameter
                double *Q;   							// surface parameter
                double *Phi;
                double *Theta;
                double *EffRad; 						// effective ionic radii
                double **U;   							// interaction energies
                double **dU;   							// first derivative
                double **d2U;   						// second derivative
                double **Psi;
                double **dPsi;
                double **d2Psi;
                double **TR; 							// TR interpolation parameter array
                double **WEps;							// indices for electrolyte specific permittivity calculation

                double* dRdP;
                double* dRdT;
                double* d2RdT2;
                double* dQdP;
                double* dQdT;
                double* d2QdT2;

                double** dUdP;
                double** dUdT;
                double** d2UdT2;

                double* ELVIS_lnGam_DH;
                double* ELVIS_lnGam_Born;
                double* ELVIS_OsmCoeff_DH;
                double* ELVIS_lnGam_UNIQUAC;
#endif

                double IS;  							// ionic strength
                double molT;  							// total molality of aqueous species (except water solvent)
                double molZ;  							// total molality of charged species


                // objects needed for debugging output
                double gammaDH[200];
                double gammaBorn[200];
                double gammaQUAC[200];
                double gammaC[200];
                double gammaR[200];

                // internal functions
                void alloc_internal();
                void free_internal();

                long int IonicStrength();

                // activity coefficient contributions
                void ELVIS_DH(double* ELVIS_lnGam_DH, double* ELVIS_OsmCoeff_DH);
                void ELVIS_Born(double* ELVIS_lnGam_Born);
                void ELVIS_UNIQUAC(double* ELVIS_lnGam_UNIQUAC);

                // Osmotic coefficient
                double Int_OsmCoeff();
                void molfrac_update();
                double FinDiff( double m_j, int j  ); 	// Finite Difference of lnGam of electrolyte 'j' with respect to its molality 'm[j]';
                double CalcWaterAct();

                // Apparent molar volume
//                double FinDiffVol( double m_j, void* params ); 					// Finite differences of mean lnGam with respect to pressure

                double trapzd( const double lower_bound, const double upper_bound, int& n, long int& species, int select ); 	// from Numerical Recipes in C, 2nd Ed.
                double qsimp( const double lower_bound, const double upper_bound, long int& species, int select ); 			// from Numerical Recipes in C, 2nd Ed.


        public:
                // Constructor
                TELVIS( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW );

                // Destructor
                ~TELVIS();

                // calculates T,P corrected interaction parameters
                long int PTparam();

                // calculates activity coefficients and osmotic coefficient by
                // numerical Integration of Bjerrum Relation
                long int MixMod();
                long int CalcAct();

                // Compute apparent molar volume of electrolyte
                double App_molar_volume();

                // calculates excess properties
                long int ExcessProp( double *Zex );

                // calculates ideal mixing properties
                long int IdealProp( double *Zid );

                // plot debug results
                void TELVIS_test_out( const char *path, const double M ) const;

                // ELVIS_FIT: get lnGamma array
                void get_lnGamma( std::vector<double>& ln_gamma );

                double FinDiffVol( double m_j, int j ); 					// Finite differences of mean lnGam with respect to pressure

};


// -------------------------------------------------------------------------------------
/// Extended Debye-Hueckel (EDH) model for aqueous electrolyte solutions, Helgesons variant.
/// References: Helgeson et al. (1981); Oelkers and Helgeson (1990); Pokrovskii and Helgeson (1995; 1997a; 1997b)
/// (c) TW July 2009
class THelgeson: public TSolMod
{
	private:

        // status flags copied from MULTI
        long int flagH2O;  ///< flag for water
        long int flagNeut;  ///< flag for neutral species
        long int flagElect;  ///< flag for selection of background electrolyte model

        // data objects copied from MULTI
        double *z;   ///< species charges
        double *m;   ///< species molalities
        double *RhoW;  ///< water density properties
        double *EpsW;  ///< water dielectrical properties
        double *an;  ///< individual ion size-parameters
        double *bg;  ///< individual extended-term parameters
        double ac;  ///< common ion size parameters
        double bc;  ///< common extended-term parameter

        // internal work objects
        double ao, daodT, d2aodT2, daodP;  ///< ion-size parameter (TP corrected)
        double bgam, dbgdT, d2bgdT2, dbgdP;  ///< extended-term parameter (TP corrected)
        double *LnG;  ///< activity coefficient
        double *dLnGdT;  ///< derivatives
		double *d2LnGdT2;
		double *dLnGdP;
        double IS;  ///< ionic strength
        double molT;  ///< total molality of aqueous species (except water solvent)
        double molZ;  ///< total molality of charged species
        double A, dAdT, d2AdT2, dAdP;  ///< A term of DH equation (and derivatives)
        double B, dBdT, d2BdT2, dBdP;  ///< B term of DH equation (and derivatives)
        double Gf, dGfdT, d2GfdT2, dGfdP;  ///< g function (and derivatives)

        // internal functions
		void alloc_internal();
		void free_internal();
		long int IonicStrength();
		long int BgammaTP();
		long int IonsizeTP();
		long int Gfunction();
		long int GShok2( double T, double P, double D, double beta,
				double alpha, double daldT, double &g, double &dgdP,
				double &dgdT, double &d2gdT2 );

	public:

        /// Constructor
                THelgeson( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW );

        /// Destructor
		~THelgeson();

        /// Calculates T,P corrected interaction parameters
		long int PTparam();

        /// Calculates activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
/// Extended Debye-Hueckel (EDH) model for aqueous electrolyte solutions, Davies variant.
/// References: Langmuir (1997)
/// (c) TW July 2009
class TDavies: public TSolMod
{
	private:

        // status flags copied from MULTI
        long int flagH2O;  ///< flag for water
        long int flagNeut;  ///< flag for neutral species
        long int flagMol;  ///< flag for molality correction

        // data objects copied from MULTI
        double *z;   ///< species charges
        double *m;   ///< species molalities
        double *RhoW;  ///< water density properties
        double *EpsW;  ///< water dielectrical properties

        // internal work objects
        double *LnG;  ///< activity coefficient
        double *dLnGdT;  ///< derivatives
		double *d2LnGdT2;
		double *dLnGdP;
        double IS;  ///< ionic strength
        double molT;  ///< total molality of aqueous species (except water solvent)
        double A, dAdT, d2AdT2, dAdP;  ///< A term of DH equation (and derivatives)

        // internal functions
		void alloc_internal();
		void free_internal();
		long int IonicStrength();

	public:

        /// Constructor
                TDavies( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW );

        /// Destructor
		~TDavies();

        /// Calculates T,P corrected interaction parameters
		long int PTparam();

        /// Calculates activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
/// Debye-Hueckel (DH) limiting law for aqueous electrolyte solutions.
/// References: Langmuir (1997)
/// (c) TW July 2009
class TLimitingLaw: public TSolMod
{
	private:

        // status flags copied from MULTI
        long int flagH2O;  ///< flag for water
        long int flagNeut;  ///< flag for neutral species

        // data objects copied from MULTI
        double *z;   ///< species charges
        double *m;   ///< species molalities
        double *RhoW;  ///< water density properties
        double *EpsW;  ///< water dielectrical properties

        // internal work objects
        double *LnG;  ///< activity coefficient
        double *dLnGdT;  ///< derivatives
		double *d2LnGdT2;
		double *dLnGdP;
        double IS;  ///< ionic strength
        double molT;  ///< total molality of aqueous species (except water solvent)
        double A, dAdT, d2AdT2, dAdP;  ///< A term of DH equation (and derivatives)

        // internal functions
		void alloc_internal();
		void free_internal();
		long int IonicStrength();

	public:

        /// Constructor
                TLimitingLaw( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW );

        /// Destructor
		~TLimitingLaw();

        /// calculates T,P corrected interaction parameters
		long int PTparam();

        /// Calculates activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
/// Two-term Debye-Hueckel (DH) model for aqueous electrolyte solutions.
/// References: Helgeson et al. (1981)
/// uses individual ion-size parameters, optionally individual salting-out coefficients
/// (c) TW July 2009
class TDebyeHueckel: public TSolMod
{
	private:

        // status flags copied from MULTI
        long int flagH2O;  ///< flag for water
        long int flagNeut;  ///< flag for neutral species

        // data objects copied from MULTI
        double *z;   ///< species charges
        double *m;   ///< species molalities
        double *RhoW;  ///< water density properties
        double *EpsW;  ///< water dielectrical properties
        double *an;  ///< individual ion size-parameters
        double *bg;  ///< individual extended-term parameters
        double ac;  ///< common ion size parameters
        double bc;  ///< common extended-term parameter

        // internal work objects
        double ao;  ///< average ion-size parameter
        double *LnG;  ///< activity coefficient
        double *dLnGdT;  ///< derivatives
		double *d2LnGdT2;
		double *dLnGdP;
        double IS;  ///< ionic strength
        double molT;  ///< total molality of aqueous species (except water solvent)
        double A, dAdT, d2AdT2, dAdP;  ///< A term of DH equation (and derivatives)
        double B, dBdT, d2BdT2, dBdP;  ///< B term of DH equation (and derivatives)

        // internal functions
		void alloc_internal();
		void free_internal();
		long int IonicStrength();

	public:

        /// Constructor
                TDebyeHueckel( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW );

        /// Destructor
		~TDebyeHueckel();

        /// Calculates T,P corrected interaction parameters
		long int PTparam();

        /// Calculates activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
/// Extended Debye-Hueckel (EDH) model for aqueous electrolyte solutions, Karpovs variant.
/// References: Karpov et al. (1997); Helgeson et al. (1981); Oelkers and Helgeson (1990);
/// Pokrovskii and Helgeson (1995; 1997a; 1997b)
/// (c) TW July 2009
class TKarpov: public TSolMod
{
	private:

        // status flags copied from MULTI
        long int flagH2O;  ///< flag for water
        long int flagNeut;  ///< flag for neutral species
        long int flagElect;  ///< flag for selection of background electrolyte model

        // data objects copied from MULTI
        double *z;   ///< species charges
        double *m;   ///< species molalities
        double *RhoW;  ///< water density properties
        double *EpsW;  ///< water dielectrical properties
        double *an;  ///< individual ion size-parameters at T,P
        double *bg;  ///< individual extended-term parameters
        double ac;  ///< common ion size parameters
        double bc;  ///< common extended-term parameter

        // internal work objects
        double ao;  ///< average ion-size parameter
        double bgam, dbgdT, d2bgdT2, dbgdP;  ///< extended-term parameter (TP corrected)
        double *LnG;  ///< activity coefficient
        double *dLnGdT;  ///< derivatives
		double *d2LnGdT2;
		double *dLnGdP;
        double IS;  ///< ionic strength
        double molT;  ///< total molality of aqueous species (except water solvent)
        double molZ;  ///< total molality of charged species
        double A, dAdT, d2AdT2, dAdP;  ///< A term of DH equation (and derivatives)
        double B, dBdT, d2BdT2, dBdP;  ///< B term of DH equation (and derivatives)
        double Gf, dGfdT, d2GfdT2, dGfdP;  ///< g function (and derivatives)

        // internal functions
		void alloc_internal();
		void free_internal();
		long int IonicStrength();
		long int BgammaTP();
		long int IonsizeTP();
		long int Gfunction();
		long int GShok2( double T, double P, double D, double beta,
				double alpha, double daldT, double &g, double &dgdP,
				double &dgdT, double &d2gdT2 );

	public:

        /// Constructor
                TKarpov( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW );

        /// Destructor
		~TKarpov();

        /// Calculates T,P corrected interaction parameters
		long int PTparam();

        /// Calculates activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
/// Extended Debye-Hueckel (EDH) model for aqueous electrolyte solutions, Shvarov variant.
/// References: Shvarov (2007); Oelkers and Helgeson (1990);
/// Pokrovskii and Helgeson (1995; 1997a; 1997b)
/// (c) TW July 2009
class TShvarov: public TSolMod
{
	private:

        // status flags copied from MULTI
        long int flagH2O;  ///< new flag for water
        long int flagNeut;  ///< new flag for neutral species
        long int flagElect;  ///< flag for selection of background electrolyte model

        // data objects copied from MULTI
        double *z;   ///< species charges
        double *m;   ///< species molalities
        double *RhoW;  ///< water density properties
        double *EpsW;  ///< water dielectrical properties
        double *bj;  ///< individual ion parameters
        double ac;  ///< common ion size parameters
        double bc;  ///< common extended-term parameter

        // internal work objects
        double ao, daodT, d2aodT2, daodP;  ///< ion-size parameter (TP corrected)
        double bgam, dbgdT, d2bgdT2, dbgdP;  ///< extended-term parameter (TP corrected)
        double *LnG;  ///< activity coefficient
        double *dLnGdT;  ///< derivatives
		double *d2LnGdT2;
		double *dLnGdP;
        double IS;  ///< ionic strength
        double molT;  ///< total molality of aqueous species (except water solvent)
        double A, dAdT, d2AdT2, dAdP;  ///< A term of DH equation (and derivatives)
        double B, dBdT, d2BdT2, dBdP;  ///< B term of DH equation (and derivatives)
        double Gf, dGfdT, d2GfdT2, dGfdP;  ///< g function (and derivatives)

        // internal functions
		void alloc_internal();
		void free_internal();
		long int IonicStrength();
		long int BgammaTP();
		long int IonsizeTP();
		long int Gfunction();
		long int GShok2( double T, double P, double D, double beta,
				double alpha, double daldT, double &g, double &dgdP,
				double &dgdT, double &d2gdT2 );

	public:

        /// Constructor
                TShvarov( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW );

        /// Destructor
		~TShvarov();

        /// Calculates T,P corrected interaction parameters
		long int PTparam();

        /// Calculates activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
/// Class for hardcoded models for solid solutions.
/// (c) TW January 2009
class TModOther: public TSolMod
{
	private:
        double PhVol;   ///< phase volume in cm3
                    // double *Pparc;  /// DC partial pressures/ pure fugacities, bar (Pc by default) [0:L-1]
                    // double *aGEX;   /// Increments to molar G0 values of DCs from pure fugacities or DQF terms, normalized [L]
                    // double *aVol;   /// DC molar volumes, cm3/mol [L]
        double *Gdqf;	///< DQF correction terms
		double *Hdqf;
		double *Sdqf;
		double *CPdqf;
		double *Vdqf;

		void alloc_internal();
		void free_internal();

	public:

        /// Constructor
                TModOther( SolutionData *sd, double *dW, double *eW );

        /// Destructor
		~TModOther();

        /// Calculates pure species properties (pure fugacities, DQF corrections)
		long int PureSpecies();

        /// Calculates T,P corrected interaction parameters
		long int PTparam();

        /// Calculates activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

                // functions for individual models (under construction)
                long int Amphibole1();
                long int Biotite1();
                long int Chlorite1();
                long int Clinopyroxene1();
                long int Feldspar1();
                long int Feldspar2();
                long int Garnet1();
                long int Muscovite1();
                long int Orthopyroxene1();
                long int Staurolite1();
                long int Talc();

};



// -------------------------------------------------------------------------------------
/// Ternary Margules (regular) model for solid solutions.
/// References: Anderson and Crerar (1993); Anderson (2006)
/// (c) TW/DK June 2009
class TMargules: public TSolMod
{
	private:

		double WU12, WS12, WV12, WG12;
		double WU13, WS13, WV13, WG13;
		double WU23, WS23, WV23, WG23;
		double WU123, WS123, WV123, WG123;

	public:

        /// Constructor
                TMargules( SolutionData *sd );

        /// Destructor
		~TMargules();

        /// Calculates T,P corrected interaction parameters
		long int PTparam( );

        /// Calculates of activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
/// Binary Margules (subregular) model for solid solutions.
/// References: Anderson and Crerar (1993); Anderson (2006)
/// (c) TW/DK June 2009
class TSubregular: public TSolMod
{
	private:

		double WU12, WS12, WV12, WG12;
		double WU21, WS21, WV21, WG21;

	public:

        /// Constructor
                TSubregular( SolutionData *sd );

        /// Destructor
		~TSubregular();

        /// Calculates T,P corrected interaction parameters
		long int PTparam( );

        /// Calculates of activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



// -------------------------------------------------------------------------------------
/// Binary Guggenheim (Redlich-Kister) model for solid solutions.
/// References: Anderson and Crerar (1993); Anderson (2006)
/// uses normalized (by RT) interaction parameters
/// (c) TW/DK June 2009
class TGuggenheim: public TSolMod
{
	private:

		double a0, a1, a2;

	public:

        /// Constructor
                TGuggenheim( SolutionData *sd );

        /// Destructor
		~TGuggenheim();

        /// Calculates T,P corrected interaction parameters
		long int PTparam( );

        /// Calculates of activity coefficients
		long int MixMod();

        /// Calculates excess properties
		long int ExcessProp( double *Zex );

        /// Calculates ideal mixing properties
		long int IdealProp( double *Zid );

};



#endif

/// _s_solmod_h
