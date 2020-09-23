//-------------------------------------------------------------------
// $Id: m_param.h 887 2013-10-23 08:34:43Z kulik $
//
/// \file m_param.h
/// Declaration of TProfil class, config and calculation functions
//
// Copyright (c) 1995-2012 S.Dmytriyeva, D.Kulik
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

#ifndef _m_param_h_
#define _m_param_h_

#include <cmath>

// Physical constants - see m_param.cpp or ms_param.cpp
extern const double R_CONSTANT, NA_CONSTANT, F_CONSTANT,
    e_CONSTANT,k_CONSTANT, cal_to_J, C_to_K, lg_to_ln, ln_to_lg, H2O_mol_to_kg, Min_phys_amount;

#ifdef IPMGEMPLUGIN

#include "gdatastream.h"
#include "ms_multi.h"
#include "verror.h"

struct BASE_PARAM /// Flags and thresholds for numeric modules
{
   long int
           PC,   ///< Mode of PhaseSelect() operation ( 0 1 2 ... ) { 1 }
           PD,   ///< abs(PD): Mode of execution of CalculateActivityCoefficients() functions { 2 }.
                 ///< Modes: 0-invoke, 1-at MBR only, 2-every MBR it, every IPM it. 3-not MBR, every IPM it.
                 ///< if PD < 0 then use test qd_real accuracy mode
           PRD,  ///< Since r1583/r409: Disable (0) or activate (-5 or less) the SpeciationCleanup() procedure { -5 }
           PSM,  ///< Level of diagnostic messages: 0- disabled (no ipmlog file); 1- errors; 2- also warnings 3- uDD trace { 1 }
           DP,   ///< Maximum allowed number of iterations in the MassBalanceRefinement() procedure {  30 }
           DW,   ///< Since r1583: Activate (1) or disable (0) error condition when DP was exceeded { 1 }
           DT,   ///< Since r1583/r409: DHB is relative for all (0) or absolute (-6 or less ) cutoff for major ICs { 0 }
           PLLG, ///< IPM tolerance for detecting divergence in dual solution { 10; range 1 to 1000; 0 disables the detection }
           PE,   ///< Flag for using electroneutrality condition in GEM IPM calculations { 0 1 }
           IIM   ///< Maximum allowed number of iterations in the MainIPM_Descent() procedure up to 9999 { 1000 }
           ;
         double DG,   ///< Standart total moles { 1e5 }
           DHB,  ///< Maximum allowed relative mass balance residual for Independent Components ( 1e-9 to 1e-15 ) { 1e-10 }
           DS,   ///< Cutoff minimum mole amount of stable Phase present in the IPM primal solution { 1e-12 }
           DK,   ///< IPM-2 convergence threshold for the Dikin criterion (may be set in the interval 1e-6 < DK < 1e-4) { 1e-5 }
           DF,   ///< Threshold for the application of the Karpov phase stability criterion: (Fa > DF) for a lost stable phase { 0.01 }
           DFM,  ///< Threshold for Karpov stability criterion f_a for insertion of a phase (Fa < -DFM) for a present unstable phase { 0.1 }
           DFYw, ///< Insertion mole amount for water-solvent { 1e-6 }
           DFYaq,///< Insertion mole amount for aqueous species { 1e-6 }
           DFYid,///< Insertion mole amount for ideal solution components { 1e-6 }
           DFYr, ///< Insertion mole amount for major solution components { 1e-6 }
           DFYh, ///< Insertion mole amount for minor solution components { 1e-6 }
           DFYc, ///< Insertion mole amount for single-component phase { 1e-6 }
           DFYs, ///< Insertion mole amount used in PhaseSelect() for a condensed phase component  { 1e-7 }
           DB,   ///< Minimum amount of Independent Component in the bulk system composition (except charge "Zz") (moles) (1e-17)
           AG,   ///< Smoothing parameter for non-ideal increments to primal chemical potentials between IPM descent iterations { -1 }
           DGC,  ///< Exponent in the sigmoidal smoothing function, or minimal smoothing factor in new functions { -0.99 }
           GAR,  ///< Initial activity coefficient value for major (M) species in a solution phase before LPP approximation { 1 }
           GAH,  ///< Initial activity coefficient value for minor (J) species in a solution phase before LPP approximation { 1000 }
           GAS,  ///< Since r1583/r409: threshold for primal-dual chem.pot.difference (mol/mol) used in SpeciationCleanup() { 1e-3 }.
                 ///< before: Obsolete IPM-2 balance accuracy control ratio DHBM[i]/b[i], for minor ICs { 1e-3 }
           DNS,  ///< Standard surface density (nm-2) for calculating activity of surface species (12.05)
           XwMin,///< Cutoff mole amount for elimination of water-solvent { 1e-9 }
           ScMin,///< Cutoff mole amount for elimination of solid sorbent {1e-7}
           DcMin,///< Cutoff mole amount for elimination of solution- or surface species { 1e-30 }
           PhMin,///< Cutoff mole amount for elimination of  non-electrolyte solution phase with all its components { 1e-10 }
           ICmin,///< Minimal effective ionic strength (molal), below which the activity coefficients for aqueous species are set to 1. { 3e-5 }
           EPS,  ///< Precision criterion of the SolveSimplex() procedure to obtain the AIA ( 1e-6 to 1e-14 ) { 1e-10 }
           IEPS, ///< Convergence parameter of SACT calculation in sorption/surface complexation models { 0.01 to 0.000001, default 0.001 }
           DKIN; ///< Tolerance on the amount of DC with two-side metastability constraints  { 1e-7 }
    char *tprn;       ///< internal

    void write(fstream& oss);
};

struct SPP_SETTING /// Base Parametres of SP
{
    char ver[TDBVERSION]; ///< Version & Copyright 64
    BASE_PARAM p; //
    void write(fstream& oss);
};

extern SPP_SETTING pa_;

/// Module contains a Flags and thresholds for numeric modules
class TProfil
{

public:

    //static TProfil* pm;

    TMulti* multi;
    MULTI *pmp;
    SPP_SETTING pa;

    TProfil( TMulti* amulti );

    const char* GetName() const
    {
        return "Project";
    }

   void outMulti( GemDataStream& ff, gstring& path  );
   void outMultiTxt( const char *path, bool append=false  );
   void readMulti( GemDataStream& ff );
   void readMulti( const char* path,  DATACH  *dCH );

   double ComputeEquilibriumState( long int& PrecLoops_, long int& NumIterFIA_, long int& NumIterIPM_ );

   inline double HelmholtzEnergy( double x )
   {
       return multi->HelmholtzEnergy(x);
   }

   inline double InternalEnergy( double TC, double P )
   {
       return multi->InternalEnergy( TC, P );
   }

   //void test_G0_V0_H0_Cp0_DD_arrays( long int nT, long int nP );
};

enum QpQdSizes {   // see m_phase.h
   QPSIZE = 180,    // This enum is for GEMS3K only!
   QDSIZE = 60
};

#else

#include "v_mod.h"
#include "ms_rmults.h"
#include "ms_mtparm.h"
#include "ms_system.h"
#include "ms_multi.h"
#include "ms_calc.h"

class GemDataStream;
class QWidget;

extern long int showMss;

struct BASE_PARAM
{ // Flags and thresholds for numeric modules
  short
    PC,   // Mode of PhaseSelect() operation ( 0 1 2 ... ) { 1 }
    PD,   // abs(PD): Mode of execution of CalculateActivityCoefficients() functions { 2 }
          // Mode of CalculateActivityCoefficients(): 0-invoke, 1-at EFD only, 2-every EFD it, every IPM it. 3-not EFD, every IPM it.
          // if PD < 0 then use test qd_real accuracy mode
    PRD,  // Since r1583/r409: Disable (0) or activate (-5 or less) the SpeciationCleanup() procedure { -5 }
    PSM,  // Level of diagnostic messages: 0- disabled (no ipmlog file); 1- errors; 2- also warnings 3- uDD trace { 1 }
    DP,   // Maximum allowed number of iterations in the MassBalanceRefinement() procedure {  30 }
    DW,   // Since r1583: Activate (1) or disable (0) error condition when DP was exceeded { 1 }
    DT,   // Since r1583/r409: DHB is relative for all (0) or absolute (-6 or less ) cutoff for major ICs { 0 }
    PLLG, // IPM tolerance for detecting divergence in dual solution { 10; range 1 to 1000; 0 disables the detection }
    PE,   // Flag for using electroneutrality condition in GEM IPM calculations { 0 1 }
    IIM   // Maximum allowed number of iterations in the MainIPM_Descent() procedure up to 9999 { 1000 }
    ;
  double DG,   // Standart total moles { 1e5 }
    DHB,  // Maximum allowed relative mass balance residual for Independent Components ( 1e-9 to 1e-15 ) { 1e-10 }
    DS,   // Cutoff minimum mole amount of stable Phase present in the IPM primal solution { 1e-12 }
    DK,   // IPM-2 convergence threshold for the Dikin criterion (may be set in the interval 1e-6 < DK < 1e-4) { 1e-5 }
    DF,   // Threshold for the application of the Karpov phase stability criterion: (Fa > DF) for a lost stable phase { 0.01 }
    DFM,  // Threshold for Karpov stability criterion f_a for insertion of a phase (Fa < -DFM) for a present unstable phase { 0.1 }
    DFYw, // Insertion mole amount for water-solvent { 1e-6 }
    DFYaq,// Insertion mole amount for aqueous species { 1e-6 }
    DFYid,// Insertion mole amount for ideal solution components { 1e-6 }
    DFYr, // Insertion mole amount for major solution components { 1e-6 }
    DFYh, // Insertion mole amount for minor solution components { 1e-6 }
    DFYc, // Insertion mole amount for single-component phase { 1e-6 }
    DFYs, // Insertion mole amount used in PhaseSelect() for a condensed phase component  { 1e-7 }
    DB,   // Minimum amount of Independent Component in the bulk system composition (except charge "Zz") (moles) (1e-17)
    AG,   // Smoothing parameter for non-ideal increments to primal chemical potentials between IPM descent iterations { -1 }
    DGC,  // Exponent in the sigmoidal smoothing function, or minimal smoothing factor in new functions { -0.99 }
    GAR,  // Initial activity coefficient value for major (M) species in a solution phase before LPP approximation { 1 }
    GAH,  // Initial activity coefficient value for minor (J) species in a solution phase before LPP approximation { 1000 }
    GAS,  // Since r1583: threshold for primal-dual chem.pot.difference (mol/mol) used in SpeciationCleanup() { 1e-3 }
          // before: Obsolete IPM-2 balance accuracy control ratio DHBM[i]/b[i], for minor ICs { 1e-3 }
    DNS,  // Standard surface density (nm-2) for calculating activity of surface species (12.05)
    XwMin,// Cutoff mole amount for elimination of water-solvent { 1e-9 }
    ScMin,// Cutoff mole amount for elimination of solid sorbent {1e-7}
    DcMin,// Cutoff mole amount for elimination of solution- or surface species { 1e-30 }
    PhMin,// Cutoff mole amount for elimination of  non-electrolyte solution phase with all its components { 1e-10 }
    ICmin,// Minimal effective ionic strength (molal), below which the activity coefficients for aqueous species are set to 1. { 3e-5 }
    EPS,  // Precision criterion of the SolveSimplex() procedure to obtain the AIA ( 1e-6 to 1e-14 ) { 1e-10 }
    IEPS, // Convergence parameter of SACT calculation in sorption/surface complexation models { 0.01 to 0.000001, default 0.001 }
    DKIN; // Tolerance on the amount of DC with two-side metastability constraints  { 1e-7 }
    char *tprn;       // internal
    void write(GemDataStream& oss);
    void read(GemDataStream& oss);
};

struct SPP_SETTING
{   // Base Parametres of SP
    char ver[TDBVERSION]; // Version & Copyright 64
    BASE_PARAM p; // Flags and thresholds for numeric modules
    char           // default codes of values
    DCpct[7],      // Default DCOMP flags and codes
    DCpdc[10],     // Default DCOMP class and units
    BCpc[7],       // Default COMPOS configuration
    REpct[7],      // Default REACDC flags and codes

    REpdc[7],      // Default REACDC class and units
    REpvc[9],      // Default REACDC configuration
    RPpdc[11],      // Default RTPARM flags and codes
    RPpvc[33],     // Default RTPARM configuration  reserved
    PHsol_t[7],    // Default PHASE model codes
    PHpvc[7],      // Default PHASE configuration
    MUpmv[11],     // Default RMULTS configuration
    TPpdc[9],      // Default class and units for MTPARM
    TPpvc[21],     // Default MTPARM configuration
    SYppc[11],     // Default class and flags for SYSTEM
    SYpvc[29],     // Default SYSTEM confifuration
    UTppc[11],     // Default DUTERM class and flags
    PEpsc[13],     // Default PROCES class and flags
    PEpvc[13],     // Default PROCES configuration
    GDcode[2][20], // Default names of screen and graphs in GTDEMO ????
    GDpsc[7],      // Default names of lines on GTDEMO graphs
    GDpcc[2][9],   // Default axis names for GTDEMO
    GDptc[7],      // Default GTDEMO configuration
    GDpgw[7],      // Default setup of graphs in GTDEMO
    SDrefKey[32],  // sdref key
    Reserv[50-32]    // (reserved) objects later
    ;
    // for RTPARM
    short NP,NT,  // Default N of points (RTPARM): P, T
    NV,       // reserved
    Mode,     // Default indexation mode RTPARM
    ResShort[5];
    float        // RTPARM
    Pi[3],    // Default interval for pressure
    Ti[3],    // Default interval for temperature, C
    Vi[3],    // Default interval for volume, cm3
    DRpst, DRtcst,   // Default Pr, Tr for DCOMP & REACDC
    lowPosNum, // MULTI Cutoff moles of DC (Ls set) { 1e-19 };
    logXw,     // log(1e-16)
    logYFk,    // log(1e-9)
    aqPar[5];  // b_gamma, a0, NeutPolicy, GamH2O, b_gam_T_dep for auto aq phase model
//    ResFloat;   // one parameter for auto gas/fluid phase

    void write(GemDataStream& oss);
    void read(GemDataStream& oss);
};

extern SPP_SETTING pa_;

struct CompItem
{
    int line;
    short delta; // RpnItemType
    CompItem(int nLine, short nDelta):
            line(nLine), delta(nDelta)
    {}

};

// Module TParam (RMULTS+MTPARM+SYSTEM+MULTY see SubModules)
class TProfil : public TCModule
{
    TRMults* rmults;
    TMTparm* mtparm;
    TSyst* syst;
    TMulti* multi;

    bool newRecord;
    // new 12.12.12 data to Project record
    char *internalBufer;  // text bufer for internal Project settings
                          // (now only for built-in TDB)

    // to compare with old Project
    bool comp_change_all;
    char (*SFold)[PH_RKLEN];// List of PHASE definition keys [0:Fi-1]             DB
    char (*SMold)[DC_RKLEN];// List of DC definition keys (DCOMP, REACDC) [0:L-1] DB
    char (*SAold)[BC_RKLEN];// List of COMPOS definition keys [0:La-1]            DB
    char (*SBold)[IC_RKLEN];// List of ICOMP record keys (stoichiometry basis)[0:N-1] DB
    short *Llold;// L1 vector, shows a number of DC included to each phase [0:Fi-1] DB
    short Nold,     // N of IC, incl. Zz (charge) and Vol (volume)
    Lold,       // L   - of DC - total for all phases
    Fiold,      // FI  - total number of phases
    Fisold,     // FIs - total number of multi-component phases
    Laold,      // La  - of references to COMPOS records
    Lsold;      // Ls  - total number of DC in multi-component phases

    // data to load SysEq <project>:G:z_cp_config:0:0:1:25:0
    bool  isSysEq;
    TCIntArray DCon;
    TCIntArray PHon;
    TCIntArray DCoff;
    TCIntArray PHoff;

protected:

    void InitFN( const char * prfName, const char* prfTemplate );
    void RenameFN( const char * prfName, const char* prfTemplate );
    bool GetFN( const char * prfName, bool show_dlg=true );
    void SetFN();
    bool rCopyFilterProfile( const char * prfName );

    void OpenProfileMode(
        const char* key, bool changeAqGas, bool addFile,  bool remakeRec );
    bool NewProfileMode( bool remakeRec, gstring& key_templ );
    bool NewProfileModeElements( bool remakeRec, gstring& key_templ );
    void CalcAllSystems(int makeDump);
    void SaveOldList();
    void DeleteOldList();
    void TestChangeProfile();
    void Push( TIArray<CompItem>& aList, int aLine,
               short aDelta, const char* dbKeywd, gstring aKey );
    void ICcompare( TIArray<CompItem>& aIComp);
    void COMPcompare( TIArray<CompItem>& aCompos);
    void DCcompare( TIArray<CompItem>& aList, int& i,int& j, int nI, int nJ);
    void PHcompare( TIArray<CompItem>& aPhase, TIArray<CompItem>& aDComp);

    // multi load
    short BAL_compare();

public:

    static TProfil* pm;

   bool userCancel;
   bool stepWise;
   bool calcFinished;
   const char * status;

    bool useAqPhase;
    bool useGasPhase;

// temporary !Use_mt_mode
    bool fStopCalc;

    //RMULTS* mup;
    //MTPARM *tpp;
    //SYSTEM *syp;
    //MULTI *pmp;
//TMulti *pmulti;

    SPP_SETTING pa;

    TProfil( int nrt );
    void InitSubModules();

    const char* GetName() const
    {
        return "Project";
    }

    void ods_link(int i=0);
    void dyn_set(int i=0);
    void dyn_kill(int i=0);
    void dyn_new(int i=0);
    void set_def(int i=0);
    void DeleteRecord( const char *key, bool errinNo=true );
    void MakeQuery();
    const char* GetHtml();

    // Setup one of 5 default IPM numerical settings
    void ChangeSettings(int nSettings);

    // work with Project
    bool initCalcMode( const char * profileKey );
    void loadSystat( const char *key=0 );
    void newSystat( int mode );
    void deriveSystat();
    void PackSystat();
//    double CalcEqstat( bool show_progress=true );
    double CalcEqstat(double &kdTime, const long kTimeStep = 0, const double kTime = 0.  );
    int  indDC( int );
    int  indPH( int );
    void  deleteAutoGenerated();
    void  systbcInput( QWidget* par, const char * p_key );
    int PhIndexforDC( int xdc, bool system );
    gstring PhNameforDC( int xdc, bool system );
    gstring PhNameforDC( int xdc, int& xph, bool system );
    TCStringArray DCNamesforPh( const char *PhName, bool system );
    void DCNamesforPh( int xph, bool system, vector<int>& xdc, vector<gstring>& dcnames);

    // Multi make functions
    void PMtest( const char *key );
    void CheckMtparam();

   void LoadFromMtparm( QWidget* par, DATACH *CSD, bool no_interpolat );
    void CalcBcc(); // Calc bulk composition
    void ShowPhaseWindow(QWidget* par, const char *objName, int nLine);
    //void ShowEqPhaseWindow();
    void ShowDBWindow( const char *objName, int nLine=0 );
    void Clear_XeA_XeD_Phm_BIun();

    // Proces make functions
    void ET_translate( int nOet, int nOpex, int JB, int JE, int jb, int je,
     tget_ndx *get_ndx = 0 )
     { multi->ET_translate( nOet, nOpex, JB, JE, jb, je, get_ndx); }
    void getNamesList( int nO, TCStringArray& lst )
     { multi->getNamesList(nO, lst); }
    double MolWeight( int N, double *ICaw, double *Smline )
     { 	return syst->MolWeight( N, ICaw, Smline ); }
    void SyTestSizes()
     { 	syst->SyTestSizes(); }

    inline double HelmholtzEnergy( double x )
    {
        return multi->HelmholtzEnergy(x);
    }
    inline double InternalEnergy( double TC, double P )
    {
        return multi->InternalEnergy( TC, P );
    }


   //test
   void outMulti( GemDataStream& ff, gstring& path  );
   // brief_mode - Do not write data items that contain only default values
   // with_comments -Write files with comments for all data entries ( in text mode)
   // addMui - Print internal indices in RMULTS to IPM file for reading into Gems back
   void outMulti( gstring& path, bool addMui, bool with_comments = true, bool brief_mode = false );
   void makeGEM2MTFiles(QWidget* par);
   void outMultiTxt( const char *path, bool append=false  );
   void readMulti( GemDataStream& ff );
   void readMulti( const char* path,  DATACH  *dCH );
   void CmReadMulti( QWidget* par, const char* path );
   double ComputeEquilibriumState( long int& NumPrecLoops, long int& NumIterFIA, long int& NumIterIPM );
   //long int testMulti( );
   bool CompareProjectName( const char* SysKey );
   void ChangeTPinKey( double T, double P );
   void SetSysSwitchesFromMulti( );
};

/* Work codes of surface site types in pm->AtNdx vector (compatibility with old-style SCMs *
enum SurTypeNdx {
    AT_SA0=0, AT_SB0=0, AT_SA1=1, AT_SB1=1, AT_SA2=2, AT_SB2=2, AT_SA3=3,
    AT_SB3=3, AT_SA4=4, AT_SB4=4, AT_SA5=5, AT_SB5=5,
    MSPN = 2   * Max number of EDL planes considered in old-style SCMs *
}; */


#endif // #ifndef IPMGEMPLUGIN


// Work DC classifier codes  pm->DCCW
enum SolDCodes {

    DC_SINGLE = 'U',        // This DC is a single-component phase
    DC_SYMMETRIC = 'I',     // This DC is in symmetric solution phase
    DC_ASYM_SPECIES = 'S',  // This is DC-solute(sorbate) in asymmetric phase
    DC_ASYM_CARRIER = 'W'   // This is carrier(solvent) DC in asymmetric phase
};


#endif  // _m_param_h


