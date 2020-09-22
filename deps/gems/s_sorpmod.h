//-------------------------------------------------------------------
// $Id: s_sorption.h 725 2012-10-02 15:43:37Z kulik $
//
// Stub declaration of new versions of SCMs (TSorpMod class)
//  Template: s_fgl.h (declarations of TSolMod class)
//
// Copyright (C) 2010-2012 D.Kulik, T.Wagner
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
#ifndef S_SORPTION_H
#define S_SORPTION_H

//#include "s_fgl.h"

// const int   MAXDCNAME = 16, MAXPHASENAME = 16, MST =   6,

const int MAXEILAYERS = 4;

struct SorptionSiteData {
    char  SiteT_;        // Site type code, see SITETYPECODES
    char  SACTC_;        // SACT equation code, see SACTCODES
    long int NSpec_;     // Number of surface species that can bind to this site (>=1; 0 if the site ignored)

    double  qC_;        //  site capacity parameter in mol/kg (CEC in eq/kg)
    double  GamC_;      //  site density parameter in mol/m2  (CEC in eq/m2)
    double  AlphaF_;    //  Frumkin interaction parameter;
    double  BETp_;      //  BET isotherm parameter p
    double  BETq_;      //  BET isotherm parameter q

    double  *arDent;    // Species denticity or coordination number on this site type [NSpecies]
    double  *arSpDUL;   // temporary upper constraint on species amount for SAT calculations [NSpecies]
    double  *ar_nxs;    // moles of sites taken by a surface species [NSpecies]
    long int *ar_xst;   // Index of surface species on surface tile (phase)  [NSpecies]
    double *arlnSACT;   // ln of SACT [NSpecies] - output (direct access, incremental) ?????
    double *arlnGamF;   // ln of Frumkin or BET term [NSpecies] - output (direct access, incremental) ????

    double *OccTot;     // Total amount of occupied sites (output scalar)
    double *FreTot;     // Total amount of free sites (output scalar)
};

class TSurfSiteMod
{
protected:

char  SiteT;   // Site type code, see SITETYPECODES
char  SACTC;   // SACT equation code, see SACTCODES

long int NSpec;   // Number of surface species on this site (at least 1; 0 if site to be ignored)

double  qC,          //  site capacity parameter in mol/kg (CEC in eq/kg)
        GamC,        //  site density parameter in mol/m2  (CEC in eq/m2)
        AlphaF,      //  Frumkin interaction parameter;
        BETp,        //  BET isotherm parameter
        BETq;        //  BET isotherm parameter

double  XSsI,        // total mole amount of surface species on this site
        XSsM,        // maximum (limiting) amount of sites
        OcTot,       // Total amount of occupied sites
        FrTot;       // Total amount of free sites

// double MASDJ[DFCN]; Parameters of surface species in surface complexation models
// enum {
//	//[0] - max site density in mkmol/(g sorbent); [1] - species charge allocated to 0 plane;
//	//[2] - surface species charge allocated to beta -or third plane; [3] - Frumkin interaction parameter;
//	//[4] species denticity or coordination number; [5]  - reserved parameter (e.g. species charge on 3rd EIL plane)]
//   XL_ST = 0, XL_EM, XL_SI, XL_SP
// };

   double *Dent;     //  Species denticity or coordination number [NSpec]
   double *SpDUL;    // temporary upper constraint on species amount for SAT calculations [NSpec]
   double *nxs;      // moles of surface species on this site (picked up) [NSpec]
   long int *xst;    // Index of surface species on surface tile (phase)  [NSpec]
// results
//   double (*D)[MST];  // Reserved; new work array for calc. surface act.coeff.
// double lnSAC[4];   // former lnSAT ln surface activity coeff and Coulomb's term  [Lads][4]
   double *ISAT;     // ISAT for each species (in SACT calculations) [NSpec]
   double *eF;       // Frumkin exponent acting on species [NSpec]
   double *SACT;     // Surface activity coefficient terms (config. entropy terms) [NSpec]
   double *lnSACT;   // ln of SACT on this site type [NSpec] - output
   double *lnGamF;   // ln activity coefficient due to Frumkin or BET isotherm - output

public:

    // Generic constructor
    TSurfSiteMod( SorptionSiteData ssd );

    // Destructor
    virtual ~TSurfSiteMod();

    // Other methods

};

struct SorptionData {

char Mod_Code_;    // Code of the sorption phase model - see SORPPHASECODES
char EIL_Code_;	  // Code for specific EIL model- see EILMODCODES  (before: SCMC)
char  PhasNam_[MAXPHASENAME];      // Phase name (for specific built-in models)

long int NSpec_;  // Total number of species assigned to this surface tile or Donnan phase
long int nlPh_;   // number of linked phases (cf. lPh), default 1 (the sorbent phase)
long int nlPhC_;  // number of linked phase parameter coefficient per link (default 0)

long int nEIml_; // number of EIL model layers (default 0, maximum MAXEILAYERS=4);
long int nEIpl_; // number of EIL params per layer (default 0);
long int nCDcf_; // number of CD coefs per DC (default 0);
long int nEIres_; // reserved

long int NsiteTs_;  // Number of surface site types per this surface patch (min 1 max 6), for Donnan 1
                    //   (if 0 then this sorption phase model is ignored)
long int nISTcf_;  // number of isotherm coeffs per site (default 1)
long int nISDcf_;  // number of isotherm coeffs per DC (default 1)
long int DelMax_;  // max.denticity of DC (default 1)

long int kSorPh_;  // Index of the sorbent phase in GEM IPM work structure (MULTI)
                   // if -1 then this is a site-balance based approach; SSa_ or sVp_ must be provided here

long int *xsMd_;    // denticity of surface species per surface site (site allocation) [NSpec][NsiteTs]

// long int *arPsDiS_;  // array of DC denticity and indexes of binding sites [NSpec*(DelMax+1)]

char *DC_Cods_;    // DC class codes for species [NSpec_]
char *IsoCt_;      // isotherm and SATC codes for surface site types [2*NsiteTs_]

long int *arPhLin_;  // indexes of linked (sorbent) phase(s) and link type code(s) [nlPh*2] read-only

    double T_k_;         // Temperature, K (initial)
    double P_bar_;       // Pressure, bar (initial)
    double IS_;          // Effective molal ionic strength of aqueous electrolyte
    double pH_;         // pH of aqueous solution
    double pe_;         // pe of aqueous solution
    double Eh_;         // Eh of aqueous solution, V

// This is taken over from the aSorMc[16] piece from MULTI
    double *arSorMc;
//    double sSA_,  // Specific surface area, m2/g, default: 0.
//           sgw_, // Standard mean surface energy of solid-aqueous interface, J/m2
//           sgg_, // Standard mean surface energy of gas-aqueous interface, J/m2
//           rX0_,    // Mean radius r0 for (spherical or cylindrical) particles, nm (reserved)
//           hX0_,    // Mean thickness h0 for cylindrical or 0 for spherical particles, nm (reserved)
    //
//           sVp_,  // Specific pore volume of phase, m3/g (default: 0)
//           frSA_,   // reactive fraction of surface area (def. 1)

//           nPh_,  // current amount of this phase, mol (read-only)
//           mPh_,  // current mass of this phase, g (read-only)
//           vPh_,  // current volume of this phase, cm3 (read-only)
//           sAPh_,  // current surface of this phase, m2
//           OmPh_,  // phase stability index (lg scale), input
        //
//           *nPul_, // pointer to upper restriction to this phase amount, mol (calculated here)
//           *nPll_, // pointer to lower restriction to this phase amount, mol (calculated here)
//           *sGP_;  // pointer to surface free energy of the phase, J (YOF*PhM)

double N_Cst_;    // Standard surface number density, 1/nm2
double G_Cst_;    // Standard surface density, mol/m2
double q_Cst_;    // Standard sorption capacity, mol/kg(sorbent) or eq/kg(sorbent)

    double Nfsp_,     // Fraction of the sorbent specific surface area or volume allocated to surface type (>0 <10000)
           MASDT_,    // Total sorption capacity for this surface type (mol/kg), before was (mkmol/g)
           XetaC_,    // Total permanent charge capacity CEC, mol/kg
           VetaP_,    // Total permanent volume charge density (eq/m3)
           XetaP_,    // Total density of permanent charge (eq/m2), before mkeq/m2
           ValP_,     // Porosity of the Donnan sorbent (d/less)
           ParD1_,    // Donnan model parameter 1
           ParD2_,    // Donnan model parameter 2
           XlamA_;    // Factor of EDL discretness  A < 1, reserved

    double *Xcap_,    // Capacitance density of EIL layers, F/m2 [MAXEILAYERS]
           *Xdl_;     // Effective thickness of EIL layers, nm, reserved
    double (*CD)[MAXEILAYERS];    // Species charges allocated to EIL planes [NspT]

double *arlPhc;  // array of phase link parameters (sum(LsPhl[k][1] over Fi)

double *arEImc;  // EIL model coefficients table [nEIml_*nEIpl]
double *armCDc;  // CD EIL model coefficients table [NSpec_*nCDcf_]
double *IsoPc;   // Isotherm equation coefficients table [NSpec_*nISDcf_]
double *IsoSc;   // Isotherm equation coeffs per surface site type [NsiteTs*nISTcf_]

double *arPparc;    // Partial pressures -> NSpecies
double *arWx;       // Species mole fractions ->NSpec_ read-only
double *arnx;     // Pointer to mole amounts of phase components (provided) [NSpec_] read-only
// output
double *arlnScalT;   // Surface/volume scaling activity correction terms [NSpec_]  direct access
double *arlnSACT;    // Pointer to ln SACT for surface species [NSpec_]  direct access
double *arlnGammaF;  // Pointer to Frumkin or BET non-electrostatic activity coefficients [NSpec_] direct access
double *arCTerms;    // Pointer to Coulombic correction terms (electrostatic activity coefficients) [NSpec_] direct access
double *arlnGamma;   // Pointer to ln activity coefficients of sorption phase components mixing [NSpec_] direct access
                   // (memory under pointers must be provided from the calling program)

double *arVol;      // molar volumes of end-members (species) cm3/mol [NSpec_] read-inly
// char  (*arSM)[MAXDCNAME];  // pointer to the list of DC names in the phase [NSpec_] read-only
char  *arDCC;   // pointer to the classifier of DCs involved in this phase [NSpec_] read-only

};

class TSorpMod
{     // Treatment of surface tile (patch) or Donnan volume phases in sorption models
protected:

    char Mod_Code;    // Code of the sorption phase model - see SORPPHASECODES
    char EIL_Code;	  // Code for specific EIL model- see EILMODCODES  (before: SCMC)
    char  PhasNam_[MAXPHASENAME];      // Phase name (for specific built-in models)

    long int NSpec;  // Total number of species assigned to this surface tile or Donnan phase
    long int nlPh;   // number of linked phases (cf. lPh), default 1 (the sorbent phase)
    long int nlPhC;  // number of linked phase parameter coefficient per link (default 0)

    long int nEIml;  // number of EIL model layers (default 0);
    long int nEIpl;  // number of EIL params per layer (default 0);
    long int nCDcf;  // number of CD coefs per DC (default 0);
    long int nEIres; // reserved

    long int NsiteTs; // Number of surface site types per surface patch type (min 1 max 6), for Donnan 1
                      //   (if 0 then this sorption phase model is ignored)
    long int nISTcf;  // number of isotherm coeffs per site (default 1)
    long int nISDcf;  // number of isotherm coeffs per DC (default 1)
    long int DelMax;  // max.denticity of DC (default 1)

    long int kSorPh;  // Index of the sorbent phase in GEM IPM work structure (MULTI)
                      // if -1 then this is a site-balance based approach; SSa_ or sVp_ must be provided

    long int **xsMd;      // denticity of surface species per surface site (site allocation) [NSpec][NsiteTs]
//    long int **PsDS;    // array of DC denticity and indexes of binding sites [NSpec][DelMax+1]
    long int (*PhLin)[2]; // indexes of linked (sorbent) phase(s) and link type code(s) [nlPh][2] read-only
    char  *DCCs;         // DC class codes for species [NSpec_]
    char (*IsoCt)[2];     // isotherm and SATC codes for surface site types NsiteTs][2]

    double T_k;        // Temperature, K (initial)
    double P_bar;      // Pressure, bar (initial)
    double IS;         // Effective molal ionic strength of aqueous electrolyte
    double pH;         // pH of aqueous solution
    double pe;         // pe of aqueous solution
    double Eh;         // Eh of aqueous solution, V

// This is taken over from the aSorMc[16] piece from MULTI
    double sSA,  // Specific surface area, m2/g, default: 0.
           sgw, // Standard mean surface energy of solid-aqueous interface, J/m2
           sgg, // Standard mean surface energy of gas-aqueous interface, J/m2
//           rX0,    // Mean radius r0 for (spherical or cylindrical) particles, nm (reserved)
//           hX0,    // Mean thickness h0 for cylindrical or 0 for spherical particles, nm (reserved)
    //
           sVp,  // Specific pore volume of phase, m3/g (default: 0)
           frSA,   // reactive fraction of surface area (def. 1)

//           nPh,  // current amount of this phase, mol (read-only)
//           mPh,  // current mass of this phase, g (read-only)
//           vPh,  // current volume of this phase, cm3 (read-only)
//           sAPh,  // current surface of this phase, m2
//           OmPh,  // phase stability index (lg scale), input
        //
           sGP_;  // surface free energy of the phase, J (YOF*PhM)

double N_Cst;    // Standard surface number density, 1/nm2
double G_Cst;    // Standard surface density, mol/m2
double q_Cst;    // Standard sorption capacity, mol/kg(sorbent) or eq/kg(sorbent)

    double R_CONST;   // R constant
    double F_CONST;   // F (Faraday's) constant

    // model parameters
    double Nfsp,     // Fraction of the sorbent specific surface area or volume allocated to surface type (>0 <10000)
           MASDT,    // Total sorption capacity for this surface type (mol/kg), before was (mkmol/g)
           XetaC,    // Total permanent charge capacity CEC, mol/kg
           VetaP,    // Total permanent volume charge density (eq/m3)
           XetaP,    // Total density of surface permanent charge (eq/m2), before mkeq/m2
           ValP,     // Porosity of the Donnan sorbent (d/less)
           ParD1,    // Donnan model parameter 1
           ParD2,    // Donnan model parameter 2
           *Xcap,    // Capacitance density of EIL layers 0 1 2 3 ... , F/m2 [nEIml]
           *Xdl,     // Effective thicknesses of EIL layers, nm, reserved
           XdlD,     // Effective thickness of diffuse layer, nm, reserved
           XlamA;    // Factor of EDL discretness  A < 1, reserved
    double (*CD)[MAXEILAYERS];    // Species charges allocated to 0, 1 and 2 planes [NspT]

    TSurfSiteMod* sitMod; // Pointer to array of TSurfSiteMod instances - [NsiT]

 // current values
    double XsTs,     // Total number of moles of species in this tile or Donnan volume
           XaTs,     // Total moles of 'solvent' (e.g. >OH group or H2O in Donnan phase)
           Sarea,    // Current area occupied by this surface tile, (m2)
           Volum,    // Current volume (if this is Donnan electrolyte), (m3)
         Xeta[MAXEILAYERS],    // Total charge of surface species on EIL planes 0,1,2, ..., equiv
         Xpsi[MAXEILAYERS],    // Relative potential at EIL planes, V
           VetaD,    // Total charge in the Donnan volume, moles
           VPsiD,    // Relative Donnan potential, V
           ISD;      // Ionic strength (for Donnan electrolyte)

    double *nx;     // pointer to moles of surface species on this surface tile (read-only) [Nspec]

    double **lnSACTs;  // Work array of ln SACT for surface species on sites [Nspec][NsiteTs]
    double *(XTS[2]);  // Total number of moles of surface DC and 'solvent' DC at surface site [2][NsiteTs]

    // current results
    double Gex, Hex, Sex, CPex, Vex, Aex, Uex;   // molar electrostatic excess properties for surface species
//    double Gid, Hid, Sid, CPid, Vid, Aid, Uid;   // molar ideal mixing properties for surface species

    double *lnScalT;   // Surface/volume scaling activity correction terms [Nspec]
    double *lnSACT;    // Pointer to ln SACT for surface species [Nspec]
    double *lnGammaF;  // Pointer to Frumkin or BET non-electrostatic activity coefficients [Nspec]
    double *CTerms;    // Pointer to Coulombic correction terms (electrostatic activity coefficients) [Nspec]
    double *lnGamma;   // Pointer to ln activity coefficients of sorption phase components mixing [Nspec]
                       // (memory under pointers must be provided from the calling program)

     public:
    // Generic constructor
    TSorpMod( SorptionData sds  );

    // Destructor
    virtual ~TSorpMod();


    virtual long int SorptionSpecies()
    {
            return 0;
    };

    virtual long int PTparam()
    {
            return 0;
    };

    virtual long int IsothermMod()
    {
            return 0;
    };

    virtual long int ElstatMod()
    {
            return 0;
    };
/*
    virtual long int ExcessProp( double *Zex )
    {
            return 0;
    };

    virtual long int IsothermProp( double *Zid )
    {
            return 0;
    };
*/
    // set new system state
    long int UpdatePT ( double T_k, double P_bar );

    bool testSizes( long int NSpecies, long int NSurSpecies, long int NSorbentEMs,
                    long int NSurfTypes, char Mod_Code, char EIL_Code );

    // getting phase name
    void GetPhaseName( const char *PhName );

};


class TNEMcalc: public TSorpMod  // Non-electrostatic sorption phase model without site balances
{
    protected:

        double Xcond, 	// conductivity of phase carrier, sm/m2, reserved
               Xeps,  	// rel. diel. permeability of phase carrier (solvent), reserved
               Aalp,  	// Specific surface area of phase carrier (sorbent) (m2/g) !!! m2/kg
               Sigw,  	// Specific surface free energy for phase-water interface (J/m2)
               Sigg,  	// Specific surface free energy for phase-gas interface (J/m2) (not yet used)
               Xr0,   // Mean radius r0 for (spherical or cylindrical) particles, nm (reserved) Xr0h0[2]
               Xh0;   // Mean thickness h0 for cylindrical or 0 for spherical particles, nm (reserved)

//        TSurfPatchMod* patchMod[]; // Pointer to array of TSurfPatch instances - size NsurT

        // internal functions
        void alloc_internal();
        void free_internal();


public:
        // Constructor
        TNEMcalc( long int NCmp, double Pp, double Tkp );
        TNEMcalc( long int NSpecies, long int NSurSpecies, long int NSorbentEMs, long int NSurfTypes,
                  char Mod_Code, char *Phase_Name, char **SM3_, char *DCC_, char **SCM_, char **SATT_,
                  double Aalp_, double *arnx_, double *arlnGam, double *arlnSACT, double T_k, double P_bar  );
        // Destructor
        ~TNEMcalc();

        // Calculates pure species properties (pure fugacities)
        long int PureSpecies();

        // Calculates T,P corrected interaction parameters
        long int PTparam();

        // Calculates activity coefficients
        long int MixMod();

        // calculates excess properties
        long int ExcessProp( double *Zex );

        // calculates ideal mixing properties
        long int IdealProp( double *Zid );

        // Calculates pure species properties (called from DCthermo)
        long int SRKCalcFugPure( double Tmin, float *Cpg, double *FugProps );

};


#endif // S_SORPTION_H
