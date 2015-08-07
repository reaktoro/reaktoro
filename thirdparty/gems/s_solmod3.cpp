//-------------------------------------------------------------------
// $Id: s_solmod3.cpp 724 2012-10-02 14:25:25Z kulik $
//
/// \file s_solmod3.cpp
/// Implementation of TSolMod derived classes
/// for activity models of mixing in condensed (solid and liquid) phases
///  (TVanLaar, TRegular, TRedlichKister, TNRTL, TWilson, TMargulesTernary,
///  TMargulesBinary, TGuggenheim, TIdeal multi-site, TBerman, TCEFmod)
//
// Copyright (c) 2007-2014  T.Wagner, D.Kulik, S.Dmitrieva
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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------

#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include "s_solmod.h"


//=============================================================================================
// Implementation of ideal mixing model for multicomponent solid solutions
// References: Price 1989
// also used with scripted models to provide ideal mixing term in the multi-site case
// (c) DK/TW November 2010
//=============================================================================================


// Generic constructor for the TIdeal class
TIdeal::TIdeal( SolutionData *sd ):
                TSolMod( sd )
{
}

TIdeal::~TIdeal()
{
}

long int TIdeal::PTparam()
{
   return 0;
}


/// Calculates ideal configurational terms in case of multi-site mixing
/// to preserve values computed in Phase scripts.
/// Only increments lnGamma[j] - may need to be cleaned before running MixMod
long int TIdeal::MixMod()
{
   long int retCode, j;
   retCode = IdealMixing();

   if(!retCode)
   {
      for(j=0; j<NComp; j++)
          lnGamma[j] += lnGamConf[j];
   }
   return 0;
}


/// calculates bulk phase excess properties
long int TIdeal::ExcessProp( double *Zex )
{

        // assignments (excess properties)
        Zex[0] = 0.;
        Zex[1] = 0.;
        Zex[2] = 0.;
        Zex[3] = 0.;
        Zex[4] = 0.;
        Zex[5] = 0.;
        Zex[6] = 0.;

        return 0;
}


/// calculates ideal mixing properties
long int TIdeal::IdealProp( double *Zid )
{
        Hid = 0.0;
        CPid = 0.0;
        Vid = 0.0;
        Sid = ideal_conf_entropy();
        Gid = Hid - Sid*Tk;
        Aid = Gid - Vid*Pbar;
        Uid = Hid - Vid*Pbar;

        // assignments (ideal mixing properties)
        Zid[0] = Gid;
        Zid[1] = Hid;
        Zid[2] = Sid;
        Zid[3] = CPid;
        Zid[4] = Vid;
        Zid[5] = Aid;
        Zid[6] = Uid;

        return 0;
}




//=============================================================================================
// Van Laar model for solid solutions
// References: Holland and Powell (2003)
// (c) TW March 2007, added sublattice ideal TW-DK in Dec 2011
//=============================================================================================


// Generic constructor for the TVanLaar class
TVanLaar::TVanLaar( SolutionData *sd ):
                TSolMod( sd )
{
    alloc_internal();
}


TVanLaar::~TVanLaar()
{
    free_internal();
}


void TVanLaar::alloc_internal()
{
	Wu = new double [NPar];
	Ws = new double [NPar];
	Wv = new double [NPar];
	Wpt = new double [NPar];
	Phi = new double [NComp];
	PsVol = new double [NComp];
}


void TVanLaar::free_internal()
{
	if(Wu) delete[]Wu;
	if(Ws) delete[]Ws;
	if(Wv) delete[]Wv;
	if(Wpt) delete[]Wpt;
	if(Phi) delete[]Phi;
	if(PsVol) delete[]PsVol;
}


/// Calculates T,P corrected binary interaction parameters
long int TVanLaar::PTparam()
{
	long int j, ip;

    if ( NPcoef < 3 || NPar < 1 )
       return 1;

    for (j=0; j<NComp; j++)
    {
    	PsVol[j] = aDCc[NP_DC*j];  // reading pseudo-volumes
    }

    for (ip=0; ip<NPar; ip++)
	{
           Wu[ip] = aIPc[NPcoef*ip];
           Ws[ip] = aIPc[NPcoef*ip+1];
           Wv[ip] = aIPc[NPcoef*ip+2];
           Wpt[ip] = Wu[ip]+ Ws[ip]*Tk + Wv[ip]*Pbar;
           aIP[ip] = Wpt[ip];
	}
    return 0;
}


/// Calculates activity coefficients
long int TVanLaar::MixMod()
{
	long int ip, j, i1, i2;
        double dj, dk, sumPhi, lnGamRT, lnGam;

	if ( NPcoef < 3 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
		return 1;

    // Trying sublattice ideal mixing model
    long int retCode;
    retCode = IdealMixing();
    if(!retCode)
    {
       for(j=0; j<NComp; j++)
           lnGamma[j] += lnGamConf[j];
    }

	// calculating Phi values
	sumPhi = 0.;
	for (j=0; j<NComp; j++)
	{
		sumPhi +=  x[j]*PsVol[j];
	}

	if( fabs(sumPhi) < 1e-30 )
		return 2;    // to prevent zerodivide

	for (j=0; j<NComp; j++)
		Phi[j] = x[j]*PsVol[j]/sumPhi;

	// calculate activity coefficients
	for (j=0; j<NComp; j++)
	{
		lnGamRT = 0.;
		for (ip=0; ip<NPar; ip++)  // inter.parameters indexed with ip
		{
			i1 = aIPx[MaxOrd*ip];
			i2 = aIPx[MaxOrd*ip+1];

			if( j == i1 )
				dj = 1.;
			else
				dj = 0.;
			if( j == i2 )
				dk = 1.;
			else
				dk = 0.;
			lnGamRT -= (dj-Phi[i1])*(dk-Phi[i2])*Wpt[ip]
			             *2.*PsVol[j]/(PsVol[i1]+PsVol[i2]);
		}
		lnGam = lnGamRT/(R_CONST*Tk);
        lnGamma[j] += lnGam;
	}
	return 0;
}


/// calculates bulk phase excess properties
long int TVanLaar::ExcessProp( double *Zex )
{
	long int ip, j, i1, i2;
	double sumPhi, g, v, s, u;

	if ( NPcoef < 3 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
		return 1;

	// calculating Phi values
	sumPhi = 0.;
	for (j=0; j<NComp; j++)
	{
		PsVol[j] = aDCc[NP_DC*j];  // reading pseudo-volumes
	    sumPhi +=  x[j]*PsVol[j];
	}

	if( fabs(sumPhi) < 1e-30 )
		return 2;    // to prevent zerodivide!

	for (j=0; j<NComp; j++)
	    Phi[j] = x[j]*PsVol[j]/sumPhi;

	// calculate bulk phase excess properties
	g = 0.0; s = 0.0; v = 0.0; u = 0.0;

	for (ip=0; ip<NPar; ip++)
	{
		i1 = aIPx[MaxOrd*ip];
	    i2 = aIPx[MaxOrd*ip+1];
	    g += Phi[i1]*Phi[i2]*2.*sumPhi/(PsVol[i1]+PsVol[i2])*Wpt[ip];
	    v += Phi[i1]*Phi[i2]*2.*sumPhi/(PsVol[i1]+PsVol[i2])*Wv[ip];
	    u += Phi[i1]*Phi[i2]*2.*sumPhi/(PsVol[i1]+PsVol[i2])*Wu[ip];
	    s -= Phi[i1]*Phi[i2]*2.*sumPhi/(PsVol[i1]+PsVol[i2])*Ws[ip];
	 }

	 Gex = g;
	 Sex = s;
	 CPex = 0.0;
	 Vex = v;
	 Uex = u;
	 Hex = Uex + Vex*Pbar;
	 Aex = Gex - Vex*Pbar;
	 Uex = Hex - Vex*Pbar;

	 // assignments (excess properties)
	 Zex[0] = Gex;
	 Zex[1] = Hex;
	 Zex[2] = Sex;
	 Zex[3] = CPex;
	 Zex[4] = Vex;
	 Zex[5] = Aex;
	 Zex[6] = Uex;

	 return 0;
}


/// calculates ideal mixing properties
long int TVanLaar::IdealProp( double *Zid )
{
//	long int j;
//	double si;
//	si = 0.0;
//	for (j=0; j<NComp; j++)
//	{
//		if ( x[j] > 1.0e-32 )
//			si += x[j]*log(x[j]);
//	}
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
    Sid = ideal_conf_entropy();
//	Sid = (-1.)*R_CONST*si;
	Gid = Hid - Sid*Tk;
	Aid = Gid - Vid*Pbar;
	Uid = Hid - Vid*Pbar;

	// assignments (ideal mixing properties)
	Zid[0] = Gid;
	Zid[1] = Hid;
	Zid[2] = Sid;
	Zid[3] = CPid;
	Zid[4] = Vid;
	Zid[5] = Aid;
	Zid[6] = Uid;

	return 0;
}


//=============================================================================================
// Regular model for multicomponent solid solutions
// References:  Holland and Powell (1993)
// (c) TW March 2007, added sublattice ideal TW-DK in Dec 2011
//=============================================================================================


// Generic constructor for the TRegular class
TRegular::TRegular( SolutionData *sd ):
                TSolMod( sd )
{
    alloc_internal();
}


TRegular::~TRegular()
{
    free_internal();
}


void TRegular::alloc_internal()
{
        Wu = new double [NPar];
        Ws = new double [NPar];
        Wv = new double [NPar];
        Wpt = new double [NPar];
}


void TRegular::free_internal()
{
        if(Wu) delete[]Wu;
        if(Ws) delete[]Ws;
        if(Wv) delete[]Wv;
        if(Wpt) delete[]Wpt;
}


/// Calculates T,P corrected binary interaction parameters
long int TRegular::PTparam()
{
        long int ip;

        if ( NPcoef < 3 || NPar < 1 )
                   return 1;

        for (ip=0; ip<NPar; ip++)
        {
            Wu[ip] = aIPc[NPcoef*ip];
            Ws[ip] = aIPc[NPcoef*ip+1];
            Wv[ip] = aIPc[NPcoef*ip+2];
            Wpt[ip] = Wu[ip]+ Ws[ip]*Tk + Wv[ip]*Pbar;
            aIP[ip] = Wpt[ip];
        }
        return 0;
}


/// Calculates activity coefficients
long int TRegular::MixMod()
{
        long int ip, j, i1, i2;
        double dj, dk, lnGamRT, lnGam;

        if ( NPcoef < 3 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
                return 1;

        // Trying sublattice ideal mixing model
        long int retCode;
        retCode = IdealMixing();
        if(!retCode)
        {
           for(j=0; j<NComp; j++)
               lnGamma[j] += lnGamConf[j];
        }

        // calculate activity coefficients
        for (j=0; j<NComp; j++)
        {
                lnGamRT = 0.;

                for (ip=0; ip<NPar; ip++)  // inter.parameters indexed with ip
                {
                        i1 = aIPx[MaxOrd*ip];
                        i2 = aIPx[MaxOrd*ip+1];

                        if( j == i1 )
                                dj = 1.;
                        else
                                dj = 0.;
                        if( j == i2 )
                                dk = 1.;
                        else
                                dk = 0.;
                        lnGamRT -= (dj-x[i1])*(dk-x[i2])*Wpt[ip];
                }

                lnGam = lnGamRT/(R_CONST*Tk);
                lnGamma[j] += lnGam;
        }
        return 0;
}


/// calculates bulk phase excess properties
long int TRegular::ExcessProp( double *Zex )
{
	long int ip, i1, i2;
	double g, v, s, u;

	if ( NPcoef < 3 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
		return 1;

	// calculate bulk phase excess properties
	g = 0.0; s = 0.0; v = 0.0; u = 0.0;

	for (ip=0; ip<NPar; ip++)
	{
		i1 = aIPx[MaxOrd*ip];
		i2 = aIPx[MaxOrd*ip+1];
		g += x[i1]*x[i2]*Wpt[ip];
		v += x[i1]*x[i2]*Wv[ip];
		u += x[i1]*x[i2]*Wu[ip];
		s -= x[i1]*x[i2]*Ws[ip];
	}

	Gex = g;
	Sex = s;
	CPex = 0.0;
	Vex = v;
	Uex = u;
	Hex = Uex + Vex*Pbar;
	Aex = Gex - Vex*Pbar;
	Uex = Hex - Vex*Pbar;

	// assignments (excess properties)
	Zex[0] = Gex;
	Zex[1] = Hex;
	Zex[2] = Sex;
	Zex[3] = CPex;
	Zex[4] = Vex;
	Zex[5] = Aex;
	Zex[6] = Uex;

	return 0;
}


/// calculates ideal mixing properties
long int TRegular::IdealProp( double *Zid )
{
//	long int j;
//	double si;
//	si = 0.0;
//	for (j=0; j<NComp; j++)
//	{
//		if ( x[j] > 1.0e-32 )
//			si += x[j]*log(x[j]);
//	}
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
    Sid = ideal_conf_entropy();
//	Sid = (-1.)*R_CONST*si;
	Gid = Hid - Sid*Tk;
	Aid = Gid - Vid*Pbar;
	Uid = Hid - Vid*Pbar;

	// assignments (ideal mixing properties)
	Zid[0] = Gid;
	Zid[1] = Hid;
	Zid[2] = Sid;
	Zid[3] = CPid;
	Zid[4] = Vid;
	Zid[5] = Aid;
	Zid[6] = Uid;

	return 0;
}




//=============================================================================================
// Redlich-Kister model for multicomponent solid solutions
// References: Hillert (1998)
// (c) TW March 2007
//=============================================================================================


// Generic constructor for the TRedlichKister class
TRedlichKister::TRedlichKister( SolutionData *sd ):
                TSolMod( sd )
{
    alloc_internal();
}


TRedlichKister::~TRedlichKister()
{
    free_internal();
}


void TRedlichKister::alloc_internal()
{
	Lu = new double [NPar][4];
	Ls = new double [NPar][4];
	Lcp = new double [NPar][4];
	Lv = new double [NPar][4];
	Lpt = new double [NPar][4];
}


void TRedlichKister::free_internal()
{
	if(Lu) delete[]Lu;
	if(Ls) delete[]Ls;
	if(Lv) delete[]Lv;
	if(Lpt) delete[]Lpt;
	if(Lcp) delete[]Lcp;
}


///   Calculates T,P corrected binary interaction parameters
long int TRedlichKister::PTparam()
{
	long int ip;

	if ( NPcoef < 16 || NPar < 1 )
		return 1;

	// read in interaction parameters
	for (ip=0; ip<NPar; ip++)
	{
		Lu[ip][0] = aIPc[NPcoef*ip+0];
	   	Ls[ip][0] = aIPc[NPcoef*ip+1];
	   	Lcp[ip][0] = aIPc[NPcoef*ip+2];
	   	Lv[ip][0] = aIPc[NPcoef*ip+3];
	   	Lpt[ip][0] = Lu[ip][0] + Ls[ip][0]*Tk + Lcp[ip][0]*Tk*log(Tk) + Lv[ip][0]*Pbar;
                aIP[ip] = Lpt[ip][0];

	   	Lu[ip][1] = aIPc[NPcoef*ip+4];
	   	Ls[ip][1] = aIPc[NPcoef*ip+5];
	   	Lcp[ip][1] = aIPc[NPcoef*ip+6];
	   	Lv[ip][1] = aIPc[NPcoef*ip+7];
	   	Lpt[ip][1] = Lu[ip][1] + Ls[ip][1]*Tk + Lcp[ip][1]*Tk*log(Tk) + Lv[ip][1]*Pbar;

	   	Lu[ip][2] = aIPc[NPcoef*ip+8];
	   	Ls[ip][2] = aIPc[NPcoef*ip+9];
	   	Lcp[ip][2] = aIPc[NPcoef*ip+10];
	   	Lv[ip][2] = aIPc[NPcoef*ip+11];
	   	Lpt[ip][2] = Lu[ip][2] + Ls[ip][2]*Tk + Lcp[ip][2]*Tk*log(Tk) + Lv[ip][2]*Pbar;

	   	Lu[ip][3] = aIPc[NPcoef*ip+12];
	   	Ls[ip][3] = aIPc[NPcoef*ip+13];
	   	Lcp[ip][3] = aIPc[NPcoef*ip+14];
	   	Lv[ip][3] = aIPc[NPcoef*ip+15];
	   	Lpt[ip][3] = Lu[ip][3] + Ls[ip][3]*Tk + Lcp[ip][3]*Tk*log(Tk) + Lv[ip][3]*Pbar;

	}
	return 0;
}


/// Calculates activity coefficients
long int TRedlichKister::MixMod()
{
	long int ip, j, i1, i2, L, I, J;
        double L0, L1, L2, L3, lnGamRT, lnGam;

	if ( NPcoef < 16 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
		return 1;

	// loop over species
	for (j=0; j<NComp; j++)
	{
		lnGamRT = 0.;
		for (ip=0; ip<NPar; ip++)  // inter.parameters indexed with ip
		{
			i1 = aIPx[MaxOrd*ip];
			i2 = aIPx[MaxOrd*ip+1];

			if ( j == i1 || j == i2) // interaction terms with j
			{
				if ( i1 == j ) // check order of idexes
				{
					L = i1;
					I = i2;
					L0 = Lpt[ip][0];
					L1 = Lpt[ip][1];
					L2 = Lpt[ip][2];
					L3 = Lpt[ip][3];
				}
				else
				{
					L = i2;
					I = i1;
					L0 = Lpt[ip][0];
					L1 = -Lpt[ip][1];
					L2 = Lpt[ip][2];
					L3 = -Lpt[ip][3];
				}

				lnGamRT += L0*x[I]*(1.-x[L])
					+ L1*x[I]*(2.*(1.-x[L])*(x[L]-x[I])+x[I])
					+ L2*x[I]*(x[L]-x[I])*(3.*(1.-x[L])*(x[L]-x[I])+2.*x[I])
					+ L3*x[I]*pow((x[L]-x[I]),2.)*(4.*(1.-x[L])*(x[L]-x[I])+3.*x[I]);
			}

			else // interaction terms without j
			{
				I = i1;
				J = i2;
				L0 = Lpt[ip][0];
				L1 = Lpt[ip][1];
				L2 = Lpt[ip][2];
				L3 = Lpt[ip][3];

				lnGamRT -= x[I]*x[J]*( L0 + L1*2.*(x[I]-x[J])
					+ L2*3.*pow((x[I]-x[J]),2.)
					+ L3*4.*pow((x[I]-x[J]),3.) );
			}
		}

		lnGam = lnGamRT/(R_CONST*Tk);
		lnGamma[j] = lnGam;
	} // j

   	return 0;
}


/// calculates bulk phase excess properties
long int TRedlichKister::ExcessProp( double *Zex )
{
	long int ip, i1, i2;
	double LU, LS, LCP, LV, LPT, g, v, s, cp, u;

	if ( NPcoef < 16 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
		return 1;

   	// calculate bulk phase excess properties
   	g = 0.0; s = 0.0; cp = 0.0; v = 0.0; u = 0.0;

   	for (ip=0; ip<NPar; ip++)
   	{
   		i1 = aIPx[MaxOrd*ip];
   	   	i2 = aIPx[MaxOrd*ip+1];

      	LPT = Lpt[ip][0] + Lpt[ip][1]*(x[i1]-x[i2])
      			+ Lpt[ip][2]*pow((x[i1]-x[i2]),2.)
      			+ Lpt[ip][3]*pow((x[i1]-x[i2]),3.);

      	LV = Lv[ip][0] + Lv[ip][1]*(x[i1]-x[i2])
      			+ Lv[ip][2]*pow((x[i1]-x[i2]),2.)
      			+ Lv[ip][3]*pow((x[i1]-x[i2]),3.);

   	   	LU = (Lu[ip][0]-Lcp[ip][0]*Tk)
   	  			+ (Lu[ip][1]-Lcp[ip][1]*Tk)*(x[i1]-x[i2])
      			+ (Lu[ip][2]-Lcp[ip][2]*Tk)*pow((x[i1]-x[i2]),2.)
      			+ (Lu[ip][3]-Lcp[ip][3]*Tk)*pow((x[i1]-x[i2]),3.);

   	   	LS = (-Ls[ip][0]-Lcp[ip][0]*(1.+log(Tk)))
      			+ (-Ls[ip][1]-Lcp[ip][1]*(1.+log(Tk)))*(x[i1]-x[i2])
      			+ (-Ls[ip][2]-Lcp[ip][2]*(1.+log(Tk)))*pow((x[i1]-x[i2]),2.)
      			+ (-Ls[ip][3]-Lcp[ip][3]*(1.+log(Tk)))*pow((x[i1]-x[i2]),3.);

   	   	LCP = (-Lcp[ip][0]) + (-Lcp[ip][1])*(x[i1]-x[i2])
      			+ (-Lcp[ip][2])*pow((x[i1]-x[i2]),2.)
      			+ (-Lcp[ip][3])*pow((x[i1]-x[i2]),3.);

      	g += x[i1]*x[i2]*LPT;
      	v += x[i1]*x[i2]*LV;
      	u += x[i1]*x[i2]*LU;
      	s += x[i1]*x[i2]*LS;
      	cp += x[i1]*x[i2]*LCP;
  	}

	Gex = g;
	Sex = s;
	CPex = cp;
	Vex = v;
	Uex = u;
	Hex = Uex + Vex*Pbar;
	Aex = Gex - Vex*Pbar;
	Uex = Hex - Vex*Pbar;

	// assignments (excess properties)
	Zex[0] = Gex;
	Zex[1] = Hex;
	Zex[2] = Sex;
	Zex[3] = CPex;
	Zex[4] = Vex;
	Zex[5] = Aex;
	Zex[6] = Uex;

	return 0;
}


/// calculates ideal mixing properties
long int TRedlichKister::IdealProp( double *Zid )
{

	long int j;
	double si;
	si = 0.0;
	for (j=0; j<NComp; j++)
	{
		if ( x[j] > 1.0e-32 )
			si += x[j]*log(x[j]);
	}
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
	Sid = (-1.)*R_CONST*si;
	Gid = Hid - Sid*Tk;
	Aid = Gid - Vid*Pbar;
	Uid = Hid - Vid*Pbar;

	// assignments (ideal mixing properties)
	Zid[0] = Gid;
	Zid[1] = Hid;
	Zid[2] = Sid;
	Zid[3] = CPid;
	Zid[4] = Vid;
	Zid[5] = Aid;
	Zid[6] = Uid;

	return 0;
}




//=============================================================================================
// NRTL model for liquid solutions
// References: Renon and Prausnitz (1968), Prausnitz et al. (1997)
// (c) TW June 2008
//=============================================================================================


// Generic constructor for the TNRTL class
TNRTL::TNRTL( SolutionData *sd ):
                TSolMod( sd )
{
    alloc_internal();
}


TNRTL::~TNRTL()
{
    free_internal();
}


void TNRTL::alloc_internal()
{
	Tau = new double *[NComp];
	dTau = new double *[NComp];
	d2Tau = new double *[NComp];
	Alp = new double *[NComp];
	dAlp = new double *[NComp];
	d2Alp = new double *[NComp];
	G = new double *[NComp];
	dG = new double *[NComp];
	d2G = new double *[NComp];

    for (long int j=0; j<NComp; j++)
    {
    	Tau[j] = new double [NComp];
    	dTau[j] = new double [NComp];
    	d2Tau[j] = new double [NComp];
    	Alp[j] = new double [NComp];
    	dAlp[j] = new double [NComp];
    	d2Alp[j] = new double [NComp];
		G[j] = new double [NComp];
		dG[j] = new double [NComp];
		d2G[j] = new double [NComp];
	}
}


void TNRTL::free_internal()
{
	// cleaning memory
	for (long int j=0; j<NComp; j++)
	{
		delete[]Tau[j];
	   	delete[]dTau[j];
	   	delete[]d2Tau[j];
	   	delete[]Alp[j];
	   	delete[]dAlp[j];
	   	delete[]d2Alp[j];
		delete[]G[j];
		delete[]dG[j];
		delete[]d2G[j];
	}
	delete[]Tau;
	delete[]dTau;
	delete[]d2Tau;
	delete[]Alp;
	delete[]dAlp;
	delete[]d2Alp;
	delete[]G;
	delete[]dG;
	delete[]d2G;
}


///   Calculates T,P corrected binary interaction parameters
long int TNRTL::PTparam()
{
	long int ip, i, j, i1, i2;
	double A, B, C, D, E, F, tau, dtau, d2tau, alp, dalp, d2alp;

    if ( NPcoef < 6 || NPar < 1 )
       return 1;

	// fill internal arrays of interaction parameters with standard value
	for (j=0; j<NComp; j++)
	{
		for ( i=0; i<NComp; i++ )
		{
			Tau[j][i] = 0.0;
			dTau[j][i] = 0.0;
			d2Tau[j][i] = 0.0;
			Alp[j][i] = 0.0;
			dAlp[j][i] = 0.0;
			d2Alp[j][i] = 0.0;
			G[j][i] = 1.0;
			dG[j][i] = 0.0;
			d2G[j][i] = 0.0;
		}
	}

	// read and convert parameters that have non-standard value
	for (ip=0; ip<NPar; ip++)
	{
		i1 = aIPx[MaxOrd*ip];
		i2 = aIPx[MaxOrd*ip+1];
		A = aIPc[NPcoef*ip+0];
		B = aIPc[NPcoef*ip+1];
		C = aIPc[NPcoef*ip+2];
		D = aIPc[NPcoef*ip+3];
		E = aIPc[NPcoef*ip+4];
		F = aIPc[NPcoef*ip+5];

		tau = A + B/Tk + C*Tk + D*log(Tk);	// partial derivatives of tau and alp
		dtau = - B/pow(Tk,2.) + C + D/Tk;
		d2tau = 2.*B/pow(Tk,3.) - D/pow(Tk,2.);
		alp = E + F*(Tk-273.15);
		dalp = F;
		d2alp = 0.0;

		Tau[i1][i2] = tau;
		dTau[i1][i2] =  dtau;
		d2Tau[i1][i2] = d2tau;
		Alp[i1][i2] = alp;
		dAlp[i1][i2] = dalp;
		d2Alp[i1][i2] =  d2alp;

		G[i1][i2] = exp(-Alp[i1][i2]*Tau[i1][i2]);
		dG[i1][i2] = - ( dAlp[i1][i2]*Tau[i1][i2] + Alp[i1][i2]*dTau[i1][i2] )
				* exp(-Alp[i1][i2]*Tau[i1][i2]);
		d2G[i1][i2] = - ( (d2Alp[i1][i2]*Tau[i1][i2] + 2.*dAlp[i1][i2]*dTau[i1][i2]
				+ Alp[i1][i2]*d2Tau[i1][i2])*G[i1][i2]
				+ (dAlp[i1][i2]*Tau[i1][i2] + Alp[i1][i2]*dTau[i1][i2])*dG[i1][i2] );

		// old version with constant Alp
		// dG[i1][i2] = -Alp[i1][i2] * exp( -Alp[i1][i2]*Tau[i1][i2] ) * dTau[i1][i2];
		// d2G[i1][i2] = -Alp[i1][i2]*(-Alp[i1][i2]*exp(-Alp[i1][i2]*Tau[i1][i2])*dTau[i1][i2]*dTau[i1][i2]
		//		+ exp(-Alp[i1][i2]*Tau[i1][i2])*d2Tau[i1][i2]);
	}
	return 0;
}


/// Calculates activity coefficients
long int TNRTL::MixMod()
{
	long int  j, i, k;
	double K, L, M, N, O, lnGam;

	if ( NPcoef < 6 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
	        return 1;

	// loop over species
	for (j=0; j<NComp; j++)
	{
		lnGam = 0.0;
		K = 0.0;
		L = 0.0;
		M = 0.0;
		for (i=0; i<NComp; i++)
		{
			N = 0.0;
			O = 0.0;
			K += ( x[i]*Tau[i][j]*G[i][j] );
			L += ( x[i]*G[i][j] );
			for (k=0; k<NComp; k++)
			{
				N += ( x[k]*G[k][i] );
				O += ( x[k]*Tau[k][i]*G[k][i] );
			}
			M += ( x[i]*G[j][i]/N * ( Tau[j][i] - O/N ) );
		}
		lnGam = K/L + M;
		lnGamma[j] = lnGam;
	}
	return 0;
}


/// calculates bulk phase excess properties
long int TNRTL::ExcessProp( double *Zex )
{
	long int  j, i;
	double U, dU, V, dV, d2U, d2V, g, dg, d2g;

	if ( NPcoef < 6 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
	        return 1;

	// calculate bulk phase excess properties
   	Gex = 0.0; Sex = 0.0; Hex = 0.0; CPex = 0.0; Vex = 0.0;
   	g = 0.0; dg = 0.0; d2g = 0.0;

   	for (j=0; j<NComp; j++)
   	{
		U = 0.0;
		V = 0.0;
		dU = 0.0;
		dV = 0.0;
		d2U = 0.0;
		d2V = 0.0;
		for (i=0; i<NComp; i++)
		{
			U += x[i]*Tau[i][j]*G[i][j];
			V += x[i]*G[i][j];
			dU += x[i] * ( dTau[i][j]*G[i][j] + Tau[i][j]*dG[i][j] );
			dV += x[i]*dG[i][j];
			d2U += x[i] * ( d2Tau[i][j]*G[i][j] + 2.*dTau[i][j]*dG[i][j]
					+ Tau[i][j]*d2G[i][j] );
			d2V += x[i]*d2G[i][j];
		}
		g += x[j]*U/V;
		dg += x[j] * (dU*V-U*dV)/pow(V,2.);
		d2g += x[j] * ( (d2U*V+dU*dV)/pow(V,2.) - (dU*V)*(2.*dV)/pow(V,3.)
				- (dU*dV+U*d2V)/pow(V,2.) + (U*dV)*(2.*dV)/pow(V,3.) );
	}

   	// final calculations
	Gex = g*R_CONST*Tk;
	Hex = -R_CONST*pow(Tk,2.)*dg;
	Sex = (Hex-Gex)/Tk;
	CPex = -R_CONST * ( 2.*Tk*dg + pow(Tk,2.)*d2g );
	Aex = Gex - Vex*Pbar;
	Uex = Hex - Vex*Pbar;

	// assignments (excess properties)
	Zex[0] = Gex;
	Zex[1] = Hex;
	Zex[2] = Sex;
	Zex[3] = CPex;
	Zex[4] = Vex;
	Zex[5] = Aex;
	Zex[6] = Uex;

	return 0;
}


/// calculates ideal mixing properties
long int TNRTL::IdealProp( double *Zid )
{
	long int j;
	double si;
	si = 0.0;
	for (j=0; j<NComp; j++)
	{
		if ( x[j] > 1.0e-32 )
			si += x[j]*log(x[j]);
	}
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
	Sid = (-1.)*R_CONST*si;
	Gid = Hid - Sid*Tk;
	Aid = Gid - Vid*Pbar;
	Uid = Hid - Vid*Pbar;

	// assignments (ideal mixing properties)
	Zid[0] = Gid;
	Zid[1] = Hid;
	Zid[2] = Sid;
	Zid[3] = CPid;
	Zid[4] = Vid;
	Zid[5] = Aid;
	Zid[6] = Uid;

	return 0;
}




//=============================================================================================
// Wilson model for liquid solutions
// References: Prausnitz et al. (1997)
// (c) TW June 2008
//=============================================================================================


// Generic constructor for the TWilson class
TWilson::TWilson( SolutionData *sd ):
                TSolMod( sd )
{
    alloc_internal();
}


TWilson::~TWilson()
{
    free_internal();
}


void TWilson::alloc_internal()
{
	Lam = new double *[NComp];
	dLam = new double *[NComp];
	d2Lam = new double *[NComp];

    for (long int j=0; j<NComp; j++)
    {
    	Lam[j] = new double [NComp];
		dLam[j] = new double [NComp];
		d2Lam[j] = new double [NComp];
	}
}


void TWilson::free_internal()
{
   	// cleaning memory
   	for (long int j=0; j<NComp; j++)
   	{
   		delete[]Lam[j];
		delete[]dLam[j];
		delete[]d2Lam[j];
	}
	delete[]Lam;
	delete[]dLam;
	delete[]d2Lam;
}


/// Calculates T-corrected interaction parameters
long int TWilson::PTparam()
{
	long int ip, i, j, i1, i2;
	double A, B, C, D, lam, dlam, d2lam;

    if ( NPcoef < 4 || NPar < 1 )
           return 1;

	// fill internal arrays of interaction parameters with standard value
	for (j=0; j<NComp; j++)
	{
		for ( i=0; i<NComp; i++ )
		{
			Lam[j][i] = 1.0;
			dLam[j][i] = 0.0;
			d2Lam[j][i] = 0.0;
		}
	}
	// read and convert parameters that have non-standard value
	for (ip=0; ip<NPar; ip++)
	{
		A = aIPc[NPcoef*ip+0];
		B = aIPc[NPcoef*ip+1];
		C = aIPc[NPcoef*ip+2];
		D = aIPc[NPcoef*ip+3];
		lam = exp( A + B/Tk + C*Tk + D*log(Tk) );
		dlam = lam*( - B/pow(Tk,2.) + C + D/Tk );
		d2lam = dlam*( - B/pow(Tk,2.) + C + D/Tk ) + lam*( 2.*B/pow(Tk,3.) - D/pow(Tk,2.) );

		i1 = aIPx[MaxOrd*ip];
		i2 = aIPx[MaxOrd*ip+1];
		Lam[i1][i2] = lam;
		dLam[i1][i2] = dlam;
		d2Lam[i1][i2] = d2lam;
	}
	return 0;
}


/// Calculates activity coefficients
long int TWilson::MixMod()
{
	long int  j, i, k;
	double K, L, M, lnGam;

	if ( NPcoef < 4 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
	        return 1;

	// loop over species
	for (j=0; j<NComp; j++)
	{
		lnGam = 0.0;
		K = 0.0;
		L = 0.0;
		for (i=0; i<NComp; i++)
		{
			M = 0.0;
			K += x[i]*Lam[j][i];
			for (k=0; k<NComp; k++)
			{
				M += x[k]*Lam[i][k];
			}
			L += x[i]*Lam[i][j]/M;
		}
		lnGam = 1.-log(K)-L;
		lnGamma[j] = lnGam;
	}

	return 0;
}


/// calculates bulk phase excess properties
long int TWilson::ExcessProp( double *Zex )
{
	long int  j, i;
	double U, dU, d2U, g, dg, d2g;

	if ( NPcoef < 4 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
	        return 1;

	// calculate bulk phase excess properties
	Gex = 0.0; Sex = 0.0; Hex = 0.0; CPex = 0.0; Vex = 0.0;
	g = 0.0; dg = 0.0; d2g = 0.0;

	// loop over species
	for (j=0; j<NComp; j++)
	{
		U = 0.0;
		dU = 0.0;
		d2U = 0.0;
		for (i=0; i<NComp; i++)
		{
			U += x[i]*Lam[j][i];
			dU += x[i]*dLam[j][i];
			d2U += x[i]*d2Lam[j][i];
		}
		g -= x[j]*log(U);
		dg -= x[j]*(1./U)*dU;
		d2g -= x[j] * ( (-1./pow(U,2.))*dU*dU + (1./U)*d2U );  // corrected 11.06.2008 (TW)
	}

	// final calculations
	Gex = g*R_CONST*Tk;
	Hex = -R_CONST*pow(Tk,2.)*dg;
	Sex = (Hex-Gex)/Tk;
	CPex = -R_CONST * ( 2.*Tk*dg + pow(Tk,2.)*d2g );

	// assignments (excess properties)
	Aex = Gex - Vex*Pbar;
	Uex = Hex - Vex*Pbar;
	Zex[0] = Gex;
	Zex[1] = Hex;
	Zex[2] = Sex;
	Zex[3] = CPex;
	Zex[4] = Vex;
	Zex[5] = Aex;
	Zex[6] = Uex;

	return 0;
}


/// calculates ideal mixing properties
long int TWilson::IdealProp( double *Zid )
{
	long int j;
	double si;
	si = 0.0;
	for (j=0; j<NComp; j++)
	{
		if ( x[j] > 1.0e-32 )
			si += x[j]*log(x[j]);
	}
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
	Sid = (-1.)*R_CONST*si;
	Gid = Hid - Sid*Tk;
	Aid = Gid - Vid*Pbar;
	Uid = Hid - Vid*Pbar;

	// assignments (ideal mixing properties)
	Zid[0] = Gid;
	Zid[1] = Hid;
	Zid[2] = Sid;
	Zid[3] = CPid;
	Zid[4] = Vid;
	Zid[5] = Aid;
	Zid[6] = Uid;

	return 0;
}




//=============================================================================================
// Ternary Margules (regular) model for solid solutions
// References: Anderson and Crerar (1993), Anderson (2006)
// (c) TW June 2009
//=============================================================================================


// Generic constructor for the TRegular class
TMargules::TMargules( SolutionData *sd ):
                TSolMod( sd )
{
    // empty constructor
}


TMargules::~TMargules()
{
    // empty destructor
}


/// Calculates T,P corrected binary interaction parameters
long int TMargules::PTparam()
{
	if ( NPcoef < 3 || NPar < 4 )
	           return 1;

	// load parameters
	WU12 = aIPc[0];
	WS12 = aIPc[1];
	WV12 = aIPc[2];
	WU13 = aIPc[3];
	WS13 = aIPc[4];
	WV13 = aIPc[5];
	WU23 = aIPc[6];
	WS23 = aIPc[7];
	WV23 = aIPc[8];
	WU123 = aIPc[9];
	WS123 = aIPc[10];
	WV123 = aIPc[11];

	// calculate parameters at (T,P)
	WG12 = WU12 - Tk*WS12 + Pbar*WV12;
	WG13 = WU13 - Tk*WS13 + Pbar*WV13;
	WG23 = WU23 - Tk*WS23 + Pbar*WV23;
	WG123 = WU123 - Tk*WS123 + Pbar*WV123;

	return 0;
}


/// Calculates activity coefficients
long int TMargules::MixMod()
{
	double a12, a13, a23, a123, lnGam1, lnGam2, lnGam3, X1, X2, X3;

	if ( NPcoef < 3 || NPar < 4 || NComp < 3 || !x || !lnGamma )
		return 1;

	a12 = WG12 / (R_CONST*Tk);
	a13 = WG13 / (R_CONST*Tk);
	a23 = WG23 / (R_CONST*Tk);
	a123 = WG123 / (R_CONST*Tk);

	X1 = x[0];
	X2 = x[1];
	X3 = x[2];

	// activity coefficients (normalized by RT)
	lnGam1 = a12*X2*(1.-X1) + a13*X3*(1.-X1) - a23*X2*X3
				+ a123*X2*X3*(1.-2.*X1);
	lnGam2 = a23*X3*(1.-X2) + a12*X1*(1.-X2) - a13*X1*X3
				+ a123*X1*X3*(1.-2.*X2);
	lnGam3 = a13*X1*(1.-X3) + a23*X2*(1.-X3) - a12*X1*X2
				+ a123*X1*X2*(1.-2.*X3);

	// assignments
	lnGamma[0] = lnGam1;
	lnGamma[1] = lnGam2;
	lnGamma[2] = lnGam3;

	return 0;
}


/// calculates bulk phase excess properties
long int TMargules::ExcessProp( double *Zex )
{
	double X1, X2, X3;

	X1 = x[0];
	X2 = x[1];
	X3 = x[2];

	// excess properties
	Gex = X1*X2*WG12 + X1*X3*WG13 + X2*X3*WG23 + X1*X2*X3*WG123;
	Vex = X1*X2*WV12 + X1*X3*WV13 + X2*X3*WV23 + X1*X2*X3*WV123;
	Uex = X1*X2*WU12 + X1*X3*WU13 + X2*X3*WU23 + X1*X2*X3*WU123;
	Sex = X1*X2*WS12 + X1*X3*WS13 + X2*X3*WS23 + X1*X2*X3*WS123;
	CPex = 0.0;
	Hex = Uex + Vex*Pbar;
	Aex = Gex - Vex*Pbar;

	// assignments (excess properties)
	Zex[0] = Gex;
	Zex[1] = Hex;
	Zex[2] = Sex;
	Zex[3] = CPex;
	Zex[4] = Vex;
	Zex[5] = Aex;
	Zex[6] = Uex;

	return 0;
}


/// calculates ideal mixing properties
long int TMargules::IdealProp( double *Zid )
{
	long int j;
	double si;
	si = 0.0;
	for (j=0; j<NComp; j++)
	{
		if ( x[j] > 1.0e-32 )
			si += x[j]*log(x[j]);
	}
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
	Sid = (-1.)*R_CONST*si;
	Gid = Hid - Sid*Tk;
	Aid = Gid - Vid*Pbar;
	Uid = Hid - Vid*Pbar;

	// assignments (ideal mixing properties)
	Zid[0] = Gid;
	Zid[1] = Hid;
	Zid[2] = Sid;
	Zid[3] = CPid;
	Zid[4] = Vid;
	Zid[5] = Aid;
	Zid[6] = Uid;

	return 0;
}




//=============================================================================================
// Binary Margules (subregular) model for solid solutions
// References: Anderson and Crerar (1993), Anderson (2006)
// (c) TW June 2009
//=============================================================================================


// Generic constructor for the TRegular class
TSubregular::TSubregular( SolutionData *sd ):
                TSolMod( sd )
{
    // empty constructor
}


TSubregular::~TSubregular()
{
    // empty destructor
}


/// Calculates T,P corrected binary interaction parameters
long int TSubregular::PTparam()
{
	if ( NPcoef < 3 || NPar < 2 )
	           return 1;

	// load parameters
	WU12 = aIPc[0];
	WS12 = aIPc[1];
	WV12 = aIPc[2];
	WU21 = aIPc[3];
	WS21 = aIPc[4];
	WV21 = aIPc[5];

	// calculate parameters at (T,P)
	WG12 = WU12 - Tk*WS12 + Pbar*WV12;
	WG21 = WU21 - Tk*WS21 + Pbar*WV21;

	return 0;
}


/// Calculates activity coefficients
long int TSubregular::MixMod()
{
	double a1, a2, lnGam1, lnGam2, X1, X2;

	if ( NPcoef < 3 || NPar < 2 || NComp < 2 || !x || !lnGamma )
		return 1;

	a1 = WG12 / (R_CONST*Tk);
	a2 = WG21 / (R_CONST*Tk);

	X1 = x[0];
	X2 = x[1];

	// activity coefficients (normalized by RT)
	lnGam1 = (2.*a2-a1)*X2*X2 + 2.*(a1-a2)*X2*X2*X2;
	lnGam2 = (2.*a1-a2)*X1*X1 + 2.*(a2-a1)*X1*X1*X1;

	// assignments
	lnGamma[0] = lnGam1;
	lnGamma[1] = lnGam2;

	return 0;
}


/// calculates bulk phase excess properties
long int TSubregular::ExcessProp( double *Zex )
{
	double X1, X2;

	X1 = x[0];
	X2 = x[1];

	// excess properties
	Gex = X1*X2*( X2*WG12 + X1*WG21 );
	Vex = X1*X2*( X2*WV12 + X1*WV21 );
	Uex = X1*X2*( X2*WU12 + X1*WU21 );
	Sex = X1*X2*( X2*WS12 + X1*WS21 );
	CPex = 0.0;
	Hex = Uex + Vex*Pbar;
	Aex = Gex - Vex*Pbar;

	// assignments
	Zex[0] = Gex;
	Zex[1] = Hex;
	Zex[2] = Sex;
	Zex[3] = CPex;
	Zex[4] = Vex;
	Zex[5] = Aex;
	Zex[6] = Uex;

	return 0;
}


/// calculates ideal mixing properties
long int TSubregular::IdealProp( double *Zid )
{
	long int j;
	double si;
	si = 0.0;
	for (j=0; j<NComp; j++)
	{
		if ( x[j] > 1.0e-32 )
			si += x[j]*log(x[j]);
	}
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
	Sid = (-1.)*R_CONST*si;
	Gid = Hid - Sid*Tk;
	Aid = Gid - Vid*Pbar;
	Uid = Hid - Vid*Pbar;

	// assignments (ideal mixing properties)
	Zid[0] = Gid;
	Zid[1] = Hid;
	Zid[2] = Sid;
	Zid[3] = CPid;
	Zid[4] = Vid;
	Zid[5] = Aid;
	Zid[6] = Uid;

	return 0;
}




//=============================================================================================
// Binary Guggenheim (Redlich-Kister) model for solid solutions
// References: Anderson and Crerar (1993), Anderson (2006)
// (c) TW June 2009
//=============================================================================================


// Generic constructor for the TBinaryGuggenheim class
TGuggenheim::TGuggenheim( SolutionData *sd ):
                TSolMod( sd )
{
    // empty constructor
}


TGuggenheim::~TGuggenheim()
{
    // empty destructor
}


/// Calculates T,P corrected binary interaction parameters
long int TGuggenheim::PTparam()
{
	if ( NPcoef < 3 || NPar < 1 )
	           return 1;

	// load parameters
	a0 = aIPc[0];
	a1 = aIPc[1];
        a2 = aIPc[2]; // Bugfix was a1 = aIPc[2];
        return 0;
}


/// Calculates activity coefficients
long int TGuggenheim::MixMod()
{
	double lnGam1, lnGam2, X1, X2;

	if ( NPcoef < 3 || NPar < 1 || NComp < 2 || !x || !lnGamma )
		return 1;

	X1 = x[0];
	X2 = x[1];

	// activity coefficients (normalized by RT)
	lnGam1 = X2*X2*( a0 + a1*(3.*X1-X2) + a2*(X1-X2)*(5.*X1-X2) );
	lnGam2 = X1*X1*( a0 - a1*(3.*X2-X1) + a2*(X2-X1)*(5.*X2-X1) );

	// assignments
	lnGamma[0] = lnGam1;
	lnGamma[1] = lnGam2;
	return 0;
}


/// calculates bulk phase excess properties
long int TGuggenheim::ExcessProp( double *Zex )
{
	double X1, X2;

	X1 = x[0];
	X2 = x[1];

	// excess properties
	Gex = (X1*X2*( a0 + a1*(X1-X2) + a2*pow((X1-X2),2.) ))* (R_CONST*Tk);
	Vex = 0.0;
	Uex = (X1*X2*( a0 + a1*(X1-X2) + a2*pow((X1-X2),2.) ))* (R_CONST*Tk);
	Sex = 0.0;
	CPex = 0.0;
	Hex = Uex + Vex*Pbar;
	Aex = Gex - Vex*Pbar;

	// assignments
	Zex[0] = Gex;
	Zex[1] = Hex;
	Zex[2] = Sex;
	Zex[3] = CPex;
	Zex[4] = Vex;
	Zex[5] = Aex;
	Zex[6] = Uex;

	return 0;
}


/// calculates ideal mixing properties
long int TGuggenheim::IdealProp( double *Zid )
{
	long int j;
	double si;
	si = 0.0;
	for (j=0; j<NComp; j++)
	{
		if ( x[j] > 1.0e-32 )
			si += x[j]*log(x[j]);
	}
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
	Sid = (-1.)*R_CONST*si;
	Gid = Hid - Sid*Tk;
	Aid = Gid - Vid*Pbar;
	Uid = Hid - Vid*Pbar;

	// assignments (ideal mixing properties)
	Zid[0] = Gid;
	Zid[1] = Hid;
	Zid[2] = Sid;
	Zid[3] = CPid;
	Zid[4] = Vid;
	Zid[5] = Aid;
	Zid[6] = Uid;

	return 0;
}


//=============================================================================================
// Berman model for multi-component sublattice solid solutions extended with reciprocal terms
// References: Wood & Nicholls (1978); Berman (1990); Price (1985); Aranovich (1991).
// (c) DK/TW December 2010, June 2011; DK July 2014 (added reciprocal terms)
//=============================================================================================

// Generic constructor for the TBerman class
TBerman::TBerman( SolutionData *sd, double *G0 ):
                TSolMod( sd )
{
    alloc_internal();
    G0f = G0;
}

TBerman::~TBerman()
{
    free_internal();
}

// n choose k  from http://stackoverflow.com/questions/15301885/calculate-value-of-n-choose-k
long int TBerman::choose( const long n, const long k )
{
    if( k == 0 ) return 1;
    return (n * choose(n - 1, k - 1)) / k;
}

void TBerman::alloc_internal()
{
    long int j, jk, jx, s, sk, sx, m, mk, mx, r, em;
    long int emx[4], si[4];  // pairwise recip. reactions; max 6 sublattices
    double mnn;

    if( !NSub || !NMoi || NComp < 2 || NSub > 6 )
    {
        NrcR = Nrc = 0;
        return;   // This is not a multi-site model or < 2 EMs or >6 sublattices
    }

    Wu = new double [NPar];
    Ws = new double [NPar];
    Wv = new double [NPar];
    Wpt = new double [NPar];
    NmoS = new long int [NSub];

    fjs = new double *[NComp];
    for( j=0; j<NComp; j++)
    {
       fjs[j] = new double[NSub];
    }

    pyp = new double [NComp];
//    pyn = new double [NComp];
    Grc = new double [NComp];
    oGf =  new double [NComp];

    // Count the number of different moieties on each sublattice
    for( s=0; s< NSub; s++ ) // Cleaning
       NmoS[s] = 0L;
    for( m=0; m<NMoi; m++ ) // Looking through moieties
    {
        bool mf=false; long int mem[8];
        for( s=0; s< NSub; s++ ) // looking through sublattices
           mem[s] = 0;
        for( j=0; j<NComp; j++) // looking through end members
        {
            for( s=0; s< NSub; s++ ) // looking through sublattices
            {
               if( mn[j][s][m] != 0. )
               {
                 mf=true;
                 NmoS[s]++;
                 mem[s]++;
                 break;
               }
            }
            if( mf == true )
               break;
        } // end j
        for( s=0; s< NSub; s++ ) // looking through sublattices
            if( mem[s] == NComp )
               NmoS[s]--;  // don't count a moiety which is the same in all end members
    } // m
    // Count the maximum number of minals L and the number of independent minals M
    // see (Aranovich, 1991 p.27)
    long int L=1, M=-1, C=1, ns=0;
    for( s=0; s< NSub; s++ )
    {
        if( NmoS[s] < 2L )
           continue; // no reciprocal contribution from sublattices with one moiety
        L *= NmoS[s];
        M += NmoS[s];
        ns++;   // counting how many sites have two or more different moieties
    }
    // Computing the number of choices by M from L
    C = choose( L, M );
    NrcR = C; // maximum possible number of reciprocal reactions
    Nrc = 0;
// cout << "NComp=" << NComp << " L=" << L << " M=" << M << " NrcR= " << NrcR << endl;

    if( ns == 2 && MixCode == MR_B_RCPT_ )  // using CEF reciprocal part in Berman model
    // NSub == 2 )
    {
       // Allocate memory for reciprocal reactions DeltaG and indexation
       DGrc = new double [NrcR];
       XrcM = new long int **[NrcR];
       for(r=0; r<NrcR; r++)
       {
          XrcM[r]   = new long int *[4];
          for(s=0; s<4; s++)
          {
             XrcM[r][s] = new long int [2];
          }
       }
       // initializing arrays
       for(r=0; r<NrcR; r++)
           DGrc[r] = 0.;
       for( r=0; r<NrcR; r++)
           for( em=0; em<4; em++)
              for( s=0; s<2; s++)
                 XrcM[r][em][s] = -1;

       Nrc = CollectReciprocalReactions2();
    }
// cout << "Nrc=" << Nrc << " NrcR= " << NrcR << endl;
}

// Collects indexes of end members involved in reciprocal reactions (2 sublattices case)
// as well as indexes of sublattices and substituted moieties on them.
// Returns the total number of reciprocal reactions actually found in this system.
long int TBerman::CollectReciprocalReactions2( void )
{
    long int rn = 0, jb1, jb2, jb3;    // max 3 sublattices
    long int jf0, jf1, jf2, jf3, s0=0, s1=0, s2=1, s3=1;

    for( jf0=0; jf0<NComp; jf0++ ) // looking through all end members
    {
        // Is there any other end member with the same moieties on s2-th sublattice?
        jb2 = 0;
        while( jb2 < NComp )
        {
            jf2 = FindIdenticalSublatticeRow( s2, jf0, jf0, jb2, NComp );
            if( jf2 < 0 )
                break; // no suitable end members found
            // found another end member involved on the right side
            jb1 = 0;
            while( jb1 < NComp )
            {
               jf1 = FindIdenticalSublatticeRow( s1, jf2, jf0, jb1, NComp );
               if( jf1 < 0 )
                      break; // no suitable end members found
               jb3 = 0;
               while( jb3 < NComp )
               {
                  jf3 = FindIdenticalSublatticeRow( s3, jf1, jf0, jb3, NComp);
                  if( jf3 < 0 )
                      break; // no suitable end members found
                  XrcM[rn][0][0] = jf0; XrcM[rn][0][1] = s0;
                  XrcM[rn][2][0] = jf2; XrcM[rn][2][1] = s2;
                  XrcM[rn][1][0] = jf1; XrcM[rn][1][1] = s1;
                  XrcM[rn][3][0] = jf3; XrcM[rn][3][1] = s3;
// cout << "rn=" << rn << " | j0=" << XrcM[rn][0][0] << " s0=" << XrcM[rn][0][1]
//                  << "  j1=" << XrcM[rn][1][0] << " s1=" << XrcM[rn][1][1]
//                  << "  j2=" << XrcM[rn][2][0] << " s2=" << XrcM[rn][2][1]
//                  << "  j3=" << XrcM[rn][3][0] << " s3=" << XrcM[rn][3][1] << endl;
                  rn++;  // next reaction
                  if( rn > NrcR )
                  { return rn-1; } // indexation error
                  jb3 = jf3+1;
               }   // while jb3
               jb1 = jf1+1;
            } // while jb1
            jb2 = jf2+1;
        }   // while jb2
    } // for j0
//    Nrc = rn;
    return rn;  // actual number of processed reciprocal reactions
}

// Looks for an identical row for the sublattice si in end member ji skipping end member ji
// and (another) end member jp among other end members in the interval of end-member indexes
// from jb until je.
// Returns the index of end member in which the identical moieties on this sublattice exists,
// or -1L otherwise (in which case the ji-th end member is not involved in any reciprocal reaction).
long int TBerman::FindIdenticalSublatticeRow(const long si, const long ji,
                                             const long jp, const long jb, const long je )
{
    long int m, j;
    bool match=true;

    for( j = jb; j < je; j++ )
    {
       if( j == ji || j == jp )
         continue;
       match = true;
       for( m=0; m < NMoi; m++ ) // Looking through moieties
       {
          if( mn[ji][si][m] == mn[j][si][m] )
              continue;
          match = false;
          break;
       }
       if(!match)
           continue;  // this site  row in j-th end member does not fit
       return j;      // match found
    }
    return -1L; // no matching end-member/sublattice row was found
}

void TBerman::free_internal()
{
    long int j,r,s;

    delete[]Wu;
    delete[]Ws;
    delete[]Wv;
    delete[]Wpt;

    for( j=0; j<NComp; j++)
    {
       delete[]fjs[j];
    }
    delete[]fjs;

    delete[]Grc;
    delete[]oGf;
    delete[]NmoS;
    delete[]pyp;
//    delete[]pyn;

    if( NSub == 2 && MixCode == MR_B_RCPT_ )  // if using CEF reciprocal part in Berman model
    {
        delete[]DGrc;
        for(r=0; r<NrcR; r++)
        {
            for(s=0; s<4; s++)
            {
                delete[]XrcM[r][s];
            }
        }
        for(r=0; r<NrcR; r++)
        {
           delete[]XrcM[r];
        }
        delete[]XrcM;
    }
}

/// Calculates T-corrected interaction parameters
long int TBerman::PTparam( )
{
    long int ip, j, r, j0, j1, j2, j3;

    if ( NPcoef < 3 || NPar < 1 )
               return 1;

    for (ip=0; ip<NPar; ip++)  // interaction parameters
    {
        Wu[ip] = aIPc[NPcoef*ip];
        Ws[ip] = aIPc[NPcoef*ip+1];
        Wv[ip] = aIPc[NPcoef*ip+2];
        Wpt[ip] = Wu[ip] - Ws[ip]*Tk + Wv[ip]*Pbar;  // This minus is a future problem...
        aIP[ip] = Wpt[ip];
    }
/*
    if( NrcPpc >= 3 && rcpcf != NULL ) // reciprocal parameters provided
    {
        long int TklnTk = Tk * log(Tk);
//cout << endl << " j" << "\trc(j,0)" << "\trc(j,1)" << "\trc(j,2)" << "\tGrc[j]" << "\tG0f[j]" << "\t|oGf[j]" << endl;
        for (j=0; j<NComp; j++)  // Reciprocal and standard Gibbs energy terms
        {
//cout << " " << j << "\t" << rcpcf[NrcPpc*j] << "\t" << rcpcf[NrcPpc*j+1] << "\t" << rcpcf[NrcPpc*j+2];
             Grc[j] = rcpcf[NrcPpc*j] + rcpcf[NrcPpc*j+1]/Tk + rcpcf[NrcPpc*j+2]*TklnTk;
                     // in J/mol iGrc[j] = a + b/T + c*T*lnT  ;
             oGf[j] = G0f[j] + Grc[j]/(R_CONST*Tk); // normalized
             aGEX[j] = Grc[j]/(R_CONST*Tk);  //this sets the respective correction in pm.fDQF[]
//cout << "\t" << Grc[j] << "\t" << G0f[j] << "\t" << oGf[j] << endl;
        }
    }
*/
    for (j=0; j<NComp; j++)  // Reciprocal and standard Gibbs energy terms
         oGf[j] = G0f[j];
    if( !Nrc )
      return 0;

    if( MixCode == MR_B_RCPT_ )  // blocking CEF reciprocal part in Berman model
     {
        // Calculation of DeltaG of reciprocal reactions (NSub==2 only)
        double dGrc; long int i;
        for( r=0; r< Nrc; r++ ) // looking through reactions
        {
           j0 = XrcM[r][0][0]; j1 = XrcM[r][1][0]; j2 = XrcM[r][2][0]; j3 = XrcM[r][3][0];
           dGrc = oGf[j0] + oGf[j1] - oGf[j2] - oGf[j3];  // Aranovich 1991, eq 1.92
// cout << "r=" << r << " : dGrc=" << dGrc << " oGF: " << oGf[j0] << " " << oGf[j1] << " " << oGf[j2] << " " << oGf[j3] << endl;
          DGrc[r] = dGrc;  // normalized!
        }
    }
    return 0;
}


// Calculates ideal config. term and activity coefficients
long int TBerman::MixMod()
{
    long int retCode, j;

    retCode = IdealMixing();
    if(!retCode)
    {
       for(j=0; j<NComp; j++)
           lnGamma[j] += lnGamConf[j];
    }
    retCode = ReciprocalPart();
    if(!retCode)
    {
       for(j=0; j<NComp; j++)
           lnGamma[j] += lnGamRecip[j];
    }
    retCode = ExcessPart();
    if(!retCode)
    {
// cout << endl << " j" << "\tlnGamC" << "\tlnGamR" << "\tlnGamE" << "\t|lnGam" << endl;
        for(j=0; j<NComp; j++)
       {
           lnGamma[j] += lnGamEx[j];
// cout << " " << j << "\t" << lnGamConf[j] << "\t" << lnGamRecip[j] << "\t" << lnGamEx[j] << "\t| " << lnGamma[j] << endl;
       }
    }
    return retCode;
}


// calculates bulk phase excess properties - to be done yet!
long int TBerman::ExcessProp( double *Zex )
{

    // check and add calculation of excess properties here
    long int ip, i1, i2;
    double g, v, s, u;

    if ( NPcoef < 3 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
            return 1;

    // calculate bulk phase excess properties
    g = 0.0; s = 0.0; v = 0.0; u = 0.0;

    for (ip=0; ip<NPar; ip++)
    {
            i1 = aIPx[MaxOrd*ip];
            i2 = aIPx[MaxOrd*ip+1];
            g += x[i1]*x[i2]*Wpt[ip];
            v += x[i1]*x[i2]*Wv[ip];
            u += x[i1]*x[i2]*Wu[ip];
            s -= x[i1]*x[i2]*Ws[ip];
    }

    Gex = g;
    Sex = s;
    CPex = 0.0;
    Vex = v;
    Uex = u;
    Hex = Uex + Vex*Pbar;
    Aex = Gex - Vex*Pbar;
    Uex = Hex - Vex*Pbar;

    // assignments (excess properties)
    Zex[0] = Gex;
    Zex[1] = Hex;
    Zex[2] = Sex;
    Zex[3] = CPex;
    Zex[4] = Vex;
    Zex[5] = Aex;
    Zex[6] = Uex;

    return 0;
}


// calculates ideal mixing properties
long int TBerman::IdealProp( double *Zid )
{
        Hid = 0.0;
        CPid = 0.0;
        Vid = 0.0;
        Sid = ideal_conf_entropy();
        Gid = Hid - Sid*Tk;
        Aid = Gid - Vid*Pbar;
        Uid = Hid - Vid*Pbar;

        // assignments (ideal mixing properties)
        Zid[0] = Gid;
        Zid[1] = Hid;
        Zid[2] = Sid;
        Zid[3] = CPid;
        Zid[4] = Vid;
        Zid[5] = Aid;
        Zid[6] = Uid;

        return 0;
}

double TBerman::KronDelta( const long int j, const long s, const long m )
{
    double krond = 0.;
    if( mn[j][s][m] != 0 )
       krond = 1.;
    return krond;
}

// CEF calculations - computing pyp[j] (product of site fractions for j-th end member) eq 42
double TBerman::PYproduct( const long int j )
{
    double pyp_j = 1., ys;
    long int s, m;

    for(s = 0; s < NSub; s++)
    {
        if( NmoS[s] < 2L )
           continue; // no reciprocal contribution from sublattices with one moiety
        ys = ysigma( j, s );
        pyp_j *= ys;
    } // s
// cout << "j=" << j << " pyp_j=" << pyp_j << endl;
    return pyp_j;
}

// calculates y_sigma(j,s) see eq 42
//
double TBerman::ysigma( const long int j, const long s )
{
    long int m;
    double ys = 0.;
    for( m=0; m < NMoi; m++ ) // Looking through moieties
    {
       if( mn[j][s][m] != 0. ) // This site is occupied in j-th end member
           ys += mn[j][s][m] * y[s][m];
    }
    if( ys > 0. )
       ys /= mns[s];
    return ys;
}

// calculates dGref/d_ysigma, eq 45
//
double TBerman::dGref_dysigma( const long int j, const long int s, const long int ex_j )
{
    long int l, m, nmem;
    double krond=0., ys, dst, dsum=0.;
    bool kron;
    for( l=0; l<NComp; l++ )
    {
       if( l == ex_j )   // this end member is marked to be skipped
           continue;
       ys = 0.; kron = false;
       for( m=0; m < NMoi; m++ ) // Looking through moieties
       {
          nmem = em_howmany( s, m );
          if( nmem == NComp )
              continue;  // ignoring the moiety which is present in all end members
          krond = KronDelta( j, s, m );
          if( krond != 0. )
          {
              ys += mn[l][s][m] * y[s][m];
              kron = true;
          }
       } // m
       if( kron == false )
           continue;  // no moieties belonging to l-th end member that are also
                      // present in j-th end member found in this sublattice
       if( ys > 0. )
       {
           ys /= mns[s];
           dst = oGf[l] * ( pyp[l] / ys );
           dsum += dst;
       }
    } // l
    return dsum;
}

// Temporary: calculates dGref/d_ysm, eq 46 (needs consideration for general case)
//
double TBerman::dGref_dysm( const long int s, const long m, const long int ex_j )
{
    long int l;
    double krond=0., ys, dst, dsum=0.;
    bool kron;
    for( l=0; l<NComp; l++ )
    {
       if( l == ex_j )   // this end member is marked to be skipped
           continue;
       ys = 0.; kron = false;
       krond = KronDelta( l, s, m );
       if( krond != 0. )
       {
           ys = y[s][m];
           kron = true;
       }
/*       for( m=0; m < NMoi; m++ ) // Looking through moieties
       {
          krond = KronDelta( l, s, m );  // j ?
          if( krond != 0. )
          {
              ys += mn[l][s][m] * y[s][m];
              kron = true;
          }
       } // m
*/     if( kron == false )
           continue;  // no moieties belonging to l-th end member found in this sublattice
       if( ys != 0. )
       {
//           ys /= mns[s];
           dst = oGf[l] * ( pyp[l] / ys );
           dsum += dst;
       }
    } // l
    return dsum;
}

// find index of end-member which has moiety m on site s
// returns end-member index or -1 if no end member has this moiety
// on this site
long int TBerman::em_which(const long int s, const long int m, const long int jb, const long int je )
{
    long int l;
    for( l=jb; l<=je; l++ )
    {
       if( mn[l][s][m] != 0 )
           return l;
    }
    return -1L;
}

// find and return the number of end-members that have the moiety m on site s
//
long int TBerman::em_howmany( long int s, long int m )
{
    long int l, jc=0;
    for( l=0; l<NComp; l++ )
    {
       if( mn[l][s][m] != 0 )
           jc++;
    }
    return jc;
}

// calculates ref.frame term (modified CEF, see eq 43)
//
double TBerman::RefFrameTerm( const long int j, const double G_ref )
{
    long int l, m, s, nmem;
    double reftj, sum_s, sum_m, ys, ydp, dgr_dys;

    sum_s = 0.;
    for( s = 0; s < NSub; s++ )  // sublattices
    {
        if( NmoS[s] < 2L )
            continue; // no reciprocal contribution from sublattices with one moiety
        dgr_dys = dGref_dysigma( j, s, -1L );
        sum_s += dgr_dys;
        sum_m = 0.;
        for( m=0; m < NMoi; m++ ) // Looking through moieties
        {
            nmem = em_howmany( s, m );
            if( nmem == NComp )
                continue;  // ignoring the moiety which is present in all end members
//           l = which_em( s, m );
//           if( l < 0 )  // moiety m does not exist on site s
//               continue;
//           if( l == j )    // not sure
//               continue;
//           ys = ysigma( l, s );
           ys = y[s][m];
           dgr_dys = dGref_dysm( s, m, -1L );
           ydp = ys * dgr_dys;
           sum_m += ydp;
        } // m
        sum_s -= sum_m;
    } // s
    reftj = G_ref + sum_s;
    return reftj;
}
/* // This is to experiment with end members with >1 moiety per site
double TBerman::RefFrameTerm( const long int j, const double G_ref )
{
    long int l, s;
    double reftj, sum_s, sum_m, ys, ydp, dgr_dys;

    sum_s = 0.;
    for(s = 0; s < NSub; s++)  // sublattices
    {
        dgr_dys = dGref_dysigma( j, s, -1L );
        sum_s += dgr_dys;
        sum_m = 0.;
        for( l=0; l<NComp; l++ )
        {
           if( l == j )    // not sure
              continue;
           ys = ysigma( l, s );
           dgr_dys = dGref_dysigma( l, s, -1L ); // j
           ydp = ys * dgr_dys;
           sum_m += ydp;
        }  // l
        sum_s -= sum_m;
    } // s
    reftj = G_ref + sum_s;
    return reftj;
}
*/
// calculates part of activity coefficients related to reciprocal energies.
// (interactions between moieties on different sublattices)
// For 2 sublattices also using reciprocal reactions
long int TBerman::ReciprocalPart()
{
    long int j, r, s, m;
    long int xm[4];  // max 4 sublattices
    bool skip;
    double rcSum, rft, yss, ysn;

    for( j=0; j<NComp; j++)
         lnGamRecip[j] = 0.;

    if( MixCode != MR_B_RCPT_ )  // blocking CEF reciprocal part in Berman model
        return 0;

    if( NSub == 0L || NMoi == 0L )
        return 1;   // this is not a multi-site model - bailing out
    if ( NSub > 6L )
        return 2;  //  too many sublattices - bailing out
    // Tables of site fractions y and end-member multiplicities mn, mns have been
    // already calculated in the IdealMixing() - here we just use them.

    // CEF calculations - computing pyp[j] (site fraction products for end members) eq 42
    // and G_ref - total reference Gibbs energy
    double G_ref=0.;
    for( j=0; j<NComp; j++)
    {
        pyp[j] = PYproduct( j );
        G_ref += pyp[j] * oGf[j];
    }
// cout << "G_ref= " << G_ref << endl;
    // Calculation of reciprocal activity terms (modified from CEF, Sundman & Agren, 1981)
    for( j=0; j<NComp; j++)
    {
       rft = RefFrameTerm( j, G_ref );
       lnGamRecip[j] = rft - oGf[j];
// cout << "j=" << j  << " rft=" << rft << " lnGam=" << lnGamRecip[j]
//     << " pyp=" << pyp[j] << endl;
    }
//    if( NSub != 2 )
        return 0;

    if( NSub == 2L && Nrc > 0 )
    {  // using reciprocal reactions for 2-site case
        // We use DGrc - normamized G effects of reciprocal reactions
        // and XrcM - the array of indexes of end members involved in recip. reactions

    for( j=0; j<NComp; j++)
    {
      lnGamRecip[j] = 0.;
      for( r=0; r< Nrc; r++)
      {
         // summation of DeltaGrc contributions,
         // skipping those that involve the moieties in this EM
         // also identify moieties on sublattices of this minal (EM)
         skip = CheckThisReciprocalReaction( r, j, xm );
         if( skip )
             continue;
         rcSum = DGrc[r];
         for(s = 0; s < NSub; s++)
         {
             if( NmoS[s] < 2L )
                continue; // no reciprocal contribution from sublattices with one moiety
             if( xm[s] < 0 )
                 continue;
             rcSum *= y[s][xm[s]];
         }
         lnGamRecip[j] -= rcSum;
      }  // r
cout << " GexRc=" << lnGamRecip[j] << endl;
    }  // j
    }
    return 0;
}

// Checks if this reaction (index r) should be skipped from summation of reciprocal
// reaction excess energy terms for the end member with index j
// Returns in xm the moiety indexes for each sublattice for picking up their site fractions
// (max. 4 sublattices can be considered)
bool TBerman::CheckThisReciprocalReaction( const long int r, const long int j, long int *xm )
{
    return true; // this reaction to be skipped
}

/// calculates part of activity coefficients related to interaction energies
/// between moieties on the same sublattice.
/// (DK/TW June 2011)
long int TBerman::ExcessPart()
{
    long int ip, sp, j, s, m, d, e, f;
    double y0jsm, yWo, Qsm, Qy, ipo, lnGamRT, lnGam;

    if( NSub < 1 || NMoi < 2 || NPar < 1 || NComp < 2 || MaxOrd < 4
        || NPcoef < 3 || !x || !lnGamma )
    {
        for( j=0; j<NComp; j++)
             lnGamEx[j] = 0.;
        return 1;   // this is not a multi-site mixing model - bailing out
    }
    // Cleaning up the fjs array
    for (j=0; j<NComp; j++)
      for( s=0; s<NSub; s++)
          fjs[j][s] =0.;

    // calculating activity coefficients
    for (j=0; j<NComp; j++)
    {
       lnGamRT = 0.;
    // tables of site fractions and end-member multiplicities have been
    // already prepared in the IdealMixing() - here we just use them
       for( s=0; s<NSub; s++)
       {
          for( m=0; m<NMoi; m++)
          {
            // Retrieving the moiety occupancy number in end member y0_j,s,m
            y0jsm = mn[j][s][m] / mns[s];
            if( y0jsm <= 0. )
               continue; // skip - this moiety is not present on s site in this end member

            // looking through the parameters list
            for (ip=0; ip<NPar; ip++)  // interaction parameters indexed with ip
            {
               sp = aIPx[MaxOrd*ip];
               if( sp != s )
                 continue;   // skip - this parameter refers to another sublattice

               d = aIPx[MaxOrd*ip+1];
               e = aIPx[MaxOrd*ip+2];
               f = aIPx[MaxOrd*ip+3];

               // Determining Q_sm
               Qsm = 0.;
               if( d == m )
                   Qsm += 1.;
               if( e == m )
                   Qsm += 1.;

               if( f < 0L )
               {  ipo = 1.; // this is symmetric interaction parameter W_de,s - eq (5.2-5)
                  yWo = y[s][d] * y[s][e] * Wpt[ip];
               }
               else { ipo = 2.; // this is asymmetric interaction parameter W_def,s - eq (5.2-6)
                  if( f == m )
                       Qsm += 1.;
                  yWo = y[s][d] * y[s][e] * y[s][f] * Wpt[ip];
               }
               Qy = Qsm * y0jsm / y[s][m] - ipo;  // eq (5.2-3)
               fjs[j][s] += y0jsm * yWo * Qy;   // fixed 29.06.2011 - was  fjs[j][s] += yWo * Qy;

               // Attention - may still be a problem with the site multiplicity factor!
               // also a problem with accounting of W_de or W_ed ( W_dee or W_eed )
               // More research is needed!  DK 08.07.2011

            }  // ip
         } // m
         lnGamRT += fjs[j][s];
      } // s

       lnGam = lnGamRT/(R_CONST*Tk);
      lnGamEx[j] = lnGam;
   } // j

   return 0;

}

//=============================================================================================
// CEF (Calphad) model for multi-component sublattice solid solutions extended with reciprocal terms
// References: Sundman & Agren (1981); Lucas et al. (2006); Hillert (1998).
// (c) DK/SN since August 2014 (still to change the excess Gibbs energy terms).
//=============================================================================================

// Generic constructor for the TCEFmod class
TCEFmod::TCEFmod( SolutionData *sd, double *G0 ):
                TSolMod( sd )
{
    alloc_internal();
    G0f = G0;
}

TCEFmod::~TCEFmod()
{
    free_internal();
}

void TCEFmod::alloc_internal()
{
    long int j, jk, jx, s, sk, sx, m, mk, mx, r, em;

    if( !NSub || !NMoi || NComp < 2 || NSub > 6 )
        return;   // This is not a multi-site model or < 2 EMs or >6 sublattices

    Wu = new double [NPar];
    Ws = new double [NPar];
    Wc = new double [NPar];
    Wv = new double [NPar];
    Wpt = new double [NPar];
    NmoS = new long int [NSub];

    fjs = new double *[NComp];
    for( j=0; j<NComp; j++)
    {
       fjs[j] = new double[NSub];
    }

    pyp = new double [NComp];
    oGf =  new double [NComp];
    Grc = new double [NComp];

    // Count the number of different moieties on each sublattice
    for( s=0; s< NSub; s++ ) // Cleaning
       NmoS[s] = 0L;
    for( m=0; m<NMoi; m++ ) // Looking through moieties
    {
        bool mf=false; long int mem[8];
        for( s=0; s< NSub; s++ ) // looking through sublattices
           mem[s] = 0;
        for( j=0; j<NComp; j++) // looking through end members
        {
            for( s=0; s< NSub; s++ ) // looking through sublattices
            {
               if( mn[j][s][m] != 0. )
               {
                 mf=true;
                 NmoS[s]++;
                 mem[s]++;
                 break;
               }
            }
            if( mf == true )
               break;
        } // end j
        for( s=0; s< NSub; s++ ) // looking through sublattices
            if( mem[s] == NComp )
               NmoS[s]--;  // don't count a moiety which is the same in all end members
    } // m
}

void TCEFmod::free_internal()
{
    long int j,r,s;

    delete[]Wu;
    delete[]Ws;
    delete[]Wc;
    delete[]Wv;
    delete[]Wpt;

    for( j=0; j<NComp; j++)
    {
       delete[]fjs[j];
    }
    delete[]fjs;

    delete[]Grc;
    delete[]oGf;
    delete[]NmoS;
    delete[]pyp;
//    delete[]pyn;
}


/// Calculates T-corrected interaction parameters
long int TCEFmod::PTparam( )
{
    long int ip, j, r, j0, j1, j2, j3;

    if ( NPcoef < 4 || NPar < 1 )
               return 1;

    for (ip=0; ip<NPar; ip++)  // interaction parameters
    {
        Wu[ip] = aIPc[NPcoef*ip];
        Ws[ip] = aIPc[NPcoef*ip+1];
        Wc[ip] = aIPc[NPcoef*ip+2];
        Wv[ip] = aIPc[NPcoef*ip+3];
        Wpt[ip] = Wu[ip] + Ws[ip]*Tk + Wc[ip]*Tk*log(Tk) + Wv[ip]*Pbar;
        // Lucas 2007 eq 5.66, Wv and Wc terms added for consistency with petrology
        aIP[ip] = Wpt[ip];
    }
/*
    if( NrcPpc >= 3 && rcpcf != NULL ) // reciprocal parameters provided
    {
        long int TklnTk = Tk * log(Tk);
//cout << endl << " j" << "\trc(j,0)" << "\trc(j,1)" << "\trc(j,2)" << "\tGrc[j]" << "\tG0f[j]" << "\t|oGf[j]" << endl;
        for (j=0; j<NComp; j++)  // Reciprocal and standard Gibbs energy terms
        {
//cout << " " << j << "\t" << rcpcf[NrcPpc*j] << "\t" << rcpcf[NrcPpc*j+1] << "\t" << rcpcf[NrcPpc*j+2];
             Grc[j] = rcpcf[NrcPpc*j] + rcpcf[NrcPpc*j+1]/Tk + rcpcf[NrcPpc*j+2]*TklnTk;
                     // in J/mol iGrc[j] = a + b/T + c*T*lnT  ;
             oGf[j] = G0f[j] + Grc[j]/(R_CONST*Tk); // normalized
             aGEX[j] = Grc[j]/(R_CONST*Tk);  //this sets the respective correction in pm.fDQF[]
//cout << "\t" << Grc[j] << "\t" << G0f[j] << "\t" << oGf[j] << endl;
        }
    }
*/
    for (j=0; j<NComp; j++)  // Reciprocal and standard Gibbs energy terms
       oGf[j] = G0f[j];
//    else
    { // no separate reciprocal free energy terms provided
// cout << "NP_DC=" << NP_DC << endl;
        for (j=0; j<NComp; j++)
        {
// cout << " j=" << j;
            if(NP_DC > 0) // use the first DCc coefficient (to be checked!)
           {
                aGEX[j] = aDCc[NP_DC*j]/(R_CONST*Tk);
// cout << " aDCc[j][0]=" << aDCc[NP_DC*j] << " aGEX[j]=" << aGEX[j];
           }
           Grc[j] = 0.;  // in J/mol
           oGf[j] = G0f[j]+aGEX[j]; // normalized
// cout << " G0f[j]=" << G0f[j] << " | oGf[j]=" << oGf[j] << endl;
        }
    }
    return 0;
}


// Calculates ideal config. term and activity coefficients
long int TCEFmod::MixMod()
{
    long int retCode, j;

    retCode = IdealMixing();
    if(!retCode)
    {
       for(j=0; j<NComp; j++)
           lnGamma[j] += lnGamConf[j];
    }
    retCode = ReciprocalPart();
    if(!retCode)
    {
       for(j=0; j<NComp; j++)
           lnGamma[j] += lnGamRecip[j];
    }
    retCode = ExcessPart();
    if(!retCode)
    {
       for(j=0; j<NComp; j++)
           lnGamma[j] += lnGamEx[j];
    }
    return retCode;
}


// calculates bulk phase excess properties - to be done yet!
long int TCEFmod::ExcessProp( double *Zex )
{

    // check and add calculation of excess properties here
    long int ip, i1, i2;
    double g, v, s, u;

    if ( NPcoef < 3 || NPar < 1 || NComp < 2 || MaxOrd < 2 || !x || !lnGamma )
            return 1;

    // calculate bulk phase excess properties
    g = 0.0; s = 0.0; v = 0.0; u = 0.0;

    for (ip=0; ip<NPar; ip++)
    {
            i1 = aIPx[MaxOrd*ip];
            i2 = aIPx[MaxOrd*ip+1];
            g += x[i1]*x[i2]*Wpt[ip];
            v += x[i1]*x[i2]*Wv[ip];
            u += x[i1]*x[i2]*Wu[ip];
            s -= x[i1]*x[i2]*Ws[ip];
    }

    Gex = g;
    Sex = s;
    CPex = 0.0;
    Vex = v;
    Uex = u;
    Hex = Uex + Vex*Pbar;
    Aex = Gex - Vex*Pbar;
    Uex = Hex - Vex*Pbar;

    // assignments (excess properties)
    Zex[0] = Gex;
    Zex[1] = Hex;
    Zex[2] = Sex;
    Zex[3] = CPex;
    Zex[4] = Vex;
    Zex[5] = Aex;
    Zex[6] = Uex;

    return 0;
}


// calculates ideal mixing properties
long int TCEFmod::IdealProp( double *Zid )
{
        Hid = 0.0;
        CPid = 0.0;
        Vid = 0.0;
        Sid = ideal_conf_entropy();
        Gid = Hid - Sid*Tk;
        Aid = Gid - Vid*Pbar;
        Uid = Hid - Vid*Pbar;

        // assignments (ideal mixing properties)
        Zid[0] = Gid;
        Zid[1] = Hid;
        Zid[2] = Sid;
        Zid[3] = CPid;
        Zid[4] = Vid;
        Zid[5] = Aid;
        Zid[6] = Uid;

        return 0;
}

double TCEFmod::KronDelta( const long int j, const long s, const long m )
{
    double krond = 0.;
    if( mn[j][s][m] != 0 )
       krond = 1.;
    return krond;
}

/// CEF - computing pyp[j] (product of site fractions for j-th end member) eq 42
//
double TCEFmod::PYproduct( const long int j )
{
    double pyp_j = 1., ys;
    long int s;

    for(s = 0; s < NSub; s++)
    {
        if( NmoS[s] < 2L )
           continue; // no reciprocal contribution from sublattices with one moiety - no substitution
        ys = ysm( j, s );
        pyp_j *= ys;
    } // s
// cout << "j=" << j << " pyp_j=" << pyp_j << endl;
    return pyp_j;
}

/// picks up y(s,m_j) under CEF (fully substituted end members only)
//
double TCEFmod::ysm( const long int j, const long s )
{
    long int m, nmem;
    double ys = 0.;
    for( m=0; m < NMoi; m++ ) // Looking through moieties
    {
       nmem = em_howmany( s, m );
       if( nmem == NComp )
           continue;  // ignoring the moiety which is present in all end members
       if( mn[j][s][m] != 0. ) // This site is occupied in j-th end member
       {
           ys = y[s][m];
           break;
       }
    }
    return ys;
}

/// calculates dGref/d_ysm, eq 45
//
double TCEFmod::dGref_dysigma( const long int j, const long int s )
{
    long int l, m, nmem;
    double krond=0., ys, dst, dsum=0.;
    bool kron;
    for( l=0; l<NComp; l++ )
    {
       ys = 0.; kron = false;
       for( m=0; m < NMoi; m++ ) // Looking through moieties
       {
          if( em_howmany( s, m ) == NComp )
              continue;  // ignoring the moiety which is present in all end members
          if( KronDelta( j, s, m ) != 0. )
          {
              ys += mn[l][s][m] * y[s][m];
              kron = true;
          }
       } // m
       if( kron == false || ys == 0. )
           continue;  // no moieties belonging to l-th end member that are also
                      // present in j-th end member found in this sublattice
       ys /= mns[s];
       dst = oGf[l] * ( pyp[l] / ys );
       dsum += dst;
    } // l
    return dsum;
}

/// Calculates dGref/d_ysm, eq 46
//
double TCEFmod::dGref_dysm( const long int s, const long int m )
{
    long int l;
    double ys, dst, dsum=0.;
    bool kron;
    for( l=0; l<NComp; l++ )
    {
       ys = 0.; kron = false;
       if( KronDelta( l, s, m ) != 0. )
       {
           ys = y[s][m];
           kron = true;
       }
       if( kron == false || ys == 0. )
           continue;  // no moieties belonging to l-th end member found in this sublattice
       dst = oGf[l] * ( pyp[l] / ys );
       dsum += dst;
    } // l
    return dsum;
}

/// find index of end-member which has moiety m on site s
/// returns end-member index or -1 if no end member has this moiety on this site
//
long int TCEFmod::em_which(const long int s, const long int m, const long int jb, const long int je )
{
    long int l;
    for( l=jb; l<=je; l++ )
    {
       if( mn[l][s][m] != 0 )
           return l;
    }
    return -1L;
}

/// find and return the number of end-members that have the moiety m on site s
//
long int TCEFmod::em_howmany( long int s, long int m )
{
    long int l, jc=0;
    for( l=0; l<NComp; l++ )
    {
       if( mn[l][s][m] != 0 )
           jc++;
    }
    return jc;
}

/// calculates reference frame term (CEF, Sundman and Agren 1981)
//
double TCEFmod::RefFrameTerm( const long int j, const double G_ref )
{
    long int l, m, s, nmem;
    double reftj, sum_s, sum_m, ys, ydp, dgr_dys;

    sum_s = 0.;
    for( s = 0; s < NSub; s++ )  // sublattices
    {
        if( NmoS[s] < 2L )
            continue; // no reciprocal contribution from sublattices with one moiety
        dgr_dys = dGref_dysigma( j, s );
        sum_s += dgr_dys;
        sum_m = 0.;
        for( m=0; m < NMoi; m++ ) // Looking through moieties
        {
            nmem = em_howmany( s, m );
            if( nmem == NComp )
                continue;  // ignoring the moiety which is present in all end members
            ys = y[s][m];
            dgr_dys = dGref_dysm( s, m );
            ydp = ys * dgr_dys;
            sum_m += ydp;
        } // m
        sum_s -= sum_m;
    } // s
    reftj = G_ref + sum_s;
    return reftj;
}

/// CEF: calculates part of activity coefficients related to reciprocal energies
/// (interactions between moieties on different sublattices)
///
long int TCEFmod::ReciprocalPart()
{
    long int j, r, s, m;
    long int xm[4];  // max 4 sublattices
    bool skip;
    double rcSum, rft, yss, ysn;

    for( j=0; j<NComp; j++)
         lnGamRecip[j] = 0.;
    if( NSub == 0L || NMoi == 0L )
        return 1;   // this is not a multi-site model - bailing out
    if ( NSub > 6L )
        return 2;  //  too many sublattices - bailing out
    // Tables of site fractions y and end-member multiplicities mn, mns have been
    // already calculated in the IdealMixing() - here we just use them.

    // CEF calculations - computing pyp[j] (site fraction products for end members) eq 42
    // and G_ref - total reference Gibbs energy
    double G_ref=0.;
    for( j=0; j<NComp; j++)
    {
        pyp[j] = PYproduct( j );
        G_ref += pyp[j] * oGf[j];
    }
// cout << "G_ref= " << G_ref << endl;
    // Calculation of reciprocal activity terms (modified from CEF, Sundman & Agren, 1981)
    for( j=0; j<NComp; j++)
    {
       rft = RefFrameTerm( j, G_ref );
       lnGamRecip[j] = rft - oGf[j];
// cout << "j=" << j  << " rft=" << rft << " lnGam=" << lnGamRecip[j]
//     << " pyp=" << pyp[j] << endl;
    }
//    if( NSub != 2 )
        return 0;
}

/// calculates part of activity coefficients related to interaction energies
/// between moieties on the same sublattice.
/// After Berman (DK/TW June 2011), needs to be changed to CEF (Lukas ea 2007)
//
long int TCEFmod::ExcessPart()
{
    long int ip, sp, j, s, m, d, e, f;
    double y0jsm, yWo, Qsm, Qy, ipo, lnGamRT, lnGam;

    if( NSub < 1 || NMoi < 2 || NPar < 1 || NComp < 2 || MaxOrd < 4
        || NPcoef < 3 || !x || !lnGamma )
    {
        for( j=0; j<NComp; j++)
             lnGamEx[j] = 0.;
        return 1;   // this is not a multi-site mixing model - bailing out
    }
    // Cleaning up the fjs array
    for (j=0; j<NComp; j++)
      for( s=0; s<NSub; s++)
          fjs[j][s] =0.;

    // calculating activity coefficients
    for (j=0; j<NComp; j++)
    {
       lnGamRT = 0.;
    // tables of site fractions and end-member multiplicities have been
    // already prepared in the IdealMixing() - here we just use them
       for( s=0; s<NSub; s++)
       {
          for( m=0; m<NMoi; m++)
          {
            // Retrieving the moiety occupancy number in end member y0_j,s,m
            y0jsm = mn[j][s][m] / mns[s];
            if( y0jsm <= 0. )
               continue; // skip - this moiety is not present on s site in this end member

            // looking through the parameters list
            for (ip=0; ip<NPar; ip++)  // interaction parameters indexed with ip
            {
               sp = aIPx[MaxOrd*ip];
               if( sp != s )
                 continue;   // skip - this parameter refers to another sublattice

               d = aIPx[MaxOrd*ip+1];
               e = aIPx[MaxOrd*ip+2];
               f = aIPx[MaxOrd*ip+3];

               // Determining Q_sm
               Qsm = 0.;
               if( d == m )
                   Qsm += 1.;
               if( e == m )
                   Qsm += 1.;

               if( f < 0L )
               {  ipo = 1.; // this is symmetric interaction parameter W_de,s - eq (5.2-5)
                  yWo = y[s][d] * y[s][e] * Wpt[ip];
               }
               else { ipo = 2.; // this is asymmetric interaction parameter W_def,s - eq (5.2-6)
                  if( f == m )
                       Qsm += 1.;
                  yWo = y[s][d] * y[s][e] * y[s][f] * Wpt[ip];
               }
               Qy = Qsm * y0jsm / y[s][m] - ipo;  // eq (5.2-3)
               fjs[j][s] += y0jsm * yWo * Qy;   // fixed 29.06.2011 - was  fjs[j][s] += yWo * Qy;

               // Attention - may still be a problem with the site multiplicity factor!
               // also a problem with accounting of W_de or W_ed ( W_dee or W_eed )
               // More research is needed!  DK 08.07.2011

            }  // ip
         } // m
         lnGamRT += fjs[j][s];
      } // s

       lnGam = lnGamRT/(R_CONST*Tk);
      lnGamEx[j] = lnGam;
   } // j

   return 0;

}

//--------------------- End of s_solmod3.cpp ----------------------------------------







