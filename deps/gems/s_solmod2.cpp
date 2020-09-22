//-------------------------------------------------------------------
// $Id: s_fgl.cpp 724 2012-10-02 14:25:25Z kulik $
//
/// \file s_solmod2.cpp
/// Implementation of TSolMod derived classes for fluid phase models
/// (TPRSVcalc, TCGFcalc, TSRKcalc, TPR78calc, TCORKcalc and TSTPcalc classes)
//
// Copyright (c) 2004-2012  T.Wagner, D.Kulik, S. Dmitrieva, S.Churakov
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
#include "verror.h"


//=======================================================================================================
// Peng-Robinson-Stryjek-Vera (PRSV) model for fluid mixtures
// References: Stryjek and Vera (1986), Proust and Vera (1989)
// (c) TW July 2006
//=======================================================================================================


// generic constructor
TPRSVcalc::TPRSVcalc( long int NCmp, double Pp, double Tkp ):
        TSolMod( NCmp, 'P', Tkp, Pp )

{
    // aGEX = 0;
    // aVol = 0;
    Pparc = 0;
    alloc_internal();
}



TPRSVcalc::TPRSVcalc( SolutionData *sd ):
                TSolMod( sd )
{
    // aGEX = arGEX;
    // aVol = arVol;
    Pparc = aPparc;
    alloc_internal();
}



TPRSVcalc::~TPRSVcalc()
{
	free_internal();
}



/// allocate work arrays for pure fluid and fluid mixture properties
void TPRSVcalc::alloc_internal()
{
	Eosparm = new double [NComp][6];
	Pureparm = new double [NComp][4];
	Fugpure = new double [NComp][6];
	Fugci = new double [NComp][4];
	a = new double *[NComp];
	b = new double *[NComp];
	KK = new double *[NComp];
	dKK = new double *[NComp];
	d2KK = new double *[NComp];
	AA = new double *[NComp];

    for (long int i=0; i<NComp; i++)
    {
    	a[i] = new double[NComp];
    	b[i] = new double[NComp];
    	KK[i] = new double[NComp];
    	dKK[i] = new double[NComp];
    	d2KK[i] = new double[NComp];
    	AA[i] = new double[NComp];
    }
}


void TPRSVcalc::free_internal()
{
	long int i;

	for (i=0; i<NComp; i++)
	{
		delete[]a[i];
		delete[]b[i];
		delete[]KK[i];
		delete[]dKK[i];
		delete[]d2KK[i];
		delete[]AA[i];
	}

	delete[]Eosparm;
	delete[]Pureparm;
	delete[]Fugpure;
	delete[]Fugci;
	delete[]a;
	delete[]b;
	delete[]KK;
	delete[]dKK;
	delete[]d2KK;
	delete[]AA;
}



/// high-level method to retrieve pure fluid fugacities
long int TPRSVcalc::PureSpecies()
{
    long int j, retCode = 0;

    for( j=0; j<NComp; j++)
    {
        // Calling PRSV EoS for pure fugacity
        retCode =  FugacityPT( j, aDCc+j*NP_DC );
        aGEX[j] = log( Fugpure[j][0] );
        Pparc[j] = Fugpure[j][0]*Pbar;  // fure fluid fugacity (required for performance)
        aVol[j] = Fugpure[j][4]*10.;  // molar volume of pure fluid component, J/bar to cm3
    } // j

    if ( retCode )
    {
        char buf[150];
        sprintf(buf, "PRSV fluid: calculation of pure fugacity failed");
                Error( "E71IPM IPMgamma: ",  buf );
    }
    return 0;
}



/// high-level method to calculate T,P corrected binary interaction parameters
long int TPRSVcalc::PTparam()
{
	long int j, i;

	PureSpecies();

	// set all interaction parameters zero
	for( j=0; j<NComp; j++ )
	{
		for( i=0; i<NComp; i++ )
		{
			KK[j][i] = 0.;
			dKK[j][i] = 0.;
			d2KK[j][i] = 0.;
		}
	}

	switch ( MixCode )
	{
        case MR_UNDEF_:
        case MR_WAAL_:
			MixingWaals();
			break;
		case MR_CONST_:
			MixingConst();
			break;
		case MR_TEMP_:
			MixingTemp();
			break;
		default:
			break;
	}

	return 0;
}



/// high-level method to retrieve activity coefficients of the fluid mixture
long int TPRSVcalc::MixMod()
{
	long int j, iRet;

    iRet = FugacitySpec( aPparc );

    phVOL[0] = PhVol * 10.;

    for(j=0; j<NComp; j++)
    {
    	if( Fugci[j][3] > 1e-23 )
    		lnGamma[j] = log( Fugci[j][3] );
        else
        	lnGamma[j] = 0;
    }
    if ( iRet )
    {
    	char buf[150];
    	sprintf(buf, "PRSV fluid: calculation failed");
    	Error( "E71IPM IPMgamma: ",  buf );
    }
    return iRet;
}



/// high-level method to retrieve residual functions of the fluid mixture
long int TPRSVcalc::ExcessProp( double *Zex )
{
	long int iRet;

        iRet = ResidualFunct( aPparc );

    if ( iRet )
    {
    	char buf[150];
    	sprintf(buf, "PRSV fluid: calculation failed");
    	Error( "E71IPM IPMgamma: ",  buf );
    }

    Ars = Grs - Vrs*Pbar;
    Urs = Hrs - Vrs*Pbar;

	// assignments (residual functions)
	Zex[0] = Grs;
	Zex[1] = Hrs;
	Zex[2] = Srs;
	Zex[3] = CPrs;
	Zex[4] = Vrs;
	Zex[5] = Ars;
	Zex[6] = Urs;

	return iRet;
}



/// calculates ideal mixing properties
long int TPRSVcalc::IdealProp( double *Zid )
{
	long int j;
	double s, sc, sp;

	s = 0.0;
	for (j=0; j<NComp; j++)
	{
		if ( x[j] > 1.0e-32 )
			s += x[j]*log(x[j]);
	}
	sc = (-1.)*R_CONST*s;
	sp = (-1.)*R_CONST*log(Pbar);
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
	Sid = sc + sp;
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



/// basic van der waals mixing rule
long int TPRSVcalc::MixingWaals()
{
	// currently no calculations

	return 0;
}



/// constant one-term interaction parameter
long int TPRSVcalc::MixingConst()
{
	long int ip, i1, i2;
	double k, dk, d2k;

	if( NPcoef > 0 )
	{
		// transfer interaction parameters that have non-standard value
		for ( ip=0; ip<NPar; ip++ )
		{
			i1 = aIPx[MaxOrd*ip];
			i2 = aIPx[MaxOrd*ip+1];
			k = aIPc[NPcoef*ip];
			dk = 0.;
			d2k = 0.;
			KK[i1][i2] = k;
			dKK[i1][i2] = dk;
			d2KK[i1][i2] = d2k;
			KK[i2][i1] = k;   // symmetric case
			dKK[i2][i1] = dk;
			d2KK[i2][i1] = d2k;
		}
	}

	return 0;
}



/// temperature dependent one-term interaction parameter
long int TPRSVcalc::MixingTemp()
{
	long int i, j, ip, i1, i2;
	double ai, aj, bi, bj, di, dj, dai, daj, d2ai, d2aj, ddi, ddj, d2di, d2dj,
				U, V, dU, dV, d2U, d2V, tmp, k, dk, d2k, C;

	// set model specific interaction parameters zero
	for( j=0; j<NComp; j++ )
	{
		for( i=0; i<NComp; i++ )
		{
			a[j][i] = 0.;
			b[j][i] = 0.;
		}
	}

	if( NPcoef > 0 )
	{
		// transfer parameters that have non-standard value
		for ( ip=0; ip<NPar; ip++ )
		{
			i1 = aIPx[MaxOrd*ip];
			i2 = aIPx[MaxOrd*ip+1];
			a[i1][i2] = aIPc[NPcoef*ip];
			b[i1][i2] = aIPc[NPcoef*ip+1];
			a[i2][i1] = aIPc[NPcoef*ip];  // symmetric case
			b[i2][i1] = aIPc[NPcoef*ip+1];
		}
	}

	// calculate binary interaction parameters
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
			if ( a[i][j] == 0.0 )
				tmp = 1.0;
			else
				tmp = a[i][j];

			// read a, b, da, d2a
			ai = Pureparm[i][0];
			aj = Pureparm[j][0];
			bi = Pureparm[i][1];
			bj = Pureparm[j][1];
			dai = Pureparm[i][2];
			daj = Pureparm[j][2];
			d2ai = Pureparm[i][3];
			d2aj = Pureparm[j][3];

			// calculate k and derivatives
			di = sqrt(ai)/bi;
			dj = sqrt(aj)/bj;
			ddi = (0.5/bi) * pow(ai,-0.5) * dai;
			ddj = (0.5/bj) * pow(aj,-0.5) * daj;
			d2di = (0.5/bi) * ( (-0.5)*pow(ai,-1.5)*dai*dai + pow(ai,-0.5)*d2ai );
			d2dj = (0.5/bj) * ( (-0.5)*pow(aj,-1.5)*daj*daj + pow(aj,-0.5)*d2aj );

			C = ( b[i][j]/tmp - 1. );
			U = a[i][j]*pow((298.15/Tk),C) - pow((di-dj),2.);
			V = (2.*di*dj);
			dU = - ( a[i][j]*C*pow((298.15/Tk),C) ) / Tk - 2.*(di-dj)*(ddi-ddj);
			dV = 2.*( ddi*dj + di*ddj );
			d2U = ( a[i][j]*pow(C,2.)*pow((298.15/Tk),C) ) / pow(Tk,2.)
						+ ( a[i][j]*C*pow((298.15/Tk),C) ) / pow(Tk,2.)
						- 2.*( pow((ddi-ddj),2.) + (di-dj)*(d2di-d2dj) );
			d2V = 2.*( d2di*dj + 2.*ddi*ddj + di*d2dj );
			k = U/V;
			dk = (dU*V-U*dV)/pow(V,2.);
			d2k = (d2U*V+dU*dV)*pow(V,2.)/pow(V,4.) - (dU*V)*(2.*V*dV)/pow(V,4.)
						- (dU*dV+U*d2V)*pow(V,2.)/pow(V,4.) + (U*dV)*(2.*V*dV)/pow(V,4.);

			// assignments
			KK[i][j] = k;
			dKK[i][j] = dk;
			d2KK[i][j] = d2k;
		}
	}

	return 0;
}



/// retrieve pure fluid properties
long int TPRSVcalc::FugacityPT( long int i, double *EoSparam )
{
	long int iRet = 0;
    double Tcrit, Pcrit, omg, k1, k2, k3, apure, bpure, da, d2a;

    // reads EoS parameters from database into work array
    if( !EoSparam )
    	return -1;  // Memory alloc error

    Eosparm[i][0] = EoSparam[0];   // critical temperature in K
    Eosparm[i][1] = EoSparam[1];   // critical pressure in bar
    Eosparm[i][2] = EoSparam[2];   // Pitzer acentric factor omega
    Eosparm[i][3] = EoSparam[3];   // empirical EoS parameter k1
    Eosparm[i][4] = EoSparam[4];   // empirical EoS parameter k2
    Eosparm[i][5] = EoSparam[5];   // empirical EoS parameter k3
    Tcrit = Eosparm[i][0];
    Pcrit = Eosparm[i][1];
    omg = Eosparm[i][2];
    k1 = Eosparm[i][3];
	k2 = Eosparm[i][4];
	k3 = Eosparm[i][5];

	AB( Tcrit, Pcrit, omg, k1, k2, k3, apure, bpure, da, d2a );

	Pureparm[i][0] = apure;  // a parameter
	Pureparm[i][1] = bpure;  // b parameter
	Pureparm[i][2] = da;  // da/dT
	Pureparm[i][3] = d2a;  // d2a/dT2

	iRet = FugacityPure( i );
	if( iRet)
		return iRet;

	return iRet;
}



/// calculates attractive (a) and repulsive (b) parameter of PRSV equation of state
/// and partial derivatives of alpha function
long int TPRSVcalc::AB( double Tcrit, double Pcrit, double omg, double k1, double k2, double k3,
		double &apure, double &bpure, double &da, double &d2a )
{
	double Tred, k0, k, alph, ac, sqa, dsqa, d2sqa;

	Tred = Tk/Tcrit;
	k0 = 0.378893 + 1.4897153*omg - 0.17131848*pow(omg,2.) + 0.0196554*pow(omg,3.);
	if(Tk >= Tcrit)
	{
		k1 = 0.0;
		k2 = 0.0;
		k3 = 0.0;
	}
	k = k0 + (k1 + k2*(k3-Tred)*(1.-sqrt(Tred))) * (1.+sqrt(Tred)) * (0.7-Tred);
	alph = pow(1. + k*(1.-sqrt(Tred)), 2.);
	ac = (0.457235)*pow(R_CONST,2.)*pow(Tcrit,2.) / Pcrit;
	apure = alph*ac;
	bpure = (0.077796)*R_CONST*Tcrit/Pcrit;
	sqa = 1.+k*(1.-sqrt(Tred));
	// dsqa = (-1.)*k0/(2.*sqrt(Tk*Tcrit)) - 1.7*k1/Tcrit + 2.*k1*Tk/(pow(Tcrit,2.));  // extend dA/dT for k2, k3
	dsqa = ( k1*(0.7-Tred)/(2.*sqrt(Tred)*Tcrit) - k1*(1.+sqrt(Tred))/Tcrit ) * (1.-sqrt(Tred))
				- (k0 + k1*(1.+sqrt(Tred))*(0.7-Tred))/(2*sqrt(Tred)*Tcrit);
	da = 2.*ac*(sqa*dsqa);
	d2sqa = ( - (k1*(0.7-Tred))/(4.*pow(Tred,1.5)*pow(Tcrit,2.)) - k1/(sqrt(Tred)*pow(Tcrit,2.)) ) * (1.-sqrt(Tred))
				+ ( - k1*(0.7-Tred)/(2.*sqrt(Tred)*Tcrit) - k1*(1.+sqrt(Tred))/Tcrit ) / (sqrt(Tred)*Tcrit)
				+ ( k0 + k1*(1.+sqrt(Tred))*(0.7-Tred) ) / (4.*pow(Tred,1.5)*pow(Tcrit,2.));
	d2a = 2.*ac*(dsqa*dsqa + sqa*d2sqa);

	return 0;
}



/// calculates fugacities and residual functions of pure fluid species
long int TPRSVcalc::FugacityPure( long int i )
{
	double Tcrit, Pcrit, Tred, aprsv, bprsv, alph, da, d2a, k, A, B, a2, a1, a0,
			z1, z2, z3, vol1, vol2, vol3, lnf1, lnf2, lnf3, z, vol, lnf;
	double gig, hig, sig, cpig, fugpure, grs, hrs, srs, cprs,
			cv, dPdT, dPdV, dVdT;

	// ideal gas changes from 1 bar to P at T of interest
	hig = 0.;
	sig = (-1.)*R_CONST*log(Pbar);
	gig = hig - Tk*sig;
	cpig = 0.;

	// retrieve a and b terms of cubic EoS
	Tcrit = Eosparm[i][0];
	Pcrit = Eosparm[i][1];
	Tred = Tk/Tcrit;
	aprsv = Pureparm[i][0];
	bprsv = Pureparm[i][1];
	da = Pureparm[i][2];
	d2a = Pureparm[i][3];

	// solve cubic equation
	A = aprsv*Pbar/(pow(R_CONST,2.)*pow(Tk,2.));
	B = bprsv*Pbar/(R_CONST*Tk);
	a2 = B - 1.;
	a1 = A - 3.*pow(B,2.) - 2.*B;
	a0 = pow(B,3.) + pow(B,2.) - A*B;
	Cardano(a2, a1, a0, z1, z2, z3);

	// find stable roots
	vol1 = z1*R_CONST*Tk/Pbar;
	vol2 = z2*R_CONST*Tk/Pbar;
	vol3 = z3*R_CONST*Tk/Pbar;
	if (z1 > B)
		lnf1 = (-1.)*log(z1-B)
			- A/(B*sqrt(8.))*log((z1+(1.+sqrt(2.))*B)/(z1+(1.-sqrt(2.))*B))+z1-1.;
	else
		lnf1 = 1000.;
	if (z2 > B)
		lnf2 = (-1.)*log(z2-B)
			- A/(B*sqrt(8.))*log((z2+(1.+sqrt(2.))*B)/(z2+(1.-sqrt(2.))*B))+z2-1.;
	else
		lnf2 = 1000.;
	if (z3 > B)
		lnf3 = (-1.)*log(z3-B)
			- A/(B*sqrt(8.))*log((z3+(1.+sqrt(2.))*B)/(z3+(1.-sqrt(2.))*B))+z3-1.;
	else
		lnf3 = 1000.;

	if (lnf2 < lnf1)
	{
		z = z2; vol = vol2; lnf = lnf2;
	}
	else
	{
		z = z1; vol = vol1; lnf = lnf1;
	}
	if (lnf3 < lnf)
	{
		z = z3; vol = vol3; lnf = lnf3;
	}
	else
	{
		z = z; vol = vol; lnf = lnf;
	}

	// calculate thermodynamic properties
	alph = aprsv/((0.457235)*pow(R_CONST,2.)*pow(Tcrit,2.) / Pcrit);
	k = (sqrt(alph)-1.)/(1.-sqrt(Tred));
	grs = R_CONST*Tk*(z-1.-log(z-B)-A/(B*sqrt(8.))
				*log((z+(1+sqrt(2.))*B)/(z+(1-sqrt(2.))*B)));
	hrs = R_CONST*Tk*(z-1.-log((z+(1+sqrt(2.))*B)/(z+(1-sqrt(2.))*B))
				*A/(B*sqrt(8.))*(1+k*sqrt(Tred)/sqrt(alph)));
	srs = (hrs-grs)/Tk;

	// heat capacity part
	cv = Tk*d2a/(bprsv*sqrt(8.))
			 * log( (z+B*(1.+sqrt(2.)))/(z+B*(1.-sqrt(2.))) );
	dPdT = R_CONST/(vol-bprsv) - da/( vol*(vol+bprsv) + bprsv*(vol-bprsv) );
	dPdV = - R_CONST*Tk/pow((vol-bprsv),2.) + 2.*aprsv*(vol+bprsv)/pow((vol*(vol+bprsv)+bprsv*(vol-bprsv)),2.);
	dVdT = (-1.)*(1./dPdV)*dPdT;
	cprs = cv + Tk*dPdT*dVdT - R_CONST;

	// increment thermodynamic properties
	fugpure = exp(lnf);
	Fugpure[i][0] = fugpure;
	Fugpure[i][1] = grs;  // changed to residual functions, 31.05.2008 (TW)
	Fugpure[i][2] = hrs;
	Fugpure[i][3] = srs;
    Fugpure[i][4] = vol;
    Fugpure[i][5] = cprs;

    return 0;
}



/// cubic equation root solver based on Cardanos method
long int TPRSVcalc::Cardano( double a2, double a1, double a0, double &z1, double &z2, double &z3 )
{
	double q, rc, q3, rc2, theta, ac, bc;

	q = (pow(a2,2.) - 3.*a1)/9.;
	rc = (2.*pow(a2,3.) - 9.*a2*a1 + 27.*a0)/54.;
	q3 = pow(q,3.);
	rc2 = pow(rc,2.);
	if (rc2 < q3)  // three real roots
	{
		theta = acos(rc/sqrt(q3));
		z1 = (-2.)*sqrt(q)*cos(theta/3.)-a2/3.;
		z2 = (-2.)*sqrt(q)*cos(theta/3.+2./3.*3.1415927)-a2/3.;
		z3 = (-2.)*sqrt(q)*cos(theta/3.-2./3.*3.1415927)-a2/3.;
	}
	else  // one real root
	{
		ac = (-1.)*rc/fabs(rc)*pow(fabs(rc)+sqrt(rc2-q3), 1./3.);
		if (ac != 0.)
			bc = q/ac;
		else
			bc = 0.;
		z1 = ac+bc-a2/3.;
		z2 = ac+bc-a2/3.;
		z3 = ac+bc-a2/3.;
	}
	return 0;
}



/// calculates mixing properties of the fluid mixture
long int TPRSVcalc::MixParam( double &amix, double &bmix )
{
	long int i, j;
	double K;
	amix = 0.;
	bmix = 0.;

	// calculate binary aij parameters
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
            K = KK[i][j];
			AA[i][j] = sqrt(Pureparm[i][0]*Pureparm[j][0])*(1.-K);
		}
	}
	// find a and b of the mixture
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
			amix = amix + x[i]*x[j]*AA[i][j];
		}
	}
	for (i=0; i<NComp; i++)
	{
		bmix = bmix + x[i]*Pureparm[i][1];
	}
	return 0;
}



/// calculates fugacity of the bulk fluid mixture
long int TPRSVcalc::FugacityMix( double amix, double bmix, double &fugmix, double &zmix,
		double &vmix )
{
	double A, B, a2, a1, a0, z1, z2, z3, vol1, vol2, vol3, lnf1, lnf2, lnf3, lnf;

	// solve cubic equation
	A = amix*Pbar/(pow(R_CONST,2.)*pow(Tk,2.));
	B = bmix*Pbar/(R_CONST*Tk);
	a2 = B - 1.;
	a1 = A - 3.*pow(B,2.) - 2.*B;
	a0 = pow(B,3.) + pow(B,2.) - A*B;
	Cardano( a2, a1, a0, z1, z2, z3 );

	// find stable roots
	vol1 = z1*R_CONST*Tk/Pbar;
	vol2 = z2*R_CONST*Tk/Pbar;
	vol3 = z3*R_CONST*Tk/Pbar;
	if (z1 > B)
		lnf1 = (-1.)*log(z1-B)
			- A/(B*sqrt(8.))*log((z1+(1.+sqrt(2.))*B)/(z1+(1.-sqrt(2.))*B))+z1-1.;
	else
		lnf1 = 1000.;
	if (z2 > B)
		lnf2 = (-1.)*log(z2-B)
			- A/(B*sqrt(8.))*log((z2+(1.+sqrt(2.))*B)/(z2+(1.-sqrt(2.))*B))+z2-1.;
	else
		lnf2 = 1000.;
	if (z3 > B)
		lnf3 = (-1.)*log(z3-B)
			- A/(B*sqrt(8.))*log((z3+(1.+sqrt(2.))*B)/(z3+(1.-sqrt(2.))*B))+z3-1.;
	else
		lnf3 = 1000.;

	if (lnf2 < lnf1)
	{
		zmix = z2; vmix = vol2; lnf = lnf2;
	}
	else
	{
		zmix = z1; vmix = vol1; lnf = lnf1;
	}
	if (lnf3 < lnf)
	{
		zmix = z3; vmix = vol3; lnf = lnf3;
	}
	else
	{
		zmix = zmix; vmix = vmix; lnf = lnf;
	}
	fugmix = exp(lnf);
        PhVol = vmix;
	return 0;
}



/// calculates fugacities and activities of fluid species in the mixture,
long int TPRSVcalc::FugacitySpec( double *fugpure )
{
    long int i, j, iRet=0;
	double fugmix=0., zmix=0., vmix=0., amix=0., bmix=0., sum=0.;
	double A, B, lnfci, fci;

    // Reload params to Pureparm
    for( j=0; j<NComp; j++ )
    {
      Fugpure[j][0] = fugpure[j]/Pbar;
    }

	// retrieve properties of the mixture
	iRet = MixParam( amix, bmix );
	iRet = FugacityMix( amix, bmix, fugmix, zmix, vmix );
	A = amix*Pbar/(pow(R_CONST, 2.)*pow(Tk, 2.));
	B = bmix*Pbar/(R_CONST*Tk);

	// calculate fugacity coefficient, fugacity and activity of species i
	for (i=0; i<NComp; i++)
	{
		sum = 0.;
		for (j=0; j<NComp; j++)
		{
			sum = sum + x[j]*AA[i][j];
		}
		lnfci = Pureparm[i][1]/bmix*(zmix-1.) - log(zmix-B)
		      + A/(sqrt(8.)*B)*(2.*sum/amix-Pureparm[i][1]/bmix)
                      * log((zmix+B*(1.-sqrt(2.)))/(zmix+B*(1.+sqrt(2.))));
		fci = exp(lnfci);
		Fugci[i][0] = fci;  // fugacity coefficient using engineering convention
		Fugci[i][1] = x[i]*fci;  // fugacity coefficient using geology convention
		Fugci[i][2] = Fugci[i][1]/Fugpure[i][0];  // activity of species
		if (x[i]>1.0e-20)
			Fugci[i][3] = Fugci[i][2]/x[i];  // activity coefficient of species
		else
			Fugci[i][3] = 1.0;
	}

	return iRet;
}



/// calculates residual functions in the mixture
long int TPRSVcalc::ResidualFunct( double *fugpure )
{
    long int i, j, iRet=0;
	double fugmix=0., zmix=0., vmix=0., amix=0., bmix=0.;
	double A, B, K, dK, d2K, Q, dQ, d2Q, damix, d2amix, ai, aj, dai, daj, d2ai, d2aj,
			cv, dPdT, dPdV, dVdT;

    // Reload params to Pureparm (probably now obsolete?)
    for( j=0; j<NComp; j++ )
    {
      Fugpure[j][0] = fugpure[j]/Pbar;
    }

	// retrieve properties of the mixture
	iRet = MixParam( amix, bmix );
	iRet = FugacityMix( amix, bmix, fugmix, zmix, vmix );
	A = amix*Pbar/(pow(R_CONST, 2.)*pow(Tk, 2.));
	B = bmix*Pbar/(R_CONST*Tk);

	// calculate total state functions of the mixture
	damix = 0.;
	d2amix = 0.;
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
			// pull parameters
			ai = Pureparm[i][0];
			aj = Pureparm[j][0];
			dai = Pureparm[i][2];
			daj = Pureparm[j][2];
			d2ai = Pureparm[i][3];
			d2aj = Pureparm[j][3];
			K = KK[i][j];
			dK = dKK[i][j];
			d2K = d2KK[i][j];

			// increments to derivatives
			Q = sqrt(ai*aj);
			dQ = 0.5*( sqrt(aj/ai)*dai + sqrt(ai/aj)*daj );
			d2Q = 0.5*( dai*daj/sqrt(ai*aj) + d2ai*sqrt(aj)/sqrt(ai) + d2aj*sqrt(ai)/sqrt(aj)
					- 0.5*( pow(dai,2.)*sqrt(aj)/sqrt(pow(ai,3.))
					+ pow(daj,2.)*sqrt(ai)/sqrt(pow(aj,3.)) ) );
			damix = damix + x[i]*x[j] * ( dQ*(1.-K) - Q*dK );
			d2amix = d2amix + x[i]*x[j] * ( d2Q*(1.-K) - 2.*dQ*dK - Q*d2K );
		}
	}

	// calculate thermodynamic properties
	Grs = (amix/(R_CONST*Tk*sqrt(8.)*bmix) * log((vmix+(1.-sqrt(2.))*bmix)
		/ (vmix+(1.+sqrt(2.))*bmix))-log(zmix*(1.-bmix/vmix))+zmix-1.)*R_CONST*Tk;
	Hrs = ((amix-Tk*damix)/(R_CONST*Tk*sqrt(8.)*bmix)*log((vmix+(1.-sqrt(2.))
		*bmix)/(vmix+(1.+sqrt(2.))*bmix))+zmix-1.)*R_CONST*Tk;
	Srs = (Hrs - Grs)/Tk;

	// heat capacity part
	cv = Tk*d2amix/(bmix*sqrt(8.))
			 * log( (zmix+B*(1.+sqrt(2.)))/(zmix+B*(1.-sqrt(2.))) );
	dPdT = R_CONST/(vmix-bmix) - damix/( vmix*(vmix+bmix) + bmix*(vmix-bmix) );
	dPdV = - R_CONST*Tk/pow((vmix-bmix),2.) + 2.*amix*(vmix+bmix)/pow((vmix*(vmix+bmix)+bmix*(vmix-bmix)),2.);
	dVdT = (-1.)*(1./dPdV)*dPdT;
	CPrs = cv + Tk*dPdT*dVdT - R_CONST;
	Vrs = vmix;

	return iRet;
}



#ifndef IPMGEMPLUGIN

/// calculates properties of pure fluids when called from DCthermo
long int TPRSVcalc::PRSVCalcFugPure( double Tmin, float *Cpg, double *FugProps )
{
	long int retCode = 0;
	double Coeff[7];

	for( int ii=0; ii<7; ii++ )
		Coeff[ii] = (double)Cpg[ii];

	if( (Tk >= Tmin) && (Tk < 1e4) && (Pbar >= 1e-5) && (Pbar < 1e5) )
	{
		retCode = FugacityPT( 0, Coeff );
		for( int i=0; i<6; i++ )
			FugProps[i] = Fugpure[0][i];
		return retCode;
	}

	else
	{
		for( int i=1; i<6; i++ )
			FugProps[i] = 0.;
		FugProps[0] = 1.;
		FugProps[4] = 8.31451*Tk/Pbar;
		return -1;
	}
}

#endif





//=======================================================================================================
// Churakov-Gottschalk (CG) model for fluid mixtures
// References: Churakov and Gottschalk (2003a, 2003b)
// (c) SC June 2003
//=======================================================================================================


/// Generic constructor
TCGFcalc::TCGFcalc( long int NCmp, double Pp, double Tkp ):
    TSolMod( NCmp, 'F', Tkp, Pp )
{
    // Pparc = 0;
    phWGT = 0;
    aX = 0;
    // aGEX = 0;
    // aVol = 0;

    set_internal();
    alloc_internal();
}



TCGFcalc::TCGFcalc( SolutionData *sd, double *arphWGT, double *arX ):
                TSolMod( sd )
{
    Pparc = aPparc;
    phWGT = arphWGT;
    aX = arX;
    // aGEX = arGEX;
    // aVol = arVol;
    set_internal();
    alloc_internal();
}



// destructor
TCGFcalc::~TCGFcalc()
{
    free_internal();
}



/// set internally used parameters
void TCGFcalc::set_internal()
{
	PI_1 = 3.141592653589793120;  // pi
	TWOPI = 6.283185307179586230;  // 2.*pi
	PISIX = 0.523598775598298927;  // pi/6.
	TWOPOW1SIX = 1.12246204830937302;  // 2^ = 1/6)
	DELTA  = 0.00001;
	DELTAMOLLIM  = 0.0000001;
	R = 8.31439;  // R constant
	NA = 0.6023;
	P1 = 1.186892378996;
	PP2 = -0.4721963005527;
	P3 = 3.259515855283;
	P4 = 3.055229342609;
	P5 = 1.095409321023;
	P6 = 1.282306659774E-2;
	P7 = 9.55712461425E-2;
	P8 = 13.67807693107;
	P9 = 35.75464856619;
	P10 = 16.04724381643;
	AA1 = -0.120078459237;
	AA2 = -.808712488307;
	AA3 = .321543801337;
	A4 = 1.16965477132;
	A5 = -.410564939543;
	A6 = -.516834310691;
	BB1 = -2.18839961483;
	BB2 = 1.59897428009;
	BB3 = -.392578806128;
	B4 = -.189396607904;
	B5 = -.576898496254;
	B6 = -0.0185167641359;
	A00 = .9985937977069455;
	A01 = .5079834224407451;
	A10 = 1.021887697885469;
	A11 = -5.136619463333883;
	A12 = -5.196188074016755;
	A21 = -6.049240839050804;
	A22 = 18.67848155616692;
	A23 = 20.10652684217768;
	A31 = 9.896491419756988;
	A32 = 14.6738380473899;
	A33 = -77.44825116542995;
	A34 = -4.82871082941229;
}



void TCGFcalc::alloc_internal()
{
    paar = 0;
    paar1 = 0;
    FugCoefs =  0;
    EoSparam =  0;
    EoSparam1 = 0;
    Cf = new double [NComp][8];
}



void TCGFcalc::free_internal()
{
    if( paar ) delete paar;
        paar = 0;
    if( paar1 ) delete paar1;
        paar1 = 0;
    if( FugCoefs ) delete[]FugCoefs;
    if( EoSparam ) delete[]EoSparam;
    if( EoSparam1 ) delete[]EoSparam1;
    delete[]Cf;
}



/// high-level method to retrieve pure fluid fugacities
long int TCGFcalc::PureSpecies()
{
	double Fugacity = 0.1, Volume = 0.0;
	double X[1] = {1.};
	double Eos4parPT[4] = { 0.0, 0.0, 0.0, 0.0 },
            Eos4parPT1[4] = { 0.0, 0.0, 0.0, 0.0 } ;
	double roro;  // added, 21.06.2008 (TW)
	long int j, retCode = 0;

	for( j=0; j<NComp; j++)
	{
		// Calling CG EoS for pure fugacity
        if( Tk >= 273.15 && Tk < 1e4 && Pbar >= 1e-6 && Pbar < 1e5 )
        	retCode = CGFugacityPT( aDCc+j*NP_DC, Eos4parPT, Fugacity, Volume, Pbar, Tk, roro );
        else {
            Fugacity = Pbar;
            Volume = 8.31451*Tk/Pbar;
            Cf[j][0] = aDCc[j*NP_DC];
            if( Cf[j][0] < 1. || Cf[j][0] > 10. )
                Cf[j][0] = 1.;                 // foolproof temporary
            Cf[j][1] = aDCc[j*NP_DC+1];
            Cf[j][2] = aDCc[j*NP_DC+2];
            Cf[j][3] = aDCc[j*NP_DC+3];
            Cf[j][4] = 0.;
            Cf[j][5] = 0.;
            Cf[j][6] = 0.;
            Cf[j][7] = 0.;
            continue;
        }

        aGEX[j] = log( Fugacity / Pbar );  // now here (since 26.02.2008)  DK
        Pparc[j] = Fugacity;  // fure fluid fugacity (required for performance)
        aVol[j] = Volume*10.;  // molar volume of pure fluid component, J/bar to cm3

        // passing corrected EoS coeffs to calculation of fluid mixtures
        Cf[j][0] = Eos4parPT[0];
        if( Cf[j][0] < 1. || Cf[j][0] > 10. )
                Cf[j][0] = 1.;                            // foolproof temporary
        Cf[j][1] = Eos4parPT[1];
        Cf[j][2] = Eos4parPT[2];
        Cf[j][3] = Eos4parPT[3];

        CGFugacityPT( aDCc+j*NP_DC, Eos4parPT1, Fugacity, Volume, Pbar, Tk+Tk*GetDELTA(), roro );

        // passing corrected EoS coeffs for T+T*DELTA
        Cf[j][4] = Eos4parPT1[0];
        if( Cf[j][4] < 1. || Cf[j][4] > 10. )
                Cf[j][4] = 1.;                            // foolproof temporary
        Cf[j][5] = Eos4parPT1[1];
        Cf[j][6] = Eos4parPT1[2];
        Cf[j][7] = Eos4parPT1[3];

        // Calculation of departure functions
        CGResidualFunct( X, Eos4parPT, Eos4parPT1, 1, roro, Tk );  // changed, 21.06.2008 (TW)
    }  // j

	if ( retCode )
	{
		char buf[150];
		sprintf(buf, "CG fluid: calculation of pure fugacity failed");
		Error( "E71IPM IPMgamma: ",  buf );
	}
	return 0;
}



/// calculates T,P corrected binary interaction parameters
long int TCGFcalc::PTparam()
{
	long int i, j;

	if( FugCoefs ) delete[]FugCoefs;
	if( EoSparam ) delete[]EoSparam;
	if( EoSparam1 ) delete[]EoSparam1;

    FugCoefs = new double[ NComp ];
    EoSparam = new double[ NComp*4 ];
    EoSparam1 = new double[ NComp*4 ];

    PureSpecies();

    // Copying T,P corrected coefficients
    for( j=0; j<NComp; j++)
    {
    	for( i=0; i<4; i++)
                EoSparam[j*4+i] = Cf[j][i];
    	for( i=0; i<4; i++)
                EoSparam1[j*4+i] = Cf[j][i+4];
    }
    return 0;
}



/// high-level method to retrieve activity coefficients in the fluid mixture
long int TCGFcalc::MixMod()
{
	long int j;
	double roro; // changed, 21.06.2008 (TW)

	if( Tk >= 273.15 && Tk < 1e4 && Pbar >= 1e-6 && Pbar < 1e5 )
	{
		CGActivCoefPT( aX, EoSparam, FugCoefs, NComp, Pbar, Tk, roro );  // changed, 21.06.2008 (TW)
		if (roro <= 0. )
		{
			char buf[150];
			sprintf(buf, "CG fluid: bad calculation of density ro= %lg", roro);
			Error( "E71IPM IPMgamma: ",  buf );
		}

		// Phase volume of the fluid in cm3 (not needed any more?)
		phVOL[0] = phWGT[0] / roro;

	}

	else  // Setting Fugcoefs to 0 outside TP interval
		for( j=0; j<NComp; j++ )
			FugCoefs[ j ] = 0.0;

		for( j=0; j<NComp; j++  )
		{
			if( FugCoefs[j] > 1e-23 )
				lnGamma[j] = log(FugCoefs[j]/Pparc[j]);
			else
				lnGamma[j] = 0;
		}  // j
	return 0;
}



/// high-level method to calculate residual functions
long int TCGFcalc::ExcessProp( double *Zex )
{
	double roro; // changed, 21.06.2008 (TW)

	if( Tk >= 273.15 && Tk < 1e4 && Pbar >= 1e-6 && Pbar < 1e5 )
	{
		CGActivCoefPT( aX, EoSparam, FugCoefs, NComp, Pbar, Tk, roro );  // changed, 21.06.2008 (TW)
		if (roro <= 0. )
		{
			char buf[150];
			sprintf(buf, "CG fluid: bad calculation of density ro= %lg", roro);
			Error( "E71IPM IPMgamma: ",  buf );
		}

		// calculate residual functions
		CGResidualFunct( aX, EoSparam, EoSparam1, NComp, roro, Tk );

	}

	else  // setting residual functions to 0 outside TP interval
	{
		Grs = 0.;
		Srs = 0.;
		Hrs = 0.;
		CPrs = 0.;
		Vrs = 0.;
	}

	Ars = Grs - Vrs*Pbar;
	Urs = Hrs - Vrs*Pbar;

	// assignments (departure functions)
	Zex[0] = Grs;
	Zex[1] = Hrs;
	Zex[2] = Srs;
	Zex[3] = CPrs;
	Zex[4] = Vrs;
	Zex[5] = Ars;
	Zex[6] = Urs;

	return 0;
}



/// calculates ideal mixing properties
long int TCGFcalc::IdealProp( double *Zid )
{
	long int j;
	double s, sc, sp;

	s = 0.0;
	for (j=0; j<NComp; j++)
	{
		if ( x[j] > 1.0e-32 )
			s += x[j]*log(x[j]);
	}
	sc = (-1.)*R_CONST*s;
	sp = (-1.)*R_CONST*log(Pbar);
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
	Sid = sc + sp;
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



/// high-level method to retrieve pure fluid properties
long int TCGFcalc::CGFugacityPT( double *EoSparam, double *EoSparPT, double &Fugacity,
        double &Volume, double P, double T, double &roro )
{
	long int iRet = 0;
	// double ro;
	double X[1] = {1.};
	double FugPure[1];

	// modification to simplify CG database structure, 20.03.2007 (TW)
	EoSparPT[0] = EoSparam[0] + EoSparam[4]*exp(T*EoSparam[5]);
	EoSparPT[1] = EoSparam[1] + EoSparam[6]*exp(T*EoSparam[7]);
	EoSparPT[2] = EoSparam[2] + EoSparam[8]/(T+EoSparam[9]);
	EoSparPT[3] = EoSparam[3] + EoSparam[10]/(T+EoSparam[11]);

	// returns density
	CGActivCoefPT( X, EoSparPT, FugPure, 1, P, T, roro );  // changed, 21.06.2008 (TW)
	if( roro < 0.  )
	{
		return -1;
	};

	Fugacity = FugPure[0];
	roro = DENSITY( X, EoSparPT, 1, P, T );

	if( roro < 0 )
	{  // error - density could not be calculated
		iRet = -2;
		roro = 1.0;
	}

	Volume = 0.1/roro;  // in J/bar
	// roro = ro;  // added, 21.06.2008 (TW)

	return iRet;
}



long int TCGFcalc::CGActivCoefPT( double *X,double *param, double *act,
		   unsigned long int NN,   double Pbar, double T, double &roro )
{
	double *xtmp,*Fx;
	double P = Pbar/10.;
	xtmp = new double [NN];
	Fx = new double [NN];

	if(!paar)
		paar = new  EOSPARAM(X, param, NN);
	else
		paar->init( X, param, NN );

	double F0,Z,F1,fideal;
	double ro,delta = DELTA,ax,dx /*,tmp*/;
	long int i;

	norm(paar->XX0,paar->NCmp());
	copy(paar->XX0,xtmp,paar->NCmp());

	paar->ParamMix(xtmp);

	ro = ROTOTALMIX(P,T,paar);

	if( ro < 0.0 )  // Too low pressure, no corrections will be done
		return ( -1 );

	Z = P/(R*T*ro);
	F0 = FTOTALMIX(T,ro,paar);

	// fideal=log(R*T*ro/BARMPA);
	fideal = log(R*T*ro/0.1);
	ax = Z - 1.+fideal;

	for ( i=0;i<paar->NCmp();i++)
	{
		if ( xtmp[i]>0. )
		{
			copy(paar->XX0,xtmp,paar->NCmp());
			dx = xtmp[i]*delta;
			xtmp[i] += dx;
			norm(xtmp,paar->NCmp());
			paar->ParamMix(xtmp);
			F1 = FTOTALMIX(T,ro,paar)*(1.+dx);
			Fx[i] = (F1-F0)/(dx);
		}
		else Fx[i] = 0.;
	};

	// GMix=0.;
	for ( i=0;i<paar->NCmp();i++)
	{
		if ( xtmp[i]>0. && Fx[i]< 100. )
		{
			// tmp=log(paar.XX0[i]);
			// GMix+=tmp*paar.XX0[i];
			act[i] = exp(ax+Fx[i]);
		}
		else
		{
			act[i] = 0.;
		}
	};

	// GMix+=F0 + ax;
	// MLPutRealList(stdlink,act,paar.NCmp());
	delete[]xtmp;
	delete[]Fx;
	roro = ro;  // added, 21.06.2008 (TW)

	return 0;  // changed, 21.06.2008 (TW)
}



/// calculate residual functions through numerical derivative
long int TCGFcalc::CGResidualFunct( double *X, double *param, double *param1, unsigned long int NN,
		double ro, double T )
{
	double F0, Z, F1, vmix;
	double delta = DELTA;
	double *xtmp = new double [NN];

	if(!paar)
		paar = new  EOSPARAM(X, param, NN);
	else
		paar->init( X, param, NN );

	if(!paar1)
		paar1 = new  EOSPARAM(X, param1, NN);
	else
		paar1->init( X, param1, NN );

	norm(paar->XX0,paar->NCmp());
	norm(paar1->XX0,paar1->NCmp());
	copy(paar->XX0,xtmp,paar->NCmp());
	paar->ParamMix(xtmp);
	paar1->ParamMix(xtmp);
	Z = ZTOTALMIX(T,ro,paar);

	F0 = FTOTALMIX(T,ro,paar);
	// recalculate param1 for T+T*delta
	F1 = FTOTALMIX(T+T*delta,ro,paar1);
	// F1 = FTOTALMIX(T+T*delta,ro,paar);

	Srs = - ( (F1-F0)/(delta*Tk)*Tk + F0 ) * R_CONST;	// corrected, 20.06.2008 (TW)
	Hrs = (F0*Tk*R_CONST + Tk*Srs) + Z*R_CONST*Tk;
	Grs = Hrs - Tk*Srs;
	CPrs = 0.;
	vmix = Z*R_CONST*Tk/Pbar;
	Vrs = vmix;

    delete[]xtmp;
    return 0;

}



// void ACTDENS(double *data,long nn, double *act )
long int TCGFcalc::CGActivCoefRhoT( double *X, double *param, double *act,
		unsigned long int NN, double ro, double T )
{
	double   F0,Z,F1,GMix,fideal;
	double delta = DELTA,ax,dx,tmp;
	long int i;
	double *Fx,*xtmp;
	xtmp = new double [NN];
	Fx = new double [NN];

	if(!paar)
		paar = new EOSPARAM(X, param, NN);
		else
			paar->init( X, param, NN );

	norm(paar->XX0,paar->NCmp());
	copy(paar->XX0,xtmp,paar->NCmp());
	paar->ParamMix(xtmp);
	Z = ZTOTALMIX(T,ro,paar);
	F0 = FTOTALMIX(T,ro,paar);
	fideal = log(R*T*ro/0.1);
	ax = Z - 1.+fideal;

	for ( i=0;i<paar->NCmp();i++)
	{
		if ( xtmp[i]>0. )
		{
			copy(paar->XX0,xtmp,NN);
			if ( xtmp[i]>DELTAMOLLIM )
			{
				dx = xtmp[i]*delta;
			}
			else
			{
				dx = DELTAMOLLIM*delta;
			}

			xtmp[i] += dx;
			norm(xtmp,paar->NCmp());
			paar->ParamMix(xtmp);
			F1 = FTOTALMIX(T,ro,paar)*(1.+dx);
			Fx[i] = (F1-F0)/(dx);

		}
		else Fx[i] = 0.;
	};

	GMix = 0.;
	for ( i=0;i<paar->NCmp();i++)
	{
		if ( xtmp[i]>0. )
		{
			tmp = log(paar->XX0[i]);
			GMix += tmp*paar->XX0[i];
			act[i] = exp(ax+Fx[i]);
		}
		else
		{
			act[i] = 0.;
		}
	};

	delete[]xtmp;
	delete[]Fx;
	return 0;
    // MLPutRealList(stdlink,act,paar.NCmp());
   }



double TCGFcalc::DIntegral( double T, double ro, unsigned long int IType )
{
	static double TOld,roOld;
	static double a,b,c,d,e;
	static double data[][6]=
		{{-0.257431, 0.439229,  0.414783,  -0.457019, -0.145520,  0.299666},
		{-0.396724, 0.690721,  0.628935,  -0.652622, -0.201462, -0.23163 },
		{-0.488498, 0.863195,  0.761344,  -0.750086, -0.218562, -0.538463},
		{-0.556600, 0.995172,  0.852903,  -0.804710, -0.214736, -0.761700},
		{-0.611295, 1.103390,  0.921359,  -0.838804, -0.197999, -0.940714},
		{-0.657866, 1.196189,  0.975721,  -0.862346, -0.172526, -1.091678},
		{-0.698790, 1.278054,  1.020604,  -0.880027, -0.140749, -1.222733},
		{-0.735855, 1.351533,  1.058986,  -0.894024, -0.104174, -1.338626},
		{-0.769504, 1.418223,  1.092052,  -0.905347, -0.063730, -1.442391},
		{-0.800934, 1.479538,  1.121453,  -0.914864, -0.020150, -1.536070},
		{-0.829779, 1.535822,  1.147161,  -0.922381, 0.026157 , -1.621183},
		{-0.856655, 1.587957,  1.169885,  -0.928269, 0.074849 , -1.698853},
		{-0.881757, 1.636402,  1.190082,  -0.932668, 0.125590 , -1.769898},
		{-0.904998, 1.681421,  1.207610,  -0.935419, 0.178283 , -1.835070},
		{-0.926828, 1.723393,  1.223088,  -0.936667, 0.232649 , -1.894899},
		{-0.946773, 1.762571,  1.236007,  -0.936403, 0.288687 , -1.949858},
		{-0.965248, 1.799170,  1.246887,  -0.934650, 0.346207 , -2.000344}};

		// static double dt12[]=
		// {-2.139734,1.971553, 0.945513, -1.901492,-0.588630,-5.390941};
		// {-0.637684, 0.708107,  0.222086,  -0.481116, -0.332141, -3.492213};

	unsigned long int n;
	double *dtmp,rez;

	if ( (T!=TOld) || (ro!=roOld) )
	{
		TOld = T;
		roOld = ro;
		e = log(T);
		b = ro*ro;
		d = ro;
		c = ro*e;
		a = b*e;
	}

	// special case
	/*
	if ( IType==12 )
	{
		rez=(dt12[0]*T + dt12[1])*b +
		(dt12[2]*T + dt12[3])*ro + dt12[4]*T + dt12[5];
		return exp(rez);
	}
	*/

	n = IType-4;
	dtmp = data[n];
	rez = dtmp[0]*a + dtmp[1]*b + dtmp[2]*c + dtmp[3]*d + dtmp[4]*e + dtmp[5];
	return exp(rez);
}



double TCGFcalc::LIntegral( double T, double ro,unsigned long int IType )
{
	static double TOld,roOld;
	static double a,b,c,d,e;
	static double data[][6]=
	{{ -1.010391, 1.628552,  2.077476,  -2.30162 , -0.689931, -2.688117},
	{ -1.228611, 2.060090,  2.463396,  -2.453303, -0.573894, -3.350638},
	{ -1.354004, 2.402034,  2.718124,  -2.462814, -0.412252, -4.018632}};

	double *dtmp,rez;

	if ( (T!=TOld) || (ro!=roOld) )
	{
		TOld = T;
		roOld = ro;
		a = ro*ro*log(T);
		b = ro*ro;
		c = ro*log(T);
		d = ro;
		e = log(T);
	}

	switch ( IType )
	{
		case 662:
			dtmp = data[0];
			break;
		case 1262:
			dtmp = data[1];
			break;
		case 12122:
			dtmp = data[2];
			break;
		default:
			return 0;
	}
	rez = dtmp[0]*a + dtmp[1]*b + dtmp[2]*c + dtmp[3]*d + dtmp[4]*e + dtmp[5];
	return -exp(rez);
}



double TCGFcalc::KIntegral( double T, double ro,unsigned long int IType )
{
	static double TOld,roOld;
	static double a,b,c,d,e;
	static double data[][6]=
	{{ -1.050534, 1.747476,  1.749366,  -1.999227, -0.661046, -3.028720},
	{ -1.309550, 2.249120,  2.135877,  -2.278530, -0.773166, -3.704690},
	{ -1.490116, 2.619997,  2.404319,  -2.420706, -0.829466, -3.930928},
	{ -1.616385, 2.881007,  2.577600,  -2.484990, -0.828596, -4.175589},
	{ -1.940503, 3.552034,  2.940925,  -2.593808, -0.724353, -4.899975}};

	double *dtmp,rez;

	if ( (T!=TOld) || (ro!=roOld) )
	{
		TOld = T;
		roOld = ro;
		a = ro*ro*log(T);
		b = ro*ro;
		c = ro*log(T);
		d = ro;
		e = log(T);
	}

	switch ( IType )
	{
		case 222333:
			dtmp = data[0];
			break;
		case 233344:
			dtmp = data[1];
			break;
		case 334445:
			dtmp = data[2];
			rez = dtmp[0]*a + dtmp[1]*b + dtmp[2]*c + dtmp[3]*d + dtmp[4]*e + dtmp[5];
			return -exp(rez);
		case 444555:
			dtmp = data[3];
			break;
		case 666777:
			dtmp = data[4];
			break;
		default:
			return 0;

	}

	rez = dtmp[0]*a + dtmp[1]*b + dtmp[2]*c + dtmp[3]*d + dtmp[4]*e + dtmp[5];
	return exp(rez);
}



double TCGFcalc::K23_13( double T, double ro )
{
	static double TOld,roOld,KOLD;
	static double a,b,c,d,e;
	static double dtmp[]=
	{ -1.050534, 1.747476,  1.749366,  -1.999227, -0.661046, -3.028720};

	if ( (T!=TOld) || (ro!=roOld) )
	{
		TOld = T;
		roOld = ro;
		a = ro*ro*log(T);
		b = ro*ro;
		c = ro*log(T);
		d = ro;
		e = log(T);
	}
	else return KOLD;

	KOLD = dtmp[0]*a + dtmp[1]*b + dtmp[2]*c + dtmp[3]*d + dtmp[4]*e + dtmp[5];
	KOLD = exp(KOLD/3.);
	return KOLD;
}



double TCGFcalc::DENSITY( double *X, double *param, unsigned long NN, double Pbar, double T )
{
	double P = Pbar * 0.1;
	double *xtmp;
	double ro;

	xtmp = new double [NN];
	if( !paar1 )
		paar1 = new EOSPARAM(X,param,NN);
	else
		paar1->init( X, param, NN );

	norm(paar1->XX0,paar1->NCmp());
	copy(paar1->XX0,xtmp,paar1->NCmp());
	paar1->ParamMix(xtmp);
	ro = ROTOTALMIX(P,T,paar1);

	delete [] xtmp;
	if( ro < 0. )
		Error( ""," Error - density cannot be found at this T,P" );
	return ro;
}



double TCGFcalc::PRESSURE( double *X,double *param,
		unsigned long int NN,double ro, double T )
{
	double *xtmp;
	xtmp = new double [NN];

	if( !paar1 )
		paar1 = new EOSPARAM(X,param,NN);
	else
		paar1->init( X, param, NN );

	norm(paar1->XX0,paar1->NCmp());
	copy(paar1->XX0,xtmp,paar1->NCmp());
	paar1->ParamMix(xtmp);
	double P = PTOTALMIX(T,ro,paar1);
	delete [] xtmp;
	return P*10.;
}



void TCGFcalc::copy( double* sours,double *dest,unsigned long int num )
{
	unsigned long int i;
	for ( i=0; i<num; i++)
	{
		dest[i]=sours[i];
	};
}



void TCGFcalc::norm( double *X,unsigned long int mNum )
{
	double tmp=0.;
	unsigned long int i;
	for ( i=0; i<mNum; i++ )
	{
		tmp += X[i];
	}
	tmp = 1./tmp;
	for ( i=0; i<mNum; i++ )
	{
		X[i] *= tmp;
	}
}



double TCGFcalc::RPA( double beta,double nuw )
{
	double fi1,fi2;
	fi1 = (1.20110+(0.064890+(-76.860+(562.686+(-2280.090+(6266.840+(-11753.40+(14053.8
			+(-9491.490 +2731.030*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw;
	fi2 = (0.588890+(-7.455360+(40.57590+(-104.8970+(60.25470+(390.6310+(-1193.080
			+(1576.350+(-1045.910+283.7580*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw)*nuw*nuw;
	return  (-12.*fi1 + 192.*fi2*beta)*beta*beta/PI_1;
}



double TCGFcalc::dHS( double beta,double ro )
{
	// service constants
	double DV112 = 1./12.;
	double DV712 = 7./12.;
	// local variables
	double T12, T112, T712, B13, dB, delta, d;
	double a0, a1, a6, a3, a4, a7, a9, a12;
	double p0, p2, p6, p3, p5, p8, p11;
	double dbdl, ri6ro, ri6ro2, d3, d2, dnew, F0, F1;
	unsigned long int i;

	T12 = sqrt(beta);
	T112 = exp(DV112*log(beta));
	T712 = exp(DV712*log(beta));
	B13 = (1+beta);
	B13 = B13*B13*B13;

	dB = (P1*T112+PP2*T712+(P3+(P4+P5*beta)*beta)*beta)/B13;
	delta = (P6+P7*T12)/(1.+(P8+(P9+P10*T12)*T12)*T12);

	dbdl = dB*delta;
	ri6ro = PISIX*ro;
	ri6ro2 = ri6ro*ri6ro;

	a0 = dB+dbdl;
	a1 = -1.;
	a3 = (-1.5*dB -3.75*dbdl)*ri6ro;
	a4 = (1.5*ri6ro);
	a6 = (2.*dB + dbdl)*0.25*ri6ro2;
	a7 = -0.5*ri6ro2;
	a9 = -2.89325*ri6ro2*ri6ro*dbdl;
	a12 = -0.755*ri6ro2*ri6ro2*dbdl;

	p0 = -1.;
	p2 = a3*3.;
	p3 = a4*4.;
	p5 = a6*6.;
	p6 = a7*7.;
	p8 = a9*9.;
	p11 = a12*12.;

	d = dB;
	i = 0;

	while ( i++<21 )
	{
		d2 = d*d;
		d3 = d*d*d;
		F0 = a0+(a1+(a3+(a4+(a6+(a7+(a9+a12*d3)*d2)*d)*d2)*d)*d2)*d;
		F1 = p0+(p2+(p3+(p5+(p6+(p8+p11*d3)*d2)*d)*d2)*d)*d2;
		dnew = d-F0/F1;
		if ( fabs(dnew-d)<1.E-7 )
		{
			return dnew;
		}
		d = dnew;
	}

	if ( i>=20 )
	{
		return dB;
	}
	return dnew;
}



double TCGFcalc::FWCA( double T,double ro )
{
	static double TOld,roOld,F;
	double d,beta,nu,nuw;
	double nu1w1,nu1w2,nu1w3,nu1w4,nu1w5;
	double a0,a1,a2,a3;
	double I2;
	double I1_6,I1_12;
	double dW,dW12,dW6;
	double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7;
	double F0,F1,FA;
	double rm,rmdw1,rmdw2,rmdw3,rmdw4,rmdw5;

	if ((T==TOld) && (ro==roOld))
	{
		return F;
	}
	else
	{
		TOld = T;
		roOld = ro;
	}

	rm = TWOPOW1SIX;
	beta = 1./T;
	d = dHS( beta, ro );
	tmp2 = PISIX*d*d*d;
	nu = tmp2*ro;
	tmp1 = (1. - nu/16.);
	nuw = nu*tmp1;
	dW = d*exp(1./3.*log(tmp1));

	nu1w1 = (1.-nuw);
	nu1w2 = nu1w1*nu1w1;
	nu1w3 = nu1w2*nu1w1;
	nu1w4 = nu1w2*nu1w2;
	nu1w5 = nu1w2*nu1w3;

	tmp1 = (1-nu);
	tmp1 = tmp1*tmp1;
	F0 = ((4.-3.*nu)*nu)/tmp1;

	a0 = fa0( nuw , nu1w2);
	a1 = fa1( nuw , nu1w3);
	a2 = fa2( nuw , nu1w4);
	a3 = fa3( nuw , nu1w5);

	I1_6 = fI1_6( nuw );
	I1_12 = fI1_12( nuw );

	rmdw1 = rm/dW;
	rmdw2 = rmdw1*rmdw1;
	rmdw3 = rmdw1*rmdw2;
	rmdw4 = rmdw2*rmdw2;
	rmdw5 = rmdw3*rmdw2;

	dW6 = dW*dW*dW;
	dW6 = 1./(dW6*dW6);
	dW12 = dW6*dW6;

	tmp1 = (a0/4.+ a1/12. + a2/24. + a3/24.)*dW6;
	tmp2 = (a0/10.+ a1/90. + a2/720. + a3/5040.)*(-dW12);
	tmp3 = (a0 - a1/3. + a2/12 - a3/60)/8.;
	tmp4 = (a0 - a1 + a2/2. - a3/6.)*rmdw2*(-9.)/40.;
	tmp5 = (a1 - a2 + a3/2)*rmdw3*(-2.)/9.;
	tmp6 = (a2 - a3)*rmdw4*(-9.)/64.;
	tmp7 = a3*(-3.)/35.*rmdw5;

	I2 = tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7;

	F1 = 48.*nuw*(I1_12*dW12-I1_6*dW6 + I2)*beta;
	FA = RPA(beta,nuw);

	F = F0+F1+FA;

	return F;
}



double TCGFcalc::ZWCANum( double T,double ro )
{
	double delta = DELTA;
	double a0,a1;
	a1 = FWCA(T,ro*(1.+delta));
	a0 = FWCA(T,ro);
	return 1.+(a1-a0)/delta;
}



double TCGFcalc::UWCANum( double T,double ro )
{
	double delta = DELTA;
	double a0,a1,beta0,beta1;
	beta0 = 1./T;
	beta1 = beta0*(1.+delta);
	a1 = FWCA(1./beta1,ro);
	a0 = FWCA(T,ro);
	return (a1-a0)/(beta1-beta0);
}



double TCGFcalc::FDipPair( double T,double ro,double m2 )
{
	double kappa,Z,U,beta,F;
	kappa = m2*m2/(24.*T);
	beta = 1./T;
	Z = ZWCANum(T,ro);
	U = UWCANum(T,ro);
	F = kappa*(4.*beta*U-Z+1.);
	return F;
}



double TCGFcalc::J6LJ( double T,double ro )
{
	double kappa,Z,U,beta,F;
	beta = 1./T;
	Z = ZWCANum(T,ro);
	kappa = -16.*PI_1*ro*beta;
	U = UWCANum(T,ro);
	F = (4.*beta*U-Z+1.)/kappa;
	return F;
}



double TCGFcalc::FTOTALMIX( double T_Real,double ro_Real,EOSPARAM* param )
{
	double FF,A0,A2,A3,AP,A1;
	// unsigned iall,inopol;
	double emix,s3mix,rotmp,T2R;
	double Jind,Jdp;
    long int /*itmp,jtmp,ktmp,*/ i,j,k;
    double s3tmp,mtmp,IK /*,atmp*/;
    double imtmp,jmtmp,iatmp,jatmp;
    double m2i,m2j,m2k;
    double s3tmpij,s3tmpik,s3tmpjk;
    double IKtmpij,IKtmpik,IKtmpjk;

    // iall=param.inonzero();
    // inopol=param.inonpolar();
    emix = param->EMIX();
    s3mix = param->S3MIX();

      rotmp = NA*ro_Real;
      T2R = T_Real*T_Real;

      A0 = FWCA(T_Real/emix,s3mix*rotmp);
     // if ( inopol< iall )
      {
        // dipole part
        A2 = 0.;
        for ( i=0; i<param->NCmp()-1; i++ )
        {
          for ( j=i+1; j<param->NCmp(); j++ )
          {
            s3tmp = param->MIXS3(i,j);
            Jdp = J6LJ(T_Real*s3tmp/param->MIXES3(i,j) , s3tmp*rotmp);
            A2 += param->M2R(i)*param->M2R(j)*Jdp*
                           param->X(i)*param->X(j)/s3tmp;
          }
        }
          A2 *= 2.;
          for ( i=0; i<param->NCmp(); i++ )
          {
            // itmp=param->ind(i);
            mtmp = param->M2R(i);
            s3tmp = param->SIG3(i);
            Jdp = J6LJ(T_Real/param->EPS(i),s3tmp*rotmp);
            A2 += mtmp*mtmp*Jdp*param->X(i)*param->X(i)/s3tmp;
          }

         A2 = -A2*TWOPI*rotmp/(3.*T2R);
         // A2 done

         if ( A2!=0. )
         {

          A3 = 0.;
          for ( i=0; i<param->NCmp(); i++ )
          {
            // itmp=param->ind(i);
            m2i = param->M2R(i);

            for ( j=0; j<param->NCmp(); j++  )
            {
             // jtmp=param->ind(j);
              m2j = param->M2R(j);

              s3tmpij = param->MIXS3(i,j);
              IKtmpij = K23_13(T_Real*s3tmpij/param->MIXES3(i,j),
                                                    s3tmpij*rotmp);
              for ( k=0; k<param->NCmp(); k++  )
              {
               // ktmp=param->ind(k);
               m2k = param->M2R(k);

               s3tmpik = param->MIXS3(i,k);
               s3tmpjk = param->MIXS3(j,k);

               IKtmpik = K23_13(T_Real*s3tmpik/param->MIXES3(i,k),
                                                    s3tmpik*rotmp);
               IKtmpjk = K23_13(T_Real*s3tmpjk/param->MIXES3(j,k),
                                                    s3tmpjk*rotmp);

               IK = IKtmpij*IKtmpik*IKtmpjk;
               A3 += m2i*m2j*m2k*IK*pow(s3tmpij*s3tmpik*s3tmpjk,-1./3.)*
               param->X(i)*param->X(j)*param->X(k);
              }
            }
          }
            A3 = A3*32.*sqrt(14.*PI_1/5.)*
                  rotmp*rotmp*PI_1*PI_1*PI_1/(135.*T_Real*T2R);
            AP = A2/(1. - A3/A2);
         }
         else AP = 0.;

        // induced interaction
        A1 = 0.;
        for ( i=0; i<param->NCmp(); i++ )
        {
         // itmp=param->ind(i);
          iatmp = param->A(i);
          imtmp = param->M2R(i);
          for ( j=0; j<param->NCmp(); j++ )
          {
            // jtmp=param->ind(j);
            jatmp = param->A(j);
            jmtmp = param->M2R(j);

            s3tmp = param->MIXS3(i,j);
            Jind = J6LJ(T_Real*s3tmp/param->MIXES3(i,j),s3tmp*rotmp);

           A1 += (iatmp*jmtmp + jatmp*imtmp)
                  *Jind*param->X(i)*param->X(j)/s3tmp;
          }
        }
        A1 = -A1*TWOPI*rotmp/T_Real;
// A1=-A1*FOURPI*rotmp/T_Real;
// A1=0.;

      }  // end of polar contribution

     FF = A0 + A1 + AP;
     //printf("%g %g %g %g %g",A0,A1,A2,A3,AP);
     //exit(1) ;

    return FF;
  }



double TCGFcalc::UTOTALMIX( double T_Real,double ro_Real,EOSPARAM* param )
{
  double T /*,ro,s3 */;
  double delta = DELTA;
  double a0,a1,beta0,beta1,eps;
  eps = param->EMIX();
  T = T_Real/eps;

  beta0 = 1./T;
  beta1 = beta0*(1.+delta);
  a1 = FTOTALMIX((1./beta1)*eps,ro_Real,param);
  a0 = FTOTALMIX(T_Real,ro_Real,param);
  return (a1-a0)/(beta1-beta0);
 }



double TCGFcalc::ZTOTALMIX( double T_Real,double ro_Real,EOSPARAM* param )
 {
  double delta = DELTA;
  double a0,a1;
  a1 = FTOTALMIX(T_Real,ro_Real*(1.+delta),param);
  a0 = FTOTALMIX(T_Real,ro_Real,param);

  return 1.+(a1-a0)/delta;
 }



double TCGFcalc::PTOTALMIX( double T_Real,double ro_Real,EOSPARAM* param )
 {
  double Z;
    Z = ZTOTALMIX(T_Real,ro_Real,param);
    return Z*R*T_Real*ro_Real;
 }



/// melting density
double TCGFcalc::Melt( double T )
 {

    return T*0.+.9;

 }



double TCGFcalc::Melt2(double T)
 {
    return T*0.+3.;
 }


 #define FIRSTSEED (15)
 #define ROMIN (1.E-2)
 #define NPOINT (5)



void  TCGFcalc::choose( double *pres, double P,unsigned long int &x1,unsigned long int &x2 )
 {
  unsigned long int i;
  double deltam = -10000000.,tmp;
  double deltap = 10000000.;

  for ( i=0; i<NPOINT; i++ )
  {
    tmp = P-pres[i];

    if ( tmp>0. )
    {
       if ( tmp<deltap )
       {
        deltap=tmp;
        x1 = i;
       }
    }
    else
    {
       if ( tmp>deltam )
       {
        deltam=tmp;
        x2 = i;
       }
    }
  }
     return ;
 }



double TCGFcalc::ROTOTALMIX( double P,double TT,EOSPARAM* param )
 {
     unsigned long int i;
     double T /*,ro*/;
     double fact, fact0, romax, dro, roarr[FIRSTSEED];
     double Ptmp[FIRSTSEED], ro0, ro1, rotest, PP0, PP1 /* ,Ptest */;
     double a,b;
     double inttofloat;
     double f[4],x[4],ff,dens[5],pres[5];
     unsigned long int x1=0L,x2=0L;
// double ptmp;

     T = TT/param->EMIX();
     fact0 = 1./(param->S3MIX()*NA);
     fact = R*TT*fact0;

     romax = Melt(T);
     inttofloat = FIRSTSEED-1;
     dro = (romax-ROMIN)/inttofloat;
     roarr[0] = ROMIN;
     roarr[1] = 2.*ROMIN;

     for ( i=2; i<FIRSTSEED; i++)
     {
       inttofloat = i;
       roarr[i] = ROMIN+inttofloat*dro;
     }

     for ( i=0; i<FIRSTSEED; i++)
     {
      Ptmp[i] = ZTOTALMIX(TT,roarr[i]*fact0,param);
      Ptmp[i] *= roarr[i] * fact;
      if ( Ptmp[i] > P )
      {
        break;
      }
     }

     if ( i==0 )  // Uses aproximation of ideal gas
     {
            return P/(R*TT);
     }

     // additional high pressure inteval
     if ( i==FIRSTSEED )
     {

     // roarr[0]=romax-0.0001;
     roarr[0] = roarr[FIRSTSEED-1];
     Ptmp[0] = Ptmp[FIRSTSEED-1];

     romax = Melt2(T);
     inttofloat = FIRSTSEED-1;
     dro = (romax-ROMIN)/inttofloat;
     for ( i=1; i<FIRSTSEED; i++)
     {
       inttofloat = i;
       roarr[i] = ROMIN+inttofloat*dro;
     }

     for ( i=1; i<FIRSTSEED; i++)
     {
      Ptmp[i] = ZTOTALMIX(TT,roarr[i]*fact0,param)*roarr[i]*fact;
      if ( Ptmp[i]>P )
      {
        break;
      }
     }

     if ( i==FIRSTSEED || i==0 )
     {
         printf( "Input pressure is too high!\n" );
            // exit(1);
         return (-1.0);
     }
     }

     ro0 = roarr[i-1];
     ro1 = roarr[i];
     PP0 = Ptmp[i-1];
     PP1 = Ptmp[i];
     i = 0;

   while ( i++<20 )
   {
     // Start interp
     ff = ro0;
     dens[0] = ro0;
     dens[1] = ro1;
     pres[0] = PP0;
     pres[1] = PP1;

     // first order
     x[0] = P-pres[0];
     f[0] = (dens[1]-dens[0])/(pres[1]-pres[0]);
     ff += f[0]*x[0];

     // second order
     dens[2] = ff;
     pres[2] = ZTOTALMIX(TT,ff*fact0,param)*ff*fact;

     if ( fabs(pres[2]-P)<1E-5 )
     {
       return ff*fact0;
     }

     x[1] = x[0]*(P-pres[1]);
     f[1] = (dens[2]-dens[1])/(pres[2]-pres[1]);

     f[0] = (f[1]-f[0])/(pres[2]-pres[0]);
     ff += f[0]*x[1];

     // third order
     dens[3] = ff;
     pres[3] = ZTOTALMIX(TT,ff*fact0,param)*ff*fact;
     if ( fabs(pres[3]-P)<1E-6 )
     {
      return ff*fact0;
     }
     x[2] = x[1]*(P-pres[2]);
     f[2] = (dens[3]-dens[2])/(pres[3]-pres[2]);
     f[1] = (f[2]-f[1])/(pres[3]-pres[1]);
     f[0] = (f[1]-f[0])/(pres[3]-pres[0]);
     ff += f[0]*x[2];
     dens[4] = ff;
     pres[4] = ZTOTALMIX(TT,ff*fact0,param)*ff*fact;
     if ( fabs(pres[4]-P)<1e-6 )
     {
      return ff*fact0;
     }

     choose(pres,P,x1,x2);

     ro0 = dens[x1];
     ro1 = dens[x2];
     PP0 = pres[x1];
     PP1 = pres[x2];

      if ( fabs((ro1-ro0))<0.001 )
      {
          a = (PP1-PP0)/(ro1-ro0);
          b = PP1-a*ro1;
          rotest = (P-b)/a;
          return rotest*(fact0);
      }
   }
        //return 10.;
         /// bad result

          a = (PP1-PP0)/(ro1-ro0);
          b = PP1-a*ro1;
          rotest = (P-b)/a;
          return rotest*(fact0);

 }

#ifndef IPMGEMPLUGIN



/// calculates properties of pure fluids when called from DCthermo
long int TCGFcalc::CGcalcFugPure( double Tmin, float *Cemp, double *FugProps )
{
	long int retCode = 0;
	double T, P, Fugacity = 0.1, Volume = 0.0;
	double X[1] = {1.};
	double roro = 1.;  // added 21.06.2008 (TW)
	double Coeff[12];  // MAXEOSPARAM = 20;
	double Eos4parPT[4] = { 0.0, 0.0, 0.0, 0.0 },
		Eos4parPT1[4] = { 0.0, 0.0, 0.0, 0.0 } ;

	T = Tk;
	P = Pbar;

	for(long int ii=0; ii<12; ii++ )
		Coeff[ii] = (double)Cemp[ii];

	// Calling CG EoS functions here
	if( (Tk >= Tmin) && (Tk < 1e4) && (Pbar >= 1e-6) && (Pbar < 1e5) )
	{
		retCode = CGFugacityPT( Coeff, Eos4parPT, Fugacity, Volume, P, T, roro );
		FugProps[0] = Fugacity/Pbar;
		FugProps[1] = 8.31451 * Tk * log( Fugacity / P );
		FugProps[4] = Volume;
		retCode = CGFugacityPT( Coeff, Eos4parPT1, Fugacity, Volume, P, T+T*DELTA, roro );
		CGResidualFunct( X, Eos4parPT, Eos4parPT1, 1, roro, T );
		FugProps[2] = Hrs;
		FugProps[3] = Srs;
		return retCode;
	}

	else
	{
		for( int i=1; i<6; i++ )
			FugProps[i] = 0.;
		FugProps[0] = 1.;
		FugProps[4] = 8.31451*Tk/Pbar;
		return -1;
	}
}

#endif




//=======================================================================================================
// Implementation of EOSPARAM class (used by TCGFcalc class)
//=======================================================================================================


void EOSPARAM::free()
{
	long int i;

	if ( NComp > 0)
	{
		for ( i=0;i<NComp;i++ )
			delete[]mixpar[i];
		delete[]mixpar;

		delete[]epspar;
		delete[]sig3par;
		delete[]XX;
		delete[]eps;
		delete[]eps05;
		delete[]sigpar;
		delete[]mpar;
		delete[]apar;
		delete[]aredpar;
		delete[]m2par;
		delete[]XX0;
		NComp = 0;
	}
}



void EOSPARAM::allocate()
{
	long int i;

	mixpar = new double*[NComp];
	for ( i=0; i<NComp; i++ )
		mixpar[i] = new double[NComp];

	epspar = new double[NComp];
	sig3par = new double[NComp];
	XX = new double[NComp];
	eps = new double[NComp];
	eps05 = new double[NComp];
	sigpar = new double[NComp];
	mpar = new double[NComp];
  	apar = new double[NComp];
	aredpar = new double[NComp];
	m2par = new double[NComp];
	XX0 = new double[NComp];
}



void EOSPARAM::init( double *Xinp, double * data, long int nn )
{
	long int i,j;
	double tmp;

	if( nn != NComp )
	{ // or error message
		free();
		NComp = nn;
		allocate();
	}

	for ( i=0;i<NComp;i++ )
	{
		XX0[i] = Xinp[i];

		sigpar[i] = data[i*4 ];
		eps[i] = data[i*4 + 1];
		mpar[i] = data[i*4 + 2];
		apar[i] = data[i*4 + 3];
	}
	for ( i=0; i<NComp; i++ )
	{
		tmp = sigpar[i];
		tmp = tmp*tmp*tmp;
		sig3par[i] = tmp;
		eps05[i] = sqrt(eps[i]);
		epspar[i] = tmp*eps[i];
		m2par[i] = mpar[i]*mpar[i]/(1.38048E-4);
		aredpar[i] = apar[i]/tmp;
	}

	// calculation of mixing properties
	for ( i=0; i<NComp-1; i++ )
	{
		for ( j=i+1; j<NComp; j++ )
		{
			tmp = (sigpar[i]+sigpar[j])*0.5;
			tmp = tmp*tmp*tmp;
			mixpar[i][j] = tmp;
			mixpar[j][i] = tmp*eps05[i]*eps05[j];
		}
	}
}


long int EOSPARAM::ParamMix( double *Xin )
  {
    long int j,i;
    double tmp,tmp1,tmp2;

    for ( i=0; i<NComp; i++ )
    	XX[i] = Xin[i];

    emix = 0.;
    s3mix = 0.;
    for ( i=0; i<NComp-1; i++ )
    {
      for ( j=i+1; j<NComp; j++ )
      {
          tmp = XX[i]*XX[j];
          tmp2 = mixpar[j][i];  //eps
          tmp1 = mixpar[i][j];  //signa
          s3mix += tmp1*tmp;
          emix += tmp2*tmp;
      }
    }
    s3mix *= 2.;
    emix *= 2.;
    for ( i=0; i<NComp; i++ )
    {
          tmp = XX[i]*XX[i];

          s3mix += sig3par[i]*tmp;
          emix += epspar[i]*tmp;
    }
    emix = emix/s3mix;
    return NComp;
  }



//=======================================================================================================
// Soave-Redlich-Kwong (SRK) model for fluid mixtures
// References: Soave (1972), Soave (1993) 
// (c) TW December 2008
//=======================================================================================================


// constructor
TSRKcalc::TSRKcalc( long int NCmp, double Pp, double Tkp ):
    TSolMod( NCmp, 'E',  Tkp, Pp )
{
    // aGEX = 0;
    // aVol = 0;
    Pparc = 0;
    alloc_internal();
}



TSRKcalc::TSRKcalc( SolutionData *sd ):
                TSolMod( sd )
{
    Pparc = aPparc;
    // aGEX = arGEX;
    // aVol = arVol;
    alloc_internal();
}



TSRKcalc::~TSRKcalc()
{
    free_internal();
}



/// allocate work arrays for pure fluid and fluid mixture properties
void TSRKcalc::alloc_internal()
{
	Eosparm = new double [NComp][4];
	Pureparm = new double [NComp][4];
	Fugpure = new double [NComp][6];
	Fugci = new double [NComp][4];
	a = new double *[NComp];
	b = new double *[NComp];
	KK = new double *[NComp];
	dKK = new double *[NComp];
	d2KK = new double *[NComp];
	AA = new double *[NComp];

	for (long int i=0; i<NComp; i++)
	{
		a[i] = new double[NComp];
		b[i] = new double[NComp];
		KK[i] = new double[NComp];
		dKK[i] = new double[NComp];
		d2KK[i] = new double[NComp];
		AA[i] = new double[NComp];
	}
}



void TSRKcalc::free_internal()
{
	long int i;

	for (i=0; i<NComp; i++)
	{
		delete[]a[i];
		delete[]b[i];
		delete[]KK[i];
		delete[]dKK[i];
		delete[]d2KK[i];
		delete[]AA[i];
	}

	delete[]Eosparm;
	delete[]Pureparm;
	delete[]Fugpure;
	delete[]Fugci;
	delete[]a;
	delete[]b;
	delete[]KK;
	delete[]dKK;
	delete[]d2KK;
	delete[]AA;

}



/// high-level method to retrieve pure fluid fugacities
long int TSRKcalc::PureSpecies()
{
    long int j, retCode = 0;

    for( j=0; j<NComp; j++)
    {
        // Calling SRK EoS for pure fugacity
        retCode =  FugacityPT( j, aDCc+j*NP_DC );
        aGEX[j] = log( Fugpure[j][0] );
        Pparc[j] = Fugpure[j][0]*Pbar;  // fure fluid fugacity (required for performance)
        aVol[j] = Fugpure[j][4]*10.;  // molar volume of pure fluid component, J/bar to cm3
    } // j

    if ( retCode )
    {
        char buf[150];
        sprintf(buf, "SRK fluid: calculation of pure fugacity failed");
                                Error( "E71IPM IPMgamma: ",  buf );
    }

    return 0;
}



/// high-level method to calculate T,P corrected binary interaction parameters
long int TSRKcalc::PTparam()
{
	long int j, i;

	PureSpecies();

	// set all interaction parameters zero
	for( j=0; j<NComp; j++ )
	{
		for( i=0; i<NComp; i++ )
		{
			KK[j][i] = 0.;
			dKK[j][i] = 0.;
			d2KK[j][i] = 0.;
		}
	}

	switch ( MixCode )
	{
        case MR_UNDEF_:
        case MR_WAAL_:
			MixingWaals();
			break;
		case MR_CONST_:
			MixingConst();
			break;
		case MR_TEMP_:
			MixingTemp();
			break;
		default:
			break;
	}

	return 0;
}



/// high-level method to retrieve activity coefficients of the fluid mixture
long int TSRKcalc::MixMod()
{
	long int j, iRet;

	iRet = FugacitySpec( Pparc );
	phVOL[0] = PhVol * 10.;

    for(j=0; j<NComp; j++)
    {
        if( Fugci[j][3] > 1e-23 )
        	lnGamma[j] = log( Fugci[j][3] );
        else
        	lnGamma[j] = 0;
    }
    if ( iRet )
    {
    	char buf[150];
    	sprintf(buf, "SRK fluid: calculation failed");
			Error( "E71IPM IPMgamma: ",  buf );
    }
    return iRet;
}



/// high-level method to retrieve residual functions of the fluid mixture
long int TSRKcalc::ExcessProp( double *Zex )
{
	long int iRet;

	iRet = ResidualFunct( Pparc );

    if ( iRet )
    {
    	char buf[150];
    	sprintf(buf, "SRK fluid: calculation failed");
    	Error( "E71IPM IPMgamma: ",  buf );
    }

	Ars = Grs - Vrs*Pbar;
	Urs = Hrs - Vrs*Pbar;

	// assignments (residual functions)
	Zex[0] = Grs;
	Zex[1] = Hrs;
	Zex[2] = Srs;
	Zex[3] = CPrs;
	Zex[4] = Vrs;
	Zex[5] = Ars;
	Zex[6] = Urs;

	return iRet;

}



/// basic van der waals mixing rule
long int TSRKcalc::MixingWaals()
{
	// currently no calculations

	return 0;
}



/// constant one-term interaction parameter
long int TSRKcalc::MixingConst()
{
	long int ip, i1, i2;
	double k, dk, d2k;

	if( NPcoef > 0 )
	{
		// transfer interaction parameters that have non-standard value
		for ( ip=0; ip<NPar; ip++ )
		{
			i1 = aIPx[MaxOrd*ip];
			i2 = aIPx[MaxOrd*ip+1];
			k = aIPc[NPcoef*ip];
			dk = 0.;
			d2k = 0.;
			KK[i1][i2] = k;
			dKK[i1][i2] = dk;
			d2KK[i1][i2] = d2k;
			KK[i2][i1] = k;   // symmetric case
			dKK[i2][i1] = dk;
			d2KK[i2][i1] = d2k;
		}
	}

	return 0;
}



/// temperature dependent one-term interaction parameter
long int TSRKcalc::MixingTemp()
{
	long int i, j, ip, i1, i2;
	double ai, aj, bi, bj, di, dj, dai, daj, d2ai, d2aj, ddi, ddj, d2di, d2dj,
				U, V, dU, dV, d2U, d2V, tmp, k, dk, d2k, C;

	// set model specific interaction parameters zero
	for( j=0; j<NComp; j++ )
	{
		for( i=0; i<NComp; i++ )
		{
			a[j][i] = 0.;
			b[j][i] = 0.;
		}
	}

	if( NPcoef > 0 )
	{
		// transfer parameters that have non-standard value
		for ( ip=0; ip<NPar; ip++ )
		{
			i1 = aIPx[MaxOrd*ip];
			i2 = aIPx[MaxOrd*ip+1];
			a[i1][i2] = aIPc[NPcoef*ip];
			b[i1][i2] = aIPc[NPcoef*ip+1];
			a[i2][i1] = aIPc[NPcoef*ip];  // symmetric case
			b[i2][i1] = aIPc[NPcoef*ip+1];
		}
	}

	// calculate binary interaction parameters
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
			if ( a[i][j] == 0.0 )
				tmp = 1.0;
			else
				tmp = a[i][j];

			// read a, b, da, d2a
			ai = Pureparm[i][0];
			aj = Pureparm[j][0];
			bi = Pureparm[i][1];
			bj = Pureparm[j][1];
			dai = Pureparm[i][2];
			daj = Pureparm[j][2];
			d2ai = Pureparm[i][3];
			d2aj = Pureparm[j][3];

			// calculate k and derivatives
			di = sqrt(ai)/bi;
			dj = sqrt(aj)/bj;
			ddi = (0.5/bi) * pow(ai,-0.5) * dai;
			ddj = (0.5/bj) * pow(aj,-0.5) * daj;
			d2di = (0.5/bi) * ( (-0.5)*pow(ai,-1.5)*dai*dai + pow(ai,-0.5)*d2ai );
			d2dj = (0.5/bj) * ( (-0.5)*pow(aj,-1.5)*daj*daj + pow(aj,-0.5)*d2aj );

			C = ( b[i][j]/tmp - 1. );
			U = a[i][j]*pow((298.15/Tk),C) - pow((di-dj),2.);
			V = (2.*di*dj);
			dU = - ( a[i][j]*C*pow((298.15/Tk),C) ) / Tk - 2.*(di-dj)*(ddi-ddj);
			dV = 2.*( ddi*dj + di*ddj );
			d2U = ( a[i][j]*pow(C,2.)*pow((298.15/Tk),C) ) / pow(Tk,2.)
						+ ( a[i][j]*C*pow((298.15/Tk),C) ) / pow(Tk,2.)
						- 2.*( pow((ddi-ddj),2.) + (di-dj)*(d2di-d2dj) );
			d2V = 2.*( d2di*dj + 2.*ddi*ddj + di*d2dj );
			k = U/V;
			dk = (dU*V-U*dV)/pow(V,2.);
			d2k = (d2U*V+dU*dV)*pow(V,2.)/pow(V,4.) - (dU*V)*(2.*V*dV)/pow(V,4.)
						- (dU*dV+U*d2V)*pow(V,2.)/pow(V,4.) + (U*dV)*(2.*V*dV)/pow(V,4.);

			// assignments
			KK[i][j] = k;
			dKK[i][j] = dk;
			d2KK[i][j] = d2k;
		}
	}

	return 0;
}



/// calculates ideal mixing properties
long int TSRKcalc::IdealProp( double *Zid )
{
	long int j;
	double s, sc, sp;

	s = 0.0;
	for (j=0; j<NComp; j++)
	{
		if ( x[j] > 1.0e-32 )
			s += x[j]*log(x[j]);
	}
	sc = (-1.)*R_CONST*s;
	sp = (-1.)*R_CONST*log(Pbar);
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
	Sid = sc + sp;
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



/// High-level method to retrieve pure fluid properties
long int TSRKcalc::FugacityPT( long int i, double *EoSparam )
{
	long int iRet = 0;
	double Tcrit, Pcrit, omg, N, apure, bpure, da, d2a;

	// reads EoS parameters from database into work array
	if( !EoSparam )
		return -1;  // Memory alloc error

	Eosparm[i][0] = EoSparam[0];  // critical temperature in K
	Eosparm[i][1] = EoSparam[1];  // critical pressure in bar
	Eosparm[i][2] = EoSparam[2];  // Pitzer acentric factor omega
	Eosparm[i][3] = EoSparam[3];  // empirical EoS parameter N
	Tcrit = Eosparm[i][0];
	Pcrit = Eosparm[i][1];
	omg = Eosparm[i][2];
	N = Eosparm[i][3];

	AB( Tcrit, Pcrit, omg, N, apure, bpure, da, d2a );

	Pureparm[i][0] = apure;  // a parameter
	Pureparm[i][1] = bpure;  // b parameter
	Pureparm[i][2] = da;  // da/dT
	Pureparm[i][3] = d2a;  // d2a/dT2

	iRet = FugacityPure( i );
    if( iRet)
    	return iRet;

    return iRet;
}



/// Calculates attractive (a) and repulsive (b) parameter of SRK equation of state
/// and partial derivatives of alpha function
long int TSRKcalc::AB( double Tcrit, double Pcrit, double omg, double N,
		double &apure, double &bpure, double &da, double &d2a )
{
	double Tred, m, alph, ac, sqa, dsqa, d2sqa;

	Tred = Tk/Tcrit;
	m = 0.48 + 1.574*omg - 0.176*pow(omg,2.);
	alph = pow(1. + m*(1-sqrt(Tred)), 2.);
	ac = (0.42747)*pow(R_CONST,2.)*pow(Tcrit,2.) / Pcrit;
	apure = alph*ac;
	bpure = (0.08664)*R_CONST*Tcrit/Pcrit;
	sqa = 1. + m*(1-sqrt(Tred));
	dsqa = -0.5*m/(sqrt(Tred)*Tcrit);
	da = 2.*ac*(sqa*dsqa);
	d2sqa = 0.25*m/(pow(Tred,1.5)*pow(Tcrit,2.));
	d2a = 2.*ac*(dsqa*dsqa + sqa*d2sqa);

	return 0;
}



/// Calculates fugacities and residual functions of pure fluid species
long int TSRKcalc::FugacityPure( long int i )
{
	double Tcrit, Pcrit, Tred, asrk, bsrk, da, d2a, A, B, a2, a1, a0,
			z1, z2, z3, vol1, vol2, vol3, lnf1, lnf2, lnf3, z, vol, lnf;
	double gig, hig, sig, cpig, fugpure, grs, hrs, srs, cprs,
			cv, dPdT, dPdV, dVdT;

	// ideal gas changes from 1 bar to P at T of interest
	hig = 0.;
	sig = (-1.)*R_CONST*log(Pbar);
	gig = hig - Tk*sig;
	cpig = 0.;

	// retrieve a and b terms of cubic EoS
	Tcrit = Eosparm[i][0];
	Pcrit = Eosparm[i][1];
	Tred = Tk/Tcrit;
	asrk = Pureparm[i][0];
	bsrk = Pureparm[i][1];
	da = Pureparm[i][2];
	d2a = Pureparm[i][3];

	// solve cubic equation
	A = asrk*Pbar/(pow(R_CONST,2.)*pow(Tk,2.));
	B = bsrk*Pbar/(R_CONST*Tk);
	a2 = (-1.);
	a1 = A-B-pow(B,2.);
	a0 = -A*B;
	Cardano(a2, a1, a0, z1, z2, z3);

	// find stable roots
	vol1 = z1*R_CONST*Tk/Pbar;
	vol2 = z2*R_CONST*Tk/Pbar;
	vol3 = z3*R_CONST*Tk/Pbar;

	if (z1 > B)
		lnf1 = z1 - 1 - log(z1-B) - A/B*log(1.+B/z1);
	else
		lnf1 = 1000.;
	if (z2 > B)
		lnf2 = z2 - 1 - log(z2-B) - A/B*log(1.+B/z2);
	else
		lnf2 = 1000.;
	if (z3 > B)
		lnf3 = z3 - 1 - log(z3-B) - A/B*log(1.+B/z3);
	else
		lnf3 = 1000.;

	if (lnf2 < lnf1)
	{
		z = z2; vol = vol2; lnf = lnf2;
	}
	else
	{
		z = z1; vol = vol1; lnf = lnf1;
	}
	if (lnf3 < lnf)
	{
		z = z3; vol = vol3; lnf = lnf3;
	}
	else
	{
		z = z; vol = vol; lnf = lnf;
	}

	// calculate thermodynamic properties
	hrs = - ( 1 - z + 1./(bsrk*R_CONST*Tk) * (asrk-Tk*da) * log(1.+ bsrk/vol) )*R_CONST*Tk;
	srs = ( log(z*(1.-bsrk/vol)) + 1./(bsrk*R_CONST) * da * log(1.+bsrk/vol) )*R_CONST;
	grs = hrs - Tk*srs;

	// heat capacity part
	cv = Tk*d2a/bsrk * log(1.+B/z);
	dPdT = R_CONST/(vol-bsrk) - da/(vol*(vol+bsrk));
	dPdV = - R_CONST*Tk/pow((vol-bsrk),2.) + asrk*(2.*vol+bsrk)/pow((vol*(vol+bsrk)),2.);
	dVdT = (-1.)*(1./dPdV)*dPdT;
	cprs = cv + Tk*dPdT*dVdT - R_CONST;

	// increment thermodynamic properties
	fugpure = exp(lnf);
	Fugpure[i][0] = fugpure;
	Fugpure[i][1] = grs;
	Fugpure[i][2] = hrs;
	Fugpure[i][3] = srs;
	Fugpure[i][4] = vol;
	Fugpure[i][5] = cprs;

	return 0;
}



/// Cubic equation root solver based on Cardanos method
long int TSRKcalc::Cardano( double a2, double a1, double a0, double &z1, double &z2, double &z3 )
{
	double q, rc, q3, rc2, theta, ac, bc;

	q = (pow(a2,2.) - 3.*a1)/9.;
	rc = (2.*pow(a2,3.) - 9.*a2*a1 + 27.*a0)/54.;
	q3 = pow(q,3.);
	rc2 = pow(rc,2.);

	if (rc2 < q3)  // three real roots
	{
		theta = acos(rc/sqrt(q3));
		z1 = (-2.)*sqrt(q)*cos(theta/3.)-a2/3.;
		z2 = (-2.)*sqrt(q)*cos(theta/3.+2./3.*3.1415927)-a2/3.;
		z3 = (-2.)*sqrt(q)*cos(theta/3.-2./3.*3.1415927)-a2/3.;
	}

	else  // one real root
	{
		ac = (-1.)*rc/fabs(rc)*pow(fabs(rc)+sqrt(rc2-q3), 1./3.);
		if (ac != 0.)
			bc = q/ac;
		else
			bc = 0.;
		z1 = ac+bc-a2/3.;
		z2 = ac+bc-a2/3.;
		z3 = ac+bc-a2/3.;
	}
	return 0;
}



/// Calculates mixing properties of the fluid mixture
long int TSRKcalc::MixParam( double &amix, double &bmix )
{
	long int i, j;
	double K;
	amix = 0.;
	bmix = 0.;

	// calculate binary aij parameters
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
			K = KK[i][j];
			AA[i][j] = sqrt(Pureparm[i][0]*Pureparm[j][0])*(1.-K);
		}
	}

	// find a and b of the mixture
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
			amix = amix + x[i]*x[j]*AA[i][j];
		}
	}
	for (i=0; i<NComp; i++)
	{
		bmix = bmix + x[i]*Pureparm[i][1];
	}
	return 0;
}



/// Calculates fugacity of the bulk fluid mixture
long int TSRKcalc::FugacityMix( double amix, double bmix,
    double &fugmix, double &zmix, double &vmix )
{
	double A, B, a2, a1, a0, z1, z2, z3, vol1, vol2, vol3, lnf1, lnf2, lnf3, lnf;

	// solve cubic equation
	A = amix*Pbar/(pow(R_CONST,2.)*pow(Tk,2.));
	B = bmix*Pbar/(R_CONST*Tk);
	a2 = (-1.);
	a1 = A-B-pow(B,2.);
	a0 = -A*B;
	Cardano(a2, a1, a0, z1, z2, z3);

	// find stable roots
	vol1 = z1*R_CONST*Tk/Pbar;
	vol2 = z2*R_CONST*Tk/Pbar;
	vol3 = z3*R_CONST*Tk/Pbar;

	if (z1 > B)
		lnf1 = z1 - 1 - log(z1-B) - A/B*log(1.+B/z1);
	else
		lnf1 = 1000.;
	if (z2 > B)
		lnf2 = z2 - 1 - log(z2-B) - A/B*log(1.+B/z2);
	else
		lnf2 = 1000.;
	if (z3 > B)
		lnf3 = z3 - 1 - log(z3-B) - A/B*log(1.+B/z3);
	else
		lnf3 = 1000.;

	if (lnf2 < lnf1)
	{
		zmix = z2; vmix = vol2; lnf = lnf2;
	}
	else
	{
		zmix = z1; vmix = vol1; lnf = lnf1;
	}
	if (lnf3 < lnf)
	{
		zmix = z3; vmix = vol3; lnf = lnf3;
	}
	else
	{
		zmix = zmix; vmix = vmix; lnf = lnf;
	}

	fugmix = exp(lnf);
	PhVol = vmix;
	return 0;
}



///  Calculates fugacities and activities of fluid species in the mixture,
long int TSRKcalc::FugacitySpec( double *fugpure )
{
	long int i, j, iRet=0;
	double fugmix=0., zmix=0., vmix=0., amix=0., bmix=0., sum=0.;
	double A, B, bi, Bi, lnfci, fci;

	// Reload params to Pureparm (possibly not required any more)
	for( j=0; j<NComp; j++ )
	{
		Fugpure[j][0] = fugpure[j]/Pbar;
	}

	// calculate properties of the mixture
	iRet = MixParam( amix, bmix);
	iRet = FugacityMix( amix, bmix, fugmix, zmix, vmix);
	A = amix*Pbar/(pow(R_CONST, 2.)*pow(Tk, 2.));
	B = bmix*Pbar/(R_CONST*Tk);

	// calculate fugacity coefficient, fugacity and activity of species i
	for (i=0; i<NComp; i++)
	{
		bi = Pureparm[i][1];
		Bi = bi*Pbar/(R_CONST*Tk);

		sum = 0.;
		for (j=0; j<NComp; j++)
		{
			sum = sum + x[j]*AA[i][j];
		}

		lnfci = Bi/B*(zmix-1.) - log(zmix-B)
			+ A/B * ( Bi/B - 2./amix*sum ) * log(1.+B/zmix);
		fci = exp(lnfci);
		Fugci[i][0] = fci;  // fugacity coefficient using engineering convention
		Fugci[i][1] = x[i]*fci;  // fugacity coefficient using geology convention
		Fugci[i][2] = Fugci[i][1]/Fugpure[i][0];  // activity of species
		if (x[i]>1.0e-20)
			Fugci[i][3] = Fugci[i][2]/x[i];  // activity coefficient of species
		else
			Fugci[i][3] = 1.0;
	}

	return iRet;
}



///  calculates residual functions in the mixture
long int TSRKcalc::ResidualFunct( double *fugpure )
{
	long int i, j, iRet=0;
	double fugmix=0., zmix=0., vmix=0., amix=0., bmix=0.;
	double A, B, K, dK, d2K, Q, dQ, d2Q, damix, d2amix, ai, aj, dai, daj, d2ai, d2aj,
				cv, dPdT, dPdV, dVdT;

	// Reload params to Pureparm (possibly not required any more)
	for( j=0; j<NComp; j++ )
	{
		Fugpure[j][0] = fugpure[j]/Pbar;
	}

	// calculate properties of the mixture
	iRet = MixParam( amix, bmix);
	iRet = FugacityMix( amix, bmix, fugmix, zmix, vmix);
	A = amix*Pbar/(pow(R_CONST, 2.)*pow(Tk, 2.));
	B = bmix*Pbar/(R_CONST*Tk);

	// calculate total state functions of the mixture
	damix = 0.;
	d2amix = 0.;
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
			// pull parameters
			ai = Pureparm[i][0];
			aj = Pureparm[j][0];
			dai = Pureparm[i][2];
			daj = Pureparm[j][2];
			d2ai = Pureparm[i][3];
			d2aj = Pureparm[j][3];
			K = KK[i][j];
			dK = dKK[i][j];
			d2K = d2KK[i][j];

			// increments to derivatives
			Q = sqrt(ai*aj);
			dQ = 0.5*( sqrt(aj/ai)*dai + sqrt(ai/aj)*daj );
			d2Q = 0.5*( dai*daj/sqrt(ai*aj) + d2ai*sqrt(aj)/sqrt(ai) + d2aj*sqrt(ai)/sqrt(aj)
					- 0.5*( pow(dai,2.)*sqrt(aj)/sqrt(pow(ai,3.))
					+ pow(daj,2.)*sqrt(ai)/sqrt(pow(aj,3.)) ) );
			damix = damix + x[i]*x[j] * ( dQ*(1.-K) - Q*dK );
			d2amix = d2amix + x[i]*x[j] * ( d2Q*(1.-K) - 2.*dQ*dK - Q*d2K );
		}
	}

	// calculate thermodynamic properties
	Hrs = - ( 1. - zmix + 1./(bmix*R_CONST*Tk) * (amix-Tk*damix )
			* log(1.+bmix/vmix) )*R_CONST*Tk;
	Srs = ( log(zmix*(1.-bmix/vmix)) + 1./(bmix*R_CONST)*damix
			* log(1.+bmix/vmix) )*R_CONST;
	Grs = Hrs - Tk*Srs;

	// heat capacity part
	cv = Tk*d2amix/bmix * log(1.+B/zmix);
	dPdT = R_CONST/(vmix-bmix) - damix/(vmix*(vmix+bmix));
	dPdV = - R_CONST*Tk/pow((vmix-bmix),2.) + amix*(2.*vmix+bmix)/pow((vmix*(vmix+bmix)),2.);
	dVdT = (-1.)*(1./dPdV)*dPdT;
	CPrs = cv + Tk*dPdT*dVdT - R_CONST;
	Vrs = vmix;

	return iRet;
}



#ifndef IPMGEMPLUGIN

/// Calculates properties of pure fluids when called from DCthermo
long int TSRKcalc::SRKCalcFugPure( double Tmin, float *Cpg, double *FugProps )
{
	long int retCode = 0;
	double Coeff[7];

	for( int ii=0; ii<7; ii++ )
		Coeff[ii] = (double)Cpg[ii];

	if( (Tk >= Tmin) && (Tk < 1e4) && (Pbar >= 1e-5) && (Pbar < 1e5) )
	{
		retCode = FugacityPT( 0, Coeff );
		for( int i=0; i<6; i++ )
			FugProps[i] = Fugpure[0][i];
		return retCode;
	}

	else
	{
		for( int i=1; i<6; i++ )
			FugProps[i] = 0.;
		FugProps[0] = 1.;
		FugProps[4] = 8.31451*Tk/Pbar;
		return -1;
	}
}

#endif





//=======================================================================================================
// Peng-Robinson (PR78) model for fluid mixtures
// References: Peng and Robinson (1976), Peng and Robinson (1978)
// (c) TW July 2009
//=======================================================================================================


// Constructor
TPR78calc::TPR78calc( long int NCmp, double Pp, double Tkp ):
    TSolMod( NCmp, '7', Tkp, Pp )
{
    // aGEX = 0;
    // aVol = 0;
    Pparc = 0;
    alloc_internal();
}



TPR78calc::TPR78calc( SolutionData *sd ):
                TSolMod( sd )
{
    Pparc = aPparc;
    // aGEX = arGEX;
    // aVol = arVol;
    alloc_internal();
}



TPR78calc::~TPR78calc()
{
    free_internal();
}



///  allocate work arrays for pure fluid and fluid mixture properties
void TPR78calc::alloc_internal()
{
	Eosparm = new double [NComp][4];
	Pureparm = new double [NComp][4];
	Fugpure = new double [NComp][6];
	Fugci = new double [NComp][4];
	a = new double *[NComp];
	b = new double *[NComp];
	KK = new double *[NComp];
	dKK = new double *[NComp];
	d2KK = new double *[NComp];
	AA = new double *[NComp];

	for (long int i=0; i<NComp; i++)
	{
		a[i] = new double[NComp];
		b[i] = new double[NComp];
		KK[i] = new double[NComp];
		dKK[i] = new double[NComp];
		d2KK[i] = new double[NComp];
		AA[i] = new double[NComp];
	}
}



void TPR78calc::free_internal()
{
	long int i;

	for (i=0; i<NComp; i++)
	{
		delete[]a[i];
		delete[]b[i];
		delete[]KK[i];
		delete[]dKK[i];
		delete[]d2KK[i];
		delete[]AA[i];
	}

	delete[]Eosparm;
	delete[]Pureparm;
	delete[]Fugpure;
	delete[]Fugci;
	delete[]a;
	delete[]b;
	delete[]KK;
	delete[]dKK;
	delete[]d2KK;
	delete[]AA;

}



///  High-level method to retrieve pure fluid fugacities
long int TPR78calc::PureSpecies()
{
    long int j, retCode = 0;

    for( j=0; j<NComp; j++)
    {
        // Calling SRK EoS for pure fugacity
        retCode =  FugacityPT( j, aDCc+j*NP_DC );
        aGEX[j] = log( Fugpure[j][0] );
        Pparc[j] = Fugpure[j][0]*Pbar;  // fure fluid fugacity (required for performance)
        aVol[j] = Fugpure[j][4]*10.;  // molar volume of pure fluid component, J/bar to cm3
    } // j

    if ( retCode )
    {
        char buf[150];
        sprintf(buf, "PR78 fluid: calculation of pure fugacity failed");
					Error( "E71IPM IPMgamma: ",  buf );
    }

    return 0;
}


///    High-level method to calculate T,P corrected binary interaction parameters
long int TPR78calc::PTparam()
{
	long int j, i;

	PureSpecies();

	// set all interaction parameters zero
	for( j=0; j<NComp; j++ )
	{
		for( i=0; i<NComp; i++ )
		{
			KK[j][i] = 0.;
			dKK[j][i] = 0.;
			d2KK[j][i] = 0.;
		}
	}

	switch ( MixCode )
	{
        case MR_UNDEF_:
        case MR_WAAL_:
			MixingWaals();
			break;
		case MR_CONST_:
			MixingConst();
			break;
		case MR_TEMP_:
			MixingTemp();
			break;
		default:
			break;
	}

	return 0;
}



/// High-level method to retrieve activity coefficients of the fluid mixture
long int TPR78calc::MixMod()
{
	long int j, iRet;

	iRet = FugacitySpec( Pparc );
	phVOL[0] = PhVol * 10.;

    for(j=0; j<NComp; j++)
    {
        if( Fugci[j][3] > 1e-23 )
        	lnGamma[j] = log( Fugci[j][3] );
        else
        	lnGamma[j] = 0;
    }
    if ( iRet )
    {
    	char buf[150];
    	sprintf(buf, "PR78 fluid: calculation failed");
			Error( "E71IPM IPMgamma: ",  buf );
    }
    return iRet;
}



/// High-level method to retrieve residual functions of the fluid mixture
long int TPR78calc::ExcessProp( double *Zex )
{
	long int iRet;

	iRet = ResidualFunct( Pparc );

    if ( iRet )
    {
    	char buf[150];
    	sprintf(buf, "PR78 fluid: calculation failed");
    	Error( "E71IPM IPMgamma: ",  buf );
    }

	Ars = Grs - Vrs*Pbar;
	Urs = Hrs - Vrs*Pbar;

	// assignments (residual functions)
	Zex[0] = Grs;
	Zex[1] = Hrs;
	Zex[2] = Srs;
	Zex[3] = CPrs;
	Zex[4] = Vrs;
	Zex[5] = Ars;
	Zex[6] = Urs;

	return iRet;

}



/// basic van der waals mixing rule
long int TPR78calc::MixingWaals()
{
	// currently no calculations

	return 0;
}



/// constant one-term interaction parameter
long int TPR78calc::MixingConst()
{
	long int ip, i1, i2;
	double k, dk, d2k;

	if( NPcoef > 0 )
	{
		// transfer interaction parameters that have non-standard value
		for ( ip=0; ip<NPar; ip++ )
		{
			i1 = aIPx[MaxOrd*ip];
			i2 = aIPx[MaxOrd*ip+1];
			k = aIPc[NPcoef*ip];
			dk = 0.;
			d2k = 0.;
			KK[i1][i2] = k;
			dKK[i1][i2] = dk;
			d2KK[i1][i2] = d2k;
			KK[i2][i1] = k;   // symmetric case
			dKK[i2][i1] = dk;
			d2KK[i2][i1] = d2k;
		}
	}

	return 0;
}



/// temperature dependent one-term interaction parameter
long int TPR78calc::MixingTemp()
{
	long int i, j, ip, i1, i2;
	double ai, aj, bi, bj, di, dj, dai, daj, d2ai, d2aj, ddi, ddj, d2di, d2dj,
				U, V, dU, dV, d2U, d2V, tmp, k, dk, d2k, C;

	// set model specific interaction parameters zero
	for( j=0; j<NComp; j++ )
	{
		for( i=0; i<NComp; i++ )
		{
			a[j][i] = 0.;
			b[j][i] = 0.;
		}
	}

	if( NPcoef > 0 )
	{
		// transfer parameters that have non-standard value
		for ( ip=0; ip<NPar; ip++ )
		{
			i1 = aIPx[MaxOrd*ip];
			i2 = aIPx[MaxOrd*ip+1];
			a[i1][i2] = aIPc[NPcoef*ip];
			b[i1][i2] = aIPc[NPcoef*ip+1];
			a[i2][i1] = aIPc[NPcoef*ip];  // symmetric case
			b[i2][i1] = aIPc[NPcoef*ip+1];
		}
	}

	// calculate binary interaction parameters
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
			if ( a[i][j] == 0.0 )
				tmp = 1.0;
			else
				tmp = a[i][j];

			// read a, b, da, d2a
			ai = Pureparm[i][0];
			aj = Pureparm[j][0];
			bi = Pureparm[i][1];
			bj = Pureparm[j][1];
			dai = Pureparm[i][2];
			daj = Pureparm[j][2];
			d2ai = Pureparm[i][3];
			d2aj = Pureparm[j][3];

			// calculate k and derivatives
			di = sqrt(ai)/bi;
			dj = sqrt(aj)/bj;
			ddi = (0.5/bi) * pow(ai,-0.5) * dai;
			ddj = (0.5/bj) * pow(aj,-0.5) * daj;
			d2di = (0.5/bi) * ( (-0.5)*pow(ai,-1.5)*dai*dai + pow(ai,-0.5)*d2ai );
			d2dj = (0.5/bj) * ( (-0.5)*pow(aj,-1.5)*daj*daj + pow(aj,-0.5)*d2aj );

			C = ( b[i][j]/tmp - 1. );
			U = a[i][j]*pow((298.15/Tk),C) - pow((di-dj),2.);
			V = (2.*di*dj);
			dU = - ( a[i][j]*C*pow((298.15/Tk),C) ) / Tk - 2.*(di-dj)*(ddi-ddj);
			dV = 2.*( ddi*dj + di*ddj );
			d2U = ( a[i][j]*pow(C,2.)*pow((298.15/Tk),C) ) / pow(Tk,2.)
						+ ( a[i][j]*C*pow((298.15/Tk),C) ) / pow(Tk,2.)
						- 2.*( pow((ddi-ddj),2.) + (di-dj)*(d2di-d2dj) );
			d2V = 2.*( d2di*dj + 2.*ddi*ddj + di*d2dj );
			k = U/V;
			dk = (dU*V-U*dV)/pow(V,2.);
			d2k = (d2U*V+dU*dV)*pow(V,2.)/pow(V,4.) - (dU*V)*(2.*V*dV)/pow(V,4.)
						- (dU*dV+U*d2V)*pow(V,2.)/pow(V,4.) + (U*dV)*(2.*V*dV)/pow(V,4.);

			// assignments
			KK[i][j] = k;
			dKK[i][j] = dk;
			d2KK[i][j] = d2k;
		}
	}

	return 0;
}



/// calculates ideal mixing properties
long int TPR78calc::IdealProp( double *Zid )
{
	long int j;
	double s, sc, sp;

	s = 0.0;
	for (j=0; j<NComp; j++)
	{
		if ( x[j] > 1.0e-32 )
			s += x[j]*log(x[j]);
	}
	sc = (-1.)*R_CONST*s;
	sp = (-1.)*R_CONST*log(Pbar);
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
	Sid = sc + sp;
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



/// High-level method to retrieve pure fluid properties
long int TPR78calc::FugacityPT( long int i, double *EoSparam )
{
	long int iRet = 0;
	double Tcrit, Pcrit, omg, N, apure, bpure, da, d2a;

	// reads EoS parameters from database into work array
	if( !EoSparam )
		return -1;  // Memory alloc error

	Eosparm[i][0] = EoSparam[0];  // critical temperature in K
	Eosparm[i][1] = EoSparam[1];  // critical pressure in bar
	Eosparm[i][2] = EoSparam[2];  // Pitzer acentric factor omega
	Eosparm[i][3] = EoSparam[3];  // empirical EoS parameter N
	Tcrit = Eosparm[i][0];
	Pcrit = Eosparm[i][1];
	omg = Eosparm[i][2];
	N = Eosparm[i][3];

	AB(Tcrit, Pcrit, omg, N, apure, bpure, da, d2a);

	Pureparm[i][0] = apure;  // a parameter
	Pureparm[i][1] = bpure;  // b parameter
	Pureparm[i][2] = da;  // da/dT
	Pureparm[i][3] = d2a;  // d2a/dT2

	iRet = FugacityPure( i );
    if( iRet)
    	return iRet;

    return iRet;
}



/// Calculates attractive (a) and repulsive (b) parameter of SRK equation of state
/// and partial derivatives of alpha function
long int TPR78calc::AB( double Tcrit, double Pcrit, double omg, double N,
		double &apure, double &bpure, double &da, double &d2a )
{
	double Tred, k, alph, ac, sqa, dsqa, d2sqa;

	Tred = Tk/Tcrit;
	if (omg <= 0.491)
			k = 0.37464 + 1.54226*omg - 0.26992*pow(omg,2.);
	else
			k = 0.379642 + 1.48503*omg - 0.164423*pow(omg,2.) + 0.0166666*pow(omg,3.);

	alph = pow(1.+k*(1.-sqrt(Tred)),2.);
	ac = (0.457235529)*pow(R_CONST,2.)*pow(Tcrit,2.) / Pcrit;
	apure = alph*ac;
	bpure = (0.0777960739)*R_CONST*Tcrit/Pcrit;
	sqa = 1.+k*(1.-sqrt(Tred));
	dsqa = -0.5*k/(sqrt(Tred)*Tcrit);
	da = 2.*ac*(sqa*dsqa);
	d2sqa = 0.25*k/(pow(Tred,1.5)*pow(Tcrit,2.));
	d2a = 2.*ac*(dsqa*dsqa + sqa*d2sqa);

	return 0;
}



/// Calculates fugacities and residual functions of pure fluid species
long int TPR78calc::FugacityPure( long int i )
{
	double Tcrit, Pcrit, Tred, apr, bpr, alph, da, d2a, k, A, B, a2, a1, a0,
			z1, z2, z3, vol1, vol2, vol3, lnf1, lnf2, lnf3, z, vol, lnf;
	double gig, hig, sig, cpig, fugpure, grs, hrs, srs, cprs,
			cv, dPdT, dPdV, dVdT;

	// ideal gas changes from 1 bar to P at T of interest
	hig = 0.;
	sig = (-1.)*R_CONST*log(Pbar);
	gig = hig - Tk*sig;
	cpig = 0.;

	// retrieve a and b terms of cubic EoS
	Tcrit = Eosparm[i][0];
	Pcrit = Eosparm[i][1];
	Tred = Tk/Tcrit;
	apr = Pureparm[i][0];
	bpr = Pureparm[i][1];
	da = Pureparm[i][2];
	d2a = Pureparm[i][3];

	// solve cubic equation
	A = apr*Pbar/(pow(R_CONST,2.)*pow(Tk,2.));
	B = bpr*Pbar/(R_CONST*Tk);
	a2 = B - 1.;
	a1 = A - 3.*pow(B,2.) - 2.*B;
	a0 = pow(B,3.) + pow(B,2.) - A*B;
	Cardano(a2, a1, a0, z1, z2, z3);

	// find stable roots
	vol1 = z1*R_CONST*Tk/Pbar;
	vol2 = z2*R_CONST*Tk/Pbar;
	vol3 = z3*R_CONST*Tk/Pbar;
	if (z1 > B)
		lnf1 = (-1.)*log(z1-B)
			- A/(B*sqrt(8.))*log((z1+(1.+sqrt(2.))*B)/(z1+(1.-sqrt(2.))*B))+z1-1.;
	else
		lnf1 = 1000.;
	if (z2 > B)
		lnf2 = (-1.)*log(z2-B)
			- A/(B*sqrt(8.))*log((z2+(1.+sqrt(2.))*B)/(z2+(1.-sqrt(2.))*B))+z2-1.;
	else
		lnf2 = 1000.;
	if (z3 > B)
		lnf3 = (-1.)*log(z3-B)
			- A/(B*sqrt(8.))*log((z3+(1.+sqrt(2.))*B)/(z3+(1.-sqrt(2.))*B))+z3-1.;
	else
		lnf3 = 1000.;

	if (lnf2 < lnf1)
	{
		z = z2; vol = vol2; lnf = lnf2;
	}
	else
	{
		z = z1; vol = vol1; lnf = lnf1;
	}
	if (lnf3 < lnf)
	{
		z = z3; vol = vol3; lnf = lnf3;
	}
	else
	{
		z = z; vol = vol; lnf = lnf;
	}

	// calculate thermodynamic properties
	alph = apr/((0.457235529)*pow(R_CONST,2.)*pow(Tcrit,2.) / Pcrit);
	k = (sqrt(alph)-1.)/(1.-sqrt(Tred));
	grs = R_CONST*Tk*(z-1.-log(z-B)-A/(B*sqrt(8.))
				*log((z+(1+sqrt(2.))*B)/(z+(1-sqrt(2.))*B)));
	hrs = R_CONST*Tk*(z-1.-log((z+(1+sqrt(2.))*B)/(z+(1-sqrt(2.))*B))
				*A/(B*sqrt(8.))*(1+k*sqrt(Tred)/sqrt(alph)));
	srs = (hrs-grs)/Tk;

	// heat capacity part
	cv = Tk*d2a/(bpr*sqrt(8.))
			 * log( (z+B*(1.+sqrt(2.)))/(z+B*(1.-sqrt(2.))) );
	dPdT = R_CONST/(vol-bpr) - da/( vol*(vol+bpr) + bpr*(vol-bpr) );
	dPdV = - R_CONST*Tk/pow((vol-bpr),2.) + 2.*apr*(vol+bpr)/pow((vol*(vol+bpr)+bpr*(vol-bpr)),2.);
	dVdT = (-1.)*(1./dPdV)*dPdT;
	cprs = cv + Tk*dPdT*dVdT - R_CONST;

	// increment thermodynamic properties
	fugpure = exp(lnf);
	Fugpure[i][0] = fugpure;
	Fugpure[i][1] = grs;
	Fugpure[i][2] = hrs;
	Fugpure[i][3] = srs;
        Fugpure[i][4] = vol;
        Fugpure[i][5] = cprs;

	return 0;
}



/// Cubic equation root solver based on Cardanos method
long int TPR78calc::Cardano( double a2, double a1, double a0, double &z1, double &z2, double &z3 )
{
	double q, rc, q3, rc2, theta, ac, bc;

	q = (pow(a2,2.) - 3.*a1)/9.;
	rc = (2.*pow(a2,3.) - 9.*a2*a1 + 27.*a0)/54.;
	q3 = pow(q,3.);
	rc2 = pow(rc,2.);
	if (rc2 < q3)  // three real roots
	{
		theta = acos(rc/sqrt(q3));
		z1 = (-2.)*sqrt(q)*cos(theta/3.)-a2/3.;
		z2 = (-2.)*sqrt(q)*cos(theta/3.+2./3.*3.1415927)-a2/3.;
		z3 = (-2.)*sqrt(q)*cos(theta/3.-2./3.*3.1415927)-a2/3.;
	}
	else  // one real root
	{
		ac = (-1.)*rc/fabs(rc)*pow(fabs(rc)+sqrt(rc2-q3), 1./3.);
		if (ac != 0.)
			bc = q/ac;
		else
			bc = 0.;
		z1 = ac+bc-a2/3.;
		z2 = ac+bc-a2/3.;
		z3 = ac+bc-a2/3.;
	}
	return 0;
}



/// Calculates mixing properties of the fluid mixture
long int TPR78calc::MixParam( double &amix, double &bmix )
{
	long int i, j;
	double K;
	amix = 0.;
	bmix = 0.;

	// calculate binary aij parameters
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
            K = KK[i][j];
			AA[i][j] = sqrt(Pureparm[i][0]*Pureparm[j][0])*(1.-K);
		}
	}

	// find a and b of the mixture
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
			amix = amix + x[i]*x[j]*AA[i][j];
		}
	}
	for (i=0; i<NComp; i++)
	{
		bmix = bmix + x[i]*Pureparm[i][1];
	}

	return 0;
}



/// Calculates fugacity of the bulk fluid mixture
long int TPR78calc::FugacityMix( double amix, double bmix,
    double &fugmix, double &zmix, double &vmix )
{
	double A, B, a2, a1, a0, z1, z2, z3, vol1, vol2, vol3, lnf1, lnf2, lnf3, lnf;

	// solve cubic equation
	A = amix*Pbar/(pow(R_CONST,2.)*pow(Tk,2.));
	B = bmix*Pbar/(R_CONST*Tk);
	a2 = B - 1.;
	a1 = A - 3.*pow(B,2.) - 2.*B;
	a0 = pow(B,3.) + pow(B,2.) - A*B;
	Cardano( a2, a1, a0, z1, z2, z3 );

	// find stable roots
	vol1 = z1*R_CONST*Tk/Pbar;
	vol2 = z2*R_CONST*Tk/Pbar;
	vol3 = z3*R_CONST*Tk/Pbar;
	if (z1 > B)
		lnf1 = (-1.)*log(z1-B)
			- A/(B*sqrt(8.))*log((z1+(1.+sqrt(2.))*B)/(z1+(1.-sqrt(2.))*B))+z1-1.;
	else
		lnf1 = 1000.;
	if (z2 > B)
		lnf2 = (-1.)*log(z2-B)
			- A/(B*sqrt(8.))*log((z2+(1.+sqrt(2.))*B)/(z2+(1.-sqrt(2.))*B))+z2-1.;
	else
		lnf2 = 1000.;
	if (z3 > B)
		lnf3 = (-1.)*log(z3-B)
			- A/(B*sqrt(8.))*log((z3+(1.+sqrt(2.))*B)/(z3+(1.-sqrt(2.))*B))+z3-1.;
	else
		lnf3 = 1000.;

	if (lnf2 < lnf1)
	{
		zmix = z2; vmix = vol2; lnf = lnf2;
	}
	else
	{
		zmix = z1; vmix = vol1; lnf = lnf1;
	}
	if (lnf3 < lnf)
	{
		zmix = z3; vmix = vol3; lnf = lnf3;
	}
	else
	{
		zmix = zmix; vmix = vmix; lnf = lnf;
	}
	fugmix = exp(lnf);
    PhVol = vmix;

	return 0;
}



/// Calculates fugacities and activities of fluid species in the mixture,
long int TPR78calc::FugacitySpec( double *fugpure )
{
    long int i, j, iRet=0;
	double fugmix=0., zmix=0., vmix=0., amix=0., bmix=0., sum=0.;
	double A, B, lnfci, fci;

    // Reload params to Pureparm
    for( j=0; j<NComp; j++ )
    {
      Fugpure[j][0] = fugpure[j]/Pbar;
    }

	// retrieve properties of the mixture
	iRet = MixParam( amix, bmix );
	iRet = FugacityMix( amix, bmix, fugmix, zmix, vmix );
	A = amix*Pbar/(pow(R_CONST, 2.)*pow(Tk, 2.));
	B = bmix*Pbar/(R_CONST*Tk);

	// calculate fugacity coefficient, fugacity and activity of species i
	for (i=0; i<NComp; i++)
	{
		sum = 0.;
		for (j=0; j<NComp; j++)
		{
			sum = sum + x[j]*AA[i][j];
		}
		lnfci = Pureparm[i][1]/bmix*(zmix-1.) - log(zmix-B)
		      + A/(sqrt(8.)*B)*(2.*sum/amix-Pureparm[i][1]/bmix)
                      * log((zmix+B*(1.-sqrt(2.)))/(zmix+B*(1.+sqrt(2.))));
		fci = exp(lnfci);
		Fugci[i][0] = fci;  // fugacity coefficient using engineering convention
		Fugci[i][1] = x[i]*fci;  // fugacity coefficient using geology convention
		Fugci[i][2] = Fugci[i][1]/Fugpure[i][0];  // activity of species
		if (x[i]>1.0e-20)
			Fugci[i][3] = Fugci[i][2]/x[i];  // activity coefficient of species
		else
			Fugci[i][3] = 1.0;
	}

	return iRet;
}



/// calculates residual functions in the mixture
long int TPR78calc::ResidualFunct( double *fugpure )
{
    long int i, j, iRet=0;
	double fugmix=0., zmix=0., vmix=0., amix=0., bmix=0.;
	double A, B, K, dK, d2K, Q, dQ, d2Q, damix, d2amix, ai, aj, dai, daj, d2ai, d2aj,
			cv, dPdT, dPdV, dVdT;

    // Reload params to Pureparm (probably now obsolete?)
    for( j=0; j<NComp; j++ )
    {
      Fugpure[j][0] = fugpure[j]/Pbar;
    }

	// retrieve properties of the mixture
	iRet = MixParam( amix, bmix );
	iRet = FugacityMix( amix, bmix, fugmix, zmix, vmix );
	A = amix*Pbar/(pow(R_CONST, 2.)*pow(Tk, 2.));
	B = bmix*Pbar/(R_CONST*Tk);

	// calculate total state functions of the mixture
	damix = 0.;
	d2amix = 0.;
	for (i=0; i<NComp; i++)
	{
		for (j=0; j<NComp; j++)
		{
			// pull parameters
			ai = Pureparm[i][0];
			aj = Pureparm[j][0];
			dai = Pureparm[i][2];
			daj = Pureparm[j][2];
			d2ai = Pureparm[i][3];
			d2aj = Pureparm[j][3];
			K = KK[i][j];
			dK = dKK[i][j];
			d2K = d2KK[i][j];

			// increments to derivatives
			Q = sqrt(ai*aj);
			dQ = 0.5*( sqrt(aj/ai)*dai + sqrt(ai/aj)*daj );
			d2Q = 0.5*( dai*daj/sqrt(ai*aj) + d2ai*sqrt(aj)/sqrt(ai) + d2aj*sqrt(ai)/sqrt(aj)
					- 0.5*( pow(dai,2.)*sqrt(aj)/sqrt(pow(ai,3.))
					+ pow(daj,2.)*sqrt(ai)/sqrt(pow(aj,3.)) ) );
			damix = damix + x[i]*x[j] * ( dQ*(1.-K) - Q*dK );
			d2amix = d2amix + x[i]*x[j] * ( d2Q*(1.-K) - 2.*dQ*dK - Q*d2K );
		}
	}

	// calculate thermodynamic properties
	Grs = (amix/(R_CONST*Tk*sqrt(8.)*bmix) * log((vmix+(1.-sqrt(2.))*bmix)
		/ (vmix+(1.+sqrt(2.))*bmix))-log(zmix*(1.-bmix/vmix))+zmix-1.)*R_CONST*Tk;
	Hrs = ((amix-Tk*damix)/(R_CONST*Tk*sqrt(8.)*bmix)*log((vmix+(1.-sqrt(2.))
		*bmix)/(vmix+(1.+sqrt(2.))*bmix))+zmix-1.)*R_CONST*Tk;
	Srs = (Hrs - Grs)/Tk;

	// heat capacity part
	cv = Tk*d2amix/(bmix*sqrt(8.))
			 * log( (zmix+B*(1.+sqrt(2.)))/(zmix+B*(1.-sqrt(2.))) );
	dPdT = R_CONST/(vmix-bmix) - damix/( vmix*(vmix+bmix) + bmix*(vmix-bmix) );
	dPdV = - R_CONST*Tk/pow((vmix-bmix),2.) + 2.*amix*(vmix+bmix)/pow((vmix*(vmix+bmix)+bmix*(vmix-bmix)),2.);
	dVdT = (-1.)*(1./dPdV)*dPdT;
	CPrs = cv + Tk*dPdT*dVdT - R_CONST;
	Vrs = vmix;

	return iRet;
}



#ifndef IPMGEMPLUGIN

/// Calculates properties of pure fluids when called from DCthermo
long int TPR78calc::PR78CalcFugPure( double Tmin, float *Cpg, double *FugProps )
{
	long int retCode = 0;
	double Coeff[7];

	for( int ii=0; ii<7; ii++ )
		Coeff[ii] = (double)Cpg[ii];

	if( (Tk >= Tmin) && (Tk < 1e4) && (Pbar >= 1e-5) && (Pbar < 1e5) )
	{
		retCode = FugacityPT( 0, Coeff );
		for( int i=0; i<6; i++ )
			FugProps[i] = Fugpure[0][i];
		return retCode;
	}

	else
	{
		for( int i=1; i<6; i++ )
			FugProps[i] = 0.;
		FugProps[0] = 1.;
		FugProps[4] = 8.31451*Tk/Pbar;
		return -1;
	}
}

#endif





//=======================================================================================================
// Compensated Redlich-Kwong (CORK) model for fluid mixtures
// References: Holland and Powell (1991)
// (c) TW May 2010
//=======================================================================================================


// Constructor
TCORKcalc::TCORKcalc( long int NCmp, double Pp, double Tkp, char Eos_Code ):
    TSolMod( NCmp, '8', Tkp, Pp )
{
    RR = 0.00831451;    // gas constant in kbar
    Pkb = Pbar/1000.;   // pressure in kbar
    // aGEX = 0;
    // aVol = 0;
    Pparc = 0;
    alloc_internal();
    EosCode[0] = Eos_Code;     // must be changed in m_dcomp.cpp line 818

}



TCORKcalc::TCORKcalc( SolutionData *sd ):
                TSolMod( sd )
{
    RR = 0.00831451;    // gas constant in kbar
    Pkb = Pbar/1000.;   // pressure in kbar
    Pparc = aPparc;
    // aGEX = arGEX;
    // aVol = arVol;
    alloc_internal();

    for( long int j=0; j<NComp; j++ )
        // EosCode[j] = sd->DC_Codes[j]; // may be change to EosCode[j] = sd->TP_Code[j][3]; SD
        EosCode[j] = sd->TP_Code[j][3];
}



TCORKcalc::~TCORKcalc()
{
    free_internal();
}



/// allocate work arrays for pure fluid and fluid mixture properties
void TCORKcalc::alloc_internal()
{
    EosCode = new char[NComp];
    phi = new double [NComp];
    dphi = new double [NComp];
    d2phi = new double [NComp];
    dphip = new double [NComp];
    Eosparm = new double [NComp][2];
    Fugpure = new double [NComp][6];
    Fugci = new double [NComp][4];
    Rho = new double [NComp][11];
    A = new double *[NComp];
    W = new double *[NComp];
    B = new double *[NComp];
    dB = new double *[NComp];
    d2B = new double *[NComp];
    dBp = new double *[NComp];

    for (long int i=0; i<NComp; i++)
    {
        A[i] = new double [NComp];
        W[i] = new double [NComp];
        B[i] = new double [NComp];
        dB[i] = new double [NComp];
        d2B[i] = new double [NComp];
        dBp[i] = new double [NComp];
    }
}



void TCORKcalc::free_internal()
{
    long int i;

    for (i=0; i<NComp; i++)
    {
        delete[]A[i];
        delete[]W[i];
        delete[]B[i];
        delete[]dB[i];
        delete[]d2B[i];
        delete[]dBp[i];
    }

    delete[]EosCode;
    delete[]phi;
    delete[]dphi;
    delete[]d2phi;
    delete[]dphip;
    delete[]Eosparm;
    delete[]Fugpure;
    delete[]Fugci;
    delete[]Rho;
    delete[]A;
    delete[]W;
    delete[]B;
    delete[]dB;
    delete[]d2B;
    delete[]dBp;
}



/// high-level method to retrieve pure fluid fugacities
long int TCORKcalc::PureSpecies()
{
    long int j, retCode = 0;

    for( j=0; j<NComp; j++ )
    {
        // Calling CORK EoS for pure fugacity
        retCode =  FugacityPT( j, aDCc+j*NP_DC );
        aGEX[j] = log( Fugpure[j][0] );
        Pparc[j] = Fugpure[j][0]*Pbar;  // fure fluid fugacity (required for performance)
        aVol[j] = Fugpure[j][4]*10.;  // molar volume of pure fluid component, J/bar to cm3
    } // j

    if ( retCode )
    {
            char buf[150];
            sprintf(buf, "CORK fluid: calculation of pure fluid fugacity failed");
                    Error( "E71IPM IPMgamma: ",  buf );
    }

    return 0;
}



/// high-level method to calculate T,P corrected binary interaction parameters
long int TCORKcalc::PTparam()
{
    long int j, i, ip, i1, i2;
    double a;

    Pkb = Pbar/1000.;

    PureSpecies();

    // set all interaction parameters zero
    for( j=0; j<NComp; j++ )
    {
        for( i=0; i<NComp; i++ )
        {
            A[j][i] = 0.;
            W[j][i] = 0.;
            B[j][i] = 0.;
            dB[j][i] = 0.;
            d2B[j][i] = 0.;
            dBp[j][i] = 0.;
        }
    }

    // transfer interaction parameters that have non-standard value
    if( NPcoef > 0 )
    {
        for ( ip=0; ip<NPar; ip++ )
        {
            i1 = aIPx[MaxOrd*ip];
            i2 = aIPx[MaxOrd*ip+1];
            a = aIPc[NPcoef*ip];
            A[i1][i2] = a;
            A[i2][i1] = a;  // symmetric case
        }
    }

    return 0;
}



/// high-level method to retrieve activity coefficients of the fluid mixture
long int TCORKcalc::MixMod()
{
    long int i, j, k;
    double dj, dk, sumphi, lnGam, Gam, vi, vj, vk;

    // calculate phi values
    sumphi = 0.;
    for (i=0; i<NComp; i++)
    {
        vi = Fugpure[i][4];
        sumphi = sumphi + x[i]*vi;
    }
    for (i=0; i<NComp; i++)
    {
        vi = Fugpure[i][4];
        phi[i] = x[i]*vi/sumphi;
    }

    // interaction parameters
    for (i=0; i<NComp; i++)
    {
        for (j=i+1; j<NComp; j++)
        {
            vi = Fugpure[i][4];
            vj = Fugpure[j][4];
            W[i][j] = A[i][j]*(vi+vj)/(vi*vj);
        }
    }

    // activity coefficients
    for (i=0; i<NComp; i++)
    {
        lnGam = 0.;
        for (j=0; j<NComp; j++)
        {
            for (k=j+1; k<NComp; k++)
            {
                vi = Fugpure[i][4];
                vj = Fugpure[j][4];
                vk = Fugpure[k][4];
                if (i==j)
                    dj = 1.;
                else
                    dj = 0.;
                if (i==k)
                    dk = 1.;
                else
                    dk = 0.;
                lnGam = lnGam - (dj-phi[j])*(dk-phi[k])*W[j][k]*2.*vi/(vj+vk);
            }
        }
        Gam = exp(lnGam/(R_CONST*Tk));
        Fugci[i][0] = Gam;
        if( Fugci[i][0] > 1e-23 )
            lnGamma[i] = log( Fugci[i][0] );
        else
            lnGamma[i] = 0.;

    }  // i

    return 0;
}



/// high-level method to retrieve residual functions of the fluid mixture
long int TCORKcalc::ExcessProp( double *Zex )
{
    long int iRet;

    iRet = ResidualFunct();

    if ( iRet )
    {
        char buf[150];
        sprintf(buf, "CORK fluid: calculation failed");
        Error( "E71IPM IPMgamma: ",  buf );
    }

    Ars = Grs - Vrs*Pbar;
    Urs = Hrs - Vrs*Pbar;

    // assignments (residual functions)
    Zex[0] = Grs;
    Zex[1] = Hrs;
    Zex[2] = Srs;
    Zex[3] = CPrs;
    Zex[4] = Vrs;
    Zex[5] = Ars;
    Zex[6] = Urs;

    return iRet;
}



/// calculates ideal mixing properties
long int TCORKcalc::IdealProp( double *Zid )
{
    long int j;
    double s, sc, sp;

    s = 0.0;
    for ( j=0; j<NComp; j++ )
    {
        if ( x[j] > 1.0e-32 )
            s += x[j]*log(x[j]);
    }
    sc = (-1.)*R_CONST*s;
    sp = (-1.)*R_CONST*log(Pbar);
    Hid = 0.0;
    CPid = 0.0;
    Vid = 0.0;
    Sid = sc + sp;
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



/// high-level method to retrieve pure fluid properties
long int TCORKcalc::FugacityPT( long int j, double *EoSparam )
{
    long int iErr = 0;

    // reads EoS parameters from database into work array
    if( !EoSparam )
        return -1;  // Memory alloc error
    Eosparm[j][0] = EoSparam[0];  // critical temperature in K
    Eosparm[j][1] = EoSparam[1];  // critical pressure in bar

    // select subroutines for different fluid types
    switch ( EosCode[j] )
    {
            // case DC_GAS_H2O_:  // H2O
            case CEM_H2O_:  // H2O
                    iErr = FugacityH2O( j );
                    break;
            // case DC_GAS_CO2_:  // CO2
            case CEM_CO2_: // CO2
                    iErr = FugacityCO2( j );
                    break;
            // case DC_GAS_COMP_:  // other fluids
            // case DC_GAS_H2_:
            // case DC_GAS_N2_:
            case CEM_GAS_:  // other fluids
            case CEM_CH4_:
            case CEM_N2_:
            case CEM_H2_:
            case CEM_O2_:
            case CEM_AR_:
            case CEM_PO_:
            case CEM_NP_:
                    iErr = FugacityCorresponding( j );
                    break;
            default:
                    iErr = 3;
                    break;
    }

    return iErr;
}



/// calculates fugacity and state functions of H2O
long int TCORKcalc::FugacityH2O( long int j )
{
    long int phState;   // 1: vapor, 2: liquid
    double a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a, da, d2a, b, c0, c1, c,
            d0, d1, d, e0, e, dc, dd, p0, Psat, vc, fc, v, g1, g2, g3, rho, z;
    double grs, hrs, srs, cprs, cvrs, dpr, dprr, dpt, dptt, dprt,
            drt, drtt, drp, drpp, drtp, dpv, dvt;

    a0 = 1113.4;
    a1 = -0.88517;
    a2 = 4.53e-3;
    a3 = -1.3183e-5;
    a4 = -0.22291;
    a5 = -3.8022e-4;
    a6 = 1.7791e-7;
    a7 = 5.8487;
    a8 = -2.1370e-2;
    a9 = 6.8133e-5;
    b = 1.465;
    c0 = -8.909e-2;
    c1 = 0.;
    d0 = 1.9853e-3;
    d1 = 0.;
    e0 = 8.0331e-2;
    c = c0 + c1*Tk;
    d = d0 + d1*Tk;
    e = e0;
    dc = c1;
    dd = d1;
    p0 = 2.;

    // molar volume and fugacity coefficient
    if (Tk > 695.)   // supercritical fluid phase
    {
        phState = 1;
        a = a0 + a4*(Tk-673.) + a5*pow((Tk-673.),2.) + a6*pow((Tk-673.),3.);
        da = -a4 - 2.*a5*(673.-Tk) - 3.*a6*pow((673.-Tk),2.);
        d2a = 2.*a5 + 6.*a6*(673.-Tk);
        VolumeFugacity(phState, Pkb, p0, a, b, c, d, e, vc, fc);
    }

    else   // subcritical region
    {
        Psat = (-13.627e-3) + (7.29395e-7)*pow(Tk,2.) - (2.34622e-9)*pow(Tk,3.) + (4.83607e-15)*pow(Tk,5.);
        if (Pkb < Psat)
        {
            phState = 1;   // vapor
            if (Tk < 673.)
            {
                a = a0 + a7*(673.-Tk) + a8*pow((673.-Tk),2.) + a9*pow((673.-Tk),3.);
                da = -a7 - 2.*a8*(673.-Tk) - 3.*a9*pow((673.-Tk),2.);
                d2a = 2.*a8 + 6.*a9*(673.-Tk);
            }
            else
            {
                a = a0 + a4*(673.-Tk) + a5*pow((673.-Tk),2.) + a6*pow((673.-Tk),3.);
                da = -a4 - 2.*a5*(673.-Tk) - 3.*a6*pow((673.-Tk),2.);
                d2a = 2.*a5 + 6.*a6*(673.-Tk);
            }
            VolumeFugacity(phState, Pkb, p0, a, b, c, d, e, vc, fc);
        }

        else
        {
            phState = 1;   // gaseous phase at Psat
            if (Tk < 673.)
            {
                a = a0 + a7*(673.-Tk) + a8*pow((673.-Tk),2.) + a9*pow((673.-Tk),3.);
                da = -a7 - 2.*a8*(673.-Tk) - 3.*a9*pow((673.-Tk),2.);
                d2a = 2.*a8 + 6.*a9*(673.-Tk);
            }
            else
            {
                a = a0 + a4*(673.-Tk) + a5*pow((673.-Tk),2.) + a6*pow((673.-Tk),3.);
                da = -a4 - 2.*a5*(673.-Tk) - 3.*a6*pow((673.-Tk),2.);
                d2a = 2.*a5 + 6.*a6*(673.-Tk);
            }
            VolumeFugacity(phState, Psat, p0, a, b, c, d, e, v, g1);

            phState = 2; // liquid phase at Psat
            if (Tk < 673.)
            {
                a = a0 + a1*(673.-Tk) + a2*pow((673.-Tk),2.) + a3*pow((673.-Tk),3.);
                da = -a1 - 2.*a2*(673.-Tk) - 3.*a3*pow((673.-Tk),2.);
                d2a = 2.*a2 + 6.*a3*(673.-Tk);
            }
            else
            {
                a = a0 + a4*(673.-Tk) + a5*pow((673.-Tk),2.) + a6*pow((673.-Tk),3.);
                da = -a4 - 2.*a5*(673.-Tk) - 3.*a6*pow((673.-Tk),2.);
                d2a = 2.*a5 + 6.*a6*(673.-Tk);
            }
            VolumeFugacity(phState, Psat, p0, a, b, c, d, e, v, g2);

            phState = 2; // fluid phase at Pkb
            if (Tk < 673.)
            {
                a = a0 + a1*(673.-Tk) + a2*pow((673.-Tk),2.) + a3*pow((673.-Tk),3.);
                da = -a1 - 2.*a2*(673.-Tk) - 3.*a3*pow((673.-Tk),2.);
                d2a = 2.*a2 + 6.*a3*(673.-Tk);
            }
            else
            {
                a = a0 + a4*(673.-Tk) + a5*pow((673.-Tk),2.) + a6*pow((673.-Tk),3.);
                da = -a4 - 2.*a5*(673.-Tk) - 3.*a6*pow((673.-Tk),2.);
                d2a = 2.*a5 + 6.*a6*(673.-Tk);
            }
            VolumeFugacity(phState, Pkb, p0, a, b, c, d, e, vc, g3);  // corrected from (phState, Pkb, psat ...)
            fc = g1/g2*g3;  // this formula needs checking
        }
    }

    // pressure and density derivatives (MRK part)
    rho = (1./vc);
    dpr = (RR*Tk)/(pow((1./rho-b),2.)*pow(rho,2.)) - a/((1./rho+b)*pow(Tk,0.5))
            - a/(rho*pow((1./rho+b),2.)*pow(Tk,0.5));
    dprr = (2.*RR*Tk)/(pow((1./rho-b),3.)*pow(rho,4.)) - (2.*RR*Tk)/(pow((1./rho-b),2.)*pow(rho,3.))
            - (2.*a)/(pow(rho,3.)*pow((1./rho+b),3.)*pow(Tk,0.5));
    dpt = RR/(1./rho-b) - (da*rho)/((1./rho+b)*pow(Tk,0.5)) + (0.5*a*rho)/((1./rho+b)*pow(Tk,1.5));
    dptt = - (d2a*rho)/((1./rho+b)*pow(Tk,0.5)) + (da*rho)/((1./rho+b)*pow(Tk,1.5))
            - (0.75*a*rho)/((1./rho+b)*pow(Tk,2.5));
    dprt = RR/(pow((1./rho-b),2.)*pow(rho,2.)) - da/((1./rho+b)*pow(Tk,0.5))
            - da/(rho*pow((1./rho+b),2.)*pow(Tk,0.5)) + (0.5*a)/((1./rho+b)*pow(Tk,1.5))
            + (0.5*a)/(rho*pow((1./rho+b),2.)*pow(Tk,1.5));
    drp = (1./dpr);
    drpp = (-1.)*dprr*pow(dpr,-3.);
    drt = (-1.)*(1./dpr)*dpt;
    drtt = (dpt - dptt + dpt*(dprt-dpr)*pow(dpr,-1.))*pow(dpr,-1.)
            + (dprt - dpt*dprr*pow(dpr,-1.))*pow(dpr,-2.)*dpt;
    drtp = - pow(dpr,-1.) - (dprt-dpr)*pow(dpr,-2.) + dprr*dpt*pow(dpr,-3.);

    // thermodynamic properties (MRK part)
    z = (Pkb*vc)/(RR*Tk);
    grs = R_CONST*Tk*log(fc);
    hrs = ( (z-1.) - 1./(b*RR*Tk)*(pow(Tk,0.5)*da-(3./2.)*a/pow(Tk,0.5))*log(vc/(vc+b)) )*R_CONST*Tk;
    srs = (hrs-grs)/Tk;
    cvrs = ( -1./(b*RR)*(-da/pow(Tk,0.5)+pow(Tk,0.5)*d2a+(3.*a)/(4.*pow(Tk,1.5))))*log(vc/(vc+b) )*R_CONST;
    dpv = (- RR*Tk/pow((vc-b),2.) + a/(pow(vc,2.)*(vc+b)*pow(Tk,0.5))
            + a/(vc*pow((vc+b),2.)*pow(Tk,0.5)))*(1000.);
    dvt = (-1.)*(1./dpv)*(dpt*1000.);
    cprs = cvrs + Tk*(dpt*1000.)*dvt - R_CONST;

    // copy results
    Fugpure[j][0] = fc;  // fugacity coefficient
    Fugpure[j][1] = grs;
    Fugpure[j][2] = hrs;
    Fugpure[j][3] = srs;
    Fugpure[j][4] = vc;
    Fugpure[j][5] = cprs;
    Rho[j][0] = rho*(0.1);       // mol cm-3
    Rho[j][1] = drt*(0.1);       // mol cm-3 K-1
    Rho[j][2] = drtt*(0.1);      // mol cm-3 K-2
    Rho[j][3] = drp*(0.001);     // mol cm-3 MPa-1
    Rho[j][4] = drpp*(0.00001);  // mol cm-3 MPa-2
    Rho[j][5] = drtp*(0.001);    // mol cm-3 K-1 MPa-1
    Rho[j][6] = dpr*(1000.);     // MPa cm3 mol-1
    Rho[j][7] = dprr*(10000.);   // MPa cm6 mol-2
    Rho[j][8] = dpt*(100.);      // MPa K-1
    Rho[j][9] = dptt*(100.);     // MPa K-2
    Rho[j][10] = dprt*(1000.);   // MPa cm3 mol-1 K-1

    return 0;
}



/// calculates fugacity and state functions of CO2
long int TCORKcalc::FugacityCO2( long int j )
{
    long int phState;
    double a0, a1, a2, a, b, c0, c1, c, d0, d1, d, e0, e, da, d2a, p0, vc, fc, rho, z;
    double grs, hrs, srs, cprs, cvrs, dpr, dprr, dpt, dptt, dprt,
            drt, drtt, drp, drpp, drtp, dpv, dvt;

    a0 = 659.8;
    a1 = 0.21078;
    a2 = -6.3976e-4;
    b = 3.057;
    c0 = -1.78198e-1;
    c1 = 2.45317e-5;
    d0 = 5.40776e-3;
    d1 = -1.59046e-6;
    e0 = 0.;
    a = a0 + a1*Tk + a2*pow(Tk,2.);
    da = a1 + 2.*a2*Tk;
    d2a = 2.*a2;
    c = c0 + c1*Tk;
    d = d0 + d1*Tk;
    e = e0;
    p0 = 5.;

    // molar volume and fugacity coefficient
    phState = 1;
    VolumeFugacity(phState, Pkb, p0, a, b, c, d, e, vc, fc);

    // pressure and density derivatives (MRK part)
    rho = (1./vc);
    dpr = (RR*Tk)/(pow((1./rho-b),2.)*pow(rho,2.)) - a/((1./rho+b)*pow(Tk,0.5))
            - a/(rho*pow((1./rho+b),2.)*pow(Tk,0.5));
    dprr = (2.*RR*Tk)/(pow((1./rho-b),3.)*pow(rho,4.)) - (2.*RR*Tk)/(pow((1./rho-b),2.)*pow(rho,3.))
            - (2.*a)/(pow(rho,3.)*pow((1./rho+b),3.)*pow(Tk,0.5));
    dpt = RR/(1./rho-b) - (da*rho)/((1./rho+b)*pow(Tk,0.5)) + (0.5*a*rho)/((1./rho+b)*pow(Tk,1.5));
    dptt = - (d2a*rho)/((1./rho+b)*pow(Tk,0.5)) + (da*rho)/((1./rho+b)*pow(Tk,1.5))
            - (0.75*a*rho)/((1./rho+b)*pow(Tk,2.5));
    dprt = RR/(pow((1./rho-b),2.)*pow(rho,2.)) - da/((1./rho+b)*pow(Tk,0.5))
            - da/(rho*pow((1./rho+b),2.)*pow(Tk,0.5)) + (0.5*a)/((1./rho+b)*pow(Tk,1.5))
            + (0.5*a)/(rho*pow((1./rho+b),2.)*pow(Tk,1.5));
    drp = (1./dpr);
    drpp = (-1.)*dprr*pow(dpr,-3.);
    drt = (-1.)*(1./dpr)*dpt;
    drtt = (dpt - dptt + dpt*(dprt-dpr)*pow(dpr,-1.))*pow(dpr,-1.)
            + (dprt - dpt*dprr*pow(dpr,-1.))*pow(dpr,-2.)*dpt;
    drtp = - pow(dpr,-1.) - (dprt-dpr)*pow(dpr,-2.) + dprr*dpt*pow(dpr,-3.);

    // thermodynamic properties (MRK part)
    z = (Pkb*vc)/(RR*Tk);
    grs = R_CONST*Tk*log(fc);
    hrs = ( (z-1.) - 1./(b*RR*Tk)*(pow(Tk,0.5)*da-(3./2.)*a/pow(Tk,0.5))*log(vc/(vc+b)) )*R_CONST*Tk;
    srs = (hrs-grs)/Tk;
    cvrs = ( -1./(b*RR)*(-da/pow(Tk,0.5)+pow(Tk,0.5)*d2a+(3.*a)/(4.*pow(Tk,1.5))))*log(vc/(vc+b) )*R_CONST;
    dpv = (- RR*Tk/pow((vc-b),2.) + a/(pow(vc,2.)*(vc+b)*pow(Tk,0.5))
            + a/(vc*pow((vc+b),2.)*pow(Tk,0.5)))*(1000.);
    dvt = (-1.)*(1./dpv)*(dpt*1000.);
    cprs = cvrs + Tk*(dpt*1000.)*dvt - R_CONST;

    // copy results
    Fugpure[j][0] = fc;     // fugacity coefficient
    Fugpure[j][1] = grs;
    Fugpure[j][2] = hrs;
    Fugpure[j][3] = srs;
    Fugpure[j][4] = vc;
    Fugpure[j][5] = cprs;
    Rho[j][0] = rho*(0.1);       // mol cm-3
    Rho[j][1] = drt*(0.1);       // mol cm-3 K-1
    Rho[j][2] = drtt*(0.1);      // mol cm-3 K-2
    Rho[j][3] = drp*(0.001);     // mol cm-3 MPa-1
    Rho[j][4] = drpp*(0.00001);  // mol cm-3 MPa-2
    Rho[j][5] = drtp*(0.001);    // mol cm-3 K-1 MPa-1
    Rho[j][6] = dpr*(1000.);     // MPa cm3 mol-1
    Rho[j][7] = dprr*(10000.);   // MPa cm6 mol-2
    Rho[j][8] = dpt*(100.);      // MPa K-1
    Rho[j][9] = dptt*(100.);     // MPa K-2
    Rho[j][10] = dprt*(1000.);   // MPa cm3 mol-1 K-1

    return 0;
}



/// calculates fugacity and state functions of fluids other than H2O and CO2
long int TCORKcalc::FugacityCorresponding( long int j )
{
    double a0, a1, a, b0, b, c0, c1, c, d0, d1, d, tcr, pcr, da, dc, dd, vc, fc, rtlnf, rho;
    double grs, hrs, srs, cprs, drt, drtt, drp, drpp, drtp,
            dvt, dvtt, dvp, dvpp, dvtp;

    a0 = 5.45963e-5;
    a1 = -8.6392e-6;
    b0 = 9.18301e-4;
    c0 = -3.30558e-5;
    c1 = 2.30524e-6;
    d0 = 6.93054e-7;
    d1 = -8.38293e-8;
    tcr = Eosparm[j][0];
    pcr = Eosparm[j][1]/1000.;  // kbar
    a = a0*pow(tcr,2.)*sqrt(tcr)/pcr + a1*tcr*sqrt(tcr)/pcr*Tk;
    b = b0*tcr/pcr;
    c = c0*tcr/(pcr*sqrt(pcr)) + c1/(pcr*sqrt(pcr))*Tk;
    d = d0*tcr/pow(pcr,2.) + d1/pow(pcr,2.)*Tk;
    da = a1*pow(tcr,1.5)/pcr;
    dc = c1/pow(pcr,1.5);
    dd = d1/pow(pcr,2.);

    // molar volume and fugacity coefficient
    vc = RR*Tk/Pkb + b - a*RR*sqrt(Tk)/((RR*Tk+b*Pkb)*(RR*Tk+2.*b*Pkb)) + c*sqrt(Pkb) + d*Pkb;
    rtlnf = RR*Tk*log(1000*Pkb) + b*Pkb + a/(b*sqrt(Tk))*(log(RR*Tk+b*Pkb) - log(RR*Tk+2*b*Pkb))
            + (2./3.)*c*Pkb*sqrt(Pkb) + (d/2.)*pow(Pkb,2.);
    fc = exp(rtlnf/(RR*Tk))/(1000*Pkb);

    // volume and density derivatives
    dvt = ( RR/Pkb - (da*RR*pow(Tk,0.5))/((RR*Tk+b*Pkb)*(RR*Tk+2.*b*Pkb))
            - (a*RR)/(2.*pow(Tk,0.5)*(RR*Tk+b*Pkb)*(RR*Tk+2.*b*Pkb))
            + (a*pow(RR,2.)*pow(Tk,0.5))/(pow((RR*Tk+b*Pkb),2.)*(RR*Tk+2.*b*Pkb))
            + (a*pow(RR,2.)*pow(Tk,0.5))/((RR*Tk+b*Pkb)*pow((RR*Tk+2.*b*Pkb),2.))
            + dc*pow(Pkb,0.5) + dd*Pkb );
    dvtt = ( - (da*RR)/(pow(Tk,0.5)*(RR*Tk+b*Pkb)*(RR*Tk+2.*b*Pkb))
            + (2.*da*pow(RR,2.)*pow(Tk,0.5))/(pow((RR*Tk+b*Pkb),2.)*(RR*Tk+2.*b*Pkb))
            + (2.*da*pow(RR,2.)*pow(Tk,0.5))/((RR*Tk+b*Pkb)*pow((RR*Tk+2.*b*Pkb),2.))
            + (a*RR)/(4.*pow(Tk,1.5)*(RR*Tk+b*Pkb)*(RR*Tk+2.*b*Pkb))
            + (a*pow(RR,2.))/(pow(Tk,0.5)*pow((RR*Tk+b*Pkb),2.)*(RR*Tk+2.*b*Pkb))
            + (a*pow(RR,2.))/(pow(Tk,0.5)*(RR*Tk+b*Pkb)*pow((RR*Tk+2.*b*Pkb),2.))
            - (2.*a*pow(RR,3.)*pow(Tk,0.5))/(pow((RR*Tk+b*Pkb),3.)*(RR*Tk+2.*b*Pkb))
            - (2.*a*pow(RR,3.)*pow(Tk,0.5))/(pow((RR*Tk+b*Pkb),2.)*pow((RR*Tk+2.*b*Pkb),2.))
            - (2.*a*pow(RR,3.)*pow(Tk,0.5))/((RR*Tk+b*Pkb)*pow((RR*Tk+2.*b*Pkb),3.)) );
    dvp = ( - RR*Tk/pow(Pkb,2.) + (a*RR*pow(Tk,0.5)*b)/(pow((RR*Tk+b*Pkb),2.)*(RR*Tk+2.*b*Pkb))
            + (2.*a*RR*pow(Tk,0.5)*b)/((RR*Tk+b*Pkb)*pow((RR*Tk+2.*b*Pkb),2.))
            + c/(2.*pow(Pkb,0.5)) + d ) / (1000.);
    dvpp = ( (2.*RR*Tk)/pow(Pkb,3.) - (2.*a*RR*pow(Tk,0.5)*pow(b,2.))/(pow((RR*Tk+b*Pkb),3.)*(RR*Tk+2.*b*Pkb))
            - (4.*a*RR*pow(Tk,0.5)*pow(b,2.))/(pow((RR*Tk+b*Pkb),2.)*pow((RR*Tk+2.*b*Pkb),2.))
            - (8.*a*RR*pow(Tk,0.5)*pow(b,2.))/((RR*Tk+b*Pkb)*pow((RR*Tk+2.*b*Pkb),3.))
            - c/(4.*pow(Pkb,1.5)) ) / (1.0e6);
    dvtp = ( - RR/pow(Pkb,2.) + (da*RR*pow(Tk,0.5)*b)/(pow((RR*Tk+b*Pkb),2.)*(RR*Tk+2.*b*Pkb))
            + (a*RR*b)/(2.*pow(Tk,0.5)*pow((RR*Tk+b*Pkb),2.)*(RR*Tk+2.*b*Pkb))
            - (2.*a*pow(RR,2.)*pow(Tk,0.5)*b)/(pow((RR*Tk+b*Pkb),3.)*(RR*Tk+2.*b*Pkb))
            - (3.*a*pow(RR,2.)*pow(Tk,0.5)*b)/(pow((RR*Tk+b*Pkb),2.)*pow((RR*Tk+2.*b*Pkb),2.))
            + (2.*da*RR*pow(Tk,0.5)*b)/((RR*Tk+b*Pkb)*pow((RR*Tk+2.*b*Pkb),2.))
            + (a*RR*b)/(pow(Tk,0.5)*(RR*Tk+b*Pkb)*pow((RR*Tk+2.*b*Pkb),2.))
            - (4.*a*pow(RR,2.)*pow(Tk,0.5)*b)/((RR*Tk+b*Pkb)*pow((RR*Tk+2.*b*Pkb),3.))
            + dc/(2.*pow(Pkb,2.)) + dd ) / (1000.);
    rho = (1./vc);
    drt = - dvt*pow(rho,2.);
    drp = - dvp*pow(rho,2.);
    drtt = - (dvtt*pow(rho,3.)-2.*pow(drt,2.))/rho;
    drpp = - (dvpp*pow(rho,3.)-2.*pow(drp,2.))/rho;
    drtp = - (dvtp*pow(rho,3.)-2.*drt*drp)/rho;

    // thermodynamic properties
    grs = R_CONST*Tk*log(fc);
    hrs = ( (pow(Tk,0.5)*da*log(RR*Tk+2.*b*Pkb))/b - (pow(Tk,0.5)*da*log(RR*Tk+b*Pkb))/b
            - (3.*a*log(RR*Tk+2.*b*Pkb))/(2.*pow(Tk,0.5)*b) + (3.*a*log(RR*Tk+b*Pkb))/(2.*pow(Tk,0.5)*b)
            - (pow(Tk,0.5)*a*RR)/(b*(RR*Tk+b*Pkb)) + (pow(Tk,0.5)*a*RR)/(b*(RR*Tk+2.*b*Pkb))
            - (2.*Tk*dc*pow(Pkb,1.5))/3. - (Tk*dd*pow(Pkb,2.))/2. + b*Pkb
            + (2.*c*pow(Pkb,1.5))/3. + (d*pow(Pkb,2.))/2. )*1000.;
    srs = (hrs-grs)/Tk;
    cprs = ( - (2.*pow(Tk,0.5)*da*RR)/(b*(RR*Tk+b*Pkb)) + (pow(Tk,0.5)*a*pow(RR,2.))/(b*pow((RR*Tk+b*Pkb),2.))
            + (2.*pow(Tk,0.5)*da*RR)/(b*(RR*Tk+2.*b*Pkb)) - (a*RR)/(pow(Tk,0.5)*b*(RR*Tk+2.*b*Pkb))
            - (da*log(RR*Tk+2.*b*Pkb))/(pow(Tk,0.5)*b) + (da*log(RR*Tk+b*Pkb))/(pow(Tk,0.5)*b)
            + (3.*a*log(RR*Tk+2.*b*Pkb))/(4.*pow(Tk,1.5)*b) - (3.*a*log(RR*Tk+b*Pkb))/(4.*pow(Tk,1.5)*b)
            - (pow(Tk,0.5)*a*pow(RR,2.))/(b*pow((RR*Tk+2.*b*Pkb),2.)) + (a*RR)/(pow(Tk,0.5)*b*(RR*Tk+b*Pkb)) )*1000.;

    // copy results
    Fugpure[j][0] = fc;     // fugacity coefficient
    Fugpure[j][1] = grs;
    Fugpure[j][2] = hrs;
    Fugpure[j][3] = srs;
    Fugpure[j][4] = vc;
    Fugpure[j][5] = cprs;
    Rho[j][0] = rho*(0.1);   // mol cm-3
    Rho[j][1] = drt*(0.1);   // mol cm-3 K-1
    Rho[j][2] = drtt*(0.1);  // mol cm-3 K-2
    Rho[j][3] = drp;         // mol cm-3 MPa-1
    Rho[j][4] = drpp*(10.);  // mol cm-3 MPa-2
    Rho[j][5] = drtp;        // mol cm-3 K-1 MPa-1
    Rho[j][6] = 1.0;         // MPa cm3 mol-1
    Rho[j][7] = 1.0;         // MPa cm6 mol-2
    Rho[j][8] = 1.0;         // MPa K-1
    Rho[j][9] = 1.0;         // MPa K-2
    Rho[j][10] = 1.0;        // MPa cm3 mol-1 K-1

    return 0;
}



/// calculate volume and fugacity coefficient
long int TCORKcalc::VolumeFugacity( long int phState, double pp, double p0, double a, double b, double c,
        double d, double e, double &vol, double &fc )
{
    double cb, cc, cd, v1, v2, v3, vmrk, vvir, lng, lnvir;

    cb = (-1.)*RR*Tk/pp;
    cc = (-1.)*(b*RR*Tk+pow(b,2.)*pp-a/sqrt(Tk))/pp;
    cd = (-1.)*a*b/(sqrt(Tk)*pp);
    Cardano(cb, cc, cd, v1, v2, v3);

    if ( phState == 1 )   // vapor root
    {
        if (v1>0.)
            vmrk = v1;
        else if (v2>0.)
            vmrk = v2;
        else
            vmrk = v3;
        if (v2>vmrk)
            vmrk = v2;
        if (v3>vmrk)
            vmrk=v3;
    }

    else   // liquid root
    {
        if (v1>0.)
            vmrk = v1;
        else if (v2>0.)
            vmrk = v2;
        else
            vmrk = v3;
        if ( (v2<vmrk) && (v2>0.) )
            vmrk = v2;
        if ( (v3<vmrk) && (v2>0.) )
            vmrk = v3;
    }

    // calculate fugacity coefficient
    lng = pp*vmrk/(RR*Tk) - 1. - log((vmrk-b)*pp/(RR*Tk)) - a/(b*RR*Tk*sqrt(Tk))*log(1.+b/vmrk);
    if (pp>p0)
    {
        vvir = c*sqrt(pp-p0) + d*(pp-p0) + e*sqrt(sqrt(pp-p0));
        lnvir = ((2./3.)*c*(pp-p0)*sqrt(pp-p0) + d/2.*pow((pp-p0),2.) + (4./5.)*e*(pp-p0)*sqrt(sqrt(pp-p0)))/(RR*Tk);
    }
    else
    {
        vvir = 0.;
        lnvir = 0.;
    }

    vol = vmrk + vvir;
    lng = lng + lnvir;
    fc = exp(lng);
    return 0;
}



/// finds roots of cubic equation
long int TCORKcalc::Cardano(double cb, double cc, double cd, double &v1, double &v2, double &v3)
{
    double co, cp, cq, cr, cp2, cq3;

    cp = (2.*pow(cb,3.)-9.*cb*cc+27.*cd)/54.;
    cq = (pow(cb,2.)-3.*cc)/9.;
    cp2 = pow(cp,2.);
    cq3 = pow(cq,3.);

    if ( (cp2-cq3)>0. )  // one real root
    {
        cr = -cp/fabs(cp)*pow(fabs(cp)+sqrt(cp2-cq3),1./3.);
        if (cr!=0)
        {
            v1 = cr + cq/cr - cb/3.;
            v2 = cr + cq/cr - cb/3.;
            v3 = cr + cq/cr - cb/3.;
        }
        else
        {
            v1 = -cb/3.;
            v2 = -cb/3.;
            v3 = -cb/3.;
        }
    }

    else  //three real roots
    {
        co = atan(sqrt(1.-cp2/(cq3))/(cp/sqrt(cq3)));
        if (co<0.)
            co = co + 3.1415927;
        v1 = (-2.)*sqrt(cq)*cos(co/3.) - cb/3.;
        v2 = (-2.)*sqrt(cq)*cos((co+2.*3.1415927)/3.) - cb/3.;
        v3 = (-2.)*sqrt(cq)*cos((co-2.*3.1415927)/3.) - cb/3.;
    }

    return 0;
}



/// calculates excess state functions of the fluid mixture
long int TCORKcalc::ResidualFunct()
{
    long int i, j;
    double sumphi, vi, vj, dvi, dvj, d2vi, d2vj, dvip, dvjp, Y, dY, d2Y, dYp,
            rhoi, drti, drtti, drpi, rhoj, drtj, drttj, drpj;
    double gex, dgex, d2gex, dgexp, sex, hex, cpex, vex, gph, sph, hph, cpph, vph;
    gex = 0.; dgex = 0.; d2gex = 0.; dgexp = 0.;
    sex = 0.; hex = 0.; vex = 0.; cpex = 0.;
    gph = 0.; sph = 0.; hph = 0.; cpph = 0.; vph = 0.;
    sumphi = 0.; Y = 0.; dY = 0.; d2Y = 0.; dYp = 0.;

    // phi values (and derivatives)
    for (i=0; i<NComp; i++)
    {
        rhoi = Rho[i][0];
        drti = Rho[i][1];
        drtti = Rho[i][2];
        drpi = Rho[i][3];
        vi = 1./rhoi;
        dvi = - pow(rhoi,-2.)*drti;
        d2vi = 2.*pow(rhoi,-3.)*pow(drti,2.) - pow(rhoi,-2.)*drtti;
        dvip = - pow(rhoi,-2.)*drpi;
        sumphi += x[i]*vi;
        Y += x[i]*vi;
        dY += x[i]*dvi;
        d2Y += x[i]*d2vi;
        dYp += x[i]*dvip;
    }

    for (i=0; i<NComp; i++)
    {
        rhoi = Rho[i][0];
        drti = Rho[i][1];
        drtti = Rho[i][2];
        drpi = Rho[i][3];
        vi = 1./rhoi;
        dvi = - pow(rhoi,-2.)*drti;
        d2vi = 2.*pow(rhoi,-3.)*pow(drti,2.) - pow(rhoi,-2.)*drtti;
        dvip = - pow(rhoi,-2.)*drpi;
        phi[i] = x[i]*vi/sumphi;
        dphi[i] = x[i]*(dvi*Y-vi*dY)/pow(Y,2.);
        dphip[i] = x[i]*(dvip*Y-vi*dYp)/pow(Y,2.);
        d2phi[i] = x[i]*(d2vi*Y+dvi*dY)/pow(Y,2.) - x[i]*(dvi*Y)*(2.*dY)/pow(Y,3.)
            - x[i]*(dvi*dY+vi*d2Y)/pow(Y,2.) + x[i]*(vi*dY)*(2.*dY)/pow(Y,3.);
    }

    // interaction parameters
    for (i=0; i<NComp; i++)
    {
        for (j=i+1; j<NComp; j++)
        {
            rhoi = Rho[i][0];
            drti = Rho[i][1];
            drtti = Rho[i][2];
            drpi = Rho[i][3];
            rhoj = Rho[j][0];
            drtj = Rho[j][1];
            drttj = Rho[j][2];
            drpj = Rho[j][3];

            vi = 1./rhoi;
            dvi = - pow(rhoi,-2.)*drti;
            d2vi = 2.*pow(rhoi,-3.)*pow(drti,2.) - pow(rhoi,-2.)*drtti;
            dvip = - pow(rhoi,-2.)*drpi;
            vj = 1./rhoj;
            dvj = - pow(rhoj,-2.)*drtj;
            d2vj = 2.*pow(rhoj,-3.)*pow(drtj,2.) - pow(rhoj,-2.)*drttj;
            dvjp = - pow(rhoj,-2.)*drpj;

            B[i][j] = A[i][j]*2*Y/(vi*vj);
            dB[i][j] = A[i][j]*(2.*dY*vi*vj-2.*Y*(dvi*vj+vi*dvj))/(pow(vi,2.)*pow(vj,2.));
            d2B[i][j] = A[i][j]*((2.*d2Y)*(vi*vj) + (2.*dY)*(dvi*vj+vi*dvj))/pow((vi*vj),2.)
                - A[i][j]*((2.*dY)*(vi*vj))*(2.*(dvi*vj+vi*dvj))/pow((vi*vj),3.)
                - A[i][j]*((2.*dY)*(dvi*vj+vi*dvj) + (2.*Y)*(d2vi*vj+2.*dvi*dvj+vi*d2vj))/pow((vi*vj),2.)
                + A[i][j]*((2.*Y)*(dvi*vj+vi*dvj))*(2.*(dvi*vj+vi*dvj))/pow((vi*vj),3.);
            dBp[i][j] = A[i][j]*(2.*dYp*vi*vj-2.*Y*(dvip*vj+vi*dvjp))/(pow(vi,2.)*pow(vj,2.));
        }
    }

    // excess properties
    for (i=0; i<NComp; i++)
    {
        for (j=i+1; j<NComp; j++)
        {
            gex += phi[i]*phi[j]*B[i][j];
            dgex += dphi[i]*phi[j]*B[i][j] + phi[i]*dphi[j]*B[i][j] + phi[i]*phi[j]*dB[i][j];
            d2gex +=d2phi[i]*phi[j]*B[i][j] + 2.*dphi[i]*dphi[j]*B[i][j] + 2.*dphi[i]*phi[j]*dB[i][j]
                + phi[i]*d2phi[j]*B[i][j] + 2.*phi[i]*dphi[j]*dB[i][j] + phi[i]*phi[j]*d2B[i][j];
            dgexp += dphip[i]*phi[j]*B[i][j] + phi[i]*dphip[j]*B[i][j] + phi[i]*phi[j]*dBp[i][j];
        }
    }

    // thermodynamic properties
    sex = - dgex;
    hex = gex + Tk*sex;
    cpex = - d2gex*Tk;
    vex = dgexp;
    Gex = gex*(10.);
    Sex = sex*(10.);
    Hex = hex*(10.);
    CPex = cpex*(10.);
    Vex = vex;
    for (j=0; j<NComp; j++)
    {
        gph += x[j]*Fugpure[j][1];
        hph += x[j]*Fugpure[j][2];
        sph += x[j]*Fugpure[j][3];
        vph += x[j]*Fugpure[j][4];
        cpph += x[j]*Fugpure[j][5];
    }
    Grs = gph + Gex;
    Hrs = hph + Hex;
    Srs = sph + Sex;
    Vrs = vph + Vex;
    CPrs = cpph + CPex;

    return 0;
}



#ifndef IPMGEMPLUGIN

/// Calculates properties of pure fluids when called from DCthermo
long int TCORKcalc::CORKCalcFugPure( double Tmin, float *Cpg, double *FugProps )
{
        long int retCode = 0;
        double Coeff[7];

        for( int ii=0; ii<7; ii++ )
                Coeff[ii] = (double)Cpg[ii];

        if( (Tk >= Tmin) && (Tk < 1e4) && (Pbar >= 1e-5) && (Pbar < 1e5) )
        {
                retCode = FugacityPT( 0, Coeff );
                for( int i=0; i<6; i++ )
                        FugProps[i] = Fugpure[0][i];
                return retCode;
        }

        else
        {
                for( int i=1; i<6; i++ )
                        FugProps[i] = 0.;
                FugProps[0] = 1.;
                FugProps[4] = 8.31451*Tk/Pbar;
                return -1;
        }
}

#endif





//=======================================================================================================
// Sterner-Pitzer (STP) model for fluid mixtures
// References: Sterner and Pitzer (1994)
// (c) TW December 2010
//=======================================================================================================


// Constructor
TSTPcalc::TSTPcalc( long int NCmp, double Pp, double Tkp, char Eos_Code ):
    TSolMod( NCmp, '6', Tkp, Pp )
{

    UpdateTauP();
    Pparc = 0;
    alloc_internal();
    set_internal();
    EosCode[0] = Eos_Code;

    // provisional
    if ( EosCode[0] == CEM_H2O_ )  // H2O
        Mw[0] = 18.015268;
    if ( EosCode[0] == CEM_CO2_ )  // CO2
        Mw[0] = 44.0098;
}



TSTPcalc::TSTPcalc( SolutionData *sd ):
                TSolMod( sd )
{
    UpdateTauP();
    Pparc = aPparc;
    alloc_internal();
    set_internal();
    for( long int j=0; j<NComp; j++ )
        EosCode[j] = sd->TP_Code[j][3];

    // provisional
    for ( long int j=0; j<NComp; j++ )
    {
        if ( EosCode[j] == CEM_H2O_ )  // H2O
            Mw[j] = 18.015268;
        if ( EosCode[j] == CEM_CO2_ )  // CO2
            Mw[j] = 44.0098;
    }

}



TSTPcalc::~TSTPcalc()
{
    free_internal();
}



/// allocate work arrays for pure fluid and fluid mixture properties
void TSTPcalc::alloc_internal()
{
    long int k;

    EosCode = new char[NComp];
    Tc = new double [NComp];
    Pc = new double [NComp];
    Psat = new double [NComp];
    Rhol = new double [NComp];
    Rhov = new double [NComp];
    Mw = new double [NComp];
    Phi = new double [NComp];
    dPhiD = new double [NComp];
    dPhiDD = new double [NComp];
    dPhiT = new double [NComp];
    dPhiTT = new double [NComp];
    dPhiDT = new double [NComp];
    dPhiDDD = new double [NComp];
    dPhiDDT = new double [NComp];
    dPhiDTT = new double [NComp];
    phi = new double [NComp];
    dphi = new double [NComp];
    d2phi = new double [NComp];
    dphip = new double [NComp];
    lng = new double [NComp];
    Fugpure = new double [NComp][7];
    Rho = new double [NComp][11];
    cfh = new double *[10];
    cfc = new double *[10];
    A = new double *[NComp];
    W = new double *[NComp];
    B = new double *[NComp];
    dB = new double *[NComp];
    d2B = new double *[NComp];
    dBp = new double *[NComp];

    for (k=0; k<10; k++)
    {
        cfh[k] = new double [6];
        cfc[k] = new double [6];
    }

    for (k=0; k<NComp; k++)
    {
        A[k] = new double [NComp];
        W[k] = new double [NComp];
        B[k] = new double [NComp];
        dB[k] = new double [NComp];
        d2B[k] = new double [NComp];
        dBp[k] = new double [NComp];
    }

}



void TSTPcalc::free_internal()
{
    long int k;

    for (k=0; k<10; k++)
    {
        delete[]cfh[k];
        delete[]cfc[k];
    }

    for (k=0; k<NComp; k++)
    {
        delete[]A[k];
        delete[]W[k];
        delete[]B[k];
        delete[]dB[k];
        delete[]d2B[k];
        delete[]dBp[k];
    }

    delete[]EosCode;
    delete[]Tc;
    delete[]Pc;
    delete[]Psat;
    delete[]Rhol;
    delete[]Rhov;
    delete[]Mw;
    delete[]Phi;
    delete[]dPhiD;
    delete[]dPhiDD;
    delete[]dPhiT;
    delete[]dPhiTT;
    delete[]dPhiDT;
    delete[]dPhiDDD;
    delete[]dPhiDDT;
    delete[]dPhiDTT;
    delete[]phi;
    delete[]dphi;
    delete[]d2phi;
    delete[]dphip;
    delete[]lng;
    delete[]Fugpure;
    delete[]Rho;
    delete[]cfh;
    delete[]cfc;
    delete[]A;
    delete[]W;
    delete[]B;
    delete[]dB;
    delete[]d2B;
    delete[]dBp;
}



/// set model specific parameters
void TSTPcalc::set_internal()
{
    long int j, k;

    RC = R_CONST;  // gas constant (bar)
    RR = 0.00831451;  // gas constant (kbar)
    TMIN = 200.;
    TMAX = 2000.;
    PMIN = 1.0e-20;
    PMAX = 10000.;

    // EoS coefficients (temperature dependence)
    double cfh_[10][6] = { { 0.0, 0.0, 0.24657688e6, 0.51359951e2, 0.0, 0.0 },
                           { 0.0, 0.0, 0.58638965, -0.28646939e-2, 0.31375577e-4, 0.0 },
                           { 0.0, 0.0, -0.62783840e1, 0.14791599e-1, 0.35779579e-3, 0.15432925e-7 },
                           { 0.0, 0.0, 0.0, -0.42719875, -0.16325155e-4, 0.0 },
                           { 0.0, 0.0, 0.56654978e4, -0.16580167e2, 0.76560762e-1, 0.0 },
                           { 0.0, 0.0, 0.0, 0.10917883, 0.0, 0.0 },
                           { 0.38878656e13, -0.13494878e9, 0.30916564e6, 0.75591105e1, 0.0, 0.0 },
                           { 0.0, 0.0, -0.65537898e5, 0.18810675e3, 0.0, 0.0 },
                           { -0.14182435e14, 0.18165390e9, -0.19769068e6, -0.23530318e2, 0.0, 0.0 },
                           { 0.0, 0.0, 0.92093375e5, 0.12246777e3, 0.0, 0.0 } };


    double cfc_[10][6] = { { 0.0, 0.0, 0.18261340e7, 0.79224365e2, 0.0, 0.0 },
                           { 0.0, 0.0, 0.0, 0.66560660e-4, 0.57152798e-5, 0.30222363e-9 },
                           { 0.0, 0.0, 0.0, 0.59957845e-2, 0.71669631e-4, 0.62416103e-8 },
                           { 0.0, 0.0, -0.13270279e1, -0.15210731, 0.53654244e-3, -0.71115142e-7 },
                           { 0.0, 0.0, 0.12456776, 0.49045367e1, 0.98220560e-2, 0.55962121e-5 },
                           { 0.0, 0.0, 0.0, 0.75522299, 0.0, 0.0 },
                           { -0.39344644e12, 0.90918237e8, 0.42776716e6, -0.22347856e2, 0.0, 0.0 },
                           { 0.0, 0.0, 0.40282608e3, 0.11971627e3, 0.0, 0.0 },
                           { 0.0, 0.22995650e8, -0.78971817e5, -0.63376456e2, 0.0, 0.0 },
                           { 0.0, 0.0, 0.95029765e5, 0.18038071e2, 0.0, 0.0 } };

    // transfer coefficients
    for (j=0; j<10; j++)
    {
            for (k=0; k<6; k++)
            {
                    cfh[j][k] = cfh_[j][k];
                    cfc[j][k] = cfc_[j][k];
            }
    }
}



long int TSTPcalc::UpdateTauP()
{
    Pkbar = Pbar/1000.;
    Pkb = Pbar/1000.;
    Pmpa = Pbar/10.;

    return 0;
}



/// Calculates pure species properties (pure fugacities)
long int TSTPcalc::PureSpecies()
{
    long int j, iErr = 0;

    // error conditions
    // 1: input temperature/pressure outside valid range
    // 2: density iteration not converged
    // 3: wrong EoS subroutine code

    for( j=0; j<NComp; j++ )
    {
        // Calling STP EoS for pure fugacity
        iErr =  FugacityPT( j, aDCc+j*NP_DC );
        aGEX[j] = log( Fugpure[j][0] );
        Pparc[j] = Fugpure[j][0]*Pbar;  // fure fluid fugacity (required for performance)
        aVol[j] = Fugpure[j][4]*10.;  // molar volume of pure fluid component, J/bar to cm3
    } // j

    if ( iErr )
    {
            char buf[150];
            sprintf(buf, "STP fluid: calculation of pure fluid fugacity failed");
                    Error( "E71IPM IPMgamma: ",  buf );
    }

    return 0;
}



/// Calculates T,P corrected interaction parameters
long int TSTPcalc::PTparam()
{
    long int j, i, ip, i1, i2;
    double a;

    UpdateTauP();

    PureSpecies();

    // set all interaction parameters zero
    for( j=0; j<NComp; j++ )
    {
        for( i=0; i<NComp; i++ )
        {
            A[j][i] = 0.;
            W[j][i] = 0.;
            B[j][i] = 0.;
            dB[j][i] = 0.;
            d2B[j][i] = 0.;
            dBp[j][i] = 0.;
        }
    }

    // transfer interaction parameters that have non-standard value
    if( NPcoef > 0 )
    {
        for ( ip=0; ip<NPar; ip++ )
        {
            i1 = aIPx[MaxOrd*ip];
            i2 = aIPx[MaxOrd*ip+1];
            a = aIPc[NPcoef*ip];
            A[i1][i2] = a;
            A[i2][i1] = a;  // symmetric case
        }
    }

    return 0;
}



/// Calculates activity coefficients
long int TSTPcalc::MixMod()
{
    long int i, j, k;
    double dj, dk;
    double sumphi, lnGam, Gam, vi, vj, vk;

    // phi values
    sumphi = 0.;
    for (i=0; i<NComp; i++)
    {
            vi = Fugpure[i][4];
            sumphi = sumphi + x[i]*vi;
    }
    for (i=0; i<NComp; i++)
    {
            phi[i] = x[i]*Fugpure[i][4]/sumphi;
    }

    // interaction parameters
    for (i=0; i<NComp; i++)
    {
            for (j=i+1; j<NComp; j++)
            {
                    vi = Fugpure[i][4];
                    vj = Fugpure[j][4];
                    W[i][j] = A[i][j]*(vi+vj)/(vi*vj);
            }
    }

    // activity coefficients
    for (i=0; i<NComp; i++)
    {
            lnGam = 0.;
            for (j=0; j<NComp; j++)
            {
                    for (k=j+1; k<NComp; k++)
                    {
                            vi = Fugpure[i][4];
                            vj = Fugpure[j][4];
                            vk = Fugpure[k][4];
                            if (i==j)
                                    dj = 1.;
                            else
                                    dj = 0.;
                            if (i==k)
                                    dk = 1.;
                            else
                                    dk = 0.;
                            lnGam = lnGam - (dj-phi[j])*(dk-phi[k])*W[j][k]*2.*vi/(vj+vk);
                    }
            }
            Gam = exp(lnGam/(RC*Tk));
            lnGamma[i] = log(Gam);
    }  // i

    return 0;
}



/// calculates excess properties
long int TSTPcalc::ExcessProp( double *Zex )
{
    long int iErr;

    iErr = ResidualFunct();

    if ( iErr )
    {
        char buf[150];
        sprintf(buf, "STP fluid: calculation failed");
        Error( "E71IPM IPMgamma: ",  buf );
    }

    Ars = Grs - Vrs*Pbar;
    Urs = Hrs - Vrs*Pbar;

    // assignments (residual functions)
    Zex[0] = Grs;
    Zex[1] = Hrs;
    Zex[2] = Srs;
    Zex[3] = CPrs;
    Zex[4] = Vrs;
    Zex[5] = Ars;
    Zex[6] = Urs;

    return iErr;
}



/// calculates ideal mixing properties
long int TSTPcalc::IdealProp( double *Zid )
{
    long int j;
    double s, sc, sp;

    s = 0.0;
    for ( j=0; j<NComp; j++ )
    {
        if ( x[j] > 1.0e-32 )
            s += x[j]*log(x[j]);
    }
    sc = (-1.)*R_CONST*s;
    sp = (-1.)*R_CONST*log(Pbar);
    Hid = 0.0;
    CPid = 0.0;
    Vid = 0.0;
    Sid = sc + sp;
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



/// high-level method to retrieve pure fluid properties
long int TSTPcalc::FugacityPT( long int j, double *EoSparam )
{
    long int iErr = 0;

    // reads EoS parameters from database into work array
    if( !EoSparam )
        return -1;  // Memory alloc error
    Tc[j] = EoSparam[0];  // critical temperature (K)
    Pc[j] = EoSparam[1]/(10.);  // critical pressure (MPa)

    // check for temperature/pressure limits
    if ( (Tk<TMIN) || (Tk>TMAX) || (Pmpa<PMIN) || (Pmpa>PMAX) )
    {
        iErr = 1;
        return iErr;
    }

    // select subroutines for different fluid types
    switch ( EosCode[j] )
    {
            case CEM_H2O_:  // H2O
                    iErr = FugacityH2O( j );
                    break;
            case CEM_CO2_: // CO2
                    iErr = FugacityCO2( j );
                    break;
            case CEM_GAS_:  // other fluids
            case CEM_CH4_:
            case CEM_N2_:
            case CEM_H2_:
            case CEM_O2_:
            case CEM_AR_:
            case CEM_PO_:
            case CEM_NP_:
                    iErr = FugacityCorresponding( j );
                    break;
            default:
                    iErr = 3;
                    break;
    }

    return iErr;
}



/// calculates fugacity and state functions of H2O
long int TSTPcalc::FugacityH2O( long int j )
{
    long int k, maxit, iErr = 0;
    double pmpa, tol, rhoguess, rho, pnew, pgrad, rhoold, rhonew, errnew, errold, step,
                    rhomin, rhomax, vol, lnfug, fug, fugc;
    double ar, dard, dardd, darddd, dart, dartt, dardt, darddt, dardtt, pig,
                    g, h, s, cp, cv, /*u,*/ a, dpr, dprr, dpt, dptt, dprt, drt, drtt, drp, drpp, drtp;

    pmpa = Pbar/10.;  // MPa
    tol = 1.e-10;
    maxit = 1000;
    rhomin = (1.e-20/1000.)/Mw[j];  // 1.e-20 kg m-3
    rhomax = (1800./1000.)/Mw[j];   // 1800 kg m-3

    // find initial density guess
    DensityGuess( j, rhoguess );
    Pressure( rhoguess, pnew, pgrad, cfh );

    rhoold = rhoguess;
    rhonew = rhoguess;
    errnew = pnew - pmpa;
    errold = pnew - pmpa;

    // find density by Newton iteration
    k = 0;
    do
    {
            k++;

            if ( k>= maxit )
            {
                    iErr = 2;
                    return iErr;
            }

            rhoold = rhonew;
            errold = errnew;
            Pressure( rhonew, pnew, pgrad, cfh );
            errnew = pnew - pmpa;
            step = - errnew/pgrad;
            rhonew = rhoold + step;

            if ( rhonew < rhomin )
                    rhonew = rhomin ;
            if ( rhonew > rhomax )
                    rhonew = rhomax;

    } while ( fabs(1.-pnew/pmpa) > tol );

    // calculate thermodynamic properties
    rho = rhonew;
    Helmholtz( j, rho, cfh );
    ar = RC*Tk*Phi[j];
    dard = RC*Tk*dPhiD[j];
    dardd = RC*Tk*dPhiDD[j];
    darddd = RC*Tk*dPhiDDD[j];
    dart = RC*(dPhiT[j]*Tk+Phi[j]);
    dartt = RC*(dPhiTT[j]*Tk+2.*dPhiT[j]);
    dardt = RC*(dPhiDT[j]*Tk+dPhiD[j]);
    darddt = RC*(dPhiDDT[j]*Tk+dPhiDD[j]);
    dardtt = RC*(dPhiDTT[j]*Tk+2.*dPhiDT[j]);

    vol = (1./rho)/10.;  // J bar-1
    lnfug = ( log(rho) + Phi[j] + pmpa/(rho*RC*Tk) ) + log(RC*Tk) - 1.;
    fug = exp(lnfug);
    fugc = fug/pmpa;
    pig = rho*RC*Tk;
    g = log(fugc)*RC*Tk;
    s = - dart + RC*log(pmpa/pig);
    h = g + Tk*s;
    a = ar - RC*Tk*log(pmpa/pig);
    cv = - Tk*dartt - RC;
    cp = cv + RC*pow((1.+(rho/RC)*dardt),2.) / (1.+(2.*rho)/(RC*Tk)*dard + pow(rho,2.)/(RC*Tk)*dardd);
    dpr = RC*Tk + 2.*rho*dard + pow(rho,2.)*dardd;
    dprr = 2.*dard + 4.*rho*dardd + pow(rho,2.)*darddd;
    dpt = rho*RC + pow(rho,2.)*dardt;
    dptt = pow(rho,2.)*dardtt;
    dprt = RC + 2.*rho*dardt + pow(rho,2.)*darddt;
    drp = (1./dpr);
    drpp = (-1.)*dprr*pow(dpr,-3.);
    drt = (-1.)*(1./dpr)*dpt;
    drtt = (dpt - dptt + dpt*(dprt-dpr)*pow(dpr,-1.))*pow(dpr,-1.)
            + (dprt - dpt*dprr*pow(dpr,-1.))*pow(dpr,-2.)*dpt;
    drtp = - pow(dpr,-1.) - (dprt-dpr)*pow(dpr,-2.) + dprr*dpt*pow(dpr,-3.);

    // copy results
    Fugpure[j][0] = fugc;  // fugacity coef.
    Fugpure[j][1] = g;     // Gres
    Fugpure[j][2] = h;     // Hres
    Fugpure[j][3] = s;     // Sres
    Fugpure[j][4] = vol;   // V
    Fugpure[j][5] = cp;    // CPres
    Fugpure[j][6] = cv;    // CVres
    Rho[j][0] = rho;      // mol cm-3
    Rho[j][1] = drt;      // mol cm-3 K-1
    Rho[j][2] = drtt;     // mol cm-3 K-2
    Rho[j][3] = drp;      // mol cm-3 MPa-1
    Rho[j][4] = drpp;     // mol cm-3 MPa-2
    Rho[j][5] = drtp;     // mol cm-3 K-1 MPa-1
    Rho[j][6] = dpr;      // MPa cm3 mol-1
    Rho[j][7] = dprr;     // MPa cm6 mol-2
    Rho[j][8] = dpt;      // MPa K-1
    Rho[j][9] = dptt;     // MPa K-2
    Rho[j][10] = dprt;    // MPa cm3 mol-1 K-1

    return iErr;
}



/// calculates fugacity and state functions of CO2
long int TSTPcalc::FugacityCO2( long int j )
{
    long int k, maxit, iErr = 0;
    double pmpa, tol, rhoguess, rho, pnew, pgrad, rhoold, rhonew, errnew, errold, step,
                    rhomin, rhomax, vol, lnfug, fug, fugc;
    double ar, dard, dardd, darddd, dart, dartt, dardt, darddt, dardtt, pig,
                    g, h, s, cp, cv, /*u,*/ a, dpr, dprr, dpt, dptt, dprt, drt, drtt, drp, drpp, drtp;

    pmpa = Pbar/10.;  // MPa
    tol = 1.e-10;
    maxit = 1000;
    rhomin = (1.e-20)/1000./Mw[j];  // 1.e-20 kg m-3
    rhomax = (2400.)/1000./Mw[j];   // 2400 kg m-3

    // find initial density guess
    DensityGuess( j, rhoguess );
    Pressure( rhoguess, pnew, pgrad, cfc );

    rhoold = rhoguess;
    rhonew = rhoguess;
    errnew = pnew - pmpa;
    errold = errnew;

    k = 0;
    do
    {
            k++;

            if ( k>= maxit )
            {
                    iErr = 2;
                    return iErr;
            }

            rhoold = rhonew;
            errold = errnew;
            Pressure( rhonew, pnew, pgrad, cfc );
            errnew = pnew - pmpa;
            step = - errnew/pgrad;
            rhonew = rhoold + step;

            if ( rhonew < rhomin )
                    rhonew = rhomin ;
            if ( rhonew > rhomax )
                    rhonew = rhomax;

    } while ( fabs(1.-pnew/pmpa) > tol );

    // calculate thermodynamic properties
    rho = rhonew;
    Helmholtz( j, rho, cfc );
    ar = RC*Tk*Phi[j];
    dard = RC*Tk*dPhiD[j];
    dardd = RC*Tk*dPhiDD[j];
    darddd = RC*Tk*dPhiDDD[j];
    dart = RC*(dPhiT[j]*Tk+Phi[j]);
    dartt = RC*(dPhiTT[j]*Tk+2.*dPhiT[j]);
    dardt = RC*(dPhiDT[j]*Tk+dPhiD[j]);
    darddt = RC*(dPhiDDT[j]*Tk+dPhiDD[j]);
    dardtt = RC*(dPhiDTT[j]*Tk+2.*dPhiDT[j]);

    vol = (1./rho)/10.;  // J bar-1
    lnfug = ( log(rho) + Phi[j] + pmpa/(rho*RC*Tk) ) + log(RC*Tk) - 1.;
    fug = exp(lnfug);
    fugc = fug/pmpa;
    pig = rho*RC*Tk;
    g = log(fugc)*RC*Tk;
    s = - dart + RC*log(pmpa/pig);
    h = g + Tk*s;
    a = ar - RC*Tk*log(pmpa/pig);
    cv = - Tk*dartt - RC;
    cp = cv + RC*pow((1.+(rho/RC)*dardt),2.) / (1.+(2.*rho)/(RC*Tk)*dard + pow(rho,2.)/(RC*Tk)*dardd);
    dpr = RC*Tk + 2.*rho*dard + pow(rho,2.)*dardd;
    dprr = 2.*dard + 4.*rho*dardd + pow(rho,2.)*darddd;
    dpt = rho*RC + pow(rho,2.)*dardt;
    dptt = pow(rho,2.)*dardtt;
    dprt = RC + 2.*rho*dardt + pow(rho,2.)*darddt;
    drp = (1./dpr);
    drpp = (-1.)*dprr*pow(dpr,-3.);
    drt = (-1.)*(1./dpr)*dpt;
    drtt = (dpt - dptt + dpt*(dprt-dpr)*pow(dpr,-1.))*pow(dpr,-1.)
            + (dprt - dpt*dprr*pow(dpr,-1.))*pow(dpr,-2.)*dpt;
    drtp = - pow(dpr,-1.) - (dprt-dpr)*pow(dpr,-2.) + dprr*dpt*pow(dpr,-3.);

    // copy results
    Fugpure[j][0] = fugc;  // fugacity coef.
    Fugpure[j][1] = g;     // Gres
    Fugpure[j][2] = h;     // Hres
    Fugpure[j][3] = s;     // Sres
    Fugpure[j][4] = vol;   // V
    Fugpure[j][5] = cp;    // CPres
    Fugpure[j][6] = cv;    // CVres
    Rho[j][0] = rho;      // mol cm-3
    Rho[j][1] = drt;      // mol cm-3 K-1
    Rho[j][2] = drtt;     // mol cm-3 K-2
    Rho[j][3] = drp;      // mol cm-3 MPa-1
    Rho[j][4] = drpp;     // mol cm-3 MPa-2
    Rho[j][5] = drtp;     // mol cm-3 K-1 MPa-1
    Rho[j][6] = dpr;      // MPa cm3 mol-1
    Rho[j][7] = dprr;     // MPa cm6 mol-2
    Rho[j][8] = dpt;      // MPa K-1
    Rho[j][9] = dptt;     // MPa K-2
    Rho[j][10] = dprt;    // MPa cm3 mol-1 K-1

    return iErr;
}



/// calculates fugacity and state functions of fluids other than H2O and CO2
/// adapted from CORK fluid model for consistency with Thermocalc
long int TSTPcalc::FugacityCorresponding( long int j )
{
    double a0, a1, a, b0, b, c0, c1, c, d0, d1, d, tcr, pcr, da, dc, dd, vc, fc, rtlnf, rho;
    double grs, hrs, srs, cprs, drt, drtt, drp, drpp, drtp,
            dvt, dvtt, dvp, dvpp, dvtp;

    a0 = 5.45963e-5;
    a1 = -8.6392e-6;
    b0 = 9.18301e-4;
    c0 = -3.30558e-5;
    c1 = 2.30524e-6;
    d0 = 6.93054e-7;
    d1 = -8.38293e-8;
    tcr = Tc[j];
    // pcr = Eosparm[j][1]/1000.;  // kbar
    pcr = Pc[j]*(10.)/(1000.);  // kbar (convert from MPa)
    a = a0*pow(tcr,2.)*sqrt(tcr)/pcr + a1*tcr*sqrt(tcr)/pcr*Tk;
    b = b0*tcr/pcr;
    c = c0*tcr/(pcr*sqrt(pcr)) + c1/(pcr*sqrt(pcr))*Tk;
    d = d0*tcr/pow(pcr,2.) + d1/pow(pcr,2.)*Tk;
    da = a1*pow(tcr,1.5)/pcr;
    dc = c1/pow(pcr,1.5);
    dd = d1/pow(pcr,2.);

    // molar volume and fugacity coefficient
    vc = RR*Tk/Pkb + b - a*RR*sqrt(Tk)/((RR*Tk+b*Pkb)*(RR*Tk+2.*b*Pkb)) + c*sqrt(Pkb) + d*Pkb;
    rtlnf = RR*Tk*log(1000*Pkb) + b*Pkb + a/(b*sqrt(Tk))*(log(RR*Tk+b*Pkb) - log(RR*Tk+2*b*Pkb))
            + (2./3.)*c*Pkb*sqrt(Pkb) + (d/2.)*pow(Pkb,2.);
    fc = exp(rtlnf/(RR*Tk))/(1000*Pkb);

    // volume and density derivatives
    dvt = ( RR/Pkb - (da*RR*pow(Tk,0.5))/((RR*Tk+b*Pkb)*(RR*Tk+2.*b*Pkb))
            - (a*RR)/(2.*pow(Tk,0.5)*(RR*Tk+b*Pkb)*(RR*Tk+2.*b*Pkb))
            + (a*pow(RR,2.)*pow(Tk,0.5))/(pow((RR*Tk+b*Pkb),2.)*(RR*Tk+2.*b*Pkb))
            + (a*pow(RR,2.)*pow(Tk,0.5))/((RR*Tk+b*Pkb)*pow((RR*Tk+2.*b*Pkb),2.))
            + dc*pow(Pkb,0.5) + dd*Pkb );
    dvtt = ( - (da*RR)/(pow(Tk,0.5)*(RR*Tk+b*Pkb)*(RR*Tk+2.*b*Pkb))
            + (2.*da*pow(RR,2.)*pow(Tk,0.5))/(pow((RR*Tk+b*Pkb),2.)*(RR*Tk+2.*b*Pkb))
            + (2.*da*pow(RR,2.)*pow(Tk,0.5))/((RR*Tk+b*Pkb)*pow((RR*Tk+2.*b*Pkb),2.))
            + (a*RR)/(4.*pow(Tk,1.5)*(RR*Tk+b*Pkb)*(RR*Tk+2.*b*Pkb))
            + (a*pow(RR,2.))/(pow(Tk,0.5)*pow((RR*Tk+b*Pkb),2.)*(RR*Tk+2.*b*Pkb))
            + (a*pow(RR,2.))/(pow(Tk,0.5)*(RR*Tk+b*Pkb)*pow((RR*Tk+2.*b*Pkb),2.))
            - (2.*a*pow(RR,3.)*pow(Tk,0.5))/(pow((RR*Tk+b*Pkb),3.)*(RR*Tk+2.*b*Pkb))
            - (2.*a*pow(RR,3.)*pow(Tk,0.5))/(pow((RR*Tk+b*Pkb),2.)*pow((RR*Tk+2.*b*Pkb),2.))
            - (2.*a*pow(RR,3.)*pow(Tk,0.5))/((RR*Tk+b*Pkb)*pow((RR*Tk+2.*b*Pkb),3.)) );
    dvp = ( - RR*Tk/pow(Pkb,2.) + (a*RR*pow(Tk,0.5)*b)/(pow((RR*Tk+b*Pkb),2.)*(RR*Tk+2.*b*Pkb))
            + (2.*a*RR*pow(Tk,0.5)*b)/((RR*Tk+b*Pkb)*pow((RR*Tk+2.*b*Pkb),2.))
            + c/(2.*pow(Pkb,0.5)) + d ) / (1000.);
    dvpp = ( (2.*RR*Tk)/pow(Pkb,3.) - (2.*a*RR*pow(Tk,0.5)*pow(b,2.))/(pow((RR*Tk+b*Pkb),3.)*(RR*Tk+2.*b*Pkb))
            - (4.*a*RR*pow(Tk,0.5)*pow(b,2.))/(pow((RR*Tk+b*Pkb),2.)*pow((RR*Tk+2.*b*Pkb),2.))
            - (8.*a*RR*pow(Tk,0.5)*pow(b,2.))/((RR*Tk+b*Pkb)*pow((RR*Tk+2.*b*Pkb),3.))
            - c/(4.*pow(Pkb,1.5)) ) / (1.0e6);
    dvtp = ( - RR/pow(Pkb,2.) + (da*RR*pow(Tk,0.5)*b)/(pow((RR*Tk+b*Pkb),2.)*(RR*Tk+2.*b*Pkb))
            + (a*RR*b)/(2.*pow(Tk,0.5)*pow((RR*Tk+b*Pkb),2.)*(RR*Tk+2.*b*Pkb))
            - (2.*a*pow(RR,2.)*pow(Tk,0.5)*b)/(pow((RR*Tk+b*Pkb),3.)*(RR*Tk+2.*b*Pkb))
            - (3.*a*pow(RR,2.)*pow(Tk,0.5)*b)/(pow((RR*Tk+b*Pkb),2.)*pow((RR*Tk+2.*b*Pkb),2.))
            + (2.*da*RR*pow(Tk,0.5)*b)/((RR*Tk+b*Pkb)*pow((RR*Tk+2.*b*Pkb),2.))
            + (a*RR*b)/(pow(Tk,0.5)*(RR*Tk+b*Pkb)*pow((RR*Tk+2.*b*Pkb),2.))
            - (4.*a*pow(RR,2.)*pow(Tk,0.5)*b)/((RR*Tk+b*Pkb)*pow((RR*Tk+2.*b*Pkb),3.))
            + dc/(2.*pow(Pkb,2.)) + dd ) / (1000.);
    rho = (1./vc);
    drt = - dvt*pow(rho,2.);
    drp = - dvp*pow(rho,2.);
    drtt = - (dvtt*pow(rho,3.)-2.*pow(drt,2.))/rho;
    drpp = - (dvpp*pow(rho,3.)-2.*pow(drp,2.))/rho;
    drtp = - (dvtp*pow(rho,3.)-2.*drt*drp)/rho;

    // thermodynamic properties
    grs = R_CONST*Tk*log(fc);
    hrs = ( (pow(Tk,0.5)*da*log(RR*Tk+2.*b*Pkb))/b - (pow(Tk,0.5)*da*log(RR*Tk+b*Pkb))/b
            - (3.*a*log(RR*Tk+2.*b*Pkb))/(2.*pow(Tk,0.5)*b) + (3.*a*log(RR*Tk+b*Pkb))/(2.*pow(Tk,0.5)*b)
            - (pow(Tk,0.5)*a*RR)/(b*(RR*Tk+b*Pkb)) + (pow(Tk,0.5)*a*RR)/(b*(RR*Tk+2.*b*Pkb))
            - (2.*Tk*dc*pow(Pkb,1.5))/3. - (Tk*dd*pow(Pkb,2.))/2. + b*Pkb
            + (2.*c*pow(Pkb,1.5))/3. + (d*pow(Pkb,2.))/2. )*1000.;
    srs = (hrs-grs)/Tk;
    cprs = ( - (2.*pow(Tk,0.5)*da*RR)/(b*(RR*Tk+b*Pkb)) + (pow(Tk,0.5)*a*pow(RR,2.))/(b*pow((RR*Tk+b*Pkb),2.))
            + (2.*pow(Tk,0.5)*da*RR)/(b*(RR*Tk+2.*b*Pkb)) - (a*RR)/(pow(Tk,0.5)*b*(RR*Tk+2.*b*Pkb))
            - (da*log(RR*Tk+2.*b*Pkb))/(pow(Tk,0.5)*b) + (da*log(RR*Tk+b*Pkb))/(pow(Tk,0.5)*b)
            + (3.*a*log(RR*Tk+2.*b*Pkb))/(4.*pow(Tk,1.5)*b) - (3.*a*log(RR*Tk+b*Pkb))/(4.*pow(Tk,1.5)*b)
            - (pow(Tk,0.5)*a*pow(RR,2.))/(b*pow((RR*Tk+2.*b*Pkb),2.)) + (a*RR)/(pow(Tk,0.5)*b*(RR*Tk+b*Pkb)) )*1000.;

    // copy results
    Fugpure[j][0] = fc;     // fugacity coefficient
    Fugpure[j][1] = grs;
    Fugpure[j][2] = hrs;
    Fugpure[j][3] = srs;
    Fugpure[j][4] = vc;
    Fugpure[j][5] = cprs;
    Rho[j][0] = rho*(0.1);   // mol cm-3
    Rho[j][1] = drt*(0.1);   // mol cm-3 K-1
    Rho[j][2] = drtt*(0.1);  // mol cm-3 K-2
    Rho[j][3] = drp;         // mol cm-3 MPa-1
    Rho[j][4] = drpp*(10.);  // mol cm-3 MPa-2
    Rho[j][5] = drtp;        // mol cm-3 K-1 MPa-1
    Rho[j][6] = 1.0;         // MPa cm3 mol-1
    Rho[j][7] = 1.0;         // MPa cm6 mol-2
    Rho[j][8] = 1.0;         // MPa K-1
    Rho[j][9] = 1.0;         // MPa K-2
    Rho[j][10] = 1.0;        // MPa cm3 mol-1 K-1

    return 0;
}



/// calculates and returns density guess for pure fluids
long int TSTPcalc::DensityGuess( long int j, double &Rhoguess )
{
    double tred, pred, rhomin, rhomax, rhoguess;

    tred = Tk/Tc[j];
    pred = Pmpa/Pc[j];

    if ( EosCode[j] == CEM_H2O_ )  // H2O
    {
            rhomin = (1.e-20/1000.)/Mw[j];  // 1.e-20 kg m-3
            rhomax = (1800./1000.)/Mw[j];   // 1800 kg m-3

            if ( Tk < Tc[j] )
            {
                    PsatH2O( j );
                    if ( Pmpa < Psat[j] )
                            rhoguess = Rhov[j];
                    else
                            rhoguess = Rhol[j];
            }

            else
            {
                    if ( Pmpa < 400. )
                            rhoguess = 1000./(Tk-273.15)*2.*Pmpa;
                    else
                            rhoguess = 800.;
            }

            Rhoguess = (rhoguess/1000.)/Mw[j];  // mol cm-3

            if ( Rhoguess < rhomin )
                    Rhoguess = rhomin;
            if ( Rhoguess > rhomax )
                    Rhoguess = rhomax;
    }

    else  // CO2 (other fluids)
    {
            rhomin = (1.e-20)/1000./Mw[j];  // 1.e-20 kg m-3
            rhomax = (2400.)/1000./Mw[j];   // 2400 kg m-3

            if ( Tk < Tc[j] )
            {
                    PsatCO2( j );
                    if ( Pmpa < Psat[j] )
                            rhoguess = Rhov[j];  // kg m-3
                    else
                            rhoguess = Rhol[j];
            }

            else
            {
                    rhoguess = (pred*Pc[j]*1000.)/((0.1889241)*tred*Tc[j]);  // kg m-3
            }

            Rhoguess = (rhoguess/1000.)/Mw[j];  // mol cm-3

            if ( Rhoguess < rhomin )
                    Rhoguess = rhomin;
            if ( Rhoguess > rhomax )
                    Rhoguess = rhomax;
    }

    return 0;
}



/// calculates pressure (P) and first density derivative (dP/dRho)
long int TSTPcalc::Pressure ( double rho, double &p, double &dpdrho, double **cf )
{
    long int k;
    double pred, dpred;
    double c[10];

    for (k=0; k<10; k++)
    {
        c[k] = cf[k][0]*pow(Tk,-4.) + cf[k][1]*pow(Tk,-2.) + cf[k][2]*pow(Tk,-1.) + cf[k][3]
            + cf[k][4]*Tk + cf[k][5]*pow(Tk,2.);
    }

    pred = rho + c[0]*pow(rho,2.) - pow(rho,2.) * ( (c[2]+2.*c[3]*rho+3.*c[4]*pow(rho,2.)+4.*c[5]*pow(rho,3.))
            / pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),2.) )
            + c[6]*pow(rho,2.)*exp(-c[7]*rho) + c[8]*pow(rho,2.)*exp(-c[9]*rho);
    dpred = 1. + 2.*c[0]*rho - 2.*rho*(c[2]+2.*c[3]*rho+3.*c[4]*pow(rho,2.)+4.*c[5]*pow(rho,3.))
            / pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),2.)
            - pow(rho,2.)*(2.*c[3]+6.*c[4]*rho+12.*c[5]*pow(rho,2.))
            / pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),2.)
            + 2.*pow(rho,2.)*pow((c[2]+2.*c[3]*rho+3.*c[4]*pow(rho,2.)+4.*c[5]*pow(rho,3.)),2.)
            / pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),3.)
            + 2.*c[6]*rho*exp(-c[7]*rho) - c[6]*pow(rho,2.)*c[7]*exp(-c[7]*rho)
            + 2.*c[8]*rho*exp(-c[9]*rho) - c[8]*pow(rho,2.)*c[9]*exp(-c[9]*rho);

    p = pred*(RC*Tk);
    dpdrho = dpred*(RC*Tk);

    return 0;
}



/// calculates reduced Helmholtz energy and derivatives
long int TSTPcalc::Helmholtz( long int j, double rho, double **cf )
{
    long int k;
    double c[10], dc[10], d2c[10];

    for (k=0; k<10; k++)
    {
            c[k] = cf[k][0]*pow(Tk,-4.) + cf[k][1]*pow(Tk,-2.) + cf[k][2]*pow(Tk,-1.) + cf[k][3]
                    + cf[k][4]*Tk + cf[k][5]*pow(Tk,2.);
            dc[k] = - 4.*cf[k][0]*pow(Tk,-5.) - 2.*cf[k][1]*pow(Tk,-3.) - cf[k][2]*pow(Tk,-2.)
                    + cf[k][4] + 2.*cf[k][5]*Tk;
            d2c[k] = 20.*cf[k][0]*pow(Tk,-6.) + 6.*cf[k][1]*pow(Tk,-4.) + 2.*cf[k][2]*pow(Tk,-3.)
                    + 2.*cf[k][5];
    }

    Phi[j] = c[0]*rho + ( 1./(c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)) - 1./c[1] )
            - (c[6]/c[7])*(exp(-c[7]*rho)-1.) - (c[8]/c[9])*(exp(-c[9]*rho)-1.);
    dPhiD[j] = c[0] - 1./pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),2.)
            * (c[2]+2.*c[3]*rho+3.*c[4]*pow(rho,2.)+4.*c[5]*pow(rho,3.))
            + c[6]*exp(-c[7]*rho) + c[8]*exp(-c[9]*rho);
    dPhiDD[j] = 2./pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),3.)
            * pow((c[2]+2.*c[3]*rho+3.*c[4]*pow(rho,2.)+4.*c[5]*pow(rho,3.)),2.)
            - 1./pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),2.)
            * (2.*c[3]+6.*c[4]*rho+12.*c[5]*pow(rho,2.))
            - c[6]*c[7]*exp(-c[7]*rho) - c[8]*c[9]*exp(-c[9]*rho);
    dPhiT[j] = dc[0]*rho - 1./pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),2.)
            * (dc[1]+dc[2]*rho+dc[3]*pow(rho,2.)+dc[4]*pow(rho,3.)+dc[5]*pow(rho,4.)) + 1./pow(c[1],2.)*dc[1]
            - dc[6]/c[7]*(exp(-c[7]*rho)-1.) + c[6]/pow(c[7],2.)*(exp(-c[7]*rho)-1.)*dc[7]
            + c[6]/c[7]*dc[7]*rho*exp(-c[7]*rho) - dc[8]/c[9]*(exp(-c[9]*rho)-1)
            + c[8]/pow(c[9],2.)*(exp(-c[9]*rho)-1.)*dc[9] + c[8]/c[9]*dc[9]*rho*exp(-c[9]*rho);
    dPhiTT[j] = - 2./pow(c[1],3.)*pow(dc[1],2.)
            + 2./pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),3.)
            * pow((dc[1]+dc[2]*rho+dc[3]*pow(rho,2.)+dc[4]*pow(rho,3.)+dc[5]*pow(rho,4.)),2.)
            - d2c[8]/c[9]*(exp(-c[9]*rho)-1.) - d2c[6]/c[7]*(exp(-c[7]*rho)-1.)
            + 2.*dc[6]/pow(c[7],2.)*(exp(-c[7]*rho)-1)*dc[7] + d2c[0]*rho
            - 1./pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),2.)
            * (d2c[1]+d2c[2]*rho+d2c[3]*pow(rho,2.)+d2c[4]*pow(rho,3.)+d2c[5]*pow(rho,4.))
            - 2.*c[8]/pow(c[9],2.)*pow(dc[9],2.)*rho*exp(-c[9]*rho) + c[8]/pow(c[9],2.)*(exp(-c[9]*rho)-1.)*d2c[9]
            - 2.*c[8]/pow(c[9],3.)*(exp(-c[9]*rho)-1.)*pow(dc[9],2.) + 2.*dc[8]/c[9]*dc[9]*rho*exp(-c[9]*rho)
            + 2.*dc[8]/pow(c[9],2.)*(exp(-c[9]*rho)-1.)*dc[9] - c[6]/c[7]*pow(dc[7],2.)*pow(rho,2.)*exp(-c[7]*rho)
            + c[6]/c[7]*d2c[7]*rho*exp(-c[7]*rho) - 2.*c[6]/pow(c[7],2.)*pow(dc[7],2.)*rho*exp(-c[7]*rho)
            + c[6]/pow(c[7],2.)*(exp(-c[7]*rho)-1.)*d2c[7] - 2.*c[6]/pow(c[7],3.)*(exp(-c[7]*rho)-1.)*pow(dc[7],2.)
            + 2.*dc[6]/c[7]*dc[7]*rho*exp(-c[7]*rho) + 1./pow(c[1],2.)*d2c[1]
            - c[8]/c[9]*pow(dc[9],2.)*pow(rho,2.)*exp(-c[9]*rho) + c[8]/c[9]*d2c[9]*rho*exp(-c[9]*rho);
    dPhiDT[j] = dc[0] + 2./pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),3.)
            * (dc[1]+dc[2]*rho+dc[3]*pow(rho,2.)+dc[4]*pow(rho,3.)+dc[5]*pow(rho,4.))
            * (c[2]+2.*c[3]*rho+3.*c[4]*pow(rho,2.)+4.*c[5]*pow(rho,3.))
            - 1./pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),2.)
            * (dc[2]+2.*dc[3]*rho+3.*dc[4]*pow(rho,2.)+4.*dc[5]*pow(rho,3.))
            + dc[6]*exp(-c[7]*rho) - c[6]*dc[7]*rho*exp(-c[7]*rho)
            + dc[8]*exp(-c[9]*rho) - c[8]*dc[9]*rho*exp(-c[9]*rho);
    dPhiDDD[j] = - 6./pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),4.)
            * pow((c[2]+2.*c[3]*rho+3.*c[4]*pow(rho,2.)+4.*c[5]*pow(rho,3.)),3.)
            + 6./pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),3.)
            * (c[2]+2.*c[3]*rho+3.*c[4]*pow(rho,2.)+4.*c[5]*pow(rho,3.))*(2.*c[3]+6.*c[4]*rho+12.*c[5]*pow(rho,2.))
            - 1./pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),2.) * (6.*c[4]+24.*c[5]*rho)
            + c[6]*pow(c[7],2.)*exp(-c[7]*rho)+c[8]*pow(c[9],2.)*exp(-c[9]*rho);
    dPhiDDT[j] = - 6./pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),4.)
            * (dc[1]+dc[2]*rho+dc[3]*pow(rho,2.)+dc[4]*pow(rho,3.)+dc[5]*pow(rho,4.))
            * pow((c[2]+2.*c[3]*rho+3.*c[4]*pow(rho,2.)+4.*c[5]*pow(rho,3.)),2.)
            + 4./pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),3.)
            * (dc[2]+2.*dc[3]*rho+3.*dc[4]*pow(rho,2.)+4.*dc[5]*pow(rho,3.))
            * (c[2]+2.*c[3]*rho+3.*c[4]*pow(rho,2.)+4.*c[5]*pow(rho,3.))
            + 2./pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),3.)
            * (dc[1]+dc[2]*rho+dc[3]*pow(rho,2.)+dc[4]*pow(rho,3.)+dc[5]*pow(rho,4.))
            * (2.*c[3]+6.*c[4]*rho+12.*c[5]*pow(rho,2.))
            - 1./pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),2.)
            * (2.*dc[3]+6.*dc[4]*rho+12.*dc[5]*pow(rho,2.))
            - dc[6]*c[7]*exp(-c[7]*rho) - c[6]*dc[7]*exp(-c[7]*rho)
            + c[6]*dc[7]*rho*c[7]*exp(-c[7]*rho) - dc[8]*c[9]*exp(-c[9]*rho)
            - c[8]*dc[9]*exp(-c[9]*rho) + c[8]*dc[9]*rho*c[9]*exp(-c[9]*rho);
    dPhiDTT[j] = 2./pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),3.)
            * (d2c[1]+d2c[2]*rho+d2c[3]*pow(rho,2.)+d2c[4]*pow(rho,3.)+d2c[5]*pow(rho,4.))
            * (c[2]+2.*c[3]*rho+3.*c[4]*pow(rho,2.)+4.*c[5]*pow(rho,3.)) + d2c[0]
            - 2.*dc[8]*dc[9]*rho*exp(-c[9]*rho) - c[6]*d2c[7]*rho*exp(-c[7]*rho)
            + c[6]*pow(dc[7],2.)*pow(rho,2.)*exp(-c[7]*rho)
            - 6./pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),4.)
            * pow((dc[1]+dc[2]*rho+dc[3]*pow(rho,2.)+dc[4]*pow(rho,3.)+dc[5]*pow(rho,4.)),2.)
            * (c[2]+2.*c[3]*rho+3.*c[4]*pow(rho,2.)+4.*c[5]*pow(rho,3.))
            + 4./pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),3.)
            * (dc[1]+dc[2]*rho+dc[3]*pow(rho,2.)+dc[4]*pow(rho,3.)+dc[5]*pow(rho,4.))
            * (dc[2]+2.*dc[3]*rho+3.*dc[4]*pow(rho,2.)+4.*dc[5]*pow(rho,3.)) + d2c[6]*exp(-c[7]*rho)
            - 1./pow((c[1]+c[2]*rho+c[3]*pow(rho,2.)+c[4]*pow(rho,3.)+c[5]*pow(rho,4.)),2.)
            * (d2c[2]+2.*d2c[3]*rho+3.*d2c[4]*pow(rho,2.)+4.*d2c[5]*pow(rho,3.))
            - 2.*dc[6]*dc[7]*rho*exp(-c[7]*rho) + c[8]*pow(dc[9],2.)*pow(rho,2.)*exp(-c[9]*rho)
            - c[8]*d2c[9]*rho*exp(-c[9]*rho) + d2c[8]*exp(-c[9]*rho);

    return 0;
}



/// calculates residual state functions of the fluid mixture
long int TSTPcalc::ResidualFunct()
{
    long int i, j;
    double sumphi, vi, vj, dvi, dvj, d2vi, d2vj, dvip, dvjp, Y, dY, d2Y, dYp,
            rhoi, drti, drtti, drpi, rhoj, drtj, drttj, drpj;
    double gex, dgex, d2gex, dgexp, sex, hex, cpex, vex, gph, sph, hph, cpph, vph;
    gex = 0.; dgex = 0.; d2gex = 0.; dgexp = 0.;
    sex = 0.; hex = 0.; vex = 0.; cpex = 0.;
    gph = 0.; sph = 0.; hph = 0.; cpph = 0.; vph = 0.;
    sumphi = 0.; Y = 0.; dY = 0.; d2Y = 0.; dYp = 0.;

    // phi values (and derivatives)
    for (i=0; i<NComp; i++)
    {
            rhoi = Rho[i][0];
            drti = Rho[i][1];
            drtti = Rho[i][2];
            drpi = Rho[i][3];
            vi = 1./rhoi;
            dvi = - pow(rhoi,-2.)*drti;
            d2vi = 2.*pow(rhoi,-3.)*pow(drti,2.) - pow(rhoi,-2.)*drtti;
            dvip = - pow(rhoi,-2.)*drpi;
            sumphi += x[i]*vi;
            Y += x[i]*vi;
            dY += x[i]*dvi;
            d2Y += x[i]*d2vi;
            dYp += x[i]*dvip;
    }

    for (i=0; i<NComp; i++)
    {
            rhoi = Rho[i][0];
            drti = Rho[i][1];
            drtti = Rho[i][2];
            drpi = Rho[i][3];
            vi = 1./rhoi;
            dvi = - pow(rhoi,-2.)*drti;
            d2vi = 2.*pow(rhoi,-3.)*pow(drti,2.) - pow(rhoi,-2.)*drtti;
            dvip = - pow(rhoi,-2.)*drpi;
            phi[i] = x[i]*vi/sumphi;
            dphi[i] = x[i]*(dvi*Y-vi*dY)/pow(Y,2.);
            dphip[i] = x[i]*(dvip*Y-vi*dYp)/pow(Y,2.);
            d2phi[i] = x[i]*(d2vi*Y+dvi*dY)/pow(Y,2.) - x[i]*(dvi*Y)*(2.*dY)/pow(Y,3.)
                    - x[i]*(dvi*dY+vi*d2Y)/pow(Y,2.) + x[i]*(vi*dY)*(2.*dY)/pow(Y,3.);
    }

    // interaction parameters
    for (i=0; i<NComp; i++)
    {
            for (j=i+1; j<NComp; j++)
            {
                    rhoi = Rho[i][0];
                    drti = Rho[i][1];
                    drtti = Rho[i][2];
                    drpi = Rho[i][3];
                    rhoj = Rho[j][0];
                    drtj = Rho[j][1];
                    drttj = Rho[j][2];
                    drpj = Rho[j][3];

                    vi = 1./rhoi;
                    dvi = - pow(rhoi,-2.)*drti;
                    d2vi = 2.*pow(rhoi,-3.)*pow(drti,2.) - pow(rhoi,-2.)*drtti;
                    dvip = - pow(rhoi,-2.)*drpi;
                    vj = 1./rhoj;
                    dvj = - pow(rhoj,-2.)*drtj;
                    d2vj = 2.*pow(rhoj,-3.)*pow(drtj,2.) - pow(rhoj,-2.)*drttj;
                    dvjp = - pow(rhoj,-2.)*drpj;

                    B[i][j] = A[i][j]*2*Y/(vi*vj);
                    dB[i][j] = A[i][j]*(2.*dY*vi*vj-2.*Y*(dvi*vj+vi*dvj))/(pow(vi,2.)*pow(vj,2.));
                    d2B[i][j] = A[i][j]*((2.*d2Y)*(vi*vj) + (2.*dY)*(dvi*vj+vi*dvj))/pow((vi*vj),2.)
                            - A[i][j]*((2.*dY)*(vi*vj))*(2.*(dvi*vj+vi*dvj))/pow((vi*vj),3.)
                            - A[i][j]*((2.*dY)*(dvi*vj+vi*dvj) + (2.*Y)*(d2vi*vj+2.*dvi*dvj+vi*d2vj))/pow((vi*vj),2.)
                            + A[i][j]*((2.*Y)*(dvi*vj+vi*dvj))*(2.*(dvi*vj+vi*dvj))/pow((vi*vj),3.);
                    dBp[i][j] = A[i][j]*(2.*dYp*vi*vj-2.*Y*(dvip*vj+vi*dvjp))/(pow(vi,2.)*pow(vj,2.));
            }
    }

    // excess properties
    for (i=0; i<NComp; i++)
    {
            for (j=i+1; j<NComp; j++)
            {
                    gex += phi[i]*phi[j]*B[i][j];
                    dgex += dphi[i]*phi[j]*B[i][j] + phi[i]*dphi[j]*B[i][j] + phi[i]*phi[j]*dB[i][j];
                    d2gex +=d2phi[i]*phi[j]*B[i][j] + 2.*dphi[i]*dphi[j]*B[i][j] + 2.*dphi[i]*phi[j]*dB[i][j]
                                    + phi[i]*d2phi[j]*B[i][j] + 2.*phi[i]*dphi[j]*dB[i][j] + phi[i]*phi[j]*d2B[i][j];
                    dgexp += dphip[i]*phi[j]*B[i][j] + phi[i]*dphip[j]*B[i][j] + phi[i]*phi[j]*dBp[i][j];
            }
    }

    // thermodynamic properties
    sex = - dgex;
    hex = gex + Tk*sex;
    cpex = - d2gex*Tk;
    vex = dgexp;
    Gex = gex*(10.);
    Sex = sex*(10.);
    Hex = hex*(10.);
    CPex = cpex*(10.);
    Vex = vex;

    // incrementing residual properties
    for (j=0; j<NComp; j++)
    {
        gph += x[j]*Fugpure[j][1];
        hph += x[j]*Fugpure[j][2];
        sph += x[j]*Fugpure[j][3];
        vph += x[j]*Fugpure[j][4];
        cpph += x[j]*Fugpure[j][5];
    }
    Grs = gph + Gex;
    Hrs = hph + Hex;
    Srs = sph + Sex;
    Vrs = vph + Vex;
    CPrs = cpph + CPex;

    return 0;
}



/// calculates saturation pressure of H2O (required only for initial density guess)
long int TSTPcalc::PsatH2O( long int j )
{
    double Rhoc, a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, b5, b6, c1, c2, c3, c4, c5, c6,
                    tau, psat, rhol, rhov, ppc, rrc;

    Rhoc = 322.0;   // kg m-3
    a1 = -7.85951783;
    a2 = 1.84408259;
    a3 = -11.7866497;
    a4 = 22.6807411;
    a5 = -15.9618719;
    a6 = 1.80122502;
    b1 = 1.99274064;
    b2 = 1.09965342;
    b3 = -0.510839303;
    b4 = -1.75493479;
    b5 = -45.5170352;
    b6 = -6.74694450e5;
    c1 = -2.03150240;
    c2 = -2.68302940;
    c3 = -5.38626492;
    c4 = -17.2991605;
    c5 = -44.7586581;
    c6 = -63.9201063;

    tau = 1 - Tk/Tc[j];
    ppc = (Tc[j]/Tk)*( a1*tau + a2*pow(tau,1.5) + a3*pow(tau,3.) + a4*pow(tau,3.5)
            + a5*pow(tau,4.) + a6*pow(tau,7.5) );
    psat = exp(ppc) * Pc[j];
    rhol = ( 1. + b1*pow(tau,(1./3.)) + b2*pow(tau,(2./3.)) + b3*pow(tau,(5./3.)) + b4*pow(tau,(16./3.))
            + b5*pow(tau,(43./3.)) + b6*pow(tau,(110./3.)) ) * Rhoc;
    rrc = c1*pow(tau,2./6.) + c2*pow(tau,4./6.) + c3*pow(tau,8./6.) + c4*pow(tau,18./6.)
            + c5*pow(tau,37./6.) + c6*pow(tau,71./6.);
    rhov = exp(rrc) * Rhoc;

    // copy results
    Psat[j] = psat;  // MPa
    Rhol[j] = rhol;
    Rhov[j] = rhov;

    return 0;
}



/// calculates saturation pressure of CO2 (required only for initial density guess)
long int TSTPcalc::PsatCO2( long int j )
{
    double Rhoc, a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, c5, tau, ppc, psat,
                    rhol, rhov, rholc, rhovc;

    Rhoc = 467.6;   // kg m-3
    a1 = -7.0602087;
    a2 = 1.9391218;
    a3 = -1.6463597;
    a4 = -3.2995634;
    b1 = 1.9245108;
    b2 = -0.62385555;
    b3 = -0.32731127;
    b4 = 0.39245142;
    c1 = -1.7074879;
    c2 = -0.82274670;
    c3 = -4.6008549;
    c4 = -10.111178;
    c5 = -29.742252;

    tau = 1. - Tk/Tc[j];
    ppc = (Tc[j]/Tk)*( a1*pow(tau,1.) + a2*pow(tau,1.5) + a3*pow(tau,2.) + a4*pow(tau,4.) );
    psat = exp(ppc)*Pc[j];
    rholc = b1*pow(tau,0.34)+ b2*pow(tau,0.5) + b3*pow(tau,(10./6.)) + b4*pow(tau,(11./6.));
    rhol = exp(rholc)*Rhoc;
    rhovc = c1*pow(tau,0.34) + c2*pow(tau,0.5) + c3*pow(tau,1.) + c4*pow(tau,(7./3.)) + c5*pow(tau,(14./3.));
    rhov = exp(rhovc)*Rhoc;

    // copy results
    Psat[j] = psat;  // MPa
    Rhol[j] = rhol;  // kg m-3
    Rhov[j] = rhov;  // kg m-3

    return 0;
}



#ifndef IPMGEMPLUGIN

/// Calculates pure species properties (called from DCthermo)
long int TSTPcalc::STPCalcFugPure( double Tmin, float *Cpg, double *FugProps )
{
    long int iErr = 0;
     double Coeff[7];

    for( int ii=0; ii<7; ii++ )
        Coeff[ii] = (double)Cpg[ii];

    if( (Tk >= Tmin) && (Tk < 1e4) && (Pbar >= 1e-5) && (Pbar < 1e5) )
    {
        iErr = FugacityPT( 0, Coeff );
        for( int i=0; i<6; i++ )
            FugProps[i] = Fugpure[0][i];
        return iErr;
    }

    else
    {
        for( int i=1; i<6; i++ )
            FugProps[i] = 0.;
        FugProps[0] = 1.;
        FugProps[4] = 8.31451*Tk/Pbar;
        return -1;
    }

    return 0;
}

#endif





//--------------------- End of s_solmod2.cpp ---------------------------


