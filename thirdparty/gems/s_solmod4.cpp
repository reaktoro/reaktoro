//-------------------------------------------------------------------
// $Id: s_fgl1.cpp 724 2012-10-02 14:25:25Z kulik $
//
/// \file s_solmod4.cpp
/// Implementation of TSolMod derived classes for ion-association aqueous
/// activity models (THelgeson, TDavies, TLimitingLaw, TDebyeHueckel,
/// TKarpov, TShvarov)
//
// Copyright (c) 2008-2012  T.Wagner, S.Dmitrieva, D.Kulik
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
//

#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include "verror.h"
#include "s_solmod.h"





//=============================================================================================
// Extended Debye-Hueckel (EDH) model for aqueous solutions, Helgeson version
// References: Helgeson et al. (1981); Oelkers and Helgeson (1990);
// Pokrovskii and Helgeson (1995; 1997a; 1997b)
// (c) TW May 2009
//=============================================================================================


// Generic constructor for the THelgeson class
THelgeson::THelgeson( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW ):
                TSolMod( sd )
{
    alloc_internal();
    m = arM;
    z = arZ;
    RhoW = dW;
    EpsW = eW;
    ac = aIPc[1];   // common ion size parameter
    bc = aIPc[0];   // common b_gamma
    flagNeut = (long int)aIPc[2];   // 0: unity, 1: calculated
    flagH2O = (long int)aIPc[3];   // 0: unity, 1: calculated
    flagElect = (long int)aIPc[4];  // 0: constant, 1: NaCl, 2: KCl, 3: NaOH, 4: KOH
}


THelgeson::~THelgeson()
{
	free_internal();
}


void THelgeson::alloc_internal()
{
	LnG = new double [NComp];
	dLnGdT = new double [NComp];
	d2LnGdT2 = new double [NComp];
	dLnGdP = new double [NComp];
}


void THelgeson::free_internal()
{
	delete[]LnG;
	delete[]dLnGdT;
	delete[]d2LnGdT2;
	delete[]dLnGdP;
}


/// Calculates T,P corrected parameters
long int THelgeson::PTparam()
{
	double alp, bet, dal, rho, eps, dedt, d2edt2, dedp;

	// read and convert rho and eps
	rho = RhoW[0];
	alp = - 1./rho*RhoW[1];
	dal = pow(alp,2.) - 1./rho*RhoW[2];
	bet = 1./rho*RhoW[3];
	eps = EpsW[0];
	dedt = 1./eps*EpsW[1];
	d2edt2 = - 1./pow(eps,2.)*pow(EpsW[1],2.) + 1./eps*EpsW[2];
	dedp = 1./eps*EpsW[3];

	// calculate A and B terms of Debye-Huckel equation (and derivatives)
	A = (1.82483e6)*sqrt(rho) / pow(Tk*eps,1.5);
	dAdT = - 3./2.*A*( dedt + 1./Tk + alp/3. );
	d2AdT2 = 1./A*pow(dAdT,2.) - 3./2.*A*( d2edt2 - 1/pow(Tk,2.) + 1/3.*dal );
	dAdP = 1./2.*A*( bet - 3.*dedp);
	B = (50.2916)*sqrt(rho) / sqrt(Tk*eps);
	dBdT = - 1./2.*B*( dedt +1./Tk + alp );
	d2BdT2 = 1./B*pow(dBdT,2.) - 1./2.*B*( d2edt2 - 1./pow(Tk,2.) + dal );
	dBdP = 1./2.*B*( bet - dedp );

 	ErrorIf( fabs(A) < 1e-9 || fabs(B) < 1e-9, "Helgeson EDH model",
 			"Error: DH parameter A or B was not calculated - no values of RhoW and EpsW !" );

	// b_gamma constant
	if ( flagElect == 0)
	{
		bgam = bc;
		dbgdT = 0.0;
		d2bgdT2 = 0.0;
		dbgdP = 0.0;
		ao = ac;
		daodT = 0.0;
		d2aodT2 = 0.0;
		daodP = 0.0;
	}

	// b_gamma TP-dependent
	else
	{
		Gfunction();
		BgammaTP();
		IonsizeTP();
		// ao = ac;
		// daodT = 0.0;
		// d2aodT2 = 0.0;
		// daodP = 0.0;
	}

	return 0;
}


/// Calculates activity coefficients
long int THelgeson::MixMod()
{
	long int j, k, w;
	double sqI, Z2, lgGam, lnGam, Nw, Lgam, lnxw, xw, Lam, sig,
			Phi, Phit, zc, za, psi, lnaw, lg_to_ln;
	zc = 1.; za = 1.; psi = 1.; lg_to_ln = 2.302585093;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molaities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		lgGam = 0.0;
		lnGam = 0.0;

		// charged species
		if ( z[j] )
		{
			lgGam = 0.0;
			Z2 = z[j]*z[j];
			lgGam = ( - A * sqI * Z2 ) / ( 1. + B * ao * sqI ) + bgam * IS ;
			lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				lgGam = 0.0;
				if ( flagNeut == 1 )
					lgGam = bgam * IS;
				else
					lgGam = 0.0;
				lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
				continue;
			}

			// water solvent
			else
			{
				lgGam = 0.0;
				lnGam = 0.0;
				if ( flagH2O == 1 )
				{
					Phit = 0.0;
					// Phi corrected using eq. (190) from Helgeson et al. (1981)
					Lam = 1. + (ao*B) * sqI;
					sig = 3./(pow(ao,3.)*pow(B,3.)*pow(IS,(3./2.))) * (Lam-1./Lam-2*log(Lam));
					// Phi = -2.3025851*(A*sqI*SigTerm/3. + Lgam/(0.0180153*2.*IS) - bgam*IS/2.);
					// Phi = - log(10.) * (molZ/molT) * ( (zc*za*A*sqI*SigTerm)/3. + (Lgam*psi)/(0.0180153*2.*IS) - bgam*IS/2. );

					for (k=0; k<(NComp-1); k++)
					{
						if ( (z[k] == 0) && (flagNeut == 0) )
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) );
						else
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) - bgam*IS/2. );
					}

					Phi = Phit/molT;
					lnaw = - Phi*molT/Nw;
					lnGam = lnaw - lnxw;
				}
				else
					lnGam = 0.0;
				lnGamma[j] = lnGam;
			}
		}
	}  // j

	return 0;
}


/// calculates excess properties
long int THelgeson::ExcessProp( double *Zex )
{
	long int j, k, w;
	double sqI, Z2, Nw, Lgam, lnxw, xw, Lam, sig, Phi, dPhidT, d2PhidT2, dPhidP,
			Phit, dPhitdT, d2PhitdT2, dPhitdP, lnaw, lg_to_ln, g, dgt, d2gt, dgp;
	double U, V, dUdT, dVdT, d2UdT2, d2VdT2, dUdP, dVdP, U1, U2, U3, V1, V2, V3,
			dU1dT, dU2dT, dU3dT, dV1dT, dV2dT, dV3dT, d2U1dT2, d2U2dT2, d2U3dT2,
			d2V1dT2, d2V2dT2, d2V3dT2, dU1dP, dU2dP, dU3dP, dV1dP, dV2dP, dV3dP,
			L, dLdT, d2LdT2, dLdP, Z, dZdT, d2ZdT2, dZdP;
	lg_to_ln = 2.302585093;
	g = 0.; dgt = 0.; d2gt = 0.; dgp = 0.;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		// charged species
		if ( z[j] )
		{
			Z2 = z[j]*z[j];
			U = - (Z2*A) * sqI;
			dUdT = - (Z2*dAdT) * sqI;
			d2UdT2 = - (Z2*d2AdT2) * sqI;
			dUdP = - (Z2*dAdP) * sqI;
			V = 1. + (ao*B) * sqI;
			dVdT = ( daodT*B + ao*dBdT ) * sqI;
			d2VdT2 = ( d2aodT2*B + 2.*daodT*dBdT + ao*d2BdT2 ) * sqI;
			dVdP = ( daodP*B + ao*dBdP ) * sqI;
			LnG[j] = ( U/V + bgam*IS ) * lg_to_ln;
			dLnGdT[j] = ( (dUdT*V - U*dVdT)/pow(V,2.) + dbgdT*IS ) * lg_to_ln;
			d2LnGdT2[j] = ( (d2UdT2*V + dUdT*dVdT)/pow(V,2.) - (dUdT*V)*(2.*dVdT)/pow(V,3.)
				- (dUdT*dVdT + U*d2VdT2)/pow(V,2.) + (U*dVdT)*(2.*dVdT)/pow(V,3.)
				+ d2bgdT2*IS ) * lg_to_ln;
			dLnGdP[j] = ( (dUdP*V - U*dVdP)/pow(V,2.) + dbgdP*IS ) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				if ( flagNeut == 1 )
				{
					LnG[j] = ( bgam*IS ) * lg_to_ln;
					dLnGdT[j] = ( dbgdT*IS ) * lg_to_ln;
					d2LnGdT2[j] = ( d2bgdT2*IS ) * lg_to_ln;
					dLnGdP[j] = ( dbgdP*IS ) * lg_to_ln;
				}

				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
				continue;
			}

			// water solvent
			else
			{
				// H2O activity coefficient calculated
				if ( flagH2O == 1 )
				{
					Phit = 0.; dPhitdT = 0.; d2PhitdT2 = 0.; dPhitdP = 0.;

					// derivatives of lambda and sigma terms
					L = 1. + (ao*B) * sqI;
					dLdT = ( daodT*B + ao*dBdT ) * sqI;
					d2LdT2 = ( d2aodT2*B + 2.*daodT*dBdT + ao*d2BdT2 ) * sqI;
					dLdP = ( daodP*B + ao*dBdP ) * sqI;

					U1 = (A*L);
					dU1dT = (dAdT*L + A*dLdT);
					d2U1dT2 = ( d2AdT2*L + 2.*dAdT*dLdT + A*d2LdT2 );
					dU1dP = ( dAdP*L + A*dLdP );
					V1 = pow(ao,3.)*pow(B,3.) * IS;
					dV1dT = ( 3.*pow(ao,2.)*daodT*pow(B,3.) + 3.*pow(ao,3.)*pow(B,2.)*dBdT ) * IS;
					d2V1dT2 = ( 6.*ao*pow(daodT,2.)*pow(B,3.) + 3.*pow(ao,2.)*d2aodT2*pow(B,3.)
								+ 18.*pow(ao,2.)*daodT*pow(B,2.)*dBdT + 6.*pow(ao,3.)*B*pow(dBdT,2.)
								+ 3.*pow(ao,3.)*pow(B,2.)*d2BdT2 ) * IS;
					dV1dP = ( 3.*pow(ao,2.)*daodP*pow(B,3.) + 3.*pow(ao,3.)*pow(B,2.)*dBdP ) * IS;

					U2 = A;
					dU2dT = dAdT;
					d2U2dT2 = d2AdT2;
					dU2dP = dAdP;
					V2 = pow(ao,3.)*pow(B,3.)*L * IS;
					dV2dT = ( 3.*pow(ao,2.)*daodT*pow(B,3.)*L + 3.*pow(ao,3.)*pow(B,2.)*dBdT*L
								+ pow(ao,3.)*pow(B,3.)*dLdT ) * IS;
					d2V2dT2 = ( 6.*ao*pow(daodT,2.)*pow(B,3.)*L + 3.*pow(ao,2.)*d2aodT2*pow(B,3.)*L
								+ 18.*pow(ao,2.)*daodT*pow(B,2.)*dBdT*L + 6.*pow(ao,2.)*daodT*pow(B,3.)*dLdT
								+ 6.*pow(ao,3.)*B*pow(dBdT,2.)*L + 3.*pow(ao,3.)*pow(B,2.)*d2BdT2*L
								+ 6.*pow(ao,3.)*pow(B,2.)*dBdT*dLdT + pow(ao,3.)*pow(B,3.)*d2LdT2 ) * IS;
					dV2dP = ( 3.*pow(ao,2.)*daodP*pow(B,3.)*L + 3.*pow(ao,3.)*pow(B,2.)*dBdP*L
								+ pow(ao,3.)*pow(B,3.)*dLdP ) * IS;

					U3 = 2.*( A*log(L) );
					dU3dT = 2.*( dAdT*log(L) + A*(1./L)*dLdT );
					d2U3dT2 = 2.*( d2AdT2*log(L) + 2.*dAdT*(1./L)*dLdT
								- A*(1./pow(L,2.))*pow(dLdT,2.) + A*(1./L)*d2LdT2 );
					dU3dP = 2.*( dAdP*log(L) + A*(1./L)*dLdP );
					V3 = pow(ao,3.)*pow(B,3.) * IS;
					dV3dT = ( 3.*pow(ao,2.)*daodT*pow(B,3.) + 3.*pow(ao,3.)*pow(B,2.)*dBdT ) * IS;
					d2V3dT2 = ( 6.*ao*pow(daodT,2.)*pow(B,3.) + 3.*pow(ao,2.)*d2aodT2*pow(B,3.)
								+ 18.*pow(ao,2.)*daodT*pow(B,2.)*dBdT + 6.*pow(ao,3.)*B*pow(dBdT,2.)
								+ 3.*pow(ao,3.)*pow(B,2.)*d2BdT2 ) * IS;
					dV3dP = ( 3.*pow(ao,2.)*daodP*pow(B,3.) + 3.*pow(ao,3.)*pow(B,2.)*dBdP ) * IS;

					Z = U1/V1 - U2/V2 - U3/V3;
					dZdT = (dU1dT*V1 - U1*dV1dT)/pow(V1,2.) - (dU2dT*V2 - U2*dV2dT)/pow(V2,2.)
								- (dU3dT*V3 - U3*dV3dT)/pow(V3,2.);
					d2ZdT2 = (d2U1dT2*V1 + dU1dT*dV1dT)/pow(V1,2.) - (dU1dT*V1)*(2.*dV1dT)/pow(V1,3.)
								- (dU1dT*dV1dT + U1*d2V1dT2)/pow(V1,2.) + (U1*dV1dT)*(2.*dV1dT)/pow(V1,3.)
								- (d2U2dT2*V2 + dU2dT*dV2dT)/pow(V2,2.) + (dU2dT*V2)*(2.*dV2dT)/pow(V2,3.)
								+ (dU2dT*dV2dT + U2*d2V2dT2)/pow(V2,2.) - (U2*dV2dT)*(2.*dV2dT)/pow(V2,3.)
								- (d2U3dT2*V3 + dU3dT*dV3dT)/pow(V3,2.) + (dU3dT*V3)*(2.*dV3dT)/pow(V3,3.)
								+ (dU3dT*dV3dT + U3*d2V3dT2)/pow(V3,2.) - (U3*dV3dT)*(2.*dV3dT)/pow(V3,3.);
					dZdP = (dU1dP*V1 - U1*dV1dP)/pow(V1,2.) - (dU2dP*V2 - U2*dV2dP)/pow(V2,2.)
								- (dU3dP*V3 - U3*dV3dP)/pow(V3,2.);

					// increments to osmotic coefficient (and derivatives)
					Lam = 1. + (ao*B) * sqI;
					sig = 3./(pow(ao,3.)*pow(B,3.)*pow(IS,(3./2.))) * (Lam-1./Lam-2*log(Lam));

					for (k=0; k<(NComp-1); k++)
					{
						if ( (z[k] == 0) && (flagNeut == 0) )
						{
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) );
							dPhitdT  += - log(10.) * m[k] * ( pow(z[k],2.)*dZdT );
							d2PhitdT2 += - log(10.) * m[k] * ( pow(z[k],2.)*d2ZdT2 );
							dPhitdP += - log(10.) * m[k] * ( pow(z[k],2.)*dZdP );
						}

						else
						{
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) - bgam*IS/2. );
							dPhitdT  += - log(10.) * m[k] * ( pow(z[k],2.)*dZdT - dbgdT*IS/2. );
							d2PhitdT2 += - log(10.) * m[k] * ( pow(z[k],2.)*d2ZdT2 - d2bgdT2*IS/2. );
							dPhitdP += - log(10.) * m[k] * ( pow(z[k],2.)*dZdP - dbgdP*IS/2. );
						}
					}

					Phi = Phit/molT;
					dPhidT = dPhitdT/molT;
					d2PhidT2 = d2PhitdT2/molT;
					dPhidP = dPhitdP/molT;

					// activity coefficient (and derivatives)
					lnaw = - Phi*molT/Nw;
					LnG[j] = lnaw - lnxw;
					dLnGdT[j] = - (molT/Nw) * dPhidT;
					d2LnGdT2[j] = - (molT/Nw) * d2PhidT2;
					dLnGdP[j] = - (molT/Nw) * dPhidP;

				}

				// H2O activity coefficient 0
				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
			}
		}

		g += x[j]*LnG[j];
		dgt += x[j]*dLnGdT[j];
		d2gt += x[j]*d2LnGdT2[j];
		dgp += x[j]*dLnGdP[j];

	} // j

	// increment thermodynamic properties
	Gex = (R_CONST*Tk) * g;
	Hex = - R_CONST*pow(Tk,2.) * dgt;
	// Sex = - R_CONST * ( g + Tk*dgt );
	Sex = (Hex-Gex)/Tk;
	CPex = - R_CONST * ( 2.*Tk*dgt + pow(Tk,2.)*d2gt );
	Vex = (R_CONST*Tk) * dgp;
	Aex = Gex - Vex*Pbar;
	Uex = Hex - Vex*Pbar;

	// assigments (excess properties)
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
long int THelgeson::IdealProp( double *Zid )
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


/// calculates true ionic strength
long int THelgeson::IonicStrength()
{
	long int j;
	double is, mt, mz;
	is = 0.0; mt = 0.0; mz = 0.0;

	// calculate ionic strength and total molalities
	// needs to be modified when nonaqueous solvents are present
	for (j=0; j<(NComp-1); j++)
	{
		is += 0.5*m[j]*z[j]*z[j];
		mt += m[j];
		if ( z[j] )
			mz += m[j];
	}

	// assignments
	IS = is;
	molT = mt;
	molZ = mz;

	return 0;
}


/// calculates TP dependence of b_gamma (and derivatives)
long int THelgeson::BgammaTP()
{
	// ni: stoichiometric number of moles of ions in one mole of electrolyte
	// rc, ra: radius of cation and anion, respectively at 298 K/1 bar
	// units are cal, kg, K, mol, bar
	double ni, nc, na, zc, za, rc, ra, a1, a2, a3, a4, a5, c1, c2, omg, bg, bs, bh, rec, rea,
			omgpt, domdt, d2omdt2, domdp, nbg, nbv, nbj, nbh;
	double eps, eta, xborn, yborn, qborn, X1, X2;

	// set parameters
	eps = EpsW[0];
	yborn = 1./pow(EpsW[0],2.)*EpsW[1];
	xborn = 1./pow(EpsW[0],2.) * ( EpsW[2] - 2./EpsW[0]*pow(EpsW[1],2.) );
	qborn = 1./pow(EpsW[0],2.)*EpsW[3];
	eta = (1.66027e5);

	switch ( flagElect )
	{
		case 1:  // NaCl
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.030056; a2 = -202.55; a3 = -2.9092; a4 = 20302;
			a5 = -0.206; c1 = -1.50; c2 = 53300.; omg = 178650.;
			bg = -174.623; bs = 2.164; rc = 0.97; ra = 1.81;
			break;
		case 2:  // KCl
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.0172; a2 = -115.36; a3 = -1.1857; a4 = 13854.2;
			a5 = -0.262; c1 = -2.53; c2 = 38628.4; omg = 164870.;
			bg = -70.0; bs = 1.727; rc = 1.33; ra = 1.81;
			break;
		case 3:  // NaOH
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.030056; a2 = -202.55; a3 = -2.9092; a4 = 20302;
			a5 = -0.206; c1 = -1.50; c2 = 53300.; omg = 205520.;
			bg = -267.4; bs = 1.836; rc = 0.97; ra = 1.40;
			break;
		case 4:  // KOH
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.0172; a2 = -115.36; a3 = -1.1857; a4 = 13854.2;
			a5 = -0.262; c1 = -2.53; c2 = 38628.4; omg = 191730.;
			bg = -335.7; bs = 1.26; rc = 1.33; ra = 1.40;
			break;
		default:  // wrong mode
			return -1;
	}

	// calculation part, extended 06.06.2009 (TW)
	bh = bg + (298.15)*bs;
	rec = rc + fabs(zc)*(0.94+Gf);
	rea = ra + fabs(za)*Gf;
	X1 = - eta*nc*( fabs(pow(zc,3.))/pow(rec,2.) - zc/pow((3.082+Gf),2.) )
			- eta*na*( fabs(pow(za,3.))/pow(rea,2.) - za/pow((3.082+Gf),2.) );
	X2 = 2.*eta*nc*( fabs(pow(zc,4.))/pow(rec,3.) - zc/pow((3.082+Gf),3.) )
			+ 2.*eta*nc * ( fabs(pow(zc,4.))/pow(rec,3.) - zc/pow((3.082+Gf),3.) );
	omgpt = eta*( nc*pow(zc,2.)/rec + na*pow(za,2.)/rea );
	// omgpt = (1.66027e5)*(1./(0.94+rc+Gf)+1./(ra+Gf));
	domdt = X1*dGfdT;
	d2omdt2 = X2*pow(dGfdT,2.) + X1*d2GfdT2;
	domdp = X1*dGfdP;
	nbg = - ni*bg/2. + ni*bs*(Tk-298.15)/2. - c1*(Tk*log(Tk/298.15)-Tk+298.15)
				+ a1*(Pbar-1.) + a2*log((2600.+Pbar)/(2600.+1.))
				- c2*((1./(Tk-228.)-1./(298.15-228.))*(228.-Tk)/228.-Tk/(228.*228.)
				* log((298.15*(Tk-228.))/(Tk*(298.15-228.))))
				+ 1./(Tk-228.)*(a3*(Pbar-1.) + a4*log((2600.+Pbar)/(2600.+1.)))
				+ a5*(omgpt*(1./eps-1.)-omg*(1./(78.24513795)-1.)+(-5.798650444e-5)*omg*(Tk-298.15));
	nbh = - ni/2.*bh + c1*(Tk-298.15) - c2*(1./(Tk-228.)-1./(298.15-228.))
				+ a1*(Pbar-1.) + a2*log((2600.+Pbar)/(2600.+1.))
				+ ((2.*Tk-228.)/pow((Tk-228.),2.))*(a3*(Pbar-1.)+a4*log((2600.+Pbar)/(2600.+1.)))
				+ a5*( omgpt*(1./eps-1.) + omgpt*Tk*yborn - Tk*(1./eps-1.)*domdt
				- omg*(1./(78.24513795)-1.) - omg*(298.15)*(-5.798650444e-5) );
	nbj = c1 + c2/pow((Tk-228.),2.)  - ( (2.*Tk)/pow((Tk-228.),3.) )
				* ( a3*(Pbar-1.) + a4*log((2600.+Pbar)/(2600.+1.)) )
				+ a5*( omgpt*Tk*xborn + 2.*Tk*yborn*domdt - Tk*(1./eps-1.)*d2omdt2 );
	nbv = a1 + a2/(2600.+Pbar) + a3/(Tk-228.) + a4/((2600.+Pbar)*(Tk-228.))
				+ a5*(-omgpt*qborn + (1./eps-1.)*domdp);
	// b_gamma = nbg/(2.*log(10.)*1.98721*Tk);
	// bgam = b_gamma;
	bgam = nbg/(2.*log(10.)*(1.98721)*Tk)*2./ni;
	dbgdT = - nbh/(2.*log(10.)*(1.98721)*pow(Tk,2.))*2./ni;
	d2bgdT2 = - nbj/( 2.*log(10.)*(1.98721)*pow(Tk,2.))*2./ni - 2./Tk*dbgdT;
	dbgdP = nbv/(2.*log(10.)*(1.98721)*Tk)*2./ni;

	return 0;
}


/// calculates TP dependence of a_not (and derivatives)
long int THelgeson::IonsizeTP()
{
	double nc, na, ni, zc, za, c;

	switch ( flagElect )
	{
		case 1:  // NaCl
			nc = 1.; na = 1.; ni = 2.;
			zc = 1.; za = -1.;
			break;
		case 2:  // KCl
			nc = 1.; na = 1.; ni = 2.;
			zc = 1.; za = -1.;
			break;
		case 3:  // NaOH
			nc = 1.; na = 1.; ni = 2.;
			zc = 1.; za = -1.;
			break;
		case 4:  // KOH
			nc = 1.; na = 1.; ni = 2.;
			zc = 1.; za = -1.;
			break;
		default:  // wrong mode
			return -1;
	}

	c = 2./ni * ( nc*fabs(zc) + na*fabs(za) );
	ao = ac + c*Gf;
	daodT = c*dGfdT;
	d2aodT2 = c*d2GfdT2;
	daodP = c*dGfdP;

	return 0;
}


/// wrapper for g-function
long int THelgeson::Gfunction()
{
	double T, P, D, beta, alpha, daldT, g, dgdP, dgdT, d2gdT2;
	double TMAX = 1000., PMAX = 5000., TOL = 1.0e-4;

	// convert parameters
	T = Tk - 273.15;
	P = Pbar;
	D = RhoW[0];
	alpha = - RhoW[1]/RhoW[0];
	daldT = pow(alpha, 2.) - RhoW[2]/RhoW[0];
	beta = RhoW[3]/RhoW[0];

	// initialize g and derivatives to zero
	g = 0.0;
	dgdP = 0.0;
	dgdT = 0.0;
	d2gdT2 = 0.0;

	if ((T > TMAX+TOL) || (P > PMAX+TOL))
		return -1;
	else
		GShok2( T, P, D, beta, alpha, daldT, g, dgdP, dgdT, d2gdT2 );

	// assignments
	Gf = g;
	dGfdT = dgdT;
	d2GfdT2 = d2gdT2;
	dGfdP = dgdP;

	return 0;
}


/// calculates g-function and derivatives
long int THelgeson::GShok2( double T, double P, double D, double beta,
		double alpha, double daldT, double &g, double &dgdP, double &dgdT, double &d2gdT2 )
{
	double a, b, dgdD, /*dgdD2,*/ dadT, dadTT, dbdT, dbdTT, dDdT, dDdP,
			dDdTT, Db, dDbdT, dDbdTT, ft, dftdT, dftdTT, fp, dfpdP,
			f, dfdP, dfdT, d2fdT2, tempy;
	double C[6]  = {-0.2037662e+01,  0.5747000e-02, -0.6557892e-05,
			0.6107361e+01, -0.1074377e-01,  0.1268348e-04 };
	double cC[3] = { 0.3666666e+02, -0.1504956e-09,  0.5017997e-13 };

	if ( D >= 1.3 )
		return -1;
	// Sveta 19/02/2000 1-D < 0 pri D>1 => pow(-number, b) = -NaN0000 error
	double pw = fabs(1.0 - D); // insert Sveta 19/02/2000

	// calculation part
	a = C[0] + C[1]*T + C[2]*pow(T,2.);
	b = C[3] + C[4]*T + C[5]*pow(T,2.);
	g = a * pow(pw, b);

	dgdD = - a*b* pow(pw, (b - 1.0));
	// dgdD2 = a * b * (b - 1.0) * pow((1.0 - D),(b - 2.0));

	dadT = C[1] + 2.0*C[2]*T;
	dadTT = 2.0*C[2];
	dbdT = C[4] + 2.0*C[5]*T;
	dbdTT = 2.0*C[5];
	dDdT = - D * alpha;
	dDdP = D * beta;
	dDdTT = - D * (daldT - pow(alpha,2.));

	// Db = pow((1.0 - D), b);  Fixed by DAK 01.11.00
	Db = pow(pw , b);
	dDbdT = - b*pow(pw,(b - 1.0))*dDdT + log(pw)*Db*dbdT;

	dDbdTT = - (b*pow(pw,(b-1.0)) *dDdTT + pow(pw,(b - 1.0)) * dDdT * dbdT
				+ b * dDdT * (-(b - 1.0) * pow(pw,(b - 2.0)) * dDdT
				+ log(pw) * pow(pw,(b - 1.0)) * dbdT))
				+ log(pw) * pow(pw,b) * dbdTT
				- pow(pw,b) * dbdT * dDdT / (1.0 - D)
				+ log(pw) * dbdT * dDbdT;

	dgdP = dgdD * dDdP;
	dgdT = a * dDbdT + Db * dadT;
	d2gdT2 = a * dDbdTT + 2.0 * dDbdT * dadT + Db * dadTT;

	if((T < 155.0) || (P > 1000.0) || (T > 355.0))
		return 0;

	tempy = ((T - 155.0) / 300.0);
	ft = pow(tempy,4.8) + cC[0] * pow(tempy,16.);

	dftdT = 4.8 / 300.0 * pow(tempy, 3.8) + 16.0 / 300.0 * cC[0] * pow(tempy, 15.);

	dftdTT = 3.8 * 4.8 / (300.0 * 300.0) * pow(tempy, 2.8)
		+ 15.0 * 16.0 / (300.0*300.0) * cC[0] * pow(tempy,14.);

	fp = cC[1] * pow((1000.0 - P),3.) + cC[2] * pow((1000.0 - P),4.);

	dfpdP = -3.0 * cC[1] * pow((1000.0 - P),2.) - 4.0 * cC[2] * pow((1000.0 - P),3.);

	f = ft * fp;
	dfdP = ft * dfpdP;
	dfdT = fp * dftdT;
	d2fdT2 = fp * dftdTT;

	g -= f;
	dgdP -= dfdP;
	dgdT -= dfdT;
	d2gdT2 -= d2fdT2;

	return 0;
}





//=============================================================================================
// Extended Debye-Hueckel (EDH) model for aqueous solutions, Davies version
// References: Langmuir (1997)
// (c) TW May 2009
//=============================================================================================


// Generic constructor for the TDavies class
TDavies::TDavies( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW ):
                TSolMod( sd )
{
    alloc_internal();
    m = arM;
    z = arZ;
    RhoW = dW;
    EpsW = eW;
    flagNeut = 0;
    flagH2O = (long int)aIPc[3];  // 0: unity, 1: calculated
    flagMol = (long int)aIPc[5];  // 0: no scale correction, 1: scale correction
}


TDavies::~TDavies()
{
	free_internal();
}


void TDavies::alloc_internal()
{
	LnG = new double [NComp];
	dLnGdT = new double [NComp];
	d2LnGdT2 = new double [NComp];
	dLnGdP = new double [NComp];
}


void TDavies::free_internal()
{
	delete[]LnG;
	delete[]dLnGdT;
	delete[]d2LnGdT2;
	delete[]dLnGdP;
}


/// Calculates T,P corrected parameters
long int TDavies::PTparam()
{
	double alp, bet, dal, rho, eps, dedt, d2edt2, dedp;

	// read and convert rho and eps
	rho = RhoW[0];
	alp = - 1./rho*RhoW[1];
	dal = pow(alp,2.) - 1./rho*RhoW[2];
	bet = 1./rho*RhoW[3];
	eps = EpsW[0];
	dedt = 1./eps*EpsW[1];
	d2edt2 = - 1./pow(eps,2.)*pow(EpsW[1],2.) + 1./eps*EpsW[2];
	dedp = 1./eps*EpsW[3];

	// calculate A term of Debye-Huckel equation (and derivatives)
	A = (1.82483e6)*sqrt(rho) / pow(Tk*eps,1.5);
	dAdT = - 3./2.*A*( dedt + 1./Tk + alp/3. );
	d2AdT2 = 1./A*pow(dAdT,2.) - 3./2.*A*( d2edt2 - 1/pow(Tk,2.) + 1/3.*dal );
	dAdP = 1./2.*A*( bet - 3.*dedp);

 	ErrorIf( fabs(A) < 1e-9, "Davies model",
 			"Error: DH parameter A was not calculated - no values of RhoW and EpsW !" );

	return 0;
}


/// Calculates activity coefficients
long int TDavies::MixMod()
{
	long int j, k, w;
	double sqI, Z2, lgGam, lnGam, Nw, Lgam, lnxw, xw, lgaw, lnaw,
				sig, zt, zz, lg_to_ln;
	lg_to_ln = 2.302585093;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		lgGam = 0.0;
		lnGam = 0.0;

		// charged species
		if ( z[j] )
		{
			lgGam = 0.0;
			Z2 = z[j]*z[j];
			lgGam = ( -A * Z2 ) * ( sqI/( 1. + sqI ) - 0.3 * IS );
			if ( flagMol == 0 )
				lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
			else
				lnGamma[j] = lgGam * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				lgGam = 0.0;
				if ( flagMol == 0 )
					lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
				else
					lnGamma[j] = lgGam * lg_to_ln;
				continue;
			}

			// water solvent
			else
			{
				lgGam = 0.0;
				lnGam = 0.0;
				if ( flagH2O == 1 )
				{
					// equation from Wolery (1990), corrected for species charges terms
					zt = 0.;
					sig = 3./sqI * ( 1. + sqI - 1./(1.+sqI) - 2.*log(1.+sqI) );

					for (k=0; k<(NComp-1); k++)
					{
						zt += m[k]*pow(z[k],2.);
					}

					zz = zt/molT;
					lgaw = (1./Nw) * ( - molT/log(10.) + 2./3.*A*sqI*sig - 2.*(-0.3/zz)*A*pow(IS,2.) );
					// lgaw = (1./Nw) * ( - molT/log(10.) + 2./3.*A*sqI*sig - 2.*(-0.3/zz)*A*pow(IS,2.) );
					if ( flagMol == 0 )
						lgaw += molT/(Nw)/log(10.) + lnxw/log(10.);
					lnaw = lgaw * lg_to_ln;
					lnGam = lnaw - lnxw;
				}
				else
					lnGam = 0.0;
				lnGamma[j] = lnGam;
			}
		}
	} // j

	return 0;
}


/// calculates excess properties
long int TDavies::ExcessProp( double *Zex )
{
	long int j, k, w;
	double sqI, Z2, Nw, Lgam, lnxw, xw, lnaw, lgGam, lgaw, sig, zt, zz,
				dawdt, d2awdt2, dawdp, lg_to_ln, g, dgt, d2gt, dgp;
	lg_to_ln = 2.302585093;
	g = 0.; dgt = 0.; d2gt = 0.; dgp = 0.;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		// charged species
		if ( z[j] )
		{
			lgGam = 0.0;
			Z2 = z[j]*z[j];
			lgGam = - ( A * Z2 ) * ( sqI/( 1. + sqI ) - 0.3 * IS );
			if ( flagMol == 0 )
				LnG[j] = lgGam * lg_to_ln;
			else
				LnG[j] = (lgGam - Lgam) * lg_to_ln;
			dLnGdT[j] = - ( dAdT * Z2 ) * ( sqI/( 1. + sqI ) - 0.3 * IS ) * lg_to_ln;
			d2LnGdT2[j] = - ( d2AdT2 * Z2 ) * ( sqI/( 1. + sqI ) - 0.3 * IS ) * lg_to_ln;
			dLnGdP[j] = - ( dAdP * Z2 ) * ( sqI/( 1. + sqI ) - 0.3 * IS ) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				lgGam = 0.0;
				if ( flagMol == 0 )
					LnG[j] = 0.0;
				else
					LnG[j] = (lgGam - Lgam) * lg_to_ln;
				dLnGdT[j] = 0.;
				d2LnGdT2[j] = 0.;
				dLnGdP[j] = 0.;
			}

			// water solvent
			else
			{
				// H2O activity coefficient calculated
				if ( flagH2O == 1 )
				{
					// equation from Wolery (1990), corrected for species charges terms
					zt = 0.;
					sig = 3./sqI * ( 1. + sqI - 1./(1.+sqI) - 2.*log(1.+sqI) );

					for (k=0; k<(NComp-1); k++)
					{
						zt += m[k]*pow(z[k],2.);
					}

					zz = zt/molT;
					// lgaw = (1./Nw) * ( - molT/log(10.) + 2./3.*A*sqI*sig - 2.*(-0.3)*A*pow(IS,2.) );
					lgaw = (1./Nw) * ( - molT/log(10.) + 2./3.*A*sqI*sig - 2.*(-0.3/zz)*A*pow(IS,2.) );
					if ( flagMol == 0 )
						lgaw += molT/(Nw)/log(10.) + lnxw/log(10.);
					lnaw = lgaw * lg_to_ln;
					dawdt = (1./Nw) * ( 2./3.*dAdT*sqI*sig - 2.*(-0.3/zz)*dAdT*pow(IS,2.) ) * lg_to_ln;
					d2awdt2 = (1./Nw) * ( 2./3.*d2AdT2*sqI*sig - 2.*(-0.3/zz)*d2AdT2*pow(IS,2.) ) * lg_to_ln;
					dawdp = (1./Nw) * ( 2./3.*dAdP*sqI*sig - 2.*(-0.3/zz)*dAdP*pow(IS,2.) ) * lg_to_ln;
					LnG[j] = lnaw - lnxw;
					dLnGdT[j] = dawdt;
					d2LnGdT2[j] = d2awdt2;
					dLnGdP[j] = dawdp;
				}

				// H2O activity coefficient 0
				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
			}
		}

		g += x[j]*LnG[j];
		dgt += x[j]*dLnGdT[j];
		d2gt += x[j]*d2LnGdT2[j];
		dgp += x[j]*dLnGdP[j];

	} // j

	// increment thermodynamic properties
	Gex = (R_CONST*Tk) * g;
	Hex = - R_CONST*pow(Tk,2.) * dgt;
	// Sex = - R_CONST * ( g + Tk*dgt );
	Sex = (Hex-Gex)/Tk;
	CPex = - R_CONST * ( 2.*Tk*dgt + pow(Tk,2.)*d2gt );
	Vex = (R_CONST*Tk) * dgp;
	Aex = Gex - Vex*Pbar;
	Uex = Hex - Vex*Pbar;

	// assigments (excess properties)
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
long int TDavies::IdealProp( double *Zid )
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


/// calculates true ionic strength
long int TDavies::IonicStrength()
{
	long int j;
	double is, mt;
	is = 0.0; mt = 0.0;

	// calculate ionic strength and total molalities
	// needs to be modified if nonaqueous solvents are present
	for (j=0; j<(NComp-1); j++)
	{
		is += 0.5*m[j]*z[j]*z[j];
		mt += m[j];
	}

	// assignments
	IS = is;
	molT = mt;

	return 0;
}





//=============================================================================================
// Debye-Hueckel (DH) limiting law for aqueous solutions
// References: Langmuir (1997)
// ((c) TW May 2009
//=============================================================================================


// Generic constructor for the TDLimitingLaw class
TLimitingLaw::TLimitingLaw( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW ):
                TSolMod( sd )
{
    alloc_internal();
    m = arM;
    z = arZ;
    RhoW = dW;
    EpsW = eW;
    flagNeut = (long int)aIPc[2];  // 0: unity
    flagH2O = (long int)aIPc[3];  // 0: unity, 1: calculated
}


TLimitingLaw::~TLimitingLaw()
{
	free_internal();
}


void TLimitingLaw::alloc_internal()
{
	LnG = new double [NComp];
	dLnGdT = new double [NComp];
	d2LnGdT2 = new double [NComp];
	dLnGdP = new double [NComp];
}


void TLimitingLaw::free_internal()
{
	delete[]LnG;
	delete[]dLnGdT;
	delete[]d2LnGdT2;
	delete[]dLnGdP;
}


/// Calculates T,P corrected parameters
long int TLimitingLaw::PTparam()
{
	double alp, bet, dal, rho, eps, dedt, d2edt2, dedp;

	// read and convert rho and eps
	rho = RhoW[0];
	alp = - 1./rho*RhoW[1];
	dal = pow(alp,2.) - 1./rho*RhoW[2];
	bet = 1./rho*RhoW[3];
	eps = EpsW[0];
	dedt = 1./eps*EpsW[1];
	d2edt2 = - 1./pow(eps,2.)*pow(EpsW[1],2.) + 1./eps*EpsW[2];
	dedp = 1./eps*EpsW[3];

	// calculate A term of Debye-Huckel equation (and derivatives)
	A = (1.82483e6)*sqrt(rho) / pow(Tk*eps,1.5);
	dAdT = - 3./2.*A*( dedt + 1./Tk + alp/3. );
	d2AdT2 = 1./A*pow(dAdT,2.) - 3./2.*A*( d2edt2 - 1/pow(Tk,2.) + 1/3.*dal );
	dAdP = 1./2.*A*( bet - 3.*dedp);

 	ErrorIf( fabs(A) < 1e-9, "DH limiting law model",
 			"Error: DH parameter A was not calculated - no values of RhoW and EpsW !" );

	return 0;
}


/// Calculates activity coefficients
long int TLimitingLaw::MixMod()
{
	long int j, k, w;
	double sqI, Z2, lgGam, lnGam, Nw, Lgam, lnxw, xw, lnaw, Phi, Phit, lg_to_ln;
	lg_to_ln = 2.302585093;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		lgGam = 0.0;
		lnGam = 0.0;

		// charged species
		if ( z[j] )
		{
			lgGam = 0.0;
			Z2 = z[j]*z[j];
			lgGam = ( -A * Z2 * sqI );
			lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				lgGam = 0.0;
				lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
				continue;
			}

			// water solvent
			else
			{
				lgGam = 0.0;
				lnGam = 0.0;
				if ( flagH2O == 1 )
				{
					Phit = 0.0;
					// Phi corrected using eq. (190) from Helgeson et al. (1981)

					for (k=0; k<(NComp-1); k++)
					{
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI)/3. + Lgam/(0.0180153*molT) );
					}

					Phi = Phit/molT;
					lnaw = - Phi*molT/Nw;
					lnGam = lnaw - lnxw;
				}
				else
					lnGam = 0.0;
				lnGamma[j] = lnGam;
			}
		}
	} // j

	return 0;
}


/// calculates excess properties
long int TLimitingLaw::ExcessProp( double *Zex )
{
	long int j, k, w;
	double sqI, Z2, Nw, Lgam, lnxw, xw, Phi, dPhidT, d2PhidT2, dPhidP,
			Phit, dPhitdT, d2PhitdT2, dPhitdP, lnaw, lg_to_ln, g, dgt, d2gt, dgp;
	lg_to_ln = 2.302585093;
	g = 0.; dgt = 0.; d2gt = 0.; dgp = 0.;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		// charged species
		if ( z[j] )
		{
			Z2 = z[j]*z[j];
			LnG[j] = - ( A * Z2 * sqI ) * lg_to_ln;
			dLnGdT[j] = - ( dAdT * Z2 * sqI ) * lg_to_ln;
			d2LnGdT2[j] = - ( d2AdT2 * Z2 * sqI ) * lg_to_ln;
			dLnGdP[j] = - ( dAdP * Z2 * sqI ) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				LnG[j] = 0.;
				dLnGdT[j] = 0.;
				d2LnGdT2[j] = 0.;
				dLnGdP[j] = 0.;
			}

			// water solvent
			else
			{
				if ( flagH2O == 1 )
				{
					Phit = 0.; dPhitdT = 0.; d2PhitdT2 = 0.; dPhitdP = 0.;
					// Phi corrected using eq. (190) from Helgeson et al. (1981)

					for (k=0; k<(NComp-1); k++)
					{
						Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI)/3. + Lgam/(0.0180153*molT) );
						dPhitdT += - log(10.) * m[k] * ( (pow(z[k],2.)*dAdT*sqI)/3. );
						d2PhitdT2 += - log(10.) * m[k] * ( (pow(z[k],2.)*d2AdT2*sqI)/3. );
						dPhitdP += - log(10.) * m[k] * ( (pow(z[k],2.)*dAdP*sqI)/3. );
					}

					Phi = Phit/molT;
					dPhidT = dPhitdT/molT;
					d2PhidT2 = d2PhitdT2/molT;
					dPhidP = dPhitdP/molT;

					// activity coefficient (and derivatives)
					lnaw = - Phi*molT/Nw;
					LnG[j] = lnaw - lnxw;
					dLnGdT[j] = - (molT/Nw) * dPhidT;
					d2LnGdT2[j] = - (molT/Nw) * d2PhidT2;
					dLnGdP[j] = - (molT/Nw) * dPhidP;
				}

				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
			}
		}

		g += x[j]*LnG[j];
		dgt += x[j]*dLnGdT[j];
		d2gt += x[j]*d2LnGdT2[j];
		dgp += x[j]*dLnGdP[j];

	} // j

	// increment thermodynamic properties
	Gex = (R_CONST*Tk) * g;
	Hex = - R_CONST*pow(Tk,2.) * dgt;
	// Sex = - R_CONST * ( g + Tk*dgt );
	Sex = (Hex-Gex)/Tk;
	CPex = - R_CONST * ( 2.*Tk*dgt + pow(Tk,2.)*d2gt );
	Vex = (R_CONST*Tk) * dgp;
	Aex = Gex - Vex*Pbar;
	Uex = Hex - Vex*Pbar;

	// assigments (excess properties)
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
long int TLimitingLaw::IdealProp( double *Zid )
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


/// calculates true ionic strength
long int TLimitingLaw::IonicStrength()
{
	long int j;
	double is, mt;
	is = 0.0; mt = 0.0;

	// calculate ionic strength and total molalities
	// needs to be modified when nonaqueous solvents are present
	for (j=0; j<(NComp-1); j++)
	{
		is += 0.5*m[j]*z[j]*z[j];
		mt += m[j];
	}

	// assignments
	IS = is;
	molT = mt;

	return 0;
}





//=============================================================================================
// Debye-Hueckel (DH) two term model for aqueous solutions
// References: Helgeson et al. (1981)
// uses individual ion-size parameters, optionally individual salting-out coefficients
// (c) TW May 2009
//=============================================================================================


// Generic constructor for the TDebyeHueckel class
TDebyeHueckel::TDebyeHueckel( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW ):
                TSolMod( sd )
{
    alloc_internal();
    m = arM;
    z = arZ;
    RhoW = dW;
    EpsW = eW;
    ac = 3.72;   // common ion size parameter
    bc = 0.064;   // common b_setch
    flagNeut = (long int)aIPc[2];   // 0: unity, 1: calculated from bg
    flagH2O = (long int)aIPc[3];   // 0: unity, 1: calculated from bg
}


TDebyeHueckel::~TDebyeHueckel()
{
	free_internal();
}


void TDebyeHueckel::alloc_internal()
{
	LnG = new double [NComp];
	dLnGdT = new double [NComp];
	d2LnGdT2 = new double [NComp];
	dLnGdP = new double [NComp];
	an = new double [NComp];
	bg = new double [NComp];
}


void TDebyeHueckel::free_internal()
{
	delete[]LnG;
	delete[]dLnGdT;
	delete[]d2LnGdT2;
	delete[]dLnGdP;
	delete[]an;
	delete[]bg;
}


/// Calculates T,P corrected parameters
long int TDebyeHueckel::PTparam()
{
	long int j;
	double alp, bet, dal, rho, eps, dedt, d2edt2, dedp;

	// read and copy individual an and bg
	for (j=0; j<NComp; j++)
	{
		an[j] = aDCc[NP_DC*j];   // individual an
		bg[j] = aDCc[NP_DC*j+1];   // individual bg
	}

	// read and convert rho and eps
	rho = RhoW[0];
	alp = - 1./rho*RhoW[1];
	dal = pow(alp,2.) - 1./rho*RhoW[2];
	bet = 1./rho*RhoW[3];
	eps = EpsW[0];
	dedt = 1./eps*EpsW[1];
	d2edt2 = - 1./pow(eps,2.)*pow(EpsW[1],2.) + 1./eps*EpsW[2];
	dedp = 1./eps*EpsW[3];

	// calculate A and B terms of Debye-Huckel equation (and derivatives)
	A = (1.82483e6)*sqrt(rho) / pow(Tk*eps,1.5);
	dAdT = - 3./2.*A*( dedt + 1./Tk + alp/3. );
	d2AdT2 = 1./A*pow(dAdT,2.) - 3./2.*A*( d2edt2 - 1/pow(Tk,2.) + 1/3.*dal );
	dAdP = 1./2.*A*( bet - 3.*dedp);
	B = (50.2916)*sqrt(rho) / sqrt(Tk*eps);
	dBdT = - 1./2.*B*( dedt +1./Tk + alp );
	d2BdT2 = 1./B*pow(dBdT,2.) - 1./2.*B*( d2edt2 - 1./pow(Tk,2.) + dal );
	dBdP = 1./2.*B*( bet - dedp );

 	ErrorIf( fabs(A) < 1e-9 || fabs(B) < 1e-9, "DH two-term model",
 			"Error: DH parameter A or B was not calculated - no values of RhoW and EpsW !" );

	return 0;
}


/// Calculates activity coefficients
long int TDebyeHueckel::MixMod()
{
	long int j, k, w;
	double sqI, Z2, lgGam, lnGam, Nw, Lgam, Lam, lnaw, lnxw, xw, lg_to_ln,
				sig, Phi, Phit;
	lg_to_ln = 2.302585093;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		lgGam = 0.0;
		lnGam = 0.0;

		// charged species
		if ( z[j] )
		{
			lgGam = 0.0;
			Z2 = z[j]*z[j];
			lgGam = ( - A * sqI * Z2 ) / ( 1. + B * an[j] * sqI ) ;
			lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				lgGam = 0.0;

				// rational Setchenow coefficient
				if ( flagNeut == 1 )
					lgGam = bg[j] * IS;

				else
					lgGam = 0.0;

				lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
				continue;
			}

			// water solvent
			else
			{
				lgGam = 0.0;
				lnGam = 0.0;

				// water activity coefficient calculated
				if ( flagH2O == 1 )
				{
					Phit = 0.0;
					// Phi corrected using eq. (190) from Helgeson et al. (1981)
					Lam = 1. + (ao*B) * sqI;
					sig = 3./(pow(ao,3.)*pow(B,3.)*pow(IS,(3./2.))) * (Lam-1./Lam-2*log(Lam));
					// Phi = -2.3025851*(A*sqI*SigTerm/3. + Lgam/(0.0180153*2.*IS) - bgam*IS/2.);
					// Phi = - log(10.) * (molZ/molT) * ( (zc*za*A*sqI*SigTerm)/3. + (Lgam*psi)/(0.0180153*2.*IS) - bgam*IS/2. );

					for (k=0; k<(NComp-1); k++)
					{
						if ( (z[k] == 0) && (flagNeut == 1) )
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) - bg[k]*IS/2. );
						else
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) );
					}

					Phi = Phit/molT;
					lnaw = - Phi*molT/Nw;
					lnGam = lnaw - lnxw;
				}
				else
					lnGam = 0.0;
				lnGamma[j] = lnGam;
			}
		}
	} // j

	return 0;
}


/// calculates excess properties
long int TDebyeHueckel::ExcessProp( double *Zex )
{
	// (under construction)
	long int j, k, w;
	double sqI, Z2, Nw, Lgam, lnxw, xw, Lam, sig, Phi, dPhidT, d2PhidT2, dPhidP,
			Phit, dPhitdT, d2PhitdT2, dPhitdP, lnaw, lg_to_ln, g, dgt, d2gt, dgp;
	double U, V, dUdT, dVdT, d2UdT2, d2VdT2, dUdP, dVdP, U1, U2, U3, V1, V2, V3,
			dU1dT, dU2dT, dU3dT, dV1dT, dV2dT, dV3dT, d2U1dT2, d2U2dT2, d2U3dT2,
			d2V1dT2, d2V2dT2, d2V3dT2, dU1dP, dU2dP, dU3dP, dV1dP, dV2dP, dV3dP,
			L, dLdT, d2LdT2, dLdP, Z, dZdT, d2ZdT2, dZdP;
	lg_to_ln = 2.302585093;
	g = 0.; dgt = 0.; d2gt = 0.; dgp = 0.;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		// charged species
		if ( z[j] )
		{
			Z2 = z[j]*z[j];
			U = - (Z2*A) * sqI;
			dUdT = - (Z2*dAdT) * sqI;
			d2UdT2 = - (Z2*d2AdT2) * sqI;
			dUdP = - (Z2*dAdP) * sqI;
			V = 1. + (an[j]*B) * sqI;
			dVdT = ( an[j]*dBdT ) * sqI;
			d2VdT2 = ( an[j]*d2BdT2 ) * sqI;
			dVdP = ( an[j]*dBdP ) * sqI;
			LnG[j] = ( ( - A * sqI * Z2 ) / ( 1. + B * an[j] * sqI ) ) * lg_to_ln;
			dLnGdT[j] = ( (dUdT*V - U*dVdT)/pow(V,2.) ) * lg_to_ln;
			d2LnGdT2[j] = ( (d2UdT2*V + dUdT*dVdT)/pow(V,2.) - (dUdT*V)*(2.*dVdT)/pow(V,3.)
				- (dUdT*dVdT + U*d2VdT2)/pow(V,2.) + (U*dVdT)*(2.*dVdT)/pow(V,3.) ) * lg_to_ln;
			dLnGdP[j] = ( (dUdP*V - U*dVdP)/pow(V,2.) ) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				// rational Setchenow coefficient
				if ( flagNeut == 1 )
				{
					LnG[j] = ( bg[j] * IS ) * lg_to_ln;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}

				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
				continue;
			}

			// water solvent
			else
			{
				// H2O activity coefficient calculated
				if ( flagH2O == 1 )
				{
					Phit = 0.; dPhitdT = 0.; d2PhitdT2 = 0.; dPhitdP = 0.;

					// derivatives of lambda and sigma terms
					L = 1. + (ao*B) * sqI;
					dLdT = ( ao*dBdT ) * sqI;
					d2LdT2 = ( ao*d2BdT2 ) * sqI;
					dLdP = ( ao*dBdP ) * sqI;

					U1 = (A*L);
					dU1dT = (dAdT*L + A*dLdT);
					d2U1dT2 = ( d2AdT2*L + 2.*dAdT*dLdT + A*d2LdT2 );
					dU1dP = ( dAdP*L + A*dLdP );
					V1 = pow(ao,3.)*pow(B,3.) * IS;
					dV1dT = ( 3.*pow(ao,3.)*pow(B,2.)*dBdT ) * IS;
					d2V1dT2 = ( 6.*pow(ao,3.)*B*pow(dBdT,2.) + 3.*pow(ao,3.)*pow(B,2.)*d2BdT2 ) * IS;
					dV1dP = ( 3.*pow(ao,3.)*pow(B,2.)*dBdP ) * IS;

					U2 = A;
					dU2dT = dAdT;
					d2U2dT2 = d2AdT2;
					dU2dP = dAdP;
					V2 = pow(ao,3.)*pow(B,3.)*L * IS;
					dV2dT = ( 3.*pow(ao,3.)*pow(B,2.)*dBdT*L + pow(ao,3.)*pow(B,3.)*dLdT ) * IS;
					d2V2dT2 = ( 6.*pow(ao,3.)*B*pow(dBdT,2.)*L + 3.*pow(ao,3.)*pow(B,2.)*d2BdT2*L
								+ 6.*pow(ao,3.)*pow(B,2.)*dBdT*dLdT + pow(ao,3.)*pow(B,3.)*d2LdT2 ) * IS;
					dV2dP = ( 3.*pow(ao,3.)*pow(B,2.)*dBdP*L + pow(ao,3.)*pow(B,3.)*dLdP ) * IS;

					U3 = 2.*( A*log(L) );
					dU3dT = 2.*( dAdT*log(L) + A*(1./L)*dLdT );
					d2U3dT2 = 2.*( d2AdT2*log(L) + 2.*dAdT*(1./L)*dLdT
								- A*(1./pow(L,2.))*pow(dLdT,2.) + A*(1./L)*d2LdT2 );
					dU3dP = 2.*( dAdP*log(L) + A*(1./L)*dLdP );
					V3 = pow(ao,3.)*pow(B,3.) * IS;
					dV3dT = ( 3.*pow(ao,3.)*pow(B,2.)*dBdT ) * IS;
					d2V3dT2 = ( 6.*pow(ao,3.)*B*pow(dBdT,2.) + 3.*pow(ao,3.)*pow(B,2.)*d2BdT2 ) * IS;
					dV3dP = ( 3.*pow(ao,3.)*pow(B,2.)*dBdP ) * IS;

					Z = U1/V1 - U2/V2 - U3/V3;
					dZdT = (dU1dT*V1 - U1*dV1dT)/pow(V1,2.) - (dU2dT*V2 - U2*dV2dT)/pow(V2,2.)
								- (dU3dT*V3 - U3*dV3dT)/pow(V3,2.);
					d2ZdT2 = (d2U1dT2*V1 + dU1dT*dV1dT)/pow(V1,2.) - (dU1dT*V1)*(2.*dV1dT)/pow(V1,3.)
								- (dU1dT*dV1dT + U1*d2V1dT2)/pow(V1,2.) + (U1*dV1dT)*(2.*dV1dT)/pow(V1,3.)
								- (d2U2dT2*V2 + dU2dT*dV2dT)/pow(V2,2.) + (dU2dT*V2)*(2.*dV2dT)/pow(V2,3.)
								+ (dU2dT*dV2dT + U2*d2V2dT2)/pow(V2,2.) - (U2*dV2dT)*(2.*dV2dT)/pow(V2,3.)
								- (d2U3dT2*V3 + dU3dT*dV3dT)/pow(V3,2.) + (dU3dT*V3)*(2.*dV3dT)/pow(V3,3.)
								+ (dU3dT*dV3dT + U3*d2V3dT2)/pow(V3,2.) - (U3*dV3dT)*(2.*dV3dT)/pow(V3,3.);
					dZdP = (dU1dP*V1 - U1*dV1dP)/pow(V1,2.) - (dU2dP*V2 - U2*dV2dP)/pow(V2,2.)
								- (dU3dP*V3 - U3*dV3dP)/pow(V3,2.);

					// increments to osmotic coefficient (and derivatives)
					Lam = 1. + (ao*B) * sqI;
					sig = 3./(pow(ao,3.)*pow(B,3.)*pow(IS,(3./2.))) * (Lam-1./Lam-2*log(Lam));

					for (k=0; k<(NComp-1); k++)
					{
						if ( (z[k] == 0) && (flagNeut == 1) )
						{
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) - bg[k]*IS/2. );
							dPhitdT  += - log(10.) * m[k] * ( pow(z[k],2.)*dZdT );
							d2PhitdT2 += - log(10.) * m[k] * ( pow(z[k],2.)*d2ZdT2 );
							dPhitdP += - log(10.) * m[k] * ( pow(z[k],2.)*dZdP );

						}

						else
						{
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) );
							dPhitdT  += - log(10.) * m[k] * ( pow(z[k],2.)*dZdT );
							d2PhitdT2 += - log(10.) * m[k] * ( pow(z[k],2.)*d2ZdT2 );
							dPhitdP += - log(10.) * m[k] * ( pow(z[k],2.)*dZdP );
						}
					}

					Phi = Phit/molT;
					dPhidT = dPhitdT/molT;
					d2PhidT2 = d2PhitdT2/molT;
					dPhidP = dPhitdP/molT;

					// activity coefficient (and derivatives)
					lnaw = - Phi*molT/Nw;
					LnG[j] = lnaw - lnxw;
					dLnGdT[j] = - (molT/Nw) * dPhidT;
					d2LnGdT2[j] = - (molT/Nw) * d2PhidT2;
					dLnGdP[j] = - (molT/Nw) * dPhidP;

				}

				// H2O activity coefficient 0
				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
			}
		}

		g += x[j]*LnG[j];
		dgt += x[j]*dLnGdT[j];
		d2gt += x[j]*d2LnGdT2[j];
		dgp += x[j]*dLnGdP[j];

	} // j

	// increment thermodynamic properties
	Gex = (R_CONST*Tk) * g;
	Hex = - R_CONST*pow(Tk,2.) * dgt;
	// Sex = - R_CONST * ( g + Tk*dgt );
	Sex = (Hex-Gex)/Tk;
	CPex = - R_CONST * ( 2.*Tk*dgt + pow(Tk,2.)*d2gt );
	Vex = (R_CONST*Tk) * dgp;
	Aex = Gex - Vex*Pbar;
	Uex = Hex - Vex*Pbar;

	// assigments (excess properties)
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
long int TDebyeHueckel::IdealProp( double *Zid )
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


/// calculates true ionic strength
long int TDebyeHueckel::IonicStrength()
{
	long int j;
	double is, mt, mz, as;
	is = 0.0; mt = 0.0; mz = 0.0; as = 0.0;

	// calculate ionic strength and total molalities
	// needs to be modified when nonaqueous solvents are present
	for (j=0; j<(NComp-1); j++)
	{
		is += 0.5*m[j]*z[j]*z[j];
		mt += m[j];
		if ( z[j] )
		{
			mz += m[j];
			as += m[j]*an[j];
		}
	}

	// assignments
	IS = is;
	molT = mt;
	ao = as/mz;

	return 0;
}





//=============================================================================================
// Extended Debye-Hueckel (EDH) model for aqueous solutions, Karpov version
// References: Karpov et al. (1997); Helgeson et al. (1981); Oelkers and Helgeson (1990);
// Pokrovskii and Helgeson (1995; 1997a; 1997b)
// (c) TW May 2009
//=============================================================================================


// Generic constructor for the TKarpov class
TKarpov::TKarpov( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW ):
                TSolMod( sd )
{
    alloc_internal();
    m = arM;
    z = arZ;
    RhoW = dW;
    EpsW = eW;
    ac = aIPc[1];   // common ion size parameter
    bc = aIPc[0];   // common b_gamma
    flagNeut = (long int)aIPc[2];   // 0: unity, 1: calculated
    flagH2O = (long int)aIPc[3];   // 0: unity, 1: calculated
    flagElect = (long int)aIPc[4];  // 0: constant, 1: NaCl, 2: KCl, 3: NaOH, 4: KOH
}


TKarpov::~TKarpov()
{
	free_internal();
}


void TKarpov::alloc_internal()
{
	LnG = new double [NComp];
	dLnGdT = new double [NComp];
	d2LnGdT2 = new double [NComp];
	dLnGdP = new double [NComp];
	an = new double [NComp];
	bg = new double [NComp];
}


void TKarpov::free_internal()
{
	delete[]LnG;
	delete[]dLnGdT;
	delete[]d2LnGdT2;
	delete[]dLnGdP;
	delete[]an;
	delete[]bg;
}


/// Calculates T,P corrected parameters
long int TKarpov::PTparam()
{
	long int j;
	double alp, bet, dal, rho, eps, dedt, d2edt2, dedp;

	// read and copy individual an and bg
	for (j=0; j<NComp; j++)
	{
		an[j] = aDCc[NP_DC*j];
		bg[j] = aDCc[NP_DC*j+1];  // individual bg (not used)
	}

	// read and convert rho and eps
	rho = RhoW[0];
	alp = - 1./rho*RhoW[1];
	dal = pow(alp,2.) - 1./rho*RhoW[2];
	bet = 1./rho*RhoW[3];
	eps = EpsW[0];
	dedt = 1./eps*EpsW[1];
	d2edt2 = - 1./pow(eps,2.)*pow(EpsW[1],2.) + 1./eps*EpsW[2];
	dedp = 1./eps*EpsW[3];

	// calculate A and B terms of Debye-Huckel equation (and derivatives)
	A = (1.82483e6)*sqrt(rho) / pow(Tk*eps,1.5);
	dAdT = - 3./2.*A*( dedt + 1./Tk + alp/3. );
	d2AdT2 = 1./A*pow(dAdT,2.) - 3./2.*A*( d2edt2 - 1/pow(Tk,2.) + 1/3.*dal );
	dAdP = 1./2.*A*( bet - 3.*dedp);
	B = (50.2916)*sqrt(rho) / sqrt(Tk*eps);
	dBdT = - 1./2.*B*( dedt +1./Tk + alp );
	d2BdT2 = 1./B*pow(dBdT,2.) - 1./2.*B*( d2edt2 - 1./pow(Tk,2.) + dal );
	dBdP = 1./2.*B*( bet - dedp );

 	ErrorIf( fabs(A) < 1e-9 || fabs(B) < 1e-9, "Karpov EDH model",
 			"Error: DH parameter A or B was not calculated - no values of RhoW and EpsW !" );

	// b_gamma constant
	if ( flagElect == 0)
	{
		bgam = bc;
		dbgdT = 0.0;
		d2bgdT2 = 0.0;
		dbgdP = 0.0;
		ao = ac;  // constant
	}

	// b_gamma TP-dependent
	else
	{
		Gfunction();
		BgammaTP();
		ao = ac;  // constant
	}

	return 0;
}


/// Calculates activity coefficients
long int TKarpov::MixMod()
{
	long int j, k, w;
	double sqI, Z2, lgGam, lnGam, Nw, Lgam, lnxw, xw, Lam, sig,
			Phi, Phit, psi, zc, za, lnaw, lg_to_ln;
	zc = 1.; za = 1.; psi = 1.; lg_to_ln = 2.302585093;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		lgGam = 0.0;
		lnGam = 0.0;

		// charged species (individual ion size parameters)
		if ( z[j] )
		{
			lgGam = 0.0;
			Z2 = z[j]*z[j];
			lgGam = ( - A * sqI * Z2 ) / ( 1. + B * an[j] * sqI ) + bgam * IS ;
			lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				lgGam = 0.0;
				if ( flagNeut == 1 )
					lgGam = bgam * IS;
				else
					lgGam = 0.0;
				lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
				continue;
			}

			// water solvent
			else
			{
				lgGam = 0.0;
				lnGam = 0.0;
				if ( flagH2O == 1 )
				{
					Phit = 0.0;
					// Phi corrected using eq. (190) from Helgeson et al. (1981)
					Lam = 1. + ao*B*sqI;
					sig = 3./(pow(ao,3.)*pow(B,3.)*pow(IS,(3./2.))) * (Lam-1./Lam-2*log(Lam));
					// Phi = -2.3025851*(A*sqI*SigTerm/3. + Lgam/(0.0180153*2.*IS) - bgam*IS/2.);
					// Phi = - log(10.) * (molZ/molT) * ( (zc*za*A*sqI*SigTerm)/3. + (psi*Lgam)/(0.0180153*2.*IS) - bgam*IS/2. );

					for (k=0; k<(NComp-1); k++)
					{
						if ( (z[k] == 0) && (flagNeut == 0) )
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) );
						else
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) - bgam*IS/2. );
					}

					Phi = Phit/molT;
					lnaw = - Phi*molT/Nw;
					lnGam = lnaw - lnxw;
				}
				else
					lnGam = 0.0;
				lnGamma[j] = lnGam;
			}
		}
	} // j

	return 0;
}


/// calculates excess properties
long int TKarpov::ExcessProp( double *Zex )
{
	long int j, k, w;
	double sqI, Z2, Nw, Lgam, lnxw, xw, Lam, sig, Phi, dPhidT, d2PhidT2, dPhidP,
			Phit, dPhitdT, d2PhitdT2, dPhitdP, lnaw, lg_to_ln, g, dgt, d2gt, dgp;
	double U, V, dUdT, dVdT, d2UdT2, d2VdT2, dUdP, dVdP, U1, U2, U3, V1, V2, V3,
			dU1dT, dU2dT, dU3dT, dV1dT, dV2dT, dV3dT, d2U1dT2, d2U2dT2, d2U3dT2,
			d2V1dT2, d2V2dT2, d2V3dT2, dU1dP, dU2dP, dU3dP, dV1dP, dV2dP, dV3dP,
			L, dLdT, d2LdT2, dLdP, Z, dZdT, d2ZdT2, dZdP;
	lg_to_ln = 2.302585093;
	g = 0.; dgt = 0.; d2gt = 0.; dgp = 0.;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		// charged species
		if ( z[j] )
		{
			Z2 = z[j]*z[j];
			U = - (Z2*A) * sqI;
			dUdT = - (Z2*dAdT) * sqI;
			d2UdT2 = - (Z2*d2AdT2) * sqI;
			dUdP = - (Z2*dAdP) * sqI;
			V = 1. + (an[j]*B) * sqI;
			dVdT = ( an[j]*dBdT ) * sqI;
			d2VdT2 = ( an[j]*d2BdT2 ) * sqI;
			dVdP = ( an[j]*dBdP ) * sqI;
			LnG[j] = ( U/V + bgam*IS ) * lg_to_ln;
			dLnGdT[j] = ( (dUdT*V - U*dVdT)/pow(V,2.) + dbgdT*IS ) * lg_to_ln;
			d2LnGdT2[j] = ( (d2UdT2*V + dUdT*dVdT)/pow(V,2.) - (dUdT*V)*(2.*dVdT)/pow(V,3.)
				- (dUdT*dVdT + U*d2VdT2)/pow(V,2.) + (U*dVdT)*(2.*dVdT)/pow(V,3.)
				+ d2bgdT2*IS ) * lg_to_ln;
			dLnGdP[j] = ( (dUdP*V - U*dVdP)/pow(V,2.) + dbgdP*IS ) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				if ( flagNeut == 1 )
				{
					LnG[j] = ( bgam*IS ) * lg_to_ln;
					dLnGdT[j] = ( dbgdT*IS ) * lg_to_ln;
					d2LnGdT2[j] = ( d2bgdT2*IS ) * lg_to_ln;
					dLnGdP[j] = ( dbgdP*IS ) * lg_to_ln;
				}

				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
				continue;
			}

			// water solvent
			else
			{
				// H2O activity coefficient calculated
				if ( flagH2O == 1 )
				{
					Phit = 0.; dPhitdT = 0.; d2PhitdT2 = 0.; dPhitdP = 0.;

					// derivatives of lambda and sigma terms
					L = 1. + (ao*B) * sqI;
					dLdT = ( ao*dBdT ) * sqI;
					d2LdT2 = ( ao*d2BdT2 ) * sqI;
					dLdP = ( ao*dBdP ) * sqI;

					U1 = (A*L);
					dU1dT = (dAdT*L + A*dLdT);
					d2U1dT2 = ( d2AdT2*L + 2.*dAdT*dLdT + A*d2LdT2 );
					dU1dP = ( dAdP*L + A*dLdP );
					V1 = pow(ao,3.)*pow(B,3.) * IS;
					dV1dT = ( 3.*pow(ao,3.)*pow(B,2.)*dBdT ) * IS;
					d2V1dT2 = ( 6.*pow(ao,3.)*B*pow(dBdT,2.) + 3.*pow(ao,3.)*pow(B,2.)*d2BdT2 ) * IS;
					dV1dP = ( 3.*pow(ao,3.)*pow(B,2.)*dBdP ) * IS;

					U2 = A;
					dU2dT = dAdT;
					d2U2dT2 = d2AdT2;
					dU2dP = dAdP;
					V2 = pow(ao,3.)*pow(B,3.)*L * IS;
					dV2dT = ( 3.*pow(ao,3.)*pow(B,2.)*dBdT*L + pow(ao,3.)*pow(B,3.)*dLdT ) * IS;
					d2V2dT2 = ( 6.*pow(ao,3.)*B*pow(dBdT,2.)*L + 3.*pow(ao,3.)*pow(B,2.)*d2BdT2*L
								+ 6.*pow(ao,3.)*pow(B,2.)*dBdT*dLdT + pow(ao,3.)*pow(B,3.)*d2LdT2 ) * IS;
					dV2dP = ( 3.*pow(ao,3.)*pow(B,2.)*dBdP*L + pow(ao,3.)*pow(B,3.)*dLdP ) * IS;

					U3 = 2.*( A*log(L) );
					dU3dT = 2.*( dAdT*log(L) + A*(1./L)*dLdT );
					d2U3dT2 = 2.*( d2AdT2*log(L) + 2.*dAdT*(1./L)*dLdT
								- A*(1./pow(L,2.))*pow(dLdT,2.) + A*(1./L)*d2LdT2 );
					dU3dP = 2.*( dAdP*log(L) + A*(1./L)*dLdP );
					V3 = pow(ao,3.)*pow(B,3.) * IS;
					dV3dT = ( 3.*pow(ao,3.)*pow(B,2.)*dBdT ) * IS;
					d2V3dT2 = ( 6.*pow(ao,3.)*B*pow(dBdT,2.) + 3.*pow(ao,3.)*pow(B,2.)*d2BdT2 ) * IS;
					dV3dP = ( 3.*pow(ao,3.)*pow(B,2.)*dBdP ) * IS;

					Z = U1/V1 - U2/V2 - U3/V3;
					dZdT = (dU1dT*V1 - U1*dV1dT)/pow(V1,2.) - (dU2dT*V2 - U2*dV2dT)/pow(V2,2.)
								- (dU3dT*V3 - U3*dV3dT)/pow(V3,2.);
					d2ZdT2 = (d2U1dT2*V1 + dU1dT*dV1dT)/pow(V1,2.) - (dU1dT*V1)*(2.*dV1dT)/pow(V1,3.)
								- (dU1dT*dV1dT + U1*d2V1dT2)/pow(V1,2.) + (U1*dV1dT)*(2.*dV1dT)/pow(V1,3.)
								- (d2U2dT2*V2 + dU2dT*dV2dT)/pow(V2,2.) + (dU2dT*V2)*(2.*dV2dT)/pow(V2,3.)
								+ (dU2dT*dV2dT + U2*d2V2dT2)/pow(V2,2.) - (U2*dV2dT)*(2.*dV2dT)/pow(V2,3.)
								- (d2U3dT2*V3 + dU3dT*dV3dT)/pow(V3,2.) + (dU3dT*V3)*(2.*dV3dT)/pow(V3,3.)
								+ (dU3dT*dV3dT + U3*d2V3dT2)/pow(V3,2.) - (U3*dV3dT)*(2.*dV3dT)/pow(V3,3.);
					dZdP = (dU1dP*V1 - U1*dV1dP)/pow(V1,2.) - (dU2dP*V2 - U2*dV2dP)/pow(V2,2.)
								- (dU3dP*V3 - U3*dV3dP)/pow(V3,2.);

					// increments to osmotic coefficient (and derivatives)
					Lam = 1. + (ao*B) * sqI;
					sig = 3./(pow(ao,3.)*pow(B,3.)*pow(IS,(3./2.))) * (Lam-1./Lam-2*log(Lam));

					for (k=0; k<(NComp-1); k++)
					{
						if ( (z[k] == 0) && (flagNeut == 0) )
						{
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) );
							dPhitdT  += - log(10.) * m[k] * ( pow(z[k],2.)*dZdT );
							d2PhitdT2 += - log(10.) * m[k] * ( pow(z[k],2.)*d2ZdT2 );
							dPhitdP += - log(10.) * m[k] * ( pow(z[k],2.)*dZdP );

						}

						else
						{
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) - bgam*IS/2. );
							dPhitdT  += - log(10.) * m[k] * ( pow(z[k],2.)*dZdT - dbgdT*IS/2. );
							d2PhitdT2 += - log(10.) * m[k] * ( pow(z[k],2.)*d2ZdT2 - d2bgdT2*IS/2. );
							dPhitdP += - log(10.) * m[k] * ( pow(z[k],2.)*dZdP - dbgdP*IS/2. );
						}
					}

					Phi = Phit/molT;
					dPhidT = dPhitdT/molT;
					d2PhidT2 = d2PhitdT2/molT;
					dPhidP = dPhitdP/molT;

					// activity coefficient (and derivatives)
					lnaw = - Phi*molT/Nw;
					LnG[j] = lnaw - lnxw;
					dLnGdT[j] = - (molT/Nw) * dPhidT;
					d2LnGdT2[j] = - (molT/Nw) * d2PhidT2;
					dLnGdP[j] = - (molT/Nw) * dPhidP;
				}

				// H2O activity coefficient 0
				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
			}
		}

		g += x[j]*LnG[j];
		dgt += x[j]*dLnGdT[j];
		d2gt += x[j]*d2LnGdT2[j];
		dgp += x[j]*dLnGdP[j];

	} // j

	// increment thermodynamic properties
	Gex = (R_CONST*Tk) * g;
	Hex = - R_CONST*pow(Tk,2.) * dgt;
	// Sex = - R_CONST * ( g + Tk*dgt );
	Sex = (Hex-Gex)/Tk;
	CPex = - R_CONST * ( 2.*Tk*dgt + pow(Tk,2.)*d2gt );
	Vex = (R_CONST*Tk) * dgp;
	Aex = Gex - Vex*Pbar;
	Uex = Hex - Vex*Pbar;

	// assigments (excess properties)
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
long int TKarpov::IdealProp( double *Zid )
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


/// calculates true ionic strength
long int TKarpov::IonicStrength()
{
	long int j;
	double is, mt, mz, as;
	is = 0.0; mt = 0.0; mz = 0.0; as = 0.0;

	// calculate ionic strength and total molalities
	// needs to be modified when nonaqueous solvents are present
	for (j=0; j<(NComp-1); j++)
	{
		is += 0.5*m[j]*z[j]*z[j];
		mt += m[j];
		if ( z[j] )
		{
			mz += m[j];
			as += m[j]*an[j];
		}
	}

	// conversions and assignments
	ao = as/mz;
	IS = is;
	molT = mt;
	molZ = mz;

	return 0;
}


/// calculates TP dependence of b_gamma (and derivatives)
long int TKarpov::BgammaTP()
{
	// ni: stoichiometric number of moles of ions in one mole of electrolyte
	// rc, ra: radius of cation and anion, respectively at 298 K/1 bar
	// units are cal, kg, K, mol, bar
	double ni, nc, na, zc, za, rc, ra, a1, a2, a3, a4, a5, c1, c2, omg, bg, bs, bh, rec, rea,
			omgpt, domdt, d2omdt2, domdp, nbg, nbv, nbj, nbh;
	double eps, eta, xborn, yborn, qborn, X1, X2;

	// set parameters
	eps = EpsW[0];
	yborn = 1./pow(EpsW[0],2.)*EpsW[1];
	xborn = 1./pow(EpsW[0],2.) * ( EpsW[2] - 2./EpsW[0]*pow(EpsW[1],2.) );
	qborn = 1./pow(EpsW[0],2.)*EpsW[3];
	eta = (1.66027e5);

	switch ( flagElect )
	{
		case 1:  // NaCl
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.030056; a2 = -202.55; a3 = -2.9092; a4 = 20302;
			a5 = -0.206; c1 = -1.50; c2 = 53300.; omg = 178650.;
			bg = -174.623; bs = 2.164; rc = 0.97; ra = 1.81;
			break;
		case 2:  // KCl
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.0172; a2 = -115.36; a3 = -1.1857; a4 = 13854.2;
			a5 = -0.262; c1 = -2.53; c2 = 38628.4; omg = 164870.;
			bg = -70.0; bs = 1.727; rc = 1.33; ra = 1.81;
			break;
		case 3:  // NaOH
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.030056; a2 = -202.55; a3 = -2.9092; a4 = 20302;
			a5 = -0.206; c1 = -1.50; c2 = 53300.; omg = 205520.;
			bg = -267.4; bs = 1.836; rc = 0.97; ra = 1.40;
			break;
		case 4:  // KOH
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.0172; a2 = -115.36; a3 = -1.1857; a4 = 13854.2;
			a5 = -0.262; c1 = -2.53; c2 = 38628.4; omg = 191730.;
			bg = -335.7; bs = 1.26; rc = 1.33; ra = 1.40;
			break;
		default:  // wrong mode
			return -1;
	}

	// calculation part, extended 06.06.2009 (TW)
	bh = bg + (298.15)*bs;
	rec = rc + fabs(zc)*(0.94+Gf);
	rea = ra + fabs(za)*Gf;
	X1 = - eta*nc*( fabs(pow(zc,3.))/pow(rec,2.) - zc/pow((3.082+Gf),2.) )
			- eta*na*( fabs(pow(za,3.))/pow(rea,2.) - za/pow((3.082+Gf),2.) );
	X2 = 2.*eta*nc*( fabs(pow(zc,4.))/pow(rec,3.) - zc/pow((3.082+Gf),3.) )
			+ 2.*eta*nc * ( fabs(pow(zc,4.))/pow(rec,3.) - zc/pow((3.082+Gf),3.) );
	omgpt = eta*( nc*pow(zc,2.)/rec + na*pow(za,2.)/rea );
	// omgpt = (1.66027e5)*(1./(0.94+rc+Gf)+1./(ra+Gf));
	domdt = X1*dGfdT;
	d2omdt2 = X2*pow(dGfdT,2.) + X1*d2GfdT2;
	domdp = X1*dGfdP;
	nbg = - ni*bg/2. + ni*bs*(Tk-298.15)/2. - c1*(Tk*log(Tk/298.15)-Tk+298.15)
				+ a1*(Pbar-1.) + a2*log((2600.+Pbar)/(2600.+1.))
				- c2*((1./(Tk-228.)-1./(298.15-228.))*(228.-Tk)/228.-Tk/(228.*228.)
				* log((298.15*(Tk-228.))/(Tk*(298.15-228.))))
				+ 1./(Tk-228.)*(a3*(Pbar-1.) + a4*log((2600.+Pbar)/(2600.+1.)))
				+ a5*(omgpt*(1./eps-1.)-omg*(1./(78.24513795)-1.)+(-5.798650444e-5)*omg*(Tk-298.15));
	nbh = - ni/2.*bh + c1*(Tk-298.15) - c2*(1./(Tk-228.)-1./(298.15-228.))
				+ a1*(Pbar-1.) + a2*log((2600.+Pbar)/(2600.+1.))
				+ ((2.*Tk-228.)/pow((Tk-228.),2.))*(a3*(Pbar-1.)+a4*log((2600.+Pbar)/(2600.+1.)))
				+ a5*( omgpt*(1./eps-1.) + omgpt*Tk*yborn - Tk*(1./eps-1.)*domdt
				- omg*(1./(78.24513795)-1.) - omg*(298.15)*(-5.798650444e-5) );
	nbj = c1 + c2/pow((Tk-228.),2.)  - ( (2.*Tk)/pow((Tk-228.),3.) )
				* ( a3*(Pbar-1.) + a4*log((2600.+Pbar)/(2600.+1.)) )
				+ a5*( omgpt*Tk*xborn + 2.*Tk*yborn*domdt - Tk*(1./eps-1.)*d2omdt2 );
	nbv = a1 + a2/(2600.+Pbar) + a3/(Tk-228.) + a4/((2600.+Pbar)*(Tk-228.))
				+ a5*(-omgpt*qborn + (1./eps-1.)*domdp);
	// b_gamma = nbg/(2.*log(10.)*1.98721*Tk);
	// bgam = b_gamma;
	bgam = nbg/(2.*log(10.)*(1.98721)*Tk)*2./ni;
	dbgdT = - nbh/(2.*log(10.)*(1.98721)*pow(Tk,2.))*2./ni;
	d2bgdT2 = - nbj/( 2.*log(10.)*(1.98721)*pow(Tk,2.))*2./ni - 2./Tk*dbgdT;
	dbgdP = nbv/(2.*log(10.)*(1.98721)*Tk)*2./ni;

	return 0;
}


/// wrapper for g-function
long int TKarpov::Gfunction()
{
	double T, P, D, beta, alpha, daldT, g, dgdP, dgdT, d2gdT2;
	double TMAX = 1000., PMAX = 5000., TOL = 1.0e-4;

	// convert parameters
	T = Tk - 273.15;
	P = Pbar;
	D = RhoW[0];
	alpha = - RhoW[1]/RhoW[0];
	daldT = pow(alpha, 2.) - RhoW[2]/RhoW[0];
	beta = RhoW[3]/RhoW[0];

	// initialize g and derivatives to zero
	g = 0.0;
	dgdP = 0.0;
	dgdT = 0.0;
	d2gdT2 = 0.0;

	if ((T > TMAX+TOL) || (P > PMAX+TOL))
		return -1;
	else
		GShok2( T, P, D, beta, alpha, daldT, g, dgdP, dgdT, d2gdT2 );

	// assignments
	Gf = g;
	dGfdT = dgdT;
	d2GfdT2 = d2gdT2;
	dGfdP = dgdP;

	return 0;
}


/// calculates g-function and derivatives
long int TKarpov::GShok2( double T, double P, double D, double beta,
		double alpha, double daldT, double &g, double &dgdP, double &dgdT, double &d2gdT2 )
{
	double a, b, dgdD, /*dgdD2,*/ dadT, dadTT, dbdT, dbdTT, dDdT, dDdP,
		dDdTT, Db, dDbdT, dDbdTT, ft, dftdT, dftdTT, fp, dfpdP,
		f, dfdP, dfdT, d2fdT2, tempy;
	double C[6]  = {-0.2037662e+01,  0.5747000e-02, -0.6557892e-05,
			0.6107361e+01, -0.1074377e-01,  0.1268348e-04 };
	double cC[3] = { 0.3666666e+02, -0.1504956e-09,  0.5017997e-13 };

	if ( D >= 1.3 )
		return -1;
	// Sveta 19/02/2000 1-D < 0 pri D>1 => pow(-number, b) = -NaN0000 error
	double pw = fabs(1.0 - D); // insert Sveta 19/02/2000

	// calculation part
	a = C[0] + C[1]*T + C[2]*pow(T,2.);
	b = C[3] + C[4]*T + C[5]*pow(T,2.);
	g = a * pow(pw, b);

	dgdD = - a*b* pow(pw, (b - 1.0));
	// dgdD2 = a * b * (b - 1.0) * pow((1.0 - D),(b - 2.0));

	dadT = C[1] + 2.0*C[2]*T;
	dadTT = 2.0*C[2];
	dbdT = C[4] + 2.0*C[5]*T;
	dbdTT = 2.0*C[5];
	dDdT = - D * alpha;
	dDdP = D * beta;
	dDdTT = - D * (daldT - pow(alpha,2.));

	// Db = pow((1.0 - D), b);  Fixed by DAK 01.11.00
	Db = pow(pw , b);
	dDbdT = - b*pow(pw,(b - 1.0))*dDdT + log(pw)*Db*dbdT;

	dDbdTT = - (b*pow(pw,(b-1.0)) *dDdTT + pow(pw,(b - 1.0)) * dDdT * dbdT
				+ b * dDdT * (-(b - 1.0) * pow(pw,(b - 2.0)) * dDdT
				+ log(pw) * pow(pw,(b - 1.0)) * dbdT))
				+ log(pw) * pow(pw,b) * dbdTT
				- pow(pw,b) * dbdT * dDdT / (1.0 - D)
				+ log(pw) * dbdT * dDbdT;

	dgdP = dgdD * dDdP;
	dgdT = a * dDbdT + Db * dadT;
	d2gdT2 = a * dDbdTT + 2.0 * dDbdT * dadT + Db * dadTT;

	if((T < 155.0) || (P > 1000.0) || (T > 355.0))
		return 0;

	tempy = ((T - 155.0) / 300.0);
	ft = pow(tempy,4.8) + cC[0] * pow(tempy,16.);

	dftdT = 4.8 / 300.0 * pow(tempy, 3.8) + 16.0 / 300.0 * cC[0] * pow(tempy, 15.);

	dftdTT = 3.8 * 4.8 / (300.0 * 300.0) * pow(tempy, 2.8)
		+ 15.0 * 16.0 / (300.0*300.0) * cC[0] * pow(tempy,14.);

	fp = cC[1] * pow((1000.0 - P),3.) + cC[2] * pow((1000.0 - P),4.);

	dfpdP = -3.0 * cC[1] * pow((1000.0 - P),2.) - 4.0 * cC[2] * pow((1000.0 - P),3.);

	f = ft * fp;
	dfdP = ft * dfpdP;
	dfdT = fp * dftdT;
	d2fdT2 = fp * dftdTT;

	g -= f;
	dgdP -= dfdP;
	dgdT -= dfdT;
	d2gdT2 -= d2fdT2;

	return 0;
}





//=============================================================================================
// Extended Debye-Hueckel (EDH) model for aqueous solutions, Shvarov version
// References: Shvarov (2007); Oelkers and Helgeson (1990);
// Pokrovskii and Helgeson (1995; 1997a; 1997b)
// (c) TW June 2009
//=============================================================================================


// Generic constructor for the TShvarov class
TShvarov::TShvarov( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW ):
                TSolMod( sd )
{
    alloc_internal();
    m = arM;
    z = arZ;
    RhoW = dW;
    EpsW = eW;
    ac = aIPc[1];   // common ion size parameter
    bc = aIPc[0];   // common b_gamma
    flagNeut = (long int)aIPc[2];   // 0: unity, 1: calculated
    flagH2O = (long int)aIPc[3];   // 0: unity, 1: calculated
    flagElect = (long int)aIPc[4];  // 0: constant, 1: NaCl, 2: KCl, 3: NaOH, 4: KOH
}


TShvarov::~TShvarov()
{
	free_internal();
}


void TShvarov::alloc_internal()
{
	LnG = new double [NComp];
	dLnGdT = new double [NComp];
	d2LnGdT2 = new double [NComp];
	dLnGdP = new double [NComp];
	bj = new double [NComp];
}


void TShvarov::free_internal()
{
	delete[]LnG;
	delete[]dLnGdT;
	delete[]d2LnGdT2;
	delete[]dLnGdP;
	delete[]bj;
}


long int TShvarov::PTparam()
{
	long int j;
	double alp, bet, dal, rho, eps, dedt, d2edt2, dedp;

	// read and copy individual bj parameters
	for (j=0; j<NComp; j++)
	{
		if ( (aDCc[NP_DC*j+1]) > 0 )
			bj[j] = aDCc[NP_DC*j+1];  // individual bj
		else
			bj[j] = 1.;
	}

	// read and convert rho and eps
	rho = RhoW[0];
	alp = - 1./rho*RhoW[1];
	dal = pow(alp,2.) - 1./rho*RhoW[2];
	bet = 1./rho*RhoW[3];
	eps = EpsW[0];
	dedt = 1./eps*EpsW[1];
	d2edt2 = - 1./pow(eps,2.)*pow(EpsW[1],2.) + 1./eps*EpsW[2];  // corrected 23.05.2009 (TW)
	dedp = 1./eps*EpsW[3];

	// calculate A and B terms of Debye-Huckel equation (and derivatives)
	A = (1.82483e6)*sqrt(rho) / pow(Tk*eps,1.5);
	dAdT = - 3./2.*A*( dedt + 1./Tk + alp/3. );
	d2AdT2 = 1./A*pow(dAdT,2.) - 3./2.*A*( d2edt2 - 1/pow(Tk,2.) + 1/3.*dal );
	dAdP = 1./2.*A*( bet - 3.*dedp);
	B = (50.2916)*sqrt(rho) / sqrt(Tk*eps);
	dBdT = - 1./2.*B*( dedt +1./Tk + alp );
	d2BdT2 = 1./B*pow(dBdT,2.) - 1./2.*B*( d2edt2 - 1./pow(Tk,2.) + dal );
	dBdP = 1./2.*B*( bet - dedp );

 	ErrorIf( fabs(A) < 1e-9 || fabs(B) < 1e-9, "Shvarov EDH model",
 			"Error: DH parameter A or B was not calculated - no values of RhoW and EpsW !" );

	// b_gamma constant
	if ( flagElect == 0)
	{
		bgam = bc;
		dbgdT = 0.0;
		d2bgdT2 = 0.0;
		dbgdP = 0.0;
		ao = ac;
		daodT = 0.0;
		d2aodT2 = 0.0;
		daodP = 0.0;
	}

	// b_gamma TP-dependent
	else
	{
		Gfunction();
		BgammaTP();
		IonsizeTP();
		// ao = ac;
		// daodT = 0.0;
		// d2aodT2 = 0.0;
		// daodP = 0.0;
	}

	return 0;
}


long int TShvarov::MixMod()
{
	long int j, k, w;
	double sqI, Z2, lgGam, lnGam, Nw, xw, lnxw, Lgam,
				msum, C, lg_to_ln;
	lg_to_ln = 2.302585093;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molaities (molT and molZ)
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);
	C = (0.5*bgam);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		lgGam = 0.0;
		lnGam = 0.0;
		msum = 0.0;

		// calculate (bj*mj) sum
		for (k=0; k<(NComp-1); k++)
		{
			msum += bj[k]*m[k];
		}

		// charged species
		if ( z[j] )
		{
			lgGam = 0.0;
			Z2 = z[j]*z[j];
			lgGam = - (A*Z2*sqI) / (1.+ao*B*sqI) + C*bj[j]*msum;
			lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				lgGam = 0.0;
				if ( flagNeut == 1 )
					lgGam = C*bj[j]*msum;
				else
					lgGam = 0.0;
				lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
				continue;
			}

			// water solvent
			else
			{
				lgGam = 0.0;
				lnGam = 0.0;
				if ( flagH2O == 1 )
				{
					lgGam = - A/(ao*B) * (2./Nw) * ( IS/(1.+ao*B*sqI) - 2.*sqI/(ao*B) + 2./pow((ao*B),2.) * log(1.+ao*B*sqI) )
								 - C*pow(msum,2.)/(2.*Nw);
				}
				else
					lgGam = 0.0;
				lnGamma[j] = (lgGam) * lg_to_ln;
			}
		}
	}  // j

	return 0;
}


long int TShvarov::ExcessProp( double *Zex )
{
	long int j, k, w;
	double sqI, Z2, Nw, xw, msum, C, dCdT, d2CdT2, dCdP, lg_to_ln,
				g, dgt, d2gt, dgp;
	double U, V, dUdT, dVdT, d2UdT2, d2VdT2, dUdP, dVdP, U1, U2, U3, V1, V2, V3,
				dU1dT, dU2dT, dU3dT, dV1dT, dV2dT, dV3dT, d2U1dT2, d2U2dT2, d2U3dT2,
				d2V1dT2, d2V2dT2, d2V3dT2, dU1dP, dU2dP, dU3dP, dV1dP, dV2dP, dV3dP,
				X, dXdT, d2XdT2, dXdP, Z, dZdT, d2ZdT2, dZdP;
	lg_to_ln = 2.302585093;
	g = 0.; dgt = 0.; d2gt = 0.; dgp = 0.;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	sqI = sqrt(IS);
	C = (0.5*bgam);
	dCdT = (0.5*dbgdT);
	d2CdT2 = (0.5*d2bgdT2);
	dCdP = (0.5*dbgdP);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		msum = 0.0;

		// calculate bj*mj sum
		for (k=0; k<(NComp-1); k++)
		{
			msum += bj[k]*m[k];
		}

		// charged species
		if ( z[j] )
		{
			Z2 = z[j]*z[j];
			U = - (Z2*A) * sqI;
			dUdT = - (Z2*dAdT) * sqI;
			d2UdT2 = - (Z2*d2AdT2) * sqI;
			dUdP = - (Z2*dAdP) * sqI;
			V = 1. + (ao*B) * sqI;
			dVdT = ( daodT*B + ao*dBdT ) * sqI;
			d2VdT2 = ( d2aodT2*B + 2.*daodT*dBdT + ao*d2BdT2 ) * sqI;
			dVdP = ( daodP*B + ao*dBdP ) * sqI;
			LnG[j] = ( U/V + C*bj[j]*msum ) * lg_to_ln;
			dLnGdT[j] = ( (dUdT*V - U*dVdT)/pow(V,2.) + dCdT*bj[j]*msum ) * lg_to_ln;
			d2LnGdT2[j] = ( (d2UdT2*V + dUdT*dVdT)/pow(V,2.) - (dUdT*V)*(2.*dVdT)/pow(V,3.)
				- (dUdT*dVdT + U*d2VdT2)/pow(V,2.) + (U*dVdT)*(2.*dVdT)/pow(V,3.)
				+ d2CdT2*bj[j]*msum ) * lg_to_ln;
			dLnGdP[j] = ( (dUdP*V - U*dVdP)/pow(V,2.) + dCdP*bj[j]*msum ) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				if ( flagNeut == 1 )
				{
					LnG[j] = ( dCdT*bj[j]*msum ) * lg_to_ln;
					dLnGdT[j] = ( dCdT*bj[j]*msum ) * lg_to_ln;
					d2LnGdT2[j] = ( d2CdT2*bj[j]*msum ) * lg_to_ln;
					dLnGdP[j] = ( dCdP*bj[j]*msum ) * lg_to_ln;
				}

				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
				continue;
			}

			// water solvent
			else
			{
				// H2O activity coefficient calculated
				if ( flagH2O == 1 )
				{
					// derivatives of lambda and sigma terms
					U1 = A * IS;
					dU1dT = dAdT * IS;
					d2U1dT2 = d2AdT2 * IS;
					dU1dP = dAdP * IS;
					V1 = (ao*B) + (pow(ao,2.)*pow(B,2.)) * sqI;
					dV1dT = ( daodT*B + ao*dBdT ) + 2.*( ao*daodT*pow(B,2.) + pow(ao,2.)*B*dBdT ) * sqI;
					d2V1dT2 = ( d2aodT2*B + 2.*daodT*dBdT + ao*d2BdT2 )
								+ 2. * ( pow(daodT,2.)*pow(B,2.) + ao*d2aodT2*pow(B,2.) + 4.*ao*daodT*B*dBdT
								+ pow(ao,2.)*pow(dBdT,2.) + pow(ao,2.)*B*d2BdT2 ) * sqI;
					dV1dP = ( daodP*B + ao*dBdP ) + 2.*( ao*daodP*pow(B,2.) + pow(ao,2.)*B*dBdP ) * sqI;

					U2 = (2.*A) * sqI;
					dU2dT = (2*dAdT) * sqI;
					d2U2dT2 = (2.*d2AdT2) * sqI;
					dU2dP = (2.*dAdP) * sqI;
					V2 = pow(ao,2.)*pow(B,2.);
					dV2dT = 2.*( ao*daodT*pow(B,2.) + pow(ao,2.)*B*dBdT );
					d2V2dT2 = 2.*( pow(daodT,2.)*pow(B,2.) + ao*d2aodT2*pow(B,2.) + 4.*ao*daodT*B*dBdT
								+ pow(ao,2.)*pow(dBdT,2.) + pow(ao,2.)*B*d2BdT2 );
					dV2dP = 2.*( ao*daodP*pow(B,2.) + pow(ao,2.)*B*dBdP );

					U3 = 2.*A*log(1.+ao*B*sqI);
					X = log(1.+ao*B*sqI);
					dXdT = pow((1.+ao*B*sqI),-1.) * ( daodT*B + ao*dBdT ) * sqI;
					d2XdT2 = - pow((1.+ao*B*sqI),-2.) * pow(( daodT*B + ao*dBdT ),2.) * IS
								+ pow((1.+ao*B*sqI),-1.) * ( d2aodT2*B + 2.*daodT*dBdT + ao*d2BdT2 ) * sqI;
					dXdP = pow((1.+ao*B*sqI),-1.) * ( daodP*B + ao*dBdP ) * sqI;
					dU3dT = 2.*( dAdT*X + A*dXdT );
					d2U3dT2 = 2.*( d2AdT2*X + 2.*dAdT*dXdT + A*d2XdT2 );
					dU3dP = 2.*( dAdP*X + A*dXdP );
					V3 = pow(ao,3.)*pow(B,3.);
					dV3dT = 3.*( pow(ao,2.)*daodT*pow(B,3.) + pow(ao,3.)*pow(B,2.)*dBdT );
					d2V3dT2 = 3.*( 2.*ao*pow(daodT,2.)*pow(B,3.) + pow(ao,2.)*d2aodT2*pow(B,3.)
								+ 6.*pow(ao,2.)*daodT*pow(B,2.)*dBdT + 2.*pow(ao,3.)*B*pow(dBdT,2.)
								+ pow(ao,3.)*pow(B,2.)*d2BdT2 );
					dV3dP = 3.*( pow(ao,2.)*daodP*pow(B,3.) + pow(ao,3.)*pow(B,2.)*dBdP );

					Z = U1/V1 - U2/V2 + U3/V3;
					dZdT = (dU1dT*V1 - U1*dV1dT)/pow(V1,2.) - (dU2dT*V2 - U2*dV2dT)/pow(V2,2.)
								+ (dU3dT*V3 - U3*dV3dT)/pow(V3,2.);
					d2ZdT2 = (d2U1dT2*V1 + dU1dT*dV1dT)/pow(V1,2.) - (dU1dT*V1)*(2.*dV1dT)/pow(V1,3.)
								- (dU1dT*dV1dT + U1*d2V1dT2)/pow(V1,2.) + (U1*dV1dT)*(2.*dV1dT)/pow(V1,3.)
								- (d2U2dT2*V2 + dU2dT*dV2dT)/pow(V2,2.) + (dU2dT*V2)*(2.*dV2dT)/pow(V2,3.)
								+ (dU2dT*dV2dT + U2*d2V2dT2)/pow(V2,2.) - (U2*dV2dT)*(2.*dV2dT)/pow(V2,3.)
								+ (d2U3dT2*V3 + dU3dT*dV3dT)/pow(V3,2.) - (dU3dT*V3)*(2.*dV3dT)/pow(V3,3.)
								- (dU3dT*dV3dT + U3*d2V3dT2)/pow(V3,2.) + (U3*dV3dT)*(2.*dV3dT)/pow(V3,3.);
					dZdP = (dU1dP*V1 - U1*dV1dP)/pow(V1,2.) - (dU2dP*V2 - U2*dV2dP)/pow(V2,2.)
								+ (dU3dP*V3 - U3*dV3dP)/pow(V3,2.);

					// activity coefficient (and derivatives)
					LnG[j] = ( - (2./Nw)*Z - C*pow(msum,2.)/(2.*Nw) ) * lg_to_ln;
					dLnGdT[j] = ( - (2./Nw)*dZdT - dCdT*pow(msum,2.)/(2.*Nw) ) * lg_to_ln;
					d2LnGdT2[j] = ( - (2./Nw)*d2ZdT2 - d2CdT2*pow(msum,2.)/(2.*Nw) ) * lg_to_ln;
					dLnGdP[j] = ( - (2./Nw)*dZdP - dCdP*pow(msum,2.)/(2.*Nw) ) * lg_to_ln;
				}

				// H2O activity coefficient 0
				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
			}

		}

		g += x[j]*LnG[j];
		dgt += x[j]*dLnGdT[j];
		d2gt += x[j]*d2LnGdT2[j];
		dgp += x[j]*dLnGdP[j];

	}  // j

	// increment thermodynamic properties
	Gex = (R_CONST*Tk) * g;
	Hex = - R_CONST*pow(Tk,2.) * dgt;
	// Sex = - R_CONST * ( g + Tk*dgt );
	Sex = (Hex-Gex)/Tk;
	CPex = - R_CONST * ( 2.*Tk*dgt + pow(Tk,2.)*d2gt );
	Vex = (R_CONST*Tk) * dgp;
	Aex = Gex - Vex*Pbar;
	Uex = Hex - Vex*Pbar;

	// assigments (excess properties)
	Zex[0] = Gex;
	Zex[1] = Hex;
	Zex[2] = Sex;
	Zex[3] = CPex;
	Zex[4] = Vex;
	Zex[5] = Aex;
	Zex[6] = Uex;

	return 0;
}


long int TShvarov::IdealProp( double *Zid )
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


long int TShvarov::IonicStrength()
{
	long int j;
	double is, mt;
	is = 0.0; mt = 0.0;

	// calculate ionic strength and total molalities
	// needs to be modified when nonaqueous solvents are present
	for (j=0; j<(NComp-1); j++)
	{
		is += 0.5*m[j]*z[j]*z[j];
		mt += m[j];
	}

	// assignments
	IS = is;
	molT = mt;

	return 0;
}


long int TShvarov::BgammaTP()
{
	// ni: stoichiometric number of moles of ions in one mole of electrolyte
	// rc, ra: radius of cation and anion, respectively at 298 K/1 bar
	// units are cal, kg, K, mol, bar
	double ni, nc, na, zc, za, rc, ra, a1, a2, a3, a4, a5, c1, c2, omg, bg, bs, bh;
	double rec, rea, omgpt, domdt, d2omdt2, domdp, nbg, nbv, nbj, nbh;
	double eps, eta, xborn, yborn, qborn, X1, X2;

	// set parameters
	eps = EpsW[0];
	yborn = 1./pow(EpsW[0],2.)*EpsW[1];
	xborn = 1./pow(EpsW[0],2.) * ( EpsW[2] - 2./EpsW[0]*pow(EpsW[1],2.) );
	qborn = 1./pow(EpsW[0],2.)*EpsW[3];
	eta = (1.66027e5);

	switch ( flagElect )
	{
		case 1:  // NaCl
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.030056; a2 = -202.55; a3 = -2.9092; a4 = 20302;
			a5 = -0.206; c1 = -1.50; c2 = 53300.; omg = 178650.;
			bg = -174.623; bs = 2.164; rc = 0.97; ra = 1.81;
			break;
		case 2:  // KCl
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.0172; a2 = -115.36; a3 = -1.1857; a4 = 13854.2;
			a5 = -0.262; c1 = -2.53; c2 = 38628.4; omg = 164870.;
			bg = -70.0; bs = 1.727; rc = 1.33; ra = 1.81;
			break;
		case 3:  // NaOH
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.030056; a2 = -202.55; a3 = -2.9092; a4 = 20302;
			a5 = -0.206; c1 = -1.50; c2 = 53300.; omg = 205520.;
			bg = -267.4; bs = 1.836; rc = 0.97; ra = 1.40;
			break;
		case 4:  // KOH
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.0172; a2 = -115.36; a3 = -1.1857; a4 = 13854.2;
			a5 = -0.262; c1 = -2.53; c2 = 38628.4; omg = 191730.;
			bg = -335.7; bs = 1.26; rc = 1.33; ra = 1.40;
			break;
		default:  // wrong mode
			return -1;
	}

	// calculation part, extended 06.06.2009 (TW)
	bh = bg + (298.15)*bs;
	rec = rc + fabs(zc)*(0.94+Gf);
	rea = ra + fabs(za)*Gf;
	X1 = - eta*nc*( fabs(pow(zc,3.))/pow(rec,2.) - zc/pow((3.082+Gf),2.) )
			- eta*na*( fabs(pow(za,3.))/pow(rea,2.) - za/pow((3.082+Gf),2.) );
	X2 = 2.*eta*nc*( fabs(pow(zc,4.))/pow(rec,3.) - zc/pow((3.082+Gf),3.) )
			+ 2.*eta*nc * ( fabs(pow(zc,4.))/pow(rec,3.) - zc/pow((3.082+Gf),3.) );
	omgpt = eta*( nc*pow(zc,2.)/rec + na*pow(za,2.)/rea );
	// omgpt = (1.66027e5)*(1./(0.94+rc+Gf)+1./(ra+Gf));
	domdt = X1*dGfdT;
	d2omdt2 = X2*pow(dGfdT,2.) + X1*d2GfdT2;
	domdp = X1*dGfdP;
	nbg = - ni*bg/2. + ni*bs*(Tk-298.15)/2. - c1*(Tk*log(Tk/298.15)-Tk+298.15)
				+ a1*(Pbar-1.) + a2*log((2600.+Pbar)/(2600.+1.))
				- c2*((1./(Tk-228.)-1./(298.15-228.))*(228.-Tk)/228.-Tk/(228.*228.)
				* log((298.15*(Tk-228.))/(Tk*(298.15-228.))))
				+ 1./(Tk-228.)*(a3*(Pbar-1.) + a4*log((2600.+Pbar)/(2600.+1.)))
				+ a5*(omgpt*(1./eps-1.)-omg*(1./(78.24513795)-1.)+(-5.798650444e-5)*omg*(Tk-298.15));
	nbh = - ni/2.*bh + c1*(Tk-298.15) - c2*(1./(Tk-228.)-1./(298.15-228.))
				+ a1*(Pbar-1.) + a2*log((2600.+Pbar)/(2600.+1.))
				+ ((2.*Tk-228.)/pow((Tk-228.),2.))*(a3*(Pbar-1.)+a4*log((2600.+Pbar)/(2600.+1.)))
				+ a5*( omgpt*(1./eps-1.) + omgpt*Tk*yborn - Tk*(1./eps-1.)*domdt
				- omg*(1./(78.24513795)-1.) - omg*(298.15)*(-5.798650444e-5) );
	nbj = c1 + c2/pow((Tk-228.),2.)  - ( (2.*Tk)/pow((Tk-228.),3.) )
				* ( a3*(Pbar-1.) + a4*log((2600.+Pbar)/(2600.+1.)) )
				+ a5*( omgpt*Tk*xborn + 2.*Tk*yborn*domdt - Tk*(1./eps-1.)*d2omdt2 );
	nbv = a1 + a2/(2600.+Pbar) + a3/(Tk-228.) + a4/((2600.+Pbar)*(Tk-228.))
				+ a5*(-omgpt*qborn + (1./eps-1.)*domdp);
	// b_gamma = nbg/(2.*log(10.)*1.98721*Tk);
	// bgam = b_gamma;
	bgam = nbg/(2.*log(10.)*(1.98721)*Tk)*2./ni;
	dbgdT = - nbh/(2.*log(10.)*(1.98721)*pow(Tk,2.))*2./ni;
	d2bgdT2 = - nbj/( 2.*log(10.)*(1.98721)*pow(Tk,2.))*2./ni - 2./Tk*dbgdT;
	dbgdP = nbv/(2.*log(10.)*(1.98721)*Tk)*2./ni;

	return 0;
}


long int TShvarov::IonsizeTP()
{
	double nc, na, ni, zc, za, c;

	switch ( flagElect )
	{
		case 1:  // NaCl
			nc = 1.; na = 1.; ni = 2.;
			zc = 1.; za = -1.;
			break;
		case 2:  // KCl
			nc = 1.; na = 1.; ni = 2.;
			zc = 1.; za = -1.;
			break;
		case 3:  // NaOH
			nc = 1.; na = 1.; ni = 2.;
			zc = 1.; za = -1.;
			break;
		case 4:  // KOH
			nc = 1.; na = 1.; ni = 2.;
			zc = 1.; za = -1.;
			break;
		default:  // wrong mode
			return -1;
	}

	c = 2./ni * ( nc*fabs(zc) + na*fabs(za) );
	ao = ac + c*Gf;
	daodT = c*dGfdT;
	d2aodT2 = c*d2GfdT2;
	daodP = c*dGfdP;

	return 0;
}


long int TShvarov::Gfunction()
{
	double T, P, D, beta, alpha, daldT;
	double g, dgdP, dgdT, d2gdT2;
	double TMAX = 1000., PMAX = 5000., TOL = 1.0e-4;

	// convert parameters
	T = Tk - 273.15;
	P = Pbar;
	D = RhoW[0];
	alpha = - RhoW[1]/RhoW[0];
	daldT = pow(alpha, 2.) - RhoW[2]/RhoW[0];
	beta = RhoW[3]/RhoW[0];

	// initialize g and derivatives to zero
	g = 0.0;
	dgdP = 0.0;
	dgdT = 0.0;
	d2gdT2 = 0.0;

	if ((T > TMAX+TOL) || (P > PMAX+TOL))
		return -1;
	else
		GShok2( T, P, D, beta, alpha, daldT, g, dgdP, dgdT, d2gdT2 );

	// assignments
	Gf = g;
	dGfdT = dgdT;
	d2GfdT2 = d2gdT2;
	dGfdP = dgdP;

	return 0;
}


long int TShvarov::GShok2( double T, double P, double D, double beta,
		double alpha, double daldT, double &g, double &dgdP, double &dgdT, double &d2gdT2 )
{
	double a, b, dgdD, /*dgdD2,*/ dadT, dadTT, dbdT, dbdTT, dDdT, dDdP,
		dDdTT, Db, dDbdT, dDbdTT, ft, dftdT, dftdTT, fp, dfpdP,
		f, dfdP, dfdT, d2fdT2, tempy;
	double C[6]  = {-0.2037662e+01,  0.5747000e-02, -0.6557892e-05,
			0.6107361e+01, -0.1074377e-01,  0.1268348e-04 };
	double cC[3] = { 0.3666666e+02, -0.1504956e-09,  0.5017997e-13 };

	if ( D >= 1.3 )
		return -1;
	// Sveta 19/02/2000 1-D < 0 pri D>1 => pow(-number, b) = -NaN0000 error
	double pw = fabs(1.0 - D); // insert Sveta 19/02/2000

	// calculation part
	a = C[0] + C[1]*T + C[2]*pow(T,2.);
	b = C[3] + C[4]*T + C[5]*pow(T,2.);
	g = a * pow(pw, b);

	dgdD = - a*b* pow(pw, (b - 1.0));
		// dgdD2 = a * b * (b - 1.0) * pow((1.0 - D),(b - 2.0));

	dadT = C[1] + 2.0*C[2]*T;
	dadTT = 2.0*C[2];
	dbdT = C[4] + 2.0*C[5]*T;
	dbdTT = 2.0*C[5];
	dDdT = - D * alpha;
	dDdP = D * beta;
	dDdTT = - D * (daldT - pow(alpha,2.));

		// Db = pow((1.0 - D), b);  Fixed by DAK 01.11.00
	Db = pow(pw , b);
	dDbdT = - b*pow(pw,(b - 1.0))*dDdT + log(pw)*Db*dbdT;

	dDbdTT = - (b*pow(pw,(b-1.0)) *dDdTT + pow(pw,(b - 1.0)) * dDdT * dbdT
				+ b * dDdT * (-(b - 1.0) * pow(pw,(b - 2.0)) * dDdT
				+ log(pw) * pow(pw,(b - 1.0)) * dbdT))
				+ log(pw) * pow(pw,b) * dbdTT
				- pow(pw,b) * dbdT * dDdT / (1.0 - D)
				+ log(pw) * dbdT * dDbdT;

	dgdP = dgdD * dDdP;
	dgdT = a * dDbdT + Db * dadT;
	d2gdT2 = a * dDbdTT + 2.0 * dDbdT * dadT + Db * dadTT;

	if((T < 155.0) || (P > 1000.0) || (T > 355.0))
		return 0;

	tempy = ((T - 155.0) / 300.0);
	ft = pow(tempy,4.8) + cC[0] * pow(tempy,16.);

	dftdT = 4.8 / 300.0 * pow(tempy, 3.8) + 16.0 / 300.0 * cC[0] * pow(tempy, 15.);

	dftdTT = 3.8 * 4.8 / (300.0 * 300.0) * pow(tempy, 2.8)
		+ 15.0 * 16.0 / (300.0*300.0) * cC[0] * pow(tempy,14.);

	fp = cC[1] * pow((1000.0 - P),3.) + cC[2] * pow((1000.0 - P),4.);

	dfpdP = -3.0 * cC[1] * pow((1000.0 - P),2.) - 4.0 * cC[2] * pow((1000.0 - P),3.);

	f = ft * fp;
	dfdP = ft * dfpdP;
	dfdT = fp * dftdT;
	d2fdT2 = fp * dftdTT;

	g -= f;
	dgdP -= dfdP;
	dgdT -= dfdT;
	d2gdT2 -= d2fdT2;

	return 0;
}

//=============================================================================================
// ELVIS activity model for aqueous electrolyte solutions
// (c) FFH Aug 2010
//=============================================================================================

// Generic constructor for the TELVIS class
TELVIS::TELVIS( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW ):
                TSolMod( sd )
{
    alloc_internal();
    z = arZ;
    m = arM;
    RhoW = dW;
    EpsW = eW;

#ifdef GEMSFIT_DEBUG
cout<<" in TELVIS::TELVIS constructor "<<endl;
#endif
}


TELVIS::~TELVIS()
{
        free_internal();
}

#ifndef ELVIS_SPEED
// declaration and initialization of dynamic storage
void TELVIS::alloc_internal()
{
        long i,j;
        IS = 0.0;

        beta0	= new double *[NComp];
        beta1	= new double *[NComp];
        alpha	= new double *[NComp];
        coord	= new double *[NComp];

        RA	= new double *[NComp];
        RC	= new double *[NComp];
        QC	= new double *[NComp];
        QA	= new double *[NComp];

        for( j=0; j<NComp; j++ )
        {
            for( i=0; i<NComp; i++ )
            {
                beta0[j][i] = 0.0;
                beta1[j][i] = 0.0;
                alpha[j][i] = 0.0;
                coord[j][i] = 0.0;
                RA[j][i]    = 0.0;
                RC[j][i]    = 0.0;
                QA[j][i]    = 0.0;
                QC[j][i]    = 0.0;
            }
        }

        R = new double [NComp];
        Q = new double [NComp];
        Phi = new double [NComp];
        Theta = new double [NComp];

        dRdP = new double [NComp];
        dRdT = new double [NComp];
        d2RdT2 = new double [NComp];
        dQdP = new double [NComp];
        dQdT = new double [NComp];
        d2QdT2 = new double [NComp];


        for( j=0; j<NComp; j++ )
        {
            R[j]      = 0.0;
            Q[j]      = 0.0;
            Phi[j]    = 0.0;
            Theta[j]  = 0.0;
            dRdP[j]   = 0.0;
            dRdT[j]   = 0.0;
            d2RdT2[j] = 0.0;
            dQdP[j]   = 0.0;
            dQdT[j]   = 0.0;
            d2QdT2[j] = 0.0;
        }

        EffRad  = new double [NComp];

        ELVIS_lnGam_DH      = new double [NComp];
        ELVIS_lnGam_Born    = new double [NComp];
        ELVIS_OsmCoeff_DH   = new double [NComp];
        ELVIS_lnGam_UNIQUAC = new double [NComp];

        for( j=0; j<NComp; j++ )
        {
            ELVIS_lnGam_DH[j]      = 0.0;
            ELVIS_lnGam_Born[j]    = 0.0;
            ELVIS_OsmCoeff_DH[j]   = 0.0;
            ELVIS_lnGam_UNIQUAC[j] = 0.0;
        }

        WEps   = new double *[NComp];
        U      = new double *[NComp];
        dU     = new double *[NComp];
        d2U    = new double *[NComp];
        Psi    = new double *[NComp];
        dPsi   = new double *[NComp];
        d2Psi  = new double *[NComp];
        TR     = new double *[NComp];

        dUdP   = new double *[NComp];
        dUdT   = new double *[NComp];
        d2UdT2 = new double *[NComp];


        for( j=0; j<NComp; j++ )
        {
            EffRad[j] = 0.0;
            WEps[j]   = new double [NComp];
            U[j]      = new double [NComp];
            dU[j]     = new double [NComp];
            d2U[j]    = new double [NComp];
            Psi[j]    = new double [NComp];
            dPsi[j]   = new double [NComp];
            d2Psi[j]  = new double [NComp];
            TR[j]     = new double [4];

            dUdP[j]   = new double [NComp];
            dUdT[j]   = new double [NComp];
            d2UdT2[j] = new double [NComp];
        }
        for( j=0; j<NComp; j++ )
        {
            for( i=0; i<NComp; i++ )
            {
                WEps[j][i] = 0.0;
                U[j][i] = 0.0;
                dU[j][i] = 0.0;
                d2U[j][i] = 0.0;
                Psi[j][i] = 1.0;
                dPsi[j][i] = 0.0;
                d2Psi[j][i] = 0.0;

                dUdP[j][i]   = 0.0;
                dUdT[j][i]   = 0.0;
                d2UdT2[j][i] = 0.0;
            }
        }

        for( j=0; j<NComp;j++ )
        {
            for( i=0; i<4; i++ )
            {
                TR[j][i]     = 0.0;
            }
        }
}


void TELVIS::free_internal()
{
    long j;
        for( j=0; j<NComp; j++ )
        {
            delete[]beta0[j];
            delete[]beta1[j];
            delete[]alpha[j];
            delete[]coord[j];

            delete[]RA[j];
            delete[]RC[j];
            delete[]QA[j];
            delete[]QC[j];

            delete[]WEps[j];
            delete[]U[j];
            delete[]dU[j];
            delete[]d2U[j];
            delete[]Psi[j];
            delete[]dPsi[j];
            delete[]d2Psi[j];
            delete[]TR[j];

            delete[]dUdP[j];
            delete[]dUdT[j];
            delete[]d2UdT2[j];
        }

        delete[]beta0;
        delete[]beta1;
        delete[]alpha;
        delete[]coord;

        delete[]RA;
        delete[]RC;
        delete[]QA;
        delete[]QC;

        delete[]R;
        delete[]Q;
        delete[]Phi;
        delete[]Theta;

        delete[]dRdP;
        delete[]dRdT;
        delete[]d2RdT2;
        delete[]dQdP;
        delete[]dQdT;
        delete[]d2QdT2;

        delete[]WEps;
        delete[]TR;
        delete[]EffRad;

        delete[]ELVIS_lnGam_DH;
        delete[]ELVIS_lnGam_Born;
        delete[]ELVIS_OsmCoeff_DH;
        delete[]ELVIS_lnGam_UNIQUAC;

        delete[]U;
        delete[]dU;
        delete[]d2U;
        delete[]Psi;
        delete[]dPsi;
        delete[]d2Psi;

        delete[]dUdP;
        delete[]dUdT;
        delete[]d2UdT2;

        R=NULL;
        Q=NULL;
        Phi=NULL;
        Theta=NULL;

        dRdP=NULL;
        dRdT=NULL;
        d2RdT2=NULL;
        dQdP=NULL;
        dQdT=NULL;
        d2QdT2=NULL;

        U=NULL;
        dU=NULL;
        d2U=NULL;

        WEps=NULL;
        TR=NULL;
        EffRad=NULL;

        ELVIS_lnGam_DH=NULL;
        ELVIS_lnGam_Born=NULL;
        ELVIS_OsmCoeff_DH=NULL;
        ELVIS_lnGam_UNIQUAC=NULL;

        U=NULL;
        dU=NULL;
        d2U=NULL;
        Psi=NULL;
        dPsi=NULL;
        d2Psi=NULL;

        dUdP=NULL;
        dUdT=NULL;
        d2UdT2=NULL;

}
#endif

#ifdef ELVIS_SPEED
void TELVIS::alloc_internal()
{
        long i,j;
        for( j=0; j<ELVIS_NCOMP; j++ )
        {
            R[j]     = 0.0;
            Q[j]     = 0.0;
            Phi[j]   = 0.0;
            Theta[j] = 0.0;
            dRdP[j]   = 0.0;
            dRdT[j]   = 0.0;
            d2RdT2[j] = 0.0;
            dQdP[j]   = 0.0;
            dQdT[j]   = 0.0;
            d2QdT2[j] = 0.0;
        }

        for( j=0; j<ELVIS_NCOMP; j++ )
        {
            EffRad[j] 	     	 	= 0.0;
            ELVIS_lnGam_DH[j]       = 0.0;
            ELVIS_lnGam_Born[j]     = 0.0;
            ELVIS_OsmCoeff_DH[j]    = 0.0;
            ELVIS_lnGam_UNIQUAC[j]  = 0.0;
        }

        for( j=0; j<ELVIS_NCOMP; j++ )
        {
            for( i=0; i<ELVIS_NCOMP; i++ )
            {
                WEps[j][i] 	= 0.0;
                U[j][i] 	= 0.0;
                dU[j][i] 	= 0.0;
                d2U[j][i] 	= 0.0;
                Psi[j][i] 	= 1.0;
                dPsi[j][i] 	= 0.0;
                d2Psi[j][i] = 0.0;

                dUdP[j][i]   = 0.0;
                dUdT[j][i]   = 0.0;
                d2UdT2[j][i] = 0.0;
           }
        }

        for( j=0; j<ELVIS_NCOMP;j++ )
        {
            for( i=0; i<4; i++ )
            {
                 TR[j][i]   = 0.0;
            }
        }
}

void TELVIS::free_internal(){}
#endif




// Initialization of vectors/arrays and calculation of T,P corrected binary interaction parameters
long int TELVIS::PTparam()
{
#ifdef GEMSFIT_DEBUG
cout<<" in TELVIS::PTparam()"<<endl;
#endif

    molfrac_update();
    IonicStrength();


    #ifdef GEMSFIT_DEBUG
    cout << "ELVIS 		PTparam():	Tk = " << Tk << endl;
    cout << "ELVIS 		PTparam():	Pbar = " << Pbar << endl;
    #endif

    long j, i, ip, i1, i2;
    int ii;
    double bet0, bet1, alp, cn, ra_, rc_, qa_, qc_;
    double u, psi, v, weps, diffU, Xw = 0.;
    double dudp, dudt, d2udt2;
    double* spec_frac = new double [NComp-1];
    double spec_sum = 0.0;
    u = 0.0; psi = 0.0; weps = 0.0;

    if ( NPcoef < 1 || NPar < 1 || NP_DC < 2 )
       return 1;

    // get index of water (assumes water is last species in phase)
    Xw = x[ NComp - 1 ];

    // read and transfer species-dependent parameters
    for (j=0; j<NComp; j++)
    {
        // Temperaure and pressure correction for effective radius

        //R[j] = aDCc[NP_DC*j+2] + aDCc[NP_DC*j+3]*(1-Xw)*(1-Xw);
        //R[j] = aDCc[NP_DC*j+1] + aDCc[NP_DC*j+2]*Tk*1e-02 + aDCc[NP_DC*j+3]*Pbar*1e-03 +  aDCc[NP_DC*j+4]*Tk*Tk*1e-04 + aDCc[NP_DC*j+5]*Tk*Pbar*1e-06 +	aDCc[NP_DC*j+6]*Tk*Tk*Tk*1e-08 + aDCc[NP_DC*j+7]*Tk*Tk*Pbar*1e-08;

        //cout<<aDCc[NP_DC*j+1]<<" "<<aDCc[NP_DC*j+2]*Tk<<" "<<aDCc[NP_DC*j+3]*Pbar<<" "<< aDCc[NP_DC*j+4]*Tk*Tk<<" "<< aDCc[NP_DC*j+5]*Tk*Pbar<<" "<<	aDCc[NP_DC*j+6]*Tk*Tk*Tk<<" "<< aDCc[NP_DC*j+7]*Tk*Tk*Pbar <<endl;

        //	Hardcoded Chlorine and Sodium Ion parameters
        //R[1] = 7.036 + (-0.03452)*Tk + 0.001275*Pbar +  0.00008134*Tk*Tk + (-0.000006346)*Tk*Pbar +	(-0.00000006675)*Tk*Tk*Tk + 0.000000007849*Tk*Tk*Pbar;

        // Concentration and temperature dependence on volume parameter
        R[j] = abs( aDCc[NP_DC*j+4] + aDCc[NP_DC*j+5]*Tk + aDCc[NP_DC*j+6]*Tk*Tk + \
               (aDCc[NP_DC*j+7] + aDCc[NP_DC*j+8]*Tk + aDCc[NP_DC*j+9]*Tk*Tk) * (1.-Xw) + \
               (aDCc[NP_DC*j+10] + aDCc[NP_DC*j+11]*Tk + aDCc[NP_DC*j+12]*Tk*Tk) * (1.-Xw)*(1.-Xw) );

        if( R[j]<0. )
        {
            cout << "R["<<j<<"] = " << R[j] << endl;
            cout << "	A = " << aDCc[NP_DC*j+4] + aDCc[NP_DC*j+5]*Tk + aDCc[NP_DC*j+6]*Tk*Tk << endl;
            cout << "	B = " << (aDCc[NP_DC*j+7] + aDCc[NP_DC*j+8]*Tk + aDCc[NP_DC*j+9]*Tk*Tk) * (1.-Xw) << endl;
            cout << "	C = " << (aDCc[NP_DC*j+10] + aDCc[NP_DC*j+11]*Tk + aDCc[NP_DC*j+12]*Tk*Tk) * (1.-Xw)*(1.-Xw) << endl;

        }



        // ONLY concentration dependence on volume parameter
//		R[j] = aDCc[NP_DC*j+2] +
//		       aDCc[NP_DC*j+3] * (1-Xw) +
//		       aDCc[NP_DC*j+4] * (1-Xw)*(1-Xw);

        Q[j] = abs( aDCc[NP_DC*j] + aDCc[NP_DC*j+2]*Tk + aDCc[NP_DC*j+3]*Tk*Tk );   // surface parameter q of UNIQUAC term

        if( Q[j]<0. )
        {
            cout << "Q["<<j<<"] = " << Q[j] << endl;
        }



// !!!!!!!!!!!!!!!!!!!!!!   ATTENTION  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! //
        R[1] = R[0];
        Q[1] = Q[0];
// !!!!!!!!!!!!!!!!!!!!!!   ATTENTION  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! //



        // !!!!!!!  For single T-P point fitting !!!!!!! :
        dRdP[j]   = aDCc[NP_DC*j+2];
        dRdT[j]   = aDCc[NP_DC*j+2];
        d2RdT2[j] = aDCc[NP_DC*j+2];

        dQdP[j]   = aDCc[NP_DC*j];
        dQdT[j]   = aDCc[NP_DC*j];
        d2QdT2[j] = aDCc[NP_DC*j];

    }

    for( j=0; j<(NComp-1); j++ )
    {
        EffRad[j] = aDCc[NP_DC*j+1];
    }


    for( ip=0; ip<NPar; ip++ )
    {
        weps = 0.0; u = 0.0; bet0 = 0.; bet1 = 0.; cn = 0.;
        i1 		= aIPx[MaxOrd*ip];
        i2 		= aIPx[MaxOrd*ip+1];
        if( aIPc[NPcoef*ip]>0.99 && aIPc[NPcoef*ip]<1.01 )
        {
            // Temperature correction for interaction parameter
            // u = aIPc[NPcoef*ip+1] + aIPc[NPcoef*ip+2]*pow((Tk-298.15), aIPc[NPcoef*ip+3]);

            // u = aIPc[NPcoef*ip+1] + aIPc[NPcoef*ip+2]*Tk + aIPc[NPcoef*ip+3]*Pbar + aIPc[NPcoef*ip+4]*Tk*Tk + aIPc[NPcoef*ip+5]*Tk*Pbar + aIPc[NPcoef*ip+6]*Tk*Tk*Tk +aIPc[NPcoef*ip+7]*Tk*Tk*Pbar;

            // Temperature dependence on interaction parameter
            u = aIPc[NPcoef*ip+1] + aIPc[NPcoef*ip+2]*Tk + aIPc[NPcoef*ip+3]*Tk*Tk;


//			u = aIPc[NPcoef*ip+1];

            u = aIPc[NPcoef*ip+1];
            U[i1][i2] 	    = u;
            U[i2][i1] 	    = u;

            weps = aIPc[NPcoef*ip+2];	// for Born term
            WEps[i1][i2] 	= weps;
                WEps[i2][i1]	= weps;

            cn = aIPc[NPcoef*ip+2] + aIPc[NPcoef*ip+3] * IS;		// variable coordination number
            coord[i1][i2] 	= weps;
                coord[i2][i1]	= weps;
/*
            bet0 = aIPc[NPcoef*ip+4];	// for Pitzer term
            beta0[i1][i2]   = bet0;
            beta0[i2][i1]   = bet0;

            bet1 = aIPc[NPcoef*ip+5];	// for Pitzer term
            beta1[i1][i2]   = bet1;
            beta1[i2][i1]   = bet1;

            alp = aIPc[NPcoef*ip+4];	// for Pitzer term
            alpha[i1][i2]   = alp;
            alpha[i2][i1]   = alp;
*/

            ra_ = aIPc[NPcoef*ip+4];	// for UNIQUAC term
            RA[i1][i2]   = ra_;
            RA[i2][i1]   = ra_;

            rc_ = aIPc[NPcoef*ip+5];	// for UNIQUAC term
            RC[i1][i2]   = rc_;
            RC[i2][i1]   = rc_;

            qa_ = aIPc[NPcoef*ip+6];	// for UNIQUAC term
            QA[i1][i2]   = qa_;
            QA[i2][i1]   = qa_;

            qc_ = aIPc[NPcoef*ip+7];	// for UNIQUAC term
            QC[i1][i2]   = qc_;
            QC[i2][i1]   = qc_;


/*				// Pressure and temperature derivatives of interaction parameter
            dudp   = aIPc[NPcoef*ip+3] + aIPc[NPcoef*ip+5]*Tk + aIPc[NPcoef*ip+7]*Tk*Tk;
            dUdP[i1][i2] 	    = dudp;
            dUdP[i2][i1] 	    = dudp;

            dudt   = aIPc[NPcoef*ip+2] + aIPc[NPcoef*ip+4]*Tk + aIPc[NPcoef*ip+5]*Pbar + aIPc[NPcoef*ip+6]*Tk*Tk + aIPc[NPcoef*ip+7]*Tk*Pbar;
            dUdT[i1][i2] 	    = dudt;
            dUdT[i2][i1] 	    = dudt;

            d2udt2 = aIPc[NPcoef*ip+4] + aIPc[NPcoef*ip+6]*Tk + aIPc[NPcoef*ip+7]*Pbar;
            d2UdT2[i1][i2] 	    = d2udt2;
            d2UdT2[i2][i1] 	    = d2udt2;
*/
            // !!!!!!!  For single T-P point fitting !!!!!!! :
            dudp   = aIPc[NPcoef*ip+1];
            dUdP[i1][i2] 	    = dudp;
            dUdP[i2][i1] 	    = dudp;

            dudt   = aIPc[NPcoef*ip+1];
            dUdT[i1][i2] 	    = dudt;
            dUdT[i2][i1] 	    = dudt;

            d2udt2 = aIPc[NPcoef*ip+1];
            d2UdT2[i1][i2] 	    = d2udt2;
            d2UdT2[i2][i1] 	    = d2udt2;
        }
        else if( aIPc[NPcoef*ip]>1.99 && aIPc[NPcoef*ip]<2.01 )
        {	// Temperature correction for permittivity term

//                                weps = aIPc[NPcoef*ip+1] + aIPc[NPcoef*ip+2]*(Tk-298.15)*1e-02 + aIPc[NPcoef*ip+3]*(Pbar-1.)*1e-03 + aIPc[NPcoef*ip+4]*(Tk-298.15)*(Tk-298.15)*1e-04 + aIPc[NPcoef*ip+5]*(Tk-298.15)*(Pbar-1.)*1e-06 + aIPc[NPcoef*ip+6]*(Tk-298.15)*(Tk-298.15)*(Tk-298.15)*1e-08 +aIPc[NPcoef*ip+7]*(Tk-298.15)*(Tk-298.15)*(Pbar-1.)*1e-08;
            weps = aIPc[NPcoef*ip+1];

            WEps[i1][i2] 	= weps;
            WEps[i2][i1]	= weps;

        }
    }

#ifdef ELVIS_DEBUG
    for( ip=0; ip<NPar; ip++ )
    {
        i1 		= aIPx[MaxOrd*ip];
        i2 		= aIPx[MaxOrd*ip+1];
        cout<<"code = "<<aIPc[NPcoef*ip]<<" | WEps["<<i1<<"]["<<i2<<"] = "<<WEps[i1][i2]<<" | u["<<i1<<"]["<<i2<<"] = "<<U[i1][i2]<<endl;
        cout<<"i1 = "<<i1<<" | i2 = "<<i2<<endl;
    }
#endif

    // calculate Psi and its partial derivatives
    for( j=0; j<NComp; j++ )
    {
        for (i=0; i<NComp; i++)
        {
            diffU = U[j][i]-U[i][i];
            psi = exp( -diffU/Tk );
            v = (U[j][i]-U[i][i])/pow(Tk,2.) - (dU[j][i]-dU[i][i])/Tk;
            Psi[j][i] = psi;
        }
    }

// ------------------------------------------------------------------------------------- //
/*
    double cat_sum = 0., ani_sum = 0.;

    for( i=0;i<(NComp-1);i++ )
    {
        if( m[i] > 0. )
        {
            cat_sum += m[i];
        }
        if( m[i] < 0. )
        {
            ani_sum += m[i];
        }
    }

    for( i=0;i<(NComp-1);i++ )
    {
        R[i] = 0.;

        Q[i] = 0.;
        for( int k=0;k<(NComp-1);k++ )
        {
            if( m[k] > 0. )
            {
                R[i] = R[i] + m[k]*RA[i][k]/cat_sum;
                Q[i] = Q[i] + m[k]*QA[i][k]/cat_sum;
            }
            if( m[k] < 0. )
            {
                R[i] = R[i] + m[k]*RC[i][k]/ani_sum;
                Q[i] = Q[i] + m[k]*QC[i][k]/ani_sum;
            }
        }


    }

*/
// ------------------------------------------------------------------------------------- //

    // Debye-Huckel Term functions

    // species fractions:
    for( i=0;i<(NComp-1);i++ )
    {
        spec_sum += m[i];
    }
    for( i=0;i<(NComp-1);i++ )
    {
        spec_frac[i]=m[i]/spec_sum;
    }

    aDH = 0.0;
    for( i=0;i<(NComp-1);i++ )
    {
        // multiply by 2 to get electrolyte size (valid for 1:1, 2:2 electrolytes)
        aDH += 2 * EffRad[i]*spec_frac[i];
#ifdef ELVIS_DEBUG
        cout<<"EffRad = "<<EffRad[i]<<endl;
        cout<<"aDH = "<<aDH<<endl;
#endif
    }

    A = 1.824829238E6 * pow(RhoW[0],0.5) / pow(EpsW[0]*Tk, 3./2.);
    //dAdP	= 0.5*1.824829238E6 * Tk * (EpsW[0]*RhoW[3]-3*RhoW[0]*EpsW[3]) / ( Tk*Tk*EpsW[0]*EpsW[0] * sqrt(RhoW[0]) );

    dAdP	= A * ( 0.5*RhoW[3]/RhoW[0] - 1.5*EpsW[3]/EpsW[0] );

    dAdT	= A * ( RhoW[1] / ( 2*RhoW[0]*pow(Tk,1.5) ) - 1.5*( EpsW[0] + Tk * EpsW[1]) / (EpsW[0]*pow(Tk,2.5)) );
    d2AdT2	= - 1.5 * A * ( EpsW[0] + Tk * EpsW[1] ) * RhoW[1] /( EpsW[0] * Tk * RhoW[0] ) + 3.75 * A * ( EpsW[0] + Tk * EpsW[1] ) * ( EpsW[0] + Tk * EpsW[1] ) / (Tk * EpsW[0] * Tk * EpsW[0]) - 1.5 * ( 2*EpsW[1] + Tk*EpsW[2] ) / pow(EpsW[0]*Tk, 2.5) + A * ( RhoW[2] / (2*sqrt(RhoW[0])) - RhoW[1] / (4*pow(RhoW[0],1.5)) ) / sqrt(RhoW[0]);

    B = 50.29158649E8 * pow(RhoW[0],0.5) / pow(EpsW[0]*Tk, 0.5);
    //dBdP	= 0.5*50.29158649E8 * Tk * (EpsW[0]*RhoW[3]-RhoW[0]*EpsW[3]) / ( pow(Tk*EpsW[0],1.5) * sqrt(RhoW[0]) );

    dBdP	= 0.5 * B * ( RhoW[3] / RhoW[0] - EpsW[3] / EpsW[0] );

    dBdT	= B * ( RhoW[1] / ( 2*RhoW[0] ) - 0.5*( EpsW[0] + Tk * EpsW[1]) / (EpsW[0]*Tk) );
    d2BdT2	= - B / (Tk*EpsW[0]) * ( RhoW[1]*( EpsW[0] + Tk * EpsW[1] ) / RhoW[0] + 0.75 * ( EpsW[0] + Tk * EpsW[1] ) * ( EpsW[0] + Tk * EpsW[1] ) / ( Tk * EpsW[0] ) - 0.5 * ( 2 * EpsW[1] + Tk * EpsW[2] ) + sqrt(Tk*EpsW[0]) * ( RhoW[2] / (2*RhoW[0]) - RhoW[1] * RhoW[1] / (4 * pow(RhoW[0],1.5)) ) );

    aDH = aDH*1e-8;

    delete[]spec_frac;
    return 0;
}


// Calculates activity coefficients
long int TELVIS::MixMod()
{
        long int j;
        double osmcoeff, msum;

//#ifdef GEMSFIT_DEBUG
//cout << " TELVIS::MixMod():	m[0] = "<<m[0]<<endl;
//#endif

/*	// compute osmotic coefficient of solvent
        osmcoeff = Int_OsmCoeff();
*/
        // compute activity coefficients of solute species
        CalcAct();

//#ifdef GEMSFIT_DEBUG
//cout<<" TELVIS::MixMod():	lnGamma[0] = "<<lnGamma[0]<<endl;
//cout<<" TELVIS::MixMod()	lnGamma[1] = "<<lnGamma[1]<<endl;
//#endif

  /*      msum = 0.0;
        for (j=0; j<(NComp-2); j++)
        {
                msum = msum + m[j];
        }

        lnGamma[NComp-1] = osmcoeff * msum / 55.508435061791985;
*/
/*        for (j=0; j<(NComp-1); j++)
        {
                 cout<<"lnGamma[j] = "<<lnGamma[j]<<endl;
        }
*/
        // Penalty function
        /*if( R[0]<1e-10 || R[0]>100 || Q[0]<1e-10 || Q[0]>100 )
        {
            lnGamma[0] = 1e5;
        }
        // ONLY FOR DEBUG SET gam(H+) and gam(OH-) to 1 !!!!
        lnGamma[2] = 0.0;
        lnGamma[3] = 0.0;
        lnGamma[4] = 0.0;
        */
return 0;
}



long int TELVIS::CalcAct()
{
        int j = 0;
        // calculation of ionic strength
        IonicStrength();

        // DH, Born and UNIQUAC contributions to activity coefficients
        ELVIS_DH(ELVIS_lnGam_DH, ELVIS_OsmCoeff_DH);
// 		ELVIS_Born(ELVIS_lnGam_Born);
        ELVIS_UNIQUAC(ELVIS_lnGam_UNIQUAC);

        // Solvent activity
        molT = 0.;
        for( j=0; j<(NComp-1); j++ )
        {
          molT += m[j];
        }
        double Gamma_gamma = log(1.+0.001801*molT);
        int w = NComp - 1;

        for( j=0;j<NComp;j++ )
        {
                // Assembling of activity coefficient contributions and converting to molality scale
                //lnGamma[j]   = ELVIS_lnGam_DH[j] + ELVIS_lnGam_Born[j] + ELVIS_lnGam_UNIQUAC[j] + log(x[w]);
                lnGamma[j]   = ELVIS_lnGam_DH[j] + ELVIS_lnGam_UNIQUAC[j] + ELVIS_lnGam_Born[j];

                // write debug results
                gammaDH[j]   = ELVIS_lnGam_DH[j];
                //gammaBorn[j] = ELVIS_lnGam_Born[j];
                gammaQUAC[j] = ELVIS_lnGam_UNIQUAC[j];

#ifdef ELVIS_DEBUG
                cout<<"ELVIS	Tk = "<<Tk<<endl;
                cout<<"ELVIS	Pbar = "<<Pbar<<endl;
                cout<<"ELVIS	m["<<j<<"] = "<<m[j]<<endl;
                cout<<"ELVIS	lnGamma["<<j<<"]	= "<<lnGamma[j]<<endl;
                cout<<"ELVIS	gammaDH["<<j<<"]	= "<<gammaDH[j]<<endl;
                //cout<<"ELVIS	gammaBorn["<<j<<"]	= "<<gammaBorn[j]<<endl;
                cout<<"ELVIS	gammaQUAC["<<j<<"]	= "<<gammaQUAC[j]<<endl;
#endif

        }


#ifdef ELVIS_DEBUG
        for( j=0;j<(NComp-1);j++ )
        {
                cout<<"lnGamma["<<j<<"] = "<<lnGamma[j]<<endl;
        }
        int fieldwidth = 20;
        ofstream my_gemactcoef;
        my_gemactcoef.open("my_gemactcoef.txt",ios::app);
        my_gemactcoef.width(fieldwidth);
        my_gemactcoef.precision(12);
        my_gemactcoef 	<< right << setw(fieldwidth) << "lnGamma"   << right << setw(fieldwidth) << "gammaDH" \
                       /* << right << setw(fieldwidth) << "gammaBorn"*/ << right << setw(fieldwidth) << "gammaQUAC" \
                        << right << setw(fieldwidth) << "gammaC"    << right << setw(fieldwidth) << "gammaR" \
                        << right << setw(fieldwidth) << "m[0]" << right << setw(fieldwidth) << "m[1]" \
                        << endl;
        for(int i=0; i<(NComp-1); i++ )
        {
           my_gemactcoef << right << setw(fieldwidth) << lnGamma[i] 	 << right << setw(fieldwidth) << gammaDH[i] \
                         /*<< right << setw(fieldwidth) << gammaBorn[i]*/ << right << setw(fieldwidth) << gammaQUAC[i] \
                         << right << setw(fieldwidth) << gammaC[i] 	 << right << setw(fieldwidth) << gammaR[i] \
                         << right << setw(fieldwidth) << m[0]         << right << setw(fieldwidth) << m[1] \
                         << endl;
        }
        my_gemactcoef<<endl;
        my_gemactcoef.close();
#endif

return 0;
}


void TELVIS::get_lnGamma( vector<double> &ln_gamma )
{
        ln_gamma.resize(5);
        copy( lnGamma, lnGamma+4, ln_gamma.begin() );
}


void TELVIS::ELVIS_DH(double* ELVIS_lnGam_DH, double* ELVIS_OsmCoeff_DH)
{
    // Debye-Huckel Term
    double a0, A_gamma, B_gamma, lambda, b, xbx, rhow, epsw;
    int i,j;
    double Mw = 0.01801528;

    rhow = RhoW[0];
    epsw  = EpsW[0];


    // species fractions:
    double* spec_frac = new double [NComp-1];
    double spec_sum = 0.0;
    for( i=0;i<(NComp-1);i++ )
    {
            spec_sum += m[i];
    }
    for( i=0;i<(NComp-1);i++ )
    {
            spec_frac[i]=m[i]/spec_sum;
    #ifdef ELVIS_DEBUG
        cout<<"spec_frac = "<<spec_frac[i]<<endl;
    #endif
    }
    a0 = 0.0;
    for( i=0;i<(NComp-1);i++ )
    {
            // multiply by 2 to get electrolyte size (valid for 1:1, 2:2 electrolytes)
            a0 += 2 * EffRad[i]*spec_frac[i];
            //a0 = 1*EffRad[0] + 2*EffRad[1];
    #ifdef ELVIS_DEBUG
            cout<<"EffRad = "<<EffRad[i]<<endl;
            cout<<"a0 = "<<a0<<endl;
    #endif
    }
    // A_gamma referring to log10 Debye Huckel term
    A_gamma = 1.824829238E6 * pow(rhow,0.5) / pow(epsw*Tk, 3./2.);
    B_gamma = 50.29158649E8 * pow(rhow,0.5) / pow(epsw*Tk, 0.5);

    // Conversion from log10 to ln Debye Huckel term
    A_gamma = A_gamma * log( 10 );


    a0 = a0*1e-8;
    //lambda = (1.+a0*B_gamma*sqrt(IS));
    b = a0*B_gamma;

    // ONLY FOR DEBUGGING
    //	b = 1.5;
    //

    lambda = (1.+b*sqrt(IS));


#ifdef ELVIS_DEBUG
    cout<<"Helgeson "<<" | a0 = "<<a0<<" | A_gamma = "<<A_gamma<<" | B_gamma = "<<B_gamma<<" | b = "<<b<<endl;
#endif

    // solutes
    for( j=0; j<(NComp-1); j++ )
    {
        ELVIS_lnGam_DH[j] = ( -A_gamma * z[j]*z[j]*pow(IS,0.5)/lambda ); // / log10(exp(1.0));//* log(10); // + Gamma_gamma;
#ifdef ELVIS_DEBUG
    cout<<" ELVIS: A_gamma = "<<A_gamma<<", B_gamma = "<<B_gamma<<", b = "<<b<<", loggam_DH["<<j<<"] = "<<ELVIS_lnGam_DH[j]<<endl;
    cout<<"lggamDM Helgeson 1982 = "<< (-A_gamma * z[j]*z[j]*pow(IS,0.5)/lambda) * log(10) <<endl; // + Gamma_gamma;
#endif
    }

    // solvent water
    xbx = 1+b*pow(IS,0.5);

    // solvent: activity: formula form dissertation of Kaj Thomsen (1997)
    ELVIS_lnGam_DH[NComp-1] = (Mw*2*A_gamma/(b*b*b))*(xbx - 1/xbx - 2*log(xbx));

    // solvent: osmotic coefficient, formula form helgeson 1981, p. 1355
    for( j=0; j<(NComp-1); j++ )
    {
          ELVIS_OsmCoeff_DH[j] = z[j]*z[j]*(A_gamma * pow(IS,0.5)*(xbx- (1./(xbx)) - 2.*log(xbx))) / ( pow(IS,(3./2.)) *b*b*b);
    }

    delete[] spec_frac;
}



void TELVIS::ELVIS_UNIQUAC( double* ELVIS_lnGam_UNIQUAC )
{
        int j, i, l, k, w;
        double Mw, Xw, b, RR, QQ, K, L, M;
        double gamC = 0.0; double gamR = 0.0;
        double lnGam = 0.0; double Gam = 0.0;
        b = 1.5; Mw = 0.01801528;

        // get index of water (assumes water is last species in phase)
        w = NComp - 1;
        Xw = x[w];

        // calculation of Phi and Theta terms
        for( j=0; j<NComp; j++ )
        {
            RR = 0.0;
            QQ = 0.0;
            for( i=0; i<NComp; i++ )
            {
                RR += x[i]*R[i];
                QQ += x[i]*Q[i];
            }
            Phi[j] = x[j]*R[j]/RR;

#ifdef ELVIS_DEBUG
    cout<<"Phi["<<j<<"] = "<<Phi[j]<<endl;
    cout<<"x["<<j<<"] = "<<x[j]<<endl;
    cout<<"R["<<j<<"] = "<<R[j]<<endl;
#endif

            Theta[j] = x[j]*Q[j]/QQ;
        }


// -------------------------------- COORDINATION NUMBER --------------------------------------------- //
        // species fractions:
        double spec_sum=0.;
        CN = 0.;
        vector<double> spec_frac((NComp-1),0);

        CN = 5.;
/*
        for( i=0;i<(NComp-1);i++ )
        {
            if( z[i] != 0 )
            {
                spec_sum += m[i];
            }
        }
        for( i=0;i<(NComp-1);i++ )
        {
            spec_frac[i]=m[i]/spec_sum;
        }
        for( i=0;i<(NComp-1);i++ )
        {
            for( j=i+1;j<(NComp-1);j++ )
            {
                CN += (spec_frac[i]*spec_frac[j]) * coord[i][j];
            }
        }
*/
// -------------------------------- COORDINATION NUMBER --------------------------------------------- //


        // loop over species
        for( j=0; j<NComp; j++ )
        {
                // species other than water solvent
                if (j < w)
                {
                        K = 0.0;
                        L = 0.0;
                        for( k=0; k<NComp; k++ )
                        {
                                M = 0.0;
                                for( l=0; l<NComp; l++ )
                                {
                                        M += Theta[l]*Psi[l][k];
                                }
                                if( fabs(M) < 1e-20 )
                                {
                                        M = 1e-20;
                                }
                                K += Theta[k]*Psi[k][j];
                                L += Theta[k]*Psi[j][k]/M;
                        }
//			gamDH = - pow(z[j],2.)*A*sqrt(IS)/(1.+b*sqrt(IS));

                        const int DivideByZero_or_NegativeLogarithm = 10;
                        try
                        {
/*								if( R[j]<1e-10 || Q[j]<1e-10 )
                                {
  */                               //if( Phi[j]<0.0 || x[j]<=0.0 || R[j]<0.0 || (Phi[j]/Theta[j])<0.0 || (R[j]*Q[w]/(R[w]*Q[j]))<0.0 || K<0.0 || Psi[w][j]<0.0 ){
                                 //       cout<<"Phi["<<j<<"] = "<<Phi[j]<<", x["<<j<<"] = "<<x[j]<<", Theta["<<j<<"] = "<<Theta[j]<<endl;
                                 //       cout<<"R["<<j<<"] = "<<R[j]<<", Q["<<j<<"] = "<<Q[j]<<", Q["<<w<<"] = "<<Q[w]<<", R["<<w<<"] = "<<R[w]<<endl;
                                 //       cout<<"K = "<<K<<", Psi["<<w<<"]["<<j<<"] = "<<Psi[w][j]<<endl;
/*									if( Phi[j]<1e-20 ){ Phi[j] = 1e-20; }
                                    if( Theta[j]<1e-20 ){ Theta[j] = 1e-20; }
                                    throw DivideByZero_or_NegativeLogarithm;
                                }
                                if( Psi[w][j]<1e-20 ){ Psi[w][j] = 1e-20; }
*/

                                gamC = log(Phi[j]/x[j]) - Phi[j]/x[j] - log(R[j]/R[w]) + R[j]/R[w]
                                       //         - 5.0*Q[j] * ( log(Phi[j]/Theta[j]) - Phi[j]/Theta[j]
                                                - CN * Q[j] * ( log(Phi[j]/Theta[j]) - Phi[j]/Theta[j]
                                                - log(R[j]*Q[w]/(R[w]*Q[j])) + R[j]*Q[w]/(R[w]*Q[j]) );
                                gamR = Q[j] * ( - log(K) - L + log(Psi[w][j]) + Psi[j][w] );

                                if( R[j]<1e-10 ){
                                //        gamC = 1e40;
                                }

                        }
                        catch( int err )
                        {
                                if( err==DivideByZero_or_NegativeLogarithm )
                                {
                                    cerr<<" R["<<j<<"] = "<<R[j]<<" | Q["<<j<<"] = "<<Q[j]<<endl;
                                //    cerr<<": Careful: a zero-divide or negative-logarithm occured in the UNIQUAC part of ELVIS !!!! Check your interaction and component specific parameters !!!! "<<" R["<<j<<"] = "<<R[j]<<endl;
                                }
                        }


                        ELVIS_lnGam_UNIQUAC[j] = gamC + gamR;

                        gammaC[j] = gamC;
                        gammaR[j] = gamR;

#ifdef ELVIS_DEBUG
                        if( j==0 || j==1 )
                        {
                            cout<<"m[0] = "<<m[0]<<", m[1] = "<<m[1]<<endl;
                            cout<<"x[0] = "<<x[0]<<", x[1] = "<<x[1]<<endl;
                            cout<<"IS = "<<IS<<endl;
                            cout<<"z["<<j<<"] = "<<z[j]<<endl;
                            cout<<"RhoW[0] = "<<RhoW[0]<<endl;
                            cout<<"EpsW[0] = "<<EpsW[0]<<endl;
                            cout<<"A = "<<A<<endl;
                            cout<<"R["<<j<<"] = "<<R[j]<<endl;
                            cout<<"Q["<<j<<"] = "<<Q[j]<<endl;
                            cout<<"x["<<j<<"] = "<<x[j]<<endl;
                            cout<<"Phi["<<j<<"] = "<<Phi[j]<<endl;
                            cout<<"Theta["<<j<<"] = "<<Theta[j]<<endl;
                            cout<<"gammaC["<<j<<"] = "<<gammaC[j]<<endl;
                            cout<<"gammaR["<<j<<"] = "<<gammaR[j]<<endl;
                        }
#endif

                }

                // water solvent
                else
                {
                        K = 0.0;
                        L = 0.0;
                        for (k=0; k<NComp; k++)
                        {
                                M = 0.0;
                                for (l=0; l<NComp; l++)
                                {
                                        M += Theta[l]*Psi[l][k];
                                }

                                K += Theta[k]*Psi[k][j];
                                L += Theta[k]*Psi[j][k]/M;
                        }

//			gamDH = Mw*2.*A/pow(b,3.) * ( 1. + b*sqrt(IS) - 1./(1.+b*sqrt(IS)) - 2*log(1.+b*sqrt(IS)) );
//                        gamC = log(Phi[j]/x[j]) + 1. - Phi[j]/x[j] - 5.0*Q[j] * ( log(Phi[j]/Theta[j]) + 1. - Phi[j]/Theta[j] );
                        gamC = log(Phi[j]/x[j]) + 1. - Phi[j]/x[j] - CN * Q[j] * ( log(Phi[j]/Theta[j]) + 1. - Phi[j]/Theta[j] );
                        gamR = Q[j] * (1. - log(K) - L );
                        lnGam = gamC + gamR;
#ifdef ELVIS_DEBUG
                            cout<<"z["<<j<<"] = "<<z[j]<<endl;
                            cout<<"R["<<j<<"] = "<<R[j]<<endl;
                            cout<<"Q["<<j<<"] = "<<Q[j]<<endl;
                            cout<<"x["<<j<<"] = "<<x[j]<<endl;
                            cout<<"Phi["<<j<<"] = "<<Phi[j]<<endl;
                            cout<<"Theta["<<j<<"] = "<<Theta[j]<<endl;
                            cout<<"gammaC["<<j<<"] = "<<gammaC[j]<<endl;
                            cout<<"gammaR["<<j<<"] = "<<gammaR[j]<<endl;
#endif


                        // Add combinatorial and residual terms without infinite dilution terms (!!!!)
                        ELVIS_lnGam_UNIQUAC[j] = gamC + gamR;

                        // write debug results
                        Gam 	  = exp(lnGam);
                }
        }

}


long int TELVIS::ExcessProp( double *Zex )
{
    long int j, i, k, w;
    double Mw, Xw;
    double gDH, gC, gR, hR, cpR, gCI, gRI, gCX, gRX, dg=0., d2g=0.;
    double DHTv, CTv, RTv, rtv1;
    double DHTg, CTg, RTg, rtg;
    double SRI = 0.0, xr = 0.0, xq = 0.0, xDrDp = 0.0, xDqDp = 0.0;
    Gex = 0.0; Hex = 0.0; Sex = 0.0; CPex = 0.0; Vex = 0.0;
    gC  = 0.0; gR  = 0.0; hR  = 0.0; cpR  = 0.0;
    DHTv = 0.0; CTv  = 0.0; RTv  = 0.0; rtv1 = 0.0;
    DHTg = 0.0; CTg  = 0.0; RTg  = 0.0; rtg  = 0.0;

    vector<double> thetapsi;
    thetapsi.resize(NComp);

    Mw = 0.01801528;

    // get index of water (assumes water is last species in phase)
    w = NComp -1;
    Xw = x[w];

    // calculation of ionic strength
    IonicStrength();

    SRI = sqrt(IS);

    // calculation of Phi and Theta terms
    for( j=0; j<NComp; j++ )
    {
        xr = 0.0;
        xq = 0.0;
        xDrDp = 0.0;
        xDqDp = 0.0;
        for (i=0; i<NComp; i++)
        {
            xr += x[i]*R[i];
            xq += x[i]*Q[i];
            xDrDp += x[i]*dRdP[i];
            xDqDp += x[i]*dQdP[i];
        }
        Phi[j] = x[j]*R[j]/xr;
        Theta[j] = x[j]*Q[j]/xq;
    }

    for( i=0; i<NComp; i++ )
    {
        thetapsi[i] = 0.0;
        for( j=0; j<NComp; j++ )
        {
            thetapsi[i] += Theta[j] * Psi[j][i];
        }
    }


    // Bulk excess volume and enthalpy

/*
cout << "ExcessProp() Rhow[0] = " << RhoW[0] << endl;
cout << "ExcessProp() Rhow[1] = " << RhoW[1] << endl;
cout << "ExcessProp() Rhow[2] = " << RhoW[2] << endl;
cout << "ExcessProp() Rhow[3] = " << RhoW[3] << endl;
*/

    // Debye Huckel term
    DHTg = - Xw * Mw * 4 * A / ( aDH*aDH*aDH*B*B*B ) * ( log(1 + aDH*B*SRI ) - aDH*B*SRI + aDH*aDH*B*B*IS/2 ) ;
    DHTv = ( 4 * Xw * Mw /pow((aDH*B),3) ) * ( ( - dAdP + 3 * A * dBdP / B ) * ( -aDH*B*SRI + 0.5*aDH*aDH*B*B*IS + log( 1 + aDH*B*SRI) ) - A * ( -aDH*dBdP*SRI + aDH*aDH*B*dBdP*IS + aDH*SRI*dBdP/(1 + aDH*B*SRI) ) );

    // outer species loop
    for( i=0; i<NComp; i++ )
    {
        //	inner species loop
        for( k=0; k<NComp; k++ )
        {
            rtv1 += Theta[k] * Psi[k][i] * dQdP[k] / Q[k] - Theta[k] * Psi[k][i] * xDqDp / xq + Theta[k] * Psi[k][i] * ( dUdP[i][i] - dUdP[k][i] ) / Tk;
            rtg  += Theta[k] * Psi[k][i];
        }

        // Combinatorial term
        CTv += x[i] * ( dRdP[i] - R[i] * xDrDp / xr ) - 5 * ( x[i] * dQdP[i] * log( Phi[i] / Theta[i] ) + Theta[i]*x[i]*Q[i]/Phi[i] * ( - Phi[i]*dQdP[i]/(Theta[i]*Q[i]) + Phi[i]*xDqDp/(x[i]*Q[i]) + x[i]*dRdP[i]/(Theta[i]*xr) - Phi[i]*xDrDp/(Theta[i]*xr) ) );
        CTg += x[i] * log( Phi[i]/x[i] ) - 5 * Q[i]*x[i]*log( Phi[i]/Theta[i] );

        // Residual term
        RTv += - ( x[i] * dQdP[i] * log( thetapsi[i] ) + x[i] * Q[i] * rtv1 / thetapsi[i] );
        RTg += Q[i]*x[i]*log( rtg );

    }
    Gex = ( DHTg + CTg - RTg ) * R_CONST * Tk;

    Vex = DHTv /*+ CTv + RTv*/;

cout << "Vex = " << Vex << endl;

/*
cout << "Gex = " << Gex << endl;
cout << "Vex = " << Vex << endl;
*/

    // increment thermodynamic properties
    //Gex = ( gDH + gRX + gCX - gRI - gCI ) * R_CONST * Tk;
    Hex = dg * pow(Tk,2.) * R_CONST;
    CPex = ( 2.*Tk*dg + pow(Tk,2.)*d2g ) * R_CONST;
    Sex = (Hex-Gex)/Tk;
    Aex = Gex - Vex*Pbar;
    Uex = Hex - Vex*Pbar;

    // assigments (excess properties)
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
long int TELVIS::IdealProp( double *Zid )
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


// Calculate ionic strength
long int TELVIS::IonicStrength()
{
        long int j;
        double Mw, Xw;
        IS = 0.0; Mw = 0.01801528;
        Xw = x[NComp-1];

        for (j=0; j<(NComp-1); j++)
        {
                IS += 0.5*x[j]*pow(z[j],2.)/(Xw*Mw);
        }

        return 0;
}


void TELVIS::ELVIS_Born(double* ELVIS_lnGam_Born)
{
    int i=0,j=0;

/*	double x=0.;
    double part1=0., part2=0.;

    // from Pitzer F-term
    for( i=0; i<(NComp-1); i++ )
    {
        if( m[i] > 0. )
        {
            for( j=0; j<(NComp-1); j++ )
            {
                if( m[j] < 0. )
                {
                    x = alpha[i][j] * sqrt(IS);
                    part1 += m[i]*m[j]*( -beta1[i][j] * 2 * (1 - (1+x+x*x/2.) * exp(-x) ) / (x*x*IS) );
                }
            }
        }
    }

    //
    for( i=0; i<(NComp-1); i++ )
    {
        part2 = 0.;
        for( j=0; j<(NComp-1); j++ )
        {
            if( m[i]<0. && m[j]>0. || m[i]>0. && m[j]<0. )
            {
                x = alpha[i][j] * sqrt(IS);
                part2 += 2 * m[j] * ( beta0[i][j] * beta1[i][j] * 2*(1-(1+x)*exp(-x))/(x*x) );
            }
        }
        ELVIS_lnGam_Born[ i ] = z[i] * z[i] * part1 + part2;	}
*/

        // species fractions:
        double spec_sum=0.;
        double EP = 0.;
        vector<double> spec_frac((NComp-1),0);

        for( i=0;i<(NComp-1);i++ )
        {
            if( z[i] != 0 )
            {
                spec_sum += m[i];
            }
        }
        for( i=0;i<(NComp-1);i++ )
        {
            spec_frac[i]=m[i]/spec_sum;
        }
        for( i=0;i<(NComp-1);i++ )
        {
            for( j=i+1;j<(NComp-1);j++ )
            {
                EP += (spec_frac[i]*spec_frac[j]) * ( beta0[i][j] + beta1[i][j]*IS*IS );
            }
        }


        ELVIS_lnGam_Born[ i ] = EP;

}


/*void TELVIS::ELVIS_Born(double* ELVIS_lnGam_Born)
{
        // Solvation Term from Helgeson Model (HKF, 1981; Shock 1992)

        double eta = 1.66027; //e5;  		// [A cal mol^-1] from Helgeson 1981
        //double R   = 1.9872;     		// gas constant /  cal/mol/�K
        double* omega1 = new double [(NComp-1)];
        memset(omega1, 0.0, sizeof omega1);
        double* omega = new double [(NComp-1)];  		// Born parameter
        memset(omega, 0.0, sizeof omega);
        double* spec_frac = new double [NComp-1];
        int i,j;

        for( j=0; j<(NComp-1); j++ )
        {
            omega1[j] = (eta*z[j]*z[j]);          // Helgeson 1981. formula 130 (only enumerator)
        }
        for( j=0; j<(NComp-1); j++ )
        {
            omega[j] = omega1[j]/(EffRad[j]);     // Helgeson 1981, formula 130
        }

#ifdef ELVIS_DEBUG
        ofstream myomega;
        myomega.open ("myomega.txt",ios::app);
        for( j=0; j<(NComp-1); j++ )
        {
            myomega<<"omega["<<j<<"] = "<<omega[j]<<endl;
        }
        myomega<<endl;

        for( j=0; j<(NComp-1); j++ )
        {
            myomega<<"omega1["<<j<<"] = "<<omega1[j]<<endl;
        }
        myomega<<endl;

        for( j=0; j<(NComp-1); j++ )
        {
            myomega<<"EffRad["<<j<<"] = "<<EffRad[j]<<endl;
        }
        myomega<<endl;

        myomega.close();
#endif

        double epsilon_sum = 0.0;
        // species fractions:
        double spec_sum=0.;
        for( i=0;i<(NComp-1);i++ )
        {
            if( z[i] != 0 )
            {
                spec_sum += m[i];
            }
        }
        for( i=0;i<(NComp-1);i++ )
        {
            spec_frac[i]=m[i]/spec_sum;
        }
        int nWEps=0;
        for( i=0;i<(NComp-1);i++ )
        {
            for( j=i+1;j<(NComp-1);j++ )
            {
                epsilon_sum += (spec_frac[i]+spec_frac[j]) * (*(*(WEps+i)+j));
                if( WEps[i][j]>0.0000001 ){ nWEps++; };
            }
        }

        epsilon_sum = epsilon_sum * IS / nWEps;
//	epsilon_sum = epsilon_sum / (2.302585092994046*R*Tk) *IS;

        for( j=0; j<(NComp-1); j++ )
        {
            ELVIS_lnGam_Born[j] = (omega[j]*epsilon_sum);
        }


        delete[] spec_frac;
        delete[] omega;
        delete[] omega1;
}*/



double TELVIS::Int_OsmCoeff()
{
        double osm_coeff   = 0.0;
        double m_infdil    = 1.1e-6;
        double bjerrum     = 0.0;
        // calc osmotic coefficient via Bjerrum relation
        long int j = 0;

        // bjerrum : \int_{m_k=0}^{m_k=m[j]} m_k  d log_gam(m_k)
        bjerrum = qsimp( m_infdil, m[j], j, 0 );
//		cout<<" integral between "<<m_infdil<<" and "<<m[j]<<" is : "<<bjerrum<<endl;
        osm_coeff = 1 + bjerrum/m[j];
        cout<<" Osmotic coefficient = "<<osm_coeff<<endl;
return osm_coeff;
}


double TELVIS::App_molar_volume()
{
        // add to app_molar_vol_part the stst molar volume of the electrolyte (outside ELVIS class) !!!!
        long int j      = 0;
        double m_infdil = 1.1e-6;
        double m_el 	= m[j]; 	// concentration of electrolyte
        double result, error;

        double app_molar_vol_part = R_CONST*Tk* qsimp( m_infdil, m[j], j, 1 )/m_el;

        // partial molar excess volume
        // double part_molar_excess_vol = R_CONST*Tk*FinDiffVol( m[j], j );
        //cout<<"part_molar_excess_vol = "<<part_molar_excess_vol<<endl;

return app_molar_vol_part; // To get the apparent molar volume of the solute, add the stst partial molar volume to 'app_molar_vol_part'
}



double TELVIS::FinDiff( double m_j, int j )
{
        double h = 1e-4;
        double DactDm;
        double m_old1 = m[j];
        double m_old2 = m[j+1];
/*	// Central Finite Difference
cout<<"m["<<j<<"] = "<<m[j]<<endl;

        m[j] = m_j - h;
cout<<"m["<<j<<"] = "<<m[j]<<endl;
        CalcAct();
        act_low = lnGamma[j];
        m[j] = m_j + h;

cout<<"m["<<j<<"] = "<<m[j]<<endl;
        CalcAct();
        DactDm = (lnGamma[j]-act_low)/(2*h);
*/
        // Forward Finite Difference
        double gam_1,gam_2,gam_3;
//cout<<"FD base 	m["<<j<<"] = "<<m[j]<<endl;
        m[j] 	= m_j;
        m[j+1] 	= m_j;
        molfrac_update();
//cout<<"FD 1 	m["<<j<<"] = "<<m[j]<<endl;
        lnGamma[j] = 0.0;
        CalcAct();
        gam_1 = 0.5*(lnGamma[j]+lnGamma[j+1]);

        m[j] 	= m_j + h;
        m[j+1]	= m_j + h;
        molfrac_update();
        lnGamma[j] = 0.0;
        lnGamma[j] = 0.0;
//cout<<"FD 2		m["<<j<<"] = "<<m[j]<<endl;
        CalcAct();
        gam_2 = 0.5*(lnGamma[j]+lnGamma[j+1]);

        m[j] 	= m_j + h + h;
        m[j+1]	= m_j + h + h;
        molfrac_update();
        lnGamma[j] 		= 0.0;
        lnGamma[j+1] 	= 0.0;
//cout<<"FD 3 	m["<<j<<"] = "<<m[j]<<endl;
        CalcAct();
        gam_3 = 0.5*(lnGamma[j]+lnGamma[j+1]);

        DactDm = (-3*gam_1 + 4*gam_2 - gam_3)/(2*h);

        m[j] 	= m_old1;
        m[j+1] 	= m_old2;
return DactDm * m_j; // for the Bjerrum relation, multiply the derivative with m_j (integrand)
}

// partial molar excess volume of solute
double TELVIS::FinDiffVol( double m_j, int j )
{
    //int j = *(int *) params;

    double lnGam_high, lnGam_low, FinDiff_cation, FinDiff_anion, FinDiffVolume = 0.0;
    double P_diff = Pbar*0.1;	// pressure in bar
    double P_old  = Pbar;
    int stoic_cation = 1;               // stoichiometric coefficient of cation in electrolyte
    int stoic_anion  = 2;               // stoichiometric coefficient of anion in electrolyte

    m[j] 	= stoic_cation*m_j;
    m[j+1] 	= stoic_anion*m_j;
    molfrac_update();


    // Cation contribution
    Pbar 			= P_old + P_diff;
    PTparam();
    CalcAct();
    lnGam_high = lnGamma[j];

    Pbar			= P_old - P_diff;
    PTparam();
    CalcAct();
    lnGam_low	= lnGamma[j];

    FinDiff_cation 		= (lnGam_high - lnGam_low)/(2*P_diff);
        //cout<<"lnGam_high = "<<lnGam_high<<" | lnGam_low = "<<lnGam_low<<endl;


    // Anion contribution
    Pbar 			= P_old + P_diff;
    PTparam();
    CalcAct();
    lnGam_high = lnGamma[j+1];

    Pbar			= P_old - P_diff;
    PTparam();
    CalcAct();
    lnGam_low	= lnGamma[j+1];

    FinDiff_anion 		= (lnGam_high - lnGam_low)/(2*P_diff);
        //cout<<"FinDiff_anion = "<<FinDiff_anion<<endl;

    // Sum
        //FinDiff = (stoic_cation+stoic_anion)*(stoic_cation*FinDiff_cation + stoic_anion*FinDiff_anion);
    FinDiffVolume = (stoic_cation*FinDiff_cation + stoic_anion*FinDiff_anion);

    Pbar  = P_old;
    PTparam();

    return FinDiffVolume;
}


void TELVIS::molfrac_update()
{
    double sum=0.0; int i;
    for( i=0; i<NComp; i++ )
    {    sum += m[i];   }

    for( i=0; i<NComp; i++ )
    {    x[i] = m[i]/sum;  }
}

// Numerical Integration code from Numerical recipes in C (4th edition).
double TELVIS::trapzd( const double m_infdil, const double m_j, int& n, long int& species, int select)
{
        double x,tnm,sum,del = 0.0;
        static double s;
        int it,j = 0;

        // select = 0: Finite diefferences over molality of electrolyte -> compute osmotic coefficient
        if (select == 0)
        {
                if( n == 1 )
                {
                        return (s=0.5*(m_j-m_infdil)*(FinDiff(m_infdil,species)+FinDiff(m_j,species)));
                }
                else
                {
                        for( it=1,j=1;j<n-1;j++ ) it <<= 1;
                        tnm=it;
                        del=(m_j-m_infdil)/tnm;
        //		This is the spacing of the points to be added.
                        x=m_infdil+0.5*del;
                        for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FinDiff(x,species);
                        s=0.5*(s+(m_j-m_infdil)*sum/tnm);
        //		This replaces s by its refined value.
                                return s;
                }
        }
        // select = 1: Finite differences over pressure -> compute apparent molar volume
        else if (select == 1)
        {

                if (n == 1)
                {
                        double lower_bound = FinDiffVol(m_infdil,species);
                        double upper_bound = FinDiffVol(m_j,species);
                        return (s=0.5*(m_j-m_infdil)*(lower_bound+upper_bound));
                }
                else
                {
                        for (it=1,j=1;j<n-1;j++) it <<= 1;
                        tnm=it;
                        del=(m_j-m_infdil)/tnm;
        //		This is the spacing of the points to be added.
                        x=m_infdil+0.5*del;
                        for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FinDiffVol(x,species);
                        s=0.5*(s+(m_j-m_infdil)*sum/tnm);
        //		This replaces s by its refined value.
                                return s;
                }
        }
        return 0.;
}


//	Returns the integral of the function func from a to b. The parameters EPS can be set to the
//	desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum allowed
//	number of steps. Integration is performed by Simpson’s rule.
// Numerical Integration code from Numerical recipes in C (4th edition).
double TELVIS::qsimp(const double m_infdil, const double m_j, long int& species, int select)
{
        int k;
        double s,st,ost,os,EPS,JMAX;
        EPS = 1.0e-6;
        JMAX = 15;
        ost = os = -1.0e30;

        for( k=1;k<=JMAX;k++ )
        {
                st=trapzd(m_infdil,m_j,k,species,select);
                if( std::isnan(st) ) break;

        s=(4.0*st-ost)/3.0;
                        if( k > 6 )
                //		Avoid spurious early convergence.
                                if( fabs(s-os) < EPS*fabs(os) || (s == 0.0 && os == 0.0) ) return s;
                os=s;
                ost=st;
        }
        cout<<"Too many steps in routine qsimp"<<endl;

        return 77777777777777777777777.0;
}



// Output of test results into text file (standalone variant only)
void TELVIS::TELVIS_test_out( const char *path, const double M ) const
{
        long int ii;//, c, a, n;

        cout << "Entered TELVIS_test_out() ... " <<endl;

        ofstream fo("ELVIS_gam.dat", ios::app );
        ErrorIf( !fo.good() , "ELVIS_gam.dat", "Fileopen error");
        //fo << "Gamma	lngamDH 	lngamQuac	lngamC	lngamR	lngamBorn"<<endl;
        fo << M <<" "<<(lnGamma[0]+lnGamma[1])/2<<endl;
        fo << M <<" "<<exp((lnGamma[0]+lnGamma[1])/2)<<endl;
        for( ii=0; ii<NComp; ii++ ){
                fo << M <<"  "<< lnGamma[ii] <<"  "<< gammaDH[ii] << "  "<< gammaQUAC[ii] <<"  "<< gammaC[ii] << "  "<< gammaR[ii] << "  "<<gammaBorn[ii] << endl;
        }
        fo.close();

        //const ios::open_mode OFSMODE = ios::out � ios::app;
        ofstream ff(path, ios::app );
        ErrorIf( !ff.good() , path, "Fileopen error");


        ff << endl << "Debye-H�ckel contribution to Activity Coefficients" << endl;
        for( ii=0; ii<NComp; ii++ )
                ff << gammaDH[ii] << "  ";

        ff << endl << "Born contribution to Activity Coefficients" << endl;
        for( ii=0; ii<(NComp-1); ii++ )
                ff << gammaBorn[ii] << "  ";

        ff << endl << "QUAC contribution to Activity Coefficients" << endl;
        for( ii=0; ii<NComp; ii++ )
                ff << gammaQUAC[ii] << "  ";

        ff << endl << "Contribution of Combinatorial Term to Activity Coefficients" << endl;
        for( ii=0; ii<NComp; ii++ )
                ff << gammaC[ii] << "  ";

        ff << endl << "Contribution of Residual Term to Activity Coefficients" << endl;
        for( ii=0; ii<NComp; ii++ )
                ff << gammaR[ii] << "  ";

        ff << endl << "ln activity coefficients of end members" << endl;
        for( ii=0; ii<NComp; ii++ )
                ff << lnGamma[ii] << "  ";

        ff << endl << "Activity coefficients of end members" << endl;
        for( ii=0; ii<NComp; ii++ )
                ff << exp(lnGamma[ii]) << "  ";
        ff << endl;

        ff.close();
}

//--------------------- End of s_solmod4.cpp ---------------------------


