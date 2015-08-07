//-------------------------------------------------------------------
// $Id: s_fgl4.cpp 724 2012-10-02 14:25:25Z kulik $
//
/// \file s_solmod5.cpp
/// Implementation of TSolMod derived classes for specific ion interaction
///  aqueous activity models (TSIT, TPitzer, TEUNIQUAC)
//
// Copyright (c) 2008-2012  F.Hingerl, T.Wagner, D.Kulik
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
// SIT model (NEA version) reimplementation for aqueous electrolyte solutions
// References:
// (c) DK/TW June 2009
//=============================================================================================


// Generic constructor for the TSIT class
TSIT::TSIT( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW ):
                TSolMod( sd )
{
    alloc_internal();
    z = arZ;
    m = arM;
    RhoW = dW;
    EpsW = eW;
}


TSIT::~TSIT()
{
    free_internal();
}


void TSIT::alloc_internal()
{
    LnG = new double [NComp];
    dLnGdT = new double [NComp];
    d2LnGdT2 = new double [NComp];
    dLnGdP = new double [NComp];
    E0 = new double *[NComp];
    E1 = new double *[NComp];
    dE0 = new double *[NComp];
    dE1 = new double *[NComp];
    d2E0 = new double *[NComp];
    d2E1 = new double *[NComp];

    for (long int j=0; j<NComp; j++)
    {
        E0[j] = new double [NComp];
        E1[j] = new double [NComp];
        dE0[j] = new double [NComp];
        dE1[j] = new double [NComp];
        d2E0[j] = new double [NComp];
        d2E1[j] = new double [NComp];
    }
}


void TSIT::free_internal()
{
    for (long int j=0; j<NComp; j++)
    {
        delete[]E0[j];
        delete[]E1[j];
        delete[]dE0[j];
        delete[]dE1[j];
        delete[]d2E0[j];
        delete[]d2E1[j];
    }
    delete[]LnG;
    delete[]dLnGdT;
    delete[]d2LnGdT2;
    delete[]dLnGdP;
    delete[]E0;
    delete[]E1;
    delete[]dE0;
    delete[]dE1;
    delete[]d2E0;
    delete[]d2E1;
}


long int TSIT::PTparam()
{
    long int j, i, ip, i1, i2;
    double p0, p1, alp, bet, dal, rho, eps, dedt, d2edt2, dedp;

    // fill internal arrays of interaction parameters with standard value
    for (j=0; j<NComp; j++)
    {
        for (i=0; i<NComp; i++)
        {
            E0[j][i] = 0.0;
            E1[j][i] = 0.0;
            dE0[j][i] = 0.0;
            dE1[j][i] = 0.0;
            d2E0[j][i] = 0.0;
            d2E1[j][i] = 0.0;
        }
    }

    // read and convert interaction parameters that have non-standard value
    for (ip=0; ip<NPar; ip++)
    {
        i1 = aIPx[MaxOrd*ip];
        i2 = aIPx[MaxOrd*ip+1];
        p0 = aIPc[NPcoef*ip];
        p1 = aIPc[NPcoef*ip+1];
        E0[i1][i2] = p0;
        E1[i1][i2] = p1;
        E0[i2][i1] = p0;
        E1[i2][i1] = p1;
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

    // calculate A term of Debye-Huckel equation (and derivatives)
    A = (1.82483e6)*sqrt(rho) / pow(Tk*eps,1.5);
    dAdT = - 3./2.*A*( dedt + 1./Tk + alp/3. );
    d2AdT2 = 1./A*pow(dAdT,2.) - 3./2.*A*( d2edt2 - 1/pow(Tk,2.) + 1/3.*dal );
    dAdP = 1./2.*A*( bet - 3.*dedp);

    ErrorIf( fabs(A) < 1e-9, "SIT model",
            "Error: DH parameter A was not calculated - no values of RhoW and EpsW !" );

    return 0;
}


/// Calculates activity coefficients in SIT (NEA) model
long int TSIT::MixMod()
{
    long int j, i1, i2, ip;
    double sqI, lgI, Z2, lgGam, SumSIT, lg_to_ln;
    lg_to_ln = 2.302585093;

    I = IonicStrength();
    sqI = sqrt(I);
    lgI = log10(I);

    // this check was already performed in CalculateActivityCoefficients()
    if( I < 1e-6 )
    {
        for( j=0; j<NComp; j++)
            lnGamma[j] = 0.;
        return 0;
    }

    // loop over species
    for( j=0; j<NComp; j++ )
    {
        lgGam = 0.;

        // Calculation of the SIT sum (new variant)
        // Corrected to 2-coeff SIT parameter and extended to neutral species by DK on 13.05.2009
        SumSIT = 0.;
        if( j != NComp-1 )  // not for water solvent
        {
            for( ip=0; ip<NPar; ip++ )
            {
                i1 = aIPx[ip*MaxOrd];  // order of indexes for binary parameters plays no role
                i2 = aIPx[ip*MaxOrd+1];

                if( i1 == i2 )
                    continue;

                if( i1 == j )
                {
                    // SumSIT += ( aIPc[ip*NPcoef] + aIPc[ip*NPcoef+1]*lgI ) * m[i2]; // epsilon
                    SumSIT += ( E0[i1][i2] + E1[i1][i2]*lgI ) * m[i2];  // epsilon
                }

                else if( i2 == j )
                {
                    // SumSIT += ( aIPc[ip*NPcoef] + aIPc[ip*NPcoef+1]*lgI ) * m[i1]; // epsilon
                    SumSIT += ( E0[i1][i2] + E1[i1][i2]*lgI ) * m[i1];  // epsilon
                }
            }
        }

        // Charged species
        if( z[j] )
        {
            lgGam = 0.;
            Z2 = z[j]*z[j];
            lgGam = ( - A * sqI * Z2 ) / ( 1. + 1.5 * sqI );  // DH part for charged species
            lgGam += SumSIT;
            lnGamma[j] = lgGam * lg_to_ln;
        }

        else // neutral species and water solvent
        {
            // neutral species
            if ( j != (NComp-1) )
            {
                lgGam = 0.;
                Z2 = 0.;
                lgGam += SumSIT;
                lnGamma[j] = lgGam * lg_to_ln;
            }

            // water solvent (osmotic coefficient not yet considered)
            else
            {
                lgGam = 0.;
                lnGamma[j] = lgGam * lg_to_ln;
            }
        }
    } // j

    return 0;
}


long int TSIT::ExcessProp( double *Zex )
{
    // (under construction)
    long int j, i1, i2, ip;
    double sqI, lgI, Z2, SumSIT, lg_to_ln, g, dgt, d2gt, dgp;
    lg_to_ln = 2.302585093;
    g = 0.; dgt = 0.; d2gt = 0.; dgp = 0.;

    I = IonicStrength();
    sqI = sqrt(I);
    lgI = log10(I);

    // loop over species
    for( j=0; j<NComp; j++ )
    {
        // Calculation of the SIT sum (new variant)
        // Corrected to 2-coeff SIT parameter and extended to neutral species by DK on 13.05.2009
        SumSIT = 0.;
        if( j != NComp-1 )  // not for water solvent
        {
            for( ip=0; ip<NPar; ip++ )
            {
                i1 = aIPx[ip*MaxOrd];  // order of indexes for binary parameters plays no role
                i2 = aIPx[ip*MaxOrd+1];

                if( i1 == i2 )
                    continue;

                if( i1 == j )
                {
                    // SumSIT += ( aIPc[ip*NPcoef] + aIPc[ip*NPcoef+1]*lgI ) * m[index2]; // epsilon
                    SumSIT += ( E0[i1][i2] + E1[i1][i2]*lgI ) * m[i2]; // epsilon
                }

                else if( i2 == j )
                {
                    // SumSIT += ( aIPc[ip*NPcoef] + aIPc[ip*NPcoef+1]*lgI ) * m[index1]; // epsilon
                    SumSIT += ( E0[i1][i2] + E1[i1][i2]*lgI ) * m[i1]; // epsilon
                }
            }
        }

        // Charged species (no TP dependence of SIT parameters)
        if( z[j] )
        {
            Z2 = z[j]*z[j];
            LnG[j] = - ( ( A * sqI * Z2 ) / ( 1. + 1.5 * sqI ) + SumSIT ) * lg_to_ln;
            dLnGdT[j] = - ( ( dAdT * sqI * Z2 ) / ( 1. + 1.5 * sqI ) ) * lg_to_ln;
            d2LnGdT2[j] = - ( ( d2AdT2 * sqI * Z2 ) / ( 1. + 1.5 * sqI ) ) * lg_to_ln;
            dLnGdP[j] = - ( ( dAdP * sqI * Z2 ) / ( 1. + 1.5 * sqI ) ) * lg_to_ln;
        }

        // neutral species and water solvent
        else
        {
            // neutral species
            if ( j != (NComp-1) )
            {
                LnG[j] = ( SumSIT ) * lg_to_ln;
                dLnGdT[j] = 0.;
                d2LnGdT2[j] = 0.;
                dLnGdP[j] = 0.;
            }

            // water solvent (osmotic coefficient not yet considered)
            else
            {
                LnG[j] = 0.;
                dLnGdT[j] = 0.;
                d2LnGdT2[j] = 0.;
                dLnGdP[j] = 0.;
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
long int TSIT::IdealProp( double *Zid )
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
double TSIT::IonicStrength()
{
    double Is = 0.;

    for( long int ii=0; ii<(NComp-1); ii++ )
        Is += 0.5*( z[ii]*z[ii]*m[ii] );

    return Is;
}



//=============================================================================================
// Pitzer model for aqueous electrolyte solutions, Harvie-Moller-Weare (HMW) version
// References: Zhang et al. (2006)
// (c) FH/SD December 2008
//=============================================================================================


// Generic constructor for the TPitzer class
TPitzer::TPitzer( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW ):
                TSolMod( sd )
{
    aZ = arZ;
    aM = arM;
    RhoW = dW;
    EpsW = eW;

    Aphi = dAphidT = d2AphidT2 = I = Is = Ffac = Zfac = 0.0;

    // calculate sizes Nc, Na, Nn, Ns
    calcSizes();

    // reallocate internal arrays and set zeros
    alloc_internal();
}


TPitzer::~TPitzer()
{
    free_internal();
}


void TPitzer::alloc_internal()
{
    long int i,j;
    long int ic, ia, in, ic1, ia1;

    // Input parameter arrays
    // 1D objects
    xcx = new long int[Nc];
    xax = new long int[Na];
    xnx = new long int[Nn];
    mc = new double[Nc];
    ma = new double[Na];
    mn = new double[Nn];
    zc = new double[Nc];
    za = new double[Na];

    //2D objects
    Bet0	= new double *[Nc];		// c,a
    Bet1	= new double *[Nc];		// c,a
    Bet2	= new double *[Nc];		// c,a
    Cphi	= new double *[Nc];		// c,a
    Theta	= new double *[Nc];		// c,c1
    Theta1	= new double *[Na];		// a,a1
    Lam		= new double *[Nn];	    // n,c
    Lam1	= new double *[Nn];		// n,a

    // coefficient[Nc][] allocations
    for(i=0; i<Nc ; i++)
    {
        Bet0[i]   = new double[Na];
        Bet1[i]   = new double[Na];
        Bet2[i]   = new double[Na];
        Cphi[i]   = new double[Na];
        Theta[i]  = new double[Nc];
    }

    // coefficient[Nn][] allocations
    for(i=0; i<Nn; i++)
    {
        Lam[i]   = new double[Nc];
        Lam1[i]  = new double[Na];
    }

    // coefficient[Na][] allocations
    for( i=0; i<Na ; i++)
    {
        Theta1[i] = new double[Na];
    }

    //3D objects
    Psi  = new double **[Nc];		// c,c1,a
    Psi1 = new double **[Na];		// a,a1,c
    Zeta = new double **[Nn];		// n,c,a

    //coefficient[Nc][Nc][Na] allocations
    for(i=0; i<Nc; i++)
    {
        Psi[i]   = new double *[Nc];
            for(j=0; j<Nc; j++)
            {
                Psi[i][j] = new double [Na];
            }
    }

    //coefficient[Na][Na][Nc] allocations
    for(i=0; i<Na; i++)
    {
        Psi1[i] = new double *[Na];
        for(j=0; j<Na; j++)
        {
            Psi1[i][j] = new double [Nc];
        }
    }

    //coefficient[Nn][Nc][Na] allocations
    for(i=0; i<Nn; i++)
    {
        Zeta[i]  = new double *[Nc];
        for(j=0; j<Nc; j++)
        {
            Zeta[i][j] = new double[Na];
        }
    }

    //MacInnes Parameter Array for binary KCl solution
    McI_PT_array = new double[13];

    //MacInnes Activity coefficient array
    long int alle = Ns + 1;
    GammaMcI = new double[alle];


    // initialize work objects with standard values

    for(ic=0; ic<Nc; ic++)
    {
        for(ic1=0; ic1<Nc; ic1++)
        {
            Theta[ic][ic1]	= 0.0;
            for(ia=0; ia<Na; ia++)
            {
                Psi[ic][ic1][ia] = 0.0;
            }
        }
    }

    for(ia=0; ia<Na; ia++)
    {
        for(ia1=0; ia1<Na; ia1++)
        {
            Theta1[ia][ia1]	= 0.0;
            for(ic=0; ic<Nc; ic++)
            {
                Psi1[ia][ia1][ic] = 0.0;
            }
        }
    }

    for(in=0; in<Nn; in++)
    {
        for(ic=0; ic<Nc; ic++)
        {
            Lam[in][ic] = 0.0;
            for(ia=0; ia<Na; ia++)
            {
                Bet0[ic][ia]     = 0.0;
                Bet1[ic][ia]     = 0.0;
                Bet2[ic][ia]     = 0.0;
                Cphi[ic][ia]     = 0.0;
                Lam1[in][ia]     = 0.0;
                Zeta[in][ic][ia] = 0.0;
            }
        }
    }

    for(ic=0; ic<13; ic++)
    {
        McI_PT_array[ic] = 0.0;
    }

    for(ic=0; ic<alle; ic++)
    {
        GammaMcI[ic] = 0.0;
    }
}


void TPitzer::free_internal()
{
    long int i, j;

    // delete internal work objects
    for(i=0; i<Nc ; i++)
    {
        delete[]Bet0[i];
        delete[]Bet1[i];
        delete[]Bet2[i];
        delete[]Cphi[i];
        delete[]Theta[i];
    }

    for( i=0; i<Na ; i++)
    {
        delete[]Theta1[i];
    }

    delete[]Bet0;
    delete[]Bet1;
    delete[]Bet2;
    delete[]Cphi;
    delete[]Theta;
    delete[]Theta1;

    for(i=0; i<Na; i++)
    {
        for(j=0; j<Na; j++)
        {
            delete[]Psi1[i][j];
        }
    }
    for(i=0; i<Na; i++)
    {
        delete[]Psi1[i];
    }
    delete[]Psi1;

    for(i=0; i<Nc; i++)
    {
        for(j=0; j<Nc; j++)
        {
            delete[]Psi[i][j];
        }
    }
    for(i=0; i<Nc; i++)
    {
            delete[]Psi[i];
    }
    delete[]Psi;



    if(Nn != 0)
    {

        for(i=0; i<Nn; i++)
        {
            for(j=0; j<Nc; j++)
            {
                delete[]Zeta[i][j];
            }
        }
        for(i=0; i<Nn; i++)
        {
            delete[]Zeta[i];
        }
        delete[]Zeta;


        for( i=0; i<Nn ; i++)
        {
            delete[]Lam[i];
        }
        delete[]Lam;

        for( i=0; i<Nn ; i++)
        {
            delete[]Lam1[i];
        }
        delete[]Lam1;


    }

    delete[]xcx;
    delete[]xax;
    delete[]xnx;
    delete[]mc;
    delete[]ma;
    delete[]mn;
    delete[]zc;
    delete[]za;
    delete[]McI_PT_array;
    delete[]GammaMcI;
}


/// Output of test results into text file (standalone variant only)
void TPitzer::Pitzer_test_out( const char *path, double Y )
{

//	long int ii, c, a, n;

    ofstream ff(path, ios::app );
    ErrorIf( !ff.good() , path, "Fileopen error");
    int number_after_decimalpoint = 4;
    ff.setf(ios::fixed);
    ff.setf(ios::showpoint);
    ff.setf(ios::showpos);
    ff.precision(number_after_decimalpoint);
    ff <<" "<<exp(lnGamma[0])<<" "<<exp(lnGamma[1])<<" "<<exp(lnGamma[2])<<" "<<exp(lnGamma[3])<<" "<<exp(lnGamma[4])<<" "<<Y<<endl;
    ff << endl;
}


long int TPitzer::PTparam( )
{
    // build conversion of species indexes between aqueous phase and Pitzer parameter tables
    setIndexes();

    // calculate vector of interaction parameters corrected to T,P of interest
    PTcalc( 0 );

    return 0;
}


long int TPitzer::MixMod()
{
    // Refreshing molalities of cations, anions, and neutral species
    // Added by DK on 01.Dec.2009
    long int ic, ia, in, iRet;

    for( ic=0; ic<Nc; ic++){
        mc[ic] = aM[xcx[ic]];
        zc[ic] = aZ[xcx[ic]];
      }

    for(ia=0; ia<Na; ia++){
        ma[ia] = aM[xax[ia]];
        za[ia] = aZ[xax[ia]];
     }

    for(in=0; in<Nn; in++){   //Important bugfix (was in>Nn) 17/12/2012 SD
        mn[in] = aM[xnx[in]];
    }

    iRet = Pitzer_calc_Gamma();
    return iRet;
}


/// calculates ideal mixing properties
long int TPitzer::IdealProp( double *Zid )
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


/// Calculation of activity coefficients
long int TPitzer::Pitzer_calc_Gamma( )
{
    long int M, N, X;

    // Ionic Strength
    Is = IonicStr( I );

    // F-Factor, Pitzer-Toughreact Report 2006 equation (A6)
    Ffac = F_Factor( Aphi );

    // Z- Term, Pitzer-Toughreact Report 2006 equation (A8)
    Zfac = Z_Term();

    lnGamma[Ns] = lnGammaH2O( Aphi );

    for( M=0; M<Nc; M++ )
    {
        lnGamma[xcx[M]] = lnGammaM( M, Aphi );
    }

    for( X=0; X<Na; X++ )
    {
        lnGamma[xax[X]] = lnGammaX( X, Aphi );
    }

    if( Nn > 0 )
    {
        for( N=0; N<Nn; N++ )
        {
                    lnGamma[xnx[N]] = lnGammaN( N );
            }
    }
    //Pitzer_test_out( "test111.dat ");

    return 0;
}


/// calculation of activity coefficient of KCl
long int TPitzer::Pitzer_McInnes_KCl( )
{
    McInnes_KCl( );
    return 0;
}


/// Build conversion of species indexes.
/// list of indexes of Nc cations in aqueous phase
/// list of indexes of Na anions in aq phase
/// list of indexes of Nn neutral species in aq phase
void TPitzer::setIndexes()
{
    long int jj, ic, ia, in;

    ic = ia = in = 0;

    for( jj=0; jj<NComp-1; jj++ ) // -1 = no check H2O
    {
        if( aZ[jj] > 0 )
        {
            xcx[ic] = jj;
            ic++;
        }

        else
            if( aZ[jj] < 0 )
            {
                xax[ia] = jj;
                ia++;
            }
            else
            {
                xnx[in] = jj;
                in++;
            }
    }

    for( ic=0; ic<Nc; ic++){
    mc[ic] = aM[xcx[ic]];
    zc[ic] = aZ[xcx[ic]];
    }

    for(ia=0; ia<Na; ia++){
    ma[ia] = aM[xax[ia]];
    za[ia] = aZ[xax[ia]];
    }

    for(in=0; in>Nn; in++){
    mn[in] = aM[xnx[in]];
    }

}


/// Pitzer parameters at T of interest (5 term)
double TPitzer::G_ex_par5(long int ii)
{
    double Tr, G_ex_par;

    Tr = 298.15;
    G_ex_par = aIPc[ii * NPcoef + 0] + aIPc[ii * NPcoef + 1]*(1./Tk-1./Tr) +
                    aIPc[ii * NPcoef + 2]*log(Tk/Tr) + aIPc[ii * NPcoef + 3]*(Tk-Tr) +
                    aIPc[ii * NPcoef + 4]*(Tk*Tk-Tr*Tr);
    return G_ex_par;
}


/// Pitzer parameters at T of interest (5 term)
double TPitzer::G_ex_par8(long int ii)
{
    double G_ex_par;

    G_ex_par = aIPc[ii * NPcoef + 0] + aIPc[ii * NPcoef + 1]*Tk + aIPc[ii * NPcoef + 2]/Tk +
                    aIPc[ii * NPcoef + 3]*log(Tk) + aIPc[ii * NPcoef + 4]/(Tk-263.) +
                    aIPc[ii * NPcoef + 5]*Tk*Tk + aIPc[ii * NPcoef + 6]/(680.-Tk) +
                    aIPc[ii * NPcoef + 7]/(Tk-227.);
    return G_ex_par;
}


/// first T derivative of Pitzer parameters (5-term)
double TPitzer::S_ex_par5(long int ii)
{
    double S_ex_par;

    S_ex_par = (-1.) * aIPc[ii * NPcoef + 1]*pow(Tk, (-2.)) + aIPc[ii * NPcoef + 2]*pow(Tk, (-1.)) +
                    aIPc[ii * NPcoef + 3] + 2. * aIPc[ii * NPcoef + 4]*Tk;

    return S_ex_par;
}


/// first T derivative of Pitzer parameters (8-term)
double TPitzer::S_ex_par8(long int ii)
{
    double S_ex_par;

    S_ex_par = aIPc[ii * NPcoef + 1] - aIPc[ii * NPcoef + 2] * pow(Tk, (-2.)) + aIPc[ii * NPcoef + 3] * pow(Tk, (-1.)) -
                    aIPc[ii * NPcoef + 4] * pow((Tk -263.), (-2.)) + 2. * aIPc[ii * NPcoef + 5] * Tk -
                    aIPc[ii * NPcoef + 6] * pow((680. - Tk), (-2.)) - aIPc[ii * NPcoef + 7] * pow((Tk - 227.), (-2.));

    return S_ex_par;
}


/// second T derivative of Pitzer parameters (5-term)
double TPitzer::CP_ex_par5(long int ii)
{
    double CP_ex_par;

    CP_ex_par = 2. * aIPc[ii * NPcoef + 1] *pow(Tk, (-3.)) - aIPc[ii * NPcoef + 2] *pow(Tk, (-2.)) +
                        2. * aIPc[ii * NPcoef + 4];

    return CP_ex_par;
}


/// second T derivative of Pitzer parameters (8-term)
double TPitzer::CP_ex_par8(long int ii)
{
    double CP_ex_par;

    CP_ex_par = 2. * aIPc[ii * NPcoef + 2] * pow(Tk, (-3.)) - aIPc[ii * NPcoef + 3] * pow(Tk, (-2.))
                    + 2. * aIPc[ii * NPcoef + 4] * pow((Tk - 263.), (-3.)) + 2. * aIPc[ii * NPcoef + 5] +
                    2. * aIPc[ii * NPcoef + 6] * pow((680. -Tk), (-3.)) + 2. * aIPc[ii * NPcoef + 7] * pow((Tk - 227.),(-3.));

    return CP_ex_par;
}


double TPitzer::setvalue(long int ii, int Gex_or_Sex)
{
    double value;
    //   int Gex_or_Sex == 0 -> Gex (unchanged formula for PT correction of interaction params used)
    //	  int Gex_or_Sex == 1 -> Sex (first temperature derivative)
    //    int Gex_or_Sex == 2 -> CPex (second temperature derivative)
    if(Gex_or_Sex==0)
    {
        if( NPcoef == 5 )
        {
            value = G_ex_par5(ii);
        }
        else if( NPcoef == 8 ){
            value = G_ex_par8(ii);
        }
        else Error( "", "PitzerHMW: Invalid number of coefficients to describe T dependence");
    }
    else if(Gex_or_Sex==1)
    {
        if ( NPcoef == 5)
        {
            value = S_ex_par5(ii);
        }
        else if( NPcoef == 8)
        {
            value = S_ex_par8(ii);
        }
        else Error( "", "PitzerHMW: Invalid number of coefficients to describe T dependence");
    }
    else if(Gex_or_Sex==2)
    {
        if ( NPcoef == 5)
        {
            value = CP_ex_par5(ii);
        }
        else if( NPcoef == 8)
        {
        value = CP_ex_par8(ii);
        }
        else Error( "", "PitzerHMW: Invalid number of coefficients to describe T dependence");
    }
    else Error ( "", "PitzerHMW: Only first and second temperature derivatives of activity function are implemented" );

    return value;
}


/// Copy data from arIPx, arIPc, arDCc to internal structures
void TPitzer::PTcalc( int Gex_or_Sex )
{
    // if G_or_ex == 1, then the parameters used in the activity coefficient calculation
    // are given the values of aIP (-> for scaled and unscaled activity coeffs). If G_or_ex
    // not equal 1,then the parameters for excess property calculation will be written into the
    // the interaction parameters
    long int ii, ic, ia, in, i;
    double Tr = 298.15;
    double N0 = 6.0221415e23;      // Avogadro
    double k = 1.3806505e-23;      // Boltzmann
    double el = 1.60217635e-19;    // Coulomb charge unit
    double eps0 = 8.854187817e-12; // Dielectricity constant vacuum
    double pi = 3.141592654;
    double Rho, Eps, const_term;

    double drho_term1, drho_term2, drho_term3, drho_term4,
                deps_term1, deps_term2, deps_term3, deps_term4,
                dTk_term1, dTk_term2, dTk_term3;


    // PT correction of interaction parameters
    for( ii=0; ii<NPar; ii++ )
    {

        switch( aIPx[ii * MaxOrd + 3] )  // type of table
        {
            case bet0_:
                ic = getIc( aIPx[ii * MaxOrd + 0] );
                if( ic < 0 )
                {
                    ic = getIc( aIPx[ii * MaxOrd + 1] );
                    ia = getIa( aIPx[ii * MaxOrd + 0] );
                }
                else
                    ia = getIa( aIPx[ii * MaxOrd + 1] );
                ErrorIf( ia<0||ic<0, "", "Cation and anion index needed here"  );
                Bet0[ic][ia] = setvalue(ii, Gex_or_Sex);
                break;

            case bet1_:
                ic = getIc( aIPx[ii * MaxOrd + 0] );
                if( ic < 0 )
                {
                    ic = getIc( aIPx[ii * MaxOrd + 1] );
                    ia = getIa( aIPx[ii * MaxOrd + 0] );
                }
                else
                    ia = getIa( aIPx[ii * MaxOrd + 1] );
                ErrorIf( ia<0||ic<0, "", "Cation and anion index needed here"  );
                Bet1[ic][ia] = setvalue(ii, Gex_or_Sex);
                break;

            case bet2_:
                ic = getIc( aIPx[ii * MaxOrd + 0] );
                if( ic < 0 )
                {
                    ic = getIc( aIPx[ii * MaxOrd + 1] );
                    ia = getIa( aIPx[ii * MaxOrd + 0] );
                }
                else
                    ia = getIa( aIPx[ii * MaxOrd + 1] );
                ErrorIf( ia<0||ic<0, "", "Cation and anion indexes needed here"  );

                Bet2[ic][ia] = setvalue(ii, Gex_or_Sex);
                break;

            case Cphi_:
                ic = getIc( aIPx[ii * MaxOrd + 0] );
                if( ic < 0 )
                {
                    ic = getIc( aIPx[ii * MaxOrd + 1] );
                    ia = getIa( aIPx[ii * MaxOrd + 0] );
                }
                else
                    ia = getIa( aIPx[ii * MaxOrd + 1] );
                ErrorIf( ia<0||ic<0, "", "Cation and anion indexes needed here"  );

                Cphi[ic][ia] = setvalue(ii, Gex_or_Sex);
                break;

            case Lam_:
                in = getIn( aIPx[ii * MaxOrd + 0] );
                if( in < 0 )
                {
                    in = getIn( aIPx[ii * MaxOrd + 1] );
                    ic = getIc( aIPx[ii * MaxOrd + 0] );
                }
                else
                    ic = getIc( aIPx[ii * MaxOrd + 1] );
                ErrorIf( in<0||ic<0, "", "Cation and neutral species indexes needed here"  );

                Lam[in][ic] = setvalue(ii, Gex_or_Sex);
                break;

            case Lam1_:
                in = getIn( aIPx[ii * MaxOrd + 0] );
                if( in < 0 )
                {
                    in = getIn( aIPx[ii * MaxOrd + 1] );
                    ia = getIa( aIPx[ii * MaxOrd + 0] );
                }
                else
                    ia = getIa( aIPx[ii * MaxOrd + 1] );
                ErrorIf( in<0||ia<0, "", "Parameters must be anion and neutral species index"  );

                Lam1[in][ia] = setvalue(ii, Gex_or_Sex);
                break;

            case Theta_:
                ic = getIc( aIPx[ii * MaxOrd + 0] );
                i = getIc( aIPx[ii * MaxOrd + 1] );
                ErrorIf( i<0||ic<0, "", "Only indexes of cations needed here"  );

                Theta[ic][i] = setvalue(ii, Gex_or_Sex);
                break;

            case Theta1_:
                ia = getIa( aIPx[ii * MaxOrd + 0] );
                i = getIa( aIPx[ii * MaxOrd + 1] );
                ErrorIf( i<0||ia<0, "", "Only indexes of anions needed here"  );

                Theta1[ia][i] = setvalue(ii, Gex_or_Sex);
                break;

            case Psi_:
                ic = getIc( aIPx[ii * MaxOrd + 0] );
                if( ic < 0 )
                {
                    ic = getIc( aIPx[ii * MaxOrd + 1] );
                    ia = getIa( aIPx[ii * MaxOrd + 0] );
                    i =  getIc( aIPx[ii * MaxOrd + 2] );
                }
                else
                {
                    i =  getIc( aIPx[ii * MaxOrd + 1] );
                    if( i<0 )
                    {
                        ia = getIa( aIPx[ii * MaxOrd + 1] );
                        i =  getIc( aIPx[ii * MaxOrd + 2] );
                    }
                    else
                        ia = getIa( aIPx[ii * MaxOrd + 2] );
                }
                ErrorIf( ic<0||ia<0||i<0, "", "Index of anion and 2 indexes of cations needed here"  );

                Psi[ic][i][ia] = setvalue(ii, Gex_or_Sex);
                break;

            case Psi1_:
                ia = getIa( aIPx[ii * MaxOrd + 0] );
                if( ia < 0 )
                {
                    ia = getIa( aIPx[ii * MaxOrd + 1] );
                    ic = getIc( aIPx[ii * MaxOrd + 0] );
                    i =  getIa( aIPx[ii * MaxOrd + 2] );
                }
                else
                {
                    i =  getIa( aIPx[ii * MaxOrd + 1] );
                    if( i<0 )
                    {
                        ic = getIc( aIPx[ii * MaxOrd + 1] );
                        i =  getIa( aIPx[ii * MaxOrd + 2] );
                    }
                    else
                        ic = getIc( aIPx[ii * MaxOrd + 2] );
                }
                ErrorIf( ic<0||ia<0||i<0, "", "Indexes of 2 anions and one cation needed here"  );

                Psi1[ia][i][ic] = setvalue(ii, Gex_or_Sex);
                break;

            case Zeta_:
                in = getIn( aIPx[ii * MaxOrd + 0] );
                if( in < 0 )
                {
                    in = getIn( aIPx[ii * MaxOrd + 1] );
                    if( in < 0 )
                        in = getIn( aIPx[ii * MaxOrd + 2] );
                }
                ic = getIc( aIPx[ii * MaxOrd + 1] );
                if( ic < 0 )
                {
                    ic = getIc( aIPx[ii * MaxOrd + 2] );
                    if( ic < 0 )
                        ic = getIc( aIPx[ii * MaxOrd + 0] );
                }
                ia = getIa( aIPx[ii * MaxOrd + 2] );
                if( ia < 0 )
                {
                    ia = getIa( aIPx[ii * MaxOrd + 0] );
                    if( ia < 0 )
                        ia = getIa( aIPx[ii * MaxOrd + 1] );
                }
                ErrorIf( ic<0||ia<0||in<0, "",
                        "Index of neutral species, index of cation and index of anion needed here"  );
                Zeta[in][ic][ia] = setvalue(ii, Gex_or_Sex);
                break;
            default:
                break;
        }
    }

    // KCl Parameter correction for MacInnes convention scaling
    // parameters from PHREEQC database pitzer.dat
    double McI_par_array[13][5] =
    {	{0.04835,	0,	0,	5.79E-04,	0},		//PT_B0_KCl
        {0.2122,	0,	0,	1.07E-03,	0},		//PT_B1_KCl
        {-0.00084,	0,	0,	-5.10E-05,	0},		//PT_Cphi_KCl
        {0.1298,	0,	0,	0,	0},				//PT_B0_KOH
        {0.32,	0,	0,	0,	0},					//PT_B1_KOH
        {0.0041,	0,	0,	0,	0},				//PT_Cphi_KOH
        {0.1775,	0,	0,	-3.08E-04,	0},		//PT_B0_HCl
        {0.2945,	0,	0,	1.42E-04,	0},		//PT_B1_HCl
        {0.0008,	0,	0,	6.21E-05,	0},		//PT_Cphi_HCl
        {0.005,	0,	0,	0,	0},					//PT_Theta_KH
        {-0.011,	0,	0,	0,	0},				//PT_Psi_KHCl
        {-0.050,	0,	0,	0,	0},				//PT_Theta_ClOH
        {-0.006,	0,	0,	0,	0}				//PT_Psi_ClOHK
    };

    for( int ii=0; ii<13; ii++ ) {
        McI_PT_array[ii] = McI_par_array[ii][0] + McI_par_array[ii][1]*(1./Tk-1./Tr) + McI_par_array[ii][2]*log(Tk/Tr) +
                        McI_par_array[ii][3]*(Tk-Tr) + McI_par_array[ii][4]*(Tk*Tk-Tr*Tr);
    }

    // Calculation of Debye-H�ckel Parameter and its derivatives
    Rho = RhoW[0];
    Eps = EpsW[0];

    drho_term1 = pow(RhoW[0], 0.5);
    drho_term2 = (0.5* pow(RhoW[0], (-0.5))* RhoW[1]);
    drho_term3 = -0.25 * pow(RhoW[0], (-1.5)) * pow(RhoW[1], 2.);
    drho_term4 = (0.5* pow(RhoW[0], -(0.5))* RhoW[2]);
    deps_term1 = pow(EpsW[0], (-3./2.) );
    deps_term2 = (-3./2. * pow(EpsW[0], (-5./2.)) *EpsW[1]);
    deps_term3 = 15./4. * pow(EpsW[0], (-7./2.)) * pow(EpsW[1], 2.);
    deps_term4 = (-3./2. * pow(EpsW[0], (-5./2.)) *EpsW[2]);
    dTk_term1 = pow(Tk, -3./2.);
    dTk_term2 = (-3./2. * pow(Tk, (-5./2.)));
    dTk_term3 = 15./4. * pow(Tk, (-7./2.));

    // Aphi
    Aphi = (1./3.) * pow((2.*pi*N0*Rho*1000.),0.5) * pow((el*el)/(Eps*4.*pi*eps0*k*Tk),1.5);
    ErrorIf( fabs(Aphi) < 1e-9, "Pitzer model",
            "Error: DH parameter A was not calculated - no values of RhoW and EpsW !" );

    // dAphidT
    const_term = 1./3. * pow((2.*pi*N0*1000.),0.5) * pow(el*el/(4.*pi*eps0*k),1.5);
    dAphidT = const_term*( (0.5* pow(RhoW[0], 0.5)* RhoW[1]) * pow(EpsW[0], -3./2. )                * pow(Tk, -3./2.) +
                            pow(RhoW[0], 0.5)              * (-3./2. * pow(EpsW[0], -5./2.) *EpsW[1]) * pow(Tk, -3./2.) +
                            pow(RhoW[0], 0.5)              * pow(EpsW[0], -3./2.)                * (-3./2. *pow(Tk, -5./2.)) );

    // d2AphidT2
    d2AphidT2 = const_term * (
            ( (drho_term3 + drho_term4) * deps_term1 * dTk_term1 +
                drho_term2 * deps_term2 * dTk_term1 +
                drho_term2 * deps_term1 * dTk_term2 ) 			+
            ( drho_term2 * deps_term2 * dTk_term1 +
                drho_term1 * (deps_term3 + deps_term4) * dTk_term1 +
                drho_term1 * deps_term2 * dTk_term2 )			+
            ( drho_term2 * deps_term1 * dTk_term2 +
                drho_term1 * deps_term2 * dTk_term2 +
                drho_term1 * deps_term1 * dTk_term3 )  				);

}


void TPitzer::calcSizes()
{
    long int jj;

    Nc = Na = Nn = 0;
    for( jj=0; jj<NComp-1; jj++ ) // -1 = no check H2O
    {
        if( aZ[jj] > 0)
            Nc++;
        else if( aZ[jj] < 0 )
            Na++;
        else
            Nn++;
    }
    Ns = NComp-1;  // index of water-solvent
}



/// Calculate Etheta and Ethetap factors.
/// Reference: Anderson (2005), p. 610
void TPitzer::Ecalc( double z, double z1, double I, double DH_term,
        double& Etheta, double& Ethetap)
{
    long int k, m;
    double xMN, xMM, xNN,  x;
    double zet=0., dzdx=0.;
    double bk[23], dk[23];
    double JMN=0., JpMN=0., JMM=0., JpMM=0., JNN=0., JpNN=0.;

    // parameters for ak1 and ak2 values from Pitzer 1991
    static double ak1[21] = {  1.925154014814667, -0.060076477753119, -0.029779077456514,
                      -0.007299499690937,  0.000388260636404,  0.000636874599598,
                       0.000036583601823, -0.000045036975204, -0.000004537895710,
                       0.000002937706971,  0.000000396566462, -0.000000202099617,
                      -0.000000025267769,  0.000000013522610,  0.000000001229405,
                      -0.000000000821969, -0.000000000050847,  0.000000000046333,
                       0.000000000001943, -0.000000000002563, -0.000000000010991 };

    static double ak2[23] = {  0.628023320520852,  0.462762985338493,  0.150044637187895,
                      -0.028796057604906, -0.036552745910311, -0.001668087945272,
                       0.006519840398744,  0.001130378079086, -0.000887171310131,
                      -0.000242107641309,  0.000087294451594,  0.000034682122751,
                      -0.000004583768938, -0.000003548684306, -0.000000250453880,
                       0.000000216991779,  0.000000080779570,  0.000000004558555,
                      -0.000000006944757, -0.000000002849257,  0.000000000237816,
                       0.0, 0.0 };

    xMN= 6. * z*z1 * DH_term * pow(I,0.5);
    xMM= 6. * z1*z1 * DH_term * pow(I,0.5);
    xNN= 6. * z*z * DH_term * pow(I,0.5);

    for( k=1; k<=3; k++ )
    {
        if( k==1)
            x=xMN;
        else if( k==2 )
            x=xMM;
        else
          x=xNN;

        if( x <= 1 )
        {
            int sign = (x<0)?(-1):(1);
            zet = sign * 4.0 * pow(fabs(x), 0.2) - 2.0;
            dzdx = sign * 0.8 * 1./pow(fabs(x),(0.8));
                // zet=4.0 * pow(x, 0.2) - 2.0;
                // dzdx=0.8 * pow(x,(-0.8));
            bk[22] = 0.; bk[21] = 0.;
            dk[21] = 0.; dk[22] = 0.;
            for( m=20; m>=0; m--)
            {
                bk[m]= zet * bk[m+1] - bk[m+2] + ak1[m];
                dk[m]= bk[m+1] + zet * dk[m+1]- dk[m+2];
            }
        }

        else if( x > 1)
        {
            zet=-22.0/9.0 + (40.0/9.0) * 1./pow(x,(0.1));
            dzdx= (-40.0/90.) * 1./pow(x,(11./10.));
            bk[22]=0.; bk[21]=0.;
            dk[21]=0.; dk[22]=0.;
            for( m=20; m>=0; m--)
            {
                bk[m] = zet *bk[m+1] - bk[m+2] + ak2[m];
                dk[m]=  bk[m+1] + zet *dk[m+1] - dk[m+2];
            }
        }

        if( k == 1 )
        {
            JMN=0.25*x -1. + 0.5* (bk[0]-bk[2]);
            JpMN=0.25 + 0.5*dzdx*(dk[0]-dk[2]);
        }

        else if( k==2 )
        {
            JMM=0.25*x -1. + 0.5*(bk[0]-bk[2]);
            JpMM=0.25 + 0.5*dzdx*(dk[0]-dk[2]);
        }

        else
        {
            JNN=0.25*x -1. +0.5*(bk[0]-bk[2]);
            JpNN=0.25 +0.5*dzdx*(dk[0]-dk[2]);
        }
    } //k

    Etheta=((z*z1) /(4.0*I)) * (JMN - 0.5*JMM - 0.5*JNN);
    Ethetap= - (Etheta/I) +((z*z1)/(8.0*I*I)) *(xMN*JpMN - 0.5*xMM*JpMM - 0.5*xNN*JpNN);
}


/// Calculate Z-Term, Pitzer-Toughreact Report 2006, equation (A8)
double TPitzer::Z_Term()
{
    double Zan=0., Zca=0., Z;
    long int a, c;

    for( a=0; a<Na; a++)
        Zan += za[a]*ma[a];

    for( c=0; c<Nc; c++)
        Zca +=zc[c]*mc[c];

    Z = fabs(Zan) + Zca;

    return Z;
}


/// Calculate Ionic Strength
double TPitzer::IonicStr( double& I )
{
    double Ia=0., Ic=0., IS;
    long int a, c;

    for( a=0; a<Na; a++ )
        Ia += za[a]*za[a]* ma[a];

    for( c=0; c<Nc; c++ )
        Ic += zc[c]* zc[c]* mc[c];

    IS =0.5*(Ia+Ic);
    I=IS;

    return pow(IS,0.5);
}


/// Calculate osmotic coefficient, activity, and activity coefficient of water-solvent
double TPitzer::lnGammaH2O( double DH_term )
{
    double OC1, OC2, alp, alp1, C, h1, h2, B3, OC3, OC3a, z, z1, Phiphi,
            OC4, OC4a, Phiphi1, OC5, OC5a, OC5b, OC6, OCges, OCmol, OC, Lna;
    double Etheta=0., Ethetap=0.;
    long int a, c, n, c1, a1;

    // Term OC1, Pitzer-Toughreact Report 2006, equation (A2)
    OC1 = (-(DH_term*pow(I,1.5)) / (1.+1.2*Is) );

    // Term OC2
    OC2 = alp = alp1 = C = h1 = h2 = B3 = 0.;
    for( c=0; c<Nc; c++)
    {
        for( a=0; a<Na; a++)
        {
            getAlp(  c,  a, alp, alp1 );
            C = Cphi[c][a] / (2.*sqrt(fabs(za[a]*zc[c])));	// Pitzer-Toughreact Report 2006, equation (A7)
            h1=alp*Is;
            h2=alp1*Is;
            if (alp1==0)
                B3 = Bet0[c][a]+ Bet1[c][a]*exp(-h1);
            else
                B3 = Bet0[c][a]+ Bet1[c][a]*exp(-h1)+(Bet2[c][a]*exp(-h2)); //Pitzer-Toughreact Report 2006 equation (A9)
            OC2 +=(mc[c]*ma[a]*(B3+Zfac*(C) ));
        }
    }

    // Term OC3
    OC3 = 0.; OC3a = 0.;
    for( c=0; c<(Nc-1); c++ )
    {
        for( c1=c+1; c1<Nc; c1++ )
        {
            for( a=0; a<Na; a++)
            {
                OC3a += (ma[a]*Psi[c][c1][a]);
            }
            z=zc[c];
            z1=zc[c1];
            Ecalc( z, z1, I, Aphi, Etheta,Ethetap);
            Theta[c1][c]=Theta[c][c1];
            Phiphi = Theta[c][c1] + Etheta + Ethetap * I;// * sqrt(I);	 Pitzer-Toughreact Report 2006, equation (A14)
            OC3 += (mc[c]*mc[c1]*(Phiphi + OC3a));
        }
    }

    // Term OC4
    OC4 = OC4a =0.;
    for( a=0; a<(Na-1); a++)
    {
        for( a1=a+1; a1<Na; a1++)
        {
            for( c=0; c<Nc; c++)
            {
                OC4a += (mc[c]*Psi1[a][a1][c]);
            }
            z=za[a];
            z1=za[a1];
            Ecalc(z,z1,I,Aphi, Etheta,Ethetap);
            Theta1[a1][a]=Theta1[a][a1];
            Phiphi1 = Theta1[a][a1] + Etheta + Ethetap * I;	// Pitzer-Toughreact Report, 2006 equation (A14)
            OC4 += (ma[a]*ma[a1]*(Phiphi1 + OC4a));
        }
    }

    // Term OC5
    OC5 = OC5a = OC5b = 0.;
    for(  n=0; n<Nn; n++)
    {
        for( c=0; c<Nc; c++)
        {
            OC5a +=(mn[n]*mc[c]*Lam[n][c]);
        }
    }

    for(  n=0; n<Nn; n++)
    {
        for( a=0; a<Na; a++)
        {
            OC5b +=(mn[n]*ma[a]*Lam1[n][a]);
        }
    }
    OC5=OC5a+OC5b;

    // Term OC6
    OC6 = 0.;
    for(  n=0; n<Nn; n++)
    {
        for( c=0; c<Nc; c++)
        {
            for( a=0; a<Na; a++)
            {
                OC6 +=(mn[n]*mc[c]*ma[a]*Zeta[n][c][a]);
            }
        }
    }

    OCges=2.*(OC1+OC2+OC3+OC4+OC5+OC6);
    OCmol= p_sum(aM, xcx, Nc)+ p_sum(aM, xax, Na)+ p_sum(aM, xnx, Nn);
    OC = 1.+ (OCges / OCmol);

    // Activity of Water
    Lna =(-18.1/1000.)*OC*OCmol;
        // double activityH2O = exp(Lna);

    return Lna-log(x[Ns]);
}


/// retrieve Alpha parameter
void TPitzer::getAlp( long int c, long int a, double& alp, double& alp1 )
{
    if( zc[c] == 1. || za[a] == -1. )
    {
        alp=2.;
            // alp1=12.;
        alp1=0.;
    }
    else if( zc[c] == 2. && za[a] == -2. )
    {
        alp=1.4;
        alp1=12.;
    }
    else
    {
        alp=2.0;
        alp1=50.;
     }
}


/// calculate g
double TPitzer::get_g( double x_alp )
{
    double g;

    g = 2.*(1.-(1.+x_alp)*exp(-x_alp))/(x_alp*x_alp);

    if(x_alp==0)
    {
        g = 0.;
    }
    return g;
}


/// calculate gp
double TPitzer::get_gp( double x_alp )
{
    double gp;

    gp = -2.*(1.-(1.+x_alp+x_alp*x_alp*0.5)*exp(-x_alp))/(x_alp*x_alp);
    if(x_alp==0){
        gp = 0.;
    }
    return gp;
}


/// Calculate F-Factor, Pitzer-Toughreact Report 2006, equation (A6)
double TPitzer::F_Factor( double DH_term )
{
    long int c, c1, a, a1;
    double z=0., z1=0., Etheta=0., Ethetap=0.;
    double F1, F2, Phip, F3, Phip1, F4, F;
    double alp, alp1, g1, g2, B1, x_alp;

    // Term F1
    F1=-DH_term*( (Is/(1.+1.2*Is)) + 2.*log(1.+1.2*Is)/1.2);

    // Term F2
    F2 = Phip=0.;
    for( c=0; c<(Nc-1); c++ )
    {
        for( c1=c+1; c1<Nc; c1++ )
        {
            z=zc[c];
            z1=zc[c1];
            Ecalc(z,z1,I,DH_term, Etheta,Ethetap);
            Phip = Ethetap;					//Pitzer-Toughreact Report 2006, equation (A16)
            F2 +=(mc[c]*mc[c1]*(Phip));
        }
    }


    // Term F3
    F3 = Phip1 = 0.;
    for( a=0; a<(Na-1); a++)
    {
        for( a1=a+1; a1<Na; a1++)
        {
            z=za[a];
            z1=za[a1];
            Ecalc(z,z1,I,DH_term, Etheta,Ethetap);
            Phip1=Ethetap;      				//Pitzer-Toughreact Report 2006, equation (A16)
            F3 +=(ma[a]*ma[a1]*(Phip1));
        }
    }


    // Term F4
    F4 = alp = alp1 = B1 = 0.;
    for( c=0; c<Nc; c++)
    {
        for( a=0; a<Na; a++)
        {
            getAlp(  c,  a, alp, alp1 );
            x_alp = alp*Is;
            g1 = get_gp( x_alp );
            x_alp = alp1*Is;
            g2 = get_gp( x_alp );
            B1= (Bet1[c][a]*g1)/I+ (Bet2[c][a]*g2)/I;
            F4 = F4+ (mc[c]*ma[a]*B1);
        }
    }
    F = F1+F2+F3+F4;
    return F;
}


/// Calculate lnGammaM - activity coefficient of a cation with index X
double TPitzer::lnGammaM( long int M, double DH_term  )
{
    double Etheta=0., Ethetap=0.;
    long int a, n, c1, a1;
    double GM1, GM2, alp, alp1, g1, g2, B2, C, x_alp,
                GM3, GM3a, Phi, z, z1, Q, GM4, GM5a, GM5, GM6a, GM6, GM;
    double actcoeffM;

    // Calculate GM1
    GM1 = (zc[M]*zc[M])*Ffac;

    // Term GM2
    GM2 = alp = alp1 = 0.;
    for( a=0; a<Na; a++)
    {
        getAlp(  M,  a, alp, alp1 );
        C = Cphi[M][a]/(2.*sqrt(fabs(za[a]*zc[M])));	// Pitzer-Toughreact Report 2006, equation (A7)
        x_alp = alp*Is;
        g1 = get_g( x_alp );
        x_alp = alp1*Is;
        g2 = get_g( x_alp );
        B2= Bet0[M][a]+(Bet1[M][a]*g1)+ (Bet2[M][a]*g2); // Pitzer-Toughreact Report 2006, equation (A10)
        GM2 = GM2+(ma[a]*(2.*B2+Zfac*(C)));
    }


    // Term GM3
    GM3 = GM3a = 0.;
    for( c1=0; c1<Nc; c1++)
    {
        for( a=0; a<Na; a++)
        {
            Psi[c1][M][a] = Psi[M][c1][a];
            GM3a += ma[a]*Psi[M][c1][a];
        }
        z = zc[M];
        z1 = zc[c1];
        Ecalc(z,z1,I,DH_term ,Etheta,Ethetap);
        Theta[c1][M] = Theta[M][c1];
        Phi = Theta[M][c1]+Etheta;  					// Pitzer-Toughreact Report 2006, equation (A15)

        if(c1==M)
        {
            Q=0;
        }
        else
        {
            Q=1;
        }
        GM3=GM3+Q*mc[c1]*(2.*Phi+ GM3a);
     }

    // Term GM4
    GM4 = 0.;
    for( a=0; a<(Na-1); a++)
    {
        for( a1=a+1; a1<Na; a1++)
        {
            Psi1[a1][a][M]=Psi1[a][a1][M];
            GM4=GM4+(ma[a]*ma[a1]*Psi1[a][a1][M]);
        }
    }



    // Term GM5
    GM5a = 0.;
    for( c1=0; c1<Nc; c1++)
    {
        for( a=0; a<Na; a++)
        {
            C = Cphi[c1][a]/(2.*sqrt(fabs(za[a]*zc[c1])));			// Pitzer-Toughreact Report 2006, equation (A7)
            GM5a = GM5a+(mc[c1]*ma[a]* (C) );
        }
    }
    GM5 = zc[M]*GM5a;


    // Term GM6
    GM6a = 0;
    for( n=0; n<Nn; n++)
    {
        GM6a += mn[n]*Lam[n][M];
    }
    GM6 = 2*GM6a;

    // Term GM
    GM = GM1+GM2+GM3+GM4+GM5+GM6;

    actcoeffM = exp(GM);
    return GM;
}


/// Calculate lnGammaX - activity coefficient of an anion with index X
double TPitzer::lnGammaX( long int X, double DH_term )
{
    double Etheta=0., Ethetap=0.;
    long int c, n, c1, a1;
    double GX1;
    double GX2, C, g1, g2, B2, alp, alp1, x_alp;
    double GX3, GX3a, z, z1, Phi1, Q ;
    double GX4;
    double GX5a, GX5;
    double GX6a, GX6;
    double GX, actcoeffX;

    // Term GX1 (Pitzer-Toughreact Report 2006, equation A4)
    GX1=(za[X]*za[X])*Ffac;

    // Term GX2
    GX2 = C = 0.;
    for( c=0; c<Nc; c++)
    {
        getAlp(  c,  X, alp, alp1 );
        C = Cphi[c][X]/(2.*sqrt(fabs(za[X]*zc[c])));
        x_alp = alp*Is;
        g1 = get_g( x_alp );
        x_alp = alp1*Is;
        g2 = get_g( x_alp );
        B2= Bet0[c][X]+Bet1[c][X]*g1+ (Bet2[c][X]*g2);
        GX2=GX2+(mc[c]*(2.*B2+Zfac*(C)));
    }

    // Term GX3
    GX3 = 0.; GX3a = 0.;
    for( a1=0; a1<Na; a1++)
    {
        for( c=0; c<Nc; c++)
        {
            Psi1[a1][X][c]=Psi1[X][a1][c];
            GX3a += mc[c]*Psi1[X][a1][c];
        }
        z = za[X];
        z1 = za[a1];
        Ecalc(z,z1,I,DH_term , Etheta,Ethetap);
        Theta1[a1][X] = Theta1[X][a1];
        Phi1 = Theta1[X][a1]+Etheta;

        if(a1 == X)
        {
            Q = 0;
        }
        else
        {
            Q = 1;
        }
        GX3 = GX3+Q*ma[a1]*(2.*Phi1+GX3a);
    }


    // Term GX4
    GX4 = 0.;
    for( c=0; c<(Nc-1); c++)
    {
        for( c1 = c+1; c1<Nc; c1++)
        {
            Psi[c1][c][X] = Psi[c][c1][X];
            GX4 = GX4+(mc[c]*mc[c1]*Psi[c][c1][X]);
        }
    }

    // Term GX5
    GX5a = 0.;
    for( c=0; c<Nc; c++)
    {
        for( a1=0; a1<Na; a1++)
        {
            C = Cphi[c][a1]/(2.*sqrt(fabs(za[a1]*zc[c])));	 // Pitzer-Toughreact Report 2006, equation (A7)
            GX5a = GX5a+(mc[c]*ma[a1]* C);
        }
    }
    GX5 = fabs(za[X])*GX5a;

    // Term GX6
    GX6a = 0.;
    for( n=0; n<Nn; n++)
        GX6a += mn[n]*Lam1[n][X];
    GX6 = 2.*GX6a;

    // Term GX
    GX = GX1+GX2+GX3+GX4+GX5+GX6;

    actcoeffX = exp(GX);

    return GX;
}


/// Calculate lngammaN - activity coefficient of a neutral species with index N
double TPitzer::lnGammaN( long int N )
{
    long int c, a;
    double GN1, GN2, GN3, GN, actcoeffN;

    // Term GN1
    GN1 = 0.;
    for( a=0; a<Na; a++)
        GN1 = GN1+(ma[a]*2.*Lam1[N][a]);

    // Term GN2
    GN2 = 0.;
    for( c=0; c<Nc; c++)
        GN2 = GN2+(mc[c]*2.*Lam[N][c]);

    // Term GN3
    GN3 = 0.;
    for( c=0; c<Nc; c++)
    {
        for( a=0; a<Na; a++)
        {
            GN3 = GN3+(mc[c]*ma[a]*Zeta[N][c][a]);
        }
    }

    // Term GN
    GN = GN1+GN2+GN3;

    actcoeffN = exp(GN);

  return GN;
}


/// mean activity coefficient of KCl in binary system KCl-H2O at system ionic Strength
///	 	and temperature
double TPitzer::McInnes_KCl( )
{
    double gammaK, gammaKCl, ln_gammaK; // ln_gammaCl, gammaCl;
    double mK, mCl, mH, mOH;
    double term1, term2, term3, term4, term5, F_McInnes, term6; // term7, term8, term9, term10;
    double v, Z;

    double B0_KCl = McI_PT_array[0];
    double B1_KCl = McI_PT_array[1];
    double Cphi_KCl = McI_PT_array[2];
    double B0_KOH = McI_PT_array[3];
    double B1_KOH = McI_PT_array[4];
    double Cphi_KOH = McI_PT_array[5];
    double B0_HCl = McI_PT_array[6];
    double B1_HCl = McI_PT_array[7];
    double Cphi_HCl = McI_PT_array[8];
    double Theta_KH = McI_PT_array[9];
    double Psi_KHCl = McI_PT_array[10];
//	double Theta_ClOH = McI_PT_array[11];
    double Psi_ClOHK = McI_PT_array[12];

    double C_KCl, C_KOH, C_HCl;
    double alp = 2.0;
    double B_KCl, B_KOH, B_HCl;
    double g;
    long int M, N, X;

    // Computation of C from Cphi parameters
    C_KCl = Cphi_KCl /2.;
    C_KOH = Cphi_KOH /2.;
    C_HCl = Cphi_HCl /2.;

    // Internal alpha parameter
    v = alp*sqrt(I);

    // Internal B parameters
    g = 2.*(1.-(1.+v)*exp(-v))/(v*v);
    B_KCl = B0_KCl + B1_KCl*g;
    B_KOH = B0_KOH + B1_KOH*g;
    B_HCl = B0_HCl + B1_HCl*g;

    // Calculation of molalities from Ionic Strength of Solution (pH fixed at 7)
    mH = 1.0e-7;
    mOH = 1.0e-7;
    mK = I - mH;
    mCl = I - mH;
    Z = mH+mOH+mK+mCl;

    // F_Term
    term1 = -Aphi*( (Is/(1.+1.2*Is)) + 2.*log(1.+1.2*Is)/1.2);
    term2 =	mK*mCl*B1_KCl*get_gp(v)/I + mK*mOH*B1_KOH*get_gp(v)/I + mH*mCl*B1_HCl*get_gp(v)/I ;
    F_McInnes = term1 + term2;

    // activity coeficient of K+
    term3 = mCl*( 2.* B_KCl + Z * C_KCl ) +
    mOH*( 2.* B_KOH + Z * C_KOH );

    term4 = mH*(2.*Theta_KH + mCl*Psi_KHCl);  // + mOH*Psi_KHOH
    term5 = mCl * mOH * Psi_ClOHK;
    term6 = mK*mCl*C_KCl + mK*mOH*C_KOH + mH*mCl*C_HCl;
    ln_gammaK = F_McInnes + term3 + term4 + term5 + term6;
    gammaK	= exp(ln_gammaK);
    gammaKCl = gammaK;


    //Calculation of unscaled Pitzer Activity Coefficients
    Pitzer_calc_Gamma( );

    // Scaling of unscaled Pitzer activity coefficients according to McInnes convention

    // this should be changed to automatic retrieval of chlorine location instead of manual input
    long int Q;
    Q = 0;
    GammaMcI[Ns] = exp( lnGamma[Ns] );

    for( M=0; M<Nc; M++ )
    {
        GammaMcI[xcx[M]] = (exp( lnGamma[xcx[M]] )) * (pow ( ( (exp( lnGamma[xax[Q]] ))
                / gammaKCl ), zc[ M ]) );
    }

    for( X=0; X<Na; X++ )
    {
        GammaMcI[xax[X]] = (exp( lnGamma[xax[X]] )) * (pow ( ( (exp( lnGamma[xax[Q]] ))
                / gammaKCl ), za[ X ]) );
    }

    if( Nn > 0 )
        for( N=0; N<Nn; N++ )
            GammaMcI[xnx[N]] = exp( lnGamma[xnx[N]] );

    return 0;
}


/// Calculation of bulk Excess Gibbs Energy per kilogram of water
long int TPitzer::ExcessProp( double *Zex )
{
    long int M, N, X;
    double OsmCoeffGex, catGex, aniGex, neutGex;
    double OsmCoeffSex, catSex, aniSex, neutSex;
    double OsmCoeffCPex, catCPex, aniCPex, neutCPex;
    OsmCoeffGex=0.0, catGex=0.0, aniGex=0.0, neutGex=0.0;

    Is = IonicStr( I );
    Ffac = F_Factor( Aphi );
    Zfac = Z_Term();
    OsmCoeffGex = lnGammaH2O( Aphi );

    for( M=0; M<Nc; M++ )
    {
        catGex += ( mc[M]* ( 1. - OsmCoeffGex + lnGammaM( M , Aphi  ) ) );
    }

    for( X=0; X<Na; X++ )
    {
        aniGex += ( ma[X]* ( 1. - OsmCoeffGex + lnGammaX( X , Aphi  ) ) );
    }

    if( Nn > 0 )
    {
        for( N=0; N<Nn; N++ )
        {
            neutGex += ( mn[N]* ( 1. - OsmCoeffGex + lnGammaN( N ) ) );
        }
    }

    Gex = R_CONST * Tk * ( catGex + aniGex + neutGex );


    // Calculation of bulk Excess Entropy per kilogram of water
    OsmCoeffSex=0.0, catSex=0.0, aniSex=0.0, neutSex=0.0;
    PTcalc( 1 );
    Ffac = F_Factor( dAphidT );
    OsmCoeffSex = lnGammaH2O( dAphidT );

    for( M=0; M<Nc; M++ )
    {
        catSex += ( mc[M]* ( lnGammaM( M , dAphidT ) - OsmCoeffSex  ) );
    }

    for( X=0; X<Na; X++ )
    {
        aniSex += ( ma[X]* ( lnGammaX( X , dAphidT ) - OsmCoeffSex ) );
    }

    if( Nn > 0 )
    {
        for( N=0; N<Nn; N++ )
        {
            neutSex += ( mc[M]* ( lnGammaN( N ) - OsmCoeffSex ) );
        }
    }

    // Excess entropy
    Sex = -( Gex / Tk + R_CONST * Tk * ( catSex + aniSex + neutSex ) );

    // Excess Enthalpy per kilogram of water
    Hex = Gex + Tk * Sex;

    // Calculation of bulk Excess Heat Capacity per kilogram of water
    OsmCoeffCPex=0.0, catCPex=0.0, aniCPex=0.0, neutCPex=0.0;

    // correct interaction parameters and Debye H�ckel term to T (and P) of system
    PTcalc( 2 );
    Ffac = F_Factor( d2AphidT2 );

    // activity and osmotic coefficients
    OsmCoeffCPex = lnGammaH2O( d2AphidT2 );

    for( M=0; M<Nc; M++ )
    {
        catCPex += ( mc[M]* ( lnGammaM( M , d2AphidT2 ) - OsmCoeffCPex  ) );
    }

    for( X=0; X<Na; X++ )
    {
        aniCPex += ( ma[X]* ( lnGammaX( X , d2AphidT2 ) - OsmCoeffCPex ) );
    }

    if( Nn > 0 )
    {
        for( N=0; N<Nn; N++ )
        {
            neutCPex += ( mc[M]* ( lnGammaN( N ) - OsmCoeffCPex ) );
        }
    }

    // Excess heat capacity
        CPex = (-4.)*Sex -2.*Gex/Tk + R_CONST*Tk*Tk*(catCPex + aniCPex + neutCPex);

    // Assignments
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





//=============================================================================================
// Extended universal quasi-chemical (EUNIQUAC) model for aqueous electrolyte solutions
// References: Nicolaisen et al. (1993), Thomsen et al. (1996), Thomsen (2005)
// (c) TW/FH February 2009
//=============================================================================================


// Generic constructor for the TEUNIQUAC class
TEUNIQUAC::TEUNIQUAC( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW ):
                TSolMod( sd )
{
    alloc_internal();
    z = arZ;
    m = arM;
    RhoW = dW;
    EpsW = eW;
}


TEUNIQUAC::~TEUNIQUAC()
{
    free_internal();
}


void TEUNIQUAC::alloc_internal()
{
    R = new double [NComp];
    Q = new double [NComp];
    Phi = new double [NComp];
    Theta = new double [NComp];
    U = new double *[NComp];
    dU = new double *[NComp];
    d2U = new double *[NComp];
    Psi = new double *[NComp];
    dPsi = new double *[NComp];
    d2Psi = new double *[NComp];

    for (long int j=0; j<NComp; j++)
    {
        U[j] = new double [NComp];
        dU[j] = new double [NComp];
        d2U[j] = new double [NComp];
        Psi[j] = new double [NComp];
        dPsi[j] = new double [NComp];
        d2Psi[j] = new double [NComp];
    }
}


void TEUNIQUAC::free_internal()
{
    for (long int j=0; j<NComp; j++)
    {
        delete[]U[j];
        delete[]dU[j];
        delete[]d2U[j];
        delete[]Psi[j];
        delete[]dPsi[j];
        delete[]d2Psi[j];
    }
    delete[]R;
    delete[]Q;
    delete[]Phi;
    delete[]Theta;
    delete[]U;
    delete[]dU;
    delete[]d2U;
    delete[]Psi;
    delete[]dPsi;
    delete[]d2Psi;
}


///  Calculates T,P corrected binary interaction parameters
long int TEUNIQUAC::PTparam()
{
    long int j, i, ip, i1, i2;
    double u0, u1, u, du, d2u, psi, dpsi, d2psi, v, dv;
    double c, alp, bet, dal, rho, eps, dedt, d2edt2, dedp;

    if ( NPcoef < 2 || NPar < 1 || NP_DC < 2 )
       return 1;

    // read and transfer species-dependent parameters
    for (j=0; j<NComp; j++)
    {
        R[j] = aDCc[NP_DC*j];   // volume parameter r
        Q[j] = aDCc[NP_DC*j+1];   // surface parameter q
    }

    // fill internal arrays of interaction parameters with standard value
    for (j=0; j<NComp; j++)
    {
        for (i=0; i<NComp; i++)
        {
            U[j][i] = 0.0;
            dU[j][i] = 0.0;
            d2U[j][i] = 0.0;
            Psi[j][i] = 1.0;
            dPsi[j][i] = 0.0;
            d2Psi[j][i] = 0.0;
        }
    }

    // read and convert interaction energies (uji) that have non-standard value
    for (ip=0; ip<NPar; ip++)
    {
        i1 = aIPx[MaxOrd*ip];
        i2 = aIPx[MaxOrd*ip+1];
        u0 = aIPc[NPcoef*ip+0];
        u1 = aIPc[NPcoef*ip+1];
        u = u0 + u1*(Tk-298.15);
        du = u1;
        d2u = 0.0;
        U[i1][i2] = u;
        dU[i1][i2] = du;
        d2U[i1][i2] = d2u;
        U[i2][i1] = u;
        dU[i2][i1] = du;
        d2U[i2][i1] = d2u;   // uij identical to uji
    }

    // calculate Psi and its partial derivatives
    for (j=0; j<NComp; j++)
    {
        for (i=0; i<NComp; i++)
        {
            psi = exp( -(U[j][i]-U[i][i])/Tk );
            v = (U[j][i]-U[i][i])/pow(Tk,2.) - (dU[j][i]-dU[i][i])/Tk;
            dv = (-2.)*(U[j][i]-U[i][i])/pow(Tk,3.) + 2.*(dU[j][i]-dU[i][i])/pow(Tk,2.)
                    - (d2U[j][i]-d2U[i][i])/Tk;
            dpsi = psi * v;
            d2psi = dpsi*v + psi*dv;
            Psi[j][i] = psi;
            dPsi[j][i] = dpsi;
            d2Psi[j][i] = d2psi;
        }
    }

    // read and convert rho and eps (check density units)
    c = (1.3287e+5);
    rho = RhoW[0]*1000.;
    alp = - 1./rho*RhoW[1]*1000.;
    dal = pow(alp,2.) - 1./rho*RhoW[2];
    bet = 1./rho*RhoW[3]*1000.;
    eps = EpsW[0];
    dedt = 1./eps*EpsW[1];
    d2edt2 = - 1./pow(eps,2.)*pow(EpsW[1],2.) + 1./eps*EpsW[2];
    dedp = 1./eps*EpsW[3];

    // calculate A term of Debye-Huckel equation (and derivatives)
    A = c*sqrt(rho)/pow((eps*Tk),1.5);
    dAdT = - 3./2.*A*( dedt + 1./Tk + alp/3. );
    d2AdT2 = 1./A*pow(dAdT,2.) - 3./2.*A*( d2edt2 - 1/pow(Tk,2.) + 1/3.*dal );
    dAdP = 1./2.*A*( bet - 3.*dedp);

    // approximation valid only for temperatures below 200 deg. C and Psat
    A = 1.131 + (1.335e-3)*(Tk-273.15) + (1.164e-5)*pow( (Tk-273.15), 2.);
    dAdT = (1.335e-3) + 2.*(1.164e-5)*(Tk-273.15);
    d2AdT2 = 2.*(1.164e-5);
    dAdP = 0.;

    ErrorIf( fabs(A) < 1e-9, "EUNIQUAC model",
            "Error: DH parameter A was not calculated - no values of RhoW and EpsW !" );

    return 0;
}


/// Calculates activity coefficients
long int TEUNIQUAC::MixMod()
{
    long int j, i, l, k, w;
    double Mw, Xw, b, RR, QQ, K, L, M, gamDH, gamC, gamR, lnGam, Gam;
    b = 1.5; Mw = 0.01801528;

    // get index of water (assumes water is last species in phase)
    w = NComp - 1;
    Xw = x[w];

    // calculation of ionic strength
    IonicStrength();

    // calculation of Phi and Theta terms
    for (j=0; j<NComp; j++)
    {
        RR = 0.0;
        QQ = 0.0;
        for (i=0; i<NComp; i++)
        {
            RR += x[i]*R[i];
            QQ += x[i]*Q[i];
        }
        Phi[j] = x[j]*R[j]/RR;
        Theta[j] = x[j]*Q[j]/QQ;
    }

    // loop over species
    for (j=0; j<NComp; j++)
    {
        // species other than water solvent
        if (j < w)
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

            gamDH = - pow(z[j],2.)*A*sqrt(IS)/(1.+b*sqrt(IS));
            gamC = log(Phi[j]/x[j]) - Phi[j]/x[j] - log(R[j]/R[w]) + R[j]/R[w]
                    - 5.0*Q[j] * ( log(Phi[j]/Theta[j]) - Phi[j]/Theta[j]
                    - log(R[j]*Q[w]/(R[w]*Q[j])) + R[j]*Q[w]/(R[w]*Q[j]) );
            gamR = Q[j] * ( - log(K) - L + log(Psi[w][j]) + Psi[j][w] );
            lnGam = gamDH + gamC + gamR;

            // convert activity coefficient to molality scale
            lnGam = lnGam + log(x[w]);
            lnGamma[j] = lnGam;

            // write debug results
            Gam = exp(lnGam);
            gammaDH[j] = gamDH;
            gammaC[j] = gamC;
            gammaR[j] = gamR;

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

            gamDH = Mw*2.*A/pow(b,3.) * ( 1. + b*sqrt(IS) - 1./(1.+b*sqrt(IS)) - 2*log(1.+b*sqrt(IS)) );
            gamC = log(Phi[j]/x[j]) + 1. - Phi[j]/x[j] - 5.0*Q[j] * ( log(Phi[j]/Theta[j]) + 1. - Phi[j]/Theta[j] );
            gamR = Q[j] * (1. - log(K) - L );
            lnGam = gamDH + gamC + gamR;
            lnGamma[j] = lnGam;

            // write debug results
            Gam = exp(lnGam);
            gammaDH[j]=gamDH;
            gammaC[j] = gamC;
            gammaR[j] = gamR;
        }
    }

    return 0;
}


long int TEUNIQUAC::ExcessProp( double *Zex )
{
    long int j, i, w;
    double Mw, Xw, b, phiti, phthi, RR, QQ, N, TPI, tpx, TPX,
                dtpx, DTPX, con;
    double gDH, gC, gR, hR, cpR, gCI, gRI, gCX, gRX, dg, d2g, dgRI, d2gRI,
                dgRX, d2gRX, dgDH, d2gDH, dgDHdP;
    b = 1.5; Mw = 0.01801528;

    // get index of water (assumes water is last species in phase)
    w = NComp - 1;
    Xw = x[w];

    // calculation of ionic strength
    IonicStrength();

    // calculation of Phi and Theta terms
    for (j=0; j<NComp; j++)
    {
        RR = 0.0;
        QQ = 0.0;
        for (i=0; i<NComp; i++)
        {
            RR += x[i]*R[i];
            QQ += x[i]*Q[i];
        }
        Phi[j] = x[j]*R[j]/RR;
        Theta[j] = x[j]*Q[j]/QQ;
    }

    // calculation of bulk phase excess properties
    Gex = 0.0; Hex = 0.0; Sex = 0.0; CPex = 0.0; Vex = 0.0;
    gC = 0.0; gR = 0.0; hR = 0.0; cpR = 0.0;

    // infinite dilution part
    gCI = 0.; gRI = 0.; dgRI = 0.; d2gRI = 0.;

    for (j=0; j<NComp; j++) // loop over species
    {
        if (j == w)
        {
            gCI += 0.;
            gRI += 0.;
            dgRI += 0.;
            d2gRI += 0.;
        }

        else
        {
            phiti = R[j]/R[w];
            phthi = R[j]*Q[w]/(Q[j]*R[w]);
            gCI += x[j]*( 1.0 + log(phiti) - phiti - 5.*Q[j]*(1.+log(phthi)-phthi) );
            gRI += x[j]*Q[j]*(1.-log(Psi[w][j])-Psi[j][w]);
            dgRI += - x[j]*Q[j]*( pow(Psi[w][j],(-1.))*dPsi[w][j] + dPsi[j][w] ) ;
            d2gRI += x[j]*Q[j]*( pow(Psi[w][j],(-2.))*dPsi[w][j]*dPsi[w][j]
                       - pow(Psi[w][j],(-1.))*d2Psi[w][j] - d2Psi[j][w] );
        }
    }

    // combinatorial and residual part
    gCX = 0.; gRX = 0.; dgRX = 0.; d2gRX = 0.;

    for (j=0; j<NComp; j++)  // loop over species
    {
        N = 0.0;
        TPI = 0.0;
        tpx = 0.0;
        TPX = 0.0;
        dtpx = 0.0;
        DTPX = 0.0;

        for (i=0; i<NComp; i++)
        {
            N += Theta[i]*Psi[i][j];
            tpx += Theta[i]*dPsi[i][j];
            dtpx += Theta[i]*d2Psi[i][j];
        }

        TPI = 1./N;
        TPX = tpx*TPI;
        DTPX = - TPX*TPX + dtpx*TPI;
        gCX += x[j]*log(Phi[j]/x[j]) - 5.0*(Q[j]*x[j]*log(Phi[j]/Theta[j]));
        gRX += ( - x[j]*Q[j]*log(N) );
        dgRX += x[j]*Q[j]*TPX;
        d2gRX += x[j]*Q[j]*DTPX;
    }

    // DH part
    con = - x[w]*Mw*4./pow(b,3.) * ( log(1.+b*sqrt(IS)) - b*sqrt(IS) + 0.5*pow(b,2.)*IS );
    gDH = con*A;
    dgDH = - con*dAdT;
    d2gDH = - con*d2AdT2;
    dgDHdP = con*dAdP;

    // increment thermodynamic properties
    dg = ( dgDH + dgRX + dgRI );
    d2g = ( d2gDH + d2gRX + d2gRI );
    Gex = ( gDH + gRX + gCX - gRI - gCI ) * R_CONST * Tk;
    Hex = dg * pow(Tk,2.) * R_CONST;
    CPex = ( 2.*Tk*dg + pow(Tk,2.)*d2g ) * R_CONST;
    Sex = (Hex-Gex)/Tk;
    Vex = dgDHdP;
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
long int TEUNIQUAC::IdealProp( double *Zid )
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


/// Calculate ionic strength
long int TEUNIQUAC::IonicStrength()
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


/// Output of test results into text file (standalone variant only)
void TEUNIQUAC::Euniquac_test_out( const char *path )
{
    long int ii;//, c, a, n;

    // const ios::open_mode OFSMODE = ios::out � ios::app;
    ofstream ff(path, ios::app );
    ErrorIf( !ff.good() , path, "Fileopen error");

    ff << endl << "Vector of interaction parameters corrected to T,P of interest" << endl;
    for( ii=0; ii<NPar; ii++ )
        ff << aIP[ii] << "  ";

    ff << endl << "Debye-H�ckel contribution to Activity Coefficients" << endl;
    for( ii=0; ii<NComp; ii++ )
        ff << gammaDH[ii] << "  ";

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
}


//--------------------- End of s_solmod5.cpp ---------------------------
