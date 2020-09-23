//-------------------------------------------------------------------
// $Id: s_fgl3.cpp 724 2012-10-02 14:25:25Z kulik $
//
// Stub implementation of TModOther derived class of TSolMod
//
// Copyright (C) 2008-2012  T.Wagner, D.Kulik, S.Dmitrieva
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
// Customized solid-solution and fluid models
// References:
// (c) TW March 2009
//=============================================================================================


// Generic constructor for the TModOther class
TModOther::TModOther( SolutionData *sd, double *dW, double *eW ):
                TSolMod( sd )
{
    alloc_internal();
}


TModOther::~TModOther()
{
    free_internal();
}


void TModOther::alloc_internal()
{
    Gdqf = new double [NComp];
    Hdqf = new double [NComp];
    Sdqf = new double [NComp];
    CPdqf = new double [NComp];
    Vdqf = new double [NComp];
}


void TModOther::free_internal()
{
	delete[]Gdqf;
	delete[]Hdqf;
	delete[]Sdqf;
	delete[]CPdqf;
	delete[]Vdqf;
}


/// calculates pure species properties (pure fugacities, DQF corrections)
long int TModOther::PureSpecies()
{
    /*
	if( !strncmp( PhaseName, "FELDSPAR", 8 ) )
	{
		TModOther::PT_Feldspar();
	}

	else if{ !strncmp( PhaseName, "Garnet", 6 ) )
	{
		TModOther::PT_Garnet();
	}
    */
	return 0;
}


/// calculates T,P corrected binary interaction parameters
long int TModOther::PTparam()
{
    /*
	if( !strncmp( PhaseName, "FELDSPAR", 8 ) )
	{
		TModOther::PT_Feldspar();
	}

	else if{ !strncmp( PhaseName, "Garnet", 6 ) )
	{
		TModOther::PT_Garnet();
	}
    */
	return 0;
}


/// calculates activity coefficients
long int TModOther::MixMod()
{
    /*
	if( !strncmp( PhaseName, "FELDSPAR", 8 ) )
	{
		TModOther::PT_Feldspar();
	}

	else if{ !strncmp( PhaseName, "Garnet", 6 ) )
	{
		TModOther::PT_Garnet();
	}
    */
	return 0;
}


/// calculates excess properties
long int TModOther::ExcessProp( double *Zex )
{
    /*
	if( !strncmp( PhaseName, "FELDSPAR", 8 ) )
	{
		TModOther::PT_Feldspar();
	}

	else if{ !strncmp( PhaseName, "Garnet", 6 ) )
	{
		TModOther::PT_Garnet();
	}
    */

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
long int TModOther::IdealProp( double *Zid )
{
    /*
	if( !strncmp( PhaseName, "FELDSPAR", 8 ) )
	{
		TModOther::PT_Feldspar();
	}

	else if{ !strncmp( PhaseName, "Garnet", 6 ) )
	{
		TModOther::PT_Garnet();
	}
    */

	// assignments (ideal mixing properties)
	Gid = Hid - Sid*Tk;
	Aid = Gid - Vid*Pbar;
	Uid = Hid - Vid*Pbar;
	Zid[0] = Gid;
	Zid[1] = Hid;
	Zid[2] = Sid;
	Zid[3] = CPid;
	Zid[4] = Vid;
	Zid[5] = Aid;
	Zid[6] = Uid;

	return 0;
}



//--------------------- End of s_solmod6.cpp ---------------------------

