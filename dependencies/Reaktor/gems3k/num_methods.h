//-------------------------------------------------------------------
// $Id: num_methods.h 771 2012-12-13 13:07:43Z kulik $
//
/// \file num_methods.h
/// Declarations of C/C++ Numerical Methods used in GEMS3K code.
//
// Copyright (C) 2006-2012 S.Dmytriyeva, D.Kulik
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

#ifndef _num_methods_h_
#define _num_methods_h_

#include <cmath>

// Calculate number of points from iterators
long int  getNpoints( double Tai[4] );
double    getStep( double *Tai, int nPoints );


// Lagrangian interpolation functions
double LagranInterp(float *y, float *x, double *d, float yoi,
		float xoi, int M, int N, int pp );
double LagranInterp(float *y, float *x, float *d, float yoi,
		float xoi, int M, int N, int pp );
double LagranInterp(double *y, double *x, double *d, double yoi,
		double xoi, long int M, long int N, long int pp );


// generic functions for calculating partial derivatives
double quot( double u, double v, double du, double dv );

double quot( double u, double v, double du, double dv, double d2u, double d2v );


double prod2( double u, double v, double du, double dv );

double prod2 ( double u, double v, double du, double dv, double d2u, double d2v );


double prod3 ( double u, double v, double w, double du, double dv, double dw );

double prod3 ( double u, double v, double w, double du, double dv, double dw,
		double d2u, double d2v, double d2w );


typedef double (*minFunction)(double x, double y );

/// Data for minimization of convex one parameter function ( f(x)=>0 )
struct GoldenSectionData
{
    double Fa;
    double Fb;
    double Fx1;
    double Fx2;
    double a;
    double b;
    double x1;
    double x2;
    double Xtol;


    GoldenSectionData( double xstart, double xend, double Xtol_)
    {
      Xtol = Xtol_;
      x1 = xstart;
      x2 = xend;
      a = min( x1, x2 );
      b = max( x1, x2 );
   }

 };


/// Class for minimization of convex one parameter function ( f(x)=>0 )
/// Method of Gold Selection
class GoldenSection
{
protected:

  GoldenSectionData dat1;
  double Ftol;
  minFunction minF;

public:

  /// Golden Selection in interval x1 to x2, to minimize function f_proc
  /// xtol, ftol tolerance for the parameter and function
  GoldenSection( double x1, double x2, double xtol, double ftol,
                   double (f_proc)(double val, double val2 )):
   dat1(x1,x2,xtol),Ftol(ftol)
  {
    minF = f_proc;
  }

  virtual double calcFunction( double x, double y )
  {
      return minF( x, y );
  }

  virtual double getMinimum( double val2=0 )
  {
     return getMinimumDat( dat1, val2 );
  }

  virtual double getMinimumDat( GoldenSectionData dat, double val2 );

};

/// Class for minimization of convex two parameter function ( f(x,y)=>0 )
/// Method of Gold Selection
class GoldenSectionTwo : public GoldenSection
{
  GoldenSectionData dat2;
  int nOperand;   // minimization for first or second parameter
  double minX;
  double minY;

public:

  /// Golden Selection in intervals x1 to x2, y1 to y2 to minimize function f_proc
  /// xtol, ytol, ftol tolerance for the parameters and function
  GoldenSectionTwo( double x1, double x2, double xtol,
                      double y1, double y2, double ytol,
                      double ftol, double (f_proc)(double val, double val2 )):
  GoldenSection( x1, x2, xtol, ftol, f_proc),  dat2(y1,y2,ytol)
  {}

  double calcFunction( double x, double y );
  double getMinimum( double val2=0 );

  double getMinX()
  { return minX; }
  double getMinY()
  { return minY; }

};


// Golden Selection in interval x1 to x2, to minimize function f_proc
//double GoldenSection( double x1, double x2, double xtol, double ftol,
//                        double val2, double (f_proc)(double val, double val2 ));
// Method of Gold Selection for two parameters
//   ystart, yend, ytol,  parameter y    - start, end, tolerance
//   xstart, xend, xtol,  parameter x    - start, end, tolerance
//   ftol  function tolerance
//   ymin and xmin return values
//   f_proc function to minimize ( f( y, x)=>0 )
//void GoldenSection2( double ystart, double yend, double Ytol,
//                       double xstart, double xend, double Xtol,
//                       double Ftol, double& ymin, double& xmin,
//                       double (f_proc)(double valy, double valx ));



#endif   // _num_methods_h_

//-----------------------End of num_methods.h--------------------------

