// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "WaterHelmholtzStateHGK.hpp"

// C++ includes
#include <cmath>
using std::log;
using std::pow;

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzState.hpp>

namespace Reaktoro {
namespace {

// Reference temperature of water in units of K
const auto referenceTemperature = 647.27;

// Reference density of water in units of kg/m3
const auto referenceDensity = 317.763;

// Reference volume of water in units of m3/kg
const auto referenceVolume = 1.0/referenceDensity;

// Reference pressure of water in units of Pa
const auto referencePressure = 22.115e+06;

// Reference viscosity of water in units of Pa*s
const auto referenceViscosity = 55.071e-06;

// Reference thermal conductivity of water in units of W/(K*m)
const auto referenceConductivity = 0.49450;

// Reference surface tension of water in units of N/m
const auto referenceSurfTension = 235.8e-03;

// Reference constant for Helmholtz function in units of J/kg
const auto referenceHelmholtz = 69595.89;

// Reference constant for entropy specific heat in units of J/(kg*K)
const auto referenceEntropy = 107.5222;

// Reference constant for sound speed in units of m/s
const auto referenceSoundSpeed = 263.810;

const double A0[] =
{
	-0.130840393653E+2,
	-0.857020420940E+2,
	 0.765192919131E-2,
	-0.620600116069E+0,
	-0.106924329402E+2,
	-0.280671377296E+1,
	 0.119843634845E+3,
	-0.823907389256E+2,
	 0.555864146443E+2,
	-0.310698122980E+2,
	 0.136200239305E+2,
	-0.457116129409E+1,
	 0.115382128188E+1,
	-0.214242224683E+0,
	 0.282800597384E-1,
	-0.250384152737E-2,
	 0.132952679669E-3,
	-0.319277411208E-5
};

const double A1[] =
{
	 0.15383053E+1,
	-0.81048367E+0,
	-0.68305748E+1,
	 0.00000000E+0,
	 0.86756271E+0
};

const double A20 = 0.42923415E+1;

const double yc[] =
{
	 0.59402227E-1,
	-0.28128238E-1,
	 0.56826674E-3,
	-0.27987451E-3
};

const double z0 = 0.317763E+0;

const int ki[] =
{
	1,	1,	1,	1,	2,	2,	2,	2,
	3,	3,	3,	3,	4,	4,	4,	4,
	5,	5,	5,	5,	6,	6,	6,	6,
	7,	7,	7,	7,	9,	9,	9,	9,
	3,	3,	1,	5
};

const int li[] =
{
	1,	2,	4,	6,	1,	2,	4,	6,
	1,	2,	4,	6,	1,	2,	4,	6,
	1,	2,	4,	6,	1,	2,	4,	6,
	1,	2,	4,	6,	1,	2,	4,	6,
	0,	3,	3,	3
};

const double A3[] =
{
	-0.76221190138079E+1,
	 0.32661493707555E+2,
	 0.11305763156821E+2,
	-0.10015404767712E+1,
	 0.12830064355028E+3,
	-0.28371416789846E+3,
	 0.24256279839182E+3,
	-0.99357645626725E+2,
	-0.12275453013171E+4,
	 0.23077622506234E+4,
	-0.16352219929859E+4,
	 0.58436648297764E+3,
	 0.42365441415641E+4,
	-0.78027526961828E+4,
	 0.38855645739589E+4,
	-0.91225112529381E+3,
	-0.90143895703666E+4,
	 0.15196214817734E+5,
	-0.39616651358508E+4,
	-0.72027511617558E+3,
	 0.11147126705990E+5,
	-0.17412065252210E+5,
	 0.99918281207782E+3,
	 0.33504807153854E+4,
	-0.64752644922631E+4,
	 0.98323730907847E+4,
	 0.83877854108422E+3,
	-0.27919349903103E+4,
	 0.11112410081192E+4,
	-0.17287587261807E+4,
	-0.36233262795423E+3,
	 0.61139429010144E+3,
	 0.32968064728562E+2,
	 0.10411239605066E+3,
	-0.38225874712590E+2,
	-0.20307478607599E+3
};

const int mi[] = { 2, 2, 2, 4 };

const int ni[] = { 0, 2, 0, 0 };

const double alpha[] = { 34, 40, 30, 1050 };

const double beta[] = { 20000, 20000, 40000, 25 };

const double ri[] =
{
	0.10038928E+1,
	0.10038928E+1,
	0.10038928E+1,
	0.48778492E+1
};

const double ti[] =
{
	0.98876821E+0,
	0.98876821E+0,
	0.99124013E+0,
	0.41713659E+0
};

const double A4[] =
{
	-0.32329494E-2,
	-0.24139355E-1,
	 0.79027651E-3,
	-0.13362857E+1
};

auto calculateWaterHelmholtzStateHGK0(real t, real d) -> WaterHelmholtzState
{
	WaterHelmholtzState s;

	const auto ln_t = log(t);

	s.helmholtz    = (A0[0] + A0[1] * t) * ln_t;
	s.helmholtzT   =  A0[0]/t + A0[1]*(ln_t + 1);
	s.helmholtzTT  = -A0[0]/(t*t) + A0[1]/t;
	s.helmholtzTTT =  2*A0[0]/(t*t*t) - A0[1]/(t*t);

	for(int i = 2; i <= 17; ++i)
	{
		const auto aux = A0[i] * pow(t, i - 4);

		s.helmholtz    += aux;
		s.helmholtzT   += aux * (i - 4)/t;
		s.helmholtzTT  += aux * (i - 4)*(i - 5)/(t*t);
		s.helmholtzTTT += aux * (i - 4)*(i - 5)*(i - 6)/(t*t*t);
	}

	return s;
}

auto calculateWaterHelmholtzStateHGK1(real t, real d) -> WaterHelmholtzState
{
	WaterHelmholtzState s;

	for(int i = 0; i <= 4; ++i)
	{
		const auto aux = d * A1[i] * pow(t, 1 - i);

		s.helmholtz    += aux;
		s.helmholtzT   -= aux * (i - 1)/t;
		s.helmholtzTT  += aux * (i - 1)*i/(t*t);
		s.helmholtzTTT -= aux * (i - 1)*i*(i + 1)/(t*t*t);
	}

	s.helmholtzD   = s.helmholtz/d;
	s.helmholtzTD  = s.helmholtzT/d;
	s.helmholtzTTD = s.helmholtzTT/d;

	return s;
}

auto calculateWaterHelmholtzStateHGK2(real t, real d) -> WaterHelmholtzState
{
	WaterHelmholtzState s;

	const auto t3   = pow(t, -3);
	const auto t5   = pow(t, -5);
	const auto ln_t = log(t);

	const auto y     = d * (yc[0] + yc[1]*ln_t + yc[2]*t3 + yc[3]*t5);
	const auto y_r   = y/d;
	const auto y_t   = d * (yc[1] - 3.0*yc[2]*t3 - 5.0*yc[3]*t5)/t;
	const auto y_rr  = 0.0;
	const auto y_tt  = d * (-yc[1] + 12.0*yc[2]*t3 + 30.0*yc[3]*t5)/(t*t);
	const auto y_rt  = y_t/d;
	const auto y_rrr = 0.0;
	const auto y_rrt = y_rt/d - y_t/(d*d);
	const auto y_rtt = y_tt/d;
	const auto y_ttt = d * (2*yc[1] - 60*yc[2]*t3 - 210*yc[3]*t5)/(t*t*t);

	const auto x    = 1.0/(1.0 - y);
	const auto x2   = x * x;
	const auto x_r  = y_r * x2;
	const auto x_t  = y_t * x2;
	const auto x_rr = y_rr * x2 + 2.0 * y_r * x_r * x;
	const auto x_tt = y_tt * x2 + 2.0 * y_t * x_t * x;
	const auto x_rt = y_rt * x2 + 2.0 * y_r * x_t * x;
	const auto x_rrr = y_rrr*x2 + 4*y_rr*x_r*x + 2*y_r*(x_rr*x + x_r*x_r);
	const auto x_rrt = y_rrt*x2 + 2*(y_rt*x_r + y_rr*x_t) + 2*y_r*(x_rt*x + x_r*x_t);
	const auto x_rtt = y_rtt*x2 + 4*y_rt*x_t*x + 2*y_r*(x_tt*x + x_t*x_t);
	const auto x_ttt = y_ttt*x2 + 4*y_tt*x_t*x + 2*y_t*(x_tt*x + x_t*x_t);

	const auto u     = log(d * x);
	const auto u_r   = x_r/x + 1.0/d;
	const auto u_t   = x_t/x;
	const auto u_rr  = x_rr/x - x_r*x_r/(x*x) - 1.0/(d*d);
	const auto u_rt  = x_rt/x - x_r*x_t/(x*x);
	const auto u_tt  = x_tt/x - x_t*x_t/(x*x);
	const auto u_rrr = x_rrr/x - 3*x_rr*x_r/(x*x) + 2*x_r*x_r*x_r/(x*x*x) + 2/(d*d*d);
	const auto u_rrt = x_rrt/x - (2*x_rt*x_r + x_rr*x_t)/(x*x) + 2*x_r*x_r*x_t/(x*x*x);
	const auto u_rtt = x_rtt/x - (2*x_rt*x_t + x_tt*x_r)/(x*x) + 2*x_t*x_t*x_r/(x*x*x);
	const auto u_ttt = x_ttt/x - 3*x_tt*x_t/(x*x) + 2*x_t*x_t*x_t/(x*x*x);

	const auto c1 = -130.0/3.0;
	const auto c2 =  169.0/6.0;
	const auto c3 = -14.0;

	s.helmholtz    = A20 * t * (u + c1*x  + c2*x*x + c3*y);
	s.helmholtzD   = A20 * t * (u_r + c1*x_r + 2*c2*x*x_r + c3*y_r);
	s.helmholtzT   = A20 * t * (u_t + c1*x_t + 2*c2*x*x_t + c3*y_t) + s.helmholtz/t;
	s.helmholtzDD  = A20 * t * (u_rr + c1*x_rr + 2*c2*(x*x_rr + x_r*x_r) + c3*y_rr);
	s.helmholtzTD  = A20 * t * (u_rt + c1*x_rt + 2*c2*(x*x_rt + x_r*x_t) + c3*y_rt) + s.helmholtzD/t;
	s.helmholtzTT  = A20 * t * (u_tt + c1*x_tt + 2*c2*(x*x_tt + x_t*x_t) + c3*y_tt) + 2*(s.helmholtzT/t - s.helmholtz/(t*t));
	s.helmholtzDDD = A20 * t * (u_rrr + c1*x_rrr + 2*c2*(3*x_r*x_rr + x*x_rrr) + c3*y_rrr);
	s.helmholtzTDD = A20 * t * (u_rrt + c1*x_rrt + 2*c2*(x_t*x_rr + 2*x_r*x_rt + x*x_rrt) + c3*y_rrt) + s.helmholtzDD/t;
	s.helmholtzTTD = A20 * t * (u_rtt + c1*x_rtt + 2*c2*(x_r*x_tt + 2*x_t*x_rt + x*x_rtt) + c3*y_rtt) + 2*(s.helmholtzTD - s.helmholtzD/t)/t;
	s.helmholtzTTT = A20 * t * (u_ttt + c1*x_ttt + 2*c2*(3*x_t*x_tt + x*x_ttt) + c3*y_ttt) + 3*(s.helmholtzTT - 2*s.helmholtzT/t + s.helmholtz/(t*t))/t;

	return s;
}

auto calculateWaterHelmholtzStateHGK3(real t, real d) -> WaterHelmholtzState
{
	WaterHelmholtzState s;

	const auto z     =  1 - exp(-z0 * d);
	const auto z_r   =  z0 * (1 - z);
	const auto z_rr  = -z0 * z_r;
	const auto z_rrr = -z0 * z_rr;

	for(int i = 0; i <= 35; ++i)
	{
		const auto lambda     =  A3[i] * pow(t, -li[i]) * pow(z, ki[i]);
		const auto lambda_r   =  ki[i]*z_r*lambda/z;
		const auto lambda_t   = -li[i]*lambda/t;
		const auto lambda_rr  =  lambda_r*(z_rr/z_r + lambda_r/lambda - z_r/z);
		const auto lambda_rt  =  lambda_r*lambda_t/lambda;
		const auto lambda_tt  =  lambda_t*(lambda_t/lambda - 1.0/t);
		const auto lambda_rrr =  lambda_rr*(z_rr/z_r + lambda_r/lambda - z_r/z) + lambda_r*(z_rrr/z_r - pow(z_rr/z_r, 2) + lambda_rr/lambda - pow(lambda_r/lambda, 2) - z_rr/z + pow(z_r/z, 2));
		const auto lambda_rrt = -pow(lambda_r/lambda, 2)*lambda_t + (lambda_rr*lambda_t + lambda_rt*lambda_r)/lambda;
		const auto lambda_rtt = -pow(lambda_t/lambda, 2)*lambda_r + (lambda_tt*lambda_r + lambda_rt*lambda_t)/lambda;
		const auto lambda_ttt =  lambda_tt * (lambda_t/lambda - 1.0/t) + lambda_t*(lambda_tt/lambda - pow(lambda_t/lambda, 2) + 1.0/(t*t));

		s.helmholtz    += lambda;
		s.helmholtzD   += lambda_r;
		s.helmholtzT   += lambda_t;
		s.helmholtzDD  += lambda_rr;
		s.helmholtzTD  += lambda_rt;
		s.helmholtzTT  += lambda_tt;
		s.helmholtzDDD += lambda_rrr;
		s.helmholtzTDD += lambda_rrt;
		s.helmholtzTTD += lambda_rtt;
		s.helmholtzTTT += lambda_ttt;
	}

	return s;
}

auto calculateWaterHelmholtzStateHGK4(real t, real d) -> WaterHelmholtzState
{
	WaterHelmholtzState s;

	for(int i = 0; i <= 3; ++i)
	{
		const auto delta   = (d - ri[i])/ri[i];
		const auto tau     = (t - ti[i])/ti[i];
		const auto delta_r = 1.0/ri[i];
		const auto tau_t   = 1.0/ti[i];

		const auto delta_m = pow(delta, mi[i]);
		const auto delta_n = pow(delta, ni[i]);

		const auto psi    = (ni[i] - alpha[i]*mi[i]*delta_m)*delta_r/delta;
		const auto psi_r  = -(ni[i] + alpha[i]*mi[i]*(mi[i] - 1)*delta_m)*pow(delta_r/delta, 2);
		const auto psi_rr = (2*ni[i] - alpha[i]*mi[i]*(mi[i] - 1)*(mi[i] - 2)*delta_m)*pow(delta_r/delta, 3);

		const auto theta     =  A4[i]*delta_n*exp(-alpha[i]*delta_m - beta[i]*tau*tau);
		const auto theta_r   =  psi*theta;
		const auto theta_t   = -2*beta[i]*tau*tau_t*theta;
		const auto theta_rr  =  psi_r*theta + psi*theta_r;
		const auto theta_tt  =  2*beta[i]*(2*beta[i]*tau*tau - 1)*tau_t*tau_t*theta;
		const auto theta_rt  = -2*beta[i]*tau*tau_t*theta_r;
		const auto theta_rrr =  psi_rr*theta + 2*psi_r*theta_r + psi*theta_rr;
		const auto theta_rrt =  psi_r*theta_r + psi*theta_rt;
		const auto theta_rtt =  psi*theta_tt;
		const auto theta_ttt = -2*beta[i]*(2*tau_t*tau_t*theta_t + tau*tau_t*theta_tt);

		s.helmholtz    += theta;
		s.helmholtzD   += theta_r;
		s.helmholtzT   += theta_t;
		s.helmholtzDD  += theta_rr;
		s.helmholtzTD  += theta_rt;
		s.helmholtzTT  += theta_tt;
		s.helmholtzDDD += theta_rrr;
		s.helmholtzTDD += theta_rrt;
		s.helmholtzTTD += theta_rtt;
		s.helmholtzTTT += theta_ttt;
	}

	return s;
}

} // namespace

auto waterHelmholtzStateHGK(real T, real D) -> WaterHelmholtzState
{
	// The dimensionless temperature and density
	const auto t = T/referenceTemperature;
	const auto r = D/referenceDensity;

	// Compute the contributions from each auxiliary Helmholtz state
	WaterHelmholtzState aux0, aux1, aux2, aux3, aux4, res;
	aux0 = calculateWaterHelmholtzStateHGK0(t, r);
	aux1 = calculateWaterHelmholtzStateHGK1(t, r);
	aux2 = calculateWaterHelmholtzStateHGK2(t, r);
	aux3 = calculateWaterHelmholtzStateHGK3(t, r);
	aux4 = calculateWaterHelmholtzStateHGK4(t, r);

	// Assemble the contributions from each auxiliary Helmholtz state
	res.helmholtz    = aux0.helmholtz    + aux1.helmholtz    + aux2.helmholtz    + aux3.helmholtz    + aux4.helmholtz;
	res.helmholtzD   = aux0.helmholtzD   + aux1.helmholtzD   + aux2.helmholtzD   + aux3.helmholtzD   + aux4.helmholtzD;
	res.helmholtzT   = aux0.helmholtzT   + aux1.helmholtzT   + aux2.helmholtzT   + aux3.helmholtzT   + aux4.helmholtzT;
	res.helmholtzDD  = aux0.helmholtzDD  + aux1.helmholtzDD  + aux2.helmholtzDD  + aux3.helmholtzDD  + aux4.helmholtzDD;
	res.helmholtzTD  = aux0.helmholtzTD  + aux1.helmholtzTD  + aux2.helmholtzTD  + aux3.helmholtzTD  + aux4.helmholtzTD;
	res.helmholtzTT  = aux0.helmholtzTT  + aux1.helmholtzTT  + aux2.helmholtzTT  + aux3.helmholtzTT  + aux4.helmholtzTT;
	res.helmholtzDDD = aux0.helmholtzDDD + aux1.helmholtzDDD + aux2.helmholtzDDD + aux3.helmholtzDDD + aux4.helmholtzDDD;
	res.helmholtzTDD = aux0.helmholtzTDD + aux1.helmholtzTDD + aux2.helmholtzTDD + aux3.helmholtzTDD + aux4.helmholtzTDD;
	res.helmholtzTTD = aux0.helmholtzTTD + aux1.helmholtzTTD + aux2.helmholtzTTD + aux3.helmholtzTTD + aux4.helmholtzTTD;
	res.helmholtzTTT = aux0.helmholtzTTT + aux1.helmholtzTTT + aux2.helmholtzTTT + aux3.helmholtzTTT + aux4.helmholtzTTT;

	// Convert the Helmholtz free energy of water and its derivatives to dimensioned form
	res.helmholtz    *= referenceHelmholtz;
	res.helmholtzD   *= referenceHelmholtz/referenceDensity;
	res.helmholtzT   *= referenceHelmholtz/referenceTemperature;
	res.helmholtzDD  *= referenceHelmholtz/referenceDensity/referenceDensity;
	res.helmholtzTD  *= referenceHelmholtz/referenceDensity/referenceTemperature;
	res.helmholtzTT  *= referenceHelmholtz/referenceTemperature/referenceTemperature;
	res.helmholtzDDD *= referenceHelmholtz/referenceDensity/referenceDensity/referenceDensity;
	res.helmholtzTDD *= referenceHelmholtz/referenceDensity/referenceDensity/referenceTemperature;
	res.helmholtzTTD *= referenceHelmholtz/referenceDensity/referenceTemperature/referenceTemperature;
	res.helmholtzTTT *= referenceHelmholtz/referenceTemperature/referenceTemperature/referenceTemperature;

	return res;
}

} // namespace Reaktoro
