// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "WaterThermoStateWagnerPruss.hpp"

// C++ includes
#include <cmath>
using std::pow;
using std::log;

// Reaktor includes
#include <Reaktor/Thermodynamics/WaterHelmholtzState.hpp>
#include <Reaktor/Thermodynamics/WaterConstants.hpp>

namespace Reaktor {
namespace internal {

const double no[] =
{
	0, -8.32044648201, 6.6832105268, 3.00632, 0.012436, 0.97315, 1.27950, 0.96956, 0.24873
};

const double gammao[] =
{
    1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105
};

const double c[] =
{
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 6, 6, 6, 6, 0, 0, 0
};

const double d[] =
{
    0, 1, 1, 1, 2, 2, 3, 4, 1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9,
    10, 11, 13, 15, 1, 2, 2, 2, 3, 4, 4, 4, 5, 6, 6, 7, 9, 9, 9,
    9, 9, 10, 10, 12, 3, 4, 4, 5, 14, 3, 6, 6, 6, 3, 3, 3
};


const double t[] =
{
    0, -0.5, 0.875, 1, 0.5, 0.75, 0.375, 1, 4, 6, 12, 1, 5, 4, 2, 13, 9, 3,
    4, 11, 4, 13, 1, 7, 1, 9, 10, 10, 3, 7, 10, 10, 6, 10, 10, 1, 2, 3, 4,
    8, 6, 9, 8, 16, 22, 23, 23, 10, 50, 44, 46, 50, 0, 1, 4
};

const double n[] =
{
	 0,
	 0.12533547935523e-01,
	 0.78957634722828e+01,
	-0.87803203303561e+01,
	 0.31802509345418,
	-0.26145533859358,
	-0.78199751687981e-02,
	 0.88089493102134e-02,
	-0.66856572307965,
	 0.20433810950965,
	-0.66212605039687e-04,
	-0.19232721156002,
	-0.25709043003438,
	 0.16074868486251,
	-0.40092828925807e-01,
	 0.39343422603254e-06,
	-0.75941377088144e-05,
	 0.56250979351888e-03,
	-0.15608652257135e-04,
	 0.11537996422951e-08,
	 0.36582165144204e-06,
	-0.13251180074668e-11,
	-0.62639586912454e-09,
	-0.10793600908932,
	 0.17611491008752e-01,
	 0.22132295167546,
	-0.40247669763528,
	 0.58083399985759,
	 0.49969146990806e-02,
	-0.31358700712549e-01,
	-0.74315929710341,
	 0.47807329915480,
	 0.20527940895948e-01,
	-0.13636435110343,
	 0.14180634400617e-01,
	 0.83326504880713e-02,
	-0.29052336009585e-01,
	 0.38615085574206e-01,
	-0.20393486513704e-01,
	-0.16554050063734e-02,
	 0.19955571979541e-02,
	 0.15870308324157e-03,
	-0.16388568342530e-04,
	 0.43613615723811e-01,
	 0.34994005463765e-01,
	-0.76788197844621e-01,
	 0.22446277332006e-01,
	-0.62689710414685e-04,
	-0.55711118565645e-09,
	-0.19905718354408,
	 0.31777497330738,
	-0.11841182425981,
	-0.31306260323435e+02,
	 0.31546140237781e+02,
	-0.25213154341695e+04,
	-0.14874640856724,
	 0.31806110878444
};

const double alpha[] = { 20, 20, 20 };

const double beta[] = { 150, 150, 250 };

const double gamma[] = { 1.21, 1.21, 1.25 };

const double epsilon[] = { 1, 1, 1 };

const double a[] = { 3.5, 3.5 };

const double b[] = { 0.85, 0.95 };

const double A[] = { 0.32, 0.32 };

const double B[] = { 0.2, 0.2 };

const double C[] = { 28, 32 };

const double F[] = { 700, 800 }; // D has been replaced by F to avoid conflicts

const double E[] = { 0.3, 0.3 };

double delta(double D)
{
	return D/waterCriticalDensity;
}

double tau(double T)
{
	return waterCriticalTemperature/T;
}

} /* namespace internal */

using namespace internal;

auto waterHelmholtzWagnerPruss(double T, double D) -> WaterHelmholtzState
{
	const double tau   = internal::tau(T);
	const double delta = internal::delta(D);

	double phio     =  log(delta) + no[1] + no[2]*tau + no[3]*log(tau);
	double phio_d   =  1.0/delta;
	double phio_t   =  no[2] + no[3]/tau;
	double phio_dd  = -1.0/pow(delta, 2);
	double phio_tt  = -no[3]/pow(tau, 2);
	double phio_dt  =  0.0;
	double phio_ddd =  2.0/pow(delta, 3);
	double phio_ttt =  2.0*no[3]/pow(tau, 3);
	double phio_dtt =  0.0;
	double phio_ddt =  0.0;

	for(int i = 4; i <= 8; ++i)
	{
		const int j = i - 4;

		const double ee = exp(gammao[j] * tau);

		phio     += no[i] * log(1.0 - 1.0/ee);
		phio_t   += no[i] * (gammao[j]/(ee - 1));
		phio_tt  -= no[i] * ee * pow((gammao[j]/(ee - 1)), 2);
		phio_ttt += no[i] * ee * (1 + ee) * pow((gammao[j]/(ee - 1)), 3);
	}

	double phir     = 0.0;
	double phir_d   = 0.0;
	double phir_t   = 0.0;
	double phir_dd  = 0.0;
	double phir_tt  = 0.0;
	double phir_dt  = 0.0;
	double phir_ddd = 0.0;
	double phir_ttt = 0.0;
	double phir_dtt = 0.0;
	double phir_ddt = 0.0;

	for(int i = 1; i <= 7; ++i)
	{
		const double A     = n[i]*pow(delta, d[i])*pow(tau, t[i]);
		const double A_d   = d[i]/delta * A;
		const double A_t   = t[i]/tau * A;
		const double A_dd  = (d[i] - 1)/delta * A_d;
		const double A_tt  = (t[i] - 1)/tau * A_t;
		const double A_dt  = t[i]*d[i]/(tau*delta) * A;
		const double A_ddd = (d[i] - 2)/delta * A_dd;
		const double A_ttt = (t[i] - 2)/tau * A_tt;
		const double A_dtt = d[i]/delta * A_tt;
		const double A_ddt = t[i]/tau * A_dd;

		phir     += A;
		phir_d   += A_d;
		phir_t   += A_t;
		phir_dd  += A_dd;
		phir_tt  += A_tt;
		phir_dt  += A_dt;
		phir_ddd += A_ddd;
		phir_ttt += A_ttt;
		phir_dtt += A_dtt;
		phir_ddt += A_ddt;
	}

	for(int i = 8; i <= 51; ++i)
	{
		const double dci = pow(delta, c[i]);

		const double B     =  n[i]*pow(delta, d[i])*pow(tau, t[i])*exp(-dci);
		const double B_d   = (d[i] - c[i]*dci)/delta * B;
		const double B_t   =  t[i]/tau * B;
		const double B_dd  = (d[i] - c[i]*dci - 1)/delta * B_d - dci*pow(c[i]/delta, 2) * B;
		const double B_tt  = (t[i] - 1)/tau * B_t;
		const double B_dt  =  t[i]/tau * B_d;
		const double B_ddd = (d[i] - c[i]*dci - 1)/delta * B_dd - ((d[i] - c[i]*dci - 1) + 2*c[i]*c[i]*dci)/pow(delta, 2) * B_d - c[i]*c[i]*dci*(c[i] - 2)/pow(delta, 3) * B;
		const double B_ttt = (t[i] - 2)/tau * B_tt;
		const double B_dtt = (t[i] - 1)/tau * B_dt;
		const double B_ddt = (d[i] - c[i]*dci - 1)/delta * B_dt - c[i]*c[i]*dci/pow(delta, 2) * B_t;

		phir     += B;
		phir_d   += B_d;
		phir_t   += B_t;
		phir_dd  += B_dd;
		phir_tt  += B_tt;
		phir_dt  += B_dt;
		phir_ddd += B_ddd;
		phir_ttt += B_ttt;
		phir_dtt += B_dtt;
		phir_ddt += B_ddt;
	}

	for(int i = 52; i <= 54; ++i)
	{
		const int j = i - 52;

		const double aux1d = (d[i]/delta - 2*alpha[j]*(delta - epsilon[j]));
		const double aux1t = (t[i]/tau - 2*beta[j]*(tau - gamma[j]));

		const double aux2d = (d[i]/pow(delta, 2) + 2*alpha[j]);
		const double aux2t = (t[i]/pow(tau, 2) + 2*beta[j]);

		const double C     = n[i]*pow(delta, d[i])*pow(tau, t[i])*exp(-alpha[j]*pow(delta - epsilon[j], 2) - beta[j]*pow(tau - gamma[j], 2));
		const double C_d   = aux1d * C;
		const double C_t   = aux1t * C;
		const double C_dd  = aux1d * C_d - aux2d * C;
		const double C_tt  = aux1t * C_t - aux2t * C;
		const double C_dt  = aux1d * aux1t * C;
		const double C_ddd = aux1d * C_dd - 2*aux2d * C_d + 2*d[i]/pow(delta, 3) * C;
		const double C_ttt = aux1t * C_tt - 2*aux2t * C_t + 2*t[i]/pow(tau, 3) * C;
		const double C_dtt = aux1t * C_dt - aux2t * C_d;
		const double C_ddt = aux1d * C_dt - aux2d * C_t;

		phir     += C;
		phir_d   += C_d;
		phir_t   += C_t;
		phir_dd  += C_dd;
		phir_tt  += C_tt;
		phir_dt  += C_dt;
		phir_ddd += C_ddd;
		phir_ttt += C_ttt;
		phir_dtt += C_dtt;
		phir_ddt += C_ddt;
	}

	for(int i = 55; i <= 56; ++i)
	{
		const int j = i - 55;

		const double dd = pow(delta - 1, 2);
		const double tt = pow(tau - 1, 2);

		const double theta     = (1 - tau) + A[j]*pow(dd, 0.5/E[j]);
		const double theta_d   = (theta + tau - 1)/(delta - 1)/E[j];
		const double theta_dd  = (1.0/E[j] - 1) * theta_d/(delta - 1);
		const double theta_ddd = (1.0/E[j] - 1) * (theta_dd/(delta - 1) - theta_d/dd);

		const double psi     = exp(-C[j]*dd - F[j]*tt);
		const double psi_d   = -2*C[j]*(delta - 1) * psi;
		const double psi_t   = -2*F[j]*(tau - 1) * psi;
		const double psi_dd  = -2*C[j]*(psi + (delta - 1) * psi_d);
		const double psi_tt  = -2*F[j]*(psi + (tau - 1) * psi_t);
		const double psi_dt  =  4*C[j]*F[j]*(delta - 1)*(tau - 1) * psi;
		const double psi_ddd = -2*C[j]*(2*psi_d + (delta - 1) * psi_dd);
		const double psi_ttt = -2*F[j]*(2*psi_t + (tau - 1) * psi_tt);
		const double psi_dtt = -2*F[j]*(psi_d + (tau - 1) * psi_dt);
		const double psi_ddt = -2*C[j]*(psi_t + (delta - 1) * psi_dt);

		const double Delta     = theta*theta + B[j]*pow(dd, a[j]);
		const double Delta_d   = 2*(theta*theta_d + a[j]*(Delta - theta*theta)/(delta - 1));
		const double Delta_t   = -2*theta;
		const double Delta_dd  = 2*(theta_d*theta_d + theta*theta_dd + a[j] * ((Delta_d - 2*theta*theta_d)/(delta - 1) - (Delta - theta*theta)/pow(delta - 1, 2)));
		const double Delta_tt  = 2;
		const double Delta_dt  = -2*theta_d;
		const double Delta_ddd = 2*(3*theta_d*theta_dd + theta*theta_ddd + a[j] * ((Delta_dd - 2*theta_d*theta_d - 2*theta*theta_dd)/(delta - 1) - 2*(Delta_d - 2*theta*theta_d)/pow(delta - 1, 2) + 2*(Delta - theta*theta)/pow(delta - 1, 3)));
		const double Delta_ttt = 0;
		const double Delta_dtt = 0;
		const double Delta_ddt = -2*theta_dd;

		const double DeltaPow     =  pow(Delta, b[j]);
		const double DeltaPow_d   =  b[j]*Delta_d/Delta * DeltaPow;
		const double DeltaPow_t   =  b[j]*Delta_t/Delta * DeltaPow;
		const double DeltaPow_dd  = (b[j]*Delta_dd/Delta + b[j]*(b[j] - 1)*pow(Delta_d/Delta, 2)) * DeltaPow;
		const double DeltaPow_tt  = (b[j]*Delta_tt/Delta + b[j]*(b[j] - 1)*pow(Delta_t/Delta, 2)) * DeltaPow;
		const double DeltaPow_dt  = (b[j]*Delta_dt/Delta + b[j]*(b[j] - 1)*Delta_d*Delta_t/Delta/Delta) * DeltaPow;
		const double DeltaPow_ddd = (b[j]*Delta_ddd/Delta + 3*b[j]*(b[j] - 1)*Delta_d*Delta_dd/Delta/Delta + b[j]*(b[j] - 1)*(b[j] - 2)*pow(Delta_d/Delta, 3)) * DeltaPow;
		const double DeltaPow_ttt = (b[j]*Delta_ttt/Delta + 3*b[j]*(b[j] - 1)*Delta_t*Delta_tt/Delta/Delta + b[j]*(b[j] - 1)*(b[j] - 2)*pow(Delta_t/Delta, 3)) * DeltaPow;
		const double DeltaPow_dtt = (b[j]*Delta_dtt/Delta + b[j]*(b[j] - 1)*(Delta_d*Delta_tt + 2*Delta_t*Delta_dt)/Delta/Delta + b[j]*(b[j] - 1)*(b[j] - 2)*Delta_t*Delta_t*Delta_d/pow(Delta, 3)) * DeltaPow;
		const double DeltaPow_ddt = (b[j]*Delta_ddt/Delta + b[j]*(b[j] - 1)*(Delta_t*Delta_dd + 2*Delta_d*Delta_dt)/Delta/Delta + b[j]*(b[j] - 1)*(b[j] - 2)*Delta_d*Delta_d*Delta_t/pow(Delta, 3)) * DeltaPow;

		const double D     = n[i]*DeltaPow*delta*psi;
		const double D_d   = n[i]*(DeltaPow*(psi + delta*psi_d) + DeltaPow_d*delta*psi);
		const double D_t   = n[i]*delta*(DeltaPow_t*psi + DeltaPow*psi_t);
		const double D_dd  = n[i]*(DeltaPow*(2*psi_d + delta*psi_dd) + 2*DeltaPow_d*(psi + delta*psi_d) + DeltaPow_dd*delta*psi);
		const double D_tt  = n[i]*delta*(DeltaPow_tt*psi + 2*DeltaPow_t*psi_t + DeltaPow*psi_tt);
		const double D_dt  = n[i]*(DeltaPow*(psi_t + delta*psi_dt) + delta*DeltaPow_d*psi_t + DeltaPow_t*(psi + delta*psi_d) + DeltaPow_dt*delta*psi);
		const double D_ddd = n[i]*(DeltaPow_ddd*delta*psi + 3*DeltaPow_dd*(psi + delta*psi_d) + 3*DeltaPow_d*(2*psi_d + delta*psi_dd) + DeltaPow*(3*psi_dd + delta*psi_ddd));
		const double D_ttt = n[i]*delta*(DeltaPow_ttt*psi + 3*DeltaPow_tt*psi_t + 3*DeltaPow_t*psi_tt + DeltaPow*psi_ttt);
		const double D_dtt = n[i]*(DeltaPow_tt*psi + 2*DeltaPow_t*psi_t + DeltaPow*psi_tt) + n[i]*delta*(DeltaPow_dtt*psi + DeltaPow_tt*psi_d + 2*DeltaPow_dt*psi_t + 2*DeltaPow_t*psi_dt + DeltaPow_d*psi_tt + DeltaPow*psi_dtt);
		const double D_ddt = n[i]*(DeltaPow_ddt*delta*psi + 2*DeltaPow_dt*(psi + delta*psi_d) + DeltaPow_dd*delta*psi_t + DeltaPow_t*(2*psi_d + delta*psi_dd) + 2*DeltaPow_d*(psi_t + delta*psi_dt) + DeltaPow*(2*psi_dt + delta*psi_ddt));

		phir     += D;
		phir_d   += D_d;
		phir_t   += D_t;
		phir_dd  += D_dd;
		phir_tt  += D_tt;
		phir_dt  += D_dt;
		phir_ddd += D_ddd;
		phir_ttt += D_ttt;
		phir_dtt += D_dtt;
		phir_ddt += D_ddt;
	}

	const double phi     = phio     + phir    ;
	const double phi_d   = phio_d   + phir_d  ;
	const double phi_t   = phio_t   + phir_t  ;
	const double phi_dd  = phio_dd  + phir_dd ;
	const double phi_tt  = phio_tt  + phir_tt ;
	const double phi_dt  = phio_dt  + phir_dt ;
	const double phi_ddd = phio_ddd + phir_ddd;
	const double phi_ttt = phio_ttt + phir_ttt;
	const double phi_dtt = phio_dtt + phir_dtt;
	const double phi_ddt = phio_ddt + phir_ddt;

	const double Tcr = waterCriticalTemperature;
	const double Dcr = waterCriticalDensity;

	const double tT   = -Tcr/(T*T);
	const double tTT  =  2*Tcr/(T*T*T);
	const double tTTT = -6*Tcr/(T*T*T*T);
	const double dD   =  1/Dcr;

	const double phiT   = phi_t*tT;
	const double phiD   = phi_d*dD;
	const double phiTT  = phi_tt*tT*tT + phi_t*tTT;
	const double phiTD  = phi_dt*tT*dD;
	const double phiDD  = phi_dd*dD*dD;
	const double phiTTT = phi_ttt*tT*tT*tau*T + 3*phi_tt*tT*tTT + phi_t*tTTT;
	const double phiTTD = phi_dtt*tT*tT*dD + phi_dt*tTT*dD;
	const double phiTDD = phi_ddt*tT*dD*dD;
	const double phiDDD = phi_ddd*dD*dD*dD;

	// The specific gas constant in units of J/(kg*K)
	const double R = 461.51805;

	WaterHelmholtzState hs;

	hs.helmholtz    = R*T*phi;
	hs.helmholtzT   = R*T*phiT + R*phi;
	hs.helmholtzD   = R*T*phiD;
	hs.helmholtzTT  = R*T*phiTT + 2*R*phiT;
	hs.helmholtzTD  = R*T*phiTD + R*phiD;
	hs.helmholtzDD  = R*T*phiDD;
	hs.helmholtzTTT = R*T*phiTTT + 3*R*phiTT;
	hs.helmholtzTTD = R*T*phiTTD + 2*R*phiTD;
	hs.helmholtzTDD = R*T*phiTDD + R*phiDD;
	hs.helmholtzDDD = R*T*phiDDD;

	return hs;
}

auto saturatedPressureWaterWagnerPruss(double T) -> double
{
	const double a1 = -7.85951783;
	const double a2 =  1.84408259;
	const double a3 = -11.7866497;
	const double a4 =  22.6807411;
	const double a5 = -15.9618719;
	const double a6 =  1.80122502;

	const double Tcr = waterCriticalTemperature;
	const double Pcr = waterCriticalPressure;

	const double t   = 1 - T/Tcr;
	const double t15 = pow(t, 1.5);
	const double t30 = t15 * t15;
	const double t35 = t15 * t * t;
	const double t40 = t30 * t;
	const double t75 = t35 * t40;

	return Pcr * exp(Tcr/T * (a1*t + a2*t15 + a3*t30 + a4*t35 + a5*t40 + a6*t75));
}

auto saturatedLiquidDensityWaterWagnerPruss(double T) -> double
{
	const double b1 =  1.99274064;
	const double b2 =  1.09965342;
	const double b3 = -0.510839303;
	const double b4 = -1.75493479;
	const double b5 = -45.5170352;
	const double b6 = -6.74694450e+05;

	const double Tcr = waterCriticalTemperature;
	const double Dcr = waterCriticalDensity;

	const double t     = 1 - T/Tcr;
	const double t13   = pow(t, 1./3);
	const double t23   = t13 * t13;
	const double t53   = t13 * t23 * t23;
	const double t163  = t13 * t53 * t53 * t53;
	const double t433  = t163 * t163 * t53 * t * t;
	const double t1103 = t433 * t433 * t163 * t53 * t;

	return Dcr * (1 + b1*t13 + b2*t23 + b3*t53 + b4*t163 + b5*t433 + b6*t1103);
}

auto saturatedVapourDensityWaterWagnerPruss(double T) -> double
{
	const double c1 = -2.03150240;
	const double c2 = -2.68302940;
	const double c3 = -5.38626492;
	const double c4 = -17.2991605;
	const double c5 = -44.7586581;
	const double c6 = -63.9201063;

	const double Tcr = waterCriticalTemperature;
	const double Dcr = waterCriticalDensity;

	const double t    = 1 - T/Tcr;
	const double t16  = pow(t, 1./6);
	const double t26  = t16 * t16;
	const double t46  = t26 * t26;
	const double t86  = t46 * t46;
	const double t186 = t86 * t86 * t26;
	const double t376 = t186 * t186 * t16;
	const double t716 = t376 * t186 * t86 * t86;

	return Dcr * exp(Tcr/T * (c1*t26 + c2*t46 + c3*t86 + c4*t186 + c5*t376 + c6*t716));
}

} // namespace Reaktor
