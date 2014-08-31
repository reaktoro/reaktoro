/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

//
//#include "GaseousActivityDuanMollerWeare.hpp"
//
//// C++ includes
//#include <cmath>
//#include <stdexcept>
//#include <iostream>
//
//namespace Reaktor {
//
///**
// * References:
// *   1) Duan et al. (1992), An equation of state for the CH4-CO2-H2O system: I. Pure systems from 0 to 1000°C and 0 to 8000 bar
// *   2) Span and Wagner (1994), A New Equation of State for Carbon Dioxide Covering the Fluid Region from the Triple-Point Temperature to 1100 K at Pressures up to 800 MPa
// *   3) Wagner and Pruss (2002), The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use
// */
//
//const double TcrCO2 = 304.20;      // the critical temperature of carbon dioxide in units of kelvin
//const double TcrH2O = 647.25;      // the critical temperature of water in units of kelvin
//const double TcrCH4 = 190.60;      // the critical temperature of methane in units of kelvin
//
//const double PcrCO2 = 73.825e+05;  // the critical pressure of carbon dioxide in units of pascal
//const double PcrH2O = 221.19e+05;  // the critical pressure of water in units of pascal
//const double PcrCH4 = 46.41e+05;   // the critical pressure of methane in units of pascal
//
//const double DcrCO2 = 467.6;       // the critical density of carbon dioxide in units of kg/m3
//const double DcrH2O = 322.0;       // the critical density of water in units of kg/m3
//const double DcrCH4 = 162.0;       // the critical density of methane in units of kg/m3
//
//const double mwCO2  = 44.0100e-03; // the molecular weight of carbon dioxide in units of kg/mol
//const double mwH2O  = 18.0153e-03; // the molecular weight of water in units of kg/mol
//const double mwCH4  = 16.0400e-03; // the molecular weight of methane in units of kg/mol
//
//const double aH2O[] =
//{
//	 8.64449220e-02,
//	-3.96918955e-01,
//	-5.73334886e-02,
//	-2.93893000e-04,
//	-4.15775512e-03,
//	 1.99496791e-02,
//	 1.18901426e-04,
//	 1.55212063e-04,
//	-1.06855859e-04,
//	-4.93197687e-06,
//	-2.73739155e-06,
//	 2.65571238e-06,
//	 8.96079018e-03,
//	 4.02000000e+00,
//	 2.57000000e-02
//};
//
//const double aCO2[] =
//{
//	 8.99288497e-02,
//	-4.94783127e-01,
//	 4.77922245e-02,
//	 1.03808883e-02,
//	-2.82516861e-02,
//	 9.49887563e-02,
//	 5.20600880e-04,
//	-2.93540971e-04,
//	-1.77265112e-03,
//	-2.51101973e-05,
//	 8.93353441e-05,
//	 7.88998563e-05,
//	-1.66727022e-02,
//	 1.39800000e+00,
//	 2.96000000e-02
//};
//
//const double aCH4[] =
//{
//	 8.72553928e-02,
//	-7.52599476e-01,
//	 3.75419887e-01,
//	 1.07291342e-02,
//	 5.49626360e-03,
//	-1.84772802e-02,
//	 3.18993183e-04,
//	 2.11079375e-04,
//	 2.01682801e-05,
//	-1.65606189e-05,
//	 1.19614546e-04,
//	-1.08087289e-04,
//	 4.48262295e-02,
//	 7.53970000e-01,
//	 7.71670000e-02
//};
//
///// Calculate the saturated pressure of water in units of pascal
//double CalculateSaturatedPressureWater(double T)
//{
//	const double Tcr = 647.096;    // critical temperature of water in units of kelvin
//	const double Pcr = 22.064e+06; // critical pressure of water in units of pascal
//
//	const double a1 = -7.85951783;
//	const double a2 =  1.84408259;
//	const double a3 = -11.7866497;
//	const double a4 =  22.6807411;
//	const double a5 = -15.9618719;
//	const double a6 =  1.80122502;
//
//	const double t   = 1 - T/Tcr;
//	const double t15 = std::pow(t, 1.5);
//	const double t30 = t15 * t15;
//	const double t35 = t15 * t * t;
//	const double t40 = t30 * t;
//	const double t75 = t35 * t40;
//
//	return Pcr * std::exp(Tcr/T * (a1*t + a2*t15 +
//		a3*t30 + a4*t35 + a5*t40 + a6*t75));
//}
//
///// Calculate the saturated liquid density of water in units of kg/m3
//double CalculateSaturatedLiquidDensityWater(double T)
//{
//	const double Tcr = 647.096; // critical temperature of water in units of kelvin
//	const double Dcr = 322.0;   // critical density of water in units of kg/m3
//
//	const double b1 =  1.99274064;
//	const double b2 =  1.09965342;
//	const double b3 = -0.510839303;
//	const double b4 = -1.75493479;
//	const double b5 = -45.5170352;
//	const double b6 = -6.74694450e+05;
//
//	const double t     = 1 - T/Tcr;
//	const double t13   = std::pow(t, 1./3);
//	const double t23   = t13 * t13;
//	const double t53   = t13 * t23 * t23;
//	const double t163  = t13 * t53 * t53 * t53;
//	const double t433  = t163 * t163 * t53 * t * t;
//	const double t1103 = t433 * t433 * t163 * t53 * t;
//
//	return Dcr * (1 + b1*t13 + b2*t23 + b3*t53 +
//		b4*t163 + b5*t433 + b6*t1103);
//}
//
///// Calculate the saturated vapour density of water in units of kg/m3
//double CalculateSaturatedVapourDensityWater(double T)
//{
//	const double Tcr = 647.096; // critical temperature of water in units of kelvin
//	const double Dcr = 322.0;   // critical density of water in units of kg/m3
//
//	const double c1 = -2.03150240;
//	const double c2 = -2.68302940;
//	const double c3 = -5.38626492;
//	const double c4 = -17.2991605;
//	const double c5 = -44.7586581;
//	const double c6 = -63.9201063;
//
//	const double t    = 1 - T/Tcr;
//	const double t16  = std::pow(t, 1./6);
//	const double t26  = t16 * t16;
//	const double t46  = t26 * t26;
//	const double t86  = t46 * t46;
//	const double t186 = t86 * t86 * t26;
//	const double t376 = t186 * t186 * t16;
//	const double t716 = t376 * t186 * t86 * t86;
//
//	return Dcr * std::exp(Tcr/T * (c1*t26 + c2*t46 +
//		c3*t86 + c4*t186 + c5*t376 + c6*t716));
//}
//
///// Calculate the saturated pressure of carbon dioxide in units of pascal
//double CalculateSaturatedPressureCarbonDioxide(double T)
//{
//	const double Tcr = 304.1282;   // critical temperature of carbon dioxide in units of kelvin
//	const double Pcr = 7.3773e+06; // critical pressure of carbon dioxide in units of pascal
//
//	const double a1 = -7.0602087;
//	const double a2 =  1.9391218;
//	const double a3 = -1.6463597;
//	const double a4 = -3.2995634;
//
//	const double tau  = 1 - T/Tcr;
//	const double tau1 = tau;
//	const double tau2 = std::pow(tau1, 1.5);
//	const double tau3 = tau1 * tau1;
//	const double tau4 = tau3 * tau3;
//
//	return Pcr * std::exp(Tcr/T*(a1*tau1 + a2*tau2 + a3*tau3 + a4*tau4));
//}
//
///// Calculate the saturated liquid density of carbon dioxide in units of kg/m3
//double CalculateSaturatedLiquidDensityCarbonDioxide(double T)
//{
//	const double Tcr = 304.1282; // critical temperature of carbon dioxide in units of kelvin
//	const double Dcr = 467.6;    // critical density of carbon dioxide in units of kg/m3
//
//	const double a1 =  1.92451080;
//	const double a2 = -0.62385555;
//	const double a3 = -0.32731127;
//	const double a4 =  0.39245142;
//
//	const double tau  = 1 - T/Tcr;
//	const double tau1 = std::pow(tau, 0.34);
//	const double tau2 = std::pow(tau, 0.50);
//	const double tau3 = std::pow(tau, 10.0/6);
//	const double tau4 = std::pow(tau, 11.0/6);
//
//	return Dcr * std::exp(a1*tau1 + a2*tau2 + a3*tau3 + a4*tau4);
//}
//
///// Calculate the saturated vapour density of carbon dioxide in units of kg/m3
//double CalculateSaturatedVapourDensityCarbonDioxide(double T)
//{
//	const double Tcr = 304.1282; // critical temperature of carbon dioxide in units of kelvin
//	const double Dcr = 467.6;    // critical density of carbon dioxide in units of kg/m3
//
//	const double a1 = -1.70748790;
//	const double a2 = -0.82274670;
//	const double a3 = -4.60085490;
//	const double a4 = -10.1111780;
//	const double a5 = -29.7422520;
//
//	const double tau  = 1 - T/Tcr;
//	const double tau1 = std::pow(tau, 0.34);
//	const double tau2 = std::pow(tau, 0.50);
//	const double tau3 = tau;
//	const double tau4 = std::pow(tau, 7.0/4);
//	const double tau5 = tau4 * tau4;
//
//	return Dcr * std::exp(a1*tau1 + a2*tau2 + a3*tau3 + a4*tau4 + a5*tau5);
//}
//
///// Calculate an approximated molar volume of water to use as initial guess
//double CalculateApproximateVolumeWater(double T, double P)
//{
//	if(T > TcrH2O)
//	{
//		return mwH2O/DcrH2O;
//	}
//	else
//	{
//		const double Psat = CalculateSaturatedPressureWater(T);
//
//		if(P <= Psat)
//			return mwH2O/CalculateSaturatedVapourDensityWater(T);
//		else
//			return mwH2O/CalculateSaturatedLiquidDensityWater(T);
//	}
//}
//
///// Calculate an approximated molar volume of carbon dioxide to use as initial guess
//double CalculateApproximateVolumeCarbonDioxide(double T, double P)
//{
//	if(T > TcrCO2)
//	{
//		return mwCO2/DcrCO2;
//	}
//	else
//	{
//		const double Psat = CalculateSaturatedPressureCarbonDioxide(T);
//
//		if(P <= Psat)
//			return mwCO2/CalculateSaturatedVapourDensityCarbonDioxide(T);
//		else
//			return mwCO2/CalculateSaturatedLiquidDensityCarbonDioxide(T);
//	}
//}
//
///// Calculate an approximated molar volume of carbon dioxide to use as initial guess
//double CalculateApproximateVolumeMethane(double T, double P)
//{
//	return mwCH4/DcrCH4;
//}
//
//double CalculateReducedVolume(double Tr, double Pr, double Vr0, const double a[])
//{
//	const double a1  = a[0];
//	const double a2  = a[1];
//	const double a3  = a[2];
//	const double a4  = a[3];
//	const double a5  = a[4];
//	const double a6  = a[5];
//	const double a7  = a[6];
//	const double a8  = a[7];
//	const double a9  = a[8];
//	const double a10 = a[9];
//	const double a11 = a[10];
//	const double a12 = a[11];
//	const double a13 = a[12];
//	const double a14 = a[13];
//	const double a15 = a[14];
//
//	unsigned maxiters = 1000;
//	double tolerance = 1.0e-8;
//
//	double Vr = Vr0;
//
//	for(unsigned i = 0; i < maxiters; ++i)
//	{
//		const double iVr  = 1/Vr;
//		const double iVr2 = iVr*iVr;
//		const double iVr3 = iVr*iVr2;
//		const double iVr4 = iVr*iVr3;
//		const double iVr5 = iVr*iVr4;
//		const double iVr6 = iVr*iVr5;
//
//		const double iTr  = 1/Tr;
//		const double iTr2 = iTr*iTr;
//		const double iTr3 = iTr*iTr2;
//
//		const double B =  a1 +  a2*iTr2 +  a3*iTr3;
//		const double C =  a4 +  a5*iTr2 +  a6*iTr3;
//		const double D =  a7 +  a8*iTr2 +  a9*iTr3;
//		const double E = a10 + a11*iTr2 + a12*iTr3;
//		const double F = a13*iTr3;
//		const double G = std::exp(-a15*iVr2);
//
//		const double f = 1 + B*iVr + C*iVr2 + D*iVr4 + E*iVr5 +
//			F*iVr2*(a14 + a15*iVr2)*G - Pr*Vr*iTr;
//
//		const double df = -B*iVr2 - 2*C*iVr3 - 4*D*iVr5 - 5*E*iVr6 -
//			2*F*iVr3*((a14 + a15*iVr2)*(1 - a15*iVr2) + a15/iVr2)*G - Pr*iTr;
//
//		const double dVr = -f/df;
//
//		const double omega = 0.1;
//
//		Vr = Vr + omega*dVr;
//
////		Vr = (Vr > -dVr) ? Vr + omega*dVr : Tr/Pr*(f + Pr*Vr*iTr);
//
////		const double rhs = 1 + B*iVr + C*iVr2 + D*iVr4 + E*iVr5 +
////			F*iVr2*(a14 + a15*iVr2)*G;
////
////		Vr = (1-omega)*Vr + omega*Tr/Pr*rhs;
//
//		std::cout << i << " " << Vr << std::endl;
//
//		if(std::abs(f) < tolerance)
//			return Vr;
//	}
//
//	throw std::runtime_error("ConvergenceError: Could not calculate the "
//		"reduced volume using the Duan et al. (1992) equations of state.");
//}
//
//double CalculateFugacityCoefficient(double Tr, double Pr, double Vr, const double a[])
//{
//	const double a1  = a[0];
//	const double a2  = a[1];
//	const double a3  = a[2];
//	const double a4  = a[3];
//	const double a5  = a[4];
//	const double a6  = a[5];
//	const double a7  = a[6];
//	const double a8  = a[7];
//	const double a9  = a[8];
//	const double a10 = a[9];
//	const double a11 = a[10];
//	const double a12 = a[11];
//	const double a13 = a[12];
//	const double a14 = a[13];
//	const double a15 = a[14];
//
//	const double iVr  = 1/Vr;
//	const double iVr2 = iVr*iVr;
//	const double iVr3 = iVr*iVr2;
//	const double iVr4 = iVr*iVr3;
//	const double iVr5 = iVr*iVr4;
//
//	const double iTr  = 1/Tr;
//	const double iTr2 = iTr*iTr;
//	const double iTr3 = iTr*iTr2;
//
//	const double B =  a1 +  a2*iTr2 +  a3*iTr3;
//	const double C =  a4 +  a5*iTr2 +  a6*iTr3;
//	const double D =  a7 +  a8*iTr2 +  a9*iTr3;
//	const double E = a10 + a11*iTr2 + a12*iTr3;
//	const double F = a13*iTr3;
//	const double G = F/(2*a15)*(a14 + 1 - (a14 + 1 + a15*iVr2)*std::exp(-a15*iVr2));
//	const double Z = Pr*Vr*iTr;
//
//	return std::exp(Z - 1 - std::log(Z) + B*iVr + 0.50*C*iVr2 +
//		0.25*D*iVr4 + 0.20*E*iVr5 + G);
//}
//
//GaseousActivityDuanMollerWeare::GaseousActivityDuanMollerWeare(
//	const std::string& gas, const std::vector<std::string>& gases)
//: GaseousActivity(gas, gases)
//{
//	     if(gas == "CO2(g)") gas_type = CO2;
//	else if(gas == "H2O(g)") gas_type = H2O;
//	else if(gas == "CH4(g)") gas_type = CH4;
//	else throw std::runtime_error("UnknownError: The provided gas is not "
//		"supported by the Duan et al. (1992) model.");
//}
//
//double GaseousActivityDuanMollerWeare::ActivityCoefficient(
//	double T, double P, const Vector& n) const
//{
//	const double Tcr[] = {TcrCO2, TcrH2O, TcrCH4};
//	const double Pcr[] = {PcrCO2, PcrH2O, PcrCH4};
//	const double*  a[] = {aCO2,   aH2O,   aCH4  };
//
//	const double R  = 8.314467;
//	const double Tr = T/Tcr[gas_type];
//	const double Pr = P/Pcr[gas_type];
//	const double Vc = R*Tcr[gas_type]/Pcr[gas_type];
//
//	double Vr0;
//	switch(gas_type)
//	{
//	case CO2: Vr0 = CalculateApproximateVolumeCarbonDioxide(T, P)/Vc; break;
//	case H2O: Vr0 = CalculateApproximateVolumeWater(T, P)/Vc; break;
//	case CH4: Vr0 = CalculateApproximateVolumeMethane(T, P)/Vc; break;
//	}
//
//	const double Vr = CalculateReducedVolume(Tr, Pr, Vr0, a[gas_type]);
//
//	return CalculateFugacityCoefficient(Tr, Pr, Vr, a[gas_type]);
//}
//
//std::string GaseousActivityDuanMollerWeare::Info() const
//{
//	return "This is the Duan et al. (1992) activity model for CO2(g),"
//		"H2O(g) or CH4(g).";
//}
//
//} // namespace Reaktor
